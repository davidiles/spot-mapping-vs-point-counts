# ==============================================================================
# 01_match_bam_to_cbbc.R
#
# Loads the CBBC census tables and the BAM point-count database, converts CBBC
# plot coordinates (DD.MM.SS text) to decimal degrees, and finds every BAM
# survey within MAX_MATCH_DIST_M of each CBBC site. Distances are kept so the
# matching radius can be chosen downstream rather than fixed here.
#
# Output (OUT_PATH): named list
#   cbbc_plots    census records (one row per RECNO) with parsed coordinates
#                 and coordinate-quality flags
#   cbbc_sites    one row per SITEID: modal coordinates, plot size, years,
#                 coordinate-quality flags
#   cbbc_species  species x census-record densities (as read)
#   site_survey   CBBC site x BAM survey pairs within MAX_MATCH_DIST_M
#   bam_counts    BAM counts (long) for all matched surveys; absent rows = 0
# Matching settings are stored in attr(out, "settings").
# ==============================================================================

library(tidyverse)
library(readxl)
library(DBI)
library(duckdb)
library(sf)       # UTM -> lat/lon for the GRIDREF coordinate check

# ------------------------------------------------------------------------------
# Paths
# ------------------------------------------------------------------------------

ROOT              <- "C:/Users/IlesD/OneDrive - EC-EC/Iles"
CBBC_TABLE_DIR    <- file.path(ROOT, "Projects/Landbirds/Breeding-Bird-Census/data/tables")
CBBC_PLOTS_PATH   <- file.path(CBBC_TABLE_DIR, "PLOTS.xlsx")
CBBC_SPECIES_PATH <- file.path(CBBC_TABLE_DIR, "SPECIES.xlsx")

BAM_DB_PATH <- "C:/Users/IlesD/Documents/Data_Local/BAM_v5/BAMDataset_duckdb/03_BAMDataset_WT-2026-08-26_EBd-Jun-2026_PRISM-2026-10-09.duckdb"

OUT_DIR  <- "output"
OUT_PATH <- file.path(OUT_DIR, "01_cbbc_bam_matched.rds")

# ------------------------------------------------------------------------------
# Constants
# ------------------------------------------------------------------------------

MAX_MATCH_DIST_M <- 10000                      # keep all surveys within this distance of a site
REPORT_RADII_M   <- c(250, 500, 1000, 2000, 5000, 10000)
VISIT_FILTER_SQL <- "NOT flag_any AND NOT flag_dup"

EARTH_RADIUS_M   <- 6371008.8                  # mean radius for haversine distances
M_PER_DEG_LAT_LO <- 110000                     # conservative (low) value so the bbox pre-filter never clips

COORD_COLS     <- c("LAT_N", "LONG_W", "GRIDREF")   # read as text: preserves DD.MM.SS and
#   stops numeric-looking grid refs being coerced
CBBC_GUESS_MAX <- 1e5
LAT_RANGE      <- c(41, 84)                    # plausible Canadian extent; parsed coords outside -> NA
LON_RANGE      <- c(-142, -52)

# GRIDREF zones: Canada spans UTM zones 7-22; NAD83 / UTM zone N = EPSG 26900 + N.
# Historical sources are likely NAD27, but the ~200 m datum shift is negligible
# for a check whose purpose is to catch kilometre-scale errors.
UTM_ZONES     <- 7:22
UTM_EPSG_BASE <- 26900

# ------------------------------------------------------------------------------
# Helpers
# ------------------------------------------------------------------------------

# "DD.MM.SS" -> decimal degrees; missing seconds treated as 0
parse_dms <- function(x) {
  p  <- str_split_fixed(string = str_trim(x), pattern = fixed("."), n = 3)
  dd <- suppressWarnings(as.numeric(p[, 1]))
  mm <- coalesce(suppressWarnings(as.numeric(p[, 2])), 0)
  ss <- coalesce(suppressWarnings(as.numeric(p[, 3])), 0)
  if_else(condition = mm < 60 & ss < 60, true = dd + mm / 60 + ss / 3600, false = NA_real_)
}

# "45.2" may be an Excel number that lost a trailing zero (45.20 -> 45.2)
dms_ambiguous <- function(x) coalesce(str_detect(string = x, pattern = "^\\s*\\d+\\.\\d\\s*$"), FALSE)

# Great-circle distance (m)
haversine_m <- function(lat1, lon1, lat2, lon2) {
  to_rad <- pi / 180
  a <- sin((lat2 - lat1) * to_rad / 2)^2 +
    cos(lat1 * to_rad) * cos(lat2 * to_rad) * sin((lon2 - lon1) * to_rad / 2)^2
  2 * EARTH_RADIUS_M * asin(pmin(1, sqrt(a)))
}

# GRIDREF holds UTM coordinates in free text. This handles the common
# "zone [letter] easting northing" layout (e.g. "17 630500 4830200" or
# "17T 630500E 4830200N"); anything else (e.g. military grid references tied
# to an NTS map sheet) is left NA rather than guessed.
parse_utm <- function(x) {
  m <- str_match(string = x, pattern = "(\\d{1,2})\\s*[A-Za-z]?\\D+(\\d{6})(?:\\.\\d+)?\\D+(\\d{7})")
  tibble(utm_zone = as.integer(m[, 2]), utm_e = as.numeric(m[, 3]), utm_n = as.numeric(m[, 4]))
}

# ==============================================================================
# 1. CBBC tables
# ==============================================================================

plot_cols <- names(read_xlsx(path = CBBC_PLOTS_PATH, n_max = 0))

cbbc_plots <- read_xlsx(
  path      = CBBC_PLOTS_PATH,
  col_types = if_else(condition = plot_cols %in% COORD_COLS, true = "text", false = "guess"),
  guess_max = CBBC_GUESS_MAX
) |>
  mutate(
    latitude        =  parse_dms(LAT_N),
    longitude       = -abs(parse_dms(LONG_W)),
    coord_ambiguous = dms_ambiguous(LAT_N) | dms_ambiguous(LONG_W),
    coord_valid     = between(latitude, LAT_RANGE[1], LAT_RANGE[2]) &
      between(longitude, LON_RANGE[1], LON_RANGE[2]),
    across(.cols = c(latitude, longitude),
           .fns  = ~ if_else(condition = coord_valid %in% TRUE, true = .x, false = NA_real_))
  )

cbbc_species <- read_xlsx(path = CBBC_SPECIES_PATH, guess_max = CBBC_GUESS_MAX)

cat("CBBC records: ", nrow(cbbc_plots),
    " | unusable coordinates: ", sum(is.na(cbbc_plots$latitude)),
    " | ambiguous DD.MM: ", sum(cbbc_plots$coord_ambiguous), "\n", sep = "")
cbbc_plots |>
  filter(is.na(latitude) | coord_ambiguous) |>
  select(RECNO, SITEID, YEAR, PROVINCE, all_of(COORD_COLS), latitude, longitude, coord_ambiguous) |>
  print(n = Inf)

# ------------------------------------------------------------------------------
# Coordinate quality
#
# Per the CBBC documentation, LAT_N / LONG_W are "of a point on the plot, OR of
# the nearest identifiable geographic entity (often CITYTOWN)", and were looked
# up approximately from maps when the source gave none. Town-based coordinates
# can be kilometres from the plot, which would match BAM surveys from the wrong
# landscape. Two independent signals are recorded:
#   coord_no_seconds  coordinates given only to the minute (+/- ~1.8 km), a
#                     sign of a coarse or map-derived position
#   gridref_dist_m    distance to the point given by GRIDREF (UTM of a point
#                     within the plot, e.g. its NW corner). Large disagreement
#                     means at least one of the two is not on the plot. Small
#                     disagreement supports the coordinates, though both could
#                     still refer to the same town.
# Script 02 uses these flags to drop sites; nothing is dropped here.
# ------------------------------------------------------------------------------

utm <- bind_cols(select(cbbc_plots, RECNO), parse_utm(x = cbbc_plots$GRIDREF)) |>
  filter(utm_zone %in% UTM_ZONES)

# Project each UTM zone separately (each zone has its own CRS). Starting from
# an empty template means the join below still works if nothing parsed.
gridref_points <- utm |>
  split(f = utm$utm_zone) |>
  map(.f = function(d) {
    xy <- sf::st_as_sf(x = d, coords = c("utm_e", "utm_n"), crs = UTM_EPSG_BASE + d$utm_zone[1]) |>
      sf::st_transform(crs = 4326) |>
      sf::st_coordinates()
    tibble(RECNO = d$RECNO, grid_lon = xy[, "X"], grid_lat = xy[, "Y"])
  }) |>
  list_rbind(ptype = tibble(RECNO = cbbc_plots$RECNO[0], grid_lon = numeric(), grid_lat = numeric()))

cbbc_plots <- cbbc_plots |>
  left_join(y = gridref_points, by = "RECNO") |>
  mutate(
    coord_no_seconds = str_count(string = LAT_N, pattern = fixed(".")) < 2 |
      str_count(string = LONG_W, pattern = fixed(".")) < 2,
    gridref_dist_m   = haversine_m(lat1 = latitude, lon1 = longitude, lat2 = grid_lat, lon2 = grid_lon)
  )

cat("\nCoordinate quality (records): no seconds = ", sum(cbbc_plots$coord_no_seconds, na.rm = TRUE),
    " | GRIDREF present = ", sum(!is.na(cbbc_plots$GRIDREF)),
    " | GRIDREF parsed = ", sum(!is.na(cbbc_plots$gridref_dist_m)), "\n", sep = "")
cat("Distance between LAT/LONG and GRIDREF (m):\n")
print(quantile(x = cbbc_plots$gridref_dist_m, probs = c(0, 0.25, 0.5, 0.75, 0.9, 1), na.rm = TRUE))
cat("Unparsed GRIDREF examples (check whether another format is common):\n")
print(head(unique(cbbc_plots$GRIDREF[!is.na(cbbc_plots$GRIDREF) & is.na(cbbc_plots$grid_lat)]), 15))

# One row per site, using the most frequent coordinate pair across census years
cbbc_sites <- cbbc_plots |>
  filter(!is.na(latitude)) |>
  add_count(SITEID, latitude, longitude, name = "n_coord") |>
  group_by(SITEID) |>
  summarise(
    PROVINCE         = first(PROVINCE),
    n_records        = n(),
    year_min         = min(YEAR, na.rm = TRUE),
    year_max         = max(YEAR, na.rm = TRUE),
    plot_ha          = median(PLOTSIZE, na.rm = TRUE),
    n_coord_versions = n_distinct(latitude, longitude),
    coord_ambiguous  = coord_ambiguous[which.max(n_coord)],
    coord_no_seconds = coord_no_seconds[which.max(n_coord)],
    gridref_dist_m   = median(gridref_dist_m, na.rm = TRUE),   # NA if no record has a parsable GRIDREF
    latitude         = latitude[which.max(n_coord)],
    longitude        = longitude[which.max(n_coord)],
    .groups          = "drop"
  ) |>
  mutate(plot_radius_m = sqrt(plot_ha * 1e4 / pi))   # radius of a circle with the plot's area

cat("\nCBBC sites with usable coordinates: ", nrow(cbbc_sites),
    " | sites whose coordinates differ among years: ", sum(cbbc_sites$n_coord_versions > 1), "\n", sep = "")

# ==============================================================================
# 2. Match BAM surveys to CBBC sites (in duckdb)
# ==============================================================================

con <- dbConnect(drv = duckdb::duckdb())
dbExecute(conn      = con,
          statement = sprintf("ATTACH '%s' AS bam (READ_ONLY)",
                              normalizePath(path = BAM_DB_PATH, winslash = "/")))

# Bounding box per site for a cheap range join; exact distances filtered after
site_bbox <- cbbc_sites |>
  mutate(
    dlat = MAX_MATCH_DIST_M / M_PER_DEG_LAT_LO,
    dlon = dlat / cos(latitude * pi / 180)
  ) |>
  transmute(SITEID, latitude, longitude,
            lat_min = latitude - dlat,  lat_max = latitude + dlat,
            lon_min = longitude - dlon, lon_max = longitude + dlon)

dbWriteTable(conn = con, name = "cbbc_sites", value = site_bbox, temporary = TRUE, overwrite = TRUE)

match_sql <- sprintf("
  WITH v AS (
    SELECT survey_id, date_time, sensor, method, duration, distance, longitude, latitude
    FROM bam.main.visit
    WHERE %s
  ),
  pairs AS (
    SELECT s.SITEID, v.*,
           2 * %f * asin(least(1, sqrt(
             pow(sin(radians(v.latitude - s.latitude) / 2), 2) +
             cos(radians(s.latitude)) * cos(radians(v.latitude)) *
             pow(sin(radians(v.longitude - s.longitude) / 2), 2)))) AS dist_m
    FROM cbbc_sites s
    JOIN v ON v.latitude  BETWEEN s.lat_min AND s.lat_max
          AND v.longitude BETWEEN s.lon_min AND s.lon_max
  )
  SELECT * FROM pairs WHERE dist_m <= %f",
                     VISIT_FILTER_SQL, EARTH_RADIUS_M, MAX_MATCH_DIST_M)

site_survey <- dbGetQuery(conn = con, statement = match_sql) |>
  as_tibble() |>
  left_join(y = select(cbbc_sites, SITEID, plot_radius_m), by = "SITEID") |>
  mutate(
    year        = year(date_time),
    doy         = yday(date_time),
    within_plot = dist_m <= plot_radius_m
  ) |>
  arrange(SITEID, dist_m)

# Counts for matched surveys only; flag_err_count kept so outliers can be dropped downstream
dbWriteTable(conn = con, name = "matched_ids", value = distinct(site_survey, survey_id),
             temporary = TRUE, overwrite = TRUE)

bam_counts <- tbl(src = con, from = I("bam.main.bird")) |>
  semi_join(y = tbl(src = con, from = "matched_ids"), by = "survey_id") |>
  collect()

dbDisconnect(conn = con, shutdown = TRUE)

# ==============================================================================
# 3. Coverage summary
# ==============================================================================

coverage <- tibble(radius_m = REPORT_RADII_M) |>
  mutate(
    n_sites = map_int(.x = radius_m, .f = ~ n_distinct(site_survey$SITEID[site_survey$dist_m <= .x])),
    n_pairs = map_int(.x = radius_m, .f = ~ sum(site_survey$dist_m <= .x))
  )

cat("\nSites with >= 1 BAM survey, by radius (", nrow(cbbc_sites), " sites total):\n", sep = "")
print(coverage)

cat("\nSites with >= 1 survey inside the plot footprint: ",
    n_distinct(site_survey$SITEID[site_survey$within_plot]), "\n", sep = "")

cat("\nSite x survey pairs by distance class and sensor:\n")
site_survey |>
  mutate(dist_class = cut(x = dist_m, breaks = c(0, REPORT_RADII_M), include.lowest = TRUE)) |>
  count(dist_class, sensor) |>
  pivot_wider(names_from = sensor, values_from = n, values_fill = 0) |>
  print()

cat("\nMatched surveys by year (decade):\n")
print(table(floor(site_survey$year / 10) * 10, useNA = "ifany"))

# ==============================================================================
# 4. Save
# ==============================================================================

out <- list(
  cbbc_plots   = cbbc_plots,
  cbbc_sites   = cbbc_sites,
  cbbc_species = cbbc_species,
  site_survey  = site_survey,
  bam_counts   = bam_counts
)
attr(out, "settings") <- list(
  bam_db           = basename(BAM_DB_PATH),
  visit_filter     = VISIT_FILTER_SQL,
  max_match_dist_m = MAX_MATCH_DIST_M,
  distance_method  = "haversine, R = 6371008.8 m",
  site_coords      = "modal DD.MM.SS pair across census years",
  gridref_check    = "GRIDREF parsed as 'zone easting northing', projected as NAD83 / UTM",
  created          = Sys.time()
)

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
saveRDS(object = out, file = OUT_PATH)
message("Saved: ", normalizePath(OUT_PATH))