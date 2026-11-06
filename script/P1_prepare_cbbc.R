# ==============================================================================
# P1_prepare_cbbc.R
#
# PREPARATION STEP 1 OF 4 (no inferences; run once, output cached)
#
# Builds the Canadian Breeding Bird Census (CBBC) data used throughout the
# comparison report (cbbc_bam_comparison.Rmd):
#   1. Reads the PLOTS and SPECIES tables.
#   2. Parses plot coordinates (DD.MM.SS text) and checks them against GRIDREF.
#   3. Applies census-quality filters to each census record.
#   4. Computes species densities (territories / km^2) with current species codes.
#   5. Builds one row per retained site (coordinates, habitat, plot size, effort).
#
# All species-code rules used by the later scripts and the report (CBBC and BAM
# recodes, split taxa, Traill's Flycatcher, non-territorial species) are defined
# HERE and stored in the output settings, so every step treats species the same.
#
# Background (CBBC database documentation; Kennedy et al. 1999, CWS Technical
# Report 342):
#   - LAT_N / LONG_W are "of a point on the plot, OR of the nearest identifiable
#     geographic entity (often CITYTOWN)"; GRIDREF gives UTM coordinates of a
#     point within the plot.
#   - PLOTDENS = territories on the plot ("+" = part of a territory, coded 0.1;
#     -999 = not given). DENSKM2 = territories / 100 ha, but "+" was entered as
#     0.1 without standardizing, and some values were rounded.
#   - Species lists hold breeding (territorial) species only; TOTSPP gives the
#     number of species the list should contain.
#   - Spot-mapping censuses territorial, non-colonial species; it cannot census
#     ducks, geese, seabirds, colonial or flocking species as territories.
#
# OUTPUT (OUT_PATH), a list:
#   cbbc_qa         every census record with coordinates, flags and one keep_*
#                   column per filter (shows exactly why a record was dropped)
#   cbbc_long       retained records x species densities (territories / km^2)
#   cbbc_sites      one row per retained site
#   traill_records  retained records listing Traill's Flycatcher
#   attr "settings" filters and species rules
# ==============================================================================

library(tidyverse)
library(readxl)
library(sf)       # UTM -> lat/lon for the GRIDREF check

# ------------------------------------------------------------------------------
# Paths
# ------------------------------------------------------------------------------

ROOT              <- "C:/Users/IlesD/OneDrive - EC-EC/Iles"
CBBC_TABLE_DIR    <- file.path(ROOT, "Projects/Landbirds/spot-mapping-vs-point-counts/data/tables")
CBBC_PLOTS_PATH   <- file.path(CBBC_TABLE_DIR, "PLOTS.xlsx")
CBBC_SPECIES_PATH <- file.path(CBBC_TABLE_DIR, "SPECIES.xlsx")

# Analysis folder inside the project. here::here() starts from the project root
# (the folder holding Spot_Mapping_vs_Point_Counts.Rproj), so paths resolve the
# same whether the script is sourced from script/ or run from the project root.
ANALYSIS_DIR <- here::here()

OUT_DIR  <- file.path(ANALYSIS_DIR, "output")
OUT_PATH <- file.path(OUT_DIR, "P1_cbbc_prepared.rds")

# ------------------------------------------------------------------------------
# Constants: coordinates
# ------------------------------------------------------------------------------

TEXT_COLS      <- c("LAT_N", "LONG_W", "GRIDREF")   # read as text: keeps DD.MM.SS and
                                                    #   stops grid refs being coerced to numbers
GUESS_MAX      <- 1e5
LAT_RANGE      <- c(41, 84)        # plausible Canadian extent; parsed coordinates outside -> NA
LON_RANGE      <- c(-142, -52)
EARTH_RADIUS_M <- 6371008.8        # mean Earth radius for haversine distances

# GRIDREF: Canada spans UTM zones 7-22; NAD83 / UTM zone N = EPSG 26900 + N.
# Historical sources are likely NAD27, but the ~200 m datum shift is negligible
# for a check meant to catch kilometre-scale errors.
UTM_ZONES     <- 7:22
UTM_EPSG_BASE <- 26900

# Records whose LAT/LONG disagree with GRIDREF by more than this may be placed
# at a nearby town rather than the plot, so nearby BAM surveys and raster cells
# would describe the wrong landscape. Records without a parsable GRIDREF are
# kept. Optionally also drop coordinates given only to the minute (+/- ~1.8 km).
COORD_MAX_DISAGREE_M <- 1000
REQUIRE_SECONDS      <- FALSE

# ------------------------------------------------------------------------------
# Constants: census-record filters
# ------------------------------------------------------------------------------

# Species-list completeness. A record is kept only if its species rows reach
# this fraction of TOTSPP: zero-filling a partial list would bias densities low
# (CWS's earlier Preliminary Catalogue gave only the 5 commonest species).
LIST_COMPLETENESS_MIN <- 0.8

# Census effort. Spot-mapping needs repeated visits to separate territories;
# the international standard (Robbins 1970) is ~8-10. Records with VISITS
# missing are dropped: census quality cannot be confirmed.
MIN_VISITS <- 6

# Plot size (ha). Records with PLOTSIZE missing are dropped (density cannot be
# computed from PLOTDENS). < 10 ha: edge effects inflate densities. > 100 ha:
# Robbins (1970) recommends ~10-40 ha in forest and up to ~100 ha in open
# habitat; much larger "plots" cannot be spot-mapped in ~10 visits, and are far
# larger than the 1 km area used to match BAM data.
MIN_PLOT_HA <- 10
MAX_PLOT_HA <- 100

# Census year. No filter by default; raise (e.g. to 1965) as a sensitivity check.
MIN_CENSUS_YEAR <- 1975

# Habitat (Erskine codes in HABCODE1). Early successional plots have almost
# certainly changed since the census, and urban plots are poorly represented by
# point counts.
#   44x clear-cut <= 7 y; 434 old field; 1xy / 2xy with y = 3 or 4: young or
#   disturbed forest; 315y with y = 3 or 4: young or disturbed swamp forest
SUCCESSIONAL_PATTERN <- "^(44|434|[12]\\d[34]|315[34])"
EXCLUDE_SUCCESSIONAL <- TRUE
EXCLUDE_URBAN        <- TRUE
HAB_MAJOR_LABELS     <- c(`1` = "Broad-leafed forest", `2` = "Conifer / mixed forest",
                          `3` = "Wetland", `4` = "Open land", `5` = "Urban")

# "+" (part of a territory on the plot) is coded 0.1 in PLOTDENS and DENSKM2
PARTIAL_TERRITORY <- 0.1

# ------------------------------------------------------------------------------
# Constants: species rules (shared by all later scripts via the settings)
# ------------------------------------------------------------------------------

# Old or pre-lump CBBC codes -> current 4-letter codes used by NA-POPS and BAM.
# Forms recoded to the same species within a census are summed.
CBBC_RECODE <- c(MYWA = "YRWA", AUWA = "YRWA", GRJA = "CAJA", GRAJ = "CAJA",
                 YSFL = "NOFL", RSFL = "NOFL", SCJU = "DEJU", ORJU = "DEJU",
                 YWAR = "YEWA", BTGW = "BTNW", RODO = "ROPI", YPWA = "PAWA",
                 RPHE = "RNEP", WEFL = "PSFL")

# BAM codes that differ from the analysis codes
BAM_RECODE <- c(NHWR = "HOWR", WEFL = "PSFL")   # EDIT: confirm BAM's House Wren code

# Taxa split after the CBBC era, assigned from plot location (approximate)
TOWHEE_SPLIT_LON <- -101           # Rufous-sided Towhee (RSTO): west = SPTO, east = EATO
VIREO_SPLIT_LON  <- -114           # Solitary Vireo (SOVI): west of this and south of
VIREO_SPLIT_LAT  <- 56             #   VIREO_SPLIT_LAT = CAVI, otherwise BHVI

# Traill's Flycatcher (pre-1973) = Alder + Willow Flycatcher, which overlap in
# southern Canada; in records listing it, ALFL and WIFL are unknown (NA), not 0
TRAILL_CODE    <- "TRFL"
TRAILL_SPECIES <- c("ALFL", "WIFL")

# Species spot-mapping cannot census as territories: kept in outputs but
# flagged (territorial = FALSE) and summarized separately in the report
EXCLUDE_SPECIES <- c(
  "BHCO",                                          # brood parasite: holds no nesting territory
  "COGR", "EUST", "HOSP", "ROPI",                  # colonial / semi-colonial nesters
  "RWBL", "BRBL", "YHBL",                          # polygynous or colonial blackbirds
  "CEDW", "AMGO", "PISI", "EVGR", "RECR", "WWCR",  # flocking, weakly territorial or nomadic
  "BARS", "CLSW", "TRES", "NRWS", "BANS", "PUMA",  # swallows / martins: aerial foragers,
  "CHSW",                                          #   often colonial; counted in flight
  "AMCR", "CORA"                                   # territories far larger than census plots
)

# ------------------------------------------------------------------------------
# Helper functions
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

# GRIDREF holds UTM coordinates in free text. Handles the common
# "zone [letter] easting northing" layout (e.g. "17 630500 4830200" or
# "17T 630500E 4830200N"); anything else is left NA rather than guessed.
parse_utm <- function(x) {
  m <- str_match(string = x, pattern = "(\\d{1,2})\\s*[A-Za-z]?\\D+(\\d{6})(?:\\.\\d+)?\\D+(\\d{7})")
  tibble(utm_zone = as.integer(m[, 2]), utm_e = as.numeric(m[, 3]), utm_n = as.numeric(m[, 4]))
}

# Most frequent value (first in case of ties); NA if every value is missing
mode_value <- function(x) {
  tab <- table(x, useNA = "no")
  if (length(tab) == 0) NA_character_ else names(which.max(tab))
}

# ==============================================================================
# 1. Read tables
# ==============================================================================

plot_cols <- names(read_xlsx(path = CBBC_PLOTS_PATH, n_max = 0))

plots <- read_xlsx(path      = CBBC_PLOTS_PATH,
                   col_types = if_else(condition = plot_cols %in% TEXT_COLS, true = "text", false = "guess"),
                   guess_max = GUESS_MAX)

species_tbl <- read_xlsx(path = CBBC_SPECIES_PATH, guess_max = GUESS_MAX)

# ==============================================================================
# 2. Coordinates and coordinate quality
# ==============================================================================

# Decimal degrees (west longitudes negative); values outside Canada -> NA
plots <- plots |>
  mutate(
    latitude         =  parse_dms(LAT_N),
    longitude        = -abs(parse_dms(LONG_W)),
    coord_ambiguous  = dms_ambiguous(LAT_N) | dms_ambiguous(LONG_W),
    coord_no_seconds = str_count(string = LAT_N, pattern = fixed(".")) < 2 |
                       str_count(string = LONG_W, pattern = fixed(".")) < 2,
    coord_valid      = between(latitude, LAT_RANGE[1], LAT_RANGE[2]) &
                       between(longitude, LON_RANGE[1], LON_RANGE[2]),
    across(.cols = c(latitude, longitude),
           .fns  = ~ if_else(condition = coord_valid %in% TRUE, true = .x, false = NA_real_))
  )

# GRIDREF -> lat/lon, projecting each UTM zone separately (each zone has its own
# CRS). Starting from an empty template keeps the join working if none parse.
utm <- bind_cols(select(plots, RECNO), parse_utm(x = plots$GRIDREF)) |>
  filter(utm_zone %in% UTM_ZONES)

gridref_points <- utm |>
  split(f = utm$utm_zone) |>
  map(.f = function(d) {
    xy <- sf::st_as_sf(x = d, coords = c("utm_e", "utm_n"), crs = UTM_EPSG_BASE + d$utm_zone[1]) |>
      sf::st_transform(crs = 4326) |>
      sf::st_coordinates()
    tibble(RECNO = d$RECNO, grid_lon = xy[, "X"], grid_lat = xy[, "Y"])
  }) |>
  list_rbind(ptype = tibble(RECNO = plots$RECNO[0], grid_lon = numeric(), grid_lat = numeric()))

plots <- plots |>
  left_join(y = gridref_points, by = "RECNO") |>
  mutate(gridref_dist_m = haversine_m(lat1 = latitude, lon1 = longitude, lat2 = grid_lat, lon2 = grid_lon))

cat("Coordinates (records): unusable = ", sum(is.na(plots$latitude)),
    " | ambiguous DD.MM = ", sum(plots$coord_ambiguous),
    " | no seconds = ", sum(plots$coord_no_seconds, na.rm = TRUE),
    " | GRIDREF present = ", sum(!is.na(plots$GRIDREF)),
    " | GRIDREF parsed = ", sum(!is.na(plots$gridref_dist_m)), "\n", sep = "")
cat("Unparsed GRIDREF examples (check whether another format is common):\n")
print(head(unique(plots$GRIDREF[!is.na(plots$GRIDREF) & is.na(plots$grid_lat)]), 15))

# ==============================================================================
# 3. Census-record filters
# ==============================================================================

# One keep_* column per filter, so cbbc_qa shows exactly why a record was
# dropped. A record is retained only if every filter passes (NA = drop).
cbbc_qa <- plots |>
  left_join(y = count(species_tbl, RECNO, name = "n_species_rows"), by = "RECNO") |>
  mutate(
    hab_code     = as.character(HABCODE1),
    hab_major    = coalesce(unname(HAB_MAJOR_LABELS[str_sub(string = hab_code, start = 1, end = 1)]), "Unknown"),
    successional = coalesce(str_detect(string = hab_code, pattern = SUCCESSIONAL_PATTERN), FALSE),
    keep_coords  = !is.na(latitude),
    keep_gridref = (is.na(gridref_dist_m) | gridref_dist_m <= COORD_MAX_DISAGREE_M) &
                   !(REQUIRE_SECONDS & coord_no_seconds %in% TRUE),
    keep_list    = !is.na(n_species_rows) & (is.na(TOTSPP) | n_species_rows >= LIST_COMPLETENESS_MIN * TOTSPP),
    keep_visits  = !is.na(VISITS) & VISITS >= MIN_VISITS,
    keep_size    = !is.na(PLOTSIZE) & between(PLOTSIZE, MIN_PLOT_HA, MAX_PLOT_HA),
    keep_year    = !is.na(YEAR) & YEAR >= MIN_CENSUS_YEAR,
    keep_habitat = !(EXCLUDE_SUCCESSIONAL & successional) & !(EXCLUDE_URBAN & hab_major == "Urban"),
    keep         = keep_coords & keep_gridref & keep_list & keep_visits & keep_size & keep_year & keep_habitat
  )

cat("\nCensus records failing each filter (a record can fail several):\n")
cbbc_qa |>
  summarise(n_records = n(), across(.cols = starts_with("keep_"), .fns = ~ sum(!.x, na.rm = TRUE))) |>
  print(width = Inf)
cat("Records retained: ", sum(cbbc_qa$keep), " of ", nrow(cbbc_qa),
    " | sites retained: ", n_distinct(cbbc_qa$SITEID[cbbc_qa$keep]), "\n", sep = "")

records <- cbbc_qa |> filter(keep)

# ==============================================================================
# 4. Densities per retained census record (territories / km^2)
# ==============================================================================

cbbc_long <- species_tbl |>
  inner_join(y = select(records, RECNO, SITEID, YEAR, PLOTSIZE, latitude, longitude), by = "RECNO") |>
  mutate(
    # Current code: location-based for split taxa, lookup for renamed / lumped
    # codes, otherwise the CBBC code as is
    species  = case_when(
      CODE == "RSTO" ~ if_else(condition = longitude < TOWHEE_SPLIT_LON, true = "SPTO", false = "EATO"),
      CODE == "SOVI" ~ if_else(condition = longitude < VIREO_SPLIT_LON & latitude < VIREO_SPLIT_LAT,
                               true = "CAVI", false = "BHVI"),
      .default = coalesce(unname(CBBC_RECODE[CODE]), CODE)
    ),
    # Density, in order of preference:
    #  1. PLOTDENS / PLOTSIZE (ha) x 100: most precise, and handles "+" as 0.1
    #     territory on the plot (PLOTDENS was added to the database for this)
    #  2. DENSKM2 = 0.1 is an unstandardized "+": convert as 0.1 territory
    #  3. otherwise DENSKM2 as given (may be rounded in some sources)
    dens_km2 = case_when(
      PLOTDENS >= 0 & PLOTSIZE > 0                    ~ PLOTDENS / PLOTSIZE * 100,
      near(DENSKM2, PARTIAL_TERRITORY) & PLOTSIZE > 0 ~ PARTIAL_TERRITORY / PLOTSIZE * 100,
      DENSKM2 >= 0                                    ~ DENSKM2,
      .default                                        = NA_real_
    )
  ) |>
  group_by(RECNO, SITEID, YEAR, species) |>          # sum forms recoded to the same species
  summarise(dens_km2 = sum(dens_km2), .groups = "drop") |>
  mutate(territorial = !species %in% EXCLUDE_SPECIES)

traill_records <- species_tbl |>
  filter(CODE == TRAILL_CODE, RECNO %in% records$RECNO) |>
  pull(RECNO) |>
  unique()

# ==============================================================================
# 5. One row per retained site
# ==============================================================================

# Coordinates: the most frequent pair among the site's retained records.
# Habitat: the most frequent major habitat across years. Plot size and visits:
# median and mean across years.
cbbc_sites <- records |>
  add_count(SITEID, latitude, longitude, name = "n_coord") |>
  group_by(SITEID) |>
  summarise(
    PROVINCE         = first(PROVINCE),
    n_census         = n(),
    year_min         = min(YEAR),
    year_max         = max(YEAR),
    year_mean        = mean(YEAR),
    plot_ha          = median(PLOTSIZE),
    visits_mean      = mean(VISITS),
    hab_major        = mode_value(hab_major),
    hab_code         = mode_value(hab_code),
    n_coord_versions = n_distinct(latitude, longitude),
    coord_no_seconds = coord_no_seconds[which.max(n_coord)],
    gridref_dist_m   = suppressWarnings(median(gridref_dist_m, na.rm = TRUE)),
    latitude         = latitude[which.max(n_coord)],
    longitude        = longitude[which.max(n_coord)],
    .groups          = "drop"
  ) |>
  mutate(plot_radius_m = sqrt(plot_ha * 1e4 / pi))   # radius of a circle with the plot's area

cat("\nRetained sites by major habitat:\n")
print(count(cbbc_sites, hab_major, sort = TRUE))

# ==============================================================================
# 6. Save
# ==============================================================================

out <- list(cbbc_qa = cbbc_qa, cbbc_long = cbbc_long, cbbc_sites = cbbc_sites, traill_records = traill_records)
attr(out, "settings") <- list(
  coord_max_disagree_m  = COORD_MAX_DISAGREE_M,
  require_seconds       = REQUIRE_SECONDS,
  list_completeness_min = LIST_COMPLETENESS_MIN,
  min_visits            = MIN_VISITS,
  min_plot_ha           = MIN_PLOT_HA,
  max_plot_ha           = MAX_PLOT_HA,
  min_census_year       = MIN_CENSUS_YEAR,
  successional_pattern  = SUCCESSIONAL_PATTERN,
  exclude_successional  = EXCLUDE_SUCCESSIONAL,
  exclude_urban         = EXCLUDE_URBAN,
  cbbc_density          = "PLOTDENS / PLOTSIZE x 100; DENSKM2 = 0.1 as '+'; else DENSKM2",
  cbbc_recode           = CBBC_RECODE,
  bam_recode            = BAM_RECODE,
  split_taxa            = c(towhee_lon = TOWHEE_SPLIT_LON, vireo_lon = VIREO_SPLIT_LON, vireo_lat = VIREO_SPLIT_LAT),
  traill_code           = TRAILL_CODE,
  traill_species        = TRAILL_SPECIES,
  exclude_species       = EXCLUDE_SPECIES,
  created               = Sys.time()
)

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
saveRDS(object = out, file = OUT_PATH)
message("Saved: ", normalizePath(OUT_PATH))
