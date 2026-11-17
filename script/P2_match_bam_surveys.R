# ==============================================================================
# P2_match_bam_surveys.R
#
# PREPARATION STEP 2 OF 5 (no inferences; run once, output cached)
#
# Finds BAM point-count surveys near each retained CBBC site (from P1), applies
# the survey filters, computes time since sunrise, excludes ALL CBBC SITEIDs
# sharing identical recorded site coordinates, assigns each remaining survey
# to its nearest unambiguous CBBC site, and attaches the counts.
#
# Matching radius: 1 km. Not smaller, because CBBC coordinates are often given
# only to the minute or refer to a nearby town, and tighter radii leave too few
# surveys per site. Not larger: a paired test on the same sites found corrected
# densities ~11% lower at 1-2 km than at 0-1 km, with mature-forest species
# declining and edge / early-successional species increasing, i.e. points
# further out sample different habitat from the plot.
#
# Survey filters, each recorded in survey_funnel so the report can show how
# the sample shrinks at every step:
#   - human point counts only (sensor PC): NA-POPS availability and EDR are not
#     calibrated for ARU recordings
#   - standardized point counts only (method PC): eBird checklists excluded
#   - breeding season (DOY_RANGE)
#   - valid duration and distance limit
#   - morning counts within the NA-POPS training range (TSSR_RANGE)
#
# OUTPUT (OUT_PATH), a list:
#   site_survey    retained BAM surveys, each assigned ONCE to its nearest
#                  uniquely located CBBC site (distance, date, design, etc.)
#   excluded_matched_sites  ambiguous CBBC sites with surveys before exclusion
#   exclusion_audit sites/surveys affected by the shared-coordinate rule
#   bam_counts     counts for retained surveys (detections only; a missing
#                  survey x species row is a zero), analysis species codes
#   survey_funnel  sites and surveys remaining after each filter
#   attr "settings"
# ==============================================================================

library(tidyverse)
library(DBI)
library(duckdb)
library(suncalc)  # sunrise times
library(lutz)     # time zone of each survey location

# ------------------------------------------------------------------------------
# Paths
# ------------------------------------------------------------------------------

# Analysis folder inside the project (here::here() starts from the folder
# holding Spot_Mapping_vs_Point_Counts.Rproj)
ANALYSIS_DIR <-  here::here()

P1_PATH     <- file.path(ANALYSIS_DIR, "output", "P1_cbbc_prepared.rds")
BAM_DB_PATH <- "C:/Users/IlesD/Documents/Data_Local/BAM_v5/BAMDataset_duckdb/03_BAMDataset_WT-2026-08-26_EBd-Jun-2026_PRISM-2026-10-09.duckdb"

OUT_DIR  <- file.path(ANALYSIS_DIR, "output")
OUT_PATH <- file.path(OUT_DIR, "P2_bam_matched.rds")

# ------------------------------------------------------------------------------
# Constants
# ------------------------------------------------------------------------------

MATCH_RADIUS_M   <- 2000           # see header
VISIT_FILTER_SQL <- "NOT flag_any AND NOT flag_dup"   # BAM's own error / duplicate flags

SENSORS    <- "PC"                 # human point counts
METHODS    <- "PC"                 # standardized point counts (eBird excluded)
DOY_RANGE  <- c(145, 196)          # 25 May - 15 Jul
TSSR_RANGE <- c(-1, 6)             # hours since sunrise; drops evening counts

EARTH_RADIUS_M   <- 6371008.8      # mean Earth radius for haversine distances
M_PER_DEG_LAT_LO <- 110000         # conservative (low) value so the bounding-box pre-filter never clips

# ------------------------------------------------------------------------------
# Helper functions
# ------------------------------------------------------------------------------

clock_hours <- function(t) hour(t) + minute(t) / 60 + second(t) / 3600

# Local clock-time sunrise (decimal hours) for locations in one time zone;
# suncalc returns sunrise in `tz`, so DST is handled automatically
sunrise_hours <- function(date, lat, lon, tz) {
  sr <- suncalc::getSunlightTimes(data = data.frame(date = date, lat = lat, lon = lon),
                                  keep = "sunrise", tz = tz)$sunrise
  clock_hours(sr)
}

# One row of the survey funnel
funnel_row <- function(df, stage) tibble(stage = stage, n_sites = n_distinct(df$SITEID), n_surveys = n_distinct(df$survey_id))

# ==============================================================================
# 1. CBBC sites and species rules from P1
# ==============================================================================

p1         <- readRDS(file = P1_PATH)
BAM_RECODE <- attr(x = p1, which = "settings")$bam_recode
sites <- p1$cbbc_sites |> select(SITEID, latitude, longitude, n_sites_at_coords, shared_coords)
if (anyNA(sites$shared_coords)) stop("P1 returned missing shared-coordinate flags; re-run P1.")
if (anyDuplicated(sites$SITEID)) stop("P1 cbbc_sites must have one row per SITEID.")

# ==============================================================================
# 2. Match BAM surveys to sites (in duckdb)
# ==============================================================================

con <- dbConnect(drv = duckdb::duckdb())
dbExecute(conn = con, statement = sprintf("ATTACH '%s' AS bam (READ_ONLY)",
                                          normalizePath(path = BAM_DB_PATH, winslash = "/")))

# Bounding box per site for a cheap range join; exact distances filtered after
site_bbox <- sites |>
  mutate(dlat = MATCH_RADIUS_M / M_PER_DEG_LAT_LO,
         dlon = dlat / cos(latitude * pi / 180)) |>
  transmute(SITEID, latitude, longitude,
            lat_min = latitude - dlat,  lat_max = latitude + dlat,
            lon_min = longitude - dlon, lon_max = longitude + dlon)

dbWriteTable(conn = con, name = "cbbc_sites", value = site_bbox, temporary = TRUE, overwrite = TRUE)

match_sql <- sprintf("
  WITH v AS (
    SELECT survey_id, location_id, project_id, organization, date_time, sensor, method,
           duration, distance, longitude, latitude
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
  VISIT_FILTER_SQL, EARTH_RADIUS_M, MATCH_RADIUS_M)

matched <- dbGetQuery(conn = con, statement = match_sql) |>
  as_tibble() |>
  mutate(year = year(date_time), doy = yday(date_time))

# ==============================================================================
# 3. Survey filters, recording the sample left after each
# ==============================================================================

f_sensor   <- matched    |> filter(sensor %in% SENSORS)
f_method   <- f_sensor   |> filter(method %in% METHODS)
f_season   <- f_method   |> filter(between(doy, DOY_RANGE[1], DOY_RANGE[2]))
f_design   <- f_season   |> filter(duration > 0, distance > 0)

# Time since sunrise. date_time is local clock time (a naive TIMESTAMP that R
# reads as UTC), so its clock hours are compared with sunrise in the survey's
# own time zone.
survey_times <- f_design |>
  distinct(survey_id, date_time, latitude, longitude) |>
  mutate(tz = lutz::tz_lookup_coords(lat = latitude, lon = longitude, method = "accurate", warn = FALSE)) |>
  filter(!is.na(tz)) |>
  group_by(tz) |>
  mutate(sunrise_h = sunrise_hours(date = as.Date(date_time), lat = latitude, lon = longitude, tz = first(tz))) |>
  ungroup() |>
  mutate(tssr = clock_hours(date_time) - sunrise_h)

eligible_pairs <- f_design |>
  inner_join(y = select(survey_times, survey_id, tssr), by = "survey_id") |>
  filter(between(tssr, TSSR_RANGE[1], TSSR_RANGE[2])) |>
  arrange(SITEID, dist_m)

# Exclude the ENTIRE coordinate group if two or more retained CBBC SITEIDs
# share the same recorded latitude/longitude. The rule is established in P1,
# independently of BAM counts, species identity or the size of density ratios.
# Retain an audit of excluded sites that would otherwise have had surveys.
excluded_matched_sites <- sites |>
  filter(shared_coords) |>
  semi_join(y = distinct(eligible_pairs, SITEID), by = "SITEID") |>
  arrange(latitude, longitude, SITEID)

unambiguous_pairs <- eligible_pairs |>
  anti_join(y = select(sites, SITEID, shared_coords) |> filter(shared_coords), by = "SITEID")

# Assign each survey only once here (rather than later in README.Rmd): P3's
# species-occurrence threshold and EVERY reported comparison now use exactly
# the same final set of CBBC sites, including any future non-identical sites
# that happen to fall within the same 1-km survey matching circle.
site_survey <- unambiguous_pairs |>
  group_by(survey_id) |>
  slice_min(order_by = dist_m, n = 1, with_ties = FALSE) |>
  ungroup() |>
  arrange(SITEID, dist_m)

if (anyDuplicated(site_survey$survey_id)) stop("A BAM survey was assigned to multiple sites.")
if (any(site_survey$SITEID %in% excluded_matched_sites$SITEID)) {
  stop("A shared-coordinate CBBC site survived the exclusion.")
}
assigned_sites <- sites |> filter(SITEID %in% site_survey$SITEID)
if (anyDuplicated(select(assigned_sites, latitude, longitude))) {
  stop("Two final matched CBBC sites have identical coordinates.")
}

exclusion_audit <- tibble(
  excluded_cbbc_sites_with_surveys = nrow(excluded_matched_sites),
  excluded_coordinate_groups_with_surveys = nrow(distinct(excluded_matched_sites, latitude, longitude)),
  distinct_surveys_matching_excluded_sites = n_distinct(eligible_pairs$survey_id[eligible_pairs$SITEID %in% excluded_matched_sites$SITEID]),
  surveys_without_unambiguous_site = n_distinct(eligible_pairs$survey_id) - n_distinct(unambiguous_pairs$survey_id),
  final_cbbc_sites = n_distinct(site_survey$SITEID),
  final_surveys = n_distinct(site_survey$survey_id)
)
cat("\nShared-coordinate exclusion and final-site audit:\n")
print(exclusion_audit, width = Inf)
print(excluded_matched_sites, n = Inf)

survey_funnel <- bind_rows(
  funnel_row(df = matched,     stage = sprintf("All BAM surveys within %g km", MATCH_RADIUS_M / 1000)),
  funnel_row(df = f_sensor,    stage = "Human point counts (sensor PC; ARUs excluded)"),
  funnel_row(df = f_method,    stage = "Standardized point counts (method PC; eBird excluded)"),
  funnel_row(df = f_season,    stage = sprintf("Breeding season (day of year %d-%d)", DOY_RANGE[1], DOY_RANGE[2])),
  funnel_row(df = f_design,    stage = "Valid duration and distance limit"),
  funnel_row(df = eligible_pairs, stage = sprintf("Morning (%g to %g h from sunrise)", TSSR_RANGE[1], TSSR_RANGE[2])),
  funnel_row(df = unambiguous_pairs, stage = "Exclude all CBBC sites sharing coordinates"),
  funnel_row(df = site_survey, stage = "Assign each survey to nearest remaining CBBC site")
)
print(survey_funnel)

# ==============================================================================
# 4. Counts for the retained surveys
# ==============================================================================

dbWriteTable(conn = con, name = "kept_ids", value = distinct(site_survey, survey_id),
             temporary = TRUE, overwrite = TRUE)

# Recoding can merge two codes in one survey, so counts are re-summed per
# survey x species; a merged row is flagged if either source row was flagged
bam_counts <- tbl(src = con, from = I("bam.main.bird")) |>
  semi_join(y = tbl(src = con, from = "kept_ids"), by = "survey_id") |>
  select(survey_id, species, count, flag_err_count) |>
  collect() |>
  mutate(species = coalesce(unname(BAM_RECODE[species]), species)) |>
  group_by(survey_id, species) |>
  summarise(count = sum(count), flag_err_count = any(flag_err_count), .groups = "drop")

dbDisconnect(conn = con, shutdown = TRUE)

# ==============================================================================
# 5. Save
# ==============================================================================

out <- list(site_survey = site_survey, bam_counts = bam_counts,
            survey_funnel = survey_funnel,
            excluded_matched_sites = excluded_matched_sites,
            exclusion_audit = exclusion_audit)
attr(out, "settings") <- list(
  bam_db           = basename(BAM_DB_PATH),
  visit_filter     = VISIT_FILTER_SQL,
  match_radius_m   = MATCH_RADIUS_M,
  sensors          = SENSORS,
  methods          = METHODS,
  doy_range        = DOY_RANGE,
  tssr_range       = TSSR_RANGE,
  exclude_shared_coordinates = TRUE,
  assign_nearest_before_species_selection = TRUE,
  distance_method  = "haversine, R = 6371008.8 m",
  created          = Sys.time()
)

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
saveRDS(object = out, file = OUT_PATH)
message("Saved: ", normalizePath(OUT_PATH))
