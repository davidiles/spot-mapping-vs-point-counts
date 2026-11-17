# ==============================================================================
# P3_napops_corrections.R
#
# PREPARATION STEP 3 OF 5 (no inferences; run once, output cached)
#
# Defines the species compared at matched sites and computes NA-POPS
# detectability corrections for every retained BAM survey x species.
#
# HOW A POINT COUNT BECOMES A DENSITY (QPAD framework)
#   Expected count = density x effective area sampled, where
#     p   = availability: probability a bird present gives a cue during the
#           count; depends on duration, day of year (od) and time since sunrise
#           (tssr). From napops::avail().
#     EDR = effective detection radius (m). From napops::edr().
#   For a count truncated at radius r:
#     effective area (km^2) = p * pi * EDR^2 * (1 - exp(-(r / EDR)^2)) / 1e6
#   which reduces to p * pi * EDR^2 for unlimited-radius counts.
#
# EDR
#   Intercept-only NA-POPS distance model: one EDR per species, averaged over
#   the habitats and road types in the NA-POPS training data. Habitat-specific
#   EDRs (survey forest cover from a land-cover raster) were tried and changed
#   densities by a median of ~3%, so they are not used.
#
# OUTPUT (OUT_PATH), a list:
#   species_candidates   species with CBBC territories at >= MIN_SITES_CBBC
#                        matched sites (whether or not BAM detected them), with
#                        NA-POPS family and flags in_napops / detected_nearby /
#                        comparable (FALSE for species in P1's
#                        non_comparable_species list)
#   survey_species       retained survey x candidate species (in BAM): count
#                        (0 if not detected) and, for NA-POPS species, p, EDR
#                        and effective area
#   species_corrections  one row per NA-POPS species: EDR (with 95% CI),
#                        availability and effective-area spread across surveys,
#                        NA-POPS model support and training-range coverage
#   attr "settings"      including the NA-POPS results release date
# ==============================================================================

library(tidyverse)
library(napops)

# ------------------------------------------------------------------------------
# Paths
# ------------------------------------------------------------------------------

# Analysis folder inside the project (here::here() starts from the folder
# holding Spot_Mapping_vs_Point_Counts.Rproj)
ANALYSIS_DIR <- here::here()

P1_PATH  <- file.path(ANALYSIS_DIR, "output", "P1_cbbc_prepared.rds")
P2_PATH  <- file.path(ANALYSIS_DIR, "output", "P2_bam_matched.rds")
OUT_DIR  <- file.path(ANALYSIS_DIR, "output")
OUT_PATH <- file.path(OUT_DIR, "P3_napops_corrections.rds")

# ------------------------------------------------------------------------------
# Constants
# ------------------------------------------------------------------------------

MIN_SITES_CBBC <- 10                # species needs territories at >= 5 FINAL retained, distinct CBBC locations

AVAIL_MODEL    <- "best"           # removal model with lowest AIC (uses od and tssr)

# Distance (EDR) model
EDR_MODEL     <- 1                 # intercept-only (see header)
EDR_QUANTILES <- c(0.025, 0.975)
EDR_SAMPLES   <- 1000              # bootstrap draws for the EDR CI

# ------------------------------------------------------------------------------
# Helper functions
# ------------------------------------------------------------------------------

# Effective sampled area (km^2) of one survey; see header
effective_area_km2 <- function(p, edr_m, dist_m) p * pi * edr_m^2 * (1 - exp(-(dist_m / edr_m)^2)) / 1e6

# Run a napops call, returning NA (with a message) if it fails for a species,
# e.g. a model not fitted or a covariance matrix that cannot be bootstrapped
try_napops <- function(expr, sp, what) {
  tryCatch(expr = suppressWarnings(expr),
           error = function(e) { message(sp, " ", what, ": ", conditionMessage(e)); NA_real_ })
}

# ==============================================================================
# 1. Inputs
# ==============================================================================

p1 <- readRDS(file = P1_PATH)
p2 <- readRDS(file = P2_PATH)

# The comparable flag (section 2) comes from P1's list of non-comparable
# species. Stop if P1 was saved before that list existed: otherwise every
# species would silently be flagged comparable.
non_comparable_species <- attr(x = p1, which = "settings")$non_comparable_species
if (is.null(non_comparable_species)) {
  stop("P1 settings have no non_comparable_species: re-run P1_prepare_cbbc.R first.")
}

# P2 has excluded whole shared-coordinate groups and assigned each BAM survey
# once to its nearest remaining CBBC site. Selection is therefore based on
# the exact retained spatial units used by all subsequent report calculations.
if (!isTRUE(attr(p2, "settings")$exclude_shared_coordinates) ||
    !isTRUE(attr(p2, "settings")$assign_nearest_before_species_selection)) {
  stop("P2 cache predates coordinate exclusion/nearest assignment: re-run P2.")
}
matched_sites <- unique(p2$site_survey$SITEID)
if (any(matched_sites %in% p1$cbbc_sites$SITEID[p1$cbbc_sites$shared_coords])) {
  stop("P3 encountered an ambiguous CBBC site; re-run P1 and P2.")
}
cat("Final unambiguous CBBC sites used for species selection:", length(matched_sites), "\n")

# NA-POPS results release in use (the package version does not identify it)
napops_data_date <- readLines(con = file.path(rappdirs::app_dir(appname = "napops")$data(), "date.txt"), warn = FALSE)[1]
cat("NA-POPS results release:", napops_data_date, "\n")

# Survey attributes that determine detectability (BAM duration is in seconds)
survey_design <- p2$site_survey |>
  distinct(survey_id, doy, tssr, duration, distance)

# ==============================================================================
# 2. Species compared at matched sites
# ==============================================================================

napops_table <- napops::list_species()
napops_ok    <- napops_table |>
  filter(Removal == 1, Distance == 1) |>
  pull(Species)

# Candidates are defined from CBBC alone (territories at >= MIN_SITES_CBBC
# matched sites), not from whether BAM detected the species: a species never
# detected near the plots is kept, with zero counts, so it cannot drop out of
# the comparison because of its result. detected_nearby records whether any
# matched survey detected it. Family (from NA-POPS) is used by the report to
# define songbirds.
# The occurrence threshold counts distinct CBBC SITEIDs AFTER shared-coordinate
# exclusions and nearest-survey assignment in P2. These are the final sampled
# locations (26 in the reference run), not the original 41 raw matched plots.
# This selects SPECIES only. Every selected species is then compared over all
# final CBBC grids, including grids where spot mapping found no territories.
# (Traill's Flycatcher is a pre-split combination, not a species to compare)
species_candidates <- p1$cbbc_long |>
  filter(SITEID %in% matched_sites, dens_km2 > 0, species != attr(x = p1, which = "settings")$traill_code) |>
  distinct(SITEID, species) |>
  count(species, name = "n_sites_cbbc", sort = TRUE) |>
  filter(n_sites_cbbc >= MIN_SITES_CBBC) |>
  left_join(y = select(napops_table, species = Species, common_name = Common_Name, family = Family), by = "species") |>
  mutate(in_napops       = species %in% napops_ok,
         detected_nearby = species %in% p2$bam_counts$species,
         comparable      = !species %in% non_comparable_species)

cat("Candidate species: ", nrow(species_candidates),
    " | with NA-POPS models: ", sum(species_candidates$in_napops),
    " | detected by a matched BAM survey: ", sum(species_candidates$detected_nearby), "\n", sep = "")
cat("Candidates never detected by a matched BAM survey. Check each code exists in BAM:",
    "a code missing from BAM entirely is a code mismatch, not a true zero.\n")
print(filter(species_candidates, !detected_nearby), n = Inf)

bam_species    <- species_candidates$species
napops_species <- species_candidates |> filter(in_napops) |> pull(species)

# ==============================================================================
# 3. NA-POPS corrections per survey x species
# ==============================================================================

# For each species: availability for every survey (pairwise = TRUE, so each
# od / tssr / duration triple is one survey), the EDR with its 95% CI, and model
# support. Warnings about covariates outside the training range are suppressed
# here and recorded instead as the share of surveys outside it.
corrections <- map(.x = napops_species, .progress = TRUE, .f = function(sp) {

  # --- Availability ---
  p <- napops::avail(species = sp, model = AVAIL_MODEL, od = survey_design$doy, tssr = survey_design$tssr,
                     time = survey_design$duration / 60, pairwise = TRUE) |>
    suppressWarnings()

  # --- EDR and its 95% CI (sampling uncertainty within this NA-POPS release only) ---
  edr <- try_napops(expr = napops::edr(species = sp, model = EDR_MODEL, road = 0, forest = 0,
                                       quantiles = EDR_QUANTILES, samples = EDR_SAMPLES),
                    sp = sp, what = "EDR")

  # --- Training ranges ---
  rem_train <- napops::covariates_removal(species = sp)

  list(
    survey = tibble(survey_id = survey_design$survey_id, species = sp,
                    p = as.numeric(p[[3]])),     # columns: TSSR, OD, p
    species = tibble(
      species        = sp,
      edr_m          = if (is.data.frame(edr)) as.numeric(edr$EDR_est) else NA_real_,
      edr_lcl        = if (is.data.frame(edr)) edr[[paste0("EDR_", EDR_QUANTILES[1] * 100)]] else NA_real_,
      edr_ucl        = if (is.data.frame(edr)) edr[[paste0("EDR_", EDR_QUANTILES[2] * 100)]] else NA_real_,
      removal_model  = napops::coef_removal(species = sp, model = AVAIL_MODEL)$Model[1],
      removal_n      = napops::coef_removal(species = sp, model = AVAIL_MODEL)$N[1],
      distance_n     = napops::coef_distance(species = sp, model = EDR_MODEL)$N[1],
      share_od_out   = mean(!between(survey_design$doy, min(rem_train$OD), max(rem_train$OD))),
      share_tssr_out = mean(!between(survey_design$tssr, min(rem_train$TSSR), max(rem_train$TSSR)))
    )
  )
})

species_corr <- map(.x = corrections, .f = "species") |> list_rbind()
survey_p     <- map(.x = corrections, .f = "survey") |> list_rbind()

# ==============================================================================
# 4. Survey x species table: counts and effective areas
# ==============================================================================

# Every retained survey x candidate species found in BAM. The bird table lists
# detections only, so a missing row is a zero count. Counts flagged as errors
# (capped outliers) remove that survey for that species only. Species without
# NA-POPS models keep their counts (used for uncorrected comparisons) with NA
# corrections.
survey_species <- expand_grid(survey_id = survey_design$survey_id, species = bam_species) |>
  left_join(y = survey_design, by = "survey_id") |>
  left_join(y = select(p2$bam_counts, survey_id, species, count, flag_err_count), by = c("survey_id", "species")) |>
  filter(!coalesce(flag_err_count, FALSE)) |>
  mutate(count = coalesce(count, 0)) |>
  select(-flag_err_count) |>
  left_join(y = survey_p, by = c("survey_id", "species")) |>
  left_join(y = select(species_corr, species, edr_m), by = "species") |>
  mutate(area_km2 = effective_area_km2(p = p, edr_m = edr_m, dist_m = distance))

# Species summary of the corrections as applied to the retained surveys
applied <- survey_species |>
  filter(species %in% napops_species) |>
  group_by(species) |>
  summarise(p_median        = median(p),
            p_q10           = quantile(p, 0.10),
            p_q90           = quantile(p, 0.90),
            area_median_km2 = median(area_km2, na.rm = TRUE),
            area_q10_km2    = quantile(area_km2, 0.10, na.rm = TRUE),
            area_q90_km2    = quantile(area_km2, 0.90, na.rm = TRUE),
            .groups         = "drop")

species_corrections <- species_corr |>
  left_join(y = applied, by = "species") |>
  left_join(y = select(species_candidates, species, n_sites_cbbc, comparable), by = "species")

# ==============================================================================
# 5. Save
# ==============================================================================

out <- list(species_candidates = species_candidates, survey_species = survey_species,
            species_corrections = species_corrections)
attr(out, "settings") <- list(
  min_sites_cbbc        = MIN_SITES_CBBC,
  avail_model           = AVAIL_MODEL,
  edr_model             = EDR_MODEL,
  edr_quantiles         = EDR_QUANTILES,
  edr_samples           = EDR_SAMPLES,
  napops_version        = as.character(packageVersion("napops")),
  napops_data_date      = napops_data_date,
  created               = Sys.time()
)

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
saveRDS(object = out, file = OUT_PATH)
message("Saved: ", normalizePath(OUT_PATH))
