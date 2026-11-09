# ==============================================================================
# P4_pif_corrections.R
#
# PREPARATION STEP 4 OF 5 (no inferences; run once, output cached)
#
# Computes the effective area of every retained BAM survey x species under the
# Partners in Flight (PIF) detectability corrections used for continental
# population estimates (Blancher et al. 2013; Stanton et al. 2019), so the
# report can compare them with the NA-POPS corrections from P3.
#
# PIF VS NA-POPS
#   PIF converts a BBS count to a density with two species-level constants,
#   which NA-POPS has replaced:
#     - time-of-day adjustment (replaced by NA-POPS availability p). The PIF
#       table stores it as PIF_time_of_day: the share of birds counted at the
#       mean BBS time of day relative to the species' peak time, on a 3-min BBS
#       stop, "scaled to be comparable to availability".
#     - maximum detection distance, MDD (m) (replaced by the NA-POPS effective
#       detection radius, EDR). Despite the name, PIF uses MDD as the radius of
#       the area a count samples, and its handbook (v3.1) treats empirical EDRs
#       as the benchmark MDD should match. So MDD is treated here as an
#       expert-elicited estimate of each species' EDR.
#   The PIF pair adjustment (birds counted -> individuals) is NOT used. NA-POPS
#   did not replace it, and the report compares counts with territories, not
#   individuals, for both sets of corrections. So only the two corrections that
#   differ are swapped, and any change in agreement with CBBC comes from them.
#
# PIF EFFECTIVE AREA OF ONE SURVEY (km^2)
#   The NA-POPS formula from P3 (half-normal detection), with MDD in place of
#   EDR, so PIF and NA-POPS differ only in the values of p and the radius:
#   area = p_pif * pi * MDD^2 * (1 - exp(-(r / MDD)^2)) / 1e6
#     r      survey distance limit (m; Inf for unlimited-distance counts, which
#            gives p_pif * pi * MDD^2)
#     p_pif  PIF availability. PIF_time_of_day refers to a 3-min BBS stop, but
#            BAM surveys last 3-10 min, so by default it is rescaled to each
#            survey's duration assuming a constant cue rate (the same rescaling
#            the report uses to put NA-POPS availability on a 5-min basis):
#              p_pif = 1 - (1 - PIF_time_of_day)^(duration / 3 min)
#            With PIF_SCALE_DURATION <- FALSE, PIF_time_of_day is applied as is
#            to every survey, as PIF applied it to BBS counts.
#   PIF has no date or time-since-sunrise effect, so p_pif is the same for all
#   surveys of one duration.
#
# HYBRID CORRECTIONS
#   Two hybrids swap one component at a time, so the report can show which
#   replacement (availability or detection distance) changes agreement:
#     area_p_mdd_km2    NA-POPS availability + PIF MDD (as EDR)
#     area_tod_edr_km2  PIF time of day      + NA-POPS EDR
#
# SPECIES MATCHING
#   The PIF table has no 4-letter codes, so species are matched by common name
#   (NA-POPS common name from P3). Names that differ are fixed in PIF_NAME_FIX;
#   candidate species without a match are printed.
#
# OUTPUT (OUT_PATH), a list:
#   pif_species     one row per NA-POPS candidate species: PIF time of day, MDD
#                   and pair adjustment (NA where the species is not in the
#                   PIF table), and the NA-POPS EDR used in P3, for reference
#   survey_species  survey x species (NA-POPS species with PIF values): PIF
#                   availability and the three alternative effective areas
#   attr "settings"
# ==============================================================================

library(tidyverse)
library(readxl)

# ------------------------------------------------------------------------------
# Paths
# ------------------------------------------------------------------------------

# Analysis folder inside the project (here::here() starts from the folder
# holding the .Rproj file)
ANALYSIS_DIR <- here::here()

P3_PATH        <- file.path(ANALYSIS_DIR, "output", "P3_napops_corrections.rds")
PIF_TABLE_PATH <- file.path(ANALYSIS_DIR, "data","tables", "Population_estimates_workflow_support.xlsx")
OUT_DIR        <- file.path(ANALYSIS_DIR, "output")
OUT_PATH       <- file.path(OUT_DIR, "P4_pif_corrections.rds")

# ------------------------------------------------------------------------------
# Constants
# ------------------------------------------------------------------------------

PIF_SHEET            <- 1                             # species table (sheet 2 documents the columns)
PIF_NA_STRINGS       <- c("", "NA", "#VALUE!")
PIF_REF_DURATION_MIN <- 3                             # BBS stop duration that PIF_time_of_day refers to
PIF_SCALE_DURATION   <- TRUE                          # rescale PIF availability to survey duration (see header)

# NA-POPS common name -> PIF table common name, where the two differ. Entries
# for names not among the candidates are ignored. EDIT after checking the
# unmatched species printed in section 2.
PIF_NAME_FIX <- c(`House Wren`               = "Northern House Wren",
                  `Pacific-slope Flycatcher` = "Western Flycatcher",
                  `Gray Jay`                 = "Canada Jay")

# ------------------------------------------------------------------------------
# Helper functions
# ------------------------------------------------------------------------------

# Common names are matched in lower case with single spaces
name_key <- function(x) str_to_lower(str_squish(x))

# Effective area (km^2), as in P3: half-normal detection with radius edr_m,
# truncated at the survey's distance limit. Used for NA-POPS, PIF and the
# hybrids alike, with PIF MDD as edr_m for PIF.
effective_area_km2 <- function(p, edr_m, dist_m) p * pi * edr_m^2 * (1 - exp(-(dist_m / edr_m)^2)) / 1e6

# PIF availability for a count lasting duration_min minutes (see header)
pif_availability <- function(tod, duration_min) {
  if (PIF_SCALE_DURATION) {
    1 - (1 - tod)^(duration_min / PIF_REF_DURATION_MIN)
  } else {
    tod
  }
}

# ==============================================================================
# 1. Inputs
# ==============================================================================

p3 <- readRDS(file = P3_PATH)

# PIF table: one row per common name (the table lists two non-Canadian
# scrub-jays twice; the first row is kept)
pif_table <- read_xlsx(path = PIF_TABLE_PATH, sheet = PIF_SHEET, na = PIF_NA_STRINGS) |>
  transmute(pif_name  = common_name,
            pif_tod   = as.numeric(PIF_time_of_day),
            pif_mdd_m = as.numeric(PIF_MDD),
            pif_pair  = as.numeric(PIF_pair)) |>
  filter(!is.na(pif_tod), !is.na(pif_mdd_m)) |>
  mutate(name_key = name_key(pif_name)) |>
  distinct(name_key, .keep_all = TRUE)

# ==============================================================================
# 2. Match NA-POPS candidate species to the PIF table
# ==============================================================================

pif_species <- p3$species_candidates |>
  filter(in_napops) |>
  select(species, common_name, comparable) |>
  mutate(pif_lookup = coalesce(unname(PIF_NAME_FIX[common_name]), common_name),
         name_key   = name_key(pif_lookup)) |>
  left_join(y = select(pif_table, name_key, pif_name, pif_tod, pif_mdd_m, pif_pair), by = "name_key") |>
  left_join(y = select(p3$species_corrections, species, edr_m), by = "species") |>
  mutate(in_pif = !is.na(pif_tod)) |>
  select(-pif_lookup, -name_key)

cat("NA-POPS candidate species: ", nrow(pif_species),
    " | matched to the PIF table: ", sum(pif_species$in_pif), "\n", sep = "")
cat("Unmatched species (add to PIF_NAME_FIX if the PIF table uses another name):\n")
print(filter(pif_species, !in_pif) |> select(species, common_name), n = Inf)

# ==============================================================================
# 3. Effective areas per survey x species
# ==============================================================================

# Same surveys and species as the NA-POPS effective areas in P3 (BAM duration is
# in seconds)
survey_species <- p3$survey_species |>
  filter(species %in% pif_species$species[pif_species$in_pif], !is.na(area_km2)) |>
  select(survey_id, species, duration, distance, p, edr_m) |>
  left_join(y = select(pif_species, species, pif_tod, pif_mdd_m), by = "species") |>
  mutate(p_pif            = pif_availability(tod = pif_tod, duration_min = duration / 60),
         area_pif_km2     = effective_area_km2(p = p_pif, edr_m = pif_mdd_m, dist_m = distance),
         area_p_mdd_km2   = effective_area_km2(p = p,     edr_m = pif_mdd_m, dist_m = distance),
         area_tod_edr_km2 = effective_area_km2(p = p_pif, edr_m = edr_m,     dist_m = distance)) |>
  select(survey_id, species, p_pif, area_pif_km2, area_p_mdd_km2, area_tod_edr_km2)

# Quick check: PIF effective area relative to NA-POPS, median across surveys
cat("\nPIF / NA-POPS effective area, median across surveys (> 1: PIF assumes more birds sampled, so lower densities):\n")
survey_species |>
  left_join(y = select(p3$survey_species, survey_id, species, area_km2), by = c("survey_id", "species")) |>
  group_by(species) |>
  summarise(pif_vs_napops     = median(area_pif_km2 / area_km2),
            availability_only = median(area_tod_edr_km2 / area_km2),
            distance_only     = median(area_p_mdd_km2 / area_km2),
            .groups           = "drop") |>
  arrange(pif_vs_napops) |>
  print(n = Inf)

# ==============================================================================
# 4. Save
# ==============================================================================

out <- list(pif_species = pif_species, survey_species = survey_species)
attr(out, "settings") <- list(
  pif_table            = basename(PIF_TABLE_PATH),
  pif_ref_duration_min = PIF_REF_DURATION_MIN,
  pif_scale_duration   = PIF_SCALE_DURATION,
  pif_name_fix         = PIF_NAME_FIX,
  pif_pair_used        = FALSE,
  pif_area             = "p_pif x pi x MDD^2 x (1 - exp(-(distance limit / MDD)^2)): MDD as EDR",
  created              = Sys.time()
)

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
saveRDS(object = out, file = OUT_PATH)
message("Saved: ", normalizePath(OUT_PATH))
