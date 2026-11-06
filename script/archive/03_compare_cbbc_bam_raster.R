# ==============================================================================
# 03_compare_cbbc_bam_raster.R
#
# PURPOSE
#   Compare CBBC spot-mapping densities (territories / km^2) with BAM v5 model
#   predictions (2020 density rasters) at each CBBC plot. Script 02 tested raw
#   BAM point counts corrected with NA-POPS; this script tests the BAM density
#   models themselves.
#
# WHAT THE RASTERS ARE
#   BAM v5 national models predict density of singing males (assumed males / ha;
#   converted to / km^2 with RASTER_TO_KM2) for 2020 at ~1 km resolution. The
#   models are fit to point counts with BAM's own QPAD offsets, so they embed a
#   detectability correction, but not the NA-POPS one used in script 02.
#   Predictions represent the whole pixel (mixed habitat) in 2020, whereas a
#   CBBC plot is a small, often habitat-specific area censused decades earlier.
#
# HOW PREDICTIONS ARE EXTRACTED
#   Two extractions per site x species:
#     buffer: area-weighted mean of the raster within BUFFER_RADIUS_M (1 km) of
#             the plot coordinate. PRIMARY: it has the same spatial support as
#             the 1 km point-count matching in script 02, absorbs coordinate
#             imprecision, and is not lost when the coordinate falls in a
#             water or no-data pixel.
#     pixel:  value of the pixel containing the plot coordinate (sensitivity).
#
# CBBC DATA
#   Census filtering, density definitions, species codes, the Traill's
#   Flycatcher rule and the territorial / non-territorial species classes are
#   all inherited from script 02 (via its saved output), so both comparisons
#   use identical CBBC data. Unlike script 02, a site does not need nearby BAM
#   surveys: every retained CBBC site with a raster value is used.
#
# COMPARISONS
#   - Species-level ratio of mean densities (raster / CBBC) across all sites
#     with predictions, with site-bootstrap CIs.
#   - Three-way comparison on the sites and species used in script 02:
#     raw BAM (NA-POPS) vs raster vs CBBC. raster / raw shows how closely the
#     model reproduces the corrected counts at the same places; if it is ~1,
#     both share the same disagreement with CBBC.
#   - Diagnostics of WHY the raster differs from the raw counts (section 9):
#     the size of the NA-POPS corrections, extrapolation (distance to the
#     nearest BAM detection), whether the model pulls the raw excess toward
#     CBBC, regional fit for southern deciduous species, and site-level
#     structure.
#   - Site community ratio by habitat and census attributes (as in script 02).
#
# INPUT
#   output/01_cbbc_bam_matched.rds   (site coordinates)
#   output/02_cbbc_bam_comparison.rds (CBBC data, settings, raw BAM densities)
#   BAM v5 rasters in RASTER_DIR, named <CODE>_Canada_<YEAR>.tif
#
# OUTPUT
#   OUT_RDS  extracted predictions, site-level pairs, species summaries,
#            three-way comparison; settings in attr(out, "settings")
#   OUT_CSV  species summary table
#   FIG_DIR  species ratio plot, species-mean scatter, raster vs raw scatter,
#            site ratio by habitat, raster / raw vs distance to nearest detection
# ==============================================================================

library(tidyverse)
library(terra)

# ------------------------------------------------------------------------------
# Paths
# ------------------------------------------------------------------------------

IN_01      <- "output/01_cbbc_bam_matched.rds"
IN_02      <- "output/02_cbbc_bam_comparison.rds"
RASTER_DIR <- "C:/Users/IlesD/Documents/Data_Local/BAM_v5/rasters"

OUT_DIR <- "output"
FIG_DIR <- file.path(OUT_DIR, "figures")
OUT_RDS <- file.path(OUT_DIR, "03_cbbc_bam_raster_comparison.rds")
OUT_CSV <- file.path(OUT_DIR, "03_species_summary.csv")

# ------------------------------------------------------------------------------
# Constants: rasters and extraction
# ------------------------------------------------------------------------------

RASTER_YEAR     <- 2020
RASTER_LAYER    <- 1               # EDIT if the layer names printed in section 1 show the
#   density estimate is not the first layer
RASTER_TO_KM2   <- 100             # males / ha -> males / km^2. Section 7 warns if the
#   typical ratio suggests a different unit
UNIT_WARN_FOLD  <- 20              # warn if the median species ratio is beyond 20x either way
BUFFER_RADIUS_M <- 1000            # same spatial support as script 02's match radius
PRIMARY_EXTRACT <- "buffer"        # "buffer" or "pixel" (see header)

# ------------------------------------------------------------------------------
# Constants: raster vs raw diagnostics (section 9)
# ------------------------------------------------------------------------------

N_SPLITS  <- 200                   # random half-splits of sites for the split-sample correlation

# Extrapolation diagnostic (9.2). Each raster has a layer giving, per pixel, the
# distance to the nearest detection of that species in the BAM dataset.
DD_LAYER      <- "detection_distance"
DD_UNITS      <- "km"              # label only; EDIT if the layer uses other units
DD_BREAKS     <- c(0, 30, 100, Inf)                  # distance classes for summaries
DD_LABELS     <- c("< 30", "30-100", ">= 100")
DD_NEAR       <- 50                # within-species split: near vs far sites
MIN_RAW_SITES <- 3                 # within-species split: sites with raw detections needed per class

# Species the BAM model under-predicted most relative to raw counts in the first
# run (temperate deciduous-forest species), tested for a regional pattern
SOUTHERN_DECIDUOUS <- c("WOTH", "SCTA", "EAWP", "WBNU", "BAOR", "RBGR")
SOUTH_LAT          <- 46           # ~northern edge of the southern Ontario Mixedwood Plains

# ------------------------------------------------------------------------------
# Constants: summaries and figures (shared with script 02)
# ------------------------------------------------------------------------------

N_BOOT <- 1000                     # bootstrap resamples of sites for ratio CIs
SEED   <- 1

REF_FOLDS     <- c(2, 5, 10)
REF_LINETYPES <- c("dashed", "longdash", "dotted")   # one per REF_FOLDS
REF_CAPTION   <- "Reference lines: solid = agreement; dashed = 2x; long-dash = 5x; dotted = 10x"
RATIO_BREAKS  <- sort(c(1 / 100, 1 / REF_FOLDS, 1, REF_FOLDS, 100))
CLASS_COLOURS <- c(Territorial = "#0072B2", `Non-territorial` = "#E69F00")   # Okabe-Ito

# ------------------------------------------------------------------------------
# Helper functions
# ------------------------------------------------------------------------------

# 95% percentile CI for a ratio of means (num / den), resampling sites with
# replacement. Ratio of means = ratio of sums, so column sums are used.
boot_ratio_ci <- function(num, den, n_boot) {
  idx <- matrix(data = sample.int(n = length(num), size = length(num) * n_boot, replace = TRUE), ncol = n_boot)
  r   <- colSums(matrix(data = num[idx], ncol = n_boot)) / colSums(matrix(data = den[idx], ncol = n_boot))
  quantile(x = r[is.finite(r)], probs = c(0.025, 0.975), names = FALSE)
}

# Typical ratio and share of species within 2x, 5x and 10x
agreement_summary <- function(df) {
  df |>
    filter(is.finite(ratio)) |>     # ratio = 0 (never predicted) counts as outside 10x
    summarise(
      n_species    = n(),
      median_ratio = median(ratio),
      within_2x    = mean(abs(log(ratio)) <= log(2)),
      within_5x    = mean(abs(log(ratio)) <= log(5)),
      within_10x   = mean(abs(log(ratio)) <= log(10)),
      median_rho   = median(spearman, na.rm = TRUE)
    )
}

# Species-level comparison of a BAM density column (`bam`) with CBBC across
# sites: ratio of mean densities with a site-bootstrap CI, and Spearman rank
# correlation across sites
summarise_species <- function(df, bam) {
  df |>
    mutate(bam = {{ bam }}) |>
    filter(!is.na(bam)) |>
    group_by(species, territorial) |>
    summarise(
      n_sites       = n(),
      n_sites_cbbc  = sum(cbbc_dens_km2 > 0),
      cbbc_mean_km2 = mean(cbbc_dens_km2),
      bam_mean_km2  = mean(bam),
      ratio         = bam_mean_km2 / cbbc_mean_km2,
      ci            = list(boot_ratio_ci(num = bam, den = cbbc_dens_km2, n_boot = N_BOOT)),
      spearman      = suppressWarnings(cor(x = cbbc_dens_km2, y = bam, method = "spearman")),
      year_gap      = median(year_gap),
      .groups       = "drop"
    ) |>
    mutate(ratio_lcl     = map_dbl(.x = ci, .f = 1),
           ratio_ucl     = map_dbl(.x = ci, .f = 2),
           ci            = NULL,
           species_class = if_else(condition = territorial, true = "Territorial", false = "Non-territorial")) |>
    arrange(ratio)
}

# Reference lines for agreement (solid) and 2x / 5x / 10x either way:
# vertical for ratio axes, diagonal for log-log density scatters
ref_vlines <- function() {
  c(list(geom_vline(xintercept = 1, colour = "grey30")),
    map2(.x = REF_FOLDS, .y = REF_LINETYPES,
         .f = ~ geom_vline(xintercept = c(1 / .x, .x), linetype = .y, colour = "grey55")))
}
ref_ablines <- function() {
  c(list(geom_abline(intercept = 0, slope = 1, colour = "grey30")),
    map2(.x = REF_FOLDS, .y = REF_LINETYPES,
         .f = ~ geom_abline(intercept = log10(c(1 / .x, .x)), slope = 1, linetype = .y, colour = "grey55")))
}

# Ratio-axis tick labels as fold differences, e.g. 0.2 -> "1/5"
fold_labels <- function(x) if_else(condition = x < 1, true = paste0("1/", round(1 / x)), false = as.character(x))

# ==============================================================================
# 0. Load inputs and inherit settings from script 02
# ==============================================================================

dat   <- readRDS(file = IN_01)
res02 <- readRDS(file = IN_02)
s02   <- attr(x = res02, which = "settings")
stopifnot(all(c("cbbc_long", "cbbc_site_info", "traill_records", "bad_coord_sites") %in% names(res02)))   # re-run 02 if missing

# Species rules come from script 02 so both comparisons treat species identically
EXCLUDE_SPECIES <- s02$exclude_species   # non-territorial species (flagged, summarized separately)
BAM_RECODE      <- s02$bam_recode        # BAM codes that differ from the codes used in the analysis
MIN_SITES_CBBC  <- s02$min_sites_cbbc    # species needs CBBC territories at >= this many sites
TRAILL_SPECIES  <- s02$traill_species    # unknown (not zero) in records listing Traill's Flycatcher

# ==============================================================================
# 1. Raster inventory
# ==============================================================================

# One file per BAM species; the 4-letter code is mapped to the analysis code.
# If two BAM codes map to one analysis code, their densities are summed later.
raster_files <- tibble(path = list.files(path = RASTER_DIR, full.names = TRUE,
                                         pattern = sprintf("^[A-Z0-9]{4}_Canada_%s\\.tif$", RASTER_YEAR))) |>
  mutate(bam_code = str_sub(string = basename(path), start = 1, end = 4),
         species  = coalesce(unname(BAM_RECODE[bam_code]), bam_code))

cat("BAM rasters found: ", nrow(raster_files), "\n", sep = "")
stopifnot(nrow(raster_files) > 0)

# Inspect one raster: CRS, resolution and layer names. Check that layer
# RASTER_LAYER is the density estimate before trusting the results.
r0 <- terra::rast(x = raster_files$path[1])
print(r0)
cat("Layer names: ", paste(names(r0), collapse = ", "), " | using layer ", RASTER_LAYER, "\n", sep = "")
stopifnot(RASTER_LAYER <= terra::nlyr(r0))

# ==============================================================================
# 2. CBBC sites to extract at
# ==============================================================================

# Every CBBC site retained by script 02's census filters (present in
# cbbc_long), excluding sites with suspect coordinates. Coordinates are the
# modal plot coordinates from script 01.
eligible_sites <- dat$cbbc_sites |>
  filter(SITEID %in% res02$cbbc_long$SITEID, !SITEID %in% res02$bad_coord_sites) |>
  select(SITEID, latitude, longitude)

# Points and buffers in the raster's CRS (buffer width in metres)
site_pts <- terra::vect(x = eligible_sites, geom = c("longitude", "latitude"), crs = "EPSG:4326") |>
  terra::project(y = terra::crs(r0))
site_buf <- terra::buffer(x = site_pts, width = BUFFER_RADIUS_M)

cat("Eligible CBBC sites: ", nrow(eligible_sites), "\n", sep = "")

# ==============================================================================
# 3. Species to compare
# ==============================================================================

# Species CBBC recorded as territorial at enough eligible sites, flagged by
# whether BAM has a raster for them and whether spot-mapping can census them
species_candidates <- res02$cbbc_long |>
  filter(SITEID %in% eligible_sites$SITEID, dens_km2 > 0) |>
  distinct(SITEID, species) |>
  count(species, name = "n_sites_cbbc", sort = TRUE) |>
  filter(n_sites_cbbc >= MIN_SITES_CBBC) |>
  mutate(in_raster   = species %in% raster_files$species,
         territorial = !species %in% EXCLUDE_SPECIES)

cat("\nCandidate species (CBBC territories at >= ", MIN_SITES_CBBC, " sites) without a BAM raster:\n", sep = "")
print(filter(species_candidates, !in_raster), n = Inf)

species_set <- species_candidates |> filter(in_raster) |> pull(species)
cat("\nSpecies compared: ", length(species_set),
    " (territorial: ", sum(!species_set %in% EXCLUDE_SPECIES), ")\n", sep = "")

# ==============================================================================
# 4. Extract raster predictions (density / km^2)
# ==============================================================================

# Speed: working out which cells each buffer covers (and by how much) is the
# slow step, and it is identical for every species because all rasters share
# one grid. So the cells and coverage weights are computed ONCE, only those
# cells are read from each raster, and the weighted means are taken in R.
raster_use <- raster_files |> filter(species %in% species_set)

same_grid <- map_lgl(.x = raster_use$path,
                     .f = ~ terra::compareGeom(x = terra::rast(x = .x), y = r0, stopOnError = FALSE))
if (!all(same_grid)) stop("Rasters not on the same grid as ", basename(raster_files$path[1]), ": ",
                          paste(basename(raster_use$path[!same_grid]), collapse = ", "))

# Cells overlapping each buffer, with the fraction of each cell inside it
# (ID indexes eligible_sites), and the cell containing each plot coordinate
buf_cells <- terra::cells(x = r0, y = site_buf, exact = TRUE) |>
  as_tibble() |>
  transmute(SITEID = eligible_sites$SITEID[ID], cell, weight = weights)
pix_cells <- tibble(SITEID = eligible_sites$SITEID,
                    cell   = terra::cellFromXY(object = r0, xy = terra::crds(site_pts)))

# Values of the density layer of every species raster at just those cells,
# read in one pass (one column per BAM code), in density / km^2
needed_cells <- sort(unique(na.omit(c(buf_cells$cell, pix_cells$cell))))
r_all        <- terra::rast(x = map(.x = raster_use$path, .f = ~ terra::rast(x = .x)[[RASTER_LAYER]]))
names(r_all) <- raster_use$bam_code

cell_values <- terra::extract(x = r_all, y = needed_cells) |>
  as_tibble() |>
  mutate(cell = needed_cells) |>
  pivot_longer(cols = -cell, names_to = "bam_code", values_to = "dens_km2") |>
  mutate(dens_km2 = dens_km2 * RASTER_TO_KM2)

# Buffer: area-weighted mean over covered cells, ignoring NA cells (water,
# outside the model region); all-NA buffers give NaN, set to NA.
# Pixel: value of the cell containing the plot coordinate.
buffer_values <- buf_cells |>
  inner_join(y = cell_values, by = "cell", relationship = "many-to-many") |>
  group_by(SITEID, bam_code) |>
  summarise(buffer_km2 = weighted.mean(x = dens_km2, w = weight, na.rm = TRUE), .groups = "drop")

pixel_values <- pix_cells |>
  inner_join(y = cell_values, by = "cell", relationship = "many-to-many") |>
  select(SITEID, bam_code, pixel_km2 = dens_km2)

# Site x species predictions; BAM codes mapped to the same species are summed
raster_long <- buffer_values |>
  left_join(y = pixel_values, by = c("SITEID", "bam_code")) |>
  left_join(y = distinct(raster_use, bam_code, species), by = "bam_code") |>
  group_by(SITEID, species) |>
  summarise(pixel_km2 = sum(pixel_km2), buffer_km2 = sum(buffer_km2), .groups = "drop") |>
  mutate(buffer_km2 = if_else(condition = is.nan(buffer_km2), true = NA_real_, false = buffer_km2))

# ==============================================================================
# 5. CBBC site densities for the compared species
# ==============================================================================

# As in script 02: explicit zeros for compared species absent from a census
# list (explicit = FALSE keeps listed-but-missing densities as NA), ALFL / WIFL
# unknown in Traill's records, then the mean over the site's census years
site_cbbc <- res02$cbbc_long |>
  filter(SITEID %in% eligible_sites$SITEID) |>
  complete(nesting(RECNO, SITEID, YEAR), species = species_set,
           fill = list(dens_km2 = 0), explicit = FALSE) |>
  filter(species %in% species_set) |>
  mutate(dens_km2 = if_else(condition = RECNO %in% res02$traill_records & species %in% TRAILL_SPECIES,
                            true = NA_real_, false = dens_km2)) |>
  group_by(SITEID, species) |>
  summarise(
    n_census       = sum(!is.na(dens_km2)),
    cbbc_dens_km2  = mean(dens_km2, na.rm = TRUE),
    cbbc_year_mean = mean(YEAR),
    .groups        = "drop"
  )

# ==============================================================================
# 6. Paired site-level densities
# ==============================================================================

# bam_dens_km2 is the primary extraction. Pairs without a prediction (outside
# the model region or all-NA cells) are dropped. year_gap = years between the
# CBBC censuses and the 2020 predictions. in_02 marks sites also used in the
# raw point-count comparison.
site_raster <- raster_long |>
  mutate(bam_dens_km2 = if (PRIMARY_EXTRACT == "buffer") buffer_km2 else pixel_km2) |>
  inner_join(y = site_cbbc, by = c("SITEID", "species")) |>
  filter(!is.na(bam_dens_km2), !is.nan(cbbc_dens_km2)) |>
  left_join(y = res02$cbbc_site_info, by = "SITEID") |>
  mutate(year_gap    = RASTER_YEAR - cbbc_year_mean,
         territorial = !species %in% EXCLUDE_SPECIES,
         in_02       = SITEID %in% res02$site_density$SITEID)

cat("\nSites with raster predictions: ", n_distinct(site_raster$SITEID), " of ", nrow(eligible_sites),
    " (", n_distinct(site_raster$SITEID[site_raster$in_02]), " also in script 02)\n", sep = "")

# ==============================================================================
# 7. Species summaries: ratio of mean densities (raster / CBBC) across sites
# ==============================================================================

set.seed(SEED)

species_summary       <- summarise_species(df = site_raster, bam = bam_dens_km2)   # primary
species_summary_pixel <- summarise_species(df = site_raster, bam = pixel_km2)      # sensitivity

cat("\nAgreement across species (BAM raster / CBBC):\n")
bind_rows(
  `All species`             = agreement_summary(df = species_summary),
  `Territorial only`        = agreement_summary(df = filter(species_summary, territorial)),
  `Territorial only, pixel` = agreement_summary(df = filter(species_summary_pixel, territorial)),
  .id = "subset"
) |>
  print()

# Unit check: a typical ratio far from 1 by a round factor (e.g. ~100) suggests
# the rasters are not in males / ha and RASTER_TO_KM2 needs changing
median_ratio_terr <- median(species_summary$ratio[species_summary$territorial], na.rm = TRUE)
if (abs(log(median_ratio_terr)) > log(UNIT_WARN_FOLD)) {
  warning("Median territorial ratio is ", signif(median_ratio_terr, 3),
          ": check raster units and RASTER_TO_KM2 before interpreting results.")
}

# ==============================================================================
# 8. Three-way comparison on script 02's sites: raw BAM vs raster vs CBBC
# ==============================================================================

# Same site x species pairs as the raw comparison (1 km, NA-POPS corrected).
# Differences between raster / CBBC and raw / CBBC here come from the model
# (and its offsets and 2020 prediction year), not from a different site set.
three_way <- res02$site_density |>
  select(SITEID, species, territorial, cbbc_dens_km2, raw_km2 = bam_dens_km2) |>
  inner_join(y = select(site_raster, SITEID, species, raster_km2 = bam_dens_km2),
             by = c("SITEID", "species"))

three_way_summary <- three_way |>
  group_by(species, territorial) |>
  summarise(
    n_sites         = n(),
    cbbc_mean_km2   = mean(cbbc_dens_km2),
    raw_mean_km2    = mean(raw_km2),
    raster_mean_km2 = mean(raster_km2),
    raw_vs_cbbc     = raw_mean_km2 / cbbc_mean_km2,
    raster_vs_cbbc  = raster_mean_km2 / cbbc_mean_km2,
    raster_vs_raw   = raster_mean_km2 / raw_mean_km2,
    rho_raster_raw  = suppressWarnings(cor(x = raw_km2, y = raster_km2, method = "spearman")),
    .groups         = "drop"
  )

cat("\nThree-way comparison on script 02 sites (territorial species; medians across species):\n")
three_way_summary |>
  filter(territorial, if_all(.cols = c(raw_vs_cbbc, raster_vs_cbbc, raster_vs_raw), .fns = is.finite)) |>
  summarise(   # new names, so later lines still see the per-species columns, not the medians
    n_species             = n(),
    raster_raw_within_2x  = mean(abs(log(raster_vs_raw)) <= log(2)),
    raster_raw_within_5x  = mean(abs(log(raster_vs_raw)) <= log(5)),
    med_raw_vs_cbbc       = median(raw_vs_cbbc),
    med_raster_vs_cbbc    = median(raster_vs_cbbc),
    med_raster_vs_raw     = median(raster_vs_raw),
    median_rho_raster_raw = median(rho_raster_raw, na.rm = TRUE)
  ) |>
  print()

# ==============================================================================
# 9. Diagnostics: why do raster predictions differ from raw corrected counts?
# ==============================================================================

# Section 8 compares the model with the raw counts at the same sites and
# species. Where they differ, the possible causes are:
#   (a) different detectability offsets: BAM's QPAD offsets vs NA-POPS,
#   (b) the model smoothing over local data, extrapolating far from
#       detections, or fitting some regions poorly,
#   (c) different spatial support: a 1 km area mean vs the survey locations,
#   (d) different years: 2020 predictions vs the surveys' years.
# The checks below separate these as far as the data allow.

# Per-species table for checks 9.1-9.2: three-way ratios plus the NA-POPS EDR
# and median availability used to build the raw densities (from script 02)
species_diag <- three_way_summary |>
  filter(territorial, is.finite(raster_vs_raw), raster_vs_raw > 0) |>
  left_join(y = select(res02$species_summary, species, edr_m, p_median), by = "species") |>
  arrange(raster_vs_raw)

# ---- 9.1 Does the gap scale with the size of the NA-POPS correction? --------
# If species with small EDRs or low availability (large NA-POPS corrections)
# had the lowest raster / raw ratios, NA-POPS would be inflating the raw
# densities. No relationship means the gap is not a simple function of how
# strongly NA-POPS corrects each species.
cat("\n9.1 log(raster / raw) ~ NA-POPS EDR + availability:\n")
print(summary(lm(formula = log(raster_vs_raw) ~ log(edr_m) + qlogis(p_median), data = species_diag)))

# ---- 9.2 Extrapolation: does disagreement grow with distance to the nearest
#          BAM detection of the species? ---------------------------------
# Small distances mean a prediction is anchored by nearby data on that species;
# large distances mean the model is extrapolating (e.g. a southern species at a
# northern site, or a boreal species at its southern range edge). If raster /
# raw departs from 1 mostly where distances are large, the biggest model-data
# disagreements are extrapolation, and the rasters should be judged mainly
# where they are anchored by detections.
# Confound: species rare at these sites tend both to be far from detections and
# to have noisy raw means (few detections), which on its own inflates
# |log(raster / raw)|. Two checks guard against this: the regression includes
# log(raw mean density), and the within-species comparison contrasts near and
# far sites of the SAME species, so species identity cannot explain a
# difference. (Far sites for a species are often outside its core range, so
# the within-species contrast overlaps with the regional check in 9.4.)

r_dd <- terra::rast(x = map(.x = raster_use$path, .f = ~ terra::rast(x = .x)[[DD_LAYER]]))
names(r_dd) <- raster_use$bam_code

# Area-weighted mean distance within each site's 1 km buffer (same cells and
# weights as section 4). A species with two BAM codes takes the smaller one.
site_dd <- terra::extract(x = r_dd, y = needed_cells) |>
  as_tibble() |>
  mutate(cell = needed_cells) |>
  pivot_longer(cols = -cell, names_to = "bam_code", values_to = "dd") |>
  inner_join(y = buf_cells, by = "cell", relationship = "many-to-many") |>
  group_by(SITEID, bam_code) |>
  summarise(dd = weighted.mean(x = dd, w = weight, na.rm = TRUE), .groups = "drop") |>
  left_join(y = distinct(raster_use, bam_code, species), by = "bam_code") |>
  group_by(SITEID, species) |>
  summarise(dist_detect = min(dd), .groups = "drop")

# Typical (median) distance across a set of sites, per species, with classes
species_dist <- function(sites) {
  site_dd |>
    filter(SITEID %in% sites) |>
    group_by(species) |>
    summarise(dist_detect = median(dist_detect, na.rm = TRUE), .groups = "drop") |>
    mutate(dist_class = cut(x = dist_detect, breaks = DD_BREAKS, labels = DD_LABELS, right = FALSE))
}

# Species level, on script 02's sites
species_diag <- species_diag |>
  left_join(y = species_dist(sites = three_way$SITEID), by = "species") |>
  mutate(abs_log_ratio = abs(log(raster_vs_raw)))

cat("\n9.2 Raster / raw by distance to the nearest BAM detection (", DD_UNITS,
    "; species medians across script 02 sites):\n", sep = "")
species_diag |>
  group_by(dist_class) |>
  summarise(n_species            = n(),
            median_fold_diff     = exp(median(abs_log_ratio)),   # typical size of disagreement, either direction
            median_raster_vs_raw = median(raster_vs_raw),        # typical direction
            within_2x            = mean(abs_log_ratio <= log(2)),
            .groups              = "drop") |>
  print()

cat("\n|log(raster / raw)| ~ log(1 + distance) + log(raw mean density):\n")
print(summary(lm(formula = abs_log_ratio ~ log1p(dist_detect) + log(raw_mean_km2), data = species_diag)))

# Within species: raster / raw at sites near vs far from detections of that
# species, keeping species with enough raw detections in both classes
within_species <- three_way |>
  filter(territorial) |>
  inner_join(y = site_dd, by = c("SITEID", "species")) |>
  mutate(distance = if_else(condition = dist_detect < DD_NEAR,
                            true = paste0("near (< ", DD_NEAR, ")"), false = paste0("far (>= ", DD_NEAR, ")"))) |>
  group_by(species, distance) |>
  summarise(n_sites       = n(),
            n_sites_raw   = sum(raw_km2 > 0),
            raster_vs_raw = mean(raster_km2) / mean(raw_km2),
            .groups       = "drop") |>
  filter(n_sites_raw >= MIN_RAW_SITES, is.finite(raster_vs_raw), raster_vs_raw > 0) |>
  group_by(species) |>
  filter(n() == 2) |>
  ungroup()

cat("\nWithin species, near vs far sites (", n_distinct(within_species$species), " species with >= ",
    MIN_RAW_SITES, " detecting sites in both):\n", sep = "")
within_species |>
  group_by(distance) |>
  summarise(n_species            = n(),
            median_fold_diff     = exp(median(abs(log(raster_vs_raw)))),
            median_raster_vs_raw = median(raster_vs_raw),
            .groups              = "drop") |>
  print()

# Does the raster agree better with CBBC where it is anchored by data? Species
# grouped by their median distance across ALL sites with predictions.
species_summary_dd <- species_summary |>
  left_join(y = species_dist(sites = site_raster$SITEID), by = "species")

cat("\nRaster / CBBC agreement (all sites, territorial species) by distance to the nearest BAM detection:\n")
species_summary_dd |>
  filter(territorial) |>
  group_by(dist_class) |>
  agreement_summary() |>
  print()

# ---- 9.3 Does the model pull the raw excess over CBBC toward CBBC? ---------
# If the species that raw counts overestimate most (high raw / CBBC) are the
# ones the model lowers most (low raster / raw), the model is partly correcting
# the raw excess, and the raster's closer median agreement with CBBC is
# systematic. If unrelated, that closer agreement is incidental: species-
# specific departures that happen to offset the raw excess on average.
# Caveat: both ratios contain the same raw mean (numerator of one, denominator
# of the other), so sampling noise in the raw mean alone produces a NEGATIVE
# correlation. The split-sample version removes this by computing raw / CBBC on
# a random half of the sites and raster / raw on the other half, so the two
# ratios share no raw data. Interpret only the split-sample result.
species_ratio_cor <- function(sites_a, sites_b) {
  a <- three_way |>
    filter(territorial, SITEID %in% sites_a) |>
    group_by(species) |>
    summarise(raw_vs_cbbc = mean(raw_km2) / mean(cbbc_dens_km2), .groups = "drop")
  b <- three_way |>
    filter(territorial, SITEID %in% sites_b) |>
    group_by(species) |>
    summarise(raster_vs_raw = mean(raster_km2) / mean(raw_km2), .groups = "drop")
  ab <- inner_join(x = a, y = b, by = "species") |>
    filter(is.finite(raw_vs_cbbc), is.finite(raster_vs_raw), raw_vs_cbbc > 0, raster_vs_raw > 0)
  cor(x = log(ab$raw_vs_cbbc), y = log(ab$raster_vs_raw))
}

tw_sites  <- unique(three_way$SITEID)
cor_naive <- species_ratio_cor(sites_a = tw_sites, sites_b = tw_sites)

set.seed(SEED)
cor_split <- map_dbl(.x = seq_len(N_SPLITS), .f = function(i) {
  half <- sample(x = tw_sites, size = floor(length(tw_sites) / 2))
  species_ratio_cor(sites_a = half, sites_b = setdiff(tw_sites, half))
})

cat("\n9.3 Correlation across species of log(raw / CBBC) and log(raster / raw):\n",
    "  same sites (biased negative, do not interpret): ", round(cor_naive, 2), "\n",
    "  split-sample, median of ", N_SPLITS, " splits: ", round(median(cor_split, na.rm = TRUE), 2),
    " (95% of splits: ", paste(round(quantile(cor_split, c(0.025, 0.975), na.rm = TRUE), 2), collapse = " to "), ")\n",
    sep = "")

# ---- 9.4 Is under-prediction of southern deciduous species regional? -------
# The model under-predicted temperate deciduous-forest species most (e.g.
# WOTH, SCTA, EAWP at ~1/7 to 1/4 of raw). If their gap is concentrated south
# of SOUTH_LAT (mostly the southern Ontario Mixedwood Plains), the model fits
# the temperate south less well than the boreal core, a caveat for using v5
# rasters there. A similar gap in both regions points instead to species-level
# causes (offsets, model structure). Other territorial species are shown for
# reference: if they also drop in the south, the issue is regional, not
# specific to deciduous species.
regional_species <- three_way |>
  filter(territorial) |>
  left_join(y = select(dat$cbbc_sites, SITEID, latitude), by = "SITEID") |>
  mutate(region = if_else(condition = latitude < SOUTH_LAT,
                          true = paste0("South (lat < ", SOUTH_LAT, ")"),
                          false = paste0("North (lat >= ", SOUTH_LAT, ")")),
         group  = if_else(condition = species %in% SOUTHERN_DECIDUOUS,
                          true = "Southern deciduous", false = "Other territorial")) |>
  group_by(group, species, region) |>
  summarise(n_sites       = n(),
            n_sites_raw   = sum(raw_km2 > 0),      # sites where the raw counts detected the species
            raster_vs_raw = mean(raster_km2) / mean(raw_km2),
            .groups       = "drop")

cat("\n9.4 Sites per region:\n")
three_way |>
  distinct(SITEID) |>
  left_join(y = select(dat$cbbc_sites, SITEID, latitude), by = "SITEID") |>
  count(south_of_SOUTH_LAT = latitude < SOUTH_LAT) |>
  print()

cat("\nRaster / raw by region (median across species in each group):\n")
regional_species |>
  filter(is.finite(raster_vs_raw), raster_vs_raw > 0) |>
  group_by(group, region) |>
  summarise(n_species = n(), median_raster_vs_raw = median(raster_vs_raw), .groups = "drop") |>
  print()

cat("\nSouthern deciduous species by region:\n")
regional_species |>
  filter(group == "Southern deciduous") |>
  select(species, region, n_sites_raw, raster_vs_raw) |>
  arrange(species, region) |>
  print(n = Inf)

# ---- 9.5 Site level: habitat, ranking and survey effort --------------------
# Community totals (territorial species) per site. If raster / raw varies by
# habitat, the model departs from the counts more at some plot types. The rank
# correlation shows whether the model orders sites the way the counts do. If
# disagreement shrinks with more surveys, it is mostly noise in the raw
# estimates; if not, it is structural (spatial support, smoothing, year).
site_diag <- three_way |>
  filter(territorial) |>
  left_join(y = distinct(res02$site_density, SITEID, species, n_surveys), by = c("SITEID", "species")) |>
  group_by(SITEID) |>
  summarise(raw_total    = sum(raw_km2),
            raster_total = sum(raster_km2),
            n_surveys    = max(n_surveys),
            .groups      = "drop") |>
  left_join(y = select(res02$cbbc_site_info, SITEID, hab_major), by = "SITEID") |>
  filter(raw_total > 0) |>
  mutate(raster_vs_raw = raster_total / raw_total)

cat("\n9.5 Site raster / raw (territorial community totals) by plot habitat:\n")
site_diag |>
  group_by(hab_major) |>
  summarise(n_sites = n(), median_ratio = median(raster_vs_raw),
            q25 = quantile(raster_vs_raw, 0.25), q75 = quantile(raster_vs_raw, 0.75)) |>
  print()
cat("Spearman correlation of site totals, raster vs raw: ",
    round(cor(x = site_diag$raw_total, y = site_diag$raster_total, method = "spearman"), 2), "\n", sep = "")
cat("\n|log(raster / raw)| ~ log(surveys per site):\n")
print(summary(lm(formula = abs(log(raster_vs_raw)) ~ log(n_surveys), data = site_diag)))

# ==============================================================================
# 10. Site-level diagnostics: habitat, census effort and plot size
# ==============================================================================

# Site community ratio, as in script 02: summed raster density / summed CBBC
# density over the territorial species compared
site_community <- site_raster |>
  filter(territorial) |>
  group_by(SITEID, hab_major, plot_ha, visits_mean, cbbc_year_mean) |>
  summarise(bam_total  = sum(bam_dens_km2),
            cbbc_total = sum(cbbc_dens_km2),
            .groups    = "drop") |>
  filter(cbbc_total > 0, bam_total > 0) |>
  mutate(community_ratio = bam_total / cbbc_total)

cat("\nSite community ratio (territorial species) by plot habitat:\n")
site_community |>
  group_by(hab_major) |>
  summarise(n_sites = n(), median_ratio = median(community_ratio),
            q25 = quantile(community_ratio, 0.25), q75 = quantile(community_ratio, 0.75)) |>
  print()

cat("\nlog(site community ratio) ~ habitat + census effort + plot size + census year:\n")
print(summary(lm(formula = log(community_ratio) ~ hab_major + log(visits_mean) + log(plot_ha) + cbbc_year_mean,
                 data = site_community)))

# ==============================================================================
# 11. Figures
# ==============================================================================

dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)

# (a) Species ratios (raster / CBBC) with 95% CIs, sorted
p_ratio <- species_summary |>
  filter(is.finite(ratio)) |>
  mutate(species = fct_reorder(.f = species, .x = ratio)) |>
  ggplot(mapping = aes(x = ratio, y = species, colour = species_class, shape = species_class)) +
  ref_vlines() +
  geom_linerange(mapping = aes(xmin = ratio_lcl, xmax = ratio_ucl)) +
  geom_point(size = 1.5) +
  scale_x_log10(breaks = RATIO_BREAKS, labels = fold_labels) +
  scale_colour_manual(values = CLASS_COLOURS) +
  labs(x = "BAM v5 raster / CBBC density", y = NULL, colour = NULL, shape = NULL, caption = REF_CAPTION) +
  theme_bw() +
  theme(legend.position = "top")

ggsave(filename = file.path(FIG_DIR, "03_species_ratio.png"), plot = p_ratio,
       width = 6, height = 1.5 + 0.15 * nrow(species_summary), dpi = 200, limitsize = FALSE)

# (b) Species mean densities, raster vs CBBC, log-log
p_means <- species_summary |>
  filter(cbbc_mean_km2 > 0, bam_mean_km2 > 0) |>
  ggplot(mapping = aes(x = cbbc_mean_km2, y = bam_mean_km2, colour = species_class, shape = species_class)) +
  ref_ablines() +
  geom_point() +
  geom_text(mapping = aes(label = species), size = 2.2, vjust = -0.7, check_overlap = TRUE, show.legend = FALSE) +
  scale_x_log10() +
  scale_y_log10() +
  scale_colour_manual(values = CLASS_COLOURS) +
  labs(x = "CBBC mean density (territories / km²)", y = "BAM v5 raster mean density (males / km²)",
       colour = NULL, shape = NULL, caption = REF_CAPTION) +
  theme_bw() +
  theme(legend.position = "top")

ggsave(filename = file.path(FIG_DIR, "03_species_means.png"), plot = p_means, width = 6, height = 6, dpi = 200)

# (c) Raster vs raw NA-POPS-corrected counts, species means on script 02's
#     sites. Points on the 1:1 line: the model reproduces the corrected counts.
p_raster_raw <- three_way_summary |>
  filter(raw_mean_km2 > 0, raster_mean_km2 > 0) |>
  mutate(species_class = if_else(condition = territorial, true = "Territorial", false = "Non-territorial")) |>
  ggplot(mapping = aes(x = raw_mean_km2, y = raster_mean_km2, colour = species_class, shape = species_class)) +
  ref_ablines() +
  geom_point() +
  geom_text(mapping = aes(label = species), size = 2.2, vjust = -0.7, check_overlap = TRUE, show.legend = FALSE) +
  scale_x_log10() +
  scale_y_log10() +
  scale_colour_manual(values = CLASS_COLOURS) +
  labs(x = "Raw BAM mean density, NA-POPS corrected (birds / km²)",
       y = "BAM v5 raster mean density (males / km²)",
       colour = NULL, shape = NULL, caption = REF_CAPTION) +
  theme_bw() +
  theme(legend.position = "top")

ggsave(filename = file.path(FIG_DIR, "03_raster_vs_raw.png"), plot = p_raster_raw, width = 6, height = 6, dpi = 200)

# (d) Site community ratio by plot habitat (territorial species)
p_habitat <- site_community |>
  mutate(hab_major = fct_reorder(.f = hab_major, .x = community_ratio, .fun = median)) |>
  ggplot(mapping = aes(x = community_ratio, y = hab_major)) +
  ref_vlines() +
  geom_boxplot(outlier.shape = NA, colour = "grey30") +
  geom_jitter(height = 0.15, width = 0, alpha = 0.5, colour = CLASS_COLOURS[["Territorial"]]) +
  scale_x_log10(breaks = RATIO_BREAKS, labels = fold_labels) +
  labs(x = "Site community ratio, BAM v5 raster / CBBC (territorial species)", y = NULL,
       caption = REF_CAPTION) +
  theme_bw()

ggsave(filename = file.path(FIG_DIR, "03_site_ratio_by_habitat.png"), plot = p_habitat, width = 7, height = 4, dpi = 200)

# (e) Raster / raw against distance to the nearest BAM detection (section
#     9.2), species medians on script 02 sites. If disagreement is mostly
#     extrapolation, points spread away from 1 as distance increases.
p_dist <- species_diag |>
  filter(!is.na(dist_detect), dist_detect > 0) |>
  ggplot(mapping = aes(x = raster_vs_raw, y = dist_detect)) +
  ref_vlines() +
  geom_point(colour = CLASS_COLOURS[["Territorial"]]) +
  geom_text(mapping = aes(label = species), size = 2.2, vjust = -0.7, check_overlap = TRUE) +
  scale_x_log10(breaks = RATIO_BREAKS, labels = fold_labels) +
  scale_y_log10() +
  labs(x = "BAM v5 raster / raw BAM (NA-POPS corrected), species means",
       y = paste0("Median distance to nearest BAM detection (", DD_UNITS, ")"),
       caption = REF_CAPTION) +
  theme_bw()

ggsave(filename = file.path(FIG_DIR, "03_raster_raw_vs_distance.png"), plot = p_dist, width = 6, height = 6, dpi = 200)

# ==============================================================================
# 12. Save results and the settings that produced them
# ==============================================================================

out <- list(
  raster_long           = raster_long,             # site x species predictions (pixel and buffer)
  site_raster           = site_raster,             # paired raster and CBBC site densities
  species_summary       = species_summary,         # species ratios, primary extraction
  species_summary_pixel = species_summary_pixel,   # species ratios, pixel extraction
  three_way             = three_way,               # site x species: CBBC, raw BAM, raster
  three_way_summary     = three_way_summary,       # species-level three-way ratios
  species_diag          = species_diag,            # raster / raw with NA-POPS EDR, availability, distance (9.1-9.2)
  site_dd               = site_dd,                 # site x species distance to nearest BAM detection (9.2)
  within_species        = within_species,          # raster / raw at near vs far sites, by species (9.2)
  species_summary_dd    = species_summary_dd,      # raster / CBBC ratios with distance classes (9.2)
  cor_split             = cor_split,               # split-sample correlations (9.3)
  regional_species      = regional_species,        # raster / raw by region (9.4)
  site_diag             = site_diag,               # site-level raster / raw (9.5)
  site_community        = site_community           # site-level summed densities (territorial species)
)
attr(out, "settings") <- list(
  raster_dir      = RASTER_DIR,
  raster_year     = RASTER_YEAR,
  raster_layer    = RASTER_LAYER,
  raster_to_km2   = RASTER_TO_KM2,
  buffer_radius_m = BUFFER_RADIUS_M,
  primary_extract = PRIMARY_EXTRACT,
  n_splits        = N_SPLITS,
  southern_decid  = SOUTHERN_DECIDUOUS,
  south_lat       = SOUTH_LAT,
  dd_layer        = DD_LAYER,
  dd_breaks       = DD_BREAKS,
  dd_near         = DD_NEAR,
  min_raw_sites   = MIN_RAW_SITES,
  settings_02     = s02,               # CBBC filters and species rules inherited from script 02
  n_boot          = N_BOOT,
  seed            = SEED,
  terra_version   = as.character(packageVersion("terra")),
  created         = Sys.time()
)

saveRDS(object = out, file = OUT_RDS)
write_csv(x = species_summary, file = OUT_CSV)
message("Saved: ", normalizePath(OUT_RDS), " and ", normalizePath(OUT_CSV))