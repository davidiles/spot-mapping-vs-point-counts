# ==============================================================================
# 02_compare_cbbc_bam.R
#
# PURPOSE
#   Compare breeding-bird densities from two independent sources at the same
#   places:
#     - CBBC (Canadian Breeding Bird Census): spot-mapping censuses that count
#       territories on a plot of known area. Treated here as the "gold
#       standard" (territories / km^2).
#     - BAM point counts: raw counts of birds at points near each CBBC plot,
#       converted to density using NA-POPS detectability corrections
#       (birds / km^2).
#   If the NA-POPS corrections are reasonable, the two densities should be of
#   similar magnitude for most species.
#
# HOW A POINT COUNT BECOMES A DENSITY (QPAD framework)
#   Expected count = density x (area effectively sampled), where the effective
#   area depends on two detection processes:
#     p   = availability: probability a bird present gives a cue (sings/calls)
#           during the count. Depends on count duration, day of year (od) and
#           time since sunrise (tssr). From napops::avail().
#     EDR = effective detection radius (m): the radius at which as many birds
#           are missed inside as are detected outside. From napops::edr().
#   For a count truncated at radius r:
#     effective area (km^2) = p * pi * EDR^2 * (1 - exp(-(r / EDR)^2)) / 1e6
#   For unlimited-radius counts (r = Inf) this reduces to p * pi * EDR^2.
#
# HOW SITE-LEVEL DENSITIES ARE FORMED
#   BAM:  sum(counts) / sum(effective areas) over all retained surveys near a
#         site (a ratio estimator). Each site gets one estimate, so sites with
#         hundreds of surveys don't dominate, and zero-heavy sites with only a
#         few surveys stay well behaved.
#   CBBC: mean territory density over the site's retained census years. A
#         species missing from a (complete) census list is a zero.
#   BAM counts are compared directly with CBBC territories (no pair
#   adjustment), since point counts mostly detect singing males.
#
# CBBC DATA HANDLING (from the CBBC database documentation and CWS Technical
# Report 342, Kennedy et al. 1999)
#   - Densities are computed from PLOTDENS (territories on the plot) / PLOTSIZE,
#     because DENSKM2 was sometimes rounded and does not standardize partial
#     ("+") territories.
#   - Census records whose species list is clearly shorter than TOTSPP are
#     dropped, because absences in a partial list are not true zeros.
#   - Census records must report both VISITS and PLOTSIZE, with >= 8 visits
#     and a plot of 10-100 ha (effort and edge effects; see constants), and
#     are filtered on habitat (successional and urban plots have changed since).
#   - Sites whose coordinates disagree with GRIDREF are dropped (coordinates
#     may refer to a nearby town rather than the plot).
#   - Species that spot-mapping cannot census as territories (colonial,
#     flocking, brood-parasitic, wide-ranging) are flagged and summarized
#     separately from territorial species.
#
# INPUT
#   output/01_cbbc_bam_matched.rds, from 01_match_bam_to_cbbc.R
#
# OUTPUT
#   OUT_RDS  census-record QA, site-level densities, species summaries,
#            per-survey corrections; settings in attr(out, "settings")
#   OUT_CSV  species summary table
#   FIG_DIR  species ratio plot, species-mean scatter, site-level scatter,
#            site-level ratio by habitat
#
# WHY BAM SURVEYS ARE MATCHED WITHIN 1 KM OF EACH PLOT
#   Not smaller: CBBC coordinates are often given only to the minute, or refer
#   to a nearby town rather than the plot, so tighter radii mostly add location
#   error and leave too few surveys per site (median 4 within 1 km).
#   Not larger (2 km was tested and rejected for the primary analysis):
#   - A paired comparison on the same 72 sites found NA-POPS-corrected BAM
#     density ~11% lower at 1-2 km than at 0-1 km (geometric mean 0.89, 95% CI
#     0.80-0.99). Mean survey year was identical in both bands (2010), so this
#     is a spatial effect, not population change.
#   - Species composition shifts with distance: mature-forest species (e.g.
#     WOTH, BBWA, GCKI, WIWR, woodpeckers) decline and edge / early-successional
#     species (e.g. SOSP, INBU, AMRE, CAWA, WTSP) increase. Points 1-2 km away
#     therefore sample different habitat from the plot, the mismatch the
#     comparison must avoid. For habitat specialists, species-level ratios
#     changed up to ~5-fold between bands.
#   - The 1-2 km ring has 3x the area of the inner circle (median 9.5 vs 4
#     surveys per site), so a 2 km site estimate is dominated by the less
#     representative outer surveys.
#   Community-level results changed only ~10% between radii, so the main
#   conclusions do not depend on this choice. The better way to gain precision
#   is to use a wider radius restricted to points whose land cover matches the
#   plot's habitat, rather than all points within 2 km.
# ==============================================================================

library(tidyverse)
library(napops)   # NA-POPS availability and EDR estimates
library(suncalc)  # sunrise times, for time since sunrise
library(lutz)     # time zone of each survey location

# ------------------------------------------------------------------------------
# Paths
# ------------------------------------------------------------------------------

IN_PATH <- "output/01_cbbc_bam_matched.rds"
OUT_DIR <- "output"
FIG_DIR <- file.path(OUT_DIR, "figures")
OUT_RDS <- file.path(OUT_DIR, "02_cbbc_bam_comparison.rds")
OUT_CSV <- file.path(OUT_DIR, "02_species_summary.csv")

# ------------------------------------------------------------------------------
# Constants: CBBC census-record filters
# ------------------------------------------------------------------------------

# Species-list completeness. Some records may hold only part of the species
# list (CWS's earlier Preliminary Catalogue gave only the 5 most common species
# per plot). Zero-filling a partial list would bias CBBC densities low, so a
# record is kept only if its species rows reach this fraction of TOTSPP.
LIST_COMPLETENESS_MIN <- 0.8

# Census effort. Spot-mapping needs repeated visits to separate territories;
# the international standard (Robbins 1970) is ~8-10 visits. Under-visited
# plots undercount territories, especially of quiet species. Records with
# VISITS missing are dropped: census quality cannot be confirmed, and every
# retained site then has known effort for the effort diagnostic (section 8).
MIN_VISITS <- 8

# Plot size, in ha. Records with PLOTSIZE missing are dropped: density cannot
# be computed from PLOTDENS, and plot size is needed for the diagnostics.
#   Lower limit: plots < 10 ha (the recommended minimum) are prone to edge
#     effects, with fractional territories of wide-ranging species inflating
#     densities.
#   Upper limit: Robbins (1970) recommends ~10-40 ha in forest and up to ~100 ha
#     in open habitat. Much larger "plots" (up to thousands of ha in the
#     database) cannot be spot-mapped in ~10 visits, so their densities likely
#     come from other methods or extrapolation. A plot far larger than the
#     match radius is also poorly represented by points around one coordinate,
#     and a few such plots would dominate the log(plot size) effect in section 8.
MIN_PLOT_HA <- 10
MAX_PLOT_HA <- 100

# Census year. No filter by default; raise (e.g. to 1965) as a sensitivity
# check on long gaps between CBBC and BAM surveys.
MIN_CENSUS_YEAR <- -Inf

# Habitat (Erskine codes in HABCODE1; see HABITAT.DOC). Plots in early
# successional habitat have almost certainly changed in the 30-50 years since
# the census, and urban plots are poorly represented by point counts.
#   44x  clear-cut <= 7 years after cutting
#   434  old field reverting to shrubland
#   1xy / 2xy with y = 3 or 4: young or disturbed / cut-over forest stands
#   315y with y = 3 or 4: young or disturbed swamp forest
SUCCESSIONAL_PATTERN <- "^(44|434|[12]\\d[34]|315[34])"
EXCLUDE_SUCCESSIONAL <- TRUE
EXCLUDE_URBAN        <- TRUE
HAB_MAJOR_LABELS     <- c(`1` = "Broad-leafed forest", `2` = "Conifer / mixed forest",
                          `3` = "Wetland", `4` = "Open land", `5` = "Urban")

# Coordinates. Sites whose LAT/LONG disagree with GRIDREF by more than the
# matching radius may be placed at a town rather than the plot, so their
# matched BAM surveys may not sample the plot's landscape. Sites without a
# parsable GRIDREF are kept. Optionally also drop minute-precision coordinates.
COORD_MAX_DISAGREE_M <- 1000
REQUIRE_SECONDS      <- FALSE

# "+" (part of a territory on the plot) is coded 0.1 in PLOTDENS, and also
# entered as 0.1 in DENSKM2 without converting to a per-km^2 value
PARTIAL_TERRITORY <- 0.1

# ------------------------------------------------------------------------------
# Constants: BAM survey selection
# ------------------------------------------------------------------------------

MATCH_RADIUS_M <- 1000             # use BAM surveys within this distance of a CBBC plot
                                   #   (see header for why 1 km, not 2 km)
DOY_RANGE      <- c(145, 196)      # breeding season window: 25 May - 15 Jul
TSSR_RANGE     <- c(-1, 6)         # hours since sunrise; drops evening counts that lie
                                   #   outside the NA-POPS training data
SENSORS        <- "PC"             # human point counts only: NA-POPS cue rates and EDRs
                                   #   are not calibrated for ARU recordings
METHODS        <- "PC"             # standardized point counts only (eBird checklists excluded)
MIN_SITES_CBBC <- 10               # a species is compared only if CBBC found territories
                                   #   at >= this many matched sites

# ------------------------------------------------------------------------------
# Constants: NA-POPS models
# ------------------------------------------------------------------------------

AVAIL_MODEL <- "best"              # removal model with lowest AIC (uses od and tssr)
EDR_MODEL   <- 1                   # intercept-only distance model (no road / forest covariates)

# ------------------------------------------------------------------------------
# Constants: species
# ------------------------------------------------------------------------------

# CBBC uses some old or pre-lump codes; map them to the current 4-letter codes
# used by NA-POPS and BAM. Lumped forms within one census are summed later.
CBBC_RECODE <- c(MYWA = "YRWA", AUWA = "YRWA", GRJA = "CAJA", GRAJ = "CAJA",
                 YSFL = "NOFL", RSFL = "NOFL", SCJU = "DEJU", ORJU = "DEJU",
                 YWAR = "YEWA", BTGW = "BTNW", RODO = "ROPI", YPWA = "PAWA",
                 RPHE = "RNEP", WEFL = "PSFL")

# BAM codes that differ from NA-POPS codes
BAM_RECODE  <- c(NHWR = "HOWR", WEFL = "PSFL")   # EDIT: confirm BAM's House Wren code

# Taxa split after the CBBC era are assigned from plot location
# (approximate range boundaries; only affects a few western plots)
TOWHEE_SPLIT_LON <- -101           # Rufous-sided Towhee (RSTO): west = SPTO, east = EATO
VIREO_SPLIT_LON  <- -114           # Solitary Vireo (SOVI): west of this longitude and
VIREO_SPLIT_LAT  <- 56             #   south of this latitude = CAVI, otherwise BHVI

# Traill's Flycatcher (pre-1973) is Alder + Willow Flycatcher combined, and the
# two overlap in southern Canada, so it cannot be assigned reliably. In records
# that list it, ALFL and WIFL densities are set to unknown (NA) rather than 0.
TRAILL_CODE    <- "TRFL"
TRAILL_SPECIES <- c("ALFL", "WIFL")

# Species the mapping method cannot census as territories. The CBBC report
# notes spot-mapping applies to territorial, non-colonial species; for these,
# CBBC "territories" and point-count detections measure different things. They
# are kept in the outputs but flagged (territorial = FALSE) and summarized
# separately.
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
# Constants: summaries and figures
# ------------------------------------------------------------------------------

N_BOOT <- 1000                     # bootstrap resamples of sites for ratio CIs
SEED   <- 1

N_SITE_PLOT_SPECIES <- 12          # number of species in the site-level scatter
# Fold-difference reference lines (each drawn at x and 1/x). Line types grade
# from dashed (closest to agreement) to dotted (furthest), so they can be told
# apart without colour.
REF_FOLDS     <- c(2, 5, 10)
REF_LINETYPES <- c("dashed", "longdash", "dotted")   # one per REF_FOLDS
REF_CAPTION   <- "Reference lines: solid = agreement; dashed = 2x; long-dash = 5x; dotted = 10x"
RATIO_BREAKS  <- sort(c(1 / 100, 1 / REF_FOLDS, 1, REF_FOLDS, 100))   # log-axis ticks for ratio plots
CLASS_COLOURS       <- c(Territorial = "#0072B2", `Non-territorial` = "#E69F00")   # Okabe-Ito

# ------------------------------------------------------------------------------
# Helper functions
# ------------------------------------------------------------------------------

# Decimal hours of a clock time, e.g. 06:30 -> 6.5
clock_hours <- function(t) hour(t) + minute(t) / 60 + second(t) / 3600

# Local clock-time sunrise (decimal hours) for locations in one time zone.
# suncalc returns sunrise as a time in `tz`, so DST is handled automatically.
sunrise_hours <- function(date, lat, lon, tz) {
  sr <- suncalc::getSunlightTimes(data = data.frame(date = date, lat = lat, lon = lon),
                                  keep = "sunrise", tz = tz)$sunrise
  clock_hours(sr)
}

# Effective sampled area (km^2) of one survey for one species; see header
effective_area_km2 <- function(p, edr_m, dist_m) p * pi * edr_m^2 * (1 - exp(-(dist_m / edr_m)^2)) / 1e6

# 95% percentile CI for a ratio of means (num / den), resampling sites with
# replacement. Ratio of means = ratio of sums, so column sums are used.
boot_ratio_ci <- function(num, den, n_boot) {
  idx <- matrix(data = sample.int(n = length(num), size = length(num) * n_boot, replace = TRUE), ncol = n_boot)
  r   <- colSums(matrix(data = num[idx], ncol = n_boot)) / colSums(matrix(data = den[idx], ncol = n_boot))
  quantile(x = r[is.finite(r)], probs = c(0.025, 0.975), names = FALSE)
}

# Typical ratio and share of species within 2x and 10x of CBBC
agreement_summary <- function(df) {
  df |>
    filter(is.finite(ratio)) |>     # ratio = 0 (never detected by BAM) counts as outside 10x
    summarise(
      n_species    = n(),
      median_ratio = median(ratio),
      within_2x    = mean(abs(log(ratio)) <= log(2)),
      within_5x    = mean(abs(log(ratio)) <= log(5)),
      within_10x   = mean(abs(log(ratio)) <= log(10)),
      median_rho   = median(spearman, na.rm = TRUE)
    )
}

# ==============================================================================
# 0. Load matched data from script 01
# ==============================================================================

# List with cbbc_plots, cbbc_sites, cbbc_species, site_survey, bam_counts
dat <- readRDS(file = IN_PATH)
stopifnot(all(c("survey_id", "species", "count", "flag_err_count") %in% names(dat$bam_counts)),
          all(c("gridref_dist_m", "coord_no_seconds") %in% names(dat$cbbc_sites)))   # re-run 01 if missing

# Harmonize BAM species codes with NA-POPS. Recoding can merge two codes in one
# survey, so counts are re-summed per survey x species; a merged row is flagged
# if either source row was flagged.
dat$bam_counts <- dat$bam_counts |>
  mutate(species = coalesce(unname(BAM_RECODE[species]), species)) |>
  group_by(survey_id, species) |>
  summarise(count = sum(count), flag_err_count = any(flag_err_count), .groups = "drop")

# ==============================================================================
# 1. CBBC census records: quality checks and filters
# ==============================================================================

# One row per census record (RECNO = one plot in one year) with the fields
# needed to judge whether its densities can be trusted. Each keep_* column is
# one filter, so the table shows exactly why a record was dropped.
cbbc_qa <- dat$cbbc_plots |>
  filter(!is.na(latitude)) |>
  left_join(y = count(dat$cbbc_species, RECNO, name = "n_species_rows"), by = "RECNO") |>
  mutate(
    hab_code     = as.character(HABCODE1),
    hab_major    = coalesce(unname(HAB_MAJOR_LABELS[str_sub(string = hab_code, start = 1, end = 1)]), "Unknown"),
    successional = coalesce(str_detect(string = hab_code, pattern = SUCCESSIONAL_PATTERN), FALSE),
    # Complete list: species rows reach LIST_COMPLETENESS_MIN x TOTSPP (kept if TOTSPP is missing)
    keep_list    = !is.na(n_species_rows) & (is.na(TOTSPP) | n_species_rows >= LIST_COMPLETENESS_MIN * TOTSPP),
    keep_visits  = !is.na(VISITS) & VISITS >= MIN_VISITS,
    keep_size    = !is.na(PLOTSIZE) & between(PLOTSIZE, MIN_PLOT_HA, MAX_PLOT_HA),
    keep_year    = YEAR >= MIN_CENSUS_YEAR,
    keep_habitat = !(EXCLUDE_SUCCESSIONAL & successional) & !(EXCLUDE_URBAN & hab_major == "Urban"),
    keep         = keep_list & keep_visits & keep_size & keep_year & keep_habitat
  )

cat("CBBC census records failing each filter (a record can fail several):\n")
cbbc_qa |>
  summarise(n_records = n(), across(.cols = starts_with("keep_"), .fns = ~ sum(!.x, na.rm = TRUE))) |>
  print()
cat("Records retained: ", sum(cbbc_qa$keep, na.rm = TRUE), " of ", nrow(cbbc_qa), "\n", sep = "")

cbbc_records <- cbbc_qa |>
  filter(keep) |>                  # NA (e.g. missing YEAR) is dropped
  select(RECNO, SITEID, YEAR, PLOTSIZE, VISITS, hab_code, hab_major, latitude, longitude)

# Site-level descriptors from the retained records, used later to stratify
# the comparison (habitat) and to test for effort / plot-size effects
cbbc_site_info <- cbbc_records |>
  group_by(SITEID) |>
  summarise(
    hab_major   = names(which.max(table(hab_major))),   # most common major habitat across years
    plot_ha     = median(PLOTSIZE, na.rm = TRUE),
    visits_mean = mean(VISITS, na.rm = TRUE),
    .groups     = "drop"
  )

# Sites whose coordinates are suspect (see 01: coord_no_seconds, gridref_dist_m).
# NA comparisons (no GRIDREF) evaluate to NA and are not flagged.
bad_coord_sites <- dat$cbbc_sites |>
  filter(gridref_dist_m > COORD_MAX_DISAGREE_M | (REQUIRE_SECONDS & coord_no_seconds)) |>
  pull(SITEID)
cat("Sites dropped for suspect coordinates: ", length(bad_coord_sites), "\n", sep = "")

# ==============================================================================
# 2. CBBC densities per census record (territories / km^2)
# ==============================================================================

cbbc_long <- dat$cbbc_species |>
  inner_join(y = cbbc_records, by = "RECNO") |>
  mutate(
    # Current species code: location-based assignment for split taxa,
    # lookup table for renamed / lumped codes, otherwise the CBBC code as is
    species  = case_when(
      CODE == "RSTO" ~ if_else(condition = longitude < TOWHEE_SPLIT_LON, true = "SPTO", false = "EATO"),
      CODE == "SOVI" ~ if_else(condition = longitude < VIREO_SPLIT_LON & latitude < VIREO_SPLIT_LAT,
                               true = "CAVI", false = "BHVI"),
      .default = coalesce(unname(CBBC_RECODE[CODE]), CODE)
    ),
    # Density in territories / km^2, in order of preference:
    #  1. PLOTDENS (territories on the plot) / PLOTSIZE (ha) x 100. Most precise:
    #     PLOTDENS was added to the database for this purpose, and it handles
    #     "+" correctly (0.1 territory on the plot). -999 = not given.
    #  2. DENSKM2 = 0.1 is a "+" that was never standardized: convert it as
    #     0.1 territory on the plot.
    #  3. Otherwise DENSKM2 as given (may be rounded in some sources).
    dens_km2 = case_when(
      PLOTDENS >= 0 & PLOTSIZE > 0              ~ PLOTDENS / PLOTSIZE * 100,
      near(DENSKM2, PARTIAL_TERRITORY) & PLOTSIZE > 0 ~ PARTIAL_TERRITORY / PLOTSIZE * 100,
      DENSKM2 >= 0                              ~ DENSKM2,
      .default                                  = NA_real_
    )
  ) |>
  # Sum forms that were recoded to the same species within a census
  group_by(RECNO, SITEID, YEAR, species) |>
  summarise(dens_km2 = sum(dens_km2), .groups = "drop")

# Records listing Traill's Flycatcher: ALFL / WIFL are unknown there, not zero
traill_records <- dat$cbbc_species |>
  filter(CODE == TRAILL_CODE, RECNO %in% cbbc_records$RECNO) |>
  pull(RECNO)
cat("Records with Traill's Flycatcher (ALFL / WIFL set to unknown): ", length(traill_records), "\n", sep = "")

# ==============================================================================
# 3. BAM surveys: selection and time since sunrise
# ==============================================================================

# site_survey holds every CBBC site x BAM survey pair within 10 km (script 01).
# Keep standardized human point counts close to plots that passed the CBBC
# filters, in the breeding season, with a valid duration and distance limit.
surveys <- dat$site_survey |>
  filter(SITEID %in% cbbc_records$SITEID,
         !SITEID %in% bad_coord_sites,
         dist_m <= MATCH_RADIUS_M,
         sensor %in% SENSORS,
         method %in% METHODS,
         between(doy, DOY_RANGE[1], DOY_RANGE[2]),
         duration > 0, distance > 0)

# Time since sunrise (tssr, hours) for each survey. date_time is stored as
# local clock time (a naive TIMESTAMP that R reads as UTC), so its clock hours
# are compared with sunrise in the survey's own time zone.
survey_times <- surveys |>
  distinct(survey_id, date_time, latitude, longitude) |>
  mutate(tz = lutz::tz_lookup_coords(lat = latitude, lon = longitude, method = "accurate", warn = FALSE)) |>
  filter(!is.na(tz)) |>
  group_by(tz) |>                  # suncalc takes one time zone per call
  mutate(sunrise_h = sunrise_hours(date = as.Date(date_time), lat = latitude, lon = longitude, tz = first(tz))) |>
  ungroup() |>
  mutate(tssr = clock_hours(date_time) - sunrise_h)

# Attach tssr and keep morning counts within the NA-POPS training range
surveys <- surveys |>
  inner_join(y = select(survey_times, survey_id, tssr), by = "survey_id") |>
  filter(between(tssr, TSSR_RANGE[1], TSSR_RANGE[2]))

matched_sites <- unique(surveys$SITEID)
cat("\nRetained surveys: ", n_distinct(surveys$survey_id), " at ", length(matched_sites), " CBBC sites\n", sep = "")

# ==============================================================================
# 4. Species to compare
# ==============================================================================

# Species with both an NA-POPS removal model (availability) and distance
# model (EDR); both are needed to build the effective area
napops_ok <- napops::list_species() |>
  filter(Removal == 1, Distance == 1) |>
  pull(Species)

# Candidates: species CBBC recorded as territorial at enough matched sites to
# give a meaningful comparison. Flag whether each can be corrected (NA-POPS),
# whether it appears in the BAM counts at all, and whether spot-mapping can
# census it as territories.
species_candidates <- cbbc_long |>
  filter(SITEID %in% matched_sites, dens_km2 > 0) |>
  distinct(SITEID, species) |>
  count(species, name = "n_sites_cbbc", sort = TRUE) |>
  filter(n_sites_cbbc >= MIN_SITES_CBBC) |>
  mutate(in_napops   = species %in% napops_ok,
         in_bam      = species %in% dat$bam_counts$species,
         territorial = !species %in% EXCLUDE_SPECIES)

# Candidates that will be dropped. A common species absent from BAM usually
# means a code mismatch: add it to CBBC_RECODE or BAM_RECODE.
cat("\nCandidate species (CBBC territories at >= ", MIN_SITES_CBBC, " matched sites) missing from NA-POPS or BAM:\n", sep = "")
print(filter(species_candidates, !in_napops | !in_bam), n = Inf)

species_set <- species_candidates |> filter(in_napops, in_bam) |> pull(species)
cat("\nSpecies compared: ", length(species_set),
    " (territorial: ", sum(!species_set %in% EXCLUDE_SPECIES), ")\n", sep = "")

# ==============================================================================
# 5. NA-POPS corrections per survey x species
# ==============================================================================

# Survey attributes that determine detectability
survey_design <- surveys |> distinct(survey_id, doy, tssr, duration, distance)

# For each species, availability for every survey (vectorized with
# pairwise = TRUE, so each od / tssr / duration triple is one survey) and one
# intercept-only EDR, combined into the effective area. BAM duration is in
# seconds; napops expects minutes. Warnings about covariates outside the
# training range are suppressed (tssr is already restricted above).
effort <- map(.x = species_set, .f = function(sp) {
  survey_design |>
    mutate(
      species  = sp,
      p        = as.numeric(suppressWarnings(napops::avail(
                   species = sp, model = AVAIL_MODEL, od = doy, tssr = tssr,
                   time = duration / 60, pairwise = TRUE))[[3]]),
      edr_m    = as.numeric(suppressWarnings(napops::edr(
                   species = sp, model = EDR_MODEL, road = 0, forest = 0))$EDR_est),
      area_km2 = effective_area_km2(p = p, edr_m = edr_m, dist_m = distance)
    )
}) |>
  list_rbind()

# Attach counts. The bird table only lists species detected, so a missing
# survey x species row is a zero count. Counts flagged as errors (capped
# outliers) remove that survey for that species only.
survey_species <- effort |>
  left_join(y = select(dat$bam_counts, survey_id, species, count, flag_err_count),
            by = c("survey_id", "species")) |>
  filter(!coalesce(flag_err_count, FALSE), !is.na(area_km2)) |>
  mutate(count = coalesce(count, 0))

# ==============================================================================
# 6. Site-level densities
# ==============================================================================

# BAM: ratio estimator per site x species (see header). A survey near two CBBC
# sites contributes to both, hence the many-to-many join.
site_bam <- surveys |>
  select(SITEID, survey_id, longitude, latitude, year) |>
  inner_join(y = survey_species, by = "survey_id", relationship = "many-to-many") |>
  group_by(SITEID, species) |>
  summarise(
    n_surveys     = n(),
    n_locations   = n_distinct(round(x = longitude, digits = 5), round(x = latitude, digits = 5)),
    sum_count     = sum(count),
    sum_area_km2  = sum(area_km2),
    bam_year_mean = mean(year),
    .groups       = "drop"
  ) |>
  mutate(bam_dens_km2 = sum_count / sum_area_km2)

# CBBC: add explicit zeros for compared species absent from a census list
# (explicit = FALSE keeps listed-but-missing densities as NA rather than 0),
# mark ALFL / WIFL unknown in Traill's records, then average over the site's
# census years (NA records are skipped)
site_cbbc <- cbbc_long |>
  filter(SITEID %in% matched_sites) |>
  complete(nesting(RECNO, SITEID, YEAR), species = species_set,
           fill = list(dens_km2 = 0), explicit = FALSE) |>
  filter(species %in% species_set) |>
  mutate(dens_km2 = if_else(condition = RECNO %in% traill_records & species %in% TRAILL_SPECIES,
                            true = NA_real_, false = dens_km2)) |>
  group_by(SITEID, species) |>
  summarise(
    n_census       = sum(!is.na(dens_km2)),
    cbbc_dens_km2  = mean(dens_km2, na.rm = TRUE),
    cbbc_year_mean = mean(YEAR),
    .groups        = "drop"
  )

# Paired site-level densities with site descriptors; year_gap = years between
# the BAM surveys and the CBBC censuses (population change over this gap is
# confounded with detectability error)
site_density <- site_bam |>
  inner_join(y = site_cbbc, by = c("SITEID", "species")) |>
  filter(!is.nan(cbbc_dens_km2)) |>
  left_join(y = cbbc_site_info, by = "SITEID") |>
  mutate(year_gap    = bam_year_mean - cbbc_year_mean,
         territorial = !species %in% EXCLUDE_SPECIES)

# ==============================================================================
# 7. Species summaries: ratio of mean densities (BAM / CBBC) across sites
# ==============================================================================

# Per species: mean density across matched sites from each source and their
# ratio (> 1 = BAM higher than CBBC), with a site-bootstrap CI. The Spearman
# correlation measures whether the two sources rank sites similarly. EDR and
# median availability are attached for diagnosing the ratio.
set.seed(SEED)

species_summary <- site_density |>
  group_by(species, territorial) |>
  summarise(
    n_sites       = n(),
    n_sites_cbbc  = sum(cbbc_dens_km2 > 0),
    n_sites_bam   = sum(sum_count > 0),
    cbbc_mean_km2 = mean(cbbc_dens_km2),
    bam_mean_km2  = mean(bam_dens_km2),
    ratio         = bam_mean_km2 / cbbc_mean_km2,
    ci            = list(boot_ratio_ci(num = bam_dens_km2, den = cbbc_dens_km2, n_boot = N_BOOT)),
    spearman      = suppressWarnings(cor(x = cbbc_dens_km2, y = bam_dens_km2, method = "spearman")),
    year_gap      = median(year_gap),
    .groups       = "drop"
  ) |>
  mutate(ratio_lcl = map_dbl(.x = ci, .f = 1), ratio_ucl = map_dbl(.x = ci, .f = 2), ci = NULL,
         species_class = if_else(condition = territorial, true = "Territorial", false = "Non-territorial")) |>
  left_join(y = effort |>
              group_by(species) |>
              summarise(edr_m = first(edr_m), p_median = median(p), .groups = "drop"),
            by = "species") |>
  arrange(ratio)

# Overall agreement, for all species and for territorial species only. The
# territorial subset is the fair test of the NA-POPS corrections.
cat("\nAgreement across species (BAM / CBBC):\n")
bind_rows(
  `All species`        = agreement_summary(df = species_summary),
  `Territorial only`   = agreement_summary(df = filter(species_summary, territorial)),
  .id = "subset"
) |>
  print()

# ==============================================================================
# 8. Site-level diagnostics: habitat, census effort and plot size
# ==============================================================================

# Community ratio per site: summed BAM density / summed CBBC density over the
# territorial species compared. Summing across species gives one stable number
# per site, so habitat and census-quality effects can be examined even though
# single-species site ratios are mostly noise.
site_community <- site_density |>
  filter(territorial) |>
  group_by(SITEID, hab_major, plot_ha, visits_mean, cbbc_year_mean) |>
  summarise(bam_total  = sum(bam_dens_km2),
            cbbc_total = sum(cbbc_dens_km2),
            n_surveys  = max(n_surveys),
            .groups    = "drop") |>
  filter(cbbc_total > 0, bam_total > 0) |>
  mutate(community_ratio = bam_total / cbbc_total)

# By major habitat: mismatch between plot habitat and the surrounding
# landscape sampled by point counts should show up as habitat differences
cat("\nSite community ratio (territorial species) by plot habitat:\n")
site_community |>
  group_by(hab_major) |>
  summarise(n_sites = n(), median_ratio = median(community_ratio),
            q25 = quantile(community_ratio, 0.25), q75 = quantile(community_ratio, 0.75)) |>
  print()

# Census effort and plot size: if under-visited or small plots distort CBBC
# densities, log(community ratio) will trend with visits or plot size
cat("\nlog(site community ratio) ~ habitat + census effort + plot size + census year:\n")
print(summary(lm(formula = log(community_ratio) ~ hab_major + log(visits_mean) + log(plot_ha) + cbbc_year_mean,
                 data = site_community)))

# ==============================================================================
# 9. Figures
# ==============================================================================

dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)

# Reference lines for agreement (solid) and 2x / 5x / 10x differences either
# way, shared by all ratio and density plots. Vertical lines are for plots with
# a ratio on the x-axis; diagonal lines for log-log density scatters, where a
# line at intercept log10(k) with slope 1 marks BAM = k x CBBC.
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

# (a) Species ratios with 95% CIs, sorted, against the reference lines. Colour
#     and shape distinguish territorial species from those spot-mapping cannot
#     census as territories.
p_ratio <- species_summary |>
  filter(is.finite(ratio)) |>
  mutate(species = fct_reorder(.f = species, .x = ratio)) |>
  ggplot(mapping = aes(x = ratio, y = species, colour = species_class, shape = species_class)) +
  ref_vlines() +
  geom_linerange(mapping = aes(xmin = ratio_lcl, xmax = ratio_ucl)) +
  geom_point(size = 1.5) +
  scale_x_log10(breaks = RATIO_BREAKS, labels = fold_labels) +
  scale_colour_manual(values = CLASS_COLOURS) +
  labs(x = "BAM (NA-POPS corrected) / CBBC density", y = NULL, colour = NULL, shape = NULL,
       caption = REF_CAPTION) +
  theme_bw() +
  theme(legend.position = "top")

ggsave(filename = file.path(FIG_DIR, "02_species_ratio.png"), plot = p_ratio,
       width = 6, height = 1.5 + 0.15 * nrow(species_summary), dpi = 200, limitsize = FALSE)

# (b) Species mean densities, log-log; points on the 1:1 line agree exactly
p_means <- species_summary |>
  filter(cbbc_mean_km2 > 0, bam_mean_km2 > 0) |>
  ggplot(mapping = aes(x = cbbc_mean_km2, y = bam_mean_km2, colour = species_class, shape = species_class)) +
  ref_ablines() +
  geom_point() +
  geom_text(mapping = aes(label = species), size = 2.2, vjust = -0.7, check_overlap = TRUE, show.legend = FALSE) +
  scale_x_log10() +
  scale_y_log10() +
  scale_colour_manual(values = CLASS_COLOURS) +
  labs(x = "CBBC mean density (territories / km²)", y = "BAM mean density (birds / km²)",
       colour = NULL, shape = NULL, caption = REF_CAPTION) +
  theme_bw() +
  theme(legend.position = "top")

ggsave(filename = file.path(FIG_DIR, "02_species_means.png"), plot = p_means, width = 6, height = 6, dpi = 200)

# (c) Site-level densities for the most widespread territorial species; point
#     size shows the number of BAM surveys behind each site estimate.
#     Pseudo-log axes keep zeros visible.
top_species <- species_candidates |>
  filter(species %in% species_set, territorial) |>
  slice_head(n = N_SITE_PLOT_SPECIES) |>
  pull(species)

p_sites <- site_density |>
  filter(species %in% top_species) |>
  mutate(species = factor(species, levels = top_species)) |>
  ggplot(mapping = aes(x = cbbc_dens_km2, y = bam_dens_km2)) +
  geom_abline(intercept = 0, slope = 1, colour = "grey40") +
  geom_point(mapping = aes(size = n_surveys), alpha = 0.5, colour = CLASS_COLOURS[["Territorial"]]) +
  scale_x_continuous(trans = scales::pseudo_log_trans(sigma = 1)) +
  scale_y_continuous(trans = scales::pseudo_log_trans(sigma = 1)) +
  scale_size_area(max_size = 4) +
  facet_wrap(facets = vars(species), scales = "free") +
  labs(x = "CBBC density (territories / km²)", y = "BAM density (birds / km²)", size = "Surveys") +
  theme_bw()

ggsave(filename = file.path(FIG_DIR, "02_site_scatter.png"), plot = p_sites, width = 11, height = 8, dpi = 200)

# (d) Site community ratio by plot habitat (territorial species)
p_habitat <- site_community |>
  mutate(hab_major = fct_reorder(.f = hab_major, .x = community_ratio, .fun = median)) |>
  ggplot(mapping = aes(x = community_ratio, y = hab_major)) +
  ref_vlines() +
  geom_boxplot(outlier.shape = NA, colour = "grey30") +
  geom_jitter(mapping = aes(size = n_surveys), height = 0.15, width = 0, alpha = 0.5,
              colour = CLASS_COLOURS[["Territorial"]]) +
  scale_x_log10(breaks = RATIO_BREAKS, labels = fold_labels) +
  scale_size_area(max_size = 3) +
  labs(x = "Site community ratio, BAM / CBBC (territorial species)", y = NULL, size = "Surveys",
       caption = REF_CAPTION) +
  theme_bw()

ggsave(filename = file.path(FIG_DIR, "02_site_ratio_by_habitat.png"), plot = p_habitat, width = 7, height = 4, dpi = 200)

# ==============================================================================
# 10. Save results and the settings that produced them
# ==============================================================================

out <- list(
  cbbc_qa         = cbbc_qa,           # census records with each filter decision
  site_density    = site_density,      # site x species paired densities
  site_community  = site_community,    # site-level summed densities (territorial species)
  species_summary = species_summary,   # species-level ratios and diagnostics
  survey_species  = survey_species,    # per survey x species: p, EDR, effective area, count
  surveys         = surveys,           # retained site x survey pairs, with tssr
  # CBBC inputs reused by 03_compare_cbbc_bam_raster.R, so that comparison uses
  # identical census filtering, density definitions and species codes
  cbbc_long       = cbbc_long,         # retained census records x species densities
  cbbc_site_info  = cbbc_site_info,    # site habitat, plot size, visits (all retained sites)
  traill_records  = traill_records,    # records where ALFL / WIFL are unknown
  bad_coord_sites = bad_coord_sites    # sites dropped for suspect coordinates
)
attr(out, "settings") <- list(
  list_completeness_min = LIST_COMPLETENESS_MIN,
  min_visits            = MIN_VISITS,
  min_plot_ha           = MIN_PLOT_HA,
  max_plot_ha           = MAX_PLOT_HA,
  min_census_year       = MIN_CENSUS_YEAR,
  successional_pattern  = SUCCESSIONAL_PATTERN,
  exclude_successional  = EXCLUDE_SUCCESSIONAL,
  exclude_urban         = EXCLUDE_URBAN,
  coord_max_disagree_m  = COORD_MAX_DISAGREE_M,
  require_seconds       = REQUIRE_SECONDS,
  cbbc_density          = "PLOTDENS / PLOTSIZE x 100; DENSKM2 = 0.1 as '+'; else DENSKM2",
  match_radius_m        = MATCH_RADIUS_M,
  doy_range             = DOY_RANGE,
  tssr_range            = TSSR_RANGE,
  sensors               = SENSORS,
  methods               = METHODS,
  min_sites_cbbc        = MIN_SITES_CBBC,
  avail_model           = AVAIL_MODEL,
  edr_model             = EDR_MODEL,
  cbbc_recode           = CBBC_RECODE,
  bam_recode            = BAM_RECODE,
  split_taxa            = c(towhee_lon = TOWHEE_SPLIT_LON, vireo_lon = VIREO_SPLIT_LON, vireo_lat = VIREO_SPLIT_LAT),
  traill                = paste(TRAILL_CODE, "->", paste(TRAILL_SPECIES, collapse = " / "), "unknown"),
  traill_species        = TRAILL_SPECIES,
  exclude_species       = EXCLUDE_SPECIES,
  cbbc_zeros            = "species absent from a complete census list = 0",
  site_estimator        = "BAM: sum(count) / sum(p * A * q); CBBC: mean over retained census years",
  n_boot                = N_BOOT,
  seed                  = SEED,
  napops_version        = as.character(packageVersion("napops")),
  created               = Sys.time()
)

saveRDS(object = out, file = OUT_RDS)
write_csv(x = species_summary, file = OUT_CSV)
message("Saved: ", normalizePath(OUT_RDS), " and ", normalizePath(OUT_CSV))
