# ==============================================================================
# 00_inspect_datasets.R
#
# Loads every dataset used in the CBBC vs point-count / model comparison and
# writes a structure report (classes, dimensions, columns, value ranges and
# join checks) to REPORT_PATH. Objects are left in the global environment so
# they can be explored with View().
#
# Each load is wrapped in safe_load(), so a missing file prints a message and
# the rest of the script continues.
# ==============================================================================

library(tidyverse)
library(readxl)
library(sf)
library(terra)
library(ebirdst)
library(Matrix)   # in case BAM bird / offset matrices are stored as sparse

# ------------------------------------------------------------------------------
# Paths
# ------------------------------------------------------------------------------

ILES_ROOT      <- "C:/Users/IlesD/OneDrive - EC-EC/Iles"
CBBC_TABLE_DIR <- file.path(ILES_ROOT, "Projects/Landbirds/Breeding-Bird-Census/data/tables")

CBBC_PLOTS_PATH   <- file.path(CBBC_TABLE_DIR, "PLOTS.xlsx")
CBBC_SPECIES_PATH <- file.path(CBBC_TABLE_DIR, "SPECIES.xlsx")
SUITABILITY_PATH  <- file.path(CBBC_TABLE_DIR, "bird_survey_suitability.xlsx")
HAMEL_PATH        <- file.path(CBBC_TABLE_DIR, "Hamel_1984_Appendix.xlsx")

BAM_NM5_PATH    <- file.path(ILES_ROOT, "Projects/Landbirds/BBMP/Data/Bird_Data_Raw/BAM/04_NM5.0_data_stratify.Rdata")
BAM_DATA_PATH   <- file.path(ILES_ROOT, "Projects/Landbirds/BAMDataset/04_BAMDataset_WT-2026-03-09_EBd-Jan-2026.Rdata")
BAM_RASTER_ROOT <- file.path(ILES_ROOT, "Data/BAMv5_model_outputs")   # parent of the year folder(s)
BCR_GPKG_PATH   <- file.path(ILES_ROOT, "Data/BCR_2026/bcr_2026g.gpkg")

# EDIT: species lookup written by 00_construct_species_list.R in the hexagon project
SPECIES_LOOKUP_PATH <- "PATH/TO/HEX_PROJECT/output/cache/species_lookup_ebird_bam.rds"

# ------------------------------------------------------------------------------
# Constants
# ------------------------------------------------------------------------------

EXAMPLE_BAM_CODE   <- "OVEN"
EXAMPLE_EBIRD_CODE <- "ovenbi1"
EXAMPLE_BAM_YEAR   <- 2020

INSPECT_NM5         <- TRUE   # set FALSE to skip the large BAM objects if memory is tight
INSPECT_BAM_DATASET <- TRUE

OUT_DIR           <- "output"
REPORT_PATH       <- file.path(OUT_DIR, "00_inspect_datasets_report.txt")
CBBC_SPECIES_CSV  <- file.path(OUT_DIR, "00_cbbc_species_codes.csv")

MAX_DESCRIBE_COLS <- 80    # wider tables: print all names, describe the first N
N_EXAMPLES        <- 3     # example values shown per column
MAX_LEVELS        <- 30    # columns with <= this many values are tabulated in full
QUANTILE_PROBS    <- c(0, 0.05, 0.25, 0.5, 0.75, 0.95, 1)

# Survey-design columns to profile in BAM tables (matched case-insensitively)
BAM_PROFILE_PATTERN <- "dur|dis|method|source|project|organi|protocol|type|aru|sensor|year|off"

# ------------------------------------------------------------------------------
# Helpers
# ------------------------------------------------------------------------------

section <- function(title) {
  cat("\n", strrep("=", 78), "\n", title, "\n", strrep("=", 78), "\n", sep = "")
}

safe_load <- function(label, expr) {
  tryCatch(expr, error = function(e) {
    cat("!! Could not load ", label, ": ", conditionMessage(e), "\n", sep = "")
    NULL
  })
}

show_all <- function(x) print(as_tibble(x), n = Inf, width = Inf)

print_names <- function(x, label = "names") {
  cat(label, " (", length(x), "):\n", sep = "")
  cat(strwrap(paste(x, collapse = " "), width = 150), sep = "\n")
}

load_rdata <- function(path) {
  env <- new.env()
  load(file = path, envir = env)
  as.list(env)
}

get_col <- function(x, col) if (is.data.frame(x)) x[[col]] else x[, col]

value_range <- function(v, f) {
  ok <- (is.numeric(v) || inherits(v, c("Date", "POSIXt"))) && any(!is.na(v))
  if (ok) format(f(v, na.rm = TRUE)) else ""
}

describe_columns <- function(x, cols = names(x)) {
  if (inherits(x, "sf")) x <- sf::st_drop_geometry(x)
  x <- as.data.frame(x)[, cols, drop = FALSE]
  tibble(
    column   = cols,
    class    = map_chr(x, ~ class(.x)[1]),
    n_na     = map_int(x, ~ sum(is.na(.x))),
    n_unique = map_int(x, n_distinct),
    min      = map_chr(x, ~ value_range(v = .x, f = min)),
    max      = map_chr(x, ~ value_range(v = .x, f = max)),
    examples = map_chr(x, ~ str_trunc(
      paste(head(unique(na.omit(.x)), N_EXAMPLES), collapse = " | "),
      width = 60
    ))
  )
}

describe_object <- function(x, name, skip_cols = NULL) {
  cat("\n---- ", name, " ----\n", sep = "")
  cat("class: ", paste(class(x), collapse = ", "), "\n", sep = "")
  cat("size:  ", format(object.size(x), units = "MB"), "\n", sep = "")

  if (is.data.frame(x)) {
    cat("dim:   ", nrow(x), " rows x ", ncol(x), " cols\n", sep = "")
    if (inherits(x, "sf")) cat("crs:   ", sf::st_crs(x)$input, "\n", sep = "")
    cols <- setdiff(names(x), c(skip_cols, attr(x, "sf_column")))
    if (length(cols) > MAX_DESCRIBE_COLS) {
      print_names(x = cols, label = "All columns")
      cat("Describing the first ", MAX_DESCRIBE_COLS, " columns only\n", sep = "")
      cols <- cols[seq_len(MAX_DESCRIBE_COLS)]
    }
    show_all(describe_columns(x = x, cols = cols))

  } else if (is.matrix(x) || inherits(x, "Matrix")) {
    cat("dim:   ", nrow(x), " x ", ncol(x), "\n", sep = "")
    cat("rownames (head): ", paste(head(rownames(x), 6), collapse = ", "), "\n", sep = "")
    print_names(x = colnames(x), label = "colnames")
    print(x[seq_len(min(5, nrow(x))), seq_len(min(8, ncol(x))), drop = FALSE])

  } else {
    str(x, max.level = 1, list.len = 20)
  }
  invisible(NULL)
}

# Report an error in a summary without stopping the script (e.g. a column
# that turns out to be character rather than numeric)
guard <- function(expr) tryCatch(expr, error = function(e) cat("!! ", conditionMessage(e), "\n", sep = ""))

tab_col <- function(df, col, label = col, transform = identity) {
  if (!col %in% names(df)) return(invisible(cat("\n[", col, " not found]\n", sep = "")))
  cat("\n", label, ":\n", sep = "")
  guard(print(table(transform(df[[col]]), useNA = "ifany")))
}

quantile_col <- function(df, col, label = col) {
  if (!col %in% names(df)) return(invisible(cat("\n[", col, " not found]\n", sep = "")))
  cat("\n", label, " quantiles:\n", sep = "")
  guard(print(quantile(df[[col]], probs = QUANTILE_PROBS, na.rm = TRUE)))
}

profile_cols <- function(df, pattern, exclude = NULL) {
  cols <- setdiff(grep(pattern = pattern, x = names(df), value = TRUE, ignore.case = TRUE), exclude)
  if (length(cols) == 0) return(invisible(cat("[no columns matching '", pattern, "']\n", sep = "")))
  for (col in cols) {
    v <- df[[col]]
    cat("\n", col, " (", class(v)[1], "):\n", sep = "")
    guard(if (n_distinct(v) <= MAX_LEVELS) {
      print(sort(table(v, useNA = "ifany"), decreasing = TRUE))
    } else if (is.numeric(v)) {
      print(quantile(v, probs = QUANTILE_PROBS, na.rm = TRUE))
    } else {
      print(head(sort(table(v, useNA = "ifany"), decreasing = TRUE), MAX_LEVELS))
    })
  }
}

# ------------------------------------------------------------------------------
# Start report
# ------------------------------------------------------------------------------

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
while (sink.number() > 0) sink()
sink(file = REPORT_PATH, split = TRUE)
options(width = 250)

cat("Report generated: ", format(Sys.time()), "\n", sep = "")
cat("R ", R.version$major, ".", R.version$minor,
    " | sf ", as.character(packageVersion("sf")),
    " | terra ", as.character(packageVersion("terra")),
    " | ebirdst ", as.character(packageVersion("ebirdst")), "\n", sep = "")

# ==============================================================================
# 1. CBBC PLOTS
# ==============================================================================

section("1. CBBC PLOTS table")

cbbc_plots <- safe_load("CBBC PLOTS", read_xlsx(path = CBBC_PLOTS_PATH))

if (!is.null(cbbc_plots)) {
  describe_object(x = cbbc_plots, name = "cbbc_plots")

  if ("SITEID" %in% names(cbbc_plots)) {
    cat("\nNumber of census records per SITEID (records : sites):\n")
    print(table(table(cbbc_plots$SITEID)))
  }

  # DD.MM.SS stored as text should have 2 dots; 1 dot or a numeric class means
  # the value may have been read as a number (trailing zeros lost)
  dot_count <- function(v) str_count(string = as.character(v), pattern = fixed("."))
  tab_col(df = cbbc_plots, col = "LAT_N",  label = "LAT_N: number of '.' separators",  transform = dot_count)
  tab_col(df = cbbc_plots, col = "LONG_W", label = "LONG_W: number of '.' separators", transform = dot_count)

  quantile_col(df = cbbc_plots, col = "PLOTSIZE", label = "PLOTSIZE (ha)")
  quantile_col(df = cbbc_plots, col = "HOURS",    label = "HOURS (person-hours)")
  tab_col(df = cbbc_plots, col = "VISITS")
  tab_col(df = cbbc_plots, col = "YEAR", label = "YEAR by decade", transform = function(v) floor(v / 10) * 10)
  tab_col(df = cbbc_plots, col = "PROVINCE")
  tab_col(df = cbbc_plots, col = "HABCODE1", label = "HABCODE1 (first 2 characters)",
          transform = function(v) str_sub(string = v, start = 1, end = 2))
}

# ==============================================================================
# 2. CBBC SPECIES
# ==============================================================================

section("2. CBBC SPECIES table")

cbbc_species <- safe_load("CBBC SPECIES", read_xlsx(path = CBBC_SPECIES_PATH))

if (!is.null(cbbc_species)) {
  describe_object(x = cbbc_species, name = "cbbc_species")

  # Density columns: flag codes documented in the metadata (-999 = missing
  # PLOTDENS, 0.1 = '+' partial territory)
  for (col in intersect(c("PLOTDENS", "DENSKM2"), names(cbbc_species))) {
    v <- cbbc_species[[col]]
    if (!is.numeric(v)) {
      cat("\n", col, " is not numeric. Non-numeric values:\n", sep = "")
      print(table(v[is.na(suppressWarnings(as.numeric(v)))], useNA = "ifany"))
      next
    }
    cat("\n", col, ": n = -999: ", sum(v == -999, na.rm = TRUE),
        " | n = 0.1: ", sum(v == 0.1, na.rm = TRUE),
        " | n <= 0: ", sum(v <= 0, na.rm = TRUE),
        " | NA: ", sum(is.na(v)), "\n", sep = "")
    print(quantile(v[v > 0], probs = QUANTILE_PROBS, na.rm = TRUE))
  }

  if (all(c("RECNO", "CODE") %in% names(cbbc_species))) {
    n_dup <- cbbc_species |> count(RECNO, CODE) |> filter(n > 1) |> nrow()
    cat("\nDuplicate RECNO x CODE combinations: ", n_dup, "\n", sep = "")
  }

  if (!is.null(cbbc_plots) && "RECNO" %in% names(cbbc_plots)) {
    cat("SPECIES RECNOs absent from PLOTS: ",
        sum(!unique(cbbc_species$RECNO) %in% cbbc_plots$RECNO), "\n", sep = "")
    cat("PLOTS RECNOs with no SPECIES rows: ",
        sum(!unique(cbbc_plots$RECNO) %in% cbbc_species$RECNO), "\n", sep = "")
  }

  # Full species list written to CSV (used for the species crosswalk)
  cbbc_species_codes <- cbbc_species |>
    count(across(any_of(c("CODE", "SPPCODE", "COMNAME", "COMMONNAME", "TAXON", "RPTNAME"))),
          name = "n_records", sort = TRUE)
  write_csv(x = cbbc_species_codes, file = CBBC_SPECIES_CSV)
  cat("\nDistinct species entries: ", nrow(cbbc_species_codes),
      " (full list written to ", CBBC_SPECIES_CSV, "). Top 30:\n", sep = "")
  print(head(cbbc_species_codes, 30), width = Inf)
}

# ==============================================================================
# 3. Suitability table and Hamel (1984) appendix
# ==============================================================================

section("3. Survey suitability table and Hamel (1984) appendix")

suitability <- safe_load("suitability table", read_xlsx(path = SUITABILITY_PATH))
if (!is.null(suitability)) describe_object(x = suitability, name = "suitability")

hamel <- safe_load("Hamel appendix", read_xlsx(path = HAMEL_PATH))
if (!is.null(hamel)) describe_object(x = hamel, name = "hamel")

# ==============================================================================
# 4. Species lookup (hexagon project) and CBBC crosswalk check
# ==============================================================================

section("4. Species lookup and CBBC code matching")

species_lookup <- safe_load("species lookup", readRDS(file = SPECIES_LOOKUP_PATH))

if (!is.null(species_lookup)) {
  describe_object(x = species_lookup, name = "species_lookup")

  if (!is.null(cbbc_species) && "CODE" %in% names(cbbc_species)) {
    lookup_codes <- unique(na.omit(c(species_lookup$bbl_code,
                                     species_lookup$AOS_CODE_4L,
                                     species_lookup$bam_species_code)))

    cbbc_unmatched <- cbbc_species_codes |> filter(!CODE %in% lookup_codes)
    if ("TAXON" %in% names(cbbc_unmatched)) {
      cbbc_unmatched$taxon_in_lookup <- cbbc_unmatched$TAXON %in% species_lookup$scientific_name
    }

    cat("\nCBBC species codes not found in lookup (bbl / AOS / BAM codes): ",
        nrow(cbbc_unmatched), "\n", sep = "")
    show_all(cbbc_unmatched)
  }
}

# ==============================================================================
# 5. BAM NM5.0 modelling dataset (visit / bird / offsets)
# ==============================================================================

section("5. BAM NM5.0 dataset")

nm5 <- if (INSPECT_NM5) safe_load("BAM NM5.0 Rdata", load_rdata(path = BAM_NM5_PATH)) else NULL

if (!is.null(nm5)) {
  cat("Objects in file:\n")
  show_all(tibble(
    object  = names(nm5),
    class   = map_chr(nm5, ~ class(.x)[1]),
    size_MB = map_dbl(nm5, ~ round(as.numeric(object.size(.x)) / 1e6, 1))
  ))

  walk(names(nm5), ~ describe_object(x = nm5[[.x]], name = .x))

  if (all(c("visit", "bird", "offsets") %in% names(nm5))) {
    visit   <- nm5$visit
    bird    <- nm5$bird
    offsets <- nm5$offsets

    cat("\n-- Alignment checks --\n")
    cat("nrow visit / bird / offsets: ", nrow(visit), " / ", nrow(bird), " / ", nrow(offsets), "\n", sep = "")
    cat("head rownames(bird):    ", paste(head(rownames(bird)), collapse = ", "), "\n", sep = "")
    cat("head rownames(offsets): ", paste(head(rownames(offsets)), collapse = ", "), "\n", sep = "")
    if (is.data.frame(offsets)) {
      cat("head first column of offsets (", names(offsets)[1], "): ",
          paste(head(offsets[[1]]), collapse = ", "), "\n", sep = "")
    }

    for (id_col in intersect(c("id", "PKEY", "pkey", "surveyid", "survey_id"), names(visit))) {
      ids <- as.character(visit[[id_col]])
      cat("visit$", id_col, " identical to rownames(bird): ", identical(ids, rownames(bird)),
          " | to rownames(offsets): ", identical(ids, rownames(offsets)), "\n", sep = "")
    }

    cat("Species columns in bird: ", ncol(bird), " | offsets: ", ncol(offsets),
        " | shared: ", length(intersect(colnames(bird), colnames(offsets))), "\n", sep = "")

    if (EXAMPLE_BAM_CODE %in% colnames(offsets)) {
      cat("\nexp(offset) for ", EXAMPLE_BAM_CODE,
          " (QPAD A*p*q; expected ~1-5 if area is in ha for unlimited-distance counts):\n", sep = "")
      print(quantile(exp(get_col(x = offsets, col = EXAMPLE_BAM_CODE)), probs = QUANTILE_PROBS, na.rm = TRUE))
    }

    cat("\n-- Survey-design columns in visit --\n")
    profile_cols(df = visit, pattern = BAM_PROFILE_PATTERN)
  }
}

# ==============================================================================
# 6. BAM survey dataset (dat; WildTrax + eBird release)
# ==============================================================================

section("6. BAM survey dataset (dat)")

bam_new <- if (INSPECT_BAM_DATASET) safe_load("BAM dataset Rdata", load_rdata(path = BAM_DATA_PATH)) else NULL

if (!is.null(bam_new)) {
  cat("Objects in file:\n")
  show_all(tibble(
    object  = names(bam_new),
    class   = map_chr(bam_new, ~ class(.x)[1]),
    size_MB = map_dbl(bam_new, ~ round(as.numeric(object.size(.x)) / 1e6, 1))
  ))

  bam_dat <- bam_new$dat

  if (!is.null(bam_dat)) {
    species_cols <- if (all(c("ABBE", "YTWA") %in% names(bam_dat))) {
      names(bam_dat)[which(names(bam_dat) == "ABBE"):which(names(bam_dat) == "YTWA")]
    } else {
      character(0)
    }

    describe_object(x = bam_dat, name = "dat (non-species columns)", skip_cols = species_cols)
    print_names(x = species_cols, label = "\nSpecies columns (ABBE:YTWA)")

    if (EXAMPLE_BAM_CODE %in% species_cols) {
      cat("\nCounts of ", EXAMPLE_BAM_CODE, " per survey:\n", sep = "")
      print(table(bam_dat[[EXAMPLE_BAM_CODE]], useNA = "ifany"))
    }

    cat("\n-- Survey-design / source columns in dat --\n")
    profile_cols(df = bam_dat, pattern = BAM_PROFILE_PATTERN, exclude = species_cols)
  }

  walk(setdiff(names(bam_new), "dat"), ~ describe_object(x = bam_new[[.x]], name = .x))
}

# ==============================================================================
# 7. BAM v5 prediction rasters
# ==============================================================================

section("7. BAM v5 prediction rasters")

cat("Sub-folders of ", BAM_RASTER_ROOT, ":\n", sep = "")
print(list.dirs(path = BAM_RASTER_ROOT, recursive = FALSE, full.names = FALSE))

cat("\nFiles for ", EXAMPLE_BAM_CODE, " (all years):\n", sep = "")
print(list.files(path = BAM_RASTER_ROOT, pattern = paste0("^", EXAMPLE_BAM_CODE, "_"), recursive = TRUE))

cat("\nNumber of .tif files in the ", EXAMPLE_BAM_YEAR, " folder: ",
    length(list.files(path = file.path(BAM_RASTER_ROOT, EXAMPLE_BAM_YEAR), pattern = "\\.tif$")), "\n", sep = "")

bam_rast <- safe_load("BAM raster", terra::rast(
  file.path(BAM_RASTER_ROOT, EXAMPLE_BAM_YEAR, paste0(EXAMPLE_BAM_CODE, "_Canada_", EXAMPLE_BAM_YEAR, ".tif"))
))

if (!is.null(bam_rast)) {
  print(bam_rast)
  cat("\nLayer summaries:\n")
  guard(print(terra::global(x = bam_rast, fun = c("min", "max", "mean"), na.rm = TRUE)))
}

# ==============================================================================
# 8. eBird Status & Trends
# ==============================================================================

section("8. eBird Status & Trends")

print(safe_load("ebirdst version", ebirdst::ebirdst_version()))

cat("\nebirdst_runs row for ", EXAMPLE_EBIRD_CODE, ":\n", sep = "")
ebirdst::ebirdst_runs |>
  filter(species_code == EXAMPLE_EBIRD_CODE) |>
  print(width = Inf)

ebird_rast <- safe_load("eBird raster", ebirdst::load_raster(
  species    = EXAMPLE_EBIRD_CODE,
  product    = "abundance",
  period     = "seasonal",
  metric     = "mean",
  resolution = "3km"
))

if (!is.null(ebird_rast)) {
  print(ebird_rast)
  cat("\nLayer summaries:\n")
  guard(print(terra::global(x = ebird_rast, fun = c("min", "max", "mean"), na.rm = TRUE)))
}

# ==============================================================================
# 9. BCR boundaries
# ==============================================================================

section("9. BCR boundaries")

bcr <- safe_load("BCR gpkg", sf::st_read(dsn = BCR_GPKG_PATH, quiet = TRUE))
if (!is.null(bcr)) describe_object(x = bcr, name = "bcr")

# ------------------------------------------------------------------------------
# Close report
# ------------------------------------------------------------------------------

sink()
message("Report written to: ", normalizePath(REPORT_PATH))
