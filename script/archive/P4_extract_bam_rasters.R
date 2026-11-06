# ==============================================================================
# P4_extract_bam_rasters.R
#
# PREPARATION STEP 4 OF 4 (no inferences; run once, output cached)
#
# Extracts BAM v5 2020 model predictions at every retained CBBC site (from P1),
# whether or not the site has nearby BAM surveys.
#
# WHAT THE RASTERS ARE
#   One file per species (<CODE>_Canada_<YEAR>.tif, ~1 km cells, EPSG:3978)
#   with three layers:
#     mean                 predicted density of singing males (males / ha;
#                          converted here to / km^2)
#     standard_deviation   not used
#     detection_distance   distance from the cell to the nearest detection of
#                          the species in the BAM dataset (km): small = the
#                          prediction is anchored by nearby data, large = the
#                          model is extrapolating
#   The models embed BAM's own QPAD offsets, not the NA-POPS corrections.
#
# EXTRACTION
#   buffer: area-weighted mean within BUFFER_RADIUS_M (1 km) of the plot
#           coordinate. Primary: same spatial support as the 1 km point-count
#           matching, absorbs coordinate error, and survives coordinates that
#           fall in water / no-data cells.
#   pixel:  value of the cell containing the plot coordinate (sensitivity).
#   Speed: the cells covered by each buffer (and their coverage fractions) are
#   identical for every species because all rasters share one grid, so they
#   are computed once and only those cells are read from each raster.
#
# OUTPUT (OUT_PATH), a list:
#   raster_long   site x species: pixel_km2, buffer_km2, dist_detect (buffer
#                 mean distance to nearest detection, km)
#   raster_info   CRS, resolution and layer names of the rasters
#   attr "settings"
# ==============================================================================

library(tidyverse)
library(terra)

# ------------------------------------------------------------------------------
# Paths
# ------------------------------------------------------------------------------

# Analysis folder inside the project (here::here() starts from the folder
# holding Spot_Mapping_vs_Point_Counts.Rproj)
ANALYSIS_DIR <- here::here("CBBC_vs_BAM")

P1_PATH    <- file.path(ANALYSIS_DIR, "output", "P1_cbbc_prepared.rds")
RASTER_DIR <- "C:/Users/IlesD/Documents/Data_Local/BAM_v5/rasters"
OUT_DIR    <- file.path(ANALYSIS_DIR, "output")
OUT_PATH   <- file.path(OUT_DIR, "P4_bam_raster_extract.rds")

# ------------------------------------------------------------------------------
# Constants
# ------------------------------------------------------------------------------

RASTER_YEAR     <- 2020
MEAN_LAYER      <- "mean"
DD_LAYER        <- "detection_distance"
RASTER_TO_KM2   <- 100             # males / ha -> males / km^2
BUFFER_RADIUS_M <- 1000            # same spatial support as the 1 km survey matching

# ==============================================================================
# 1. Inputs and raster inventory
# ==============================================================================

p1         <- readRDS(file = P1_PATH)
BAM_RECODE <- attr(x = p1, which = "settings")$bam_recode

# One file per BAM species; codes mapped to analysis codes (if two BAM codes map
# to one species, their densities are summed below)
raster_files <- tibble(path = list.files(path = RASTER_DIR, full.names = TRUE,
                                         pattern = sprintf("^[A-Z0-9]{4}_Canada_%s\\.tif$", RASTER_YEAR))) |>
  mutate(bam_code = str_sub(string = basename(path), start = 1, end = 4),
         species  = coalesce(unname(BAM_RECODE[bam_code]), bam_code))
stopifnot(nrow(raster_files) > 0)

r0 <- terra::rast(x = raster_files$path[1])
print(r0)
stopifnot(all(c(MEAN_LAYER, DD_LAYER) %in% names(r0)))

same_grid <- map_lgl(.x = raster_files$path, .f = ~ terra::compareGeom(x = terra::rast(x = .x), y = r0, stopOnError = FALSE))
if (!all(same_grid)) stop("Rasters not on the same grid as ", basename(raster_files$path[1]), ": ",
                          paste(basename(raster_files$path[!same_grid]), collapse = ", "))

# ==============================================================================
# 2. Sites, buffers and the cells they cover
# ==============================================================================

sites <- p1$cbbc_sites |> select(SITEID, latitude, longitude)

site_pts <- terra::vect(x = as.data.frame(sites), geom = c("longitude", "latitude"), crs = "EPSG:4326") |>
  terra::project(y = terra::crs(r0))

# Projecting a SpatVector can silently drop points that fail to transform,
# which would shift every later ID onto the wrong SITEID
stopifnot(nrow(site_pts) == nrow(sites))
site_buf <- terra::buffer(x = site_pts, width = BUFFER_RADIUS_M)

# Cells overlapping each buffer with the fraction of each cell inside it (ID
# indexes sites), and the cell containing each plot coordinate
buf_cells <- terra::cells(x = r0, y = site_buf, exact = TRUE) |>
  as_tibble() |>
  transmute(SITEID = sites$SITEID[ID], cell, weight = weights)
pix_cells <- terra::cells(x = r0, y = site_pts) |>        # matched by ID, so a point that fails
  as_tibble() |>                                          #   to project cannot misalign the sites
  transmute(SITEID = sites$SITEID[ID], cell)

needed_cells <- sort(unique(na.omit(c(buf_cells$cell, pix_cells$cell))))

# ==============================================================================
# 3. Read the needed cells from every raster (one pass per layer type)
# ==============================================================================

read_layer <- function(layer) {
  r <- terra::rast(x = map(.x = raster_files$path, .f = ~ terra::rast(x = .x)[[layer]]))
  names(r) <- raster_files$bam_code
  terra::extract(x = r, y = needed_cells) |>
    as_tibble() |>
    mutate(cell = needed_cells) |>
    pivot_longer(cols = -cell, names_to = "bam_code", values_to = "value")
}

cell_values <- read_layer(layer = MEAN_LAYER) |>
  rename(dens_km2 = value) |>
  mutate(dens_km2 = dens_km2 * RASTER_TO_KM2) |>
  left_join(y = rename(read_layer(layer = DD_LAYER), dist_km = value), by = c("cell", "bam_code"))

# ==============================================================================
# 4. Site x species values
# ==============================================================================

# Buffer: area-weighted means over covered cells, ignoring NA cells (water,
# outside the model region); all-NA buffers give NaN, set to NA below.
buffer_values <- buf_cells |>
  inner_join(y = cell_values, by = "cell", relationship = "many-to-many") |>
  group_by(SITEID, bam_code) |>
  summarise(buffer_km2  = weighted.mean(x = dens_km2, w = weight, na.rm = TRUE),
            dist_detect = weighted.mean(x = dist_km,  w = weight, na.rm = TRUE),
            .groups     = "drop")

# Pixel: value of the cell containing the plot coordinate
pixel_values <- pix_cells |>
  inner_join(y = cell_values, by = "cell", relationship = "many-to-many") |>
  select(SITEID, bam_code, pixel_km2 = dens_km2)

# BAM codes mapped to the same species: densities summed, the smaller distance kept
raster_long <- buffer_values |>
  left_join(y = pixel_values, by = c("SITEID", "bam_code")) |>
  left_join(y = distinct(raster_files, bam_code, species), by = "bam_code") |>
  group_by(SITEID, species) |>
  summarise(pixel_km2   = sum(pixel_km2),
            buffer_km2  = sum(buffer_km2),
            dist_detect = min(dist_detect),
            .groups     = "drop") |>
  mutate(across(.cols = c(buffer_km2, dist_detect), .fns = ~ if_else(condition = is.nan(.x), true = NA_real_, false = .x)))

cat("Sites with predictions: ", n_distinct(raster_long$SITEID[!is.na(raster_long$buffer_km2)]),
    " of ", nrow(sites), " | species: ", n_distinct(raster_long$species), "\n", sep = "")

# ==============================================================================
# 5. Save
# ==============================================================================

out <- list(
  raster_long = raster_long,
  raster_info = list(crs = terra::crs(r0, describe = TRUE)$name, resolution = terra::res(r0), layers = names(r0))
)
attr(out, "settings") <- list(
  raster_dir      = RASTER_DIR,
  raster_year     = RASTER_YEAR,
  mean_layer      = MEAN_LAYER,
  dd_layer        = DD_LAYER,
  raster_to_km2   = RASTER_TO_KM2,
  buffer_radius_m = BUFFER_RADIUS_M,
  terra_version   = as.character(packageVersion("terra")),
  created         = Sys.time()
)

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
saveRDS(object = out, file = OUT_PATH)
message("Saved: ", normalizePath(OUT_PATH))
