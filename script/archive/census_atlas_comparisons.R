library(tidyverse)
library(sf)
library(rnaturalearth)
library(naturecounts)
library(ggrepel)

rm(list=ls())

setwd("C:/Users/IlesD/OneDrive - EC-EC/Iles/Projects/Landbirds/Breeding-Bird-Census/script")

# ----------------------------------------
# Functions
# ----------------------------------------

# Function to convert DD.MM.SS format to decimal degrees
dms_to_dd <- function(dms) {
  parts <- unlist(strsplit(dms, "\\."))
  degrees <- as.numeric(parts[1])
  minutes <- as.numeric(parts[2])
  seconds <- ifelse(length(parts) == 3, as.numeric(parts[3]), 0)
  sign <- ifelse(degrees < 0, -1, 1)
  abs(degrees) + (minutes / 60) + (seconds / 3600)
}

# Function to create square polygon from NW corner and area
create_square_plot <- function(point, area_m2) {
  side_length <- sqrt(area_m2)
  coords <- st_coordinates(point)
  nw_x <- coords[1]
  nw_y <- coords[2]
  
  # Define the square polygon (NW to NE to SE to SW to NW)
  square_coords <- matrix(
    c(
      nw_x, nw_y,                        # NW
      nw_x + side_length, nw_y,          # NE
      nw_x + side_length, nw_y - side_length, # SE
      nw_x, nw_y - side_length,          # SW
      nw_x, nw_y                         # Back to NW
    ), ncol = 2, byrow = TRUE
  )
  
  st_polygon(list(square_coords))
}

# ----------------------------------------
# Load spatial boundaries of study areas
# ----------------------------------------

canada <- ne_countries(country = "Canada", returnclass = "sf")

bcr <- sf::st_read("C:/Users/IlesD/OneDrive - EC-EC/Iles/Data/BCR_2025.gdb") %>%
  st_intersection(st_transform(canada,st_crs(.))) %>%
  dplyr::rename(geometry = SHAPE) %>%
  group_by(bcr_label, name_en) %>% 
  summarize(geometry = st_union(geometry), .groups = "drop") %>%
  subset(bcr_label != 23)

# ----------------------------------------
# Load BBC plot information (location / visit years)
# ----------------------------------------

bbc_plots <- readxl::read_xlsx("../data/tables/PLOTS.xlsx") %>%
  dplyr::rename(lon = LONG_W, lat = LAT_N)

# Convert lat/lon columns to decimal degrees
bbc_plots$lat <- sapply(bbc_plots$lat, dms_to_dd)
bbc_plots$lon <- -sapply(bbc_plots$lon, dms_to_dd)

# Area in m2 of each plot
bbc_plots$area_m2 <- bbc_plots$PLOTSIZE * 10000
bbc_plots$area_km2 <- bbc_plots$area_m2 / 1000000

# Convert DataFrame to SF Object
bbc_plots_sf <- st_as_sf(bbc_plots, coords = c("lon", "lat"), crs = 4326) %>% st_transform(3857)

# Only select plots surveyed since 1970
bbc_plots <- subset(bbc_plots, area_m2 > 0.04 & area_km2 < 0.25)

# ----------------------------------------
# Load species count records from BBC
# ----------------------------------------

bbc_dat <- readxl::read_xlsx("../data/tables/SPECIES.xlsx")

# ----------------------------------------
# Load BAM point count dataset
# ----------------------------------------

load("C:/Users/IlesD/OneDrive - EC-EC/Iles/Projects/Landbirds/BBMP/Data/Bird_Data_Raw/BAM/04_NM5.0_data_stratify.Rdata")

# Only use bam data from prior to 2010
#visit <- subset(visit,year <= 2010)
#offsets <- offsets[visit$id,]
#bird <- bird[visit$id,]

bam_sf <- st_as_sf(visit, coords = c("lon","lat"), crs = 4326) %>% st_transform(st_crs(bbc_plots_sf))

# ----------------------------------------
# Specify buffer within which to extract bam data around bbc locations
# ----------------------------------------

# Distance buffer around each bbc location, to extract bam data
distance_threshold <- 5000 # metres
bbc_buff <- bbc_plots_sf %>% st_buffer(distance_threshold)

# Which BAM data are within range of *any* bbc location?
bam_intersects <- st_intersects(bam_sf,bbc_buff)
bam_intersects <- which(lengths(bam_intersects) > 0)
bird2 <- bird[bam_intersects,]
visit2 <- bam_sf[bam_intersects,]
offsets2 <- offsets[bam_intersects,]

# **************************************************************
# **************************************************************
# Comparisons across species
# **************************************************************
# **************************************************************

results_df <- data.frame()

# Identify species in both BAM and BBC
species_list <- search_species_code(colnames(offsets[,-1]))
species_list <- subset(species_list, english_name %in% bbc_dat$COMMONNAME)

passerines <- readxl::read_xlsx("../species/passerines.xlsx") %>%
  subset(Passerine == "Yes")

species_list <- subset(species_list, english_name %in% passerines$Species)

for (species_code in species_list$BSCDATA){
  
  if (species_code == "YWAR") species_code <- "YEWA"
  if (species_code == "AGOL") species_code <- "AMGO"
  
  # ----------------------------------------
  # Select an example species, and compare estimated densities at each bbc plot
  # ----------------------------------------
  
  spec_name <- search_species_code(species_code)
  
  print(species_code)
  if (spec_name$english_name[1] %in% results_df$species_name) next
  
  # Process BBC data
  # Fill in species abundance data for each site in BBC dataset
  spec_data <- subset(bbc_dat, COMMONNAME == spec_name$english_name[1]) %>%
    dplyr::select(RECNO,COMMONNAME,DENSKM2)
  
  spec_data <- bbc_plots_sf %>% left_join(.,spec_data)
  
  # Fill in zeroes
  spec_data$DENSKM2[is.na(spec_data$DENSKM2)] <- 0
  
  # Calculate mean for each location
  spec_data <- spec_data %>%
    group_by(SITEID) %>%
    summarise(
      geometry = st_geometry(geometry)[1],
      bbc_mean = mean(DENSKM2)
    ) %>%
    st_as_sf()
    
  # --------
  # Process BAM data from within distance threshold at each bbc location
  # --------
  
  spec_data <- st_buffer(spec_data,distance_threshold)
  
  spec_data$has_bam <- 0
  spec_data$bam_count_sum <- NA
  spec_data$bam_dens_km2 <- NA
  
  for (i in 1:nrow(spec_data)){
    
    # bbc data for this site
    idat <- bbc_buff[i,]
    
    # bam data assocaited with this site
    ibam_intersects <- st_intersects(visit2,idat)
    ibam_intersects <- which(lengths(ibam_intersects) > 0)
    ibam_counts <- bird2[ibam_intersects,]
    ibam_offsets <- offsets[ibam_intersects,]
    
    if (length(ibam_intersects)>0){
      spec_data$has_bam[i] <- 1
      spec_data$bam_count_sum[i] <- sum(ibam_counts[,species_code])
      spec_data$bam_dens_km2[i] <- sum(ibam_counts[,species_code]/exp(ibam_offsets[,species_code]))
    }
    
  }
  
  # ----------------------------------------
  # plot data
  # ----------------------------------------
  
  ggplot(data = spec_data) +
    geom_point(aes(x = bbc_mean,y = bam_dens_km2)) +
    theme_minimal()
  
  # bbc
  ggplot() +
    geom_sf(data = bcr, col = "gray90", fill = "gray95") +
    geom_sf(data = spec_data %>% st_centroid(), aes(size = bbc_mean),col = "red", fill = "transparent") +
    theme_minimal()
  
  # bam
  ggplot() +
    geom_sf(data = bcr, col = "gray90", fill = "gray95") +
    geom_sf(data = spec_data %>% st_centroid(), aes(size = bam_dens_km2),col = "red", fill = "transparent") +
    theme_minimal()
  
  # ----------------------------------------
  # Calculate means across all plots that were sampled by both methods
  # ----------------------------------------
  
  species_means <- subset(spec_data, !is.na(bbc_mean) & !is.na(bam_dens_km2)) %>%
    as.data.frame() %>%
    summarize(bbc_mean = mean(bbc_mean),
              bam_mean = mean(bam_dens_km2),
              bam_mean_uncorrected = mean(bam_count_sum)) %>%
    mutate(species_name = spec_name$english_name[1]) %>%
    relocate(species_name)
  
  results_df <- rbind(results_df,species_means)
  
  # ----------------------------------------------------
  # Plot results
  # ----------------------------------------------------
  
  result_plot <- ggplot()+
    geom_point(data = results_df, aes(x = bbc_mean, y = bam_mean))+
    geom_text_repel(data = results_df, aes(x = bbc_mean, y = bam_mean, label = species_name), size = 2)+
    
    scale_x_continuous(trans = "log10")+
    scale_y_continuous(trans = "log10")+
    geom_abline(intercept = 0, slope = 1)+
    xlab("Breeding Bird Census (birds/km2)")+
    ylab("Point counts (birds/km2)")+
    theme_minimal()
  
  print(result_plot)
}

# Remove results with unrealistic densities (due to offsets from bam)
results_df2 <- subset(results_df, bam_mean <= 30 & bbc_mean > 0)

result_plot <- ggplot()+
  geom_point(data = results_df2, aes(x = bbc_mean, y = bam_mean))+
  geom_text_repel(data = results_df2, aes(x = bbc_mean, y = bam_mean, label = species_name), size = 2)+
  
  scale_x_continuous(limits = c(0.005,30), trans = "log10")+
  scale_y_continuous(limits = c(0.005,30), trans = "log10")+
  geom_abline(intercept = 0, slope = 1)+
  xlab("Breeding Bird Census (birds/km2)")+
  ylab("Point counts (birds/km2)")+
  ggtitle("Comparison of Breeding Bird Census to Point Counts (+offsets)")+
  theme_bw()
print(result_plot)


mean(results_df2$bam_mean/results_df2$bbc_mean)
median(results_df2$bam_mean/results_df2$bbc_mean)
mean(results_df2$bam_mean > results_df2$bbc_mean)
cor(results_df2$bam_mean,results_df2$bbc_mean, method = "spearman")

mean(results_df2$bam_mean/results_df2$bbc_mean > 2)



# Remove results with unrealistic densities (due to offsets from bam)
results_df2 <- subset(results_df, bam_mean <= 30 & bbc_mean > 2)

result_plot <- ggplot()+
  geom_point(data = results_df2, aes(x = bbc_mean, y = bam_mean))+
  geom_text_repel(data = results_df2, aes(x = bbc_mean, y = bam_mean, label = species_name), size = 2)+
  
  scale_x_continuous(limits = c(1,30), trans = "log10")+
  scale_y_continuous(limits = c(1,30), trans = "log10")+
  geom_abline(intercept = 0, slope = 1)+
  xlab("Breeding Bird Census (birds/km2)")+
  ylab("Point counts (birds/km2)")+
  ggtitle("Comparison of Breeding Bird Census to Point Counts (+offsets)")+
  theme_bw()
print(result_plot)


mean(results_df2$bam_mean/results_df2$bbc_mean)
median(results_df2$bam_mean/results_df2$bbc_mean)
mean(results_df2$bam_mean > results_df2$bbc_mean)
cor(results_df2$bam_mean,results_df2$bbc_mean, method = "spearman")

mean(results_df2$bam_mean/results_df2$bbc_mean > 2)

