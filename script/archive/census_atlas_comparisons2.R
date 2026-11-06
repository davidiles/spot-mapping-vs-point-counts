library(tidyverse)
library(sf)
library(rnaturalearth)
library(naturecounts)
library(ggrepel)
library(viridis)

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

# Only select plots with reasonable area
bbc_plots <- subset(bbc_plots, area_km2 > 0.01 & area_km2 < 0.25) #  & YEAR >= 1960

# Convert DataFrame to SF Object
bbc_plots_sf <- st_as_sf(bbc_plots, coords = c("lon", "lat"), crs = 4326) %>% st_transform(3857)

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

# ----------------------------------------
# For each unique bbc plot, identify which BAM point counts are associated with it
# ----------------------------------------

bbc_plots_buffer <- bbc_plots_sf %>% 
  group_by(SITEID) %>%
  summarize(geometry = st_geometry(geometry)[1],
            first_year = min(YEAR),
            final_year = max(YEAR),
            n_visits_bbc = n(),
            area_km2_bbc = area_km2[1]) %>%
  st_as_sf() %>%
  st_buffer(distance_threshold)

point_indices_by_polygon <- lapply(seq_len(nrow(bbc_plots_buffer)), function(i) {
  which(st_within(visit2, bbc_plots_buffer[i, ], sparse = FALSE))
})

# Number of bam point counts associated with each bbc plot
bbc_plots_buffer$n_points_per_polygon <- lengths(point_indices_by_polygon)

# Plot BAM data availability at each BBC plot
ggplot() +
  geom_sf(data = bcr, col = "gray90", fill = "white") +
  geom_sf(data = bbc_plots_buffer %>% st_centroid() %>% subset(n_points_per_polygon>0), 
          aes(size = area_km2_bbc, col = n_points_per_polygon)) +
  scale_color_gradientn(colours = c("lightblue","darkblue"), trans = "sqrt")+
  theme_minimal()

# **************************************************************
# **************************************************************
# Comparisons across species
# **************************************************************
# **************************************************************

results_df <- data.frame()

# Identify species in both BAM and BBC
species_list <- search_species_code(colnames(offsets[,-1]))
species_list <- subset(species_list, english_name %in% bbc_dat$COMMONNAME)
species_list <- subset(species_list, BSCDATA %in% colnames(bird2))

for (species_code in species_list$BSCDATA){
  
  if (species_code == "YWAR") species_code <- "YEWA"
  if (species_code == "AGOL") species_code <- "AMGO"
  
  # ----------------------------------------
  # Select an example species, and summarize available BBC data at each siteID
  # ----------------------------------------
  
  spec_name <- search_species_code(species_code)
  in_bam <- which(spec_name$BSCDATA %in% colnames(bird2) & spec_name$BSCDATA %in% species_list$BSCDATA) %>% min()
  spec_name <- spec_name[in_bam,]
  
  if (species_code %in% results_df$species_code) next
  
  # Process BBC data
  # Fill in species abundance data for each site/visit in BBC dataset
  spec_data <- subset(bbc_dat, COMMONNAME == spec_name$english_name[1]) %>%
    dplyr::select(RECNO,COMMONNAME,DENSKM2)
  
  spec_data <- bbc_plots_sf %>% left_join(.,spec_data)
  
  # Fill in zeroes
  spec_data$DENSKM2[is.na(spec_data$DENSKM2)] <- 0
  
  # Calculate mean/sum of observed densities for each location
  spec_data <- spec_data %>%
    group_by(SITEID) %>%
    summarise(
      geometry = st_geometry(geometry)[1],
      bbc_mean = mean(DENSKM2),
      bbc_sum = sum(DENSKM2)
    ) %>%
    as.data.frame()
  
  spec_data <- bbc_plots_buffer %>% left_join(spec_data[,c("SITEID","bbc_mean","bbc_sum")])
    
  # --------
  # Associate each site with counts from the BAM dataset
  # --------
  
  spec_data$bam_count_sum <- 0
  spec_data$bam_dens_km2 <- 0
  spec_data$bam_sum_offsets <- 0
  
  for (i in 1:nrow(spec_data)){
    
    if (spec_data$n_points_per_polygon[i] > 0){
      
      # Indices in the bam visit dataset
      ibam <- unlist(point_indices_by_polygon[i])
      ibam_counts <- bird2[ibam,species_code]
      spec_data$bam_count_sum[i] <- sum(ibam_counts)
      spec_data$bam_dens_km2[i] <- sum(ibam_counts/exp(offsets2[ibam,species_code]))
      spec_data$bam_sum_offsets[i] <- sum(exp(offsets2[ibam,species_code]))
      
    }
  }
  
  # ----------------------------------------
  # A SITEID must have at least 1 point count to be considered in this analysis
  # ----------------------------------------
  spec_data <- subset(spec_data, n_points_per_polygon >= 1)
  
  # ----------------------------------------
  # Calculate precision-weighted estimate of means for each survey method
  # ----------------------------------------
  
  # For bbc, precision weights are calculated as number_of_visits x plot_area
  spec_data$bbc_sum_effort <- spec_data$n_visits_bbc * spec_data$area_km2_bbc
  
  # For bam, precision weights are calculated based on sum of offsets
  # spec_data$bam_sum_offsets
  
  # prepare to fit a weighted linear model to assess weighted means of each group
  spec_data_2 <- subset(spec_data, n_points_per_polygon>0)
  
  # Calculate weighted means separately (just for sanity checking)
  weighted.mean(spec_data_2$bbc_sum, w = spec_data_2$bbc_sum_effort)
  weighted.mean(spec_data_2$bam_count_sum, w = spec_data_2$bam_sum_offsets)
  
  # Fit weighted linear model to estimate both simultaneously and conduct test of differences
  x <- c(spec_data_2$bbc_sum, spec_data_2$bam_count_sum)
  group <- factor(rep(c("BBC", "PC"), c(length(spec_data_2$bbc_sum), length(spec_data_2$bam_count_sum))))
  weights <- c(spec_data_2$bbc_sum_effort, spec_data_2$bam_sum_offsets)
  model <- lm(x ~ group, weights = weights)
  summary(model)
  
  predictions <- predict(model, newdata = data.frame(group = c("BBC","PC")), interval = "confidence") %>%
    as.data.frame()
  
  species_summary <- data.frame(species_name = spec_name$english_name,
                                species_code = species_code,
                                bbc_n_obs = sum(spec_data$bbc_mean > 0),
                                bbc_mean = predictions$fit[1],
                                bbc_lcl = predictions$lwr[1],
                                bbc_ucl = predictions$upr[1],
                                
                                pc_n_obs = sum(spec_data$bam_count_sum > 0),
                                pc_mean = predictions$fit[2],
                                pc_lcl = predictions$lwr[2],
                                pc_ucl = predictions$upr[2],
                                
                                p_val = summary(model)$coefficients[2,4]
  )
  
  results_df <- rbind(results_df,species_summary)
  
  # ----------------------------------------------------
  # Plot results
  # ----------------------------------------------------
  
  # result_plot <- ggplot()+
  #   geom_point(data = results_df, aes(x = bbc_mean, y = pc_mean))+
  #   geom_text_repel(data = results_df, aes(x = bbc_mean, y = pc_mean, label = species_name), size = 2)+
  #   
  #   scale_x_continuous(trans = "log10")+
  #   scale_y_continuous(trans = "log10")+
  #   geom_abline(intercept = 0, slope = 1)+
  #   xlab("Breeding Bird Census (birds/km2)")+
  #   ylab("Point counts (birds/km2)")+
  #   theme_minimal()
  # 
  # print(result_plot)
  
}

# ----------------------------------------------------
# Separate results into territorial songbird species that appropriate for monitoring with point counts and spot mapping
# (territorial, vocal, detectable, not flocking)
# ----------------------------------------------------

# Only consider species that were detected in at least 5 plots by both methods
results_df2 <- subset(results_df, bbc_n_obs >= 5 & pc_n_obs >= 5)

# Species list that was analyzed
analyzed_species <- paste0(results_df2$species_name, collapse = ",")


# Of those species, which are considered suitably monitored by both point counts and spot mapping
suitability <- readxl::read_xlsx("../data/tables/bird_survey_suitability.xlsx")

results_df2 <- left_join(results_df2, suitability, by = c("species_name" = "Species"))

# ----------------------------------------------------
# Examine results for all species
# ----------------------------------------------------

lim <- range(results_df2[,c("bbc_mean","pc_mean")])
result_plot2 <- ggplot(data = results_df2)+
  
  geom_abline(intercept =  c(log10(1/5),log10(5)), slope = 1, col = "orange")+
  geom_abline(intercept = 0, slope = 1, col = "black")+
  
  geom_text_repel(aes(x = bbc_mean, y = pc_mean, label = species_name), size = 1, col = "gray70")+
  geom_point(aes(x = bbc_mean, y = pc_mean))+
 
  #scale_x_continuous(limits = lim, trans = "log10")+
  #scale_y_continuous(limits = lim,trans = "log10")+
  scale_x_continuous(limits = lim)+
  scale_y_continuous(limits = lim)+
  
  coord_equal() +
  
  xlab("Breeding Bird Census (birds/km2)")+
  ylab("Point counts (birds/km2)")+
  theme_bw()

print(result_plot2)

mean(results_df2$pc_mean/results_df2$bbc_mean)   # 1.9
median(results_df2$pc_mean/results_df2$bbc_mean) # 0.88
mean(results_df2$pc_mean > results_df2$bbc_mean) # 0.47
cor(results_df2$pc_mean,results_df2$bbc_mean)    # 0.74
mean(abs(log(results_df2$pc_mean/results_df2$bbc_mean)) > log(2)) # 0.57 - over half of species are off by at least a factor of 2
mean(abs(log(results_df2$pc_mean/results_df2$bbc_mean)) > log(5)) # 0.20 - 20% of species are off by more than a factor of 5


# ----------------------------------------------------
# Only plot "suitable" species
# ----------------------------------------------------

results_df3 <- subset(results_df2, suitable == "Yes")
result_plot3 <- ggplot(data = results_df3)+
  
  geom_abline(intercept =  c(log10(1/5),log10(5)), slope = 1, col = "orange")+
  geom_abline(intercept = 0, slope = 1, col = "black")+
  
  geom_text_repel(aes(x = bbc_mean, y = pc_mean, label = species_name), size = 1, col = "gray70")+
  geom_point(aes(x = bbc_mean, y = pc_mean))+
  
  #scale_x_continuous(trans = "log10", limits = lim)+
  #scale_y_continuous(trans = "log10", limits = lim)+
  
  scale_x_continuous(limits = lim)+
  scale_y_continuous(limits = lim)+
  
  coord_equal() +
  xlab("Breeding Bird Census (birds/km2)")+
  ylab("Point counts (birds/km2)")+
  theme_bw()

print(result_plot3)

mean(results_df3$pc_mean/results_df3$bbc_mean)   # 1.3
median(results_df3$pc_mean/results_df3$bbc_mean) # 0.7
mean(results_df3$pc_mean > results_df3$bbc_mean) # 0.4
cor(results_df3$pc_mean,results_df3$bbc_mean)    # 0.74
mean(abs(log(results_df3$pc_mean/results_df3$bbc_mean)) > log(2)) # 0.53 - approx half of species are off by at least a factor of 2
mean(abs(log(results_df3$pc_mean/results_df3$bbc_mean)) > log(5)) # 0.17 - 17 % of species are off by more than a factor of 5

# ----------------------------------------------------
# Calculate precision-weighted correlation
# ----------------------------------------------------


# ----------------------------------------------------
# Compare to results from Hamel
# ----------------------------------------------------

hamel <- readxl::read_xlsx("../data/tables/Hamel_1984_Appendix.xlsx") %>% 
  left_join(suitability) %>%
  subset(suitable == "Yes") %>%
  na.omit()

hamel$Species

setdiff(hamel$Species, results_df$species_name)
hamel_lm <- lm(log_pc~log_sm, data = hamel)

lim <- range(results_df2[,"bbc_mean"],results_df2[,"pc_mean"], hamel$`Spot mapping`,hamel$VPC)
result_plot3.1 <- ggplot(data = results_df3)+
  
  geom_abline(intercept =  c(log10(1/5),log10(5)), slope = 1, col = "black", linetype = 2)+
  geom_abline(intercept = 0, slope = 1, col = "black", linewidth = 2)+
  
  geom_text_repel(aes(x = bbc_mean, y = pc_mean, label = species_name, col = "This analysis"), size = 1)+
  geom_point(aes(x = bbc_mean, y = pc_mean, col = "This analysis"))+
  stat_smooth(aes(x = bbc_mean, y = pc_mean, col = "This analysis", fill = "This analysis"), method = "lm")+
  
  # geom_point(data = hamel, aes(x = `Spot mapping`, y = `VPC`, col = "Hamel (1984)"))+
  # stat_smooth(data = hamel, aes(x = `Spot mapping`, y = `VPC`, col = "Hamel (1984)", fill = "Hamel (1984)"), method = "lm")+
  
  scale_color_manual(values = c("dodgerblue"), name = "Source")+
  scale_fill_manual(values = c("dodgerblue"), name = "Source")+
  
  
  scale_x_continuous(trans = "log10", limits = lim)+
  scale_y_continuous(trans = "log10", limits = lim)+
  coord_equal() +
  xlab("Breeding Bird Census (birds/km2)")+
  ylab("Point counts (birds/km2)")+
  theme_bw()+
  ggtitle("Results from this comparison")

result_plot3.2 <- ggplot(data = results_df3)+
  
  geom_abline(intercept =  c(log10(1/5),log10(5)), slope = 1, col = "black", linetype = 2)+
  geom_abline(intercept = 0, slope = 1, col = "black", linewidth = 2)+
  
  # geom_text_repel(aes(x = bbc_mean, y = pc_mean, label = species_name, col = "This analysis"), size = 1)+
  # geom_point(aes(x = bbc_mean, y = pc_mean, col = "This analysis"))+
  # stat_smooth(aes(x = bbc_mean, y = pc_mean, col = "This analysis", fill = "This analysis"), method = "lm")+
  
  geom_point(data = hamel, aes(x = `Spot mapping`, y = `VPC`, col = "Hamel (1984)"))+
  stat_smooth(data = hamel, aes(x = `Spot mapping`, y = `VPC`, col = "Hamel (1984)", fill = "Hamel (1984)"), method = "lm")+
  
  scale_color_manual(values = c("orangered","dodgerblue"), name = "Source")+
  scale_fill_manual(values = c("orangered","dodgerblue"), name = "Source")+
  
  scale_x_continuous(trans = "log10", limits = lim)+
  scale_y_continuous(trans = "log10", limits = lim)+
  coord_equal() +
  xlab("Breeding Bird Census (pairs / 40 ha)")+
  ylab("Point counts (pairs / 40 ha)")+
  theme_bw()+
  ggtitle("Results from Hamel 1984")

print(result_plot3.1)
print(result_plot3.2)


result_plot3 <- ggplot(data = results_df3)+
  
  geom_abline(intercept =  c(log10(1/5),log10(5)), slope = 1, col = "black", linetype = 2)+
  geom_abline(intercept = 0, slope = 1, col = "black", linewidth = 2)+
  
  geom_text_repel(aes(x = bbc_mean, y = pc_mean, label = species_name, col = "This analysis"), size = 1)+
  geom_point(aes(x = bbc_mean, y = pc_mean, col = "This analysis"))+
  stat_smooth(aes(x = bbc_mean, y = pc_mean, col = "This analysis", fill = "This analysis"), method = "lm")+
  
  geom_point(data = hamel, aes(x = `Spot mapping`, y = `VPC`, col = "Hamel (1984)"))+
  stat_smooth(data = hamel, aes(x = `Spot mapping`, y = `VPC`, col = "Hamel (1984)", fill = "Hamel (1984)"), method = "lm")+
  
  scale_color_manual(values = c("orangered","dodgerblue"), name = "Source")+
  scale_fill_manual(values = c("orangered","dodgerblue"), name = "Source")+
  
  scale_x_continuous(trans = "log10", limits = lim)+
  scale_y_continuous(trans = "log10", limits = lim)+
  coord_equal() +
  xlab("Breeding Bird Census (birds/km2)")+
  ylab("Point counts (birds/km2)")+
  theme_bw()+
  ggtitle("Results for both studies")

result_plot3
