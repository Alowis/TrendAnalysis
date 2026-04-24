# 1. Setup and Libraries
library(ncdf4)
library(sf)
library(raster)
library(dplyr)
library(ggplot2)
library(modifiedmk) # For sens.slope calculation
library(exactextractr)
library(viridis)
library(tidyverse)
library(trend)
source("D:/tilloal/Documents/LFRuns_utils/TrendAnalysis/R/functions_trends.R")

palet2=c(hcl.colors(9, palette = "Blues", alpha = NULL, rev = TRUE, fixup = TRUE))
outletname="GeoData/efas_rnet_100km_01min"
outll=outletopen(hydroDir,outletname)
cord.dec=outll[,c(2,3)]
cord.dec = SpatialPoints(cord.dec, proj4string=CRS("+proj=longlat"))
cord.UTM <- spTransform(cord.dec, CRS("+init=epsg:3035"))
nco=cord.UTM@coords
world <- ne_countries(scale = "medium", returnclass = "sf")
Europe <- world[which(world$continent == "Europe"),]
e2=st_transform(Europe,  crs=3035)
w2=st_transform(world,  crs=3035)
biobase=st_transform(biogeo, crs=3035)
tsize=12
osize=12
Impdates=seq(1950,2020,by=10)
valuenames=paste0("Y",Impdates)
catmap=cst7
basemap=w2

# Directories (Adjust paths as needed)
hydroDir <- "D:/tilloal/Documents/LFRuns_utils/data"
workDir  <- "D:/tilloal/Documents/06_Floodrivers"
selectedHybas <- "2050016510"  # Example ID from your script

# 2. Extract Specific Catchment
hybas_eu <- read_sf(paste0(hydroDir, "/Catchments/hydrosheds/hybas_eu_lev05_v1c.shp"))
target_catchment <- hybas_eu %>% filter(HYBAS_ID == selectedHybas)

hybas_bb <- read_sf(paste0(hydroDir, "/Catchments/hydrosheds/hybas_eu_lev10_v1c.shp"))

# 1. Make the target catchment valid first
target_catchment <- st_make_valid(target_catchment)

# 2. Filter with a 'validity check' on the fly
# This might take a minute because hybas_bb is large
sub_catchments <- hybas_bb %>%
  st_transform(st_crs(target_catchment)) %>%
  st_make_valid() %>% 
  st_filter(target_catchment, .predicate = st_intersects)

sub_catchments <- st_intersection(sub_catchments, st_geometry(target_catchment))


# 3. Compute Sen's Slope for Precipitation and Evapotranspiration
# This function loads yearly files and calculates trend per pixel
calc_climate_trend <- function(prefix, start_yr=1951, end_yr=2020, target_sf,dir) {
  years <- start_yr:end_yr
  file_list <- paste0(dir, "/sums/",prefix,"/", prefix, "_ysum", years, ".nc")
  
  # Load as a stack and crop immediately to the catchment to save memory
  stack_data <- stack(file_list)
  cropped_stack <- mask(crop(stack_data, target_sf), target_sf)
  
  # Sen's Slope calculation function for pixel-wise application
  sen_func <- function(x) {
    if (all(is.na(x))) return(NA)
    res <- modifiedmk::sens.slope(x)
    return(res$estimates)
  }
  
  # Apply trend calculation across the stack
  trend_raster <- calc(cropped_stack, fun = sen_func)
  return(trend_raster)
}

# 3. Fast Trend Calculation Function
# We use sens.slope() which is significantly faster than mmkh() for returning just the trend
get_slope <- function(x) {
  if (any(is.na(x)) || all(x == x[1])) return(NA)
  return(as.numeric(trend::sens.slope(x)$estimates))
}

compute_meteo_trend <- function(prefix, start_yr=1951, end_yr=2020, target_sf,dir) {
  # Build file list and load stack
  years <- start_yr:end_yr
  files <- paste0(dir, "/sums/", prefix,"/",prefix, "_ysum", years, ".nc")
  s <- stack(files)
  
  # CROP FIRST: This is the key to speed. Reducing the number of pixels before calculating.
  s_sub <- mask(crop(s, target_sf), target_sf)
  
  # Calculate trend pixel-wise
  # 'calc' is optimized for this; 'force_slope' returns a single layer
  trend_rast <- calc(s_sub, fun = get_slope)
  return(trend_rast)
}


# 1. Define the Year Lists
yrlist1 <- 1951:1980
yrlist2 <- 1991:2020

# 2. Fast Difference Function
compute_period_diff <- function(prefix, poly,dir) {
  # Generate file paths for both periods
  files1 <- paste0(dir, "/sums/",prefix,"/", prefix, "_ysum", yrlist1, ".nc")
  files2 <- paste0(dir, "/sums/",prefix,"/",prefix, "_ysum", yrlist2, ".nc")
  
  # Load as stacks and crop/mask once to the catchment
  s1 <- mask(crop(stack(files1), poly), poly)
  s2 <- mask(crop(stack(files2), poly), poly)
  
  # Calculate Means (Vectorized and much faster than Sen's Slope)
  mean_p1 <- mean(s1, na.rm = TRUE)
  mean_p2 <- mean(s2, na.rm = TRUE)
  

  mean_p2=mean_p2/4
  mean_p1=mean_p1/4

  # Return the absolute difference between the periods
  return(mean_p2 - mean_p1)
}

compute_aridity_change <- function(poly, dir) {
  # 1. Helper to load and mean a stack
  get_period_mean <- function(prefix, years, poly, dir) {
    files <- paste0(dir, "/sums/", prefix, "/", prefix, "_ysum", years, ".nc")
    s <- mask(crop(stack(files), poly), poly)
    # Applying your /4 correction factor
    return(mean(s, na.rm = TRUE) / 4)
  }
  
  # 2. Load Precipitation (P) and Potential Evapotranspiration (PET)
  # Period 1
  p_p1   <- get_period_mean("pr6", yrlist1, poly, dir)
  pet_p1 <- get_period_mean("et", yrlist1, poly, dir)
  
  # Period 2
  p_p2   <- get_period_mean("pr6", yrlist2, poly, dir)
  pet_p2 <- get_period_mean("et", yrlist2, poly, dir)
  
  # 3. Calculate Aridity Index (AI = P / PET) for each period
  ai_p1 <- p_p1 / pet_p1
  ai_p2 <- p_p2 / pet_p2
  
  # 4. Return the Change (Delta AI)
  # Positive = Getting wetter/less arid | Negative = Getting drier/more arid
  return(ai_p2 - ai_p1)
}

# Run computation for Precip and ET (1951 - 2020)
mdir="D:/tilloal/Documents/06_Floodrivers/meteo/"
# 3. Execute for Precipitation and ET
rast_pr_diff <- compute_period_diff(prefix="pr6", poly=target_catchment,dir=mdir)
rast_et_diff <- compute_period_diff("et", target_catchment,mdir) # Assuming 'et' prefix

rast_delta_ai <- compute_aridity_change(target_catchment, mdir)

# Convert to dataframe for ggplot
df_ai <- as.data.frame(rast_delta_ai, xy = TRUE, na.rm = TRUE)
colnames(df_ai) <- c("x", "y", "delta_ai")
# 4. Convert to Dataframe for ggplot
df_pr <- as.data.frame(rast_pr_diff, xy = TRUE, na.rm = TRUE)
names(df_pr)[3] <- "change"

# Run computation for Precip and ET (1951 - 2020)
mdir="D:/tilloal/Documents/06_Floodrivers/meteo/"
# rast_precip_trend <- compute_meteo_trend(prefix="pr6", start_yr=1951, end_yr=2020, target_sf=target_catchment,dir=mdir)
# rast_et_trend     <- calc_climate_trend("et", 1951, 2020, target_catchment)


# 4. Extract Land Use Change and Water Demand
# Land use: Forest and Sealed Fractions
rast_forest <- raster(paste0(workDir, "/landuse/fracforest_ch20201951.tif"))
rast_sealed <- raster(paste0(workDir, "/landuse/fracsealed_ch20201951.tif"))

forest_crop <- mask(crop(rast_forest, target_catchment), target_catchment)
sealed_crop <- mask(crop(rast_sealed, target_catchment), target_catchment)


# Water Demand: Aggregated at NUTS3 level
nuts3_shp <- read_sf(paste0(hydroDir, "/Countries/NUTS3/NUTS3_modified.shp"))
rast_wd_total <- raster(paste0(workDir, "/wateruse/all_ysum_ch20201951.tif"))
rast_wd_2020   <- raster(paste0(workDir, "/wateruse/wateruse_sums/all_demands_2020.tif"))

# 1. Ensure CRS matches
if (st_crs(nuts3_shp) != st_crs(target_catchment)) {
  nuts3_shp <- st_transform(nuts3_shp, st_crs(target_catchment))
}

# 2. Calculate centroids for ALL NUTS3 regions
# Note: Do this on the original 'nuts3_shp', not the intersected one
all_nuts3_centroids <- st_centroid(nuts3_shp)

# 3. Find which centroids are inside the target catchment
# st_intersects returns a logical vector of which points are inside the polygon
is_inside <- st_intersects(all_nuts3_centroids, target_catchment, sparse = FALSE)

# 4. Filter the ORIGINAL polygons based on that centroid check
# This keeps the full, un-chopped NUTS3 shapes
nuts3_retained <- nuts3_shp[as.vector(is_inside), ]
# 4. Now calculate your Water Demand statistics on these specific polygons
# Using 'sum' on the full NUTS3 polygons
nuts3_retained$wd_change <- exact_extract(rast_wd_total, nuts3_retained, 'sum')
nuts3_retained$total_wd_2020 <- exact_extract(rast_wd_2020, nuts3_retained, 'sum')

# 5. Create a centroid object for plotting (bubbles) for these retained regions
nuts3_centroids <- st_centroid(nuts3_retained) %>%
  mutate(lon = st_coordinates(.)[,1],
         lat = st_coordinates(.)[,2])

# 5. Extract New Reservoirs
# Using the logic from ReservoirOpen to find points within catchment

res2020=resOpen(hydroDir,"/reservoirs/reservoirs_volumes_2020_Domain2.nc")
res2020$idla=2970-res2020$idla+1
res2020$idlalo=paste(res2020$idlo,res2020$idla,sep=" ")
res1951=resOpen(hydroDir,"/reservoirs/reservoirs_volumes_1951.nc")
max(res2020$res)
matres=na.omit(match(res1951$idlalo,res2020$idlalo))
res_old=res2020[matres,]
res_new=res2020[-matres,]
res_comp=left_join(res_old,res1951,by="idlalo")

new_res_sf <- st_as_sf(res_new, coords = c("Var1", "Var2"), crs = 4326) # Use the CRS of your meteo data (usually WGS84)

# 3. Ensure CRS matches target_catchment to avoid the intersection error
if (st_crs(new_res_sf) != st_crs(target_catchment)) {
  new_res_sf <- st_transform(new_res_sf, st_crs(target_catchment))
}

# 4. Perform the intersection
new_res_points <- st_intersection(new_res_sf, target_catchment)


# 6. Final Multi-Layer Map
# Convert rasters to dataframes for ggplot compatibility
precip_df <- as.data.frame(rast_precip_trend, xy = TRUE, na.rm = TRUE)

forest_df <- as.data.frame(forest_crop, xy = TRUE, na.rm = TRUE)
forest_df$fracforest_ch20201951[which(forest_df$fracforest_ch20201951==0)]=NA
sealed_df <- as.data.frame(sealed_crop, xy = TRUE, na.rm = TRUE)
sealed_df$fracsealed_ch20201951[which(sealed_df$fracsealed_ch20201951==0)]=NA

df_forest_high <- forest_df %>% 
  filter(fracforest_ch20201951 > 0.2)

df_sealed_high <- sealed_df %>% 
  filter(fracsealed_ch20201951 > 0.2)

ggplot() +
  # Background: Precipitation Trend (Raster)
  # geom_raster() +
  # scale_fill_viridis_c(option = "mako", name = "Precip Trend") +
  
  # 2. Forest Change: Diverging Fill + Dynamic Alpha
  # Alpha is mapped to the absolute value so 0 is transparent
  geom_tile(data = forest_df, 
            aes(x = x, y = y, 
                fill = fracforest_ch20201951, 
                alpha = abs(fracforest_ch20201951)), 
            linewidth = 0) +
  scale_fill_gradient2(low = "#8c510a", mid = "white", high = "#01665e", 
                       midpoint = 0, na.value = "transparent", name = "Forest Change") +
  scale_alpha_continuous(range = c(0.3, 1), guide = "none") + # guide="none" hides alpha legend
  
  # --- SECOND FILL SCALE RESET ---
  new_scale_fill() +
  
  # 3. Sealed Change: Diverging Fill + Dynamic Alpha
  geom_tile(data = sealed_df, 
            aes(x = x, y = y, 
                fill = fracsealed_ch20201951, 
                alpha = abs(fracsealed_ch20201951)+0.4), 
            linewidth = 0) +
  scale_fill_gradient2(low = "#5e3c99", mid = "white", high = "#e66101", 
                       midpoint = 0, na.value = "transparent", name = "Sealed Change") +
  scale_alpha_continuous(range = c(0.3, 1), guide = "none") +
  # Water Demand: Circles centered on NUTS3
  # size = wd_change maps the volume of change to the radius of the circle
  # Water Demand Circles
  # SIZE = Total 2020, COLOR = Change 1951-2020
  geom_sf(data = nuts3_centroids, 
          aes(size = total_wd_2020, color = wd_change), 
          alpha = 0.8) +
  scale_size_continuous(range = c(2, 12), name = "Total Demand 2020") +
  scale_color_gradient2(low = "blue", mid = "grey", high = "red", name = "Demand Change") +
  # New Reservoirs (Points)
  geom_sf(data = new_res_points, color = "black", size = 5, shape = 18) +
  
  # Catchment Boundary
  geom_sf(data = target_catchment, fill = "transparent", color = "gray23", size = 1) +
  
  theme_minimal() +
  labs(title = paste("Catchment Changes (ID:", selectedHybas, ")"),
       subtitle = "Precip (Raster), Land Use (Pixels), Water Demand (NUTS3), Reservoirs (Points)")

# We resample 'et' to match 'pr' to fix the "different extent" error
rast_et_resampled <- resample(rast_et_diff, rast_pr_diff, method = "bilinear")

climate_stack <- stack(rast_pr_diff, rast_et_resampled)
names(climate_stack) <- c("P_change", "ET_change")


# 1. Define your manual thresholds (3 classes require 2 internal split points)
# Anything below -20 is "1", -20 to 20 is "2", above 20 is "3"
p_cuts  <- c(-20, 20)
et_cuts <- c(-20, 20)

# 2. Manually assign classes 1, 2, and 3
# findInterval returns 1, 2, or 3 based on where the value sits
df_climate$P_manual_class  <- findInterval(df_climate$P_change, p_cuts, left.open = FALSE) + 1
df_climate$ET_manual_class <- findInterval(df_climate$ET_change, et_cuts, left.open = FALSE) + 1

# 3. Constrain values to 1-3 (in case of NAs or extremes)
df_climate$P_manual_class  <- pmin(pmax(df_climate$P_manual_class, 1), 3)
df_climate$ET_manual_class <- pmin(pmax(df_climate$ET_manual_class, 1), 3)
table(df_climate$ET_manual_class)
# 4. Create the final 'bi_class' string that ggplot needs
# This creates the "1-1", "2-3", etc. labels manually
df_bi <- df_climate %>%
  mutate(bi_class = paste0(P_manual_class, "-", ET_manual_class))

# 5. Check your work
table(df_bi$bi_class)
# Create categorical rasters (1, 2, 3)
# 'findInterval' works on vectors, so we wrap it in 'calc' or use raster math
rast_p_class <- reclassify(rast_pr_diff, c(-Inf, p_cuts[1], 1,  p_cuts[1], p_cuts[2], 2,  p_cuts[2], Inf, 3))
rast_et_class <- reclassify(rast_et_final, c(-Inf, et_cuts[1], 1,  et_cuts[1], et_cuts[2], 2,  et_cuts[2], Inf, 3))

# 4. Apply your Manual Bivariate Breaks
# Adjust these cuts based on your specific climate range

sub_catchments$P_change<- exact_extract(rast_pr_diff, sub_catchments, 'mean')
sub_catchments$ET_change <- exact_extract(rast_et_diff, sub_catchments, 'mean')
sub_catchments$AI_change <- exact_extract(rast_delta_ai, sub_catchments, 'mean')

sub_catchments <- sub_catchments %>%
  mutate(
    P_class = findInterval(P_change, p_cuts) + 1,
    ET_class = findInterval(ET_change, et_cuts) + 1,
    bi_class = paste0(pmin(pmax(P_class, 1), 3), "-", pmin(pmax(ET_class, 1), 3))
  )


#load bivariate changes
hydroDir<-("D:/tilloal/Documents/LFRuns_utils/ChangingHydroExtremes/data")

load (file=paste0(hydroDir,"/databivar_clim.Rdata"))
load (file=paste0(hydroDir,"/dataSEchanges_bivariate.Rdata"))

head(databipic)

# 1. Transform the catchment to match the points (LAEA Europe)
# This is usually more accurate for spatial operations in Europe
databipic_wgs84<- st_transform(databipic, st_crs(target_catchment))

# 2. Crop the points to the catchment boundary
# This acts like a 'cookie cutter'
databipic_cropped <- st_intersection(databipic_wgs84, target_catchment)

#removal of pixel with NA direction
databipic_cropped$bi_class[which(databipic_cropped$bi_class=="NA-NA")]=NA
databipic_cropped$bi_class[which(is.na(databipic_cropped$Y2015))]=NA
databipic_cropped$bi_class[which(is.na(databipic_cropped$d2015))]=NA

# 2. Crop the points to the catchment boundary
# This acts like a 'cookie cutter'
databipise_wgs84<- st_transform(databipise, st_crs(target_catchment))
databipise_cropped <- st_intersection(databipise_wgs84, target_catchment)

databipise_cropped$bi_class[which(databipise_cropped$bi_class=="NA-NA")]=NA
databipise_cropped$bi_class[which(is.na(databipise_cropped$Y2015))]=NA
databipise_cropped$bi_class[which(is.na(databipise_cropped$d2015))]=NA

# 3. (Optional) Check how many points are left
print(paste("Points before:", nrow(databipic), "| Points after:", nrow(databipic_cropped)))



# 3. Final Multi-Layer Map
library(ggnewscale)
library(biscale)
library(cowplot)
# 4. Prepare Other Layers (Water Demand & Land Use)
# [Assuming previous logic for nuts3_centroids, forest_df, and sealed_df is loaded]




# 5. Build the Bivariate Map
ggplot() +
  # LAYER 1: Bivariate Climate (P and ET combined)
  geom_sf(data = sub_catchments, aes(fill= AI_change),color="transparent") +
  
  # ggplot() +
  #   # Background: Precipitation Trend (Raster)
  #   # 1. Background: Precipitation Trend
  #   geom_raster(data = df_pr, aes(x = x, y = y, fill = change)) +
  scale_fill_viridis_c(option = "mako", name = "Precip Trend") +
  
  # Catchment Boundary
  geom_sf(data = target_catchment, fill = "transparent", color = "gray23", size = 1) +
  
  
  # 2. Forest Change: Diverging Fill + Dynamic Alpha
  # Alpha is mapped to the absolute value so 0 is transparent
  theme_minimal() +
  labs(title = paste("Catchment Changes (ID:", selectedHybas, ")"),
       subtitle = "Precip (Raster), Land Use (Pixels), Water Demand (NUTS3), Reservoirs (Points)")

# 1. Define 5 values for 4 classes
p_breaks_4  <- c(-50, -10, 0, 10, 50)
et_breaks_4 <- c(0, 20, 50,100)

# 2. Update the findInterval logic
sub_catchments$P_class <- findInterval(sub_catchments$P_change, p_breaks_4, all.inside = TRUE)
sub_catchments$ET_class <- findInterval(sub_catchments$ET_change, et_breaks_4, all.inside = TRUE)

summary(sub_catchments$ET_change)
# 3. Create the bi_class string (now ranging from 1-1 to 4-4)
df_bi <- sub_catchments %>%
  mutate(bi_class = paste0(P_class, "-", ET_class))



colorp <- c(
  "1-1" = "#dd6a40",  # high x, low y
  "2-1" = "#d9926a",  # medium-high x, low y
  "3-1" = "#DEB887",  # medium-low x, low y
  "4-1" = "#FFD39B",  # low x, low y
  
  "1-2" = "#a36229",  # high x, medium-low y
  "2-2" = "#999999",  # medium-high x, medium-low y
  "3-2" = "#999999",  # medium-low x, medium-low y
  "4-2" = "#9cc4d2",  # low x, medium-low y
  
  "1-3" = "#635929",  # high x, medium-high y
  "2-3" = "#999999",  # medium-high x, medium-high y
  "3-3" = "#999999",  # medium-low x, medium-high y
  "4-3" = "#5fb2d1",  # low x, medium-high y
  
  "1-4" = "#174f28",  # high x, high y
  "2-4" = "#166d68",  # medium-high x, high y
  "3-4" = "#16869e",  # medium-low x, high y
  "4-4" = "#169dd0"  # low x, high y
  
)

loscolors=c("Accelerating" = "#174f28","Drying" = "#dd6a29","Stable"="gray60","Wetting" = "#169dd0","Decelerating" = "burlywood")

bi_pal(pal = colors, dim = 4)


legend <- bi_legend(pal = colorp,
                    dim = 4,
                    xlab = "  +  Drought intensity  -  ",
                    ylab = "  -  Flood intensity  +  ",
                    size = 16,
                    arrows = FALSE)

pl=ggarrange(map, legend,
             labels = c("Map", "Key"),
             ncol = 2, nrow = 1,widths = c(2,1), heights=c(1,1), vjust=-1)


ggplot(Europe) +
  geom_sf(fill="white")+
  # Background: Precipitation Trend (Raster)
  # geom_raster() +
  # scale_fill_viridis_c(option = "mako", name = "Precip Trend") +
  
  # 2. Forest Change: Diverging Fill + Dynamic Alpha
  # Alpha is mapped to the absolute value so 0 is transparent
  # LAYER 1: Bivariate Climate (P and ET combined)
  # geom_sf(data = df_bi, aes( fill = bi_class), show.legend = FALSE,color="transparent") +
  # # Catchment Boundary
  # geom_sf(data = target_catchment, fill = "transparent", color = "gray23", size = 2) +
  # 
  # 
  # bi_scale_fill(pal = "DkViolet2", dim = 4,flip_axes = TRUE) + # "BlueOr" fits your Red/Blue intent
  # new_scale_fill() +
  # ggplot() +
  #   # Background: Precipitation Trend (Raster)
  #   # 1. Background: Precipitation Trend
  #   geom_raster(data = df_pr, aes(x = x, y = y, fill = change)) +
# Layer 2: The Cropped Points (databipic)
# Using size to show the magnitude of change (d2015)
# Catchment Boundary
geom_sf(data = target_catchment, fill = "white", color = "gray23", size = 1) +
  geom_sf(data = databipise_cropped, aes(color = bi_class),size=1,alpha=1,stroke=0,shape=15) +
  
  scale_color_manual(values = colorp, guide="none") +
  
  geom_tile(data = df_forest_high, 
            aes(x = x, y = y), 
            fill="#01665e",
                alpha = .7, 
            linewidth = 0.1) +
  # scale_fill_gradient2(low = "#8c510a", mid = "white", high = "#01665e", 
  #                      midpoint = 0, na.value = "transparent", name = "Forest Change") +
  # scale_alpha_continuous(range = c(0.8, 1), guide = "none") + # guide="none" hides alpha legend
  
  # --- SECOND FILL SCALE RESET ---
  new_scale_fill() +
  
  # 3. Sealed Change: Diverging Fill + Dynamic Alpha
  geom_tile(data = df_sealed_high, 
            aes(x = x, y = y),fill="brown", 
            linewidth = 0,alpha=.7) +
  # scale_fill_gradient2(low = "#5e3c99", mid = "white", high = "#e66101", 
  #                      midpoint = 0, na.value = "transparent", name = "Sealed Change") +
  # scale_alpha_continuous(range = c(0.8, 1), guide = "none") +
  # Water Demand: Circles centered on NUTS3
  # size = wd_change maps the volume of change to the radius of the circle
  # Water Demand Circles
  # SIZE = Total 2020, COLOR = Change 1951-2020
  
  
  # --- SECOND FILL SCALE RESET ---
  new_scale_color() +
  new_scale_fill() +
  
  geom_sf(data = nuts3_centroids, 
          aes(size = total_wd_2020/1000, color = wd_change/1000), 
          alpha = 0.6) +
  scale_size_continuous(range = c(2, 12), name = "Total Demand 2020 (km3)") +
  scale_color_gradient2(low = "blue", mid = "grey", high = "red", name = "Demand Change (km3)") +
  # New Reservoirs (Points)
  
  geom_sf(data = new_res_points, color = "black", size = 5, shape = 18) +


  theme(axis.title=element_text(size=tsize),
        panel.background = element_rect(fill = "aliceblue", colour = "grey1"),
        panel.border = element_rect(linetype = "solid", fill = NA, colour="black"),
        legend.title = element_text(size=tsize),
        legend.text = element_text(size=osize),
        legend.position = "right",
        panel.grid.minor = element_line(colour = "grey90"),
        legend.key = element_rect(fill = "transparent", colour = "transparent"),
        legend.key.size = unit(.8, "cm"))+
  
  labs(title = paste("Socioeconomic changes | Rhone catchment"),
       subtitle = "Land Use (Pixels), Water Demand (NUTS3), Reservoirs (Points)")


library(ggnewscale) # Essential for multiple color/fill scales

# 1. Calculate the bounding box of your target catchment for the zoom
bbox <- st_bbox(target_catchment)

ggplot(Europe) +
  # --- BASE LAYER: Europe Background ---
  geom_sf(fill = "white", color = "grey80") +
  
  # --- LAYER 1: Catchment Area ---
  geom_sf(data = target_catchment, fill = "white", color = "gray23", size = 1) +
  
  # --- LAYER 2: Climate Points (databipise) ---
  geom_sf(data = databipise_cropped, aes(color = bi_class), size = 1, alpha = 1, shape = 15) +
  scale_color_manual(values = colorp, guide = "none") + 
  
  # --- Layer 2: Land Use Changes (High Contrast) ---
  new_scale_fill() +
  geom_tile(data = df_forest_high, aes(x = x, y = y, fill = "Forest Increase (>0.2)"), alpha = 0.8) +
  geom_tile(data = df_sealed_high, aes(x = x, y = y, fill = "Sealed Increase (>0.2)"), alpha = 0.8) +
  scale_fill_manual(
    name = "Land Use Change",
    values = c("Forest Increase (>0.2)" = "#00441b", # Deep Emerald Green
               "Sealed Increase (>0.2)" = "#762a83")  # Deep Purple/Magenta
  ) +
  
  # --- Layer 3: Water Demand (NUTS3 Bubbles) ---
  new_scale_color() +
  geom_sf(data = nuts3_centroids, 
          aes(size = total_wd_2020/1000, color = wd_change/1000), 
          alpha = 0.7) +
  scale_size_continuous(range = c(2, 12), name = "Total Demand 2020 (km³)") +
  scale_color_gradient(low = "grey30", high = "black", name = "Demand Change (km³)") +
  
  # --- Layer 4: Reservoirs (The "Water" Points) ---
  new_scale_fill() +
  geom_sf(data = new_res_points, aes(fill = "New Reservoir"), size = 4, shape = 23,color="black",stroke=1) +
  scale_fill_manual(name = "Infrastructure", values = c("New Reservoir" = "darkblue")) + # Bright Blue
  
  # --- Layer 5: Catchment Boundary ---
  geom_sf(data = target_catchment, fill = NA, color = "black", size = 0.8) +
  
  # --- Zoom and Styling ---
  coord_sf(xlim = c(bbox["xmin"], bbox["xmax"]), 
           ylim = c(bbox["ymin"], bbox["ymax"]), 
           expand = FALSE) +
  theme_minimal() +
  # --- Styling ---
  theme_minimal() +
  theme(axis.title = element_text(size = tsize),
        panel.background = element_rect(fill = "aliceblue", colour = "grey1"),
        panel.border = element_rect(linetype = "solid", fill = NA, colour = "black"),
        legend.title = element_text(size = tsize, face = "bold"),
        legend.text = element_text(size = osize),
        # This removes the boxes and borders from the legend icons
        legend.key = element_blank(), 
        plot.title = element_text(size = tsize + 2, face = "bold"),
        plot.subtitle = element_text(size = osize),
        # This ensures the legend background itself doesn't have a white box
        legend.background = element_blank(),
        legend.position = "right",
        legend.box = "vertical", # Stacks the multiple legends
        panel.grid.major = element_line(colour = "grey90"))+
  
  labs(title = "Socioeconomic changes | Rhone catchment",
       subtitle = "Land Use (Pixels), Water Demand (NUTS3), Reservoirs (Points)",
       x = "Longitude", y = "Latitude")

library(sf)
library(terra)
library(dplyr)
library(ggnewscale)

# Define Target CRS (LAEA Europe)
target_crs <- 3035

# A. Reproject Vector Layers
Europe_laea           <- st_transform(Europe, target_crs)
target_catchment_laea <- st_transform(target_catchment, target_crs)
databipise_laea       <- st_transform(databipise_cropped, target_crs)
nuts3_laea            <- st_transform(nuts3_centroids, target_crs)
res_laea              <- st_transform(new_res_points, target_crs)
databipise_claea       <- st_transform(databipise_cropped, target_crs)
sub_catchments_laea <- st_transform(sub_catchments, target_crs)
library(terra)

# 1. Convert RasterLayer to SpatRaster
rast_spat <- rast(rast_delta_ai)

# 2. Reproject using the EPSG code
rast_laea <- project(rast_spat, "EPSG:3035")

# 3. Convert to dataframe for ggplot
df_ai_laea <- as.data.frame(rast_laea, xy = TRUE, na.rm = TRUE)
colnames(df_ai_laea) <- c("x", "y", "delta_ai")
# C. Calculate the LAEA Bounding Box for the Zoom
bbox_laea <- st_bbox(target_catchment_laea)

## 1. Transform to SF and reproject
forest_sf_laea <- df_forest_high %>%
  st_as_sf(coords = c("x", "y"), crs = 4326) %>%
  st_transform(3035)

# 2. Extract the new coordinates and convert BACK to a dataframe
df_forest_laea <- forest_sf_laea %>%
  mutate(x = st_coordinates(.)[,1],
         y = st_coordinates(.)[,2]) %>%
  st_drop_geometry() # This removes the 'spatial' part so it's a regular df

# 3. Repeat for Sealed
df_sealed_laea <- df_sealed_high %>%
  st_as_sf(coords = c("x", "y"), crs = 4326) %>%
  st_transform(3035) %>%
  mutate(x = st_coordinates(.)[,1],
         y = st_coordinates(.)[,2]) %>%
  st_drop_geometry()

cities <- data.frame(
  name = c("Lyon", "Geneva", "Marseille", "Avignon", "Grenoble", 
           "Valence", "Bourg-en-Bresse", "Chambéry", "Dijon", "Alès"),
  lon = c(4.8357, 6.1432, 5.3698, 4.8055, 5.7245, 
          4.8924, 5.2272, 5.9175, 5.0222, 4.0328),
  lat = c(45.7640, 46.2044, 43.2965, 43.9488, 45.1885, 
          44.9334, 46.2052, 45.5646, 47.3220, 44.1272)
)
# Convert to LAEA (3035) and extract coordinates for ggrepel
cities_laea<- cities %>%
  st_as_sf(coords = c("lon", "lat"), crs = 4326) %>%
  st_transform(3035)
# Extract coordinates for repel
cities_coords <- cities_laea %>%
  mutate(x = st_coordinates(.)[,1],
         y = st_coordinates(.)[,2]) %>%
  st_drop_geometry()

library(ggrepel)
ggplot(Europe_laea) +
  # --- Layer 0: Background Europe ---
  geom_sf(fill = "grey95", color = "white") +
  
  # --- LAYER 2: Climate Points (databipise) ---
  geom_sf(data = databipise_claea, aes(color = bi_class), size = 1, alpha = 1, shape = 15) +
  scale_color_manual(values = colorp, guide = "none") + 
  # --- Layer: Main Cities ---
  geom_sf(data = cities_laea, color = "black", size = 2) +
  
  # --- Layer 2: Land Use Changes (High Contrast) ---
  new_scale_fill() +
  geom_tile(data = df_forest_laea, aes(x = x, y = y, fill = "Forest Increase (>0.2)"),
            width = 2000, height = 2000,alpha = 0.8) +
  geom_tile(data = df_sealed_laea, aes(x = x, y = y, fill = "Sealed Increase (>0.2)"),
            width = 2000, height = 2000, alpha = 0.8) +
  scale_fill_manual(
    name = "Land Use Change",
    values = c("Forest Increase (>0.2)" = "green3", # Deep Emerald Green
               "Sealed Increase (>0.2)" = "maroon3") 
  ) +
  
  
  # --- Layer 3: Water Demand (NUTS3 Bubbles) ---
  new_scale_color() +
  geom_sf(data = nuts3_laea, 
          aes(size = total_wd_2020/1000, color = wd_change/1000), 
          alpha = 0.7) +
  scale_size_continuous(range = c(2, 12), name = "Total Water Demand 2020 (km³)",guide = guide_legend(
    direction = "horizontal", 
    title.position = "top", 
    nrow = 1, 
    order = 3
  )) +
  scale_color_gradient(low = "grey60", high = "grey20", name = "Water Demand Change (km³)",guide = guide_colorbar(
    direction = "horizontal", 
    title.position = "top", 
    barwidth = 10,         # You can control the width of the bar here
    order = 4
  )) +
  
  # --- Layer 4: Reservoirs (The "Water" Points) ---
  # --- Layer 4: Reservoirs (The "Water" Points) ---
  new_scale_fill() +
  geom_sf(data = res_laea, aes(fill = "New Reservoir"), size = 4, shape = 23,color="black",stroke=1) +
  scale_fill_manual(name = "Infrastructure", values = c("New Reservoir" = "darkblue"),guide = guide_legend(
    direction = "horizontal", 
    title.position = "top", 
    order = 2
  )) + # Bright Blue
  

  
  # --- Layer: City Labels with Transparent Boxes ---
  geom_label_repel(data = cities_coords, 
                   aes(x = x, y = y, label = name),
                   size = 3.5, 
                   fontface = "bold",
                   # Fill: White (#FFFFFF) with ~50% transparency (88)
                   fill = "#FFFFFF88", 
                   # Label border: Set to 0 to remove it, or a small number for a thin line
                   label.size = 0.2,
                   label.padding = unit(0.15, "lines"),
                   box.padding = 0.5,
                   # Ensure the label doesn't cover the city point
                   point.padding = 0.3,
                   segment.color = "black", # Line connecting box to point
                   segment.size = 0.3) +
  
  # --- Layer 5: Catchment Boundary ---
  geom_sf(data = target_catchment_laea, fill = NA, color = "black", size = 0.8) +
  
  # --- Zoom and Styling ---
  coord_sf(xlim = c(bbox_laea["xmin"], bbox_laea["xmax"]), 
           ylim = c(bbox_laea["ymin"], bbox_laea["ymax"]), 
           expand = FALSE) +

  # --- Styling ---
  theme_minimal() +
  theme(axis.title = element_text(size = tsize),
        panel.background = element_rect(fill = "aliceblue", colour = "grey1"),
        panel.border = element_rect(linetype = "solid", fill = NA, colour = "black"),
        legend.title = element_text(size = tsize, face = "bold"),
        legend.text = element_text(size = osize),
        # This removes the boxes and borders from the legend icons
        legend.key = element_blank(), 
        plot.title = element_text(size = tsize + 2, face = "bold"),
        plot.subtitle = element_text(size = osize),
        # This ensures the legend background itself doesn't have a white box
        legend.background = element_blank(),
        legend.position = "right",
        legend.box = "vertical", # Stacks the multiple legends
        panel.grid.major = element_line(colour = "grey90"))+
  
  labs(title = "Socioeconomic changes | Rhone catchment",
       subtitle = "Land Use (Pixels), Water Demand (NUTS3), Reservoirs (Points)",
       x = "Longitude", y = "Latitude")



s_plot=ggplot(Europe_laea) +
  # --- Layer 0: Europe Background ---
  geom_sf(fill = "white", color = "grey80") +
  
  # --- Layer 1: Climate/AI Background ---
  geom_sf(data = databipise_claea, aes(color = bi_class), size = 1, alpha = 1, shape = 15) +
  scale_color_manual(values = colorp, guide = "none") + 
  
  # --- Layer: Cities & Catchment ---
  geom_sf(data = cities_laea, color = "black", size = 2) +

  
  # --- Layer 2: Land Use Change (LEGEND ORDER 1) ---
  new_scale_fill() +
  geom_tile(data = df_forest_laea, aes(x = x, y = y, fill = "Forest Increase (>0.2)"),
            width = 2000, height = 2000, alpha = 0.7) +
  geom_tile(data = df_sealed_laea, aes(x = x, y = y, fill = "Sealed Increase (>0.2)"),
            width = 2000, height = 2000, alpha = 0.7) +
  scale_fill_manual(
    name = "Land use change",
    values = c("Forest Increase (>0.2)" = "green3", 
               "Sealed Increase (>0.2)" = "maroon3"),
    guide = guide_legend(order = 1) # <--- FIRST
  ) +
  
  # --- Layer 3: Reservoirs (LEGEND ORDER 2) ---
  new_scale_fill() +
  geom_sf(data = res_laea, aes(fill = "New reservoirs"), size = 4, shape = 23, color = "black", stroke = 1) +
  scale_fill_manual(
    name = "Infrastructure", 
    values = c("New reservoirs" = "darkblue"),
    guide = guide_legend(order = 2) # <--- SECOND
  ) +
  
  # --- Layer 4: Water Demand (LEGEND ORDER 3) ---
  new_scale_color() +
  geom_sf(data = nuts3_laea, 
          aes(size = total_wd_2020/1000, color = wd_change/1000), 
          alpha = 0.7) +
  scale_size_continuous(range = c(2, 12), name = "Total water demand 2020 (km³)",
                        guide = guide_legend(
                          direction = "horizontal", 
                          title.position = "top", 
                          nrow = 1, 
                          order = 3
                        )) + # <--- THIRD (Size)
  scale_color_gradient(low = "grey70", high = "grey30", name = "Water demand change (km³)",
                       guide = guide_colorbar(
                         direction = "horizontal", 
                         title.position = "top", 
                         barwidth = 10,         # You can control the width of the bar here
                         order = 4
                       )) + # <--- FOURTH (Color)
  
  # --- Layer: Cities & Catchment ---
  # geom_sf(data = cities_laea, color = "black", size = 2) +
  geom_label_repel(data = cities_coords, aes(x = x, y = y, label = name),
                   size = 3.5, fontface = "bold", fill = "#FFFFFF88",
                   label.size = 0.2, label.padding = unit(0.15, "lines"),
                   box.padding = 0.5, point.padding = 0.3,
                   segment.color = "black", segment.size = 0.3) +
  geom_sf(data = target_catchment_laea, fill = NA, color = "black", size = 0.8) +
  
  # --- Zoom & Final Theme ---
  coord_sf(xlim = c(bbox_laea["xmin"], bbox_laea["xmax"]), 
           ylim = c(bbox_laea["ymin"], bbox_laea["ymax"]), 
           expand = FALSE) +
  theme_minimal() +
  theme(
    axis.title = element_text(size = tsize),
    panel.background = element_rect(fill = "aliceblue", colour = "grey1"),
    panel.border = element_rect(linetype = "solid", fill = NA, colour = "black"),
    legend.title = element_text(size = tsize, face = "bold"),
    legend.text = element_text(size = osize),
    legend.key = element_blank(), 
    legend.background = element_blank(),
    legend.position = "right",
    legend.box = "vertical",
    plot.title = element_text(size = tsize + 2, face = "bold"),
    plot.subtitle = element_text(size = osize),
    panel.grid.major = element_line(colour = "grey90")
  ) +
  labs(title = " (b) Socioeconomic changes",
       subtitle = "Land use (Pixels), water demand (NUTS3), reservoirs (Points)",
       x = "Longitude", y = "Latitude")

ggsave(paste0("D:/tilloal/Documents/LFRuns_utils/TrendAnalysis/plots/rhone_sechanges3.jpg"), width=20, height=20, units=c("cm"),dpi=1000)




#climate:


# 4. Plotting
# Note: You MUST use a palette that supports dim = 4
ggplot() +
  geom_sf(data = sub_catchments_laea, aes(fill = AI_change)) +
  bi_scale_fill(pal = "DkViolet2", dim = 4) 


ggplot() +
  geom_sf(data = sub_catchments, aes(fill = AI_change)) +
  bi_scale_fill(pal = "DkViolet2", dim = 4) 



# Update your legend to match
legend <- bi_legend(pal = "DkViolet2",
                    dim = 4,
                    xlab = "More Precip ",
                    ylab = "More ET ",
                    size = 8,
                    flip_axes = TRUE)

library(cowplot)
br=seq(-.5,0.3,0.1)
limi=c(-0.4,0)
palet=c(hcl.colors(11, palette = "YlOrRd", alpha = NULL, rev = F, fixup = TRUE))
main_plot <- ggplot() +
  # --- Layer 1: Bivariate Climate Background ---
  #geom_sf(data = df_bi, aes( fill = bi_class), show.legend = FALSE,color="transparent") +
  geom_sf(data = sub_catchments, aes( fill = AI_change),color="transparent") +
  # Catchment Boundary
  geom_sf(data = target_catchment, fill = "transparent", color = "gray23", size = 1) +
  
  scale_fill_gradientn(
    colors=palet, breaks=br,limits=limi,trans=scales::modulus_trans(1),
    oob = scales::squish)   +
  
  # bi_scale_fill(pal = "DkViolet2", dim = 4,flip_axes = TRUE) + 
  # scale_fill
  geom_sf(data = databipic_cropped, aes(color = bi_class),size=1,alpha=1,stroke=0,shape=15) +
  
  scale_color_manual(values = colorp, guide="none") +
  # "BlueOr" fits your Red/Blue intent
  theme(axis.title=element_text(size=tsize),
        panel.background = element_rect(fill = "aliceblue", colour = "grey1"),
        panel.border = element_rect(linetype = "solid", fill = NA, colour="black"),
        legend.title = element_text(size=tsize),
        legend.text = element_text(size=osize),
        legend.position = "bottom",
        panel.grid.major = element_line(colour = "grey70"),
        panel.grid.minor = element_line(colour = "grey90"),
        legend.key = element_rect(fill = "transparent", colour = "transparent"),
        legend.key.size = unit(.8, "cm"))+
labs(title = paste("Climatic changes | Rhone catchment"),
     subtitle = "Yearly precipitation, Yearly PET  (1951-1980 vs 1991-2020)")

main_plot
# Combine main_plot and the square legend
final_map=ggarrange(main_plot, legend,
             labels = c("Map", "Key"),
             ncol = 2, nrow = 1,widths = c(2,1), heights=c(1,1), vjust=-1)
# Display result
print(final_map)

# 1. Calculate the bounding box for consistent zooming
bbox <- st_bbox(target_catchment_laea)

# 2. Main Map Construction
main_plot <- ggplot(Europe_laea) +
  # --- Layer 0: Europe Background ---
  geom_sf(fill = "white", color = "grey80") +
  
  # --- Layer 1: Bivariate Climate Polygons (Aggregated) ---
  # geom_sf(data = df_bi, aes(fill = bi_class), color = "transparent", show.legend = FALSE) +
  # bi_scale_fill(pal = "DkViolet2", dim = 4, flip_axes = TRUE) + 
  geom_sf(data = sub_catchments_laea, aes( fill = AI_change),alpha=0.7,color="transparent") +
  # Catchment Boundary
  geom_sf(data = target_catchment_laea, fill = "transparent", color = "gray23", size = 1) +
  
  scale_fill_gradientn(
    colors=palet, breaks=br,limits=limi,trans=scales::modulus_trans(1),
    oob = scales::squish,name = "Δ AI", guide = guide_colorbar(direction = "horizontal", 
    title.position = "top", 
    barwidth = 10 ))   +
  # --- Layer 2: Climate Trend Points (Cropped databipic) ---
  # Using new_scale_color to ensure the points use your specific colorp palette
  new_scale_color() +
  geom_sf(data = databipic_cropped, aes(color = bi_class), 
          size = 1, alpha = 1, stroke = 0, shape = 15) +
  scale_color_manual(values = colorp, guide = "none") +
  
  # --- Layer 3: Catchment Boundary ---
  geom_sf(data = target_catchment_laea, fill = "transparent", color = "gray23", size = 1) +
  # --- Layer: Main Cities ---
  geom_sf(data = cities_laea, color = "black", size = 2) +
  
  
  # --- Layer: City Labels with Transparent Boxes ---
  geom_label_repel(data = cities_coords, 
                   aes(x = x, y = y, label = name),
                   size = 3.5, 
                   fontface = "bold",
                   # Fill: White (#FFFFFF) with ~50% transparency (88)
                   fill = "#FFFFFF88", 
                   # Label border: Set to 0 to remove it, or a small number for a thin line
                   label.size = 0.2,
                   label.padding = unit(0.15, "lines"),
                   box.padding = 0.5,
                   # Ensure the label doesn't cover the city point
                   point.padding = 0.3,
                   segment.color = "black", # Line connecting box to point
                   segment.size = 0.3) +
  
  # --- CRITICAL: Zoom to Catchment but keep Europe background ---
  coord_sf(xlim = c(bbox["xmin"], bbox["xmax"]), 
           ylim = c(bbox["ymin"], bbox["ymax"]), 
           expand = FALSE) +
  
  # --- Styling to match Socioeconomic Map ---
  # --- Styling ---
  theme_minimal() +
  theme(axis.title = element_text(size = tsize),
        panel.background = element_rect(fill = "aliceblue", colour = "grey1"),
        panel.border = element_rect(linetype = "solid", fill = NA, colour = "black"),
        legend.title = element_text(size = tsize, face = "bold"),
        legend.text = element_text(size = osize),
        # This removes the boxes and borders from the legend icons
        legend.key = element_blank(), 
        plot.title = element_text(size = tsize + 2, face = "bold"),
        plot.subtitle = element_text(size = osize),
        # This ensures the legend background itself doesn't have a white box
        legend.background = element_blank(),
        legend.position = "left",
        legend.box = "vertical", # Stacks the multiple legends
        panel.grid.major = element_line(colour = "grey90"))+
  
  
  labs(title = "(a) Climatic changes",
       subtitle = "Aridity Index (1951-1980 vs 1991-2020)",
       x = "Longitude", y = "Latitude")

main_plot
ggsave(paste0("D:/tilloal/Documents/LFRuns_utils/TrendAnalysis/plots/rhone_Cchanges3.jpg"),main_plot, width=24, height=24, units=c("cm"),dpi=1000)


# 3. Create/Update the Legend
# Make sure the legend background is also transparent or matches theme
legend <- bi_legend(pal = "DkViolet2",
                    dim = 4,
                    xlab = "More Precip",
                    ylab = "More ET",
                    size = 10,
                    flip_axes = TRUE)

# 4. Final Combination
final_map <- ggarrange(main_plot, s_plot,
                       labels = c("", ""), # Removed "Map/Key" labels for a cleaner look
                       ncol = 2, nrow = 1, 
                       widths = c(40, 40))

print(final_map)

library(patchwork)
combined_plot <- main_plot + s_plot + 
  plot_layout(guides = "collect") # Optional: combines legends if they are identical

cplot=combined_plot & theme(plot.margin = margin(5, 5, 5, 5))

# Your legend object
legend_biv <- bi_legend(
  pal = colorp, 
  dim = 4,
  xlab = "  +  Drought  -  ",
  ylab = "  -  Flood  +  ",
  size = 10, # Adjusted for a sidebar fit
  arrows = FALSE
) + 
  ggtitle( "Hydrological extremes") + 
  theme(
    plot.title = element_text(
      size = tsize, 
      face = "bold", 
      hjust = 0.5,           # Perfectly centered over the square
      margin = margin(b = 5) # Controls the gap between title and legend
    ),
    legend.title = element_text(size = tsize, face = "bold"),
    plot.background = element_rect(fill = "transparent", colour = NA),
    panel.background = element_rect(fill = "transparent", colour = NA),
    # Ensure the text matches your 'osize' and 'tsize'
    axis.title = element_text(size = osize - 2) 
  )


legend_wrapped <- wrap_elements(panel = legend_biv)

design <- "AAAAABBBBBCC"
combined_plot <- main_plot + s_plot + legend_wrapped + 
  plot_layout(
    design = design, 
    widths = c(10, 10, 10), # Maps are 10x wider than the legend
    guides = "collect"     # Collects standard legends (AI, Land Use, etc.) to the right
  )

# Apply global theme changes
cplot <- combined_plot & 
  theme(
    legend.position = "right",
    legend.box = "vertical",
   # plot.margin = margin(5, 5, 5, 5)
  )

ggsave(
  filename = "D:/tilloal/Documents/LFRuns_utils/TrendAnalysis/plots/Rhone_Combined_Analysis4.jpg",
  plot = cplot,
  width = 40,      # Doubled width for side-by-side
  height = 20,     # Kept height same
  units = "cm",
  dpi = 1000        # 1000 is very high, 600 is usually enough for publication
)

# 1. Combine your two maps and collect the standard legends on the right
# This creates the layout you had, but without the bivariate square
# 1. Base maps with collected guides
base_combined <- main_plot + s_plot + 
  plot_layout(guides = "collect") & 
  theme(legend.position = "right",
        legend.justification = "bottom")

# 2. Inset the legend with "Center-Aligned" coordinates
final_cplot <- base_combined + 
  inset_element(
    legend_wrapped, 
    # Adjust 'left' and 'right' to control width and centering
    # (Increase the gap between them to make the legend bigger)
    left = 0.68,   
    right = 0.95,  
    
    # Adjust 'bottom' and 'top' to control height and vertical position
    bottom = 0.6, 
    top = 1,    
    
    align_to = 'full', 
    on_top = TRUE
  )

# 3. View the result
#final_cplot

ggsave(
  filename = "D:/tilloal/Documents/LFRuns_utils/TrendAnalysis/plots/Rhone_Combined_Analysis5.jpg",
  plot = final_cplot,
  width = 40,      # Doubled width for side-by-side
  height = 20,     # Kept height same
  units = "cm",
  dpi = 1000        # 1000 is very high, 600 is usually enough for publication
)
#ggsave(paste0("D:/tilloal/Documents/LFRuns_utils/TrendAnalysis/plots/rhone_Cchanges.jpg"), width=20, height=20, units=c("cm"),dpi=1000)




#load land use change

#Library calling
suppressWarnings(suppressMessages(library(ncdf4)))
suppressWarnings(suppressMessages(library(sf)))
suppressWarnings(suppressMessages(library(rnaturalearth)))
suppressWarnings(suppressMessages(library(rnaturalearthdata)))
suppressWarnings(suppressMessages(library(rgeos)))
suppressWarnings(suppressMessages(library(dplyr)))
library(Kendall)
library(biscale)
library(cowplot)
library(ggpubr)
library(ggridges)
library(ggplot2)
library(viridis)
library(hrbrthemes)
library(tidyverse)
library(raster)
library(modifiedmk)
library(ks)
library(pracma)
library(data.table)



#Functions
outletopen=function(dir,outletname,nrspace=rep(NA,5)){
  ncbassin=paste0(dir,"/",outletname,".nc")
  ncb=nc_open(ncbassin)
  name.vb=names(ncb[['var']])
  namev=name.vb[1]
  if ("Band1"%in% name.vb)namev="Band1"
  name.lon="lon"
  name.lat="lat"
  if (!is.na(nrspace[1])){
    start=as.numeric(nrspace[c(2,4)])
    count=as.numeric(nrspace[c(3,5)])-start+1
  }else{
    londat = ncvar_get(ncb,name.lon) 
    llo=length(londat)
    latdat = ncvar_get(ncb,name.lat)
    lla=length(latdat)
    start=c(1,1)
    count=c(llo,lla)
  }
  
  londat = ncvar_get(ncb,name.lon,start=start[1],count=count[1]) 
  llo=length(londat)
  latdat = ncvar_get(ncb,name.lat,start=start[2],count=count[2])
  lla=length(latdat)
  outlets = ncvar_get(ncb,namev,start = start, count= count) 
  outlets=as.vector(outlets)
  outll=expand.grid(londat,latdat)
  lonlatloop=expand.grid(c(1:llo),c(1:lla))
  outll$idlo=lonlatloop$Var1
  outll$idla=lonlatloop$Var2
  
  outll=outll[which(!is.na(outlets)),]
  outlets=outlets[which(!is.na(outlets))]
  outll=data.frame(outlets,outll)
  return (outll)
}

UpAopen=function(dir,outletname,Sloc_final){
  ncbassin=paste0(dir,outletname)
  ncb=nc_open(ncbassin)
  name.vb=names(ncb[['var']])
  namev=name.vb[2]
  #time <- ncvar_get(ncb,"time")
  
  #timestamp corretion
  name.lon="lon"
  name.lat="lat"
  londat = ncvar_get(ncb,name.lon) 
  llo=length(londat)
  latdat = ncvar_get(ncb,name.lat)
  lla=length(latdat)
  start=c(1,1)
  count=c(llo,lla)
  
  
  londat = ncvar_get(ncb,name.lon,start=start[1],count=count[1]) 
  llo=length(londat)
  latdat = ncvar_get(ncb,name.lat,start=start[2],count=count[2])
  lla=length(latdat)
  outlets = ncvar_get(ncb,namev,start = start, count= count) 
  outlets=as.vector(outlets)/1000000
  outll=expand.grid(londat,latdat)
  lonlatloop=expand.grid(c(1:llo),c(1:lla))
  outll$upa=outlets
  outll$idlo=lonlatloop$Var1
  outll$idla=lonlatloop$Var2
  
  #outll$idlalo=paste(outll$idlo,outll$idla,sep=" ")
  outll$latlong=paste(round(outll$Var1,4),round(outll$Var2,4),sep=" ")
  outfinal=inner_join(outll, Sloc_final, by="latlong")
  return (outfinal)
}

ReservoirOpen=function(dir,outletname,Sloc_final){
  ncbassin=paste0(dir,outletname)
  ncb=nc_open(ncbassin)
  name.vb=names(ncb[['var']])
  namev=name.vb[1]
  #time <- ncvar_get(ncb,"time")
  
  #timestamp corretion
  name.lon="lon"
  name.lat="lat"
  londat = ncvar_get(ncb,name.lon) 
  llo=length(londat)
  latdat = ncvar_get(ncb,name.lat)
  lla=length(latdat)
  start=c(1,1)
  count=c(llo,lla)
  
  
  londat = ncvar_get(ncb,name.lon,start=start[1],count=count[1]) 
  llo=length(londat)
  latdat = ncvar_get(ncb,name.lat,start=start[2],count=count[2])
  lla=length(latdat)
  outlets = ncvar_get(ncb,namev,start = start, count= count) 
  outlets=as.vector(outlets)/1000000
  outll=expand.grid(londat,latdat)
  lonlatloop=expand.grid(c(1:llo),c(1:lla))
  outll$upa=outlets
  outll$idlo=lonlatloop$Var1
  outll$idla=lonlatloop$Var2
  
  #outll$idlalo=paste(outll$idlo,outll$idla,sep=" ")
  outll$latlong=paste(round(outll$Var1,4),round(outll$Var2,4),sep=" ")
  outfinal=inner_join(outll, Sloc_final, by="latlong")
  return (outfinal)
}



#Set data directory
hydroDir<-("D:/tilloal/Documents/LFRuns_utils/data")

#load outf

outf=c()
for( Nsq in 1:88){
  print(Nsq)
  rspace= read.csv(paste0(hydroDir,"/subspace_efas.csv"))
  rspace=rspace[,-1]
  nrspace=rspace[Nsq,]
  #outletname="outletsv8_hybas07_01min"
  #outletname="outlets_hybas09_01min"
  outletname="efas_rnet_100km_01min"
  
  outhybas=outletopen(hydroDir,outletname,nrspace)
  Idstart=as.numeric(Nsq)*10000
  Idstart2=as.numeric(Nsq)*100000
  if (length(outhybas$outlets)>0){
    outhybas$outlets=seq((Idstart+1),(Idstart+length(outhybas$outlets)))
    outhybas$outl2=seq((Idstart2+1),(Idstart2+length(outhybas$outlets)))
    outhybas$latlong=paste(round(outhybas$Var1,4),round(outhybas$Var2,4),sep=" ")
    #outcut=which(!is.na(match(outhybas$outlets,parlist$catchment)))
    # zebi=seq(parlist$catchment[1],parlist$catchment[length(parlist$catchment)])
    # outcut=which(!is.na(match(outhybas$outlets,zebi)))
    outhloc=outhybas
    outf=rbind(outf,outhloc)
  }
}

#Load my shapefile on which to aggregate

### Hybas07 ----
Catchmentrivers7=read.csv(paste0(hydroDir,"/Catchments/from_hybas_eu_onlyid.csv"),encoding = "UTF-8", header = T, stringsAsFactors = F)
hybas07 <- read_sf(dsn = paste0(hydroDir,"/Catchments/hydrosheds/hybas_eu_lev05_v1c.shp"))
hybasf7=fortify(hybas07)
Catamere07=inner_join(hybasf7,Catchmentrivers7,by= "HYBAS_ID")
Catamere07$llcoord=paste(round(Catamere07$POINT_X,4),round(Catamere07$POINT_Y,4),sep=" ")
cst7=right_join(Catamere07,outf,by= c("llcoord"="latlong"))
GNF=cst7
st_geometry(GNF)=NULL

UnHY=unique(GNF$HYBAS_ID)
# 
#cst7=st_transform(cst7,  crs=3035)



### HydroRegions ----

GridHR=raster( paste0(hydroDir,"/HydroRegions_raster_WGS84.tif"))
GHR=as.data.frame(GridHR,xy=T)
GHR=GHR[which(!is.na(GHR[,3])),]
GHR$llcoord=paste(round(GHR$x,4),round(GHR$y,4),sep=" ")
GHR_riv=inner_join(GHR,outf,by= c("llcoord"="latlong"))
GHshpp <- read_sf(dsn ="Z:/ClimateRun4/nahaUsers/tilloal/HydroRegions/her_all_adjusted.shp")
HydroRsf=fortify(GHshpp)

UnHY=unique(GHR_riv$HydroRegions_raster_WGS84)



### NUTS3 ----


# NUTS3 <- read_sf(dsn = paste0(hydroDir,"/Countries/NUTS3/NUTS3_Extended_domain.shp"))
# NUTS3$N3ID=c(1:length(NUTS3$NUTS_ID))
# N2ID=unique(NUTS3$NUTS2_ID)
# N2IDn=c(1:length(N2ID))
# mati=match(NUTS3$NUTS2_ID,N2ID)
# NUTS3$N2ID=N2IDn[mati]
# st_write(NUTS3, paste0(hydroDir,"/Countries/NUTS3/NUTS3_modified.shp"), driver = "ESRI Shapefile")

NUTS3 <- read_sf(dsn = paste0(hydroDir,"/Countries/NUTS3/NUTS3_modified.shp"))
GridNUTS3=raster( paste0(hydroDir,"/Countries/NUTS3/NUTS3_Raster3ID.tif"))
GN3=as.data.frame(GridNUTS3,xy=T)
GN3=GN3[which(!is.na(GN3[,3])),]
GN3$llcoord=paste(round(GN3$x,4),round(GN3$y,4),sep=" ")
GN3_riv=right_join(GN3,outf,by= c("llcoord"="latlong"))

GridNUTS2=raster( paste0(hydroDir,"/Countries/NUTS3/NUTS3_Raster2ID.tif"))
GN2=as.data.frame(GridNUTS2,xy=T)
GN2=GN2[which(!is.na(GN2[,3])),]
GN2$llcoord=paste(round(GN2$x,4),round(GN2$y,4),sep=" ")
GN2_riv=right_join(GN2,outf,by= c("llcoord"="latlong"))

GNF=right_join(GN3,GN2_riv,by="llcoord")

GNUTS3sf=fortify(NUTS3)

GNFx=GNF[which(is.na(GNF$NUTS3_Raster3ID)),]


#Catchment aggregation
Catchments <- read_sf(dsn = paste0("D:/tilloal/Documents/01_Projects/RegimeShifts/OtherCatchments_polygon.shp"))

#I choose a specific catchment, all inputs will be cropped onto that catchment

selectedHybas="2050016510"

#load data to be aggregated



### Plot parameters ----
palet2=c(hcl.colors(9, palette = "Blues", alpha = NULL, rev = TRUE, fixup = TRUE))
# Plot of ordered change by region, can be important
palet=c(hcl.colors(11, palette = "RdYlBu", alpha = NULL, rev = F, fixup = TRUE))
outletname="efas_rnet_100km_01min"
outll=outletopen(hydroDir,outletname)
cord.dec=outll[,c(2,3)]
cord.dec = SpatialPoints(cord.dec, proj4string=CRS("+proj=longlat"))
cord.UTM <- spTransform(cord.dec, CRS("+init=epsg:3035"))
nco=cord.UTM@coords
world <- ne_countries(scale = "medium", returnclass = "sf")
Europe <- world[which(world$continent == "Europe"),]
e2=st_transform(Europe,  crs=3035)
w2=st_transform(world,  crs=3035)
tsize=12
osize=12
Impdates=seq(1950,2020,by=10)
valuenames=paste0("Y",Impdates)
catmap=cst7
basemap=w2



#Climate variables

#Load every yearly file, compare mean of first 30 years and last 30 years
#use cdo instead and just load the two files

yrlist1=c(1951:1980)
yrlist2=c(1991:2020)
for (yr in yrlist1)
{
  
}

#Land use
library(exactextractr)
rastforest=raster( "D:/tilloal/Documents/06_Floodrivers/landuse/fracforest_ch20201951.tif")

rastsealed=raster( "D:/tilloal/Documents/06_Floodrivers/landuse/fracsealed_ch20201951.tif")

rastirrigated=raster( "D:/tilloal/Documents/06_Floodrivers/landuse/fracirrigated_ch20201951.tif")

rastother=raster( "D:/tilloal/Documents/06_Floodrivers/landuse/fracother_ch20201951.tif")

rastrice=raster( "D:/tilloal/Documents/06_Floodrivers/landuse/fracrice_ch20201951.tif")

rastwater=raster( "D:/tilloal/Documents/06_Floodrivers/landuse/fracwater_ch20201951.tif")

#I could do a relative change as well

forestchange<- exact_extract(rastforest, GHshpp, 'mean')
GHshpp$forestchange=forestchange
GHshpp$sealedchange <- exact_extract(rastsealed, GHshpp, 'mean')
GHshpp$irrigatedchange <- exact_extract(rastirrigated, GHshpp, 'mean')
GHshpp$otherchange <- exact_extract(rastother, GHshpp, 'mean')
GHshpp$ricechange <- exact_extract(rastrice, GHshpp, 'mean')
GHshpp$waterchange <- exact_extract(rastwater, GHshpp, 'mean')

mhh=match(UnHY,GHshpp$Id)
GHshppH=GHshpp[mhh,]
df_GHshppH=data.frame(GHshppH)
st_geometry(df_GHshppH)<-NULL


lucmap=list()
luclass=c("forest","sealed","irrigated","other","rice","water")

for (li in 1:length(luclass))
{
  lu=luclass[li]
  print(lu)
  GHshppH$fill=as.numeric(df_GHshppH[,7+li])
  
  flims=(quantile(GHshppH$fill,c(0.01,0.99),na.rm=T))
  lims=c(-round(max(abs(flims)),1),round(max(abs(flims)),1))
  if (diff(lims)==0){
    lims=c(-round(max(abs(flims)),2),round(max(abs(flims)),2))
  }
  if (diff(lims)==0){
    lims=c(-round(max(abs(flims)),4),round(max(abs(flims)),4))
  }
  #lims=c(-25,25)
  lucmap<-ggplot(basemap) +
    geom_sf(fill="gray95",color="transparent",size=0.5)+
    geom_sf(data=GHshppH,aes(fill=fill*100,geometry=geometry),alpha=1,color="transparent")+
    geom_sf(fill="transparent",color="gray30",size=0.5)+
    coord_sf(xlim = c(min(nco[,1]),max(nco[,1])), ylim = c(min(nco[,2]),max(nco[,2])))+
    scale_fill_gradientn(
      colors=palet,
      limits=lims*100,oob = scales::squish,
      name=paste0("Change in ",lu,"(%)"))   +
    labs(x="Longitude", y = "Latitude")+
    theme(axis.title=element_text(size=tsize),
          title = element_text(size=16),
          axis.text=element_text(size=osize),
          panel.background = element_rect(fill = "aliceblue", colour = "grey1"),
          panel.border = element_rect(linetype = "solid", fill = NA, colour="black"),
          legend.title = element_text(size=tsize),
          legend.text = element_text(size=osize),
          legend.position = "right",
          panel.grid.major = element_line(colour = "grey70"),
          panel.grid.minor = element_line(colour = "grey90"),
          legend.key = element_rect(fill = "transparent", colour = "transparent"),
          legend.key.size = unit(1, "cm"))+
    ggtitle(lu)
  ggsave(paste0("D:/tilloal/Documents/LFRuns_utils/TAplots/luchange_",lu,"_20201951_HR.jpg"), lucmap,width=30, height=20, units=c("cm"),dpi=300)
  
}



#Same for water demand change



#divide by area
workDir<-("D:/tilloal/Documents/06_Floodrivers/")
ncpix= paste0(workDir,"/mapscal/pixarea_European_01min.nc")
ncp=nc_open(ncpix)
t=ncp$var[[2]]
tsize<-t$varsize
pixarea=ncvar_get(ncp,names(ncp[['var']])[2]) 
plon=ncvar_get(ncp,"lon") 
plat=ncvar_get(ncp,"lat") 
pixarea=as.matrix(t(pixarea))

rast_totwd<-raster( "D:/tilloal/Documents/06_Floodrivers/wateruse/wateruse_sums/all_demands_2020.tif")
rast_totwd<-raster( "D:/tilloal/Documents/06_Floodrivers/wateruse/wateruse_sums/all_demands_1951.tif")
#Convert RasterLayer to matrix
rast_tmat <- (as.matrix(rast_totwd))
#multiply by 30.4 to have the real sum values
rast_tmat=rast_tmat*30.4
#m3/m2 to m3
rast_tmat=rast_tmat*pixarea
#convert mm to m3
rast_tmat=rast_tmat/1e3
#m3 to km3
rast_tmat=rast_tmat/1e9
rast_tmout <- raster(nrows=nrow(rast_totwd), ncols=ncol(rast_totwd), ext=extent(rast_totwd))
crs(rast_tmout) <- crs(rast_tmat)
values(rast_tmout) <- rast_tmat
rast_totwd=rast_tmout

HRwgs84 <- st_transform(GHshpp, crs = 4326)
HRwgs84$totalwd2020 <- exact_extract(rast_totwd, HRwgs84, 'sum')
sum(HRwgs84$totalwd2020)
HRwgs84$wdpkm2=HRwgs84$totalwd2020/HRwgs84$SURF_KM2*1000*1000




mhh=match(UnHY,GHshpp$HYBAS_ID)
hybas07$totalwd2020 <- exact_extract(rast_totwd, hybas07, 'sum')
hybas07H=hybas07[mhh,]
sum(hybas07H$totalwd2020,na.rm=T)
hybas07H$wdpkm2=hybas07H$totalwd2020/hybas07H$UP_AREA*1000*1000

#unit is now in mm
wd="total water demand in 1951"
tsize=12
flims=(quantile(HRwgs84$wdpkm2,c(0.1,0.95),na.rm=T))
lims=c(0,round(max(abs(flims)),1))
lims=c(0,120)
wdmap<-ggplot(basemap) +
  geom_sf(fill="gray95",color="transparent",size=0.5)+
  geom_sf(data=HRwgs84,aes(fill=wdpkm2,geometry=geometry),alpha=1,color="transparent")+
  geom_sf(fill="transparent",color="gray30",size=0.5)+
  coord_sf(xlim = c(min(nco[,1]),max(nco[,1])), ylim = c(min(nco[,2]),max(nco[,2])))+
  scale_fill_gradientn(
    colors=palet2,
    limits=lims,oob = scales::squish, trans="sqrt", breaks=c(1,5,20,50,100,200),
    name=paste0("(mm/year)"))   +
  labs(x="Longitude", y = "Latitude")+
  theme(axis.title=element_text(size=tsize),
        title = element_text(size=16),
        axis.text=element_text(size=osize),
        panel.background = element_rect(fill = "aliceblue", colour = "grey1"),
        panel.border = element_rect(linetype = "solid", fill = NA, colour="black"),
        legend.title = element_text(size=tsize),
        legend.text = element_text(size=osize),
        legend.position = "right",
        panel.grid.major = element_line(colour = "grey70"),
        panel.grid.minor = element_line(colour = "grey90"),
        legend.key = element_rect(fill = "transparent", colour = "transparent"),
        legend.key.size = unit(1, "cm"))+
  ggtitle(wd)

wdmap
ggsave(paste0("D:/tilloal/Documents/LFRuns_utils/TAplots/TotalWaterDemand_1951.jpg"), wdmap,width=30, height=20, units=c("cm"),dpi=300)



rast_ene=raster( "D:/tilloal/Documents/06_Floodrivers/wateruse/ene_ysum_ch20201951.tif")
#Convert RasterLayer to matrix
rast_enemat <- (as.matrix(rast_ene))
#multiply by 30.4 to have the real sum values
rast_enemat=rast_enemat*30.4
#m3/m2 to m3
rast_enemat=rast_enemat*pixarea
#convert mm to m3
rast_enemat=rast_enemat/1e3
#m3 to km3
rast_enemat=rast_enemat/1e9
rast_eneout <- raster(nrows=nrow(rast_ene), ncols=ncol(rast_ene), ext=extent(rast_ene))
crs(rast_eneout) <- crs(rast_ene)
values(rast_eneout) <- rast_enemat
rast_ene=rast_eneout

rast_dom=raster( "D:/tilloal/Documents/06_Floodrivers/wateruse/dom_ysum_ch20201951.tif")
#Convert RasterLayer to matrix
rast_dommat <- (as.matrix(rast_dom))
#multiply by 30.4 to have the real sum values
rast_dommat=rast_dommat*30.4
#m3/m2 to m3
rast_dommat=rast_dommat*pixarea
#convert mm to m3
rast_dommat=rast_dommat/1e3
#m3 to km3
rast_dommat=rast_dommat/1e9
rast_domout <- raster(nrows=nrow(rast_dom), ncols=ncol(rast_dom), ext=extent(rast_dom))
crs(rast_domout) <- crs(rast_dom)
values(rast_domout) <- rast_dommat
rast_dom=rast_domout
rast_liv=raster( "D:/tilloal/Documents/06_Floodrivers/wateruse/liv_ysum_ch20201951.tif")
#Convert RasterLayer to matrix
rast_livmat <- (as.matrix(rast_liv))
#multiply by 30.4 to have the real sum values
rast_livmat=rast_livmat*30.4
#m3/m2 to m3
rast_livmat=rast_livmat*pixarea
#convert mm to m3
rast_livmat=rast_livmat/1e3
#m3 to km3
rast_livmat=rast_livmat/1e9
rast_livout <- raster(nrows=nrow(rast_liv), ncols=ncol(rast_liv), ext=extent(rast_liv))
crs(rast_livout) <- crs(rast_liv)
values(rast_livout) <- rast_livmat
rast_liv=rast_livout

rast_ind=raster( "D:/tilloal/Documents/06_Floodrivers/wateruse/ind_ysum_ch20201951.tif")
#Convert RasterLayer to matrix
rast_indmat <- (as.matrix(rast_ind))
#multiply by 30.4 to have the real sum values
rast_indmat=rast_indmat*30.4
#m3/m2 to m3
rast_indmat=rast_indmat*pixarea
#convert mm to m3
rast_indmat=rast_indmat/1e3
#m3 to km3
rast_indmat=rast_indmat/1e9
rast_indout <- raster(nrows=nrow(rast_ind), ncols=ncol(rast_ind), ext=extent(rast_ind))
crs(rast_indout) <- crs(rast_ind)
values(rast_indout) <- rast_indmat
rast_ind=rast_indout

rast_total=raster( "D:/tilloal/Documents/06_Floodrivers/wateruse/all_ysum_ch20201951.tif")
#Convert RasterLayer to matrix
rast_totalmat <- (as.matrix(rast_total))
#multiply by 30.4 to have the real sum values
rast_totalmat=rast_totalmat*30.4
#m3/m2 to m3
rast_totalmat=rast_totalmat*pixarea
#convert mm to m3
rast_totalmat=rast_totalmat/1e3
#m3 to km3
rast_totalmat=rast_totalmat/1e9
rast_totalout <- raster(nrows=nrow(rast_total), ncols=ncol(rast_total), ext=extent(rast_total))
crs(rast_totalout) <- crs(rast_total)
values(rast_totalout) <- rast_totalmat
rast_total=rast_totalout


HRwgs84$enechange <- exact_extract(rast_ene, HRwgs84, 'sum')
HRwgs84$domchange <- exact_extract(rast_dom, HRwgs84, 'sum')
HRwgs84$livchange <- exact_extract(rast_liv, HRwgs84, 'sum')
HRwgs84$indchange <- exact_extract(rast_ind, HRwgs84, 'sum')
HRwgs84$totalchange <- exact_extract(rast_total, HRwgs84, 'sum')


mhh=match(UnHY,HRwgs84$Id)
HRwgs84h=HRwgs84[mhh,]
df_HRwgs84h=data.frame(HRwgs84h)
st_geometry(df_HRwgs84h)<-NULL

# mhh=match(UnHY,hybas07$HYBAS_ID)
# hybas07H=hybas07[mhh,]
# df_hybas07H=data.frame(hybas07H)
# st_geometry(df_hybas07H)<-NULL

wdclass=c("ene","dom","liv","ind","total")
wdmap=list()
for (li in 1:length(wdclass))
{
  wd=wdclass[li]
  print(wd)
  HRwgs84h$fill=as.numeric(df_HRwgs84h[,15+li])
  
  flims=(quantile(HRwgs84h$fill,c(0.01,0.99),na.rm=T))
  lims=c(-round(max(abs(flims)),1),round(max(abs(flims)),1))
  if (diff(lims)==0){
    lims=c(-round(max(abs(flims)),2),round(max(abs(flims)),2))
  }
  if (diff(lims)==0){
    lims=c(-round(max(abs(flims)),4),round(max(abs(flims)),4))
  }
  tsize=12
  wdmap[[li]]<-ggplot(basemap) +
    geom_sf(fill="gray95",color="transparent",size=0.5)+
    geom_sf(data=HRwgs84h,aes(fill=fill,geometry=geometry),alpha=1,color="transparent")+
    geom_sf(fill="transparent",color="gray30",size=0.5)+
    coord_sf(xlim = c(min(nco[,1]),max(nco[,1])), ylim = c(min(nco[,2]),max(nco[,2])))+
    scale_fill_gradientn(
      colors=palet,
      limits=lims,oob = scales::squish,
      name=paste0("Change (km3/year)"))   +
    labs(x="Longitude", y = "Latitude")+
    theme(axis.title=element_text(size=tsize),
          title = element_text(size=16),
          axis.text=element_text(size=osize),
          panel.background = element_rect(fill = "aliceblue", colour = "grey1"),
          panel.border = element_rect(linetype = "solid", fill = NA, colour="black"),
          legend.title = element_text(size=tsize),
          legend.text = element_text(size=osize),
          legend.position = "right",
          panel.grid.major = element_line(colour = "grey70"),
          panel.grid.minor = element_line(colour = "grey90"),
          legend.key = element_rect(fill = "transparent", colour = "transparent"),
          legend.key.size = unit(1, "cm"))+
    ggtitle(wd)
  ggsave(paste0("D:/tilloal/Documents/LFRuns_utils/TAplots/wdchangeHR_",wd,"_20201951.jpg"),wdmap[[li]],width=30, height=20, units=c("cm"),dpi=300)
  
}

wdmap[[1]]

###################################################################
#detailed analysis of water demand (seasonal)

#load netcdf of water demand

print(year)
df_wd=c()
for (w in wdclass[-5]){
  
  nc <- nc_open(paste0("D:/tilloal/Documents/06_Floodrivers/wateruse/wateruse_histo/",w,"_1950_2020.nc"))
  print(w)
  # Extract the dimensions of the file
  lon_dim <- nc$dim[["lon"]]
  lat_dim <- nc$dim[["lat"]]
  time_dim <- nc$dim[["time"]]
  name.vb=names(nc[['var']])
  namev=name.vb[1]
  time <- ncvar_get(nc,"time")
  timestamp=as.Date(time,origin="1950-01-01")
  lt=length(time)
  
  
  # Extract the longitude and latitude values
  lon <- ncvar_get(nc, "lon")
  lat <- ncvar_get(nc, "lat")
  
  outll=expand.grid(lon,lat)
  mvs=c()
  for (d in 1:lt){
    print(d/lt*100)
    mon_val <- ncvar_get(nc, namev, 
                         start = c(1, 1, d), 
                         count = c(lon_dim$len, lat_dim$len, 1))
    md=sum(mon_val,na.rm=T)
    
    mvs=c(mvs,md)
  }
  
  df_wd=rbind(df_wd,mvs)
}

df_wd=as.data.frame(t(df_wd))
names(df_wd)=wdclass[-5]
df_wd$sum=df_wd$ind+df_wd$liv+df_wd$dom+df_wd$ene
df_wd$time=timestamp

save(df_wd,file=paste0("D:/tilloal/Documents/06_Floodrivers/wateruse/wd_domain.Rdata"))


# Create a data frame from the first day
dx <- data.frame(lon = as.vector(outll$Var1), 
                 lat = as.vector(outll$Var2), 
                 value = as.vector(first_day))

####################################################################
#load results from drought catchment analysis 

#floods
load(file=paste0(hydroDir,"/TSEVA/output_plots/outputs_flood_year_qsp.Rdata"))
#droughts
load(file=paste0(hydroDir,"/TSEVA/output_plots/outputs_drought_nonfrost_qsp.Rdata"))

#I extract the trend at MUTS3 level fist

FloodTrends=Output_fl_year$TrendOutlets
# mhh=match(UnHY,FloodTrends$HydroR)
FloodTrends=FloodTrends[which(FloodTrends$driver=="Landuse"),]
length(which(!is.na(FloodTrends$Y2020)))
Output_dr_nfrost=Output_dr_year
DroughtTrends=Output_dr_nfrost$TrendOutlets
DroughtTrends=DroughtTrends[which(DroughtTrends$driver=="Landuse"),]


#join flood and land use changes
#Flood_xplain=inner_join(FloodTrends,df_hybas07H,by=c("HydroR"="HYBAS_ID"))
Flood_xplain=inner_join(FloodTrends,df_hybas07H,by=c("HYBAS_ID"))

plot(Flood_xplain$Rchange.Y2020,Flood_xplain$sealedchange)
corsealed=cor(Flood_xplain$Rchange.Y2020,Flood_xplain$sealedchange)

plot(Flood_xplain$Rchange.Y2020,Flood_xplain$forestchange)
corforest=cor(Flood_xplain$Rchange.Y2020,Flood_xplain$forestchange)

plot(Flood_xplain$Rchange.Y2020,Flood_xplain$sealedchange)
corsealed=cor(Flood_xplain$Rchange.Y2020,Flood_xplain$sealedchange)

plot(Flood_xplain$Rchange.Y2020,Flood_xplain$sealedchange)
corsealed=cor(Flood_xplain$Rchange.Y2020,Flood_xplain$sealedchange)

Flood_xplainLU=Flood_xplain[,c(71,86,87,88,89,90,91)]
Flood_xplainLU=Flood_xplain[,c(95,112,113,114,115,116,117)]


library(tidyverse)
library(ggpubr)



# Reshape the data from wide to long format
long_data <- Flood_xplainLU %>%
  gather(key = "variable", value = "value", -Y2020)
long_data=long_data[-which(abs(long_data$value)<=1e-3),]
# Create the base plot with facets
p <- long_data %>%
  ggplot(aes(x = value*100, y = Y2020)) +
  facet_wrap(~ variable, scales = "free") +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +  # Add a linear regression line
  # Calculate the Pearson correlation for each facet and add it as text
  geom_text(data = long_data %>% 
              group_by(variable) %>%
              summarize(correlation = cor(value*100, Y2020, use = "complete.obs")) %>%
              mutate(label = paste0("r = ", round(correlation, 2))),
            aes(label = label, x = Inf, y = Inf),
            hjust = 1.1, vjust = 1.1, inherit.aes = FALSE, check_overlap = TRUE)+
  labs(title = "change in 10y flood (l/s/km2) attributed to LUC vs changes in land use fractions (%)")

# Print the final plot
print(p)
ggsave(paste0("D:/tilloal/Documents/LFRuns_utils/TAplots/LUdriversFLOOD_20201951_outlets.jpg"), p,width=30, height=20, units=c("cm"),dpi=300)


#Same for drought


#join flood and land use changes
Drought_xplain=inner_join(DroughtTrends,df_hybas07H,by=c("HYBAS_ID"))

#Drought_xplainLU=Drought_xplain[,c(71,86,87,88,89,90,91)]
Drought_xplainLU=Drought_xplain[,c(95,112,113,114,115,116,117)]

library(tidyverse)
library(ggpubr)



# Reshape the data from wide to long format
long_data <- Drought_xplainLU %>%
  gather(key = "variable", value = "value", -Y2020)
long_data=long_data[-which(abs(long_data$value)<=1e-3),]
# Create the base plot with facets
p <- long_data %>%
  ggplot(aes(x = value*100, y = Y2020)) +
  facet_wrap(~ variable, scales = "free") +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +  # Add a linear regression line
  # Calculate the Pearson correlation for each facet and add it as text
  geom_text(data = long_data %>% 
              group_by(variable) %>%
              summarize(correlation = cor(value*100, Y2020, use = "complete.obs")) %>%
              mutate(label = paste0("r = ", round(correlation, 2))),
            aes(label = label, x = Inf, y = Inf),
            hjust = 1.1, vjust = 1.1, inherit.aes = FALSE, check_overlap = TRUE)+
  labs(title = "change in 10y drought (l/s/km2) attributed to LUC vs changes in land use fractions (%)")

# Print the final plot
print(p)
ggsave(paste0("D:/tilloal/Documents/LFRuns_utils/TAplots/LUdriversDROUGHT_20201951_outlets.jpg"), p,width=30, height=20, units=c("cm"),dpi=300)





#Water demand
DroughtTrends=Output_dr_nfrost$TrendOutlets
DroughtTrends=DroughtTrends[which(DroughtTrends$driver=="Wateruse"),]

Drought_xplain=inner_join(DroughtTrends,df_hybas07H,by=c("HYBAS_ID"))

Drought_xplainWD=Drought_xplain[,c(71,92,93,94,95,96)]
Drought_xplainWD=Drought_xplain[,c(95,118,119,120,121,122)]



# Reshape the data from wide to long format
long_data <- Drought_xplainWD %>%
  gather(key = "variable", value = "value", -Y2020)
long_data=long_data[-which(abs(long_data$value)<=1e-3),]
# Create the base plot with facets
p <- long_data %>%
  ggplot(aes(x = value, y = Y2020)) +
  facet_wrap(~ variable, scales = "free") +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +  # Add a linear regression line
  # Calculate the Pearson correlation for each facet and add it as text
  geom_text(data = long_data %>% 
              group_by(variable) %>%
              summarize(correlation = cor(value, Y2020, use = "complete.obs")) %>%
              mutate(label = paste0("R = ", round(correlation, 2))),
            aes(label = label, x = Inf, y = Inf),
            hjust = 1.1, vjust = 1.1, inherit.aes = FALSE, check_overlap = TRUE)+
  labs(title = "change in 10y drought (% of 1951 RL) attributed to water demand vs changes in water demand (km3)")

# Print the final plot
print(p)
#ggsave(paste0("D:/tilloal/Documents/LFRuns_utils/TAplots/WDdriversDROUGHT_20201951_outlets.jpg"), p,width=30, height=20, units=c("cm"),dpi=300)


#Reservoirs



#load water demand change/ Aggregate at NUTS3

uppressWarnings(suppressMessages(library(ncdf4)))
suppressWarnings(suppressMessages(library(sf)))
suppressWarnings(suppressMessages(library(rnaturalearth)))
suppressWarnings(suppressMessages(library(rnaturalearthdata)))
suppressWarnings(suppressMessages(library(rgeos)))
suppressWarnings(suppressMessages(library(dplyr)))
library(ggplot2)
library(raster)

hydroDir<-("D:/tilloal/Documents/LFRuns_utils/data")
workDir<-("D:/tilloal/Documents/06_Floodrivers/")

#continue this part
ncpix= paste0(workDir,"/mapscal/pixarea_European_01min.nc")
ncp=nc_open(ncpix)
t=ncp$var[[2]]
tsize<-t$varsize
pixarea=ncvar_get(ncp,names(ncp[['var']])[2]) 


Area_domain=raster( paste0(workDir,"/DataPaper/Revisions/GIS/area_HERA_Domain.tif"))
Area_domain=as.data.frame(Area_domain,xy=T)
Area_domain=Area_domain[which(Area_domain$area_HERA_domain==1),]
Area_domain$latlong=paste(round(Area_domain$x,3),round(Area_domain$y,3),sep=" ")

#France=countries[which(countries$NAME_FREN=="France"),]


Wd2="D:/tilloal/Documents/LFRuns_utils"

seclist=c("dom","ene","ind","liv")
swd_vr=c()
for (sector in seclist){
  print(sector)
  ncdis= paste0(workDir,"/wateruse/wateruse_sums/",sector,"_ysum.nc")
  ncd=nc_open(ncdis)
  
  nav=names(ncd[['var']])
  #Band1 is the second variable
  t=ncd$var[[2]]
  name.var=names(ncd$var)[3]
  tsize<-t$varsize
  tdims<-t$ndims
  nt1<-tsize[tdims]
  
  name.lon="lon"
  name.lat="lat"
  
  lon=ncvar_get(ncd,name.lon)
  lat=ncvar_get(ncd,name.lat)
  lon=lon+0.01666667
  lat=lat+0.01666667
  efas_lat=c(34.50,72.25)
  efas_lon=c(-25.25, 35.00)
  # lolon=c(which(round(lon,2)==efas_lon[1]),which(round(lon,2)==efas_lon[2]))
  # lolat=c(which(round(lat,2)==efas_lat[1]),which(round(lat,2)==efas_lat[2]))
  # 
  # lon[2]-lon[1]
  ll=expand.grid(lon,lat)
  ll=as.data.frame(ll)
  min(ll$Var2)
  min(Area_domain$y)
  ll$latlong=paste(round(ll$Var1,3),round(ll$Var2,3),sep=" ")
  
  matcho=match(Area_domain$latlong,ll$latlong)
  
  start <- rep(1,tdims) # begin with start=(1,1,1,...,1)
  count <- tsize # begin w/count=(nx,ny,nz,...,nt), reads entire var
  count[3]=1
  yearlist=c(1951:2020)
  swd_yr=c()
  for(id in 2:71){
    print(id)
    start[3]=id
    yr=yearlist[id]
    #Here I need to extract only some values as this is the full EFAS domain
    w_demand   = ncvar_get(ncd,name.var,start = start, count= count) 
    #multiply by 30 to have total year sum (1 value per month in initial dataset)
    mult=365.25/12
    w_demand=w_demand*mult
    
    #now multiply by area in m2
    # w_d_w=w_demand[matcho]
    # wk=sum(w_d_w)/length(w_d_w)
    w_demand_a=w_demand*pixarea
    pixarea[1]
    #set the mask for only HERA domain
    w_demand_dom=w_demand_a[matcho]
    #w_demand_dom=w_demand_a
    
    #values are in mm/yr
    #convert to m3/yr
    w_demand_dom=w_demand_dom/1e3
    swdv=sum(w_demand_dom,na.rm=T)
    
    swd1=swdv
    swd_yr=c(swd_yr,swd1)
    
  }
  
  swd_vr=cbind(swd_vr,swd_yr)
  
}


#do a plot with all water demand 



#transform m3 to Km3
swd_vrkm3=swd_vr/1e9
swd_vrkm3=as.data.frame(swd_vrkm3)
swd_vrkm3$total=swd_vrkm3[,1]+swd_vrkm3[,2]+swd_vrkm3[,3]+swd_vrkm3[,4]
plot(sum_vrkm3)
plot(swd_vrkm3[,1],type="l",ylim=c(0,1e2))
lines(swd_vrkm3[,2],col="red")
lines(swd_vrkm3[,3],col="blue")
lines(swd_vrkm3[,4],col="green")

sum_wdkm3=data.frame(year=yearlist,swd_vrkm3)
sectors=c("Domestic", "Energy", "Industrial", "Livestock")
names(sum_wdkm3)[2:5]=sectors

wdkm3_long <- reshape2::melt(sum_wdkm3, id.vars = "year")

wdp=ggplot(wdkm3_long, aes(x = year, y = value, color = variable)) +
  geom_line(lwd=2) +
  labs(title = "",
       x = "Time",
       y = expression(paste("Volume (",km^3/year,")")),
       color = "") +
  scale_x_continuous(breaks=seq(1950,2020, by=10)) +
  scale_color_manual(values = c("black", "darkred", "red", "springgreen3","royalblue"),
                     guide = guide_legend(nrow = 2, keywidth = unit(3, "cm")))+
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf,
           fill = NA, color = "black", size = 3)+
  theme(axis.title=element_text(size=24),
        axis.text = element_text(size=22),
        panel.background = element_rect(fill = "white", colour = "white"),
        panel.grid = element_blank(),
        panel.border = element_rect(linetype = "solid", fill = NA, colour="black"),
        legend.title = element_text(size=22),
        legend.text = element_text(size=22),
        legend.position="bottom",
        panel.grid.major = element_line(colour = "grey60"),
        panel.grid.minor.y = element_line(colour = "grey80",linetype="dashed"),
        legend.key = element_rect(fill = "transparent", colour = "transparent"),
        legend.key.size = unit(.8, "cm"))
#save the plot
wdp
ggsave(paste0("D:/tilloal/Documents/06_Floodrivers/wateruse/Water_demand_volumesHERA.jpg"), wdp,width=30, height=20, units=c("cm"),dpi=300)

write.csv(sum_wdkm3,file=paste0(workDir,"/Datapaper/Revisions/water_demand_yearsAD.csv"))



# METEO PART ------------------------

#oad the meteo files
ncmet= paste0(workDir,"meteo/sums/et_ysum1951.nc")
ncmet= paste0(workDir,"meteo/sums/ta6CF_ysum1951.nc")
ncmet= paste0(workDir,"meteo/sums/ta/ta6_ysum1951.nc")

ncm=nc_open(ncmet)
nav=names(ncm[['var']])
#Band1 is the second variable
t=ncm$var[[2]]
name.var=names(ncm$var)[2]
tsize<-t$varsize

name.lon="lon"
name.lat="lat"

lon=ncvar_get(ncm,name.lon)
lat=ncvar_get(ncm,name.lat)
# lon=lon+0.01666667
# lat=lat+0.01666667
efas_lat=c(34.50,72.25)
efas_lon=c(-25.25, 35.00)
# lolon=c(which(round(lon,2)==efas_lon[1]),which(round(lon,2)==efas_lon[2]))
# lolat=c(which(round(lat,2)==efas_lat[1]),which(round(lat,2)==efas_lat[2]))
# 
# lon[2]-lon[1]
ll=expand.grid(lon,lat)
ll=as.data.frame(ll)
min(ll$Var2)
min(Area_domain$y)
ll$latlong=paste(round(ll$Var1,3),round(ll$Var2,3),sep=" ")

matcho=match(Area_domain$latlong,ll$latlong)
start=c(1,1,1)
count <- tsize 

met_svar=c()
fv=c("et","ta","pr")
var=c("et","ta6","pr6")
sce=c("","CF")
yearlist=c(1951:2020)
svar=3
for (svar in c(1:3)){
  for (ss in c(1,2)){
    v=paste0(var[svar],sce[ss])
    vf=fv[svar]
    print(v)
    name.var=var[svar]
    met_vsum=c()
    for (yr in yearlist){
      print(yr)
      ncmet= paste0(workDir,"meteo/sums/",vf,"/",v,"_ysum",yr,".nc")
      ncm=nc_open(ncmet)
      # nav=names(ncm[['var']])
      # #Band1 is the second variable
      # t=ncm$var[[2]]
      # name.var=names(ncm$var)[2]
      # tsize<-t$varsize
      # 
      # name.lon="lon"
      # name.lat="lat"
      # 
      # lon=ncvar_get(ncm,name.lon)
      # lat=ncvar_get(ncm,name.lat)
      # # lon=lon+0.01666667
      # # lat=lat+0.01666667
      # efas_lat=c(34.50,72.25)
      # efas_lon=c(-25.25, 35.00)
      # # lolon=c(which(round(lon,2)==efas_lon[1]),which(round(lon,2)==efas_lon[2]))
      # # lolat=c(which(round(lat,2)==efas_lat[1]),which(round(lat,2)==efas_lat[2]))
      # # 
      # # lon[2]-lon[1]
      # ll=expand.grid(lon,lat)
      # ll=as.data.frame(ll)
      # min(ll$Var2)
      # min(Area_domain$y)
      # ll$latlong=paste(round(ll$Var1,3),round(ll$Var2,3),sep=" ")
      # 
      # matcho=match(Area_domain$latlong,ll$latlong)
      # start=c(1,1,1)
      # count <- tsize 
      meteo_w   = ncvar_get(ncm,name.var,start = start, count= count) 
      
      if (v=="pr6" | v=="pr6CF" | v=="et" |v=="etCF" ){
        meteo_w=meteo_w*pixarea
        
        #divide by 4 since there are four values per day in mm/day
        meteo_w=meteo_w/4
      }
      #set the mask for only HERA domain
      meteo_w_dom=meteo_w[matcho]
      
      smwd=mean(meteo_w_dom,na.rm=T)
      #values are in mm/yr
      #convert to m/yr
      if (v=="pr6" | v=="pr6CF" | v=="et" |v=="etCF" ){
        meteo_w_dom=meteo_w_dom/1e3
        smwd=sum(meteo_w_dom,na.rm=T)
      }
      met_vsum=c(met_vsum,smwd)
    }
    met_svar=cbind(met_svar,met_vsum)
  }
}
met_svar=as.data.frame(met_svar)
colnames(met_svar)=c("et","etCF","ta6","ta6CF","pr6","pr6CF")
#transform m3 to Km3
met_svar[,c(1,2,5,6)]=met_svar[,c(1,2,5,6)]/1e9
plot(met_svar[,1],type="l")
lines(met_svar[,2],col=2)
lines(met_svar[,5],col=3)
lines(met_svar[,6],col=4)

sum_mkm3=data.frame(year=yearlist,met_svar)
#names(sum_mkm3)[c(2,3)]=c("Potential evapotranspiration","Precipitation")

metvol_long <- reshape2::melt(sum_mkm3, id.vars = "year")


metvol_et=metvol_long[which(metvol_long$variable=="et" | metvol_long$variable=="etCF"),]
metvol_pr=metvol_long[which(metvol_long$variable=="pr6" | metvol_long$variable=="pr6CF"),]
metvol_ta=metvol_long[which(metvol_long$variable=="ta6" | metvol_long$variable=="ta6CF"),]


# tmean=data.frame(year=yearlist,met_svar)
# 
# names(tmean)[c(2,3)]=c("factual Tmean","counterfactual Tmean")

#tm_long <- reshape2::melt(tmean, id.vars = "year")
yl=expression(paste("Volume (",km^3/year,")"))

pvar=c("et","ta6","pr6")
vp=pvar[3]

if (vp=="ta6"){
  yl=expression(paste("Mean Temperature (°C)"))
  metplot=metvol_ta
  
}
if (vp=="pr6"){
  metplot=metvol_pr
}
if (vp=="et"){
  metplot=metvol_et
}

wdp=ggplot(metplot, aes(x = year, y = value, color = variable)) +
  geom_line(lwd=1.4) +
  labs(title = "",
       x = "Time",
       y = yl,
       color = "") +
  scale_x_continuous(breaks=seq(1950,2020, by=10)) +
  scale_color_manual(values = c("darkred","darkblue"),
                     guide = guide_legend(nrow = 2, keywidth = unit(3, "cm")))+
  #scale_linetype_manual(values = c("solid", "dotted","solid", "dotted")) +
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf,
           fill = NA, color = "black", size = 3)+
  guides(linetype="none")+
  theme(axis.title=element_text(size=24),
        axis.text = element_text(size=22),
        panel.background = element_rect(fill = "white", colour = "white"),
        panel.grid = element_blank(),
        panel.border = element_rect(linetype = "solid", fill = NA, colour="black"),
        legend.title = element_text(size=22),
        legend.text = element_text(size=22),
        legend.position="bottom",
        panel.grid.major = element_line(colour = "grey60"),
        panel.grid.minor.y = element_line(colour = "grey80",linetype="dashed"),
        legend.key = element_rect(fill = "transparent", colour = "transparent"),
        legend.key.size = unit(.8, "cm"))
#save the plot
wdp
ggsave(paste0(workDir,"/Datapaper/Revisions/ag",vp,".jpg"), wdp,width=20, height=20, units=c("cm"),dpi=300)


#this could be used to compute statistics for any region 


#reuse that code to extract points in different Hydroregions and see the changes in the drivers

### HydroRegions ----

GridHR=raster( paste0(hydroDir,"/HydroRegions_raster_WGS84.tif"))
GHR=as.data.frame(GridHR,xy=T)
GHR=GHR[which(!is.na(GHR[,3])),]
GHR$llcoord=paste(round(GHR$x,4),round(GHR$y,4),sep=" ") 
GHR_riv=inner_join(GHR,outf,by= c("llcoord"="latlong"))
GHshpp <- read_sf(dsn ="Z:/ClimateRun4/nahaUsers/tilloal/HydroRegions/her_all_adjusted.shp")
HydroRsf=fortify(GHshpp) 


# Load the raster package
library(raster)

#load my shitty rdata
outletname="efas_rnet_100km_01min"

outf=outletopen(hydroDir,outletname)
#outf$latlong=paste(round(outf$Var1,4),round(outf$Var2,4),sep=" ")



main_path = 'D:/tilloal/Documents/06_Floodrivers/'
valid_path = paste0(main_path,'DataPaper/')
outletname="/GIS/upArea_European_01min.nc"
dir=valid_path
outf$idlalo=paste(outf$idlo, outf$idla, sep=" ")
UpArea=UpAopen(valid_path,outletname,outf)
head(UpArea)

coord_upa = SpatialPoints(UpArea[,c(1,2)], proj4string=CRS("+proj=longlat"))
cord_upat <- spTransform(coord_upa, CRS("+init=epsg:3035"))


rast_k=m3[,c(4,5,6)]

rk2 <- rasterFromXYZ(rast_k)


values <- extract(rk2, cord_upat)




write.csv(sum_mkm3,file=paste0(workDir,"/Datapaper/Revisions/meteo_volumes_years.csv"))
#now save everything to csv





zoneNA=upArea[c(lolon[2]:count[1]),c(lolat[1]:count[2])]
upArea[c(lolon[2]:count[1]),]=NA  
upArea[,c(lolat[1]:count[2])]=NA
upAc=as.vector(upArea)
upA=as.vector(upArea[which(!is.na(upArea))])


#reload lon and lat
lon=ncvar_get(ncd,name.lon)
lat=ncvar_get(ncd,name.lat)
ll=expand.grid(lon,lat)
upArea=data.frame(upAc,ll)
upArea$llcoord=paste(round(upArea$Var1,4),round(upArea$Var2,4),sep =" ")
narm=which(is.na(upArea$upAc))
upAreal=upArea[-narm,]



#length(unique(gHybas09$hybas09_raster))
# gHybas09l=gHybas09[-narm,]
# gHybas09l$llcoord=paste(round(gHybas09l$x,4),round(gHybas09l$y,4),sep=" ")

#this is where I can modify the output file

upAreaSub=upAreal[which( upAreal$upAc>=20 &  upAreal$upAc<100),]




####Extracting only pixels belonging to a given area

#Step 1: load the shapefile containing countries

countries <- read_sf(dsn = paste0(hydroDir,"/Countries/Countries_EFAS.shp"))

#France=countries[which(countries$NAME_FREN=="France"),]

points <- st_as_sf(upAreaSub, coords = c("Var1", "Var2"), crs = 4326)
points <- st_transform(points, crs = 3035)

pixOut <- st_intersection(points, countries)
mac=which(!is.na((match(upAreaSub$llcoord,pixOut$llcoord))))
pixOut=pixOut[,c(1:3)]
st_geometry(pixOut)=NULL
upAreaSEU=upAreaSub[mac,]
pixOut=inner_join(pixOut,upAreaSub,by="llcoord")



############################################################
###########  MATCHING PIXELS WITH CATCHMENTS ################
############################################################




outMatch=which(!is.na((match(upAreal$llcoord,gHybas09l$lloc))))
library(dplyr)
Rivdat=right_join(upAreal,gHybas09l, by= "llcoord")

Rivdat100=Rivdat[which(Rivdat$upAc>=100),]
Rivdat100 %>%
  ggplot(aes(x = Var1, y = Var2, fill = hybas09_raster)) +
  geom_tile()+
  scale_fill_distiller(palette = "Spectral",
                       direction = -1,oob = scales::squish) 


Catall=inner_join(hybasf9,Rivdat,by= c("SORT"="hybas09_raster"))


filtered_data <- list()
Rivdat=Rivdat[which(!is.na(Rivdat$hybas09_raster)),]
max(unique(Rivdat$hybas09_raster),na.rm=T)
# Iterate over unique values of HYBAS_ID
indic=1
for (i in unique(Rivdat$hybas09_raster)) {
  print(indic)
  print(i)
  # i=unique(Rivdat$HYBAS_ID)[18]
  # Get the subset of rows with the current HYBAS_ID value
  subset <- Rivdat[Rivdat$hybas09_raster == i, ]
  
  # Find the row with the maximum catch_area value
  subset=subset[order(subset$upAc,decreasing = T),]
  #I take the number 3 to avoid getting a river from another catchment at the mouth
  if (length(subset$upAc)>3) {
    sel_row <- subset[3, ]
  }else{
    sel_row <- subset[1, ]
  }
  
  
  indic=indic+1
  # Append the max row to the filtered data list
  filtered_data[[i]] <- sel_row
}

# Combine the filtered data frames into a single data frame
final_data <- do.call(rbind, filtered_data)
#remove the 0
final_data <- final_data[-1,]

final_data %>%
  ggplot(aes(x = Var1, y = Var2, fill = hybas09_raster)) +
  geom_tile()+
  scale_fill_distiller(palette = "Spectral",
                       direction = -1,oob = scales::squish) 

final_data$flagsmall=0
final_data$flagsmall[which(final_data$upAc<=100)]=1
###NOw need to find a way to give neighbour's attribute is there is no point > 100km2

final_geom=inner_join(hybasf9,final_data,by= c("SORT"="hybas09_raster"))

final_df=dplyr::select(as.data.frame(final_geom), -geometry)
plot(final_df$Var1,final_df$Var2,pch=16,col=2)
points(final_df$Var1[which(final_df$flagsmall==1)],final_df$Var2[which(final_df$flagsmall==1)],pch=16,col=3)


final_small=final_df[which(final_df$flagsmall==1),]
final_big=final_df[which(final_df$flagsmall==0),]
test=aggregate(list(flag=final_df$flagsmall),
               by = list( bas=final_df$MAIN_BAS),
               FUN = function(x) c(len=length(x),sum=sum(x)))
test <- do.call(data.frame, test)    

#A loop to find closest catchment
id=unique(final_small$SORT)
bID=c()
for (idu in id){
  smally=final_small[which(final_small$SORT==idu),]
  cid=smally$MAIN_BAS
  biggy=final_big$SORT
  bigID=biggy[which.min(abs(smally$SORT-biggy))]
  if (length(biggy)<1)bigID=NA
  bID=c(bID,bigID)
}
final_small$bigID=bID

#save final_df and final_small

################################################################################
################################################################################


final_df=pixOut
names(final_df)[1]="upAc"
#create the ncdf outlet

j2 <- sapply(final_df$Var1, function(x) which.min(abs(lon-x)))
k2 <- sapply(final_df$Var2, function(x) which.min(abs(lat-x)))

j2=as.matrix(j2)
llon=length(lon)
llat=length(lat)

fillvalue <- NA
# partial loop avoidance for tmp_array3
temp_array <- array(fillvalue, dim=c(llon,llat))

temp_array[cbind(j2,k2)] <- as.matrix(final_df$upAc) 


#function to create netcdf

history = 'Created Nov 2022' #####
Conventions = 'CF-1.6'
Source_Software = 'R netCDF4'
reference = 'JRC Climate risk team'  #####
title = 'Lisflood area maps for EUROPE setting Feb. 2023'
keywords = 'Lisflood, Global'
source = 'JRC Ispra'
institution = 'European Commission - Economics of climate change Unit (JRC.C.6) : https://ec.europa.eu/jrc/en/research-topic/climate-change'
comment = 'no.'

ncdf_creator<- function(arrin, ncname ,path, lon, lat, tsize, longname,name.var){
  
  ncfname <- paste(path, ncname, ".nc", sep="")
  dname <- name.var  # note: tmp means temperature (not temporary)
  # create and write the netCDF file -- ncdf4 version
  # define dimensions
  londim <- ncdim_def("lon","degrees_east",as.double(lon)) 
  latdim <- ncdim_def("lat","degrees_north",as.double(lat)) 
  
  
  arret0=as.array(arrin)
  # define variables
  fillvalue <- 0.0
  dlname <- longname
  if (length(tsize)==3){
    timedim <- ncdim_def("time",tunits2,as.double(time2))
    tmp_def <- ncvar_def(name.var,t$units,list(londim,latdim,timedim),fillvalue,dlname,prec="float",compression=4)
    proj <-ncvar_def("wgs_1984","1",NULL,NULL,longname="wgs_1984",prec="integer")
  }else{
    tmp_def <- ncvar_def(name.var,t$units,list(londim,latdim),fillvalue,dlname,prec="float",compression=1)
    #proj <-ncvar_def("wgs_1984","1",NULL,NULL,longname="wgs_1984",prec="integer")
  }
  
  
  
  # create netCDF file and put arrays
  ncout <- nc_create(ncfname,tmp_def,force_v4=TRUE)
  
  # put variables
  ncvar_put(ncout,tmp_def,arret0)
  
  
  # put additional attributes into dimension and data variables
  ncatt_put(ncout,"lon","axis","X") #,verbose=FALSE) #,definemode=FALSE)
  ncatt_put(ncout,"lat","axis","Y")
  
  ncatt_put(ncout,tmp_def,"standard_name",longname)
  # ncatt_put(ncout,tmp_def,"grid_mapping","wgs_1984")
  # ncatt_put(ncout,tmp_def,"esri_pe_string",'GEOGCS["WGS 84",DATUM["WGS_1984",SPHEROID["WGS 84",6378137,298.257223563,AUTHORITY["EPSG","7030"]],AUTHORITY["EPSG","6326"]],PRIMEM["Greenwich",0],UNIT["degree",0.0174532925199433],AUTHORITY["EPSG","4326"]]')
  ncatt_put(ncout,tmp_def,"missing_value",as.numeric(0.0))
  
  # put the CRS attributes
  # projname <- "wgs_1984"
  # ncatt_put(ncout,proj,"name",projname)
  # ncatt_put(ncout,proj,"long_name",projname)
  # ncatt_put(ncout,proj,"grid_mapping_name","latitude_longitude")
  # ncatt_put(ncout,proj,"semi_major_axis", 6378137)
  # ncatt_put(ncout,proj,"inverse_flattening", 298.257223563)
  # ncatt_put(ncout,proj,"proj4_params", "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")
  # ncatt_put(ncout,proj,"EPSG_code","EPSG:4326")
  
  
  # add global attributes
  ncatt_put(ncout,0,"title",title)
  ncatt_put(ncout,0,"institution",institution)
  ncatt_put(ncout,0,"source",source)
  ncatt_put(ncout,0,"references",reference)
  history <- paste("A.M. Tilloy", date(), sep=", ")
  ncatt_put(ncout,0,"history",history)
  ncatt_put(ncout,0,"Conventions",Conventions)
  
  # Get a summary of the created file:
  return(ncout)
  
  #nc_close(ncout)
}

dataDir<-("D:/tilloal/Documents/LFRuns_utils/data")
ncpath <- paste0(dataDir,"/")
# ncname <- "outlets_hybas09_01min" 
ncname <- "EUstreams_20to100_01min" 
nvar="streams"
temp_arfl=as.double(temp_array)
ncfout=ncdf_creator(temp_arfl, ncname, ncpath,lon=lon, lat=lat, tsize=2, longname="outlets EU streams 20 to 100 km2 upstream area",name.var=nvar)
nc_close(ncfout)


write.csv(final_df,file=paste0(hydroDir,"/hybas09_attributes.csv"))
write.csv(final_small,file=paste0(hydroDir,"/hybas09_small_attributes.csv"))

#load changes in precipitation and et0

#load bivariate change in flood and drought




#Identification of pre and post 1951 reservoirs in EFAS-----
# Library calling --------------------------------------------------
library(rgdal)
library(raster)
library(rgdal)
library(raster)
library(ncdf4)
library(lubridate)
library(ggplot2)
library(rasterVis)
library(tidyverse)
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)
library(rgeos)

dataDir<-("D:/tilloal/Documents/LFRuns_utils/data")
# Load the CSV file into a data frame
data <- read.csv(paste0(dataDir,"/Catchments/from_hybas_eu_allATTR.csv"), header = TRUE)

# Convert the HYBAS_ID column to a factor
data$HYBAS_ID <- as.character(data$HYBAS_ID)

data$llcoord=paste(round(data$POINT_X,4),round(data$POINT_Y,4),sep =" ")
# Create an empty list to store the filtered data frames
dir=hydroDir

# Function that open netcdf outlet files
outletopen=function(dir,outletname,nrspace=rep(NA,5)){
  ncbassin=paste0(dir,"/",outletname,".nc")
  ncb=nc_open(ncbassin)
  name.vb=names(ncb[['var']])
  namev=name.vb[1]
  if ("Band1"%in% name.vb)namev="Band1"
  name.lon="lon"
  name.lat="lat"
  if (!is.na(nrspace[1])){
    start=as.numeric(nrspace[c(2,4)])
    count=as.numeric(nrspace[c(3,5)])-start+1
  }else{
    londat = ncvar_get(ncb,name.lon) 
    llo=length(londat)
    latdat = ncvar_get(ncb,name.lat)
    lla=length(latdat)
    start=c(1,1)
    count=c(llo,lla)
  }
  
  londat = ncvar_get(ncb,name.lon,start=start[1],count=count[1]) 
  llo=length(londat)
  latdat = ncvar_get(ncb,name.lat,start=start[2],count=count[2])
  lla=length(latdat)
  outlets = ncvar_get(ncb,namev,start = start, count= count) 
  outlets=as.vector(outlets)
  outll=expand.grid(londat,latdat)
  lonlatloop=expand.grid(c(1:llo),c(1:lla))
  outll$idlo=lonlatloop$Var1
  outll$idla=lonlatloop$Var2
  
  outll=outll[which(!is.na(outlets)),]
  outlets=outlets[which(!is.na(outlets))]
  outll=data.frame(outlets,outll)
  return (outll)
}
#Import reservoir locations from netcdf
resOpen=function(dir,outletname){
  ncbassin=paste0(dir,outletname)
  ncb=nc_open(ncbassin)
  name.vb=names(ncb[['var']])
  namev=name.vb[1]
  #time <- ncvar_get(ncb,"time")
  
  #timestamp corretion
  name.lon="lon"
  name.lat="lat"
  londat = ncvar_get(ncb,name.lon) 
  llo=length(londat)
  latdat = ncvar_get(ncb,name.lat)
  lla=length(latdat)
  start=c(1,1)
  count=c(llo,lla)
  
  
  londat = ncvar_get(ncb,name.lon,start=start[1],count=count[1]) 
  llo=length(londat)
  latdat = ncvar_get(ncb,name.lat,start=start[2],count=count[2])
  lla=length(latdat)
  outlets = ncvar_get(ncb,namev,start = start, count= count) 
  outlets=as.vector(outlets)
  outll=expand.grid(londat,latdat)
  lonlatloop=expand.grid(c(1:llo),c(1:lla))
  outll$res=outlets
  outll$idlo=lonlatloop$Var1
  outll$idla=lonlatloop$Var2
  
  outll$idlalo=paste(outll$idlo,outll$idla,sep=" ")
  outfinal=outll[which(!is.na(outll$res)),]
  return (outfinal)
}

hydroDir<-("D:/tilloal/Documents/LFRuns_utils/data")
#comparison file between 2020 and 1951

res2020=resOpen(hydroDir,"/reservoirs/reservoirs_volumes_2020_Domain2.nc")
res2020$idla=2970-res2020$idla+1
res2020$idlalo=paste(res2020$idlo,res2020$idla,sep=" ")
res1951=resOpen(hydroDir,"/reservoirs/reservoirs_volumes_1951.nc")

max(res2020$res)
matres=na.omit(match(res1951$idlalo,res2020$idlalo))
res_old=res2020[matres,]
res_new=res2020[-matres,]

res_comp=left_join(res_old,res1951,by="idlalo")


### [Plot] - Figure S4 - Old and new reservoirs included in HERA -----

palet2=c(hcl.colors(9, palette = "Blues", alpha = NULL, rev = TRUE, fixup = TRUE))
outletname="efas_rnet_100km_01min"
outll=outletopen(hydroDir,outletname)
cord.dec=outll[,c(2,3)]
cord.dec = SpatialPoints(cord.dec, proj4string=CRS("+proj=longlat"))
cord.UTM <- spTransform(cord.dec, CRS("+init=epsg:3035"))
nco=cord.UTM@coords
world <- ne_countries(scale = "medium", returnclass = "sf")
Europe <- world[which(world$continent == "Europe"),]
e2=st_transform(Europe,  crs=3035)
w2=st_transform(world,  crs=3035)
tsize=12
osize=12
Impdates=seq(1950,2020,by=10)
valuenames=paste0("Y",Impdates)
basemap=w2

colorz = c("new"="orange","old" ='tomato4')
lab1=c("post-1951","pre-1951")
res_new$status="new"
res_old$status="old"

res_f=rbind(res_new,res_old)
pointout <- st_as_sf(res_f, coords = c("Var1", "Var2"), crs = 4326)
pointout <- st_transform(pointout, crs = 3035)

rnetwork=st_as_sf(outll, coords = c("Var1", "Var2"), crs = 4326)
rnetwork <- st_transform(rnetwork, crs = 3035)



damm<-ggplot(basemap) +
  geom_sf(fill="gray95",color="gray10",size=0.5)+
  geom_sf(data=pointout,aes(geometry=geometry,size=res,col=status),alpha=.9,stroke=0,shape=16)+
  geom_sf(data=rnetwork,aes(geometry=geometry),col="royalblue",size=0.1,alpha=.9,stroke=0,shape=15)+
  coord_sf(xlim = c(min(nco[,1]),max(nco[,1])), ylim = c(min(nco[,2]),max(nco[,2])))+
  scale_colour_manual(values = colorz, name="Dam construction", labels=lab1) +
  # scale_fill_manual(values = colorz, name="Largest change driver",labels=lab1) +
  scale_size(range = c(0.1, 5), trans="sqrt",name= expression(paste("Reservoir volume ", (m^3),sep = " ")),
             breaks=c(1e+5,1e+6,1e+7,1e+8,1e+9,1e+10))+
  #ggtitle(title)+
  labs(x="Longitude", y = "Latitude")+
  guides(colour = guide_legend(override.aes = list(size = 10)))+
  theme(axis.title=element_text(size=tsize),
        title = element_text(size=osize),
        axis.text=element_text(size=osize),
        panel.background = element_rect(fill = "aliceblue", colour = "grey1"),
        panel.border = element_rect(linetype = "solid", fill = NA, colour="black"),
        legend.title = element_text(size=tsize),
        legend.text = element_text(size=osize),
        legend.position = "right",
        panel.grid.major = element_line(colour = "grey70"),
        panel.grid.minor = element_line(colour = "grey90"),
        legend.key = element_rect(fill = "transparent", colour = "transparent"),
        legend.key.size = unit(1, "cm"))

ggsave(paste0("D:/tilloal/Documents/LFRuns_utils/TrendAnalysis/plots/DamKontraktion.jpg"), damm, width=20, height=20, units=c("cm"),dpi=1000) 
