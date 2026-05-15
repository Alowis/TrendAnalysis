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

#load dem
dem <- ("D:/tilloal/Documents/LFRuns_utils/ChangingHydroExtremes/data/GeoData/dem.nc")



# 2. Extract Specific Catchment
hybas_eu <- read_sf(paste0(hydroDir, "/Catchments/hydrosheds/hybas_eu_lev05_v1c.shp"))
target_catchment <- hybas_eu %>% filter(HYBAS_ID == selectedHybas)

# 1. Calculate the bounding box of your target catchment for the zoom
bbox <-terra::ext(target_catchment)+1
GHshpp      <- read_sf("Z:/ClimateRun4/nahaUsers/tilloal/HydroRegions/her_all_adjusted_cor.shp")

#target_catchment  <- GHshpp[GHshpp$CODEB == 74,  ]
# Load as stacks and crop/mask once to the catchment

# ── 1. Define your projections ─────────────────────────────────
crs_wgs84 <- "EPSG:4326"
crs_laea  <- "EPSG:3035"   # ETRS89-LAEA (standard for Europe)

dem2=rast(dem)
dem_bbox_crop <- terra::crop(dem2, (bbox))
dem_laea      <- terra::project(dem_bbox_crop, crs_laea, method = "bilinear")
rm(dem2)
gc()

slope  <- terra::terrain(dem_laea, "slope",  unit = "radians")
aspect <- terra::terrain(dem_laea, "aspect", unit = "radians")
hill   <- terra::shade(slope, aspect, angle = 45, direction = 315)

hill_df <- as.data.frame(hill, xy = TRUE, na.rm = TRUE)
colnames(hill_df)[3] <- "shade"
ggplot() +
  geom_raster(data = hill_df, aes(x = x, y = y, fill = shade)) +
  scale_fill_continuous(palette = "Light Grays",
                       na.value = "transparent", guide = "none")

#dem_df <- as.data.frame(dem_bbox_crop, xy = TRUE, na.rm = TRUE)

# ── 2. Reproject the DEM to LAEA ───────────────────────────────
# Use method = "bilinear" for continuous data like elevation
dem_laea <- terra::project(dem_bbox_crop, crs_laea, method = "bilinear")

# ── 3. Convert to data frame for ggplot ────────────────────────
dem_df <- as.data.frame(dem_laea, xy = TRUE, na.rm = TRUE)
colnames(dem_df)[3] <- "elevation"   # rename the value column


# Transform to LAEA and keep only regions that intersect the catchment
her_rhone <- st_transform(GHshpp, st_crs(hybas_eu)) %>%
  st_filter(target_catchment, .predicate = st_intersects)

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
  s1 <- raster::mask(crop(stack(files1), poly), poly)
  s2 <- raster::mask(crop(stack(files2), poly), poly)
  
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
    s <- raster::mask(crop(stack(files), poly), poly)
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
rast_pr_diff <- compute_period_diff("pr6", poly=target_catchment,dir=mdir)
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

forest_crop <- raster::mask(crop(rast_forest, target_catchment), target_catchment)
sealed_crop <- raster::mask(crop(rast_sealed, target_catchment), target_catchment)


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
  geom_tile(data = dem_df, 
            aes(x = x, y = y, 
                fill = layer), 
            linewidth = 0) +
  scale_fill_gradient2(low = "#8c510a", mid = "white", high = "#01665e", 
                       midpoint = 0, na.value = "transparent", name = "Forest Change")
  # geom_tile(data = forest_df, 
  #           aes(x = x, y = y, 
  #               fill = fracforest_ch20201951, 
  #               alpha = abs(fracforest_ch20201951)), 
  #           linewidth = 0) +
  # scale_fill_gradient2(low = "#8c510a", mid = "white", high = "#01665e", 
  #                      midpoint = 0, na.value = "transparent", name = "Forest Change") +
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

# 
# # 1. Define your manual thresholds (3 classes require 2 internal split points)
# # Anything below -20 is "1", -20 to 20 is "2", above 20 is "3"
# p_cuts  <- c(-20, 20)
# et_cuts <- c(-20, 20)
# 
# # 2. Manually assign classes 1, 2, and 3
# # findInterval returns 1, 2, or 3 based on where the value sits
# df_climate$P_manual_class  <- findInterval(df_climate$P_change, p_cuts, left.open = FALSE) + 1
# df_climate$ET_manual_class <- findInterval(df_climate$ET_change, et_cuts, left.open = FALSE) + 1
# 
# # 3. Constrain values to 1-3 (in case of NAs or extremes)
# df_climate$P_manual_class  <- pmin(pmax(df_climate$P_manual_class, 1), 3)
# df_climate$ET_manual_class <- pmin(pmax(df_climate$ET_manual_class, 1), 3)
# table(df_climate$ET_manual_class)
# # 4. Create the final 'bi_class' string that ggplot needs
# # This creates the "1-1", "2-3", etc. labels manually
# df_bi <- df_climate %>%
#   mutate(bi_class = paste0(P_manual_class, "-", ET_manual_class))
# 
# # 5. Check your work
# table(df_bi$bi_class)
# # Create categorical rasters (1, 2, 3)
# # 'findInterval' works on vectors, so we wrap it in 'calc' or use raster math
# rast_p_class <- reclassify(rast_pr_diff, c(-Inf, p_cuts[1], 1,  p_cuts[1], p_cuts[2], 2,  p_cuts[2], Inf, 3))
# rast_et_class <- reclassify(rast_et_final, c(-Inf, et_cuts[1], 1,  et_cuts[1], et_cuts[2], 2,  et_cuts[2], Inf, 3))

# 4. Apply your Manual Bivariate Breaks
# Adjust these cuts based on your specific climate range
# table(st_geometry_type(sub_catchments))
# sub_catchments <- sub_catchments[st_is(sub_catchments, c("POLYGON", "MULTIPOLYGON")), ]
# sub_catchments$P_change<- exact_extract(rast_pr_diff, sub_catchments, 'mean')
# sub_catchments$ET_change <- exact_extract(rast_et_diff, sub_catchments, 'mean')

# 
# sub_catchments <- sub_catchments %>%
#   mutate(
#     P_class = findInterval(P_change, p_cuts) + 1,
#     ET_class = findInterval(ET_change, et_cuts) + 1,
#     bi_class = paste0(pmin(pmax(P_class, 1), 3), "-", pmin(pmax(ET_class, 1), 3))
#   )
# 

sub_catchments <- sub_catchments[st_is(sub_catchments, c("POLYGON", "MULTIPOLYGON")), ]
sub_catchments$AI_change <- exact_extract(rast_delta_ai, sub_catchments, 'mean')
#load bivariate changes
hydroDir<-("D:/tilloal/Documents/LFRuns_utils/ChangingHydroExtremes/data")

load (file=paste0(hydroDir,"/CLchanges_bivariate_F.Rdata"))
load (file=paste0(hydroDir,"/SEchanges_bivariate_F.Rdata"))

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
  "2-2" = "#FFFFFF",  # medium-high x, medium-low y
  "3-2" = "#FFFFFF",  # medium-low x, medium-low y
  "4-2" = "#9cc4d2",  # low x, medium-low y
  
  "1-3" = "#635929",  # high x, medium-high y
  "2-3" = "#FFFFFF",  # medium-high x, medium-high y
  "3-3" = "#FFFFFF",  # medium-low x, medium-high y
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
her_rhone_laea <-st_transform(her_rhone, target_crs)



## 1. Transform to SF and reproject

dem_sf_laea <- dem_df %>%
  st_as_sf(coords = c("x", "y"), crs = 4326) %>%
  st_transform(3035)

# 2. Extract the new coordinates and convert BACK to a dataframe
df_dem_laea <- dem_sf_laea %>%
  mutate(x = st_coordinates(.)[,1],
         y = st_coordinates(.)[,2]) %>%
  st_drop_geometry() # This removes the 'spatial' part so it's a regular df


library(terra)

her_rhone_laea <- st_collection_extract(st_make_valid(her_rhone_laea[idx, ]), "POLYGON")


# 1. Convert RasterLayer to SpatRaster
rast_spat <- rast(rast_delta_ai)

# 2. Reproject using the EPSG code
rast_laea <- project(rast_spat, "EPSG:3035")

# 3. Convert to dataframe for ggplot
df_ai_laea <- as.data.frame(rast_laea, xy = TRUE, na.rm = TRUE)
colnames(df_ai_laea) <- c("x", "y", "delta_ai")
# C. Calculate the LAEA Bounding Box for the Zoom
bbox_laea <- st_bbox(target_catchment_laea)
her_rhone_laea=st_make_valid(her_rhone_laea)
# Keep only the part of each HER that lies within the catchment
her_rhone_inside <- st_crop(
  st_make_valid(her_rhone_laea),
  (bbox_laea)
)

ggplot(Europe_laea) +
  geom_sf(fill = "white", color = "grey80") +
  geom_sf(data = her_rhone_inside, fill = "pink", color = "black", alpha = 0.7) +
  coord_sf(xlim = c(bbox_laea["xmin"], bbox_laea["xmax"]),
           ylim = c(bbox_laea["ymin"], bbox_laea["ymax"]), expand = T)

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
# cities <- data.frame(
#   name = c("Alès", "Marseille", "Nice", "Barcelona"),
#   lon  = c(4.0825, 5.3698, 7.2620, 2.1734),
#   lat  = c(44.1275, 43.2965, 43.7102, 41.3851)
# )
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
  geom_sf(fill = "grey5", color = "white") +
  geom_sf(data = her_rhone_inside, fill = "grey80", color = NA, alpha = 0.5) +
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

ggplot()+
  geom_sf(data = her_rhone, fill = "grey80", color = NA, alpha = 0.5) 




s_plot=ggplot(Europe_laea) +
  # --- Layer 0: Europe Background ---
  geom_sf(fill = "grey87", color = "grey40") +
  
  #geom_sf(data = her_rhone_inside, fill = "darkgrey", color = NA, alpha = 0.5) +
  
  
  # ---- pre-layer: DEM -------------
  geom_raster(data = hill_df, aes(x = x, y = y, fill = shade), alpha=0.5) +
  scale_fill_continuous(palette = "Light Grays",
                        na.value = "transparent", guide = "none") +
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
          aes(size = total_wd_2020/1000, color = wd_change/1000),stroke=0, 
          alpha = 0.7) +
  geom_sf(data = nuts3_laea, 
          aes(size = total_wd_2020/1000),stroke=1, shape=21, color="black",
          alpha = 0.7) +
  scale_size_continuous(range = c(2, 12), name = "Total water demand 2020 (km³)",
                        guide = guide_legend(
                          direction = "horizontal", 
                          title.position = "top", 
                          nrow = 1, 
                          override.aes = list(shape = 21) ,
                          order = 3
                        )) + # <--- THIRD (Size)
  scale_color_gradient(low = "grey70", high = "mediumpurple4", name = "Water demand change (km³)",
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
           expand = F, clip="on") +
  theme_minimal() +
  theme(
    axis.title = element_blank(),
    axis.text = element_blank(),
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
s_plot
ggsave(paste0("D:/tilloal/Documents/LFRuns_utils/TrendAnalysis/plots/rhone_sechanges_FINAL.jpg"), width=20, height=20, units=c("cm"),dpi=1000)




#climate:



# 
# # Update your legend to match
# legend <- bi_legend(pal = "DkViolet2",
#                     dim = 4,
#                     xlab = "More Precip ",
#                     ylab = "More ET ",
#                     size = 8,
#                     flip_axes = TRUE)
# 
# library(cowplot)
# main_plot <- ggplot() +
#   # --- Layer 1: Bivariate Climate Background ---
#   #geom_sf(data = df_bi, aes( fill = bi_class), show.legend = FALSE,color="transparent") +
#   geom_sf(data = sub_catchments, aes( fill = AI_change),color="transparent") +
#   # Catchment Boundary
#   geom_sf(data = target_catchment, fill = "transparent", color = "gray23", size = 1) +
#   
#   scale_fill_gradientn(
#     colors=palet, breaks=br,limits=limi,trans=scales::modulus_trans(1),
#     oob = scales::squish)   +
#   
#   # bi_scale_fill(pal = "DkViolet2", dim = 4,flip_axes = TRUE) + 
#   # scale_fill
#   geom_sf(data = databipic_cropped, aes(color = bi_class),size=1,alpha=1,stroke=0,shape=15) +
#   
#   scale_color_manual(values = colorp, guide="none") +
#   # "BlueOr" fits your Red/Blue intent
#   theme(axis.title=element_text(size=tsize),
#         panel.background = element_rect(fill = "aliceblue", colour = "grey1"),
#         panel.border = element_rect(linetype = "solid", fill = NA, colour="black"),
#         legend.title = element_text(size=tsize),
#         legend.text = element_text(size=osize),
#         legend.position = "bottom",
#         panel.grid.major = element_line(colour = "grey70"),
#         panel.grid.minor = element_line(colour = "grey90"),
#         legend.key = element_rect(fill = "transparent", colour = "transparent"),
#         legend.key.size = unit(.8, "cm"))+
# labs(title = paste("Climatic changes | Rhone catchment"),
#      subtitle = "Yearly precipitation, Yearly PET  (1951-1980 vs 1991-2020)")
# 
# main_plot
# # Combine main_plot and the square legend
# final_map=ggarrange(main_plot, legend,
#              labels = c("Map", "Key"),
#              ncol = 2, nrow = 1,widths = c(2,1), heights=c(1,1), vjust=-1)
# # Display result
# print(final_map)

# 1. Calculate the bounding box for consistent zooming
#bbox <- st_bbox(target_catchment_laea)

# 2. Main Map Construction
br=seq(-.5,0.3,0.1)
limi=c(-0.4,0)
palet=c(hcl.colors(11, palette = "YlOrRd", alpha = NULL, rev = F, fixup = TRUE))
main_plot <- ggplot(Europe_laea) +
  # --- Layer 0: Europe Background ---

  geom_sf(fill = "grey87", color = NA) +
  geom_raster(data = hill_df, aes(x = x, y = y, fill = shade), alpha=0.5) +
  scale_fill_continuous(palette = "Light Grays",
                        na.value = "transparent", guide = "none") +
  new_scale_fill()+
  #geom_sf(data = her_rhone_inside, fill = "darkgrey", color = NA, alpha = 0.5) +
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
  coord_sf(xlim = c(bbox_laea["xmin"], bbox_laea["xmax"]), 
           ylim = c(bbox_laea["ymin"], bbox_laea["ymax"]), 
           expand = FALSE) +
  
  # --- Styling to match Socioeconomic Map ---
  # --- Styling ---
  theme_minimal() +
  theme(axis.title = element_blank(),
        axis.text = element_blank(),
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
ggsave(paste0("D:/tilloal/Documents/LFRuns_utils/TrendAnalysis/plots/rhone_Cchanges_FINAL.jpg"),main_plot, width=24, height=24, units=c("cm"),dpi=1000)


# # 3. Create/Update the Legend
# # Make sure the legend background is also transparent or matches theme
# legend <- bi_legend(pal = "DkViolet2",
#                     dim = 4,
#                     xlab = "More Precip",
#                     ylab = "More ET",
#                     size = 10,
#                     flip_axes = TRUE)
# 
# # 4. Final Combination
# final_map <- ggarrange(main_plot, s_plot,
#                        labels = c("", ""), # Removed "Map/Key" labels for a cleaner look
#                        ncol = 2, nrow = 1, 
#                        widths = c(40, 40))
# 
# print(final_map)

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


# Your original legend (without xlab/ylab because we'll use corner labels)
legend_biv <- bi_legend(
  pal = colorp, 
  dim = 4,
  xlab = " ",        # remove original axis titles
  ylab = " ",
  size = 10,
  arrows = FALSE
) + 
  ggtitle("Hydrological trajectories") + 
  theme(
    plot.title = element_text(
      size = tsize, 
      face = "bold", 
      hjust = 0.5,
      margin = margin(b = 5)
    ),
    legend.title = element_text(size = tsize, face = "bold"),
    plot.background = element_rect(fill = "transparent", colour = NA),
    panel.background = element_rect(fill = "transparent", colour = NA),
    axis.title = element_text(size = osize - 2),
    # Optional: remove axis ticks/text if you only want the corner labels
    axis.text = element_blank(),
    axis.ticks = element_blank()
  )

# Add the four corner labels
legend_biv <- legend_biv +
  annotate(
    "text",
    x = c(-0.2, 5.2, -.2, 5.2),        # bottom-left, bottom-right, top-right, top-left
    y = c(1, 1, 4, 4),
    label = c("DRY", "DEC",  "ACC","WET"),
    size = osize / 3,        # adjust size as needed (e.g., relative to osize)
    fontface = "bold",
    hjust=.5
  )+
  coord_equal(clip = "off")  

# Display
print(legend_biv)

legend_wrapped <- wrap_elements(panel = legend_biv)
# 
# design <- "AAAAABBBBBCC"
# combined_plot <- main_plot + s_plot + legend_wrapped + 
#   plot_layout(
#     design = design, 
#     widths = c(10, 10, 10), # Maps are 10x wider than the legend
#     guides = "collect"     # Collects standard legends (AI, Land Use, etc.) to the right
#   )
# 
# # Apply global theme changes
# cplot <- combined_plot & 
#   theme(
#     legend.position = "right",
#     legend.box = "vertical",
#    # plot.margin = margin(5, 5, 5, 5)
#   )
# 
# ggsave(
#   filename = "D:/tilloal/Documents/LFRuns_utils/TrendAnalysis/plots/Rhone_Combined_Analysis5.jpg",
#   plot = cplot,
#   width = 40,      # Doubled width for side-by-side
#   height = 20,     # Kept height same
#   units = "cm",
#   dpi = 1000        # 1000 is very high, 600 is usually enough for publication
# )

# 1. Combine your two maps and collect the standard legends on the right
# This creates the layout you had, but without the bivariate square
# 1. Base maps with collected guides
library(patchwork)
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
    right = 0.98,  
    
    # Adjust 'bottom' and 'top' to control height and vertical position
    bottom = 0.6, 
    top = 1,    
    
    align_to = 'full', 
    on_top = TRUE
  )

# 3. View the result
#final_cplot

ggsave(
  filename = "D:/tilloal/Documents/LFRuns_utils/TrendAnalysis/plots/Rhone_Combined_Analysis_FINAL.jpg",
  plot = final_cplot,
  width = 40,      # Doubled width for side-by-side
  height = 20,     # Kept height same
  units = "cm",
  dpi = 500        # 1000 is very high, 600 is usually enough for publication
)

