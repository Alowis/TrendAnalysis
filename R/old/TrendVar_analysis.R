################################################################################
#  TrendVar_analysis.R
#  Structured analysis of hazard trend & variability across scenarios
#
#  Workflow:
#    1. Setup & spatial data
#    2. Load scenario data with process_hazard_data()
#    3. Merge all scenarios into df_Main
#    4. Analyse 4 drivers x 3 variables with analyze_hazard_symmetric()
#    5. Save driver maps
#    6. Aggregate global trends with aggregate_global_trends()
#    7. Save aggregation plots
#    8. Save RL100mat (RP = 10 yr) for all scenarios
################################################################################


# ==============================================================================
# 0. LIBRARIES & SOURCE
# ==============================================================================

library(dplyr)
library(tidyr)
library(lubridate)
library(data.table)
library(ggplot2)
library(ggnewscale)
library(sf)
library(raster)
library(sp)
library(rnaturalearth)
library(reshape2)
library(scales)
library(Hmisc)

setwd("D:/tilloal/Documents/LFRuns_utils/TrendAnalysis/R/")
source("functions_trends.R")


analyze_hazard_symmetric <- function(param_df, 
                                     data_save_list, 
                                     var_type = "RP", 
                                     scenario = "WStat", 
                                     haz= "drought",
                                     base_year = 1955, 
                                     return_period = 100,
                                     years_to_agg = seq(1955, 2015, by = 10),
                                     eps = 1e-1) {
  require(dplyr)
  require(tidyr)
  print(scenario)
  # 1. Configuration: Scenario -> (Current Suffix, Baseline Suffix, Type)
  # 1. Configuration: Use if/else instead of case_when for lists
  if (scenario == "SCF") {
    cfg <- list(curr = "",   base = NA,   type = "climatic")
  } else if (scenario == "RWStat") {
    cfg <- list(curr = "RW", base = "",   type = "driver")
  } else if (scenario == "WStat") {
    cfg <- list(curr = "W",  base = "RW", type = "driver")
  } else if (scenario == "Histo") {
    cfg <- list(curr = "H",  base = "W",  type = "driver")
  } else {
    stop(paste("Unknown scenario provided:", scenario))
  }
  
  curr_sfx      <- cfg$curr
  base_sfx      <- cfg$base
  analysis_type <- cfg$type
  
  # Helper to get dynamic column names
  get_col <- function(v, sfx) {
    if (sfx == "") return(if(v == "Sigma") "sigmaGPD" else if(v == "threshold") "thresholdGPD" else v)
    return(paste0(v, sfx))
  }
  
  # 2. Logic for Driver Attribution (Comparison against a specific baseline)
  if (analysis_type == "driver") {
    
    curr_col_name <- get_col(var_type, curr_sfx)
    base_col_name <- get_col(var_type, base_sfx)
    
    if (var_type == "threshold") {
      # SIMPLE DIFFERENCE for Driver Thresholds
      if (haz=="flood"){
        param_df$calc_val <- (param_df[[curr_col_name]] - param_df[[base_col_name]])/
          ((param_df[[curr_col_name]] + param_df[[base_col_name]]) / 2 + eps) * 100
        # walla=param_df$calc_val[which(param_df$catchment==4200005)]
        # plot(walla)
      }
      if(haz=="drought"){
        param_df$calc_val <- (-param_df[[curr_col_name]] - (-param_df[[base_col_name]]))/
          ((-param_df[[curr_col_name]] - param_df[[base_col_name]]) / 2 + eps) * 100
      }
      
    } else if (var_type == "Sigma") {
      # RATIO for Driver Sigma
      param_df$calc_val <- param_df[[curr_col_name]] / (param_df[[base_col_name]])
      param_df$calc_val [param_df$calc_val > 10]=10
      
    } else if (var_type == "RP") {
      # Symmetric change for RP (Keeping this as requested previously, 
      # or you can change to ratio by removing the symmetric formula)
      val_curr <- calcGPDReturnLevel_Single(param_df$epsilonGPD, param_df[[get_col("Sigma", curr_sfx)]], 
                                            param_df[[get_col("threshold", curr_sfx)]], param_df$nPeaks, 
                                            70, return_period)
      val_base <- calcGPDReturnLevel_Single(param_df$epsilonGPD, param_df[[get_col("Sigma", base_sfx)]], 
                                            param_df[[get_col("threshold", base_sfx) ]], param_df$nPeaks, 
                                            70, return_period)
      if (haz=="drought"){
        val_curr<- -val_curr
        val_base<- -val_base
        mean(val_curr,na.rm=T)
        mean(val_base,na.rm=T)
        print((length(which(val_curr<0))/length(val_curr))*100)
        print((length(which(val_base<0))/length(val_curr))*100)
        val_curr[which(val_curr<0)]=0
        val_base[which(val_base<0)]=0
        
      }
      
      param_df$calc_val <- (val_curr - val_base) / ((val_curr + val_base) / 2 + eps) * 100
      
    } else {
      # Default: Simple subtraction for 'variability' or 'trend'
      param_df$calc_val <- param_df[[curr_col_name]] - param_df[[base_col_name]]
    }
    active_var <- "calc_val"
    
    # 3. Logic for Climatic Trend (SCF internal temporal change)
  } else {
    active_var <- if(var_type == "Sigma") "sigmaGPD" else if(var_type == "threshold") "thresholdGPD" else var_type
    
    if (var_type == "RP") {
      if (haz=="drought"){
        param_df$RP_val <- -calcGPDReturnLevel_Single(param_df$epsilonGPD, param_df$sigmaGPD, 
                                                      param_df$thresholdGPD, param_df$nPeaks, 
                                                      70, return_period)
        print(length(which(param_df$RP_val<0))/length(param_df$RP_val)*100)
        param_df$RP_val[which(param_df$RP_val<0)]=0
      }else{
        param_df$RP_val <- calcGPDReturnLevel_Single(param_df$epsilonGPD, param_df$sigmaGPD, 
                                                     param_df$thresholdGPD, param_df$nPeaks, 
                                                     70, return_period)
      }
      
      active_var <- "RP_val"
    }
  }
  
  if(haz=="Drought") param_df$thresholdGPD=-param_df$thresholdGPD
  # 4. Pivot and Aggregate

  
  df_wide <- param_df %>%
    dplyr::select(catchment, Year, !!sym(active_var)) %>%
    pivot_wider(names_from = Year, values_from = !!sym(active_var)) %>%
    as.data.frame()

  
  data_cols <- setdiff(names(df_wide), "catchment")
  
  if (analysis_type == "climatic") {
    base_yr_str <- as.character(base_year)
    # Applying the symmetric logic to the internal climatic trend as well (Relative to 1955)
    row_init <- df_wide[[base_yr_str]]
    
    if (var_type == "threshold" | var_type == "RP") {
      # (Year_X - Year_1955) / Mean(Year_X, Year_1955) * 100
      
      
      vals <- (df_wide[, data_cols] - (row_init)) / ((df_wide[, data_cols] + (row_init)) / 2 + eps) * 100
      #capping super high values
      vals[vals > 1e3]=1e3
      vals[vals < -1e3]=-1e3
      
    } else if (var_type == "Sigma") {
      vals <- df_wide[, data_cols] / row_init
      #capping super high values
      # max(vals,na.rm=T)
      # length(which(vals[,70]>10))
      vals[vals > 10]=10
    } else {
      vals <- df_wide[, data_cols]
    }
    DataV_vals <- vals
  } else {
    # Driver results are already symmetric percentages from step 2
    DataV_vals <- df_wide[, data_cols]
  }
  
  DataV_vals[DataV_vals > 1e4] <- NA
  # 5. Final Assembly
  DataV <- data.frame(outl2 = df_wide$catchment, DataV_vals)
  DataO <- data_save_list[[1]][, -c(12:81)]
  DataV <- inner_join(DataO, DataV, by = "outl2")
  
  yr_cols <- intersect(paste0("X", years_to_agg), names(DataV))
  if(length(yr_cols) == 0) yr_cols <- as.character(years_to_agg)
  
  trendClim <- aggregate(DataV[, yr_cols], by = list(HydroR = DataV$HER), 
                         FUN = function(x) mean(x, na.rm = TRUE)) %>% do.call(data.frame, .)
  
  col_2015 <- if("X2015" %in% names(DataV)) "X2015" else "2015"
  pointClim <- aggregate(list(Rchange_rel = DataV[[col_2015]]), by = list(HydroR = DataV$HER), 
                         FUN = function(x) {
                           c(mean = mean(x, na.rm = TRUE), med = median(x, na.rm = TRUE))
                         }) %>% do.call(data.frame, .)
  
  return(list(DataV = DataV, point = pointClim, trendClim = trendClim))
}


aggregate_global_trends <- function(trend_df,weight,
                                    decades = c(1955, 1965, 1975, 1985, 1995, 2005, 2015)) {
  require(reshape2)
  require(dplyr)
  
  # 1. Melt the regional trend data
  tr_melted <- suppressWarnings(
    melt(trend_df, id.vars = "HydroR", variable.name = "variable", value.name = "value")
  )
  
  # 2. Extract Year from variable names (handles 'X1955' or '1955')
  tr_melted$yr <- as.numeric(gsub("X", "", as.character(tr_melted$variable)))
  
  # 3. Filter for specified decades/years
  tr_filtered <- tr_melted %>%
    filter(yr %in% decades) %>%
    mutate(decad = yr)
  
  hmatct=match(tr_filtered$HydroR,weight$reg)
  tr_filtered$weight=weight$value[hmatct]
  
  # 4. Global aggregation with statistics
  # tr_global <- aggregate(
  #   list(value = tr_filtered$value),
  #   by = list(yr = tr_filtered$decad),
  #   FUN = function(x) {
  #     c(mean = mean(x, na.rm = TRUE),
  #       l    = length(x),
  #       med  = median(x, na.rm = TRUE),
  #       ql   = quantile(x, 0.25, na.rm = TRUE),
  #       qh   = quantile(x, 0.75, na.rm = TRUE),
  #       w1   = quantile(x, 0.025, na.rm = TRUE),
  #       w2   = quantile(x, 0.975, na.rm = TRUE))
  #   }
  # )
  
  tr_global <- tr_filtered %>%
    group_by(decad) %>%
    summarise(
      mean_w = weighted.mean(value, w = weight, na.rm = TRUE),
      l      = n(),
      med = wtd.quantile(value, weights = weight, probs = 0.5, na.rm = TRUE),
      ql   = wtd.quantile(value, weights = weight, probs = 0.25, na.rm = TRUE),
      qh   = wtd.quantile(value, weights = weight, probs = 0.75, na.rm = TRUE),
      w1   = wtd.quantile(value, weights = weight, probs = 0.025, na.rm = TRUE),
      w2   = wtd.quantile(value, weights = weight, probs = 0.975, na.rm = TRUE)
    )
  
  # 5. Convert matrix output to a clean data frame
  tDataHuman <- do.call(data.frame, tr_global)
  
  return(tDataHuman)
}


calcGPDReturnLevel_Single <- function(epsilon, sigma, threshold, nPeaks, sampleTimeHorizon, returnPeriod) {
  
  # Calculate expected number of events in T years
  XX <- (nPeaks / sampleTimeHorizon) * returnPeriod
  
  # Initialize result vector with NAs
  returnLevel <- rep(NA_real_, length(sigma))
  
  # 1. Handle GPD Case: epsilon is NOT NA and NOT zero
  idx_gpd <- !is.na(epsilon) & abs(epsilon) > 1e-9
  if (any(idx_gpd)) {
    returnLevel[idx_gpd] <- threshold[idx_gpd] + 
      (sigma[idx_gpd] / epsilon[idx_gpd]) * (XX[idx_gpd]^epsilon[idx_gpd] - 1)
  }
  
  # 2. Handle Gumbel Case: epsilon is NOT NA and IS zero
  idx_gum <- !is.na(epsilon) & abs(epsilon) <= 1e-9
  if (any(idx_gum)) {
    returnLevel[idx_gum] <- threshold[idx_gum] + 
      sigma[idx_gum] * log(XX[idx_gum])
  }
  
  # Note: Rows where epsilon/sigma/threshold were originally NA will remain NA
  return(returnLevel)
}


process_hazard_data <- function(sub_dir = "SCFX", 
                                hazard = "Flood", 
                                base_data_dir = dataDir, 
                                fileVar="Var_TrendX_agg.Rdata",
                                hydro_base_dir = "D:/tilloal/Documents/LFRuns_utils/ChangingHydroExtremes/data" ) {
  
  require(dplyr)
  require(tidyr)
  require(lubridate)
  require(data.table)
  
  # 1. Determine the filename suffix based on sub_dir
  # This handles cases like SCF -> SocCF, Histo -> Histo, WStat -> WStat
  file_suffix <- case_when(
    sub_dir == "SCFX"   ~ "SocCF2_revfF",
    sub_dir == "HistoX" ~ "Histo",
    sub_dir == "WStatX" ~ "WCF",
    sub_dir == "RWStatX" ~ "RWCF",
    TRUE               ~ sub_dir  # For WStat and RWStat
  )
  
  if (hazard == "Drought") {
    namefile <- paste0("Drought.nonfrost.", file_suffix)
  } else {
    namefile <- paste0("Flood.year.", file_suffix)
  }
  
  # 2. Construct paths and Load Files
  res_path <- file.path(base_data_dir, sub_dir, fileVar)
  param_path <- file.path(hydro_base_dir, hazard, paste0("params.", namefile, ".Rdata"))
  
  message("Loading: ", res_path)
  load(res_path, envir = .GlobalEnv) 
  
  message("Loading: ", param_path)
  load(param_path, envir = .GlobalEnv)
  
  #fixing errors
  #p43x=Paramsfl[which(Paramsfl$catchment==4304271),]
  # 3. Process Parameters (Paramsfl)
  # Selecting columns by index as per your original script
  rmcol=c("epsilonGEV", "sigmaGEV", "muGEV" ,"epsilonStdErrGEV",  "sigmaStdErrGEV"  ,  "muStdErrGEV"  ,"SampleTimeHorizon")
  #keep_cols <- setdiff(colnames(Paramsfl), rmcol)
  params_clean <- Paramsfl[, !..rmcol] %>% as.data.table()
  rm(Paramsfl, envir = .GlobalEnv)
  
  # 4. Process Variability (ResTV[[2]])
  var_long <- ResTV[[2]] %>%
    rename(time = 1) %>%
    pivot_longer(cols = -time, names_to = "catchment", values_to = "variability") %>%
    mutate(
      Year = year(time), 
      catchment = as.numeric(catchment)
    )
  
  # 5. Process Trend (ResTV[[1]])
  trend_long <- ResTV[[1]] %>%
    rename(time = 1) %>%
    pivot_longer(cols = -time, names_to = "catchment", values_to = "trend") %>%
    mutate(
      Year = year(time), 
      catchment = as.numeric(catchment)
    )
  if (hazard == "Drought") {
    trend_long$trend=-trend_long$trend
  }
  
  # 6. Merge and Calculate
  final_df <- var_long %>%
    full_join(params_clean, by = c("catchment", "Year")) %>%
    full_join(trend_long, by = c("catchment", "Year", "time")) %>%
    mutate(
      Sigma0 = sigmaGPD / variability,
      threshold0 = (thresholdGPD - trend) / variability
    )
  
  f2=final_df[which(final_df$Year==1951),]
  gc() # Clear memory after large loads
  return(final_df)
}

attach_scenario_cols <- function(base_df, source_df, suffix) {
  # Select only the joining keys and the columns we need to transfer
  to_add <- source_df %>%
    dplyr::select(catchment, Year, time, 
                  source_variability = variability, 
                  source_trend = trend)
  
  # Join and calculate the new scenario-specific parameters
  res <- base_df %>%
    inner_join(to_add, by = c("catchment", "Year", "time")) %>%
    mutate(
      !!paste0("variability", suffix) := source_variability,
      !!paste0("trend", suffix)       := source_trend,
      !!paste0("Sigma", suffix)       := Sigma0 * source_variability,
      !!paste0("threshold", suffix)   := threshold0 * source_variability + source_trend
    ) %>%
    dplyr::select(-source_variability, -source_trend) # Clean up temp cols
  
  return(res)
}



# ==============================================================================
# 1. PATHS & GLOBAL SETTINGS
# ==============================================================================

hazard   <- "Flood"    # "Flood" or "Drought"
haz      <- tolower(hazard)

# Root data directories
hydroDir <- "D:/tilloal/Documents/LFRuns_utils/ChangingHydroExtremes/data"
dataDir  <- switch(hazard,
  Flood   = "D:/tilloal/Documents/LFRuns_utils/data/Flood/HPC/Calibrated/revision/TrendVar/",
  Drought = "D:/tilloal/Documents/LFRuns_utils/data/Drought/HPC/Calibrated/revision/TrendVar/"
)
plotDir  <- "D:/tilloal/Documents/LFRuns_utils/TAplots/"

# Return period used throughout
RETURN_PERIOD <- 10

# Aggregation decades
DECADES <- seq(1955, 2015, by = 10)

# Driver labels (used in plots and file names)
DRIVER_LABELS  <- c("Climate", "Land use", "Reservoirs", "Water demand")
DRIVER_COLOURS <- c(
  "Clim" = "royalblue",
  "LUC"  = "orange",
  "Res"  = "tomato4",
  "WU"   = "limegreen"
)


# ==============================================================================
# 2. SPATIAL DATA
# ==============================================================================

# --- 2.1  River network outlets --------------------------------------------

rspace <- read.csv(paste0(hydroDir, "/subspace_efas.csv"))[, -1]

if (!exists("outf")) {
  outf <- c()
  for (Nsq in 1:88) {
    message("Reading outlet square ", Nsq)
    nrspace    <- rspace[Nsq, ]
    outletname <- "GeoData/efas_rnet_100km_01min"
    outhybas   <- outletopen(hydroDir, outletname, nrspace)

    Idstart  <- as.numeric(Nsq) * 10000
    Idstart2 <- as.numeric(Nsq) * 100000

    if (length(outhybas$outlets) > 0) {
      outhybas$outlets <- seq((Idstart + 1), (Idstart + length(outhybas$outlets)))
      outhybas$outl2   <- seq((Idstart2 + 1), (Idstart2 + length(outhybas$outlets)))
      outhybas$latlong <- paste(round(outhybas$Var1, 4), round(outhybas$Var2, 4), sep = " ")
      outf <- rbind(outf, outhybas)
    }
  }
}

# --- 2.2  HydroBASINS level 7 catchments -----------------------------------

Catchmentrivers7 <- read.csv(
  paste0(hydroDir, "/GeoData/HYBAS07/from_hybas_eu_onlyid.csv"),
  encoding = "UTF-8", header = TRUE, stringsAsFactors = FALSE
)
hybas07    <- read_sf(dsn = paste0(hydroDir, "/GeoData/HYBAS07/hybas_eu_lev07_v1c.shp"))
hybasf7    <- fortify(hybas07)
Catamere07 <- inner_join(hybasf7, Catchmentrivers7, by = "HYBAS_ID")
Catamere07$llcoord <- paste(round(Catamere07$POINT_X, 4), round(Catamere07$POINT_Y, 4), sep = " ")

catmap <- right_join(Catamere07, outf, by = c("llcoord" = "latlong"))
GNF    <- catmap
st_geometry(GNF) <- NULL
rm(Catamere07)

# Match HYBAS07 outlets to pixel IDs
outhybas07 <- outletopen(hydroDir, "/GeoData/HYBAS07/outletsv8_hybas07_01min")
outhybas07$latlong <- paste(round(outhybas07$Var1, 4), round(outhybas07$Var2, 4), sep = " ")
mhy <- match(outhybas07$latlong, outf$latlong)
outhybas07$outID <- outf$outl2[mhy]

# --- 2.3  European biogeographic regions -----------------------------------

biogeo      <- read_sf(dsn = paste0(hydroDir, "/GeoData/eea_3035_biogeo-regions_2016/BiogeoRegions2016_wag84.shp"))
biogeof     <- fortify(biogeo)
st_geometry(biogeof) <- NULL
biogeoregions <- raster(paste0(hydroDir, "/GeoData/eea_3035_biogeo-regions_2016/Biogeo_rasterized_wsg84.tif"))
Gbiogeoregions <- as.data.frame(biogeoregions, xy = TRUE)

biogeomatch <- inner_join(biogeof, Gbiogeoregions, by = c("PK_UID" = "Biogeo_rasterized_wsg84"))
biogeomatch$latlong <- paste(round(biogeomatch$x, 4), round(biogeomatch$y, 4), sep = " ")
biogeo_rivers <- right_join(biogeomatch, outf, by = "latlong")

# --- 2.4  HydroEcoRegions --------------------------------------------------

GridHR  <- raster(paste0(hydroDir, "/GeoData/HER/HydroRegions_raster_WGS84.tif"))
GHR     <- as.data.frame(GridHR, xy = TRUE)
GHR     <- GHR[!is.na(GHR[, 3]), ]
GHR$llcoord <- paste(round(GHR$x, 4), round(GHR$y, 4), sep = " ")
GHR_riv <- inner_join(GHR, outf, by = c("llcoord" = "latlong"))
GHshpp  <- read_sf(dsn = paste0(hydroDir, "/GeoData/HER/her_all_adjusted.shp"))
HydroRsf <- fortify(GHshpp)

# --- 2.5  Base-map & plot parameters ---------------------------------------

outletname <- "/GeoData/efas_rnet_100km_01min"
outll      <- outletopen(hydroDir, outletname)
cord.dec   <- SpatialPoints(outll[, c(2, 3)], proj4string = CRS("+proj=longlat"))
cord.UTM   <- spTransform(cord.dec, CRS("+init=epsg:3035"))
nco        <- cord.UTM@coords

world   <- ne_countries(scale = "medium", returnclass = "sf")
Europe  <- world[world$continent == "Europe", ]
e2      <- st_transform(Europe, crs = 3035)
w2      <- st_transform(world,  crs = 3035)
basemap <- w2

tsize <- 16
osize <- 12

# --- 2.6  Upstream area ----------------------------------------------------

outf$idlalo <- paste(outf$idlo, outf$idla, sep = " ")
UpArea <- UpAopen(hydroDir, "/GeoData/upArea_European_01min.nc", outf)


# ==============================================================================
# 3. LOAD SCENARIO DATA  (process_hazard_data — unchanged)
# ==============================================================================

message("\n=== Loading scenario data ===")

df_SCF    <- process_hazard_data(sub_dir = "SCFX",   hazard = hazard,
                                 base_data_dir = dataDir,
                                 fileVar = "Var_TrendX_agg.Rdata")

df_Histo  <- process_hazard_data(sub_dir = "HistoX", hazard = hazard,
                                 base_data_dir = dataDir,
                                 fileVar = "Var_TrendX_agg.Rdata")

df_WStat  <- process_hazard_data(sub_dir = "WStatX",  hazard = hazard,
                                 base_data_dir = dataDir,
                                 fileVar = "Var_TrendX_agg.Rdata")

df_RWStat <- process_hazard_data(sub_dir = "RWStatX", hazard = hazard,
                                 base_data_dir = dataDir,
                                 fileVar = "Var_TrendX_agg.Rdata")

message("All scenario data loaded.")


# ==============================================================================
# 4. MERGE SCENARIOS INTO df_Main
#
#  Symmetric attribution chain:
#    SCF  = full socio-climatic forcing (base columns: sigmaGPD, thresholdGPD)
#    RWStat -> suffix "RW"  (Reservoirs + Water demand static)
#    WStat  -> suffix "W"   (Water demand static only)
#    Histo  -> suffix "H"   (historical)
# ==============================================================================


message("\n=== Merging scenarios into df_Main ===")

df_Main <- attach_scenario_cols(df_SCF,   df_RWStat, "RW")
rm(df_RWStat); gc()

df_Main <- attach_scenario_cols(df_Main,  df_WStat,  "W")
rm(df_WStat);  gc()

df_Main <- attach_scenario_cols(df_Main,  df_Histo,  "H")
rm(df_Histo, df_SCF); gc()

# Convert to data.table once for speed
setDT(df_Main)
df_Main$nsq=floor(df_Main$catchment/100000)
s43=which(!is.na(match(df_Main$catchment,4304271)))
p43=df_Main[s43,]
duplicates <- p43 %>%
  group_by(catchment, Year) %>%
  summarise(n = n(), .groups = "drop") %>%
  filter(n > 1)

print(duplicates)

message("df_Main ready: ", nrow(df_Main), " rows, ",
        length(unique(df_Main$catchment)), " catchments.")


# ==============================================================================
# 5. DRIVER ANALYSIS WITH analyze_hazard_symmetric()
#
#  Drivers  | Scenario pair (curr vs base) | Meaning
#  ---------|------------------------------|-------------------------------
#  Climate  | SCF  (internal time trend)   | full climatic change
#  Reservoirs | RWStat vs SCF (RW vs "")   | land use effect
#  Land use   | WStat  vs RWStat  (W vs RW)| reservoir effect
#  WaterDemand| Histo  vs WStat   (H vs W) | water-demand effect
#
#  Variable types: "threshold", "Sigma", "RP"
# ==============================================================================

message("\n=== Driver analysis ===")

# Load DataSave (reference spatial metadata object used inside the function)
load(file = paste0("D:/tilloal/Documents/LFRuns_utils/data/TSEVA/output_plots/Flood_pixChange_RL100_v1.Rdata"))

#length of each HER for weight
Dtweight=DataSave$Total
Dtweight=aggregate(
  list(value = Dtweight$outlets),
  by = list(reg = Dtweight$HER),
  FUN = function(x) {
    c(l = length(x))}
)


# Map scenario labels to driver names
scenario_driver_map <- list(
  SCF    = "climate",
  RWStat = "laduse",
  WStat  = "reservoir",
  Histo  = "waterdemand"
)

var_types <- c("threshold", "Sigma", "RP")

# Storage for aggregation results (used in Section 6)
agg_results <- list()   # agg_results[[var_type]][[scenario]]

for (var_type in var_types) {

  message("\n--- Variable type: ", var_type, " ---")
  agg_results[[var_type]] <- list()

  # --- 5.1  Configure plot aesthetics per variable -------------------------

  if (var_type == "Sigma") {
    br      <- seq(0.5, 1.5, by = 0.1)
    limi    <- c(0.5, 1.5)
    legend2 <- "Sigma ratio"
    palet   <- hcl.colors(11, palette = "RdYlBu", rev = FALSE)
    paletf  <- hcl.colors(11, palette = "RdBu",   rev = FALSE)
  } else {                              # threshold or RP
    br      <- seq(-50, 50, by = 5)
    limi    <- c(-50, 50)
    legend2 <- if (var_type == "threshold") "Threshold change (%)" else "RL change (%)"
    palet   <- hcl.colors(11, palette = "RdYlBu", rev = FALSE)
    paletf  <- hcl.colors(11, palette = "RdBu",   rev = FALSE)
  }

  # --- 5.2  Loop over scenarios / drivers ----------------------------------

  for (sce in names(scenario_driver_map)) {
    driver <- scenario_driver_map[[sce]]
    message("  Scenario: ", sce, " | Driver: ", driver)

    res <- analyze_hazard_symmetric(
      param_df      = df_Main ,
      data_save_list = DataSave ,
      var_type      = var_type ,
      scenario      = sce,
      haz           = haz,
      base_year     = 1955 ,
      return_period = RETURN_PERIOD ,
      years_to_agg  = DECADES
    )

    # Store trend data for aggregation in Section 6
    agg_results[[var_type]][[sce]] <- res$trendClim

    # --- 5.3  Build map ----------------------------------------------------

    # Point layer
    points_sf <- res$DataV[!is.na(res$DataV$Var1), ]
    points_sf <- st_as_sf(points_sf, coords = c("Var1", "Var2"), crs = 4326)
    points_sf <- st_transform(points_sf, crs = 3035)

    # Add upstream area
    # mp_upa <- match(points_sf$outl2, outf$outl2)
    # points_sf$upa <- UpArea$upa[mp_upa]

    # Regional polygon layer
    pag <- inner_join(HydroRsf, res$point, by = c("CODEB" = "HydroR"))

    fmap <- ggplot(basemap) +
      geom_sf(fill = "white", color = "darkgrey", linewidth = 0.3) +
      geom_sf(data = pag,
              aes(fill = Rchange_rel.mean, geometry = geometry),
              alpha = 0.2, color = "transparent") +
      geom_sf(data = points_sf,
              aes(col = X2015, geometry = geometry, size = upa),
              alpha = 0.9, stroke = 0, shape = 15) +
      scale_size(range = c(0.08, 0.4), trans = "sqrt", guide = "none") +
      scale_fill_gradientn(
        colors = paletf, breaks = br, limits = limi,
        trans = scales::modulus_trans(0.3),
        oob = scales::squish, name = legend2) +
      scale_color_gradientn(
        colors = palet, breaks = br, limits = limi,
        trans = scales::modulus_trans(0.3),
        oob = scales::squish, name = legend2) +
      coord_sf(xlim = range(nco[, 1]), ylim = range(nco[, 2])) +
      guides(colour = guide_colourbar(barwidth = 22, barheight = 1), fill = "none") +
      labs(x = "Longitude", y = "Latitude",
           title = paste0(var_type, " – ", driver, " driver (", hazard, ")")) +
      theme(
        axis.title       = element_text(size = tsize),
        axis.text        = element_text(size = osize),
        title            = element_text(size = osize),
        panel.background = element_rect(fill = "aliceblue", colour = "grey1"),
        panel.border     = element_rect(linetype = "solid", fill = NA, colour = "black"),
        legend.title     = element_text(size = tsize),
        legend.text      = element_text(size = osize),
        legend.position  = "bottom",
        legend.box       = "vertical",
        panel.grid.major = element_line(colour = "grey70"),
        panel.grid.minor = element_line(colour = "grey90"),
        legend.key       = element_rect(fill = "transparent", colour = "transparent"),
        legend.key.size  = unit(1, "cm")
      )

    # --- 5.4  Save map -----------------------------------------------------

    map_file <- paste0(plotDir, "map_", var_type, "_", driver, "_", haz, "_FINAL.jpg")
    ggsave(map_file, fmap, width = 22, height = 20, units = "cm", dpi = 800)
    message("  Map saved: ", map_file)
  }
}


# ==============================================================================
# 6. GLOBAL TREND AGGREGATION  (aggregate_global_trends — unchanged)
# ==============================================================================

message("\n=== Global trend aggregation ===")

# Driver metadata for plotting
driver_meta <- data.frame(
  scenario = c("SCF", "RWStat", "WStat", "Histo"),
  driver   = c("Clim", "LUC", "Res", "WU"),
  label    = DRIVER_LABELS,
  stringsAsFactors = FALSE
)

for (var_type in var_types) {
 #var_type="threshold"
  message("\n--- Aggregating: ", var_type, " ---")

  if (var_type == "Sigma") {
    name  <- "scale"
    nplot <- "Scale parameter (ratio)"
    br_p  <- seq(0, 10, by = 0.25)
  } else if (var_type == "threshold") {
    name  <- "location"
    nplot <- "Location parameter – relative difference (%)"
    br_p  <- seq(-50, 50, by = 10)
  } else {
    name  <- "RL"
    nplot <- paste0("Return level (RP", RETURN_PERIOD, ") – relative difference (%)")
    br_p  <- seq(-50, 50, by = 10)
  }

  popo=agg_results[[var_type]][[sce]]
  # Aggregate each scenario and bind
  trtF <- NULL
  for (i in seq_len(nrow(driver_meta))) {
    sce    <- driver_meta$scenario[i]
    drv    <- driver_meta$driver[i]
    tx     <- aggregate_global_trends(agg_results[[var_type]][[sce]],weight=Dtweight, decades = DECADES)
    tx$driver <- drv
    trtF <- rbind(trtF, tx)
  }

  # Rename columns to consistent names
  names(trtF)[c(2, 4, 5, 6, 7, 8)] <- c("changeC", "med", "cq1", "cq2", "w1", "w2")
  trtF$year <- rep(DECADES, nrow(driver_meta))

  colorn <- c("WU" = "limegreen", "Res" = "tomato4", "LUC" = "orange", "Clim" = "royalblue")
  colorz <- c("Clim" = "dodgerblue4", "LUC" = "gold4", "Res" = "firebrick4", "WU" = "olivedrab")

  # --- 6.1  Build aggregation boxplot-in-time plot -------------------------

  fbox <- ggplot() +
    geom_linerange(data = trtF,
                   aes(x = year, ymin = w1, ymax = w2,
                       color = factor(driver), group = factor(driver)),
                   position = position_dodge2(width = 9), lwd = 1, alpha = 0.6) +
    scale_color_manual(values = colorn, name = "Drivers", labels = DRIVER_LABELS) +
    new_scale_color() +
    geom_rect(data = trtF,
              aes(xmin = year - 4.5, xmax = year + 4.5,
                  ymin = cq1, ymax = cq2,
                  fill = factor(driver), group = factor(driver)),
              alpha = 0.5, position = position_dodge(width = 9)) +
    geom_point(data = trtF,
               aes(x = year, y = changeC,
                   color = factor(driver), group = factor(driver)),
               position = position_dodge(width = 9), size = 3) +
    geom_rect(data = trtF,
              aes(xmin = year - 4.5, xmax = year + 4.5,
                  ymin = med - 1e-3, ymax = med + 1e-3,
                  fill = factor(driver), group = factor(driver)),
              alpha = 1, position = position_dodge(width = 9)) +
    scale_y_continuous(name = nplot, breaks = br_p,
                       trans = scales::modulus_trans(0.6)) +
    scale_x_continuous(name = "Decade",
                       breaks = DECADES, labels = DECADES,
                       minor_breaks = DECADES, expand = c(0.01, 0.01)) +
    scale_fill_manual(values = colorn, name = "Drivers", labels = DRIVER_LABELS) +
    scale_color_manual(values = colorz, name = "Drivers", labels = DRIVER_LABELS) +
    guides(color = guide_legend(override.aes = list(color = colorz))) +
    ggtitle(paste0("Europe – ", var_type)) +
    theme(
      axis.title       = element_text(size = 18, face = "bold"),
      title            = element_text(size = 22, face = "bold"),
      axis.text        = element_text(size = 16),
      axis.text.x      = element_text(size = 16, face = "bold"),
      panel.background = element_rect(fill = "white", colour = "white"),
      panel.border     = element_rect(linetype = "solid", fill = NA, colour = "black"),
      panel.grid       = element_blank(),
      panel.grid.major.y = element_line(color = "lightgray", linetype = "dashed"),
      panel.grid.minor.x = element_line(colour = "grey23",   linetype = "dashed"),
      legend.title     = element_text(size = 20, face = "bold"),
      legend.text      = element_text(size = 16),
      legend.position  = "right",
      legend.key       = element_rect(fill = "transparent", colour = "transparent"),
      legend.key.size  = unit(0.8, "cm"),
      axis.ticks.y     = element_blank()
    )

  # --- 6.2  Save aggregation plot ------------------------------------------

  agg_file <- paste0(plotDir, "bxp_EU_", haz, "_", name, "_", var_type, ".jpg")
  ggsave(agg_file, fbox, width = 30, height = 20, units = "cm", dpi = 1000)
  message("  Aggregation plot saved: ", agg_file)
}


# ==============================================================================
# 7. SAVE RL100mat  (return period = 10 yr) FOR ALL SCENARIOS
#
#  Scenarios and their sigma / threshold column mapping in df_Main:
#    SCF    -> sigmaGPD    / thresholdGPD   (base columns)
#    RWStat -> SigmaRW     / thresholdRW
#    WStat  -> SigmaW      / thresholdW
#    Histo  -> SigmaH      / thresholdH
# ==============================================================================

message("\n=== Saving RL100mat (RP = ", RETURN_PERIOD, " yr) for all scenarios ===")

rl_scenarios <- list(
  SCF    = list(sigma_col = "sigmaGPD",    threshold_col = "thresholdGPD"),
  RWStat = list(sigma_col = "SigmaRW",     threshold_col = "thresholdRW"),
  WStat  = list(sigma_col = "SigmaW",      threshold_col = "thresholdW"),
  Histo  = list(sigma_col = "SigmaH",      threshold_col = "thresholdH")
)

# Pre-split df_Main by year once (shared across all scenarios)
params_split   <- split(df_Main, df_Main$Year)
Years          <- names(params_split)
num_catchments <- length(unique(df_Main$catchment))

for (s_name in names(rl_scenarios)) {
  sigma_col     <- rl_scenarios[[s_name]]$sigma_col
  threshold_col <- rl_scenarios[[s_name]]$threshold_col

  cat("\n--- RL100mat for scenario:", s_name,
      " | sigma:", sigma_col, "| threshold:", threshold_col, "---\n")

  # Pre-allocate result matrix
  RL100mat <- matrix(NA_real_, nrow = num_catchments, ncol = length(Years))
  colnames(RL100mat) <- Years

  for (i in seq_along(Years)) {
    yr_label <- Years[i]
    cat("  Year:", yr_label, "\r")
    df_year <- params_split[[yr_label]]

    RL100mat[, i] <- calcGPDReturnLevel_Single(
      epsilon           = df_year$epsilonGPD,
      sigma             = df_year[[sigma_col]],
      threshold         = df_year[[threshold_col]],
      nPeaks            = df_year$nPeaks,
      sampleTimeHorizon = 70,
      returnPeriod      = RETURN_PERIOD
    )
  }

  # colnames(RL100mat)[1] <- "1951"
  RL100mat <- cbind(RL100mat, df_year$catchment)
  colnames(RL100mat) <- paste0("Y", colnames(RL100mat))
  colnames(RL100mat)[71] <- "unikout"
  RLmat <- data.frame(RL100mat)
  print(RLmat[1,])
  # Build output filename following existing convention
  out_name <- if (hazard == "Drought") {
    paste0("Drought.nonfrost.", s_name)
  } else {
    paste0("Flood.year.", s_name)
  }

  out_file <- paste0(hydroDir, "/", hazard, "/RL100xx.", out_name, ".Rdata")
  save(RLmat, file = out_file)
  cat("\n  Saved:", out_file, "\n")
}

message("\n=== All done ===")
