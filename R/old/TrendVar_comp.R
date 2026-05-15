setwd("D:/tilloal/Documents/LFRuns_utils/TrendAnalysis/R/")
source("functions_trends.R")
hydroDir<-("D:/tilloal/Documents/LFRuns_utils/data")


rspace= read.csv(paste0(hydroDir,"/subspace_efas.csv"))
rspace=rspace[,-1]
hazard="Flood"
if (hazard=="Flood"){
  dataDir="D:/tilloal/Documents/LFRuns_utils/data/Flood/HPC/Calibrated/revision/TrendVar/"
  }

if (hazard=="Drought"){
  dataDir="D:/tilloal/Documents/LFRuns_utils/data/Drought/HPC/Calibrated/revision/TrendVar/"
}
#outlets file outf
hydroDir<-("D:/tilloal/Documents/LFRuns_utils/ChangingHydroExtremes/data")
if (!exists("outf")){
  outf=c()
  for( Nsq in 1:88){
    print(Nsq)
    nrspace=rspace[Nsq,]
    outletname="GeoData/efas_rnet_100km_01min"
    
    outhybas=outletopen(hydroDir,outletname,nrspace)
    Idstart=as.numeric(Nsq)*10000
    Idstart2=as.numeric(Nsq)*100000
    if (length(outhybas$outlets)>0){
      outhybas$outlets=seq((Idstart+1),(Idstart+length(outhybas$outlets)))
      outhybas$outl2=seq((Idstart2+1),(Idstart2+length(outhybas$outlets)))
      outhybas$latlong=paste(round(outhybas$Var1,4),round(outhybas$Var2,4),sep=" ")
      outhloc=outhybas
      outf=rbind(outf,outhloc)
    }
  }
}

##2.1 Spatial data for catchments ----
#


### Hybas07 ----
Catchmentrivers7=read.csv(paste0(hydroDir,"/GeoData/HYBAS07/from_hybas_eu_onlyid.csv"),encoding = "UTF-8", header = T, stringsAsFactors = F)
hybas07 <- read_sf(dsn = paste0(hydroDir,"/GeoData/HYBAS07/hybas_eu_lev07_v1c.shp"))
hybasf7=fortify(hybas07)
Catamere07=inner_join(hybasf7,Catchmentrivers7,by= "HYBAS_ID")
Catamere07$llcoord=paste(round(Catamere07$POINT_X,4),round(Catamere07$POINT_Y,4),sep=" ")
cst7=right_join(Catamere07,outf,by= c("llcoord"="latlong"))
GNF=cst7
length(unique(GNF$HYBAS_ID))
st_geometry(GNF)=NULL
rm(Catamere07)
outlethybas07="/GeoData/HYBAS07/outletsv8_hybas07_01min"
outhybas07=outletopen(hydroDir,outlethybas07)

#matching outlets with pixel Ids
outhybas07$latlong=paste(round(outhybas07$Var1,4),round(outhybas07$Var2,4),sep=" ")
mhy=match(outhybas07$latlong,outf$latlong)
outhybas07$outID=outf$outl2[mhy]

### European Biogeo regions ----
biogeo <- read_sf(dsn = paste0(hydroDir,"/GeoData/eea_3035_biogeo-regions_2016/BiogeoRegions2016_wag84.shp"))
biogeof=fortify(biogeo)
st_geometry(biogeof)<-NULL
biogeoregions=raster( paste0(hydroDir,"/GeoData/eea_3035_biogeo-regions_2016/Biogeo_rasterized_wsg84.tif"))
Gbiogeoregions=as.data.frame(biogeoregions,xy=T)
biogeomatch=inner_join(biogeof,Gbiogeoregions,by= c("PK_UID"="Biogeo_rasterized_wsg84"))
biogeomatch$latlong=paste(round(biogeomatch$x,4),round(biogeomatch$y,4),sep=" ")
biogeo_rivers=right_join(biogeomatch,outf, by="latlong")

### HydroRegions ----

GridHR=raster( paste0(hydroDir,"/GeoData/HER/HydroRegions_raster_WGS84.tif"))
GHR=as.data.frame(GridHR,xy=T)
GHR2=GHR
GHR2$llcoord=paste(round(GHR2$x,4),round(GHR2$y,4),sep=" ")
GHR=GHR[which(!is.na(GHR[,3])),]
GHR$llcoord=paste(round(GHR$x,4),round(GHR$y,4),sep=" ")
GHR_riv=inner_join(GHR,outf,by= c("llcoord"="latlong"))
GHshpp <- read_sf(dsn =paste0(hydroDir,"/GeoData/HER/her_all_adjusted.shp"))
HydroRsf=fortify(GHshpp)

### Plot parameters ----
palet2=c(hcl.colors(9, palette = "Blues", alpha = NULL, rev = TRUE, fixup = TRUE))
outletname="/GeoData/efas_rnet_100km_01min"
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
rm(cst7)
basemap=w2

##2.2 Loading saved results in .Rdata ---------------------------

###load UpArea -----
#load upstream area
# main_path = 'D:/tilloal/Documents/06_Floodrivers/'
# valid_path = paste0(main_path,'DataPaper/')
outletname="/GeoData/upArea_European_01min.nc"
#dir=valid_path
outf$idlalo=paste(outf$idlo, outf$idla, sep=" ")
UpArea=UpAopen(hydroDir,outletname,outf)
head(UpArea)

#now load parameters


# load(file=paste0(dataDir,"/SCF/Var_Trend_agg.Rdata"))
# ###load Socio-CF run -----
# hydroDir<-("D:/tilloal/Documents/LFRuns_utils/ChangingHydroExtremes/data")
# haz="Flood"
# if (haz == "Drought") namefile="Drought.nonfrost.SocCF"
# if (haz == "Flood") namefile="Flood.year.socCF"
# 
# load(file=paste0(hydroDir,"/",haz,"/params.",namefile,".Rdata"))
# 
# Paramsfl=Paramsfl[,-c(4:9,17)]
# ParamsflSCF=data.table(Paramsfl)
# rm(Paramsfl)
# gc()

# rm(ParamVarTrend,ParamVarTrendHx,ParamVarTrendH,catmap)
# gc()

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
                                fileVar="Var_Trend_agg.Rdata",
                                hydro_base_dir = "D:/tilloal/Documents/LFRuns_utils/ChangingHydroExtremes/data" ) {
  
  require(dplyr)
  require(tidyr)
  require(lubridate)
  require(data.table)
  
  # 1. Determine the filename suffix based on sub_dir
  # This handles cases like SCF -> SocCF, Histo -> Histo, WStat -> WStat
  file_suffix <- case_when(
    sub_dir == "SCFX"   ~ "SocCF2",
    sub_dir == "HistoX" ~ "Histo",
    sub_dir == "WStat" ~ "WCF",
    sub_dir == "RWStat" ~ "RWCF",
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

hazard="Flood"
# Process all four versions
df_SCF    <- process_hazard_data(sub_dir = "SCFX",hazard = hazard,
                                 base_data_dir = dataDir,fileVar="Var_TrendX_agg.Rdata" )

df_Histo  <- process_hazard_data(sub_dir = "HistoX",hazard = hazard,base_data_dir = dataDir,
                                 fileVar="Var_TrendX_agg.Rdata")
df_WStat  <- process_hazard_data(sub_dir = "WStat",hazard = hazard,base_data_dir = dataDir,
                                 fileVar="Var_TrendX_agg.Rdata")
df_RWStat <- process_hazard_data(sub_dir = "RWStat",hazard = hazard,base_data_dir = dataDir,
                                 fileVar="Var_TrendX_agg.Rdata")


# 2. Add Histo columns to SCF (Suffix "H")
df_Main <- attach_scenario_cols(df_SCF, df_RWStat, "RW")
rm(df_RWStat)
df_Main <- attach_scenario_cols(df_Main, df_WStat, "W")
rm(df_WStat)
df_Main <- attach_scenario_cols(df_Main, df_Histo, "H")
rm(df_Histo)
gc()

df_Main <- attach_scenario_cols(df_SCF, df_Histo, "H")

df_Main$sq=floor(df_Main$catchment/100000)
df_51=df_Main[which(df_Main$sq==42),]
df_51y=df_51[which(df_51$Year==2015),]
merde=df_51y[which(abs((df_51y$trendRW-df_51y$trendW))>100),]

tc=4200002
dfm_test=df_Main[which(df_Main$catchment==tc),]

dfscf_test=df_SCF[which(df_SCF$catchment==tc),]



#load parameters form previous method
haz="Drought"
if (haz == "Drought") namefile="Drought.nonfrost.Histo"
if (haz == "Flood") namefile="flood.year.Histo"

load(file=paste0(hydroDir,"/",haz,"/params.",namefile,".Rdata"))
gc()
Paramsfl=(Paramsfl[,-c(4:9,17)])
ParamsflH=Paramsfl
rm(Paramsfl)
gc()
###load Socio-CF run -----
if (haz == "Drought") namefile="Drought.nonfrost.SocCF"
if (haz == "Flood") namefile="Flood.year.socCF"

load(file=paste0(hydroDir,"/",haz,"/params.",namefile,".Rdata"))

Paramsfl=Paramsfl[,-c(4:9,17)]
ParamsflSCF=data.table(Paramsfl)
rm(Paramsfl)
gc()


ParamHl=ParamsflH[which(ParamsflH$catchment==tc),]
ParamSl=ParamsflSCF[which(ParamsflSCF$catchment==tc),]

plot(dfm_test$time, dfm_test$variability)
points(dfm_test$time,dfm_test$variabilityH,col=2)


plot(dfm_test$time, -dfm_test$trend)
points(dfm_test$time,-dfm_test$trendH,col=2)

#intrusion from other code
# plot(transformed_data$timeStamps, transformed_data$trendSeries,col=1)
# lines(trasfData$timeStamps, trasfData$trendSeries,col=2)



hist(df_Main$SigmaW[which(df_Main$Year==1995)]/df_Main$SigmaRW[which(df_Main$Year==1995)],breaks=100000, xlim=c(0.8,1.2))
mean((df_Main$thresholdW[which(df_Main$Year==1955)]-df_Main$thresholdRW[which(df_Main$Year==1955)]),na.rm=T)
#load result file for format
load(file=paste0("D:/tilloal/Documents/LFRuns_utils/data/TSEVA/output_plots/Flood_pixChange_RL100_v1.Rdata"))

analyze_trends <- function(param_df, 
                                   data_save_list, 
                                   var_type = "Sigma", 
                                   scenario = "SCF",
                                   base_year = 1955, 
                                   return_period = 10,
                                   years_to_agg = seq(1955, 2015, by = 10)) {
  require(dplyr)
  require(tidyr)
  
  # 1. Map Scenario to Suffix
  # SCF is the base (no suffix), others use H, W, RW
  suffix <- case_when(
    scenario == "Histo"  ~ "H",
    scenario == "WStat"  ~ "W",
    scenario == "RWStat" ~ "RW",
    TRUE                 ~ ""    # Default for SCF
  )
  
  # 2. Identify the correct column name based on var_type and scenario
  # Example: var_type "Sigma" + scenario "Histo" -> "SigmaH"
  active_var <- paste0(var_type, suffix)
  
  # 3. Handle Special Case: RP (Return Period)
  # If user wants RP, we calculate it using the scenario-specific params
  if (var_type == "RP") {
    param_df$RP_val <- calcGPDReturnLevel_Single(
      epsilon = param_df$shapeGPD, 
      sigma = param_df[[paste0("Sigma", suffix)]], 
      threshold = param_df[[paste0("threshold", suffix)]], 
      nPeaks = param_df$nPeaks, 
      sampleTimeHorizon = param_df$sampleTimeHorizon, 
      returnPeriod = return_period
    )
    active_var <- "RP_val"
  }

  # 4. Pivot Data Wide
  df_wide <- param_df %>%
    dplyr::select(catchment, Year, !!sym(active_var)) %>%
    pivot_wider(names_from = Year, values_from = !!sym(active_var)) %>%
    as.data.frame()
  
  data_cols <- setdiff(names(df_wide), "catchment")
  base_year_col <- as.character(base_year)
  
  # 5. Variable-Specific Transformations (DataV)
  if (var_type == "threshold") {
    # Relative change (%) vs base year, normalized by row mean
    row_means <- rowMeans(df_wide[, data_cols], na.rm = TRUE)
    vals <- (df_wide[, data_cols] - df_wide[[base_year_col]]) / row_means * 100
    DataV <- data.frame(outl2 = df_wide$catchment, vals)
    
  } else if (var_type == "Sigma") {
    # Ratio vs base year
    vals <- df_wide[, data_cols] / df_wide[[base_year_col]]
    DataV <- data.frame(outl2 = df_wide$catchment, vals)
    
  } else {
    # For variability, trend, and RP
    DataV <- data.frame(outl2 = df_wide$catchment, df_wide[, data_cols])
  }
  
  # 6. Join with Metadata (DataO)
  DataO <- data_save_list[[1]][, -c(12:81)]
  DataV <- inner_join(DataO, DataV, by = "outl2")
  
  # 7. Aggregate Trends (trendClim)
  yr_cols <- paste0("X", years_to_agg)
  valid_cols <- intersect(yr_cols, names(DataV))
  if(length(valid_cols) == 0) valid_cols <- as.character(years_to_agg)
  
  trendagg <- aggregate(DataV[, valid_cols], 
                        by = list(HydroR = DataV$HER), 
                        FUN = function(x) mean(x, na.rm = TRUE))
  trend <- do.call(data.frame, trendagg)
  
  # 8. Point Statistics for 2015 (pointClim)
  col_2015 <- if("X2015" %in% names(DataV)) "X2015" else "2015"
  
  pointagg <- aggregate(list(Rchange_rel = DataV[[col_2015]]), 
                        by = list(HydroR = DataV$HER), 
                        FUN = function(x) {
                          c(mean = mean(x, na.rm = TRUE), 
                            dev = sd(x, na.rm = TRUE), 
                            len = length(x), 
                            med = median(x, na.rm = TRUE), 
                            q1 = quantile(x, 0.025, na.rm = TRUE), 
                            q3 = quantile(x, 0.975, na.rm = TRUE))
                        })
  point <- do.call(data.frame, pointagg)
  
  return(list(
    DataV = DataV,
    point = point,
    trend = trend,
    variable_plotted = active_var
  ))
}



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
        param_df$thresholdGPD=-param_df$thresholdGPD
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
        
    } else if (var_type == "Sigma") {
      vals <- df_wide[, data_cols] / row_init
    } else {
      vals <- df_wide[, data_cols]
    }
    DataV_vals <- vals
  } else {
    # Driver results are already symmetric percentages from step 2
    DataV_vals <- df_wide[, data_cols]
  }
  
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


#histogram of n peaks


# 1. Look at Histo Sigma

driver = "climate"
haz="flood"

var_type = "threshold"
res_driver <- analyze_hazard_symmetric(
  param_df = df_Main, 
  data_save_list = DataSave, 
  var_type = var_type, 
  scenario = "SCF",
  haz= haz,
  return_period=10
)




peakH=df_Main[which(df_Main$Year==2015),]


  psds<-ggplot(peakH, aes(x=nPeaks/70)) + 
  geom_histogram(color="steelblue", fill="slategray1",bins=200,alpha=0.9,lwd=1)+
  scale_y_continuous(breaks=seq(0,200000, by=10000),name="Number of pixels")+
  scale_x_continuous(breaks=seq(1,3, by=.2),limits=c(0,3))+
  # scale_x_sqrt(name=expression(paste("Upstream area ", (km^2),sep = " ")),
  #               breaks=c(100,1000,10000,100000), minor_breaks = log10_minor_break(),
  #               labels=c("100","1 000","10 000","100 000")) +
  theme(axis.title=element_text(size=16, face="bold"),
        axis.text = element_text(size=16),
        panel.background = element_rect(fill = "white", colour = "white"),
        panel.grid = element_blank(),
        panel.border = element_rect(linetype = "solid", fill = NA, colour="black"),
        legend.title = element_text(size=14),
        legend.text = element_text(size=12),
        panel.grid.major = element_line(colour = "grey80"),
        panel.grid.minor.x = element_line(colour = "grey90",linetype="dashed"),
        legend.key = element_rect(fill = "transparent", colour = "transparent"),
        legend.key.size = unit(.8, "cm"))







craps=res_driver$DataV[-which(is.na(res_driver$DataV$X2015)),]
points <- st_as_sf(res_driver$DataV[-which(is.na(res_driver$DataV$Var1)),], coords = c("Var1", "Var2"), crs = 4326)
#points=points[-which(is.na(points$X1955)),]
points <- st_transform(points, crs = 3035)

min(points$X2015,na.rm=T)
Regio=HydroRsf
pag <- inner_join(Regio, res_driver$point, by = c("CODEB" = "HydroR"))

if (var_type=="Sigma"){
  br=seq(0.5,1.5,by=0.1)
  labels=br
  limi=c(0.5,1.5)
  legend2="Sigma ratio"
}else if (var_type=="RL"){
  br=seq(-50,50,by=05)
  labels=br
  limi=c(-50,50)
  legend2="RL change"
}else if (var_type=="threshold"){
  br=seq(-50,50,by=05)
  labels=br
  limi=c(-50,50)
  legend2="Threshold change"
}

# br=seq(-1,1,by=.1)
# labels=br
# limi=c(-1,1)

tsize=16
osize=12

palet=c(hcl.colors(11, palette = "RdYlBu", alpha = NULL, rev = F, fixup = TRUE))
paletf=c(hcl.colors(11, palette = "RdBu", alpha = NULL, rev = F, fixup = TRUE))






####[Plot] - Figure 2b - Map of changes in flood 10-Y RL driven by climatic changes ----

#titleX=paste0("Change in 10-year ",haz," attributed \nto ",driver," changes (% of  100y flood) -  1955-2015")
fmap<-ggplot(basemap) +
  geom_sf(fill="white",color="darkgrey",size=0.5)+
  geom_sf(data=pag,aes(fill=Rchange_rel.mean,geometry=geometry),alpha=0.2,color="transparent")+
  geom_sf(data=points,aes(col=X2015,geometry=geometry,size=upa),alpha=.9,stroke=0,shape=15)+ 
  
  scale_size(range = c(0.08, 0.4), trans="sqrt",name= expression(paste("Upstream area ", (km^2),
                                                                       sep = " ")),
             breaks=c(101,1000,10000,100000,500000), labels=c("100","1000", "10 000", "100 000", "500 000"),
             guide = "none")+
  scale_fill_gradientn(
    colors=paletf,
    breaks=br,limits=limi,trans=scales::modulus_trans(.3),
    oob = scales::squish, name=legend2)   +
  guides(colour = guide_colourbar(barwidth = 22, barheight = 1), fill = "none")+
  # new_scale_fill()+
  # geom_sf(data=pointsInside,aes(geometry=geometry, fill=sign),alpha=.6,size=.5,stroke=0,shape=21,color="black")+
  # scale_fill_manual(values=c("tomato4","steelblue4"), name="Significant trends")+
  coord_sf(xlim = c(min(nco[,1]),max(nco[,1])), ylim = c(min(nco[,2]),max(nco[,2])))+
  scale_color_gradientn(
    colors=palet,
    breaks=br,limits=limi,trans=scales::modulus_trans(.3),
    oob = scales::squish, name=legend2)   +
  labs(x="Longitude", y = "Latitude")+
  # guides(colour = guide_colourbar(barwidth = 22, barheight = 1),
  #        fill = guide_legend(override.aes = list(size = 10)))+
  theme(axis.title=element_text(size=tsize),
        title = element_text(size=osize),
        axis.text=element_text(size=osize),
        panel.background = element_rect(fill = "aliceblue", colour = "grey1"),
        panel.border = element_rect(linetype = "solid", fill = NA, colour="black"),
        legend.title = element_text(size=tsize),
        legend.text = element_text(size=osize),
        legend.position = "bottom",
        legend.box = "vertical",  # Stack legends vertically
        panel.grid.major = element_line(colour = "grey70"),
        panel.grid.minor = element_line(colour = "grey90"),
        legend.key = element_rect(fill = "transparent", colour = "transparent"),
        legend.key.size = unit(1, "cm"))

ggsave(paste0("D:/tilloal/Documents/LFRuns_utils/TrendAnalysis/plots/map_",var_type,"_",driver,"_",haz,"X.jpg"), fmap,
       width=22, height=20, units=c("cm"),dpi=800) 




mp=match(points$outl2,peakH$catchment)
points$npeaks=peakH$nPeaks[mp]
palet=c(hcl.colors(11, palette = "RdYlBu", alpha = NULL, rev = T, fixup = TRUE))
cazzo=ggplot(basemap) +
  geom_sf(fill="gray95",color="gray10",size=0.5)+
  geom_sf(data=points,aes(col=npeaks,geometry=geometry,size=upa),alpha=.9,stroke=0,shape=15)+
  coord_sf(xlim = c(min(nco[,1]),max(nco[,1])), ylim = c(min(nco[,2]),max(nco[,2])))+
  # scale_colour_manual(values = colIR, name="IR", labels=c("0"="perennial","1"="casi-perennial","2"="IRES")) +
  scale_size(range = c(0.08, 0.4), trans="sqrt",name= expression(paste("Upstream area ", (km^2),
                                                                       sep = " ")),
             breaks=c(101,1000,10000,100000,500000), labels=c("100","1000", "10 000", "100 000", "500 000"),
             guide = "none")+
  scale_color_gradientn(
    colors=palet,
    breaks=seq(0,200,by=2), limits=c(0,200),
    oob = scales::squish)   +
  labs(x="Longitude", y = "Latitude")+
  # guides(colour = guide_legend(override.aes = list(size = 10)))+
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

ggsave(cazzo,filename=paste0("D:/tilloal/Documents/LFRuns_utils/TrendAnalysis/plots/maps_npeaks",haz,".jpg"), width=23, height=20, units=c("cm"),dpi=1000) 













aggregate_global_trends <- function(trend_df, 
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
  
  # 4. Global aggregation with statistics
  tr_global <- aggregate(
    list(value = tr_filtered$value),
    by = list(yr = tr_filtered$decad),
    FUN = function(x) {
      c(mean = mean(x, na.rm = TRUE),
        l    = length(x),
        med  = median(x, na.rm = TRUE),
        ql   = quantile(x, 0.25, na.rm = TRUE),
        qh   = quantile(x, 0.75, na.rm = TRUE),
        w1   = quantile(x, 0.025, na.rm = TRUE),
        w2   = quantile(x, 0.975, na.rm = TRUE))
    }
  )
  
  # 5. Convert matrix output to a clean data frame
  tDataHuman <- do.call(data.frame, tr_global)
  
  return(tDataHuman)
}

var_type = "Sigma"
if (var_type=="Sigma"){
  name="scale"
  br=c(seq(0,15,1))
  v="ratio"
}
if (var_type=="threshold"){
  name="location"
  br=c(seq(-100,100,5))
  v="relative difference"
}
if (var_type=="RP"){
  name="RL"
  br=c(seq(-100,100,5))
  v="relative difference"
}
lsce=c("SCF","RWStat","WStat","Histo")
trtF=c()
for (sce in lsce){
  print(paste0("scenario ",sce))
  res_driver <- analyze_hazard_symmetric(
    param_df = df_Main, 
    data_save_list = DataSave, 
    var_type = var_type, 
    scenario = sce,
    haz= haz,
    return_period=10
  )
  
  trend_df=res_driver$trendClim
  
  tx=aggregate_global_trends(trend_df)
  
  trtF<-rbind(trtF,tx)
}
bio_names=unique(biogeo$code)

#I keep only some bioregions
bio_names=bio_names[c(1,3,4,6,7,9,11)]

# combine trends from all drivers
#trtF=rbind(TdataClim$tGlobal,TdataLuse$tGlobal,TdataRes$tGlobal,TdataWuse$tGlobal)
trtF$year=rep(seq(1950,2010,10),4)
trtF$driver=c(rep("Clim",7),rep("LUC",7),rep("Res",7),rep("WU",7))

names(trtF)[c(2,4,5,6,7,8)]=c("changeC","med","cq1","cq2","w1","w2")

colorz = c("Clim" ='dodgerblue4',"LUC" ='gold4',"Res" ='firebrick4',"WU"="olivedrab")
colorn = c("WU" ='limegreen',"Res" ='tomato4',"LUC" ='orange',"Clim" ='royalblue')


#### [Plot] - Figure 1 - Boxplot in time ----

fac=1
xlabs=seq(1950,2010,10)
clabels=c("Climate","Land use","Reservoirs", "Water demand")
nplot=paste0(name," parameters ",v)

ggplot() +
  # IQR represented as rectangles
  geom_linerange(data=trtF,aes(x=year, ymin=fac*(w1),ymax=fac*w2,color = factor(driver),group=factor(driver)),
                 position = position_dodge2(width = 9),lwd=1,alpha=0.6) +
  scale_color_manual(values = colorn, name = "Drivers", labels = clabels) +
  new_scale_color()+
  
  geom_rect(data=trtF, aes(xmin = year - 4.5, xmax = year + 4.5, 
                           ymin = cq1, ymax = cq2, fill = factor(driver), group = factor(driver)), 
            alpha = 0.5, position = position_dodge(width = 9)) +
  
  # Mean as points over the IQR
  geom_point(data=trtF, aes(x = year, y = changeC, color = factor(driver), group = factor(driver)),
             position = position_dodge(width = 9), size = 3) +
  
  # Median as a horizontal line within the IQR rectangle
  geom_rect(data=trtF, aes(xmin = year - 4.5, xmax = year + 4.5, 
                           ymin = med-1e-3, ymax = med+1e-3, fill = factor(driver), group = factor(driver)), 
            alpha = 1, position = position_dodge(width = 9)) +
  
  # Y-axis settings
  scale_y_continuous(name = nplot, breaks = br, trans=scales::modulus_trans(.6)) +
  
  # X-axis settings
  scale_x_continuous(breaks = xlabs, labels = xlabs, name = "Decades",
                     minor_breaks = seq(1955,2005,10), expand = c(.01,0.01)) +
  
  # Manual fill and color scales
  scale_fill_manual(values = colorn, name = "Drivers", labels = clabels) +
  scale_color_manual(values = colorz, name = "Drivers", labels = clabels) +
  
  # Customize legend
  guides(color = guide_legend(override.aes = list(color = colorz))) +
  
  # Theme customization
  theme(
    axis.title = element_text(size = 18, face = "bold"),
    title = element_text(size = 22, face = "bold"),
    axis.text = element_text(size = 16),
    axis.text.x = element_text(size = 16, face = "bold"),
    panel.background = element_rect(fill = "white", colour = "white"),
    panel.grid = element_blank(),
    panel.border = element_rect(linetype = "solid", fill = NA, colour = "black"),
    legend.title = element_text(size = 20, face = "bold"),
    legend.text = element_text(size = 16),
    axis.ticks.y = element_blank(),
    panel.grid.major.y = element_line(color = "lightgray",linetype = "dashed"),
    #panel.grid.minor.y = element_line(color = "lightgray"),
    legend.position = "right",
    #panel.grid.major = element_line(colour = "grey80"),
    panel.grid.minor.x = element_line(colour = "grey23",linetype = "dashed"),
    legend.key = element_rect(fill = "transparent", colour = "transparent"),
    legend.key.size = unit(0.8, "cm")
  ) +
  
  # Add title
  ggtitle("Europe")

ggsave(paste0("D:/tilloal/Documents/LFRuns_utils/TrendAnalysis/plots/bxp_in_time_EU_",haz,"_",name,"_",var_type,"X.jpg"), width=30, height=20, units=c("cm"),dpi=1000) 

# yrlist=c(1951:2020)
# driverlist=c("climate","landuse","reservoirs","wateruse","all")
# driver=driverlist[2]
# Dchangelist=list()
# for (driver in driverlist){
#   print(driver)
#   if (driver=="climate"){
#     trendPlot=trendClim
#     datap=DataC
#     pointagg=pointClim
#   }
#   if (driver=="landuse"){
#     trendPlot=trendSoc
#     datap=DataL
#     pointagg=pointSoc
#   }
#   if (driver=="reservoirs"){
#     trendPlot=trendRes
#     datap=DataR
#     datap$Y2015[match(Rcrap$outl2,datap$outl2)]=0
#     pointagg=pointRes
#   }
#   if (driver=="wateruse"){
#     trendPlot=trendWu
#     datap=DataW
#     pointagg=pointWu
#   }
#   if (driver=="all"){
#     trendPlot=trendTot
#     datap=DataT
#     pointagg=pointTot
#   }
#   Pplot=calculatePoints(trendPlot, yrlist, pointagg, Regio, GHshpp, datap)
#   
#   save=T
#   
#   if (save==T){ 
#     colNA="transparent"
#     
#     if (driver=="climate"){
#       ####[Plot] - supplement- ordered change aggregated at the HER level ----
#       pointP=Pplot$PagD
#       uhi=unique(pointP$CODEB)
#       pointP=pointP[match(uhi,pointP$CODEB),]
#       pointP=pointP[order(pointP$Rchange_rel.mean),]
#       pointP$id=c(1:length(pointP$Rchange_rel.mean))
#       
#       hist(pointP$Rchange_rel.len,breaks=100)
#       hist(GHshpp$SURF_KM2,breaks=100)
#       mean(GHshpp$SURF_KM2)
#       median(pointP$Rchange_rel.len)
#       print(length(which(pointP$change==-2)))/length(pointP$Id)
#       print(length(which(pointP$change==2)))/length(pointP$Id)
#       manualcol=c("-2"="#A51122","-1"="#F1C363", "1"="#ACD2BB","2"= "#324DA0")
#       manualab=c("sig. decrease","decrease", "increase","sig. increase")
#       brl=c(-200,200)
#       by=50
#       
#       br=c(brl[1],-100,-50,-10,0,10,50,100,brl[2])
#       limi=c(-5000,5000)
#       
#       ggplot() +
#         coord_cartesian(ylim=c(brl[1],brl[2]))+
#         geom_hline(yintercept = 0,lwd=1, col="black")+
#         geom_segment(data=pointP,aes(x=id, xend=id, y=Rchange_rel.q1.2.5.,
#                                      yend=Rchange_rel.q3.97.5.,color=factor(change),size = Rchange_rel.len), alpha=.99)+
#         scale_color_manual(values = manualcol,breaks=manualab,labels=manualab,name="")+
#         geom_point(data=pointP, aes(x=id, y=Rchange_rel.mean), pch=21, fill="white",colour="gray3",size=3,stroke=1,alpha=1) + 
#         geom_text(data=pointP, aes(x=id, y=Rchange_rel.mean,label = CODEB), size=1.5, color = "black",fontface = "bold") +
#         scale_size(range = c(0.8, 4),trans="sqrt",
#                    guide = "none")+
#         scale_y_continuous(limits=limi,breaks=br,name="Change (%)",trans=scales::modulus_trans(.3))+
#         scale_x_continuous(name="HER",expand=c(.01,.01),breaks=c(-100,200))+
#         guides(colour = guide_legend(override.aes = list(size = 10)))+
#         theme(axis.title=element_text(size=20, face="bold",color="black"),
#               axis.text = element_text(size=18,color="black"),
#               panel.background = element_rect(fill = "white", colour = "white"),
#               panel.grid = element_blank(),
#               panel.border = element_rect(linetype = "solid", fill = NA, colour="black",linewidth=2),
#               legend.title = element_text(size=20),
#               legend.text = element_text(size=18,color="black"),
#               legend.position = "none",
#               panel.grid.major = element_line(colour = "transparent"),
#               panel.grid.minor.y = element_line(colour = "transparent",linetype="dashed"),
#               legend.key = element_rect(fill = "transparent", colour = "transparent"),
#               legend.key.size = unit(.8, "cm"))
#       
#       ggsave(paste0("D:/tilloal/Documents/LFRuns_utils/TrendAnalysis/plots/ordered_HR",it,"_",haz,".jpg"),width=40, height=8, units=c("cm"),dpi=400) 
#       
#     }
#     
#     #differentiated plot schemes for flood and drought
#     if (haz=="Flood"){
#       br=c(-50,-20,-10,-5,0,5,10,20,50)
#       labels=br
#       limi=c(-50,50)
#       tsize=16
#       osize=12
#       legend2="Change (%)    "
#       palet=c(hcl.colors(11, palette = "RdYlBu", alpha = NULL, rev = F, fixup = TRUE))
#       paletf=c(hcl.colors(11, palette = "RdBu", alpha = NULL, rev = F, fixup = TRUE))
#       points=Pplot$points
#       pag=Pplot$PagD
#       pointsInside=Pplot$psp
#       ####[Plot] - Figure 2b - Map of changes in flood 10-Y RL driven by climatic changes ----
#       if (driver=="climate" | driver=="all"){
#         
#         titleX=paste0("Change in 10-year ",haz," attributed \nto ",driver," changes (% of  10y flood) -  1955-2015")
#         fmap<-ggplot(basemap) +
#           geom_sf(fill="white",color="darkgrey",size=0.5)+
#           geom_sf(data=pag,aes(fill=Rchange_rel.mean,geometry=geometry),alpha=0.2,color="transparent")+
#           geom_sf(data=points,aes(col=Y2015,geometry=geometry,size=upa),alpha=.9,stroke=0,shape=15)+ 
#           
#           scale_size(range = c(0.08, 0.4), trans="sqrt",name= expression(paste("Upstream area ", (km^2),
#                                                                                sep = " ")),
#                      breaks=c(101,1000,10000,100000,500000), labels=c("100","1000", "10 000", "100 000", "500 000"),
#                      guide = "none")+
#           scale_fill_gradientn(
#             colors=paletf,
#             breaks=br,limits=limi,trans=scales::modulus_trans(.3),
#             oob = scales::squish,na.value=colNA, name=legend2)   +
#           guides(fill = "none")+
#           new_scale_fill()+
#           geom_sf(data=pointsInside,aes(geometry=geometry, fill=sign),alpha=.6,size=.5,stroke=0,shape=21,color="black")+
#           scale_fill_manual(values=c("tomato4","steelblue4"), name="Significant trends")+
#           coord_sf(xlim = c(min(nco[,1]),max(nco[,1])), ylim = c(min(nco[,2]),max(nco[,2])))+
#           scale_color_gradientn(
#             colors=palet,
#             breaks=br,limits=limi,trans=scales::modulus_trans(.3),
#             oob = scales::squish,na.value=colNA, name=legend2)   +
#           labs(x="Longitude", y = "Latitude")+
#           guides(colour = guide_colourbar(barwidth = 22, barheight = 1),
#                  fill = guide_legend(override.aes = list(size = 10)))+
#           theme(axis.title=element_text(size=tsize),
#                 title = element_text(size=osize),
#                 axis.text=element_text(size=osize),
#                 panel.background = element_rect(fill = "aliceblue", colour = "grey1"),
#                 panel.border = element_rect(linetype = "solid", fill = NA, colour="black"),
#                 legend.title = element_text(size=tsize),
#                 legend.text = element_text(size=osize),
#                 legend.position = "bottom",
#                 legend.box = "vertical",  # Stack legends vertically
#                 panel.grid.major = element_line(colour = "grey70"),
#                 panel.grid.minor = element_line(colour = "grey90"),
#                 legend.key = element_rect(fill = "transparent", colour = "transparent"),
#                 legend.key.size = unit(1, "cm"))+
#           ggtitle(titleX)
#         
#         
#         #frequency of significant change regions
#         ls=length(unique(pointsInside$Id))
#         lsp=length(unique(pointsInside$Id[which(pointsInside$sign=="increase")]))
#         lsn=length(unique(pointsInside$Id[which(pointsInside$sign=="decrease")]))
#         lx=length(pag$Id)
#         lsp/lx
#         lsn/lx
#         ggsave(paste0("D:/tilloal/Documents/LFRuns_utils/TrendAnalysis/plots/mapF_",driver,"_",haz,mmx,"rel_New2",it,".jpg"), fmap, width=22, height=20, units=c("cm"),dpi=1000) 
#         
#       }else{
#         ####[Plot] - Figure 3 - Map of changes in flood 10-Y RL driven by socioeconomic changes ----
#         titleX=paste0("Change in 10-year ",haz," attributed \nto ",driver," changes (% of 10y flood) -  1955-2015")
#         legend2="Change (%)"
#         ocrap<-ggplot(basemap) +
#           geom_sf(fill="white",color="darkgrey",size=0.5)+
#           geom_sf(data=pag,aes(fill=Rchange_rel.mean,geometry=geometry),alpha=0.6,color="transparent")+
#           geom_sf(data=points,aes(col=Y2015,geometry=geometry,size=upa),alpha=.9,stroke=0,shape=15)+ 
#           
#           scale_size(range = c(0.1, 0.5), trans="sqrt",name= expression(paste("Upstream area ", (km^2),
#                                                                               sep = " ")),
#                      breaks=c(101,1000,10000,100000,500000), labels=c("100","1000", "10 000", "100 000", "500 000"),
#                      guide = "none")+
#           scale_fill_gradientn(
#             colors=palet,
#             breaks=br,limits=limi,trans=scales::modulus_trans(.3),
#             oob = scales::squish,na.value=colNA, name="Change (%)")   +
#           coord_sf(xlim = c(min(nco[,1]),max(nco[,1])), ylim = c(min(nco[,2]),max(nco[,2])))+
#           scale_color_gradientn(
#             colors=palet,
#             breaks=br,limits=limi,trans=scales::modulus_trans(.3),
#             oob = scales::squish,na.value="transparent", name="Change (%)")   +
#           labs(x="Longitude", y = "Latitude")+
#           guides(colour = guide_colourbar(barwidth = 1.5, barheight = 14),
#                  fill = guide_colourbar(barwidth = 1.5, barheight = 14))+
#           theme(axis.title=element_text(size=tsize),
#                 title = element_text(size=osize),
#                 axis.text=element_text(size=osize),
#                 panel.background = element_rect(fill = "aliceblue", colour = "grey1"),
#                 panel.border = element_rect(linetype = "solid", fill = NA, colour="black"),
#                 legend.text = element_text(size=8),
#                 legend.title = element_text(size = osize, margin = margin(t = 2, r = 2, b = 6, l = 0)),
#                 legend.spacing.x = unit(0.2, "cm"),
#                 legend.position = "bottom",
#                 legend.box = "vertical",  # Stack legends vertically
#                 panel.grid.major = element_line(colour = "grey70"),
#                 panel.grid.minor = element_line(colour = "grey90"),
#                 legend.key = element_rect(fill = "transparent", colour = "transparent"),
#                 legend.key.size = unit(1, "cm"))+
#           ggtitle(titleX)
#         
#         ggsave(paste0("D:/tilloal/Documents/LFRuns_utils/TrendAnalysis/plots/maF_",driver,"_",haz,"_",it,".jpg"), ocrap, width=23, height=20, units=c("cm"),dpi=1000) 
#         
#       } 
#     }else if(haz=="Drought"){
#       br=c(-50,-20,-10,-5,0,5,10,20,50)
#       labels=br
#       limi=c(-50,50)
#       tsize=16
#       osize=12
#       legend2="Change (%)"
#       palet=c(hcl.colors(11, palette = "RdYlBu", alpha = NULL, rev = F, fixup = TRUE))
#       paletf=c(hcl.colors(11, palette = "RdYlBu", alpha = NULL, rev = F, fixup = TRUE))
#       points=Pplot$points
#       pag=Pplot$PagD
#       pointsInside=Pplot$psp
#       ####[Plot] - Figure 2a - Map of changes in drought 10-Y RL driven by climatic changes ----
#       if (driver=="climate" | driver=="all"){
#         #specifically designed plot for climate
#         titleX=paste0("Change in 10-year ",haz," attributed \n to ",driver, "changes (% of 10y drought) - 1955-2015")
#         points=points[-which(is.na(points$Y2015)),]
#         ocrap<-ggplot(basemap) +
#           geom_sf(fill="white",color="darkgrey",size=0.5)+
#           geom_sf(data=pag,aes(fill=Rchange_rel.mean,geometry=geometry),alpha=0.2,color="transparent")+
#           geom_sf(data=points,aes(col=Y2015,geometry=geometry,size=upa),alpha=.9,stroke=0,shape=15)+ 
#           
#           scale_size(range = c(0.08, 0.4), trans="sqrt",name= expression(paste("Upstream area ", (km^2),
#                                                                                sep = " ")),
#                      breaks=c(101,1000,10000,100000,500000), labels=c("100","1000", "10 000", "100 000", "500 000"),
#                      guide = "none")+
#           scale_fill_gradientn(
#             colors=paletf,
#             breaks=br,limits=limi,labels = labels, trans=scales::modulus_trans(.3),
#             oob = scales::squish,na.value=colNA, name=legend2)   +
#           guides(fill = "none")+
#           new_scale_fill()+
#           geom_sf(data=pointsInside,aes(geometry=geometry, fill=sign),alpha=.6,size=.5,stroke=0,shape=21,color="black")+
#           scale_fill_manual(values=c("tomato4","steelblue4"), name="Significant trends")+
#           coord_sf(xlim = c(min(nco[,1]),max(nco[,1])), ylim = c(min(nco[,2]),max(nco[,2])))+
#           scale_color_gradientn(
#             colors=palet,
#             breaks=br,limits=limi,labels = labels,trans=scales::modulus_trans(.3),
#             oob = scales::squish,na.value="transparent", name=legend2)   +
#           labs(x="Longitude", y = "Latitude")+
#           guides(colour = guide_colourbar(barwidth = 1.5, barheight = 14),
#                  fill = guide_legend(override.aes = list(size = 10)))+
#           theme(axis.title=element_text(size=tsize),
#                 title = element_text(size=osize),
#                 axis.text=element_text(size=osize),
#                 panel.background = element_rect(fill = "aliceblue", colour = "grey1"),
#                 panel.border = element_rect(linetype = "solid", fill = NA, colour="black"),
#                 legend.title = element_text(size=tsize),
#                 legend.text = element_text(size=osize),
#                 legend.position = "right",
#                 panel.grid.major = element_line(colour = "grey70"),
#                 panel.grid.minor = element_line(colour = "grey90"),
#                 legend.key = element_rect(fill = "transparent", colour = "transparent"),
#                 legend.key.size = unit(1, "cm"))+
#           ggtitle(titleX)
#         
#         ggsave(paste0("D:/tilloal/Documents/LFRuns_utils/TrendAnalysis/plots/mapF_",driver,"_",haz,"_",it,".jpg"), ocrap, width=22, height=20, units=c("cm"),dpi=1000) 
#         
#       }else{
#         ####[Plot] - Figure 3 - Map of changes in drought 10-Y RL driven by socioeconomic changes ----
#         titleX=paste0("Change in 10-year ",haz," attributed \nto ", driver," changes (% of  10y drought) - 1955-2015")
#         ocrap<-ggplot(basemap) +
#           geom_sf(fill="white",color="darkgrey",size=0.5)+
#           geom_sf(data=pag,aes(fill=Rchange_rel.mean,geometry=geometry),alpha=0.6,color="transparent")+
#           geom_sf(data=points,aes(col=Y2015,geometry=geometry,size=upa),alpha=.9,stroke=0,shape=15)+ 
#           
#           scale_size(range = c(0.1, 0.5), trans="sqrt",name= expression(paste("Upstream area ", (km^2),
#                                                                               sep = " ")),
#                      breaks=c(101,1000,10000,100000,500000), labels=c("100","1000", "10 000", "100 000", "500 000"),
#                      guide = "none")+
#           scale_fill_gradientn(
#             colors=paletf,
#             breaks=br,limits=limi,labels = labels,trans=scales::modulus_trans(.3),
#             oob = scales::squish,na.value=colNA, name="Change (%)")   +
#           coord_sf(xlim = c(min(nco[,1]),max(nco[,1])), ylim = c(min(nco[,2]),max(nco[,2])))+
#           scale_color_gradientn(
#             colors=paletf,
#             breaks=br,limits=limi,labels = labels, trans=scales::modulus_trans(.3),
#             oob = scales::squish,na.value=colNA, name="Change (%)")   +
#           labs(x="Longitude", y = "Latitude")+
#           guides(colour = guide_colourbar(barwidth = 1.5, barheight = 14),
#                  fill = guide_colourbar(barwidth = 1.5, barheight = 14))+
#           theme(axis.title=element_text(size=tsize),
#                 title = element_text(size=osize),
#                 axis.text=element_text(size=osize),
#                 panel.background = element_rect(fill = "aliceblue", colour = "grey1"),
#                 panel.border = element_rect(linetype = "solid", fill = NA, colour="black"),
#                 legend.text = element_text(size=8),
#                 legend.title = element_text(size = osize, margin = margin(t = 2, r = 2, b = 6, l = 0)),
#                 legend.spacing.x = unit(0.2, "cm"),
#                 legend.position = "right",
#                 legend.box = "horizontal",  # Stack legends vertically
#                 panel.grid.major = element_line(colour = "grey70"),
#                 panel.grid.minor = element_line(colour = "grey90"),
#                 legend.key = element_rect(fill = "transparent", colour = "transparent"),
#                 legend.key.size = unit(1, "cm"))+
#           ggtitle(titleX)
#         
#         ggsave(paste0("D:/tilloal/Documents/LFRuns_utils/TrendAnalysis/plots/mapF_",driver,"_",haz,"_",it,".jpg"), ocrap, width=23, height=20, units=c("cm"),dpi=1000) 
#         
#       }
#     }
#   }
#   Dchangelist=c(Dchangelist,list(Pplot))
# }


#NOW COMPUTE DIFFERENT rlS

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



# 1. Define the Suffixes and Scenarios
# SCF is your "Base" (Sigma0/thresholdGPD), others are the suffixes you attached
scenarios <- c("SCF" = "", "RW" = "RW", "W" = "W", "H" = "H")

# 2. Convert to data.table if not already for speed
setDT(df_Main)

# 3. Loop through each scenario suffix
for (s_name in names(scenarios)) {
  suffix <- scenarios[[s_name]]
  
  cat("\n--- Calculating Return Levels for Scenario:", s_name, "---\n")
  
  # Identify the correct column names for this scenario
  # Logic: If suffix is empty, use the base names, otherwise use SigmaRW, etc.
  sigma_col     <- if(suffix == "") "sigmaGPD"    else paste0("Sigma", suffix)
  threshold_col <- if(suffix == "") "thresholdGPD" else paste0("threshold", suffix)
  
  # 1. Pre-calculate groups by Year (Massive speed boost)
  params_split <- split(df_Main, df_Main$Year)
  Years <- names(params_split)
  
  # 2. Pre-allocate Result Matrix
  # Rows = unique catchments, Cols = Years
  num_catchments <- length(unique(df_Main$catchment))
  RL100mat <- matrix(NA_real_, nrow = num_catchments, ncol = length(Years))
  colnames(RL100mat) <- Years
  
  # 3. Yearly Loop
  for (i in seq_along(Years)) {
    year_label <- Years[i]
    cat("  Processing Year:", year_label, "\r")
    
    df_year <- params_split[[year_label]]
    
    # Vectorized Return Level Calculation
    RL100mat[, i] <- calcGPDReturnLevel_Single(
      epsilon = df_year$epsilonGPD,   # Epsilon usually stays constant across scenarios
      sigma   = df_year[[sigma_col]], 
      threshold = df_year[[threshold_col]],
      nPeaks  = df_year$nPeaks, 
      sampleTimeHorizon = 70, 
      returnPeriod = 10
    )
  }
  
  # 4. Save results using your namefile logic
  out_name <- if (hazard == "Drought") {
    paste0("Drought.nonfrost.", s_name)
  } else {
    paste0("Flood.year.", s_name)
  }
  
  save(RL100mat, file = paste0(hydroDir, "/", hazard, "/RL100x.", out_name, ".Rdata"))
  cat("\nFinished and Saved:", out_name, "\n")
}





###load historical run -----
haz="Flood"
for (haz in c("Flood","Drought")){
  
  if (haz == "Drought") namefile="Drought.nonfrost.Histo"
  if (haz == "Flood") namefile="flood.year.Histo"
  print(namefile)
  load(file=paste0(hydroDir,"/",haz,"/params.",namefile,".Rdata"))
  
  Paramsfl=(Paramsfl[,-c(4:9,17)])
  ParamsflH=Paramsfl
  rm(Paramsfl)
  gc()
  
  # 1. Pre-calculate the groups once (massive speed boost over 'which')
  # This creates a list where each element is the data for one year
  params_split <- split(ParamsflH, ParamsflH$Year)
  Years <- names(params_split)
  
  # 2. Pre-allocate the result matrix (prevents memory fragmentation)
  # Rows = number of spatial points/observations, Cols = number of years
  num_rows_per_year <- nrow(params_split[[1]]) 
  RL100mat <- matrix(NA_real_, nrow = num_rows_per_year, ncol = length(Years))
  colnames(RL100mat) <- Years
  
  # 3. The Optimized Loop
  for (i in seq_along(Years)) {
    year_label <- Years[i]
    cat("Processing Year:", year_label, "\n") # faster than print()
    
    # Extract the data for this year from our pre-split list
    df_year <- params_split[[year_label]]
    
    # Run the vectorized function (using the NA-safe version we discussed)
    # This calculates all 300k (or whatever your spatial count is) rows at once
    RL100s <- calcGPDReturnLevel_Single(
      epsilon = df_year$epsilonGPD, 
      sigma = df_year$sigmaGPD, 
      threshold = df_year$thresholdGPD,
      nPeaks = df_year$nPeaks, 
      sampleTimeHorizon = 70, 
      returnPeriod = 100
    )
    
    # Assign directly into the pre-allocated column (very fast)
    RL100mat[, i] <- RL100s
  }
  
  save(RL100mat,file=paste0(hydroDir,"/",haz,"/RL100x.",namefile,".Rdata"))
  
  
  ###load Socio-CF run -----
  if (haz == "Drought") namefile="Drought.nonfrost.SocCF"
  if (haz == "Flood") namefile="Flood.year.socCF"
  print(namefile)
  load(file=paste0(hydroDir,"/",haz,"/params.",namefile,".Rdata"))
  
  Paramsfl=Paramsfl[,-c(4:9,17)]
  ParamsflSCF=data.table(Paramsfl)
  rm(Paramsfl)
  gc()
  
  
  # 1. Pre-calculate the groups once (massive speed boost over 'which')
  # This creates a list where each element is the data for one year
  params_split <- split(ParamsflSCF, ParamsflH$Year)
  Years <- names(params_split)
  
  # 2. Pre-allocate the result matrix (prevents memory fragmentation)
  # Rows = number of spatial points/observations, Cols = number of years
  num_rows_per_year <- nrow(params_split[[1]]) 
  RL100mat <- matrix(NA_real_, nrow = num_rows_per_year, ncol = length(Years))
  colnames(RL100mat) <- Years
  
  # 3. The Optimized Loop
  for (i in seq_along(Years)) {
    year_label <- Years[i]
    cat("Processing Year:", year_label, "\n") # faster than print()
    
    # Extract the data for this year from our pre-split list
    df_year <- params_split[[year_label]]
    
    # Run the vectorized function (using the NA-safe version we discussed)
    # This calculates all 300k (or whatever your spatial count is) rows at once
    RL100s <- calcGPDReturnLevel_Single(
      epsilon = df_year$epsilonGPD, 
      sigma = df_year$sigmaGPD, 
      threshold = df_year$thresholdGPD,
      nPeaks = df_year$nPeaks, 
      sampleTimeHorizon = 70, 
      returnPeriod = 100
    )
    
    # Assign directly into the pre-allocated column (very fast)
    RL100mat[, i] <- RL100s
  }
  
  save(RL100mat,file=paste0(hydroDir,"/",haz,"/RL100x.",namefile,".Rdata"))
  
  ###load results from Res+WU CF run -----
  if (haz == "Drought") namefile="Drought.nonfrost.RWCF"
  if (haz == "Flood") namefile="flood.year.RWCF"
  print(namefile)
  load(file=paste0(hydroDir,"/",haz,"/params.",namefile,".Rdata"))
  Paramsfl=Paramsfl[,-c(4:9,17)]
  ParamsflRWCF=data.table(Paramsfl)
  rm(Paramsfl)
  gc()
  
  # 1. Pre-calculate the groups once (massive speed boost over 'which')
  # This creates a list where each element is the data for one year
  params_split <- split(ParamsflRWCF, ParamsflH$Year)
  Years <- names(params_split)
  
  # 2. Pre-allocate the result matrix (prevents memory fragmentation)
  # Rows = number of spatial points/observations, Cols = number of years
  num_rows_per_year <- nrow(params_split[[1]]) 
  RL100mat <- matrix(NA_real_, nrow = num_rows_per_year, ncol = length(Years))
  colnames(RL100mat) <- Years
  
  # 3. The Optimized Loop
  for (i in seq_along(Years)) {
    year_label <- Years[i]
    cat("Processing Year:", year_label, "\n") # faster than print()
    
    # Extract the data for this year from our pre-split list
    df_year <- params_split[[year_label]]
    
    # Run the vectorized function (using the NA-safe version we discussed)
    # This calculates all 300k (or whatever your spatial count is) rows at once
    RL100s <- calcGPDReturnLevel_Single(
      epsilon = df_year$epsilonGPD, 
      sigma = df_year$sigmaGPD, 
      threshold = df_year$thresholdGPD,
      nPeaks = df_year$nPeaks, 
      sampleTimeHorizon = 70, 
      returnPeriod = 100
    )
    
    # Assign directly into the pre-allocated column (very fast)
    RL100mat[, i] <- RL100s
  }
  
  save(RL100mat,file=paste0(hydroDir,"/",haz,"/RL100x.",namefile,".Rdata"))
  
  ###load results from Water CF run -----
  if (haz == "Drought") namefile="Drought.nonfrost.WCF"
  if (haz == "Flood") namefile="flood.year.WCF"
  print(namefile)
  load(file=paste0(hydroDir,"/",haz,"/params.",namefile,".Rdata"))
  
  
  Paramsfl=Paramsfl[,-c(4:9,17)]
  ParamsflWCF=data.table(Paramsfl)
  rm(Paramsfl)
  gc()
  
  # 1. Pre-calculate the groups once (massive speed boost over 'which')
  # This creates a list where each element is the data for one year
  params_split <- split(ParamsflWCF, ParamsflH$Year)
  Years <- names(params_split)
  
  # 2. Pre-allocate the result matrix (prevents memory fragmentation)
  # Rows = number of spatial points/observations, Cols = number of years
  num_rows_per_year <- nrow(params_split[[1]]) 
  RL100mat <- matrix(NA_real_, nrow = num_rows_per_year, ncol = length(Years))
  colnames(RL100mat) <- Years
  
  # 3. The Optimized Loop
  for (i in seq_along(Years)) {
    year_label <- Years[i]
    cat("Processing Year:", year_label, "\n") # faster than print()
    
    # Extract the data for this year from our pre-split list
    df_year <- params_split[[year_label]]
    
    # Run the vectorized function (using the NA-safe version we discussed)
    # This calculates all 300k (or whatever your spatial count is) rows at once
    RL100s <- calcGPDReturnLevel_Single(
      epsilon = df_year$epsilonGPD, 
      sigma = df_year$sigmaGPD, 
      threshold = df_year$thresholdGPD,
      nPeaks = df_year$nPeaks, 
      sampleTimeHorizon = 70, 
      returnPeriod = 100
    )
    
    # Assign directly into the pre-allocated column (very fast)
    RL100mat[, i] <- RL100s
  }
  
  save(RL100mat,file=paste0(hydroDir,"/",haz,"/RL100x.",namefile,".Rdata"))
}





# 2. Look at WStat Return Period (100yr)
res_wstat_rp <- analyze_trends(df_Main, DataSave, 
                                       var_type = "RP", 
                                       scenario = "WStat", 
                                       return_period = 10)

# 3. Look at original SCF variability
res_scf_var <- analyze_climate_trends(df_Main, DataSave, 
                                      var_type = "variability", 
                                      scenario = "SCF")

# Example of your cross-scenario calculation:
df_Histo$SigmaSCF <- df_Histo$Sigma0 * df_SCF$variability


res_path <- file.path(dataDir, "SCF/Var_Trend_agg.Rdata")

message("Loading: ", res_path)
load(res_path, envir = .GlobalEnv) 


#reobtain the original parameters
VariabilitySave=ResTV[[2]]
names(VariabilitySave)[1]="time"
Var_long <- VariabilitySave %>%
  pivot_longer(
    cols = -time,           # Keep 'time' as is, pivot everything else
    names_to = "catchment",   # Column headers become 'location'
    values_to = "variability"      # Cell values become 'value'
  )

Var_long$Year=year(Var_long$time)
Var_long$catchment=as.numeric(Var_long$catchment)

ParamVar<- Var_long %>%
  inner_join(ParamsflSCF, by = c("catchment", "Year"))

ParamVar$Sigma0=ParamVar$sigmaGPD/ParamVar$variability

TrendSave=ResTV[[1]]
names(TrendSave)[1]="time"
Trend_long <- TrendSave %>%
  pivot_longer(
    cols = -time,           # Keep 'time' as is, pivot everything else
    names_to = "catchment",   # Column headers become 'location'
    values_to = "trend"      # Cell values become 'value'
  )

Trend_long$Year=year(Trend_long$time)
Trend_long$catchment=as.numeric(Trend_long$catchment)
trendfuck=Trend_long[which(Trend_long$catchment==tc),]
ParamVarTrend<- ParamVar %>%
  inner_join(Trend_long, by = c("catchment", "Year"))

ParamVarTrend$threshold0=(ParamVarTrend$thresholdGPD-ParamVarTrend$trend)/ParamVarTrend$variability
#plot changes in trend and variability

#load result file for format
load(file=paste0("D:/tilloal/Documents/LFRuns_utils/data/TSEVA/output_plots/Flood_pixChange_RL100_v1.Rdata"))

load(file=paste0(dataDir,"/Histo/Var_Trend_agg.Rdata"))
###load Histo run -----
hydroDir<-("D:/tilloal/Documents/LFRuns_utils/ChangingHydroExtremes/data")
haz="Flood"
if (haz == "Drought") namefile="Drought.nonfrost.Histo"
if (haz == "Flood") namefile="Flood.year.Histo"
load(file=paste0(hydroDir,"/",haz,"/params.",namefile,".Rdata"))

Paramsfl=Paramsfl[,-c(4:9,17)]
ParamsflH=data.table(Paramsfl)
rm(Paramsfl)
gc()

#reobtain the original parameters
VariabilitySaveH=ResTV[[2]]
names(VariabilitySaveH)[1]="time"
Var_long <- VariabilitySaveH %>%
  pivot_longer(
    cols = -time,           # Keep 'time' as is, pivot everything else
    names_to = "catchment",   # Column headers become 'location'
    values_to = "variability"      # Cell values become 'value'
  )

Var_long$Year=year(Var_long$time)
Var_long$catchment=as.numeric(Var_long$catchment)

ParamVarH<- Var_long %>%
  inner_join(ParamsflH, by = c("catchment", "Year"))

ParamVarH$Sigma0=ParamVarH$sigmaGPD/ParamVarH$variability


#create new version of sigmaSCF
ParamVarH$SigmaSCF=ParamVarH$Sigma0*ParamVar$variability


TrendSave=ResTV[[1]]
names(TrendSave)[1]="time"
Trend_long <- TrendSave %>%
  pivot_longer(
    cols = -time,           # Keep 'time' as is, pivot everything else
    names_to = "catchment",   # Column headers become 'location'
    values_to = "trend"      # Cell values become 'value'
  )

Trend_long$Year=year(Trend_long$time)
Trend_long$catchment=as.numeric(Trend_long$catchment)

ParamVarTrendH<- ParamVarH %>%
  inner_join(Trend_long, by = c("catchment", "Year"))

process_catchment_params <- function(res_list, params_ref) {
  require(dplyr)
  require(tidyr)
  require(lubridate)
  
  # 1. Process Variability Data
  variability_raw <- res_list[[2]]
  names(variability_raw)[1] <- "time"
  
  var_long <- variability_raw %>%
    pivot_longer(
      cols = -time,
      names_to = "catchment",
      values_to = "variability"
    ) %>%
    mutate(
      Year = year(time),
      catchment = as.numeric(catchment)
    )
  
  # 2. Join with Parameters and Calculate Sigmas
  param_var <- var_long %>%
    inner_join(params_ref, by = c("catchment", "Year")) %>%
    mutate(
      Sigma0 = sigmaGPD / variability,
      threshold0=(thresholdGPD-trend)/variability,
      # Note: Using internal 'variability' col
    )
  
  # 3. Process Trend Data
  trend_raw <- res_list[[1]]
  names(trend_raw)[1] <- "time"
  
  trend_long <- trend_raw %>%
    pivot_longer(
      cols = -time,
      names_to = "catchment",
      values_to = "trend"
    ) %>%
    mutate(
      Year = year(time),
      catchment = as.numeric(catchment)
    )
  
  # 4. Final Join
  final_df <- param_var %>%
    inner_join(trend_long, by = c("catchment", "Year", "time"))
  
  return(final_df)
}

# Simply call the function and assign it to a name
ParamVarTrendHx <- process_catchment_params(ResTV, ParamsflH)




#keep improving this, include other runs

#Save parameters and recompute RLs with other script


ParamVarTrendH$threshold0=(ParamVarTrendH$thresholdGPD-ParamVarTrendH$trend)/ParamVarTrendH$variability
ParamVarTrendH$thresholdSCF=ParamVarTrendH$threshold0*ParamVarTrend$variability+ParamVarTrend$trend
ParamVarTrendH$variabilitySCF=ParamVarTrend$variability
ParamVarTrendH$trendSCF=ParamVarTrend$trend


ParamVarTrend$thresholdH=ParamVarTrend$threshold0*ParamVarTrendH$variability+ParamVarTrendH$trend
ParamVarTrend$variabilityH=ParamVarTrendH$variability
ParamVarTrend$SigmaH=ParamVarTrend$Sigma0*ParamVarTrendH$variability
ParamVarTrend$trendH=ParamVarTrendH$trend

#plot changes in trend and variability

ParamVarTrendH$dth=ParamVarTrendH$thresholdGPD-ParamVarTrendH$thresholdSCF

ParamVarTrend$dth=(ParamVarTrend$thresholdH-ParamVarTrend$thresholdGPD)/ParamVarTrend$thresholdGPD*100
ParamVarTrend$dSig=ParamVarTrend$SigmaH/ParamVarTrend$sigmaGPD


rm(ParamVarH,ParamVar,ParamsflH,ParamsflSCF)
gc()

#some tests
ParamSCF_1=ParamVarTrend[which(ParamVarTrend$Year==2020),]

ParamH_1=ParamVarTrendH[which(ParamVarTrendH$Year==2020),]

cloc=ParamVarTrendH$catchment[which.min(ParamVarTrendH$dth)]

Paramcat1=ParamVarTrendH[which(ParamVarTrendH$catchment==cloc),]
Paramcat2=ParamVarTrend[which(ParamVarTrend$catchment==cloc),]

plot(Paramcat2$SigmaH/Paramcat2$sigmaGPD)
#create stationnary params

#Now work with ParamSCF1

#Loading fitting results from the 4 runs and for all 282 000 river pixels in
# in the domain. Requires at least 20 GB of free RAM.
###load historical run -----
yname=seq(1955,2015, by=10)

haz="Flood"
library(tidyr)
vs="dSig"
vs="dth"
idv=match(vs,colnames(ParamVarTrend))
Var_long=ParamVarTrend[,c(1,2,idv,4)]
df_wide <- Var_long[,-1] %>%
  pivot_wider(
    names_from = Year,       # The values in 'year' become column headers
    values_from = vs      # The values in 'value' fill the cells
  )




cd=as.numeric(which(colnames(df_wide)=="1955"))
df_wid=data.frame(df_wide)
##4.1 change from climate--------------------

if (vs=="thresholdGPD"){
  DataV=data.frame(outl2=df_wid[,1],(df_wid[,c(2:70)]-(df_wid[,cd]))/rowMeans(df_wid[,c(2:70)],(df_wid[,cd]))*100)
}
if (vs=="sigmaGPD"){
  DataV=data.frame(outl2=df_wid[,1],df_wid[,c(2:70)]/(df_wid[,cd]))
}
if(vs=="dSig"){
  DataV=data.frame(outl2=df_wid[,1],df_wid[,c(2:70)])
}

if(vs=="dth"){
  DataV=data.frame(outl2=df_wid[,1],df_wid[,c(2:70)])
}


DataO=DataSave[[1]][,-c(12:81)]

DataV=inner_join(DataO,DataV,by="outl2")
yrange=match(paste0("X",yname),names(DataV))
trendagg <- aggregate(list(Rchange = DataV[, yrange]),
                      by = list(HydroR = DataV$HER),
                      FUN = function(x) c(mean = mean(x, na.rm = TRUE)))
trendClim <- do.call(data.frame, trendagg)


pointagg <- aggregate(list(Rchange_rel = DataV$X2015),
                      by = list(HydroR = DataV$HER),
                      FUN = function(x) c(mean = mean(x, na.rm = TRUE), 
                                          dev = sd(x, na.rm = TRUE), 
                                          len = length(x), 
                                          med = median(x, na.rm = TRUE), 
                                          q1 = quantile(x, 0.025, na.rm = TRUE), 
                                          q3 = quantile(x, 0.975, na.rm = TRUE)))
pointClim <- do.call(data.frame, pointagg)



plot(trendClim$Rchange.X2015)


#map of trendclim




TdataLuse=processTrendData(trendSoc,DataL,id_var = "HydroR")
TdataRes=processTrendData(trendRes,DataR,id_var = "HydroR")
TdataWuse=processTrendData(trendWu,DataW,id_var = "HydroR")

lm=match(pointsag,RegioRLi$HydroR)
RegioRLi$HydroR[lm]

paggC=Dchangelist[[1]]$PagD
paggL=Dchangelist[[2]]$PagD
paggR=Dchangelist[[3]]$PagD
paggW=Dchangelist[[4]]$PagD
paggC$driver="Clim"
paggL$driver="Lu"
paggR$driver="Res"
paggW$driver="Wu"
colnames(paggC)
colnames(paggR)
pointsAD=rbind(paggC,paggL,paggR,paggW)

length(which(DataC$Y2020>15))/length(DataC$Y2020)
length(which(pointSoc$Rchange_rel.mean<0))/length(pointSoc$Rchange_rel.mean)


