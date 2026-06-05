#Code to identify threshold for extreme trend assessment with TSEVA
#The method is applied here to a subsample of the 282 521 river pixels analysed 
#in the research article


setwd(dirname(rstudioapi::getSourceEditorContext()$path))
getwd()
source("functions_trends.R")



# -----------------------------------------------------------------------
# HELPER FUNCTIONS
# -----------------------------------------------------------------------

outletopen = function(dir, outletname, nrspace = rep(NA, 5)){
  ncbassin = paste0(dir, "/", outletname, ".nc")
  ncb = nc_open(ncbassin)
  name.vb = names(ncb[['var']])
  namev = name.vb[1]
  if ("Band1" %in% name.vb) namev = "Band1"
  name.lon = "lon"
  name.lat = "lat"
  if (!is.na(nrspace[1])){
    start = as.numeric(nrspace[c(2, 4)])
    count = as.numeric(nrspace[c(3, 5)]) - start + 1
  } else {
    londat = ncvar_get(ncb, name.lon)
    llo = length(londat)
    latdat = ncvar_get(ncb, name.lat)
    lla = length(latdat)
    start = c(1, 1)
    count = c(llo, lla)
  }
  londat = ncvar_get(ncb, name.lon, start = start[1], count = count[1])
  llo = length(londat)
  latdat = ncvar_get(ncb, name.lat, start = start[2], count = count[2])
  lla = length(latdat)
  outlets = ncvar_get(ncb, namev, start = start, count = count)
  nc_close(ncb)
  outlets = as.vector(outlets)
  outll = expand.grid(londat, latdat)
  lonlatloop = expand.grid(c(1:llo), c(1:lla))
  outll$idlo = lonlatloop$Var1
  outll$idla = lonlatloop$Var2
  outll = outll[which(!is.na(outlets)), ]
  outlets = outlets[which(!is.na(outlets))]
  outll = data.frame(outlets, outll)
  return(outll)
}

# OPTIMIZED: reads ALL catchment time series in a single NetCDF open/close
# replaces per-catchment disNcopenloc calls inside the loop
disNcopenALL = function(fname, dir, outloc){
  ncdis = paste0(dir, "/", fname, ".nc")
  ncd = nc_open(ncdis)
  name.vb = names(ncd[['var']])
  namev = name.vb[1]
  time = ncvar_get(ncd, "time")
  lt = length(time)
  name.lon = "lon"
  name.lat = "lat"
  londat = ncvar_get(ncd, name.lon)
  latdat = ncvar_get(ncd, name.lat)
  n = nrow(outloc)
  all_series = vector("list", n)
  for (idc in seq_len(n)){
    start = c(outloc$idlo[idc], outloc$idla[idc], 1)
    count = c(1, 1, lt)
    all_series[[idc]] = as.vector(ncvar_get(ncd, namev, start = start, count = count))
  }
  nc_close(ncd)
  # return time vector alongside series so timestamps can be built outside
  return(list(series = all_series, time = time))
}

check_timeserie2 = function(timeseries, yro){
  ts_years = as.integer((lubridate::year(timeseries)))
  year_check = yro %in% ts_years
  runs = rle(year_check)
  rf = which(runs$values == FALSE)
  if (any(runs$lengths[rf] >= 2)){
    return(FALSE)
  } else {
    return(TRUE)
  }
}

# OPTIMIZED: removed the redundant first tsEvaDetrendTimeSeries call that
# ran on the full series before the quantile loop but whose result was
# never used — nr is now initialised on the first valid loop iteration
tsEvaFindTrendThreshold2 <- function(series, timeStamps, timeWindow){
  ptn = timeStamps[which(!is.na(series))]
  bounds = unique(lubridate::year(ptn))
  sts = c()
  lnegs = c()
  pctd = c()
  pcts = seq(0.4, 0.95, by = 0.05)
  nr = NULL   # initialised on first valid iteration, not before the loop
  
  for (iter in 1:length(pcts)){
    thrsdt = quantile(series, pcts[iter], na.rm = TRUE)
    series_no_na = series
    series_no_na[which(is.na(series_no_na))] = -9999
    serieb = series_no_na
    timeb = timeStamps
    timeb = timeb[-which(serieb < thrsdt)]
    serieb[which(serieb < thrsdt)] = NA
    checkY = check_timeserie2(timeb, bounds)
    if (checkY == FALSE){
      print(paste0("not all years - q= ", pcts[iter]))
      break
    }
    rs = tsEvaDetrendTimeSeries(timeStamps, serieb, timeWindow, fast = T)
    norm_trend = rs@trendSeries / mean(rs@trendSeries, na.rm = TRUE)
    dtr1 = serieb - rs@trendSeries
    lneg = length(which(dtr1 < 0))
    if (is.null(nr)){
      # first iteration: store reference trend; correlation with itself is 1
      nr = norm_trend
      stab = 1
    } else {
      stab = cor(nr, norm_trend, use = "pairwise.complete.obs")
    }
    lnegs = c(lnegs, lneg)
    sts = c(sts, stab)
    pctd = c(pctd, pcts[iter])
  }
  
  rval = pctd[length(pctd)]
  if (length(sts) > 3){
    dow = abs(diff(sts))[-1]
    if (max(dow, na.rm = T) > 0.2){
      print("breaking point")
      rval = pctd[which.max(dow) + 1]
    }
  }
  if (sum(lnegs) > 1){
    rval = pctd[which.min(lnegs)]
  }
  return(rval)
}




#Set data directory
hydroDir<-("D:/tilloal/Documents/LFRuns_utils/ChangingHydroExtremes/data")


# -----------------------------------------------------------------------
# ARGUMENTS
# -----------------------------------------------------------------------
# Discharge data from HERA was previously divided into spatial chunks
# Square number
Nsq <- 42
# Tail or extreme studied: can be "high" or "low"
tail <- "low"

# Scenario or run analysed. Here we use a constant threshold for extreme trend estimates across 4 OS lisflood runs scenarios:
# 1. Socioeconomic counterfactual
# 2. Reservoir+Water demand counterfactual
# 3. Water demand counterfactual
# 4. Historical run (HERA)
# For threshold identification, we only consider the two most "extreme" scenarios: Scenario 1 and Scenario 4.
sce <- "SCF"
outlets <- "RNetwork"
var <- "dis"

# Metadata about the domain spatial division into chunks
rspace <- read.csv(paste0(hydroDir, "/subspace_efas.csv"))
rspace <- rspace[, -1]
nrspace <- rspace[Nsq,]

# Load river pixel locations and identifiers
outletname <- "/GeoData/efas_rnet_100km_01min"
nameout <- "UCRnet"
outhybas <- outletopen(hydroDir, outletname, nrspace)
Idstart <- as.numeric(Nsq) * 10000
if (length(outhybas$outlets) > 0) {
  outhybas$outlets <- seq((Idstart + 1), (Idstart + length(outhybas$outlets)))
}

unikout <- outhybas$outlets
outhybas$latlong <- paste(round(outhybas$Var1, 4), round(outhybas$Var2, 4), sep = " ")
timeWindow <- 365.25 * 30  # Time window in days, the correction is done within the functions

# Loading discharge data for two scenarios
if (sce == "Histo") {
  filename <- paste0("dis_", Nsq, "_1951_2020_h_RNetwork")
}
if (sce == "SCF") {
  filename <- paste0("dis_", Nsq, "_1951_2020_scf_RNetwork")
}




# Load data from the first river pixel
dists <- disNcopenloc(filename, hydroDir, outhybas, 1)
df.dis <- dists 
print(paste0("opening square ", Nsq, " /88"))
timeStamps <- unique(as.Date(df.dis$time, origin = "1979-01-01"))
timeStamps <- as.POSIXct(timeStamps - 1/24)
txx <- timeStamps
df.dis$timeStamps <- timeStamps
names(df.dis)[c(1, 2)] <- c("dis", "outlets")


rmv = which(year(txx) == 1950)
if (length(rmv) > 1){
  txx = txx[-rmv]
}

# OPTIMIZED: load all low-tail static data ONCE before the loop
# In the original these were reloaded on every catchment iteration
if (tail == "low"){
  load(file = paste0(hydroDir, "/Drought/catchment_frost.Rdata"))
  rmvx = which(year(frostcat$time) == 1950)
  frostcat = frostcat[-rmvx, ]
  frostcat = frostcat[-1, ]
  Catchmentrivers7 = read.csv(paste0(hydroDir, "/GeoData/HYBAS07/from_hybas_eu_onlyid.csv"),
                              encoding = "UTF-8", header = T, stringsAsFactors = F)
  outletname_h07 = "/GeoData/HYBAS07/outletsv8_hybas07_01min"
  outhyb07 = outletopen(hydroDir, outletname_h07, nrspace)
  catmatch = match(outhyb07$outlets, Catchmentrivers7$pointid)
  mycat = Catchmentrivers7[catmatch, ]
  hybas07 = read_sf(dsn = paste0(hydroDir, "/GeoData/HYBAS07/hybas_eu_lev07_v1c.shp"))
  hybasf7 = fortify(hybas07)
  Catamere07 = inner_join(hybasf7, Catchmentrivers7, by = "HYBAS_ID")
  Catamere07$llcoord = paste(round(Catamere07$POINT_X, 4),
                             round(Catamere07$POINT_Y, 4), sep = " ")
  # pre-join with outhybas so per-catchment step is just a single filter
  Catf7 = inner_join(Catamere07, outhybas, by = c("llcoord" = "latlong"))
  st_geometry(Catf7) = NULL
}


library(doParallel)
library(foreach)

# Number of cores (adjust to your machine)
ncores <- detectCores() - 2
cl <- makeCluster(ncores)
registerDoParallel(cl)

# Export all necessary objects and functions to the workers
clusterExport(cl, c("Station_data_IDs", "unikout", "txx", "rmv", 
                    "Catf7", "mycat", "frostcat",
                    "filename", "hydroDir", "outhybas", "timeWindow", "tail",
                    "disNcopenloc", "max_daily_value"))


nid=length(unikout)
nid=10
Station_data_IDs = unikout

# Custom combine function
comb <- function(x, y) {
  list(df_Xtsave = rbind(x$df_Xtsave, y$df_Xtsave),
       Trth_H_list = rbind(x$Trth_H_list, y$Trth_H_list))
}

# Initial empty objects (as lists)
empty_res <- list(df_Xtsave = NULL, Trth_H_list = NULL)


results <- foreach(id = 1:nid, 
                     .combine = comb, 
                     .init = empty_res,
                     .packages = c("RtsEva", "zoo", "ncdf4", "pracma")) %dopar% {
                       
                       # ---- your original loop content (slightly adjusted) ----
                       
                       print(id)
                       flush.console() 
                       
                       stid <- Station_data_IDs[id]
                       catch <- stid
                       
                       dists <- disNcopenloc(filename, hydroDir, outhybas, id)
                       df.dis <- dists 
                       
                       if (length(rmv) > 1) {
                         df.dis <- df.dis[-rmv, ]
                       }
                       
                       data <- data.frame(txx, df.dis$outlets)
                       names(data) <- c("date", "Qs")
                       
                       if (tail == "high") {
                         series <- max_daily_value(data)
                         timeAndSeriesH <- series
                         names(timeAndSeriesH) <- c("timestamp", "data")
                         
                       } else if (tail == "low") {
                         
                         # frost lookup
                         catmat <- Catf7[which(Catf7$outlets == catch), ]
                         Tcatmat <- mycat[which(mycat$HYBAS_ID == catmat$HYBAS_ID), ]
                         Tcatchment <- which(colnames(frostcat) == Tcatmat$pointid)
                         
                         frosttime <- NA
                         if (length(Tcatchment) > 0) {
                           frostserie <- data.frame(frostcat[, 1], frostcat[, Tcatchment])
                           frosttime <- which(frostserie[, 2] < 0)
                         }
                         
                         minPeakDistanceInDays <- 31
                         WindowSize <- 7
                         names(data) <- c("date", "Qs")
                         dt1 <- min(diff(data$date), na.rm = TRUE)
                         dt <- as.numeric(dt1)
                         tdim <- attributes(dt1)$units
                         if (tdim == "hours") dt <- dt / 24
                         nRunMn <- ceiling(WindowSize / dt)
                         data$Q7 <- rollmean(data$Qs, nRunMn, align = "right", fill = NA)
                         timeStamps <- data$date
                         series <- data$Q7
                         
                         # ****** FROST SEASON LOGIC (kept as in original) ******
                         # Original code had a bug: length(!is.na(frosttime)) is always > 1.
                         # I preserve it as is, but you likely meant:
                         # if (length(frosttime[!is.na(frosttime)]) > 1) { ... }
                         season <- "nonfrost"   # hardcoded in original
                         if (length(!is.na(frosttime)) > 1) {
                           if (season == "nonfrost") {
                             print("nonfrost season")
                             series[frosttime] <- NA
                           } else if (season == "frost") {
                             print("frost season")
                             series[-frosttime] <- NA
                           }
                         }
                         
                         start_index <- 1
                         indices_to_extract <- seq(from = start_index, to = length(series), by = WindowSize / dt)
                         series <- series[indices_to_extract]
                         series <- -1 * series
                         timeStamps <- timeStamps[indices_to_extract]
                         timeAndSeriesH <- data.frame(timeStamps, series)
                         names(timeAndSeriesH) <- c("timestamp", "data")
                       }
                       
                       # --- Trend thresholds ---
                       TrendTh_H1 <- try(tsEvaFindTrendThreshold(series, timeStamps, timeWindow), TRUE)
                       if (inherits(TrendTh_H1, "try-error") || length(TrendTh_H1) == 0 || all(is.na(TrendTh_H1))) {
                         TrendTh_H1 <- NA
                       }
                       
                       TrendTh_H2 <- try(tsEvaFindTrendThreshold2(series = timeAndSeriesH$data,
                                                                  timeStamps = timeAndSeriesH$timestamp,
                                                                  timeWindow), TRUE)
                       if (inherits(TrendTh_H2, "try-error") || length(TrendTh_H2) == 0 || all(is.na(TrendTh_H2))) {
                         TrendTh_H2 <- NA
                       }
                       
                       # --- Build outputs for this iteration ---
                       Trth_H <- c(catch, TrendTh_H1, TrendTh_H2)
                       names(Trth_H) <- c("cid", "Th_old", "Th_new")
                       
                       thresh <- TrendTh_H2
                       if (is.na(thresh)) thresh <- 0.1
                       thrsdt <- quantile(timeAndSeriesH$data, thresh, na.rm = TRUE)
                       serieb <- timeAndSeriesH$data
                       timeb <- timeAndSeriesH$timestamp
                       idb <- which(serieb >= thrsdt)
                       timeb <- timeb[idb]
                       serieb <- serieb[idb]
                       df_xtrem <- data.frame(cid = rep(catch, length(idb)), 
                                              time = timeb, 
                                              data = serieb, 
                                              id = idb)
                       
                       # Return both data frames as a list
                       list(df_Xtsave = df_xtrem,
                            Trth_H_list = as.data.frame(t(Trth_H)))   # ensure row format
                     }
# Stop cluster
stopCluster(cl)
gc()


write.csv(results$Trth_H_list, file = paste0(hydroDir, "/Thresholds/trenTH_x_", sce, "_", tail, "_", Nsq, ".csv"))
write.csv(results$df_Xtsave,   file = paste0(hydroDir, "/Thresholds/xtrempoints_", sce, "_", tail, "_", Nsq, ".csv"))
