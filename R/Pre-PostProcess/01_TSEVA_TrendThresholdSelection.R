# Step 1 of Large-scale TSEVA fitting: --------------------
## Identification of threshold for extreme trend estimates

####################### Libraries loading #############################
library(RtsEva)
library(lubridate)
library(zoo)
suppressWarnings(suppressMessages(library(ncdf4)))
suppressWarnings(suppressMessages(library(sf)))
suppressWarnings(suppressMessages(library(rnaturalearth)))
suppressWarnings(suppressMessages(library(rnaturalearthdata)))
suppressWarnings(suppressMessages(library(rgeos)))
suppressWarnings(suppressMessages(library(dplyr)))

# Function to open and extract data from a netCDF file containing basin information
outletopen <- function(dir, outletname, nrspace = rep(NA, 5)) {
  ncbassin <- paste0(dir, "/", outletname, ".nc")
  ncb <- nc_open(ncbassin)
  name.vb <- names(ncb[['var']])
  namev <- name.vb[1]
  if ("Band1" %in% name.vb) namev <- "Band1"
  name.lon <- "lon"
  name.lat <- "lat"
  
  if (!is.na(nrspace[1])) {
    start <- as.numeric(nrspace[c(2, 4)])
    count <- as.numeric(nrspace[c(3, 5)]) - start + 1
  } else {
    londat <- ncvar_get(ncb, name.lon) 
    llo <- length(londat)
    latdat <- ncvar_get(ncb, name.lat)
    lla <- length(latdat)
    start <- c(1, 1)
    count <- c(llo, lla)
  }
  
  londat <- ncvar_get(ncb, name.lon, start = start[1], count = count[1]) 
  llo <- length(londat)
  latdat <- ncvar_get(ncb, name.lat, start = start[2], count = count[2])
  lla <- length(latdat)
  outlets <- ncvar_get(ncb, namev, start = start, count = count) 
  outlets <- as.vector(outlets)
  outll <- expand.grid(londat, latdat)
  lonlatloop <- expand.grid(c(1:llo), c(1:lla))
  outll$idlo <- lonlatloop$Var1
  outll$idla <- lonlatloop$Var2
  
  outll <- outll[which(!is.na(outlets)),]
  outlets <- outlets[which(!is.na(outlets))]
  outll <- data.frame(outlets, outll)
  return(outll)
}

# Function to open and extract discharge data from a netCDF file for a specific location
disNcopenloc <- function(fname, dir, outloc, idc) {
  ncdis <- paste0(dir, "/", fname, ".nc")
  ncd <- nc_open(ncdis)
  name.vb <- names(ncd[['var']])
  namev <- name.vb[1]
  time <- ncvar_get(ncd, "time")
  lt <- length(time)
  
  # Timestamp correction
  name.lon <- "lon"
  name.lat <- "lat"
  londat <- ncvar_get(ncd, name.lon) 
  llo <- length(londat)
  latdat <- ncvar_get(ncd, name.lat) 
  lla <- length(latdat)
  
  idm <- 1 + (idc - 1) * lt
  start <- c(outloc$idlo[idc], outloc$idla[idc], 1)
  count <- c(1, 1, lt)
  outlets <- ncvar_get(ncd, namev, start = start, count = count)
  outlets <- as.vector(outlets)
  outid <- rep(outloc[idc, 1], length(time))
  lon <- rep(londat[start[1]], length(time))
  lat <- rep(latdat[start[2]], length(time))
  outll <- data.frame(outlets, outid, lon, lat, time)
  
  return(outll)
}

# Function to check if a time series contains data for all years in a given range
check_timeserie2 <- function(timeseries, yro) {
  ts_years <- as.integer(lubridate::year(timeseries))
  year_check <- yro %in% ts_years
  runs <- rle(year_check)
  rf <- which(runs$values == FALSE)
  if (any(runs$lengths[rf] >= 2)) {
    return(FALSE)
  } else {
    return(TRUE)
  }
}

# Function to find the trend threshold for a time series
tsEvaFindTrendThreshold2 <- function(series, timeStamps, timeWindow) {
  ptn <- timeStamps[which(!is.na(series))]
  bounds <- unique(lubridate::year(ptn))
  
  sts <- c()
  lnegs <- c()
  pctd <- c()
  pcts <- seq(0.4, 0.95, by = 0.05)
  
  rs <- tsEvaDetrendTimeSeries(timeStamps, series, timeWindow, fast = TRUE)
  
  nr <- rep(1, length(rs@trendSeries))
  nr <- nr + rnorm(length(rs@trendSeries), 0, 1e-05)
  
  for (iter in 1:length(pcts)) {
    thrsdt <- quantile(series, pcts[iter], na.rm = TRUE)
    series_no_na <- series
    series_no_na[which(is.na(series_no_na))] <- -9999
    serieb <- series_no_na
    timeb <- timeStamps
    timeb <- timeb[-which(serieb < thrsdt)]
    serieb[which(serieb < thrsdt)] <- NA
    checkY <- check_timeserie2(timeb, bounds)
    if (checkY == FALSE) {
      print(paste0("not all years - q= ", pcts[iter]))
      break
    }
    rs <- tsEvaDetrendTimeSeries(timeStamps, serieb, timeWindow, fast = TRUE)
    
    norm_trend <- rs@trendSeries / mean(rs@trendSeries, na.rm = TRUE)
    dtr1 <- serieb - rs@trendSeries
    lneg <- length(which(dtr1 < 0))
    stab <- cor(nr, norm_trend, use = "pairwise.complete.obs")
    if (iter == 1) nr <- norm_trend
    lnegs <- c(lnegs, lneg)
    sts <- c(sts, stab)
    pctd <- c(pctd, pcts[iter])
  }
  
  rval <- pctd[length(pctd)]
  if (length(sts) > 3) {
    dow <- abs(diff(sts))[-1]
    print(dow)
    if (max(dow, na.rm = TRUE) > 0.2) {
      print("breaking point")
      rval <- pctd[which.max(dow) + 1]
    }
  }
  if (sum(lnegs) > 1) {
    rval <- pctd[which.min(lnegs)]
  }
  print(pctd)
  return(rval)
}

####################### Arguments importation #############################

# Discharge data from HERA was previously divided into spatial chunks
# Square number
Nsq <- 33

# Tail or extreme studied: can be "high" or "low"
tail <- "high"

# Scenario or run analysed. Here we use a constant threshold for extreme trend estimates across 4 OS lisflood runs scenarios:
# 1. Socioeconomic counterfactual
# 2. Reservoir+Water demand counterfactual
# 3. Water demand counterfactual
# 4. Historical run (HERA)
# For threshold identification, we only consider the two most "extreme" scenarios: Scenario 1 and Scenario 4.
sce <- "Histo"
outlets <- "RNetwork"
var <- "dis"

# Loading repositories where the data is stored
workDir <- "/BGFS/CLIMEX/tilloal/HydroMeteo/TSEVA/"
hydroDir <- "/BGFS/CLIMEX/tilloal/HydroMeteo"
setwd(workDir)

# Metadata about the domain spatial division into chunks
rspace <- read.csv(paste0(hydroDir, "/subspace_efas.csv"))
rspace <- rspace[, -1]
nrspace <- rspace[Nsq,]

# Load river pixel locations and identifiers
outletname <- "efas_rnet_100km_01min"
nameout <- "UCRnet"
outhybas <- outletopen(hydroDir, outletname, nrspace)
Idstart <- as.numeric(Nsq) * 10000
if (length(outhybas$outlets) > 0) {
  outhybas$outlets <- seq((Idstart + 1), (Idstart + length(outhybas$outlets)))
}

unikout <- outhybas$outlets
outhybas$latlong <- paste(round(outhybas$Var1, 4), round(outhybas$Var2, 4), sep = " ")

# Check if file already exists before launching computation
file_out <- paste0(workDir, "out/Thresholds/trenTH_31d_", sce, "_", tail, "_", Nsq, ".csv")
fex <- file.exists(file_out)
if (fex == TRUE) {
  f.size <- as.numeric(fs::file_size(file_out))
} else {
  f.size <- 0
}

if (length(outhybas$outlets) > 0 & f.size < 500) {
  
  timeWindow <- 365.25 * 30  # Time window in days, the correction is done within the functions
  
  # Loading discharge data for two scenarios
  if (sce == "Histo") {
    filename <- paste0("/Timeseries/dis6_Histo/dis_", Nsq, "_1951_2020_h")
  }
  if (sce == "SCF") {
    filename <- paste0("/Timeseries/dis6_SCF/dis_", Nsq, "_1951_2020_scf")
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
  rmv <- which(year(txx) == 1950)
  if (length(rmv) > 1) {
    df.dis <- df.dis[-rmv,]
    txx <- txx[-rmv]
  }
  
  Station_data_IDs <- unikout
  Trth_H_list <- c()
  
  for (id in c(1:length(unikout))) {
    print(id)
    stid <- Station_data_IDs[id]
    catch <- stid
    dists <- disNcopenloc(filename, hydroDir, outhybas, id)
    df.dis <- dists 
    if (length(rmv) > 1) {
      df.dis <- df.dis[-rmv,]
    }
    data <- data.frame(txx, df.dis$outlets)
    names(data) <- c("date", "Qs")
    
    # Preparation of discharge timeseries for the desired extreme analysis
    if (tail == "high") {
      # Extract the maximum daily value
      series <- max_daily_value(data)
      timeAndSeriesH <- series
      names(timeAndSeriesH) <- c("timestamp", "data")
    } else if (tail == "low") {
      minPeakDistanceInDays <- 31
      # 7 day average flow for drought analysis
      WindowSize2 <- 7
      names(data) <- c("date", "Qs")
      dt1 <- min(diff(data$date), na.rm = TRUE)
      dt <- as.numeric(dt1)
      tdim <- attributes(dt1)$units
      if (tdim == "hours") dt <- dt / 24
      nRunMn <- ceiling(WindowSize2 / dt);
      data$Q7 <- rollmean(data$Qs, nRunMn, align = "right", fill = NA)
      timeStamps <- data$date
      series <- data$Q7
      start_index <- 1
      indices_to_extract <- seq(from = start_index, to = length(series), by = WindowSize2 / dt)
      series <- series[indices_to_extract]
      series <- -1 * series
      timeStamps <- timeStamps[indices_to_extract]
      timeAndSeriesH <- data.frame(timeStamps, series)
      names(timeAndSeriesH) <- c("timestamp", "data")
    }
    
    TrendTh_H <- tsEvaFindTrendThreshold(series = timeAndSeriesH$data, timeStamps = timeAndSeriesH$timestamp, timeWindow)
    
    if (is.null(TrendTh_H)) TrendTh_H <- NA
    if (length(TrendTh_H) == 0) TrendTh_H <- NA
    
    print(paste0("catch:", catch, " Threshold: ", TrendTh_H1))
    Trth_H <- c(catch, TrendTh_H)
    
    names(Trth_H) <- c("cid", "Th_new")
    
    Trth_H_list <- rbind(Trth_H_list, Trth_H)
  }
  
  # Saving outputs
  write.csv(Trth_H_list, file = paste0(workDir, "out/Thresholds/trenTH_31d_", sce, "_", tail, "_", Nsq, ".csv"))
}
