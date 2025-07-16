#Code to identify threshold for extreme trend assessment with TSEVA
#The method is applied here to a subsample of the 282 521 river pixels analysed 
#in the research article


setwd(dirname(rstudioapi::getSourceEditorContext()$path))
getwd()
source("functions_trends.R")

#Set data directory
hydroDir<-("D:/tilloal/Documents/LFRuns_utils/ChangingHydroExtremes/data")

# Discharge data from HERA was previously divided into spatial chunks
# Square number
Nsq <- 42

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
  filename <- paste0("/Timeseries/dis_", Nsq, "_1951_2020_h")
}
if (sce == "SCF") {
  filename <- paste0("/Timeseries/dis_", Nsq, "_1951_2020_scf")
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


Station_data_IDs <- unikout

# 
# dis_path<-paste0(main_path,'HERA/')
# load(file=paste0(dis_path,"HERA6h_H07_19502020.Rdata"))
# Station_data_IDs <- as.vector(t(Q_sim[1, ]))[-1]
# tqx=as.POSIXct((Q_sim[,1]*3600*24), origin="1979-01-01 00:00:00")
# txx=tqx[-1]
# Q_simH=Q_sim
# 
# #extract values for 1 location
# #remove locations without discharge
# Q_simH=Q_simH[,-which(is.na(Q_simH[2,]))]
# Station_data_IDs=Station_data_IDs[-wolala]
# txx=tqx[-1]
# rmv=which(year(txx)==1950)
# txx=txx[-rmv]

timeWindow = 365.25*30; #time windows in days, the correction is done within the functions

# Register parallel backend
cl <- makeCluster(4)
registerDoParallel(cl)
Trth_H_list=c()
tail="low"
Trth_H_list <- foreach(id = 1:length(Station_data_IDs), .combine = rbind, .packages=c("RtsEva","zoo")) %dopar% {
  #print(id)
  # id=1
  # stid <- Station_data_IDs[id]
  # catch <- stid
  # Q_H <- Q_simH[,id+1][-1]
  # Q_H <- Q_H[-rmv]
  
  stid <- Station_data_IDs[id]
  catch <- stid
  dists <- disNcopenloc(filename, hydroDir, outhybas, id)
  df.dis <- dists 
  if (length(rmv) > 1) {
    df.dis <- df.dis[-rmv,]
  }
  data <- data.frame(txx, df.dis$outlets)
  names(data) <- c("date", "Qs")
  # data <- data.frame(txx,Q_H)
  # names(data) <- c("date","Qs")
  if (tail=="high"){
    # Extract the maximum daily value
    series <- max_daily_value(data)
    timeAndSeriesH=series
    names(timeAndSeriesH)=c("timestamp","data")
  }else if (tail=="low"){
    minPeakDistanceInDays=30
    #7 day average flow for drought analysis
    WindowSize=7
    names(data)=c("date","Qs")
    dt1=min(diff(data$date),na.rm=T)
    dt=as.numeric(dt1)
    tdim=attributes(dt1)$units
    if (tdim=="hours") dt=dt/24
    nRunMn = ceiling(WindowSize/dt);
    colnames(data)
    data$Q7=rollmean(data$Qs,nRunMn, align = "right", fill=NA)
    timeStamps=data$date
    series=data$Q7
    start_index=1
    indices_to_extract <- seq(from = start_index, to = length(series), by = WindowSize/dt)
    series=series[indices_to_extract]
    series=-1*series
    timeStamps=timeStamps[indices_to_extract]
    timeAndSeriesH=data.frame(timeStamps,series)
    names(timeAndSeriesH)=c("timestamp","data")
  }
  
  TrendTh_H <- tsEvaFindTrendThreshold2(series=timeAndSeriesH$data, timeStamps=timeAndSeriesH$timestamp, timeWindow)
  m1=median(timeAndSeriesH$data)
  o=findpeaks(timeAndSeriesH$data, threshold = m1)
  plot(o[,1])
  mean(o[,1])
  m1=mean(timeAndSeriesH$data)
  if (is.null(TrendTh_H))TrendTh_H=NA
  Trth_H <- c(catch,TrendTh_H)
}

stopCluster(cl)
write.csv(Trth_H_list,file=paste0("trenTH_H",tail,".csv"))
rm(Q_simH)

dis_path<-paste0(main_path,'HERA_CFS/')
load(file=paste0(dis_path,"HERA_CFS6h_H07_19502020.Rdata"))
Q_simCFS=Q_sim
Q_simCFS=Q_simCFS[,-which(is.na(Q_simCFS[2,]))]

cl <- makeCluster(4)
registerDoParallel(cl)
Trth_S_list <- foreach(id = 1:length(Station_data_IDs), .combine = rbind, .packages=c("RtsEva","zoo")) %dopar% {
  
  
  #print(id)
  stid <- Station_data_IDs[id]
  catch <- stid
  
  Q_SCF <- Q_simCFS[,id+1][-1]
  Q_SCF <- Q_SCF[-rmv]

  data <- data.frame(txx,Q_SCF)
  names(data) <- c("date","Qs")
  if (tail=="high"){
    # Extract the maximum daily value
    series <- max_daily_value(data)
    timeAndSeriesS=series
    names(timeAndSeriesS)=c("timestamp","data")
  } else if (tail=="low"){
    
    minPeakDistanceInDays=30
    #7 day average flow for drought analysis
    WindowSize=7
    names(data)=c("date","Qs")
    dt1=min(diff(data$date),na.rm=T)
    dt=as.numeric(dt1)
    tdim=attributes(dt1)$units
    if (tdim=="hours") dt=dt/24
    nRunMn = ceiling(WindowSize/dt);
    colnames(data)
    data$Q7=rollmean(data$Qs,nRunMn, align = "right", fill=NA)
    timeStamps=data$date
    series=data$Q7
    start_index=1
    indices_to_extract <- seq(from = start_index, to = length(series), by = WindowSize/dt)
    series=series[indices_to_extract]
    series=-1*series
    timeStamps=timeStamps[indices_to_extract]
    timeAndSeriesS=data.frame(timeStamps,series)
    names(timeAndSeriesS)=c("timestamp","data")
  }
  
  TrendTh_S <- tsEvaFindTrendThreshold2(series=timeAndSeriesS$data, timeStamps=timeAndSeriesS$timestamp, timeWindow)
  if (is.null(TrendTh_S2))TrendTh_S=NA
  Trth_S <- c(catch,TrendTh_S)
}
stopCluster(cl)
write.csv(Trth_S_list,file=paste0("trenTH_S",tail,".csv"))

plot(Trth_H_list[,2], Trth_S_list[,2],ylim=c(0.5,1),xlim=c(0.5,1))

dx=data.frame(hist=Trth_H_list[,2],scf=Trth_S_list[,2], catch=Trth_H_list[,1])
dx2=dx[which(dx$hist!=dx$scf),]
dx$common=dx$hist
dx$common[which(dx$hist>1)]=0.5
save(dx,file=paste0("thresholds_",tail,"_catchments.rdata"))
