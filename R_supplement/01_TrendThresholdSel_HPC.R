#Code to identify trend threshold for each catchment, 
#if threshold is different, I need a flag and I check manually what is the best

getwd()
setwd("D:/tilloal/Documents/LFRuns_utils/TrendAnalysis/")
source("R/functions_trends.R")
library(RtsEva)
library(lubridate)
library(zoo)
suppressWarnings(suppressMessages(library(ncdf4)))
suppressWarnings(suppressMessages(library(sf)))
suppressWarnings(suppressMessages(library(rnaturalearth)))
suppressWarnings(suppressMessages(library(rnaturalearthdata)))
suppressWarnings(suppressMessages(library(rgeos)))
suppressWarnings(suppressMessages(library(dplyr)))
main_path = "D:/tilloal/Documents/06_Floodrivers/"
valid_path = paste0(main_path)



disNcopenloc=function(fname,dir,outloc,idc){
  ncdis=paste0(dir,"/",fname,".nc")
  ncd=nc_open(ncdis)
  name.vb=names(ncd[['var']])
  namev=name.vb[1]
  time <- ncvar_get(ncd,"time")
  lt=length(time)
  
  #timestamp corretion
  name.lon="lon"
  name.lat="lat"
  londat = ncvar_get(ncd,name.lon) 
  llo=length(londat)
  latdat = ncvar_get(ncd,name.lat) 
  lla=length(latdat)
  
  idm=1+(idc-1)*lt
  start=c(outloc$idlo[idc],outloc$idla[idc],1)
  count=c(1,1,lt)
  outlets = ncvar_get(ncd,namev,start = start, count= count)
  outlets=as.vector(outlets)
  outid=rep(outloc[idc,1],length(time))
  #lonlatloop=expand.grid(c(1:lla),c(1:lt))
  lon=rep(londat[start[1]],length(time))
  lat=rep(latdat[start[2]],length(time))
  outll=data.frame(outlets,outid,lon,lat,time)
  
  return (outll)
}
check_timeserie2=function(timeseries,yro){
  ts_years <- as.integer((lubridate::year(timeseries)))
  year_check <- yro %in% ts_years
  runs <- rle(year_check)
  rf=which(runs$values==FALSE)
  if (any(runs$lengths[rf] >= 4)) {
    return(FALSE)
  }
  else {
    return(TRUE)
  }
}
tsEvaFindTrendThreshold2<-function(series, timeStamps, timeWindow){
  ptn = timeStamps[which(!is.na(series))]
  bounds = unique(lubridate::year(ptn))
  nr <- rep(1, length(series))
  nr = nr + rnorm(length(series), 0, 1e-05)
  sts <- c()
  lnegs = c()
  pctd = c()
  pcts <- seq(0.4, 0.95, by = 0.05)
  for (iter in 1:length(pcts)) {
    thrsdt <- quantile(series, pcts[iter], na.rm = TRUE)
    series_no_na <- series
    series_no_na[which(is.na(series_no_na))] <- -9999
    serieb <- series_no_na
    timeb = timeStamps
    timeb = timeb[-which(serieb < thrsdt)]
    serieb[which(serieb < thrsdt)] <- NA
    checkY=check_timeserie2(timeb,bounds)
    if (checkY == FALSE) {
      print(paste0("not all years - q= ",pcts[iter]))
      break
    }
    rs <- tsEvaDetrendTimeSeries(timeStamps, serieb, timeWindow, 
                                 fast = T)
    # varianceSeries <- tsEvaNanRunningVariance(serieb, rs@nRunMn)
    # varianceSeries <- tsEvaNanRunningMean(varianceSeries, 
    #                                       ceiling(rs@nRunMn/2))
    # norm_trend <- rs@trendSeries/mean(rs@trendSeries, na.rm = TRUE)
    dtr1 = serieb - rs@trendSeries
    lneg = length(which(dtr1 < 0))
    # stab <- cor(nr, norm_trend, use = "pairwise.complete.obs")
    # if (iter == 1) 
    #   nr <- norm_trend
    if (lneg >= 1) 
      lnegs = c(lnegs, lneg)
    # sts <- c(sts, stab)
    pctd = c(pctd, pcts[iter])
  }
  
  rval = pctd[length(pctd)]
  if (sum(lnegs) > 1) {
    rval = pctd[which.min(lnegs)]
  }
  
  return(rval)
}

#######################  Arguments importation #############################
args <- commandArgs(TRUE)
argus=as.vector(unlist(strsplit(args, split = " ")))
Nsq = 47
tail = "high"
sce="Histo"
outlets="RNetwork"

# haz="flood"
var = "dis"

# workDir = "/BGFS/CLIMEX/tilloal/HydroMeteo/"
# setwd(workDir)
# hydroDir<-"/BGFS/CLIMEX/tilloal/HydroMeteo/Timeseries/dis6_Uncal2"

hydroDir<-("D:/tilloal/Documents/LFRuns_utils/data")
workDir<-("D:/tilloal/Documents/06_Floodrivers/dis")
setwd(workDir)
rspace= read.csv(paste0(hydroDir,"/subspace_efas.csv"))
rspace=rspace[,-1]
nrspace=rspace[Nsq,]
#outlets="Hybas09"
if (outlets=="Hybas07"){
  outletname="outletsv8_hybas07_01min"
  nameout="UCH07"
  outhybas=outletopen(hydroDir,outletname,nrspace)
}else if (outlets=="Hybas09"){
  outletname="outlets_hybas09_01min"
  nameout="UCH09"
  outhybas=outletopen(hydroDir,outletname,nrspace)
}else if (outlets=="RNetwork"){
  outletname="efas_rnet_100km_01min"
  nameout="UCRnet"
  outhybas=outletopen(hydroDir,outletname,nrspace)
  Idstart=as.numeric(Nsq)*10000
  if (length(outhybas$outlets)>0){
    outhybas$outlets=seq((Idstart+1),(Idstart+length(outhybas$outlets)))
  }
}

unikout=outhybas$outlets
outhybas$latlong=paste(round(outhybas$Var1,4),round(outhybas$Var2,4),sep=" ")


# file_out=paste0(workDir,"TSEVA/out/",sce,"/",haz,"/",season,"/parCat6h_",outlets,haz,"_",Nsq,Quarter,"_1951_2020.Rdata")
# fex=file.exists(file_out)
# if (fex==TRUE){
#   f.size=as.numeric(fs::file_size(file_out))
# }else{
#   f.size=0
# }

timeWindow = 365.25*30; #time windows in days, the correction is done within the functions


library(foreach)
library(doParallel)




#loading of the data for this scenario

#Load the file
#loading the files as netcdf (needs to be checked offline)

if (sce=="Histo"){
  filename=paste0("Timeseries/dis6_Calout/dis_",Nsq,"_1950_2020_cf")
}
if (sce=="SCF"){
  filename=paste0("Timeseries/dis6_CalCFSout/dis_",Nsq,"_1950_2020_cf")
}
if (sce=="RWCF"){
  filename=paste0("Timeseries/dis6_CalCFRWout/dis_",Nsq,"_1951_2020_rcf")
}
dists=disNcopenloc(filename,hydroDir,outhybas,1)
df.dis=dists 
print(paste0("opening square ", Nsq, " /88"))
timeStamps=unique(as.Date(df.dis$time,origin="1979-01-01"))
timeStamps=as.POSIXct(timeStamps-1/24)
txx=timeStamps
print(tail(txx))
length(txx)
df.dis$timeStamps=timeStamps

names(df.dis)[c(1,2)]=c("dis","outlets")

rmv=which(year(txx)==1950)
if (length(rmv)>1){
  df.dis=df.dis[-rmv,]
  txx=txx[-rmv]
}

# Register parallel backend
Station_data_IDs=unikout
cl <- makeCluster(4)
registerDoParallel(cl)
Trth_H_list=c()
tail="high"
Trth_H_list <- foreach(id = 1:10, .combine = rbind, .packages=c("RtsEva","zoo","ncdf4")) %dopar% {
  #print(id)
  stid <- Station_data_IDs[id]
  catch <- stid
  
  dists=disNcopenloc(filename,hydroDir,outhybas,id)
  df.dis=dists 
  
  if (length(rmv)>1){
    df.dis=df.dis[-rmv,]
  }
  data <- data.frame(txx,df.dis$outlets)
  names(data) <- c("date","Qs")
  # Q_H <- Q_simH[,id+1][-1]
  # Q_H <- Q_H[-rmv]
  
  if (tail=="high"){
    # Extract the maximum daily value
    series <- max_daily_value(data)
    timeAndSeriesH=series
    names(timeAndSeriesH)=c("timestamp","data")
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
    timeAndSeriesH=data.frame(timeStamps,series)
    names(timeAndSeriesH)=c("timestamp","data")
  }
  
  TrendTh_H1 <- tsEvaFindTrendThreshold(series=timeAndSeriesH$data, timeStamps=timeAndSeriesH$timestamp, timeWindow)
  
  if (is.null(TrendTh_H1))TrendTh_H1=NA
  
  TrendTh_H2 <- tsEvaFindTrendThreshold2(series=timeAndSeriesH$data, timeStamps=timeAndSeriesH$timestamp, timeWindow)
  if (is.null(TrendTh_H2))TrendTh_H2=NA
  Trth_H <- c(catch,TrendTh_H1,TrendTh_H2)
  names(Trth_H)=c("cid","Th_old","Th_new")
  
  Trth_H
}

stopCluster(cl)
write.csv(Trth_H_list,file=paste0("trenTH_",sce,"_",tail,"_",Nsq,".csv"))



# Trth_H_list=c()
# Trth_R_list=c()
# Trth_S_list=c()
# tail="low"
# for (id in 1:length(Station_data_IDs)){
#   
#   id=1
#   print(id)
#   stid=Station_data_IDs[id]
#   id=which(Station_data_IDs==394554)
#   catch=stid
#   Q_H=Q_simH[,id+1][-1]
#   Q_SCF=Q_simCFS[,id+1][-1]
#   Q_RCF=Q_simCFR[,id+1][-1]
#   
#   Q_H=Q_H[-rmv]
#   Q_SCF=Q_SCF[-rmv]
#   #Q_RCF=Q_RCF[-rmv]
#   
#   data=data.frame(txx,Q_H)
#   names(data)=c("date","Qs")
#   if (tail=="high"){
#     # Extract the maximum daily value
#     series <- max_daily_value(data)
#     timeAndSeriesH=series
#     names(timeAndSeriesH)=c("timestamp","data")
#   } else if (tail=="low"){
#     
#     minPeakDistanceInDays=30
#     #7 day average flow for drought analysis
#     WindowSize=7
#     names(data)=c("date","Qs")
#     dt1=min(diff(data$date),na.rm=T)
#     dt=as.numeric(dt1)
#     tdim=attributes(dt1)$units
#     if (tdim=="hours") dt=dt/24
#     nRunMn = ceiling(WindowSize/dt);
#     colnames(data)
#     data$Q7=rollmean(data$Qs,nRunMn, align = "right", fill=NA)
#     timeStamps=data$date
#     series=data$Q7
#     start_index=1
#     indices_to_extract <- seq(from = start_index, to = length(series), by = WindowSize/dt)
#     series=series[indices_to_extract]
#     series=-1*series
#     timeStamps=timeStamps[indices_to_extract]
#     timeAndSeriesH=data.frame(timeStamps,series)
#     names(timeAndSeriesH)=c("timestamp","data")
#   }
#   data=data.frame(txx,Q_RCF)
#   names(data)=c("date","Qs")
#   if (tail=="high"){
#     # Extract the maximum daily value
#     series <- max_daily_value(data)
#     timeAndSeriesR=series
#     names(timeAndSeriesR)=c("timestamp","data")
#   } else if (tail=="low"){
#     
#     minPeakDistanceInDays=30
#     #7 day average flow for drought analysis
#     WindowSize=7
#     names(data)=c("date","Qs")
#     dt1=min(diff(data$date),na.rm=T)
#     dt=as.numeric(dt1)
#     tdim=attributes(dt1)$units
#     if (tdim=="hours") dt=dt/24
#     nRunMn = ceiling(WindowSize/dt);
#     colnames(data)
#     data$Q7=rollmean(data$Qs,nRunMn, align = "right", fill=NA)
#     timeStamps=data$date
#     series=data$Q7
#     start_index=1
#     indices_to_extract <- seq(from = start_index, to = length(series), by = WindowSize/dt)
#     series=series[indices_to_extract]
#     series=-1*series
#     timeStamps=timeStamps[indices_to_extract]
#     timeAndSeriesR=data.frame(timeStamps,series)
#     names(timeAndSeriesR)=c("timestamp","data")
#   }
# 
#   data=data.frame(txx,Q_SCF)
#   names(data)=c("date","Qs")
#   if (tail=="high"){
#     # Extract the maximum daily value
#     series <- max_daily_value(data)
#     timeAndSeriesS=series
#     names(timeAndSeriesS)=c("timestamp","data")
#   } else if (tail=="low"){
#     
#     minPeakDistanceInDays=30
#     #7 day average flow for drought analysis
#     WindowSize=7
#     names(data)=c("date","Qs")
#     dt1=min(diff(data$date),na.rm=T)
#     dt=as.numeric(dt1)
#     tdim=attributes(dt1)$units
#     if (tdim=="hours") dt=dt/24
#     nRunMn = ceiling(WindowSize/dt);
#     colnames(data)
#     data$Q7=rollmean(data$Qs,nRunMn, align = "right", fill=NA)
#     timeStamps=data$date
#     series=data$Q7
#     start_index=1
#     indices_to_extract <- seq(from = start_index, to = length(series), by = WindowSize/dt)
#     series=series[indices_to_extract]
#     series=-1*series
#     timeStamps=timeStamps[indices_to_extract]
#     timeAndSeriesS=data.frame(timeStamps,series)
#     names(timeAndSeriesS)=c("timestamp","data")
#   }
# 
#   
#   plot(timeAndSeriesS)
#   TrendTh_H1=tsEvaFindTrendThreshold(series=timeAndSeriesH$data, timeStamps=timeAndSeriesH$timestamp, timeWindow)
#   TrendTh_R1=tsEvaFindTrendThreshold(series=timeAndSeriesR$data, timeStamps=timeAndSeriesR$timestamp, timeWindow)
#   TrendTh_S1=tsEvaFindTrendThreshold(series=timeAndSeriesS$data, timeStamps=timeAndSeriesS$timestamp, timeWindow)
#   
#   TrendTh_H2=tsEvaFindTrendThreshold2(series=timeAndSeriesH$data, timeStamps=timeAndSeriesH$timestamp, timeWindow)
#   TrendTh_R2=tsEvaFindTrendThreshold2(series=timeAndSeriesR$data, timeStamps=timeAndSeriesR$timestamp, timeWindow)
#   TrendTh_S2=tsEvaFindTrendThreshold2(series=timeAndSeriesS$data, timeStamps=timeAndSeriesS$timestamp, timeWindow)
#   
#   Trth_H=c(catch,TrendTh_H1,TrendTh_H2)
#   Trth_R=c(catch,TrendTh_R1,TrendTh_R2)
#   Trth_S=c(catch,TrendTh_S1,TrendTh_S2)
#   
#   Trth_H_list=rbind(Trth_H_list,Trth_H)
#   Trth_R_list=rbind(Trth_R_list,Trth_R)
#   Trth_S_list=rbind(Trth_S_list,Trth_S)
#   
# }
# 