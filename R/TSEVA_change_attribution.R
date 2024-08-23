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

#  Function declaration ---------------------------------------------------
source("D:/tilloal/Documents/LFRuns_utils/TrendAnalysis/R/function_loading.R")
# Load inputs from HPC computation ----------------------------------------




# Pre-loaded results ------------------------------------------------------



## Spatial data for catchments ----
#Set data directory
hydroDir<-("D:/tilloal/Documents/LFRuns_utils/data")
hybas07 <- read_sf(dsn = paste0(hydroDir,"/Catchments/hydrosheds/hybas_eu_lev07_v1c.shp"))
hybasf7=fortify(hybas07) 

#load my shitty rdata
load(file=paste0(hydroDir,"/dataHybas07_ids.Rdata"))
outletname="outletsv8_hybas07_01min"
nameout="UCH07"

### Hybas07 ----
Catchmentrivers7=read.csv(paste0(hydroDir,"/Catchments/from_hybas_eu_onlyid.csv"),encoding = "UTF-8", header = T, stringsAsFactors = F)
hybas07 <- read_sf(dsn = paste0(hydroDir,"/Catchments/hydrosheds/hybas_eu_lev07_v1c.shp"))
hybasf7=fortify(hybas07) 
Catamere07=inner_join(hybasf7,Catchmentrivers7,by= "HYBAS_ID")
outf=outletopen(hydroDir,outletname)
outf$latlong=paste(round(outf$Var1,4),round(outf$Var2,4),sep=" ")
Catamere07$llcoord=paste(round(Catamere07$POINT_X,4),round(Catamere07$POINT_Y,4),sep=" ") 
Catf7=inner_join(Catamere07,outf,by= c("llcoord"="latlong"))
cst7=st_transform(Catf7,  crs=3035)

# zebi=unique(parlist$catchment)
# parmerde=parlist[(match(zebi,parlist$catchment)),]
# outcut=which(!is.na(match(outf$outlets,zebi)))
# outhloc=outf[outcut,]
# Catfchier=inner_join(Catamere07,outhloc,by= c("llcoord"="latlong"))

reservoirdog=outletopen(paste0(hydroDir,"/reservoirs"),"res_ratio_diff_2020-1951")
reservoirdog$llcoord=paste(round(reservoirdog$Var1,4),round(reservoirdog$Var2,4),sep=" ") 

#other method to extract right catchments
#efas rnt
out1=outletopen(hydroDir,"efas_rnet_100km_01min")
out1$latlong=paste(round(out1$Var1,4),round(out1$Var2,4),sep=" ")
out2=inner_join(out1,outf, by="latlong")


### Plot parameters ----
palet2=c(hcl.colors(9, palette = "Blues", alpha = NULL, rev = TRUE, fixup = TRUE))
outll=outletopen(hydroDir,outletname)
cord.dec=outll[,c(2,3)]
cord.dec = SpatialPoints(cord.dec, proj4string=CRS("+proj=longlat"))
cord.UTM <- spTransform(cord.dec, CRS("+init=epsg:3035"))
nco=cord.UTM@coords
coord.lamb=spTransform(cord.dec, CRS("+proj=laea +lat_0=55 +lon_0=20 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=km +no_defs"))
nci=coord.lamb@coords
world <- ne_countries(scale = "medium", returnclass = "sf")
Europe <- world[which(world$continent == "Europe"),]
e2=st_transform(Europe,  crs=3035)
w2=st_transform(world,  crs=3035)
tsize=12
osize=12
catchmentproj=cst7
Impdates=seq(1951,2020,by=1)
valuenames=paste0("X",Impdates)
catmap=catchmentproj
basemap=w2



#-------------Experiment: loading datasets and fitting TSEVA--------------------
#load the three dataset together and fit TSEVA on same peaks?

main_path = "D:/tilloal/Documents/06_Floodrivers/"
valid_path = paste0(main_path)

dis_path<-paste0(main_path,'HERA_CFS/')
load(file=paste0(dis_path,"HERA_CFS6h_H07_19502020.Rdata"))
Station_data_IDs <- as.vector(t(Q_sim[1, ]))[-1]
tqx=as.POSIXct((Q_sim[,1]*3600*24), origin="1979-01-01 00:00:00")
txx=tqx[-1]
Q_simCFS=Q_sim

dis_path<-paste0(main_path,'HERA/')
load(file=paste0(dis_path,"HERA6h_H07_19502020.Rdata"))
Q_simH=Q_sim

dis_path<-paste0(main_path,'HERA_RWstat/')
load(file=paste0(dis_path,"HERA_CFR6h_H07_19502020.Rdata"))
Q_simCFR=Q_sim

#extract values for 1 location

#remove locations without discharge
wolala=which(is.na(Q_simCFR[2,]))-1
Q_simH=Q_simH[,-which(is.na(Q_simH[2,]))]
Q_simCFR=Q_simCFR[,-which(is.na(Q_simCFR[2,]))]
Q_simCFS=Q_simCFS[,-which(is.na(Q_simCFS[2,]))]
Station_data_IDs=Station_data_IDs[-wolala]

#match HYBAS_id with station_IDs

hybmatch=match(Station_data_IDs,cst7$outlets)

hybas_IDs=cst7$HYBAS_ID[hybmatch]

stid=382509
catch=stid

id=which(Station_data_IDs==stid)
Q_H=Q_simH[,id+1][-1]
Q_SCF=Q_simCFS[,id+1][-1]
Q_RCF=Q_simCFR[,id+1][-1]


txx=tqx[-1]
rmv=which(year(txx)==1950)
txx=txx[-rmv]
Q_H=Q_H[-rmv]
Q_SCF=Q_SCF[-rmv]


plot(Q_SCF)
#TSEVA for each scenario
catmat=Catf7[which(Catf7$pointid==stid),]
minPeakDistanceInDays=7
interflag=0
frosttime=NA
trans="rev"
tail="low"



data=data.frame(txx,Q_H)
names(data)=c("date","Qs")
# Extract the maximum daily value
series <- max_daily_value(data)
names(series)=c("timestamp","dis")


dt1=min(diff(series$timestamp),na.rm=T)
dt=as.numeric(dt1)
tdim=attributes(dt1)$units
if (tdim=="hours") dt=dt/24
if (dt==1){
  timeDays=series$timestamp
}else{
  timeDays=unique(as.Date(series$timestamp))
}

bounds=c(year(timeDays[1]),year(timeDays[length(timeDays)]))
realbound=bounds
tbound=c(as.Date(paste0(realbound[1],"-12-31")),as.Date(paste0(realbound[2],"-12-31")))
Impdates=seq(tbound[1],tbound[2],by="1 year")
tindexes=match(Impdates,timeDays)
names(series)=c("timestamp","dis")
nv=length(unique(series$dis))
timeAndSeries=series
names(timeAndSeries)=c("timestamp","data")

tsm=1/dt
series=timeAndSeries[,2]
timeWindow = 365.25*30; #time windows in days, the correction is done within the functions
windowSize=366
timeStamps=timeAndSeries$timestamp


reservoir_bg=inner_join(reservoirdog,outf,by=c("llcoord"="latlong"))
reservoir_loc=reservoir_bg[match(catch,reservoir_bg$outlets.y),]
Reservoir_alteration=reservoir_loc$outlets.x

#Keep the same trend threshold for other scenarios

#Ok restart from the beginning
data=data.frame(txx,Q_H)
names(data)=c("date","Qs")
# Extract the maximum daily value
series <- max_daily_value(data)
timeAndSeriesH=series
names(timeAndSeriesH)=c("timestamp","data")

data=data.frame(txx,Q_RCF)
names(data)=c("date","Qs")
# Extract the maximum daily value
series <- max_daily_value(data)
timeAndSeriesR=series
names(timeAndSeriesR)=c("timestamp","data")

data=data.frame(txx,Q_SCF)
names(data)=c("date","Qs")
# Extract the maximum daily value
series <- max_daily_value(data)
timeAndSeriesS=series
names(timeAndSeriesS)=c("timestamp","data")



check_timeserie2=function(timeseries,yro){
  ts_years <- as.integer((lubridate::year(timeseries)))
  year_check <- yro %in% ts_years
  runs <- rle(year_check)
  rf=which(runs$values==FALSE)
  if (any(runs$lengths[rf] >= 2)) {
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
    varianceSeries <- tsEvaNanRunningVariance(serieb, rs@nRunMn)
    varianceSeries <- tsEvaNanRunningMean(varianceSeries, 
                                          ceiling(rs@nRunMn/2))
    norm_trend <- rs@trendSeries/mean(rs@trendSeries, na.rm = TRUE)
    dtr1 = serieb - rs@trendSeries
    lneg = length(which(dtr1 < 0))
    stab <- cor(nr, norm_trend, use = "pairwise.complete.obs")
    if (iter == 1) 
      nr <- norm_trend
    if (lneg >= 1) 
      lnegs = c(lnegs, lneg)
    sts <- c(sts, stab)
    pctd = c(pctd, pcts[iter])
  }
  sts[1] <- 1
  plot(sts)
  cptm <- try(changepoint::cpt.var(sts, method = "AMOC", 
                                   penalty = "BIC", minseglen = 3), T)
  if (inherits(cptm, "try-error")) {
    rval = pctd[which.max(sts[-1])]
  } else if (cptm@cpts[1] == length(sts)) {
    message("no change point")
    cptm@cpts[1] <- length(sts)/2
    rval = pctd[cptm@cpts[1]]
  }
  if (sum(lnegs) > 1) {
    rval = pctd[which.min(lnegs)]
  }
  
  return(rval)
}

TrendTh=tsEvaFindTrendThreshold2(series=timeAndSeriesH$data, timeStamps=timeAndSeriesS$timestamp, timeWindow)

print(TrendTh)

Nonstat<-TsEvaNs(timeAndSeriesH, timeWindow, transfType="trendPeaks", minPeakDistanceInDays = minPeakDistanceInDays, tail=tail, trans=trans,TrendTh=TrendTh)
nonStationaryEvaParamsH=Nonstat[[1]]
stationaryTransformDataH=Nonstat[[2]]

TrendTh=tsEvaFindTrendThreshold2(series=timeAndSeriesR$data, timeStamps=timeAndSeriesS$timestamp, timeWindow)

Nonstat<-TsEvaNs(timeAndSeriesR, timeWindow, transfType="trendPeaks", minPeakDistanceInDays = minPeakDistanceInDays, tail=tail, trans=trans,TrendTh=TrendTh)
nonStationaryEvaParamsR=Nonstat[[1]]
stationaryTransformDataR=Nonstat[[2]]

TrendTh=tsEvaFindTrendThreshold2(series=timeAndSeriesS$data, timeStamps=timeAndSeriesS$timestamp, timeWindow)

Nonstat<-TsEvaNs(timeAndSeriesS, timeWindow, transfType="trendPeaks", minPeakDistanceInDays = minPeakDistanceInDays, tail=tail, trans=trans, TrendTh=TrendTh)
nonStationaryEvaParamsS=Nonstat[[1]]
stationaryTransformDataS=Nonstat[[2]]



stationaryTransformDataH$timeStampsDay=unique(as.Date(stationaryTransformDataH$timeStamps))
pikosH=data.frame(nonStationaryEvaParamsH$potObj$parameters$peaks,nonStationaryEvaParamsH$potObj$parameters$peakID,nonStationaryEvaParamsH$potObj$parameters$peakST, 
                  nonStationaryEvaParamsH$potObj$parameters$peakEN)
names(pikosH)=c("value","timeID","tIDstart","tIDend")
pikosH$time=stationaryTransformDataH$timeStamps[pikosH$timeID]
pikosH$catch=rep(catch,length(pikosH[,1]))

stationaryTransformDataR$timeStampsDay = unique(as.Date(stationaryTransformDataR$timeStamps))
pikosR = data.frame(nonStationaryEvaParamsR$potObj$parameters$peaks,
                    nonStationaryEvaParamsR$potObj$parameters$peakID,
                    nonStationaryEvaParamsR$potObj$parameters$peakST,
                    nonStationaryEvaParamsR$potObj$parameters$peakEN)
names(pikosR) = c("value", "timeID", "tIDstart", "tIDend")
pikosR$time = stationaryTransformDataR$timeStamps[pikosR$timeID]
pikosR$catch = rep(catch, length(pikosR[,1]))

stationaryTransformDataS$timeStampsDay = unique(as.Date(stationaryTransformDataS$timeStamps))
pikosS = data.frame(nonStationaryEvaParamsS$potObj$parameters$peaks,
                    nonStationaryEvaParamsS$potObj$parameters$peakID,
                    nonStationaryEvaParamsS$potObj$parameters$peakST,
                    nonStationaryEvaParamsS$potObj$parameters$peakEN)
names(pikosS) = c("value", "timeID", "tIDstart", "tIDend")
pikosS$time = stationaryTransformDataS$timeStamps[pikosS$timeID]
pikosS$catch = rep(catch, length(pikosS[,1]))

plot(pikosS$value)
#Here I need to convert the timeStamp to a daily one if dt is not 1
dt1=min(diff(timeStamps),na.rm=T)
dt=as.numeric(dt1)
tdim=attributes(dt1)$units
if (tdim=="hours") dt=dt/24
if (dt==1){
  timeDays=stationaryTransformDataH$timeStamps
}else{
  timeDays=stationaryTransformDataH$timeStampsDay
}

bounds=c(year(timeDays[1]),year(timeDays[length(timeDays)]))
#realbound=10*ceiling(bounds/10)
realbound=bounds
tbound=c(as.Date(paste0(realbound[1],"-12-31")),as.Date(paste0(realbound[2],"-12-31")))
Impdates=seq(tbound[1],tbound[2],by="1 years")
datex=yday(timeDays)
dtect=c(diff(datex),-1)
last_days <- timeDays[which(dtect<0)]
tindexes=match(last_days,timeDays)

# Compute return periods and levels
RPgoal=10
timeIndex=tindexes[1]
RLevs100H=ComputeReturnLevels(nonStationaryEvaParamsH, RPgoal, timeIndex)
RLevs100R=ComputeReturnLevels(nonStationaryEvaParamsR, RPgoal, timeIndex)
RLevs100S=ComputeReturnLevels(nonStationaryEvaParamsS, RPgoal, timeIndex)

RLevs100H
RLevs100R
RLevs100S

paramsS=c(catch,year(Impdates[1]),timeIndex,RLevs100S$Params)
names(paramsS)[1:3]=c("catchment","Year","timeIndex")

RLgevH=RLevs100H$ReturnLevels[2]
RLgpdH=RLevs100H$ReturnLevels[3]
ERgevH=RLevs100H$ReturnLevels[4]
ERgpdH=RLevs100H$ReturnLevels[5]
nRPgevH=nRPgpdH=10
paramsH=c()
for (t in 2:length(Impdates)){
  timeIndex=tindexes[t]
  RLevs100i=ComputeReturnLevels(nonStationaryEvaParamsH, RPgoal, timeIndex)
  params=c(catch,year(Impdates[t]),timeIndex,RLevs100i$Params)
  names(params)[1:3]=c("catchment","Year","timeIndex")
  
  Rper=RPcalc(params,RPiGEV=RLevs100H$ReturnLevels[2],RPiGPD=RLevs100H$ReturnLevels[3])
  nRPgpdH=c(nRPgpdH,Rper[2])
  nRPgevH=c(nRPgevH,Rper[1])
  RLgevH=cbind(RLgevH,RLevs100i$ReturnLevels[2])
  RLgpdH=cbind(RLgpdH,RLevs100i$ReturnLevels[3])
  ERgevH=cbind(ERgevH,RLevs100i$ReturnLevels[4])
  ERgpdH=cbind(ERgpdH,RLevs100i$ReturnLevels[5])
  if (length(parlist)>1) colnames(parlist)=names(params)
  parlist=rbind(parlist,params)
}

RLgevH=as.data.frame(RLgevH)
names(RLgevH)=year(Impdates)
rownames(RLgevH)=RPgoal

RLgpdH=as.data.frame(RLgpdH)
names(RLgpdH)=year(Impdates)
rownames(RLgpdH)=RPgoal

nRPgevH=as.data.frame(t(nRPgevH))
names(nRPgevH)=year(Impdates)

nRPgpdH=as.data.frame(t(nRPgpdH))
names(nRPgpdH)=year(Impdates)




RLgevR=RLevs100R$ReturnLevels[2]
RLgpdR=RLevs100R$ReturnLevels[3]
ERgevR=RLevs100R$ReturnLevels[4]
ERgpdR=RLevs100R$ReturnLevels[5]
nRPgevR=nRPgpdR=10
paramsR=c()
for (t in 2:length(Impdates)){
  timeIndex=tindexes[t]
  RLevs100i=ComputeReturnLevels(nonStationaryEvaParamsR, RPgoal, timeIndex)
  params=c(catch,year(Impdates[t]),timeIndex,RLevs100i$Params)
  names(params)[1:3]=c("catchment","Year","timeIndex")
  
  Rper=RPcalc(params,RPiGEV=RLevs100R$ReturnLevels[2],RPiGPD=RLevs100R$ReturnLevels[3])
  nRPgpdR=c(nRPgpdR,Rper[2])
  nRPgevR=c(nRPgevR,Rper[1])
  RLgevR=cbind(RLgevR,RLevs100i$ReturnLevels[2])
  RLgpdR=cbind(RLgpdR,RLevs100i$ReturnLevels[3])
  ERgevR=cbind(ERgevR,RLevs100i$ReturnLevels[4])
  ERgpdR=cbind(ERgpdR,RLevs100i$ReturnLevels[5])
  if (length(parlist)>1) colnames(parlist)=names(params)
  parlist=rbind(parlist,params)
}

RLgevR=as.data.frame(RLgevR)
names(RLgevR)=year(Impdates)
rownames(RLgevR)=RPgoal

RLgpdR=as.data.frame(RLgpdR)
names(RLgpdR)=year(Impdates)
rownames(RLgpdR)=RPgoal

nRPgevR=as.data.frame(t(nRPgevR))
names(nRPgevR)=year(Impdates)

nRPgpdR=as.data.frame(t(nRPgpdR))
names(nRPgpdR)=year(Impdates)


parlistS=paramsS

RLgevS=RLevs100S$ReturnLevels[2]
RLgpdS=RLevs100S$ReturnLevels[3]
ERgevS=RLevs100S$ReturnLevels[4]
ERgpdS=RLevs100S$ReturnLevels[5]
nRPgevS=nRPgpdS=10
for (t in 2:length(Impdates)){
  timeIndex=tindexes[t]
  RLevs100i=ComputeReturnLevels(nonStationaryEvaParamsS, RPgoal, timeIndex)
  paramsS=c(catch,year(Impdates[t]),timeIndex,RLevs100i$Params)
  names(paramsS)[1:3]=c("catchment","Year","timeIndex")
  
  Rper=RPcalc(params,RPiGEV=RLevs100S$ReturnLevels[2],RPiGPD=RLevs100S$ReturnLevels[3])
  nRPgpdS=c(nRPgpdS,Rper[2])
  nRPgevS=c(nRPgevS,Rper[1])
  RLgevS=cbind(RLgevS,RLevs100i$ReturnLevels[2])
  RLgpdS=cbind(RLgpdS,RLevs100i$ReturnLevels[3])
  ERgevS=cbind(ERgevS,RLevs100i$ReturnLevels[4])
  ERgpdS=cbind(ERgpdS,RLevs100i$ReturnLevels[5])
  if (length(parlist)>1) colnames(parlist)=names(params)
  parlistS=rbind(parlistS,paramsS)
}

RLgevS=as.data.frame(RLgevS)
names(RLgevS)=year(Impdates)
rownames(RLgevS)=RPgoal

RLgpdS=as.data.frame(RLgpdS)
names(RLgpdS)=year(Impdates)
rownames(RLgpdS)=RPgoal

nRPgevS=as.data.frame(t(nRPgevS))
names(nRPgevS)=year(Impdates)

nRPgpdS=as.data.frame(t(nRPgpdS))
names(nRPgpdS)=year(Impdates)


parlistS=data.frame(parlistS)
length(t(RLgpdH))
length(parlistS$timeIndex)
plot(nonStationaryEvaParamsH$gevObj$parameters$annualMaxIndx,nonStationaryEvaParamsH$gevObj$parameters$annualMax,type="p",pch=16)
lines(nonStationaryEvaParamsR$gevObj$parameters$annualMaxIndx,nonStationaryEvaParamsR$gevObj$parameters$annualMax,type="p",col="blue",pch=16)
lines(nonStationaryEvaParamsS$gevObj$parameters$annualMaxIndx,nonStationaryEvaParamsS$gevObj$parameters$annualMax,type="p",col="red",pch=16)

lines(parlistS$timeIndex , t(RLgpdH))
lines(parlistS$timeIndex,t(RLgpdR),col="blue")
lines(parlistS$timeIndex,t(RLgpdS),col="red")

RLGPDH=as.vector(t(RLgpdH))
RLGPDR=as.vector(t(RLgpdR))
RLGPDS=as.vector(t(RLgpdS))

#difference attributed to reservoir konstraktion
#Reservoir changes
DRes=(RLGPDH-RLGPDR)/RLGPDH[1]*100
plot(DRes)

#Socioeconomic changes
DSoc=(RLGPDR-RLGPDS)/RLGPDH[1]*100
plot(DSoc)

#Climate change
DClim=(RLGPDS-RLGPDS[1])/RLGPDH[1]*100


#check other approaches

#1 I remove the climate trend
RLGPDSdt=RLGPDS-(RLGPDS-RLGPDS[1])
RLGPDRdt=RLGPDR-(RLGPDS-RLGPDS[1])
RLGPDHdt=RLGPDH-(RLGPDS-RLGPDS[1])
plot(RLGPDSdt)

#Now I remove the socioeconomic trend
RLGPDRdtp=RLGPDRdt-(RLGPDRdt-RLGPDSdt)
RLGPDHdtp=RLGPDHdt-(RLGPDRdt-RLGPDSdt)
points(RLGPDRdtp,col="purple")
points(RLGPDHdtp,col="red")
#And the reservoir
RLGPDHdtpx=RLGPDHdtp-(RLGPDHdtp-RLGPDRdtp)
points(RLGPDHdtpx,col="cyan")

#change from climate
Climtrend=(RLGPDS-RLGPDS[1])
plot(Climtrend,ylim=c(-100,100))
#change from socioeco
#Soctrend=(RLGPDRdt-RLGPDRdt[1])+(RLGPDRdt-RLGPDSdt)
Soctrend=(RLGPDRdt-RLGPDSdt)
points(Soctrend,col=2)


#change from reservoirs
Restrend=(RLGPDHdtp-RLGPDRdtp)
#Restrend=(RLGPDHdtp-RLGPDHdtp[1])+(RLGPDHdtp-RLGPDRdtp)
points(Restrend,col="blue")

df <- data.frame(
  time = c(year(timeDays[1]):year(timeDays[length(timeDays)])),
  Dres = Restrend/RLGPDS[1]*100,
  Dsoc = Soctrend/RLGPDS[1]*100,
  Dclim=Climtrend/RLGPDS[1]*100
)

df_long <- reshape2::melt(df, id.vars = "time", variable.name = "driver", value.name = "value")

ggplot(df_long, aes(x = time, y = value, fill = driver)) +
  geom_area() +
  theme_minimal() +
  labs(title = "10y RL ",
       x = "Time",
       y = "Value",
       fill = "driver")

total=Restrend+Soctrend+Climtrend
TotalChange=(RLGPDH-RLGPDS[1])

df$catch=rep(stid,length(Restrend))

#OK that would be an output plot for each catchment


# Plotting the outputs ----------------------------------------------------


plot_riverchange<-function(plot_inputs,main,dates){
  
  yr.deb <-  seq(as.POSIXct("1950-01-15"), by="5 years", length=14)
  Climtrend=plot_inputs$Climtrend
  Soctrend=plot_inputs$Soctrend
  Restrend=plot_inputs$Restrend
  RLGPDH=plot_inputs$RLGPDH
  RLGPDS=plot_inputs$RLGPDS
  RLGPDR=plot_inputs$RLGPDR
  pid=plot_inputs$pid
  prid=plot_inputs$prid
  psid=plot_inputs$psid
  
  ctrend=c(Climtrend)
  lctrend=c(Climtrend+Soctrend)
  rwutrend=c(lctrend+Restrend)
  
  ptrend=c(lctrend,rev(ctrend))
  rtrend=c(rwutrend,rev(lctrend))
  pdates=c(dates,rev(dates))
  qlim=c(0, 1.2*max(c(pid$value,prid$value,psid$value)))
  
  plot(pid$time, pid$value, col=alpha("grey",.8) ,pch=16,axes=FALSE,xaxs="i",yaxs="i",ylim=qlim,
       xlab = NA, ylab="")
  points(psid$time,psid$value,col=alpha("royalblue",.8),type="p",pch=16)
  points(prid$time,prid$value,col=alpha("orange",.8),type="p",pch=16)
  
  lines(dates,RLGPDH,col="lightgrey",lwd=2)
  lines(dates,RLGPDS,col="royalblue",lwd=2)
  lines(dates,RLGPDR,col="orange",lwd=2)
  mtext(main,3,font = 2,line = 0.5,cex = 1.5)
  abline(v = yr.deb, col="lightgrey", lty=2)
  abline(h = 0, col=1)
  axis(2,cex.axis=1)
  title(ylab = expression(paste("Q (",m^3/s,")")),cex.lab=1.5,line=2, xlab="years")
  axis(1, yr.deb, label=format(yr.deb,"%Y"),cex.axis=1)
  box()
  ## Trace des erreurs absolues
  polygon(c(dates[1],dates,dates[length(dates)]),c(0,abs(Climtrend),0),
          col=alpha("royalblue",.5),border="royalblue")
  polygon(pdates,abs(ptrend),
          col=alpha("orange",.5),border="orange")
  polygon(pdates,abs(rtrend),
          col=alpha("grey",.5),border="lightgrey")
  ## Trace de la grille mensuelle
  
  ## Definition de la legende
  legend("topleft", leg=c("Historical 10y RL","RW static 10y RL","Socio cf 10y RL"),
         lwd=c(2,2,2), col=c("grey","orange","royalblue"),
         cex=1, lty=c(1,1,1), bg=alpha("white",.6))
} 


#load upstream area
main_path = 'D:/tilloal/Documents/06_Floodrivers/'
valid_path = paste0(main_path,'DataPaper/')
outletname="/GIS/upArea_European_01min.nc"
dir=valid_path
outf$idlalo=paste(outf$idlo, outf$idla, sep=" ")
UpArea=UpAopen(valid_path,outletname,outf)
head(UpArea)


## 1. selection of data to be plotted ----
#Loading saved results in .Rdata

tail="low"

if (tail=="high"){
  #load(file=paste0(hydroDir,"/TSEVA_hybas/outputs/parCat07_Histo_flood_1950_2020.Rdata"))
  load(file=paste0(hydroDir,"/TSEVA_hybas/outputs/ResCat07_Histo_floodv2_1951_2020.Rdata"))
  Flood_histo=Results
  
  #load(file=paste0(hydroDir,"/TSEVA_hybas/outputs/parCat07_SocCF_flood_1950_2020.Rdata"))
  load(file=paste0(hydroDir,"/TSEVA_hybas/outputs/ResCat07_SocCF_floodv2_1951_2020.Rdata"))
  Flood_SocCF=Results
  
  
  #load(file=paste0(hydroDir,"/TSEVA_hybas/outputs/parCat07_ResCF_flood_1950_2020.Rdata"))
  load(file=paste0(hydroDir,"/TSEVA_hybas/outputs/ResCat07_RWstat_floodv2_1951_2020.Rdata"))
  Flood_ResCF=Results
  
  RLGPD_SCF=Flood_SocCF$RetLevGPD
  RLGPD_SCF=data.frame(RLGPD_SCF,unikout=Flood_SocCF$catrest$catlist)
  Peaks_SCF=Flood_SocCF$Peaks
  
  
  RLGPD_H=Flood_histo$RetLevGPD
  RLGPD_H=data.frame(RLGPD_H,unikout=Flood_histo$catrest$catlist)
  Peaks_H=Flood_histo$Peaks
  
  
  RLGPD_RCF=Flood_ResCF$RetLevGPD
  RLGPD_RCF=data.frame(RLGPD_RCF,unikout=Flood_ResCF$catrest$catlist)
  Peaks_RCF=Flood_ResCF$Peaks
  
  Station_data_IDs=RLGPD_H$unikout
}

if (tail=="low"){
  #load(file=paste0(hydroDir,"/TSEVA_hybas/outputs/parCat07_Histo_flood_1950_2020.Rdata"))
  load(file=paste0(hydroDir,"/TSEVA_hybas/outputs/ResCat07_Histo_droughtv2_1951_2020.Rdata"))
  Drought_histo=Results
  IRES_H=Drought_histo$catrest
  
  
  #load(file=paste0(hydroDir,"/TSEVA_hybas/outputs/parCat07_SocCF_flood_1950_2020.Rdata"))
  load(file=paste0(hydroDir,"/TSEVA_hybas/outputs/ResCat07_SocCF_Droughtv2_1951_2020.Rdata"))
  Drought_SocCF=Results
  IRES_SCF=Drought_SocCF$catrest
  
  #load(file=paste0(hydroDir,"/TSEVA_hybas/outputs/parCat07_ResCF_Drought_1950_2020.Rdata"))
  load(file=paste0(hydroDir,"/TSEVA_hybas/outputs/ResCat07_RWstat_Droughtv2_1951_2020.Rdata"))
  Drought_ResCF=Results
  IRES_RCF=Drought_ResCF$catrest
  
  length(which(IRES_H$IRES==1))
  length(which(IRES_RCF$IRES==1))
  length(which(IRES_SCF$IRES==1))
  
  IRES_all=inner_join(IRES_H,IRES_RCF,by="catlist")
  IRES_all=inner_join(IRES_all,IRES_SCF,by="catlist")
  IRES_all$tor=IRES_all$IRES.y+IRES_all$IRES+IRES_all$IRES.x
  
  rmIRES=which(IRES_all$tor>0)
  IRES_all$catlist[rmIRES]
  RLGPD_SCF=Drought_SocCF$RetLevGPD
  RLGPD_SCF=data.frame(RLGPD_SCF,unikout=Drought_SocCF$catrest$catlist)
  Peaks_SCF=Drought_SocCF$Peaks
  
  
  RLGPD_H=Drought_histo$RetLevGPD
  RLGPD_H=data.frame(RLGPD_H,unikout=Drought_histo$catrest$catlist)
  Peaks_H=Drought_histo$Peaks
  
  
  RLGPD_RCF=Drought_ResCF$RetLevGPD
  RLGPD_RCF=data.frame(RLGPD_RCF,unikout=Drought_ResCF$catrest$catlist)
  Peaks_RCF=Drought_ResCF$Peaks
  
  Station_data_IDs=RLGPD_H$unikout[-rmIRES]
}

bigfish=Drought_histo$RetPerGPD
bigfish=data.frame(bigfish[,70],unikout=Drought_ResCF$catrest$catlist)
#reproduce what was done before for all catchments


a=RLGPD_SCF$X1956-RLGPD_RCF$X1956
plot(a)
#look at a pixel with high negative change

#create a time vector with years
dates <- seq.Date(ymd("1951-06-01"), ymd("2020-06-01"), by = "year")
dates= as.POSIXct(dates, format = "%Y-%m-%d")
dfall=c()
plotOut=FALSE
for (id in 1:length(Station_data_IDs)){
  #id=733
  #id=23
  print(id)
  stid=Station_data_IDs[id]
  stid=202512
  catch=stid
  
  
  RLGPDH=as.vector(t(RLGPD_H[which(RLGPD_H$unikout==stid),-71]))
  RLGPDR=as.vector(t(RLGPD_RCF[which(RLGPD_RCF$unikout==stid),-71]))
  RLGPDS=as.vector(t(RLGPD_SCF[which(RLGPD_SCF$unikout==stid),-71]))
  
  
  pid=Peaks_H[which(Peaks_H$catch==stid),]
  psid=Peaks_SCF[which(Peaks_SCF$catch==stid),]
  prid=Peaks_RCF[which(Peaks_RCF$catch==stid),]
  
  
  if (tail=="low"){
    RLGPDH[which(RLGPDH<0)]=0
    RLGPDS[which(RLGPDS<0)]=0
    RLGPDR[which(RLGPDR<0)]=0
    RLGPDH[which(is.na(RLGPDH))]= RLGPDH[which(is.na(RLGPDH))[1]-1]
    pid$value=-pid$value
    prid$value=-prid$value
    psid$value=-psid$value
  }

  
  #1 I remove the climate trend
  RLGPDSdt=RLGPDS-(RLGPDS-RLGPDS[1])
  RLGPDRdt=RLGPDR-(RLGPDS-RLGPDS[1])
  RLGPDHdt=RLGPDH-(RLGPDS-RLGPDS[1])

  #Now I remove the socioeconomic trend (land use)
  RLGPDRdtp=RLGPDRdt-(RLGPDRdt-RLGPDSdt)
  RLGPDHdtp=RLGPDHdt-(RLGPDRdt-RLGPDSdt)

  #And the reservoir + water use
  RLGPDHdtpx=RLGPDHdtp-(RLGPDHdtp-RLGPDRdtp)

  
  #change from climate
  Climtrend=(RLGPDS-RLGPDS[1])

  #change from socioeco
  #Soctrend=(RLGPDRdt-RLGPDRdt[1])+(RLGPDRdt-RLGPDSdt)
  Soctrend=(RLGPDRdt-RLGPDSdt)

  
  
  #change from reservoirs
  Restrend=(RLGPDHdtp-RLGPDRdtp)
  #Restrend=(RLGPDHdtp-RLGPDHdtp[1])+(RLGPDHdtp-RLGPDRdtp)
  
  if (plotOut==TRUE)
    {
    Cairo::Cairo(
      20, #length
      15, #width
      file = paste("D:/tilloal/Documents/LFRuns_utils/TrendAnalysis/plots/",catch, "flood.png", sep = ""),
      type = "png", #tiff
      bg = "white", #white or transparent depending on your requirement 
      dpi = 300,
      units = "cm" #you can change to pixels etc 
    )
    
  # create a function from this  
    main="Moselle @ Koblenz"
    
    
    plot_inputs=list(Climtrend=Climtrend,Soctrend=Soctrend,Restrend=Restrend,
                     RLGPDH=RLGPDH,RLGPDR=RLGPDR,RLGPDS=RLGPDS,
                     pid=pid,psid=psid,prid=prid)
    
    
    test=plot_riverchange(plot_inputs,main,dates)
      
    dev.off()
  }
  
  
  
  #  compute change in specific discharge instead of relative change
  upag=match(stid ,UpArea$outlets)
  upareaL=UpArea$upa[upag]
  
  df <- data.frame(
    time = c(1951:2020),
    Dres = Restrend*1000/upareaL,
    Dsoc = Soctrend*1000/upareaL,
    Dclim=Climtrend*1000/upareaL
  )
  
  # df_long <- reshape2::melt(df, id.vars = "time", variable.name = "driver", value.name = "value")
  # 
  # ggplot(df_long, aes(x = time, y = value, fill = driver)) +
  #   geom_area() +
  #   theme_minimal() +
  #   labs(title = "10y RL ",
  #        x = "Time",
  #        y = "Value",
  #        fill = "driver")
  # 
  total=Restrend+Soctrend+Climtrend
  TotalChange=(RLGPDH-RLGPDS[1])
  
  df$catch=rep(stid,length(Restrend))
  
  dfall=rbind(dfall,df)

}

#boxplot of the different contributions... in 2020
df_2020=dfall[which(dfall$time==2020),-5]

df_l2020 <- reshape2::melt(df_2020, id.vars = "time", variable.name = "driver", value.name = "value")


meds <- c(by(df_l2020$value, df_l2020$driver, median))
q <- c(by(df_l2020$value, df_l2020$driver, quantile))
merdecol=match(df_l2020$driver,names(meds))
df_l2020$col=meds[merdecol]
palet=c(hcl.colors(9, palette = "Blues", alpha = NULL, rev = T, fixup = TRUE))
colorz = c("Dres" ='darkgrey',"Dsoc" ='orange',"Dclim" ='royalblue')

ggplot(df_l2020, aes(x=factor(driver), y=value,fill=driver)) +
  ggdist::stat_halfeye(adjust = 1.5, width = .6, justification = -.3, .width = .1,scale=0.6,
                       trim=TRUE, point_colour = NA, normalize="groups") + 
  #ggdist::stat_dots(side = "left", dotsize = .1, justification = 1.1, binwidth = .1) +
  #geom_boxplot(width = .1, outlier.shape = NA) +
  geom_boxplot(notch=F,width = .1,position=position_dodge(.9),alpha=.8,aes(fill=driver),linewidth=0.8,outlier.alpha = 0.4)+
  scale_y_continuous(limits = c(-200,200),name="Contribution to change in 10y Flood RL(l/s/km2) (1951-2020)",breaks = seq(-100,100,by=10),minor_breaks = seq(-100,100,5))+
  scale_x_discrete(labels=c("Dres" = "Reservoirs", "Dsoc" = "Land use + water use",
                            "Dclim" = "Climate"),name="Driver")+
  coord_cartesian(ylim = c(-20,20))+
  scale_fill_manual(values = colorz, name=" ") +
  theme(axis.title=element_text(size=16, face="bold"),
        axis.text = element_text(size=12),
        panel.background = element_rect(fill = "white", colour = "white"),
        panel.grid = element_blank(),
        panel.border = element_rect(linetype = "solid", fill = NA, colour="black"),
        legend.title = element_text(size=14),
        legend.text = element_text(size=12),
        legend.position = "none",
        panel.grid.major = element_line(colour = "grey60"),
        panel.grid.minor.y = element_line(colour = "grey80",linetype="dashed"),
        legend.key = element_rect(fill = "transparent", colour = "transparent"),
        legend.key.size = unit(.8, "cm"))





#Deep comparison between 1951 for all catchments

data_fy=df_2020=dfall[which(dfall$time==2020), ]
data_fyp=data_fy[which(is.infinite(data_fy$Dclim)),]

### Hybas07 ----
Catchmentrivers7=read.csv(paste0(hydroDir,"/Catchments/from_hybas_eu_onlyid.csv"),encoding = "UTF-8", header = T, stringsAsFactors = F)
hybas07 <- read_sf(dsn = paste0(hydroDir,"/Catchments/hydrosheds/hybas_eu_lev07_v1c.shp"))
hybasf7=fortify(hybas07) 
Catamere07=inner_join(hybasf7,Catchmentrivers7,by= "HYBAS_ID")
outletname="outletsv8_hybas07_01min"
outf=outletopen(hydroDir,outletname)
outf$latlong=paste(round(outf$Var1,4),round(outf$Var2,4),sep=" ")
Catamere07$llcoord=paste(round(Catamere07$POINT_X,4),round(Catamere07$POINT_Y,4),sep=" ") 
Catf7=inner_join(Catamere07,outf,by= c("llcoord"="latlong"))
cst7=st_transform(Catf7,  crs=3035)

#load reservoir influence
reservoirdog=outletopen(paste0(hydroDir,"/reservoirs"),"res_ratio_diff_2020-1951")
reservoirdog$llcoord=paste(round(reservoirdog$Var1,4),round(reservoirdog$Var2,4),sep=" ") 
reservoir_bg=inner_join(reservoirdog,outf,by=c("llcoord"="latlong"))


hybmatch=match(data_fy$catch,cst7$outlets)
data_fy$hybasID=cst7$HYBAS_ID[hybmatch]

resmatch=match(data_fy$catch,reservoir_bg$outlets.y)

data_fy$resinf=reservoir_bg$outlets.x[resmatch]
plot(data_fy$resinf,data_fy$Dres)
data_reso=data_fy[which(data_fy$resinf>0),]
driver="socioeconomic changes"
driver="climate"
driver="reservoirs"
#catchment plot with difference
data=right_join(catmap,data_fy,by = c("outlets"="catch"))
title=paste0("10y flood RL change due to ", driver ," (1951-2020)")
legend="change(l/s/km2)"
palet=c(hcl.colors(11, palette = "RdYlBu", alpha = NULL, rev = F, fixup = TRUE))
br=c(-100,-80,-50,-30,-20,-10,-8,-6,-4,-2,0,2,4,6,8,10,20,30,50,80,100)
br=c(-100,-80,-50,-20,-10,0,10,20,50,80,100)
limi=c(-50,50)
trans=scales::modulus_trans(.8)
colNA="darkgrey"


datapl <- st_as_sf(data_fy, coords = c("Var1", "Var2"), crs = 4326)
datapl <- st_transform(datapl, crs = 3035)

ggplot(basemap) +
  geom_sf(fill="white")+
  #geom_sf(data=data,aes(fill=Dclim,geometry=geometry),color="transparent")+ 
  geom_sf(data=datapl,aes(geometry=geometry,fill=Dclim),alpha=.7,stroke=0,shape=16)+
  geom_sf(fill=NA, color="grey") +
  scale_fill_gradientn(
    colors=palet,
    breaks=br,limits=limi,
    oob = scales::squish,na.value=colNA, name=legend,trans=trans)   +
  #geom_sf(data=dt3,aes(geometry=geometry,fill=factor(plot)),color="gray12")+ 
  #scale_fill_manual(values=c("1"="purple","-1"="tan2"), name="IRES trend")+
  coord_sf(xlim = c(min(nco[,1]),max(nco[,1])), ylim = c(min(nco[,2]),max(nco[,2])))+
  labs(x="Longitude", y = "Latitude")+
  guides(fill = guide_colourbar(barwidth = 20, barheight = .8))+
  theme(axis.title=element_text(size=tsize),
        panel.background = element_rect(fill = "aliceblue", colour = "grey1"),
        panel.border = element_rect(linetype = "solid", fill = NA, colour="black"),
        legend.title = element_text(size=tsize),
        legend.text = element_text(size=osize),
        legend.position = "bottom",
        legend.title.position = "top",
        panel.grid.major = element_line(colour = "grey70"),
        panel.grid.minor = element_line(colour = "grey90"),
        legend.key = element_rect(fill = "transparent", colour = "transparent"),
        legend.key.size = unit(.8, "cm"))+
  ggtitle(title)




### Changes in RLs or RPs at pixel and catchment levels ----
datarF=RLGPD_H
datarF$IRES=0
Plot.change=plotchangemaps_qsp(basemap,catmap=catmap,datarF, UpArea,GHR_riv,HydroRsf, law="GPD",type="RLchange",period=c(1951,2020),parlist,valuenames,haz=haz)
Plot.change[[3]]

#Deep comparison between 1951 for all catchments

data_fy=df_2020=dfall[which(dfall$time==2020), ]
data_fyp=data_fy[which(is.infinite(data_fy$Dclim)),]

### Hybas07 ----
Catchmentrivers7=read.csv(paste0(hydroDir,"/Catchments/from_hybas_eu_onlyid.csv"),encoding = "UTF-8", header = T, stringsAsFactors = F)
hybas07 <- read_sf(dsn = paste0(hydroDir,"/Catchments/hydrosheds/hybas_eu_lev07_v1c.shp"))
hybasf7=fortify(hybas07) 
Catamere07=inner_join(hybasf7,Catchmentrivers7,by= "HYBAS_ID")
outletname="outletsv8_hybas07_01min"
outf=outletopen(hydroDir,outletname)
outf$latlong=paste(round(outf$Var1,4),round(outf$Var2,4),sep=" ")
Catamere07$llcoord=paste(round(Catamere07$POINT_X,4),round(Catamere07$POINT_Y,4),sep=" ") 
Catf7=inner_join(Catamere07,outf,by= c("llcoord"="latlong"))
cst7=st_transform(Catf7,  crs=3035)

#load reservoir influence
reservoirdog=outletopen(paste0(hydroDir,"/reservoirs"),"res_ratio_diff_2020-1951")
reservoirdog$llcoord=paste(round(reservoirdog$Var1,4),round(reservoirdog$Var2,4),sep=" ") 
reservoir_bg=inner_join(reservoirdog,outf,by=c("llcoord"="latlong"))


hybmatch=match(data_fy$catch,cst7$outlets)
data_fy$hybasID=cst7$HYBAS_ID[hybmatch]

resmatch=match(data_fy$catch,reservoir_bg$outlets.y)

data_fy$resinf=reservoir_bg$outlets.x[resmatch]
plot(data_fy$resinf,data_fy$Dres)
data_reso=data_fy[which(data_fy$resinf>0),]
driver="socioeconomic changes"
driver="climate"
driver="reservoirs"
#catchment plot with difference
data=right_join(catmap,data_fy,by = c("outlets"="catch"))
title=paste0("10y flood RL change due to ", driver ," (1951-2020)")
legend="change(l/s/km2)"
palet=c(hcl.colors(11, palette = "RdYlBu", alpha = NULL, rev = F, fixup = TRUE))
br=c(-100,-80,-50,-30,-20,-10,-8,-6,-4,-2,0,2,4,6,8,10,20,30,50,80,100)
br=c(-100,-80,-50,-20,-10,0,10,20,50,80,100)
limi=c(-1,1)
trans=scales::modulus_trans(.8)
colNA="darkgrey"


ggplot(basemap) +
  geom_sf(fill="white")+
  geom_sf(data=data,aes(fill=Dclim,geometry=geometry),color="transparent")+ 
  geom_sf(fill=NA, color="grey") +
  scale_fill_gradientn(
    colors=palet,
    breaks=br,limits=limi,
    oob = scales::squish,na.value=colNA, name=legend,trans=trans)   +
  #geom_sf(data=dt3,aes(geometry=geometry,fill=factor(plot)),color="gray12")+ 
  #scale_fill_manual(values=c("1"="purple","-1"="tan2"), name="IRES trend")+
  coord_sf(xlim = c(min(nco[,1]),max(nco[,1])), ylim = c(min(nco[,2]),max(nco[,2])))+
  labs(x="Longitude", y = "Latitude")+
  guides(fill = guide_colourbar(barwidth = 20, barheight = .8))+
  theme(axis.title=element_text(size=tsize),
        panel.background = element_rect(fill = "aliceblue", colour = "grey1"),
        panel.border = element_rect(linetype = "solid", fill = NA, colour="black"),
        legend.title = element_text(size=tsize),
        legend.text = element_text(size=osize),
        legend.position = "bottom",
        legend.title.position = "top",
        panel.grid.major = element_line(colour = "grey70"),
        panel.grid.minor = element_line(colour = "grey90"),
        legend.key = element_rect(fill = "transparent", colour = "transparent"),
        legend.key.size = unit(.8, "cm"))+
  ggtitle(title)




ggsave("plots/RLChange_ResImpact.jpg", width=15, height=20, units=c("cm"),dpi=500)
#maybe I should compute this change is qsp


library(exactextractr)
rastforest=raster( "D:/tilloal/Documents/06_Floodrivers/landuse/fracforest_ch20201951.tif")

rastsealed=raster( "D:/tilloal/Documents/06_Floodrivers/landuse/fracsealed_ch20201951.tif")


data$forestchange<- exact_extract(rastforest, data, 'mean')

data$sealedchange <- exact_extract(rastsealed, data, 'mean')

plot(data$forestchange,data$Dsoc)

#compare this with the chance from the SocCF scenario
haz="flood"
Plot.change=plotchangemaps(basemap,catmap=catmap,RLGPD_SCF, law="GPD",type="RLchange",period=c(1951,2020),parlost,valuenames)
Plot.change[[3]]

#success

#check individual catchments with unexpected changes to see if thresold is problematic (E.g murcia)


#Look by hydroregions as well could be cool

library(ghibli)
library(ggdist)

GridHR=raster( paste0(hydroDir,"/HydroRegions_raster_WGS84.tif"))
GHR=as.data.frame(GridHR,xy=T)
GHR=GHR[which(!is.na(GHR[,3])),]
GHR$llcoord=paste(round(GHR$x,4),round(GHR$y,4),sep=" ") 
GHR_riv=inner_join(GHR,outf,by= c("llcoord"="latlong"))
GHshpp <- read_sf(dsn ="Z:/ClimateRun4/nahaUsers/tilloal/HydroRegions/her_all_adjusted.shp")
HydroRsf=fortify(GHshpp) 

data=right_join(GHR_riv,data,by = c("outlets"))

pointagg=aggregate(list(Schange=data$Dsoc,Cchange=data$Dclim,Rchange=data$Dres),
                   by = list(HydroR=data$HydroRegions_raster_WGS84),
                   FUN = function(x) c(med=median(x,na.rm=T),dev=sd(x,na.rm=T),len=length(x),q1=quantile(x, 0.25, na.rm=T),q3=quantile(x, 0.75, na.rm=T)))
pointagg <- do.call(data.frame, pointagg)

pointaglite=pointagg[,c(1,2,7,12)]
#pointaglite$Rchange.med[which(abs(round(pointaglite$Rchange.med,3))<0.01)]=NA

pa_l <- reshape2::melt(pointaglite, id.vars = "HydroR", variable.name = "driver", value.name = "value")

meds <- c(by(pa_l$value, pa_l$driver, median))
q <- c(by(pa_l$value, pa_l$driver, quantile))
merdecol=match(pa_l$driver,names(meds))
pa_l$col=meds[merdecol]
pa_l$order="1"
pa_l$order[which(pa_l$driver=="Schange.med")]="2"
pa_l$order[which(pa_l$driver=="Rchange.med")]="3"
palet=c(hcl.colors(9, palette = "Blues", alpha = NULL, rev = T, fixup = TRUE))
colorz = c("Rchange.med" ='darkgrey',"Schange.med" ='orange',"Cchange.med" ='royalblue')

# ggplot(pa_l, aes(x=factor(driver), y=value)) +
#   ggdist::stat_halfeye(adjust = 1, width = .3, .width = .1, justification = -.3, point_colour = NA) + 
#   ggdist::stat_dots(side = "left", dotsize = .1, justification = 1.1, binwidth = .1) +
#   #geom_boxplot(width = .1, outlier.shape = NA) +
#   geom_boxplot(notch=F,position=position_dodge(.9),alpha=.8,aes(fill=driver),linewidth=0.8,outlier.alpha = 0.4)+
#   scale_y_continuous(limits = c(-200,200),name="Contribution to change in 10y Flood RL(l/s/km2) (1951-2020)",breaks = seq(-100,100,by=10),minor_breaks = seq(-100,100,5))+
#   scale_x_discrete(labels=c("Rchange.med" = "Reservoirs", "Schange.med" = "Land use + water use",
#                             "Cchange.med" = "Climate"),name="Driver")+
#   coord_cartesian(ylim = c(-15,15))+
#   scale_fill_manual(values = colorz, name=" ") +
#   theme(axis.title=element_text(size=16, face="bold"),
#         axis.text = element_text(size=12),
#         panel.background = element_rect(fill = "white", colour = "white"),
#         panel.grid = element_blank(),
#         panel.border = element_rect(linetype = "solid", fill = NA, colour="black"),
#         legend.title = element_text(size=14),
#         legend.text = element_text(size=12),
#         legend.position = "none",
#         panel.grid.major = element_line(colour = "grey60"),
#         panel.grid.minor.y = element_line(colour = "grey80",linetype="dashed"),
#         legend.key = element_rect(fill = "transparent", colour = "transparent"),
#         legend.key.size = unit(.8, "cm"))


limi=c(-1,1)


ggplot(pa_l, aes(x=order, y=value, fill=driver)) +
  #stat_dotsinterval(subguide = 'integer')+
  # Line below sets the Studio Ghibli color pallete, for the sake of nostalgia )
  geom_boxplot(width = 0.25,position = position_nudge(x = -0.2)) +
  theme(axis.title=element_text(size=12, face="bold"),
        axis.text = element_text(size=12),
        panel.background = element_rect(fill = "white", colour = "white"),
        panel.grid = element_blank(),
        panel.border = element_rect(linetype = "solid", fill = NA, colour="black"),
        legend.title = element_text(size=14),
        legend.text = element_text(size=10),
        legend.position = "none",
        panel.grid.major = element_line(colour = "grey60"),
        panel.grid.minor.y = element_line(colour = "grey80",linetype="dashed"),
        legend.key = element_rect(fill = "transparent", colour = "transparent"),
        legend.key.size = unit(.8, "cm"))+
  scale_fill_manual(values = colorz, name=" ") +
  scale_y_continuous(limits = c(-200,200),name="Contribution to change in 10y drought RL(l/s/km2) (1951-2020)",
                     breaks = seq(-100,100,by=1),minor_breaks = seq(-100,100,.25))+
  scale_x_discrete(labels=c("3" = "Reservoirs + Water use", "2" = "Land use",
                            "1" = "Climate"),name="Driver")+
  coord_cartesian(ylim = limi)+
  # Line below adds dot plots from {ggdist} package 
 # stat_dots(side = "left", justification = 1.12, color=NA,binwidth = unit(c(0.1, Inf), "mm"), overflow = "compress", alpha = 0.75,stroke=0) +
  # Line below adds half-violin from {ggdist} package
  stat_halfeye(adjust = 1, width = .6, justification = -.2, .width = .1,scale=0.5,
               trim=TRUE, point_colour = NA, interval_color=NA, normalize="groups")

ggsave("D:/tilloal/Documents/LFRuns_utils/Figures/Bxplot_droughtchange2.jpg", width=20, height=15, units=c("cm"),dpi=300)

#natch pointagg with hybasf
pointplot=inner_join(HydroRsf,pointagg,by= c("Id"="HydroR"))


pointplot$Rchange.med
drivers=c("climate","reservoirs + water demand", "land use")
title=paste0("Contribution of ",drivers[3]," changes to \ntotal changes in 10-y RL drought")
palet=c(hcl.colors(11, palette = "RdYlBu", alpha = NULL, rev = F, fixup = TRUE))

datapl <- st_as_sf(data, coords = c("Var1.y", "Var2.y"), crs = 4326)
datapl <- st_transform(datapl, crs = 3035)
limi=c(-.5,.5)

pl2=ggplot(basemap) +
  geom_sf(fill="white")+
  geom_sf(data=pointplot,aes(fill=Schange.med,geometry=geometry),color="transparent")+ 
  geom_sf(fill=NA, color="grey") +
  geom_sf(data=datapl,aes(geometry=geometry,color=Dsoc),alpha=.6,stroke=0,shape=16, size=.8)+
  geom_sf(data=datapl,aes(geometry=geometry),col="grey20",alpha=1,stroke=0.05,shape=1,size=.8)+ 
  scale_fill_gradientn(
    colors=palet, limits=limi,
    oob = scales::squish,na.value=colNA, name=legend)   +
  scale_color_gradientn(
    colors=palet, limits=limi,
    oob = scales::squish,na.value=colNA, name=legend)   +
  #geom_sf(data=dt3,aes(geometry=geometry,fill=factor(plot)),color="gray12")+ 
  #scale_fill_manual(values=c("1"="purple","-1"="tan2"), name="IRES trend")+
  coord_sf(xlim = c(min(nco[,1]),max(nco[,1])), ylim = c(min(nco[,2]),max(nco[,2])))+
  labs(x="Longitude", y = "Latitude")+
  guides(fill = guide_colourbar(barwidth = 20, barheight = .8))+
  theme(axis.title=element_text(size=tsize),
        panel.background = element_rect(fill = "aliceblue", colour = "grey1"),
        panel.border = element_rect(linetype = "solid", fill = NA, colour="black"),
        legend.title = element_text(size=tsize),
        legend.text = element_text(size=osize),
        legend.position = "bottom",
        legend.title.position = "top",
        panel.grid.major = element_line(colour = "grey70"),
        panel.grid.minor = element_line(colour = "grey90"),
        legend.key = element_rect(fill = "transparent", colour = "transparent"),
        legend.key.size = unit(.8, "cm"))+
  ggtitle(title)

pl2

ggsave("D:/tilloal/Documents/LFRuns_utils/TrendAnalysis/plots/map_drchange_soc.jpg",pl2, width=15, height=20, units=c("cm"),dpi=1200)


pointplot$forestchange<- exact_extract(rastforest, pointplot, 'mean')

pointplot$sealedchange <- exact_extract(rastsealed, pointplot, 'mean')




pl3=ggplot(basemap) +
  geom_sf(fill="white")+
  geom_sf(data=pointplot,aes(fill=sealedchange,geometry=geometry),color="transparent")+ 
  geom_sf(fill=NA, color="grey") +
  scale_fill_gradientn(
    colors=palet, limits=c(-.05,.05),
    oob = scales::squish,na.value=colNA, name=legend)   +
  #geom_sf(data=dt3,aes(geometry=geometry,fill=factor(plot)),color="gray12")+ 
  #scale_fill_manual(values=c("1"="purple","-1"="tan2"), name="IRES trend")+
  coord_sf(xlim = c(min(nco[,1]),max(nco[,1])), ylim = c(min(nco[,2]),max(nco[,2])))+
  labs(x="Longitude", y = "Latitude")+
  guides(fill = guide_colourbar(barwidth = 20, barheight = .8))+
  theme(axis.title=element_text(size=tsize),
        panel.background = element_rect(fill = "aliceblue", colour = "grey1"),
        panel.border = element_rect(linetype = "solid", fill = NA, colour="black"),
        legend.title = element_text(size=tsize),
        legend.text = element_text(size=osize),
        legend.position = "bottom",
        legend.title.position = "top",
        panel.grid.major = element_line(colour = "grey70"),
        panel.grid.minor = element_line(colour = "grey90"),
        legend.key = element_rect(fill = "transparent", colour = "transparent"),
        legend.key.size = unit(.8, "cm"))

pl3



plot(pointplot$sealedchange,pointplot$Rchange.med)



ggplot(df_l2020, aes(x=factor(driver), y=value)) +
  ggdist::stat_halfeye(adjust = 1, width = .3, .width = .1, justification = -.3, point_colour = NA) + 
  ggdist::stat_dots(side = "left", dotsize = .1, justification = 1.1, binwidth = .1) +
  geom_boxplot(width = .1, outlier.shape = NA) +
  # geom_boxplot(notch=F,position=position_dodge(.9),alpha=.8,aes(fill=driver),linewidth=0.8,outlier.alpha = 0.4)+
  scale_y_continuous(limits = c(-200,200),name="Contribution to change in 10y Flood RL(l/s/km2) (1951-2020)",breaks = seq(-100,100,by=10),minor_breaks = seq(-100,100,5))+
  scale_x_discrete(labels=c("Dres" = "Reservoirs", "Dsoc" = "Land use + water use",
                            "Dclim" = "Climate"),name="Driver")+
  coord_cartesian(ylim = c(-15,15))+
  scale_fill_manual(values = colorz, name=" ") +
  theme(axis.title=element_text(size=16, face="bold"),
        axis.text = element_text(size=12),
        panel.background = element_rect(fill = "white", colour = "white"),
        panel.grid = element_blank(),
        panel.border = element_rect(linetype = "solid", fill = NA, colour="black"),
        legend.title = element_text(size=14),
        legend.text = element_text(size=12),
        legend.position = "none",
        panel.grid.major = element_line(colour = "grey60"),
        panel.grid.minor.y = element_line(colour = "grey80",linetype="dashed"),
        legend.key = element_rect(fill = "transparent", colour = "transparent"),
        legend.key.size = unit(.8, "cm"))

#boxplot with Hydroregions instead
library(ggdist)
ggplot(pointplot, aes(x=1,Rchange.med)) + 
  ggdist::stat_halfeye(adjust = .5, width = .3, .width = 1, justification = -.3, point_colour = NA) + 
  ggdist::stat_dots(side = "left", dotsize = .9, justification = 1.1, binwidth = .1)


#Function to compare two scenarios




#I need to work on this function
plotSceComp_Clevel=function(basemap,data1,data2,date=c(2020),hybasf,valuenames,nco){
  
  
  #data=right_join(outloc,datar,by = c("outlets"="unikout"))
  data1=RLGPD_H
  data2=RLGPD_RCF
  
  
  datacomb=inner_join(data1,data2,by="unikout")
  
  
  #mcor=unique(match(data$unikout,parlist$catchment))
  date=1951
  seq=c(1951:2020)
  datacol=names(datacomb)
  aystat=c()
  for (it in 1:length(seq)){
    yr=seq[it]
    dloc1=paste0("X",yr,".x")
    dloc2=paste0("X",yr,".y")
    #this need to be changed to make advantage of the whole data
    vc1=match(dloc1,datacol)
    vc2=match(dloc2,datacol)
    tmpval=(datacomb[,vc1])/(datacomb[,vc2])*100-100
    stats=c(yr,mean(tmpval),median(abs(tmpval)),median(tmpval),sd(tmpval))
    aystat=rbind(aystat,stats)
  }
  plot(aystat[,2])
  dateloc1=paste0("X",date,".x")
  dateloc2=paste0("X",date,".y")
  datacol=names(datacomb)
  
  vcol1=match(dateloc1,datacol)
  vcol2=match(dateloc2,datacol)
  tmpval=(datacomb[,vcol1])/(datacomb[,vcol2])*100-100
  datacomb$plot=tmpval
  
  #difference in changes
  ploc1x=paste0("X",period[1],".x")
  ploc1y=paste0("X",period[1],".y")
  
  ploc2x=paste0("X",period[2],".x")
  ploc2y=paste0("X",period[2],".y")
  
  datacol=names(datacomb)
  
  vp1x=match(ploc1x,datacol)
  vp2x=match(ploc2x,datacol)
  
  vp1y=match(ploc1y,datacol)
  vp2y=match(ploc2y,datacol)
  
  changef=(datacomb[,vp1x])/(datacomb[,vp2x])*100-100
  changecf=(datacomb[,vp1y])/(datacomb[,vp2y])*100-100
  changediff=(changecf-changef)
  hist(changediff,xlim=c(-10,10),breaks=10000)
  datacomb$chdiff=changediff
  
  title=paste0("2020 Relative difference in 100-years flood RL (factual-counterfactual)")
  palet=c(hcl.colors(11, palette = "RdBu", alpha = NULL, rev = T, fixup = TRUE))
  #aggregate by hybas07 catchment
  
  # pointagg=aggregate(list(Rchange=datacomb[,156]),
  #                    by = list(HYBAS_ID=datacomb$HYBAS_ID.x),
  #                    FUN = function(x) c(mean=mean(x)))
  # pointagg <- do.call(data.frame, pointagg)
  
  trans=scales::modulus_trans(.6)
  colNA="darkgrey"
  
  
  
  
  #give a projection to pixels
  #match with hybas
  
  points=right_join(GHR_riv,datacomb,by = c("outlets"="unikout"))
  points=right_join(catmap,points,by = c("outlets"="outlets"))
  points <- st_as_sf(datacomb[,c(1,2,3,5,157)], coords = c("Var1.x", "Var2.x"), crs = 4326)
  points <- st_transform(points, crs = 3035)
  
  
  #natch pointagg with hybasf
  # pointsag=inner_join(hybasf,pointagg,by= "HYBAS_ID")
  
  
  # #TBC with the river pixels and better colors
  #   ggplot(basemap) +
  #   geom_sf(data=pointsInside,aes(geometry=geometry, fill=sign),alpha=.7,size=2,stroke=0,shape=21,color="black")+
  #   coord_sf(xlim = c(min(nco[,1]),max(nco[,1])), ylim = c(min(nco[,2]),max(nco[,2])))+
  #   scale_fill_manual(values=c("tomato4","steelblue4"), name="Significant trends")
  #   
  # ggplot(regiod)+
  # geom_sf(mapping=aes(geometry=geometry,group=name,fill=name))
  legend="Difference(%)"
  br=c(-100,-80,-50,-30,-20,-10,0,10,20,30,50,80,100)
  limi=c(-30,30)
  tsize=16
  osize=12
  
  ocrap<-ggplot(basemap) +
    geom_sf(fill="gray85",color="darkgrey",size=0.5)+
    geom_sf(data=points,aes(fill=plot,geometry=geometry),alpha=1,color="transparent")+ 
    geom_sf(fill=NA, color="grey") +
    #geom_sf(regiod,mapping=aes(geometry=geometry,group=name),fill=NA,color="orange")+
    coord_sf(xlim = c(min(nco[,1]),max(nco[,1])), ylim = c(min(nco[,2]),max(nco[,2])))+
    ggtitle(title)+
    scale_fill_gradientn(
      colors=palet,
      breaks=br,limits=limi,trans=trans,
      oob = scales::squish,na.value=colNA, name=legend)   +
    labs(x="Longitude", y = "Latitude")+
    guides(colour = guide_colourbar(barwidth = 1.5, barheight = 10),
           fill = guide_legend(override.aes = list(size = 10)))+
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
  ocrap
  
  
  #investigate catchments with strange values (basically positive values)
  
  
  
  title=paste0("Contribution of Socioeconomic changes to total changes in 100-y RL flood")
  legend="Contribution (%)"
  palet=c(hcl.colors(11, palette = "RdBu", alpha = NULL, rev = F, fixup = TRUE))
  br=c(-100,-80,-50,-30,-20,-10,-8,-6,-4,-2,0,2,4,6,8,10,20,30,50,80,100)
  limi=c(-10,10)
  ocrap2<-ggplot(basemap) +
    geom_sf(fill="gray85",color="darkgrey",size=0.5)+
    geom_sf(data=points,aes(fill=chdiff,geometry=geometry),color="transparent")+ 
    geom_sf(fill=NA, color="grey") +
    #geom_sf(regiod,mapping=aes(geometry=geometry,group=name),fill=NA,color="orange")+
    coord_sf(xlim = c(min(nco[,1]),max(nco[,1])), ylim = c(min(nco[,2]),max(nco[,2])))+
    ggtitle(title)+
    scale_fill_gradientn(
      colors=palet,
      breaks=br,limits=limi,
      oob = scales::squish,na.value=colNA, name=legend)   +
    labs(x="Longitude", y = "Latitude")+
    guides(colour = guide_colourbar(barwidth = 1.5, barheight = 10),
           fill = guide_legend(override.aes = list(size = 10)))+
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
  
  ocrap2
  #ggsave("D:/tilloal/Documents/LFRuns_utils/Figures/RLChange_SocImpact.jpg", width=20, height=15, units=c("cm"),dpi=1500)
  
  datacomb$impact=changef*changecf
  datacomb$impact[which(datacomb$impact>0)]="no change"
  datacomb$impact[which(datacomb$impact<0)]="change"
  
  length(which(datacomb$impact=="change"))/length(datacomb$impact)
  points <- st_as_sf(datacomb[,c(1,2,3,5,158)], coords = c("Var1.x", "Var2.x"), crs = 4326)
  points <- st_transform(points, crs = 3035)
  title=paste0("Socioeconomic impact on 100y RL flood changes sign")
  legend="Impact"
  palet=c(hcl.colors(11, palette = "RdBu", alpha = NULL, rev = F, fixup = TRUE))
  br=c(-100,-80,-50,-30,-20,-10,-8,-6,-4,-2,0,2,4,6,8,10,20,30,50,80,100)
  limi=c(-10,10)
  
  ocrap3<-ggplot(basemap) +
    geom_sf(fill="gray85",color="darkgrey",size=0.5)+
    geom_sf(data=points,aes(col=factor(impact),geometry=geometry),alpha=1,size=0.25,stroke=0,shape=15)+ 
    geom_sf(fill=NA, color="grey") +
    #geom_sf(regiod,mapping=aes(geometry=geometry,group=name),fill=NA,color="orange")+
    coord_sf(xlim = c(min(nco[,1]),max(nco[,1])), ylim = c(min(nco[,2]),max(nco[,2])))+
    ggtitle(title)+
    scale_color_manual(name=legend,
                       values=c("change" ="darkred","no change"="white"),
                       na.value="grey")+
    # scale_color_gradientn(
    #   colors=palet,
    #   breaks=br,limits=limi,trans=trans,
    #   oob = scales::squish,na.value=colNA, name=legend)   +
    labs(x="Longitude", y = "Latitude")+
    guides(color = guide_legend(override.aes = list(size = 10)))+
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
  
  ocrap3
  
}






