


##########################################################################################
############   THIS SCRPT IS FOR RUNNING 1 SQUARE OF THE HERA DOMAIN ON A HPC ############
##########################################################################################

setwd(dirname(rstudioapi::getSourceEditorContext()$path))
setwd("D:/tilloal/Documents/LFRuns_utils/TrendAnalysis/R")
getwd()
source("functions_trends.R")

#Set data directory
hydroDir<-("D:/tilloal/Documents/LFRuns_utils/ChangingHydroExtremes/data")

# Arguments importation #############################

#default: tail is high
tail="low"
haz = "drought"
var = "dis"
outlets="RNetwork"
outletname <- "/GeoData/efas_rnet_100km_01min"
season="nonfrost"
Nsq = 42
sce <- "Histo"



rspace= read.csv(paste0(hydroDir,"/subspace_efas.csv"))
rspace=rspace[,-1]
nrspace=rspace[Nsq,]
print(nrspace)
nameout="UCRnet"
outhybas=outletopen(hydroDir,outletname,nrspace)
Idstart=as.numeric(Nsq)*100000
if (length(outhybas$outlets)>0){
  outhybas$outlets=seq((Idstart+1),(Idstart+length(outhybas$outlets)))
}

unikout=outhybas$outlets
outhybas$latlong=paste(round(outhybas$Var1,4),round(outhybas$Var2,4),sep=" ")

UpAname="/GeoData/upArea_European_01min.nc"
#dir=valid_path

UpArea=UpAopen(hydroDir,UpAname,outhybas)
head(UpArea)
outhybas$upa=UpArea$upa



tsEvaTransformSeriesToStationaryPeakTrendX<-function (timeStamps, series, timeWindow, xid) 
{
  message("\ncomputing the trend on extremes...\n")
  series_above_threshold <- series
  series_above_threshold[-xid] <- NA
  detrend_result <- tsEvaDetrendTimeSeries(timeStamps, series_above_threshold, 
                                           timeWindow)
  detrended_series <- series - detrend_result@trendSeries
  detrended_series_above_threshold <- series_above_threshold - 
    detrend_result@trendSeries
  n_run_mean <- detrend_result@nRunMn
  #plot(detrend_result@trendSeries)
  variance_series <- tsEvaNanRunningVariance(detrended_series_above_threshold, 
                                             n_run_mean)
  variance_series <- tsEvaNanRunningMean(variance_series, ceiling(n_run_mean/2))
  std_dev_series <- sqrt(variance_series)
  avg_std_dev <- mean(std_dev_series, na.rm = TRUE)
  S <- 2
  N <- timeWindow * 4
  std_dev_error <- avg_std_dev * (2 * S^2/N^3)^(1/4)
  stationary_series <- detrended_series/std_dev_series
  stationary_series_no_na <- stats::na.omit(stationary_series)
  running_moments <- tsEvaNanRunningStatistics(stationary_series_no_na, 
                                               n_run_mean)
  stat_ser_3_mom <- tsEvaNanRunningMean(running_moments$rn3mom, 
                                        ceiling(n_run_mean))
  stat_ser_4_mom <- tsEvaNanRunningMean(running_moments$rn4mom, 
                                        ceiling(n_run_mean))
  trend_error <- avg_std_dev/sqrt(N)
  transformed_data <- list(runningStatsMulteplicity = n_run_mean, 
                           stationarySeries = stationary_series, trendSeries = detrend_result@trendSeries, 
                           trendSeriesNonSeasonal = NULL, trendError = trend_error, 
                           stdDevSeries = std_dev_series, stdDevSeriesNonSeasonal = NULL, 
                           stdDevError = std_dev_error * rep(1, length(std_dev_series)), 
                           timeStamps = timeStamps, nonStatSeries = series, statSer3Mom = stat_ser_3_mom, 
                           statSer4Mom = stat_ser_4_mom)
  return(transformed_data)
}
# file_out=file=paste0(workDir,"TSEVA/out/",sce,"/",haz,"/",season,"/ResCat6h_",outlets,haz,"_",Nsq,Quarter,"_1951_2020.Rdata")
# fex=file.exists(file_out)
# if (fex==TRUE){
#   f.size=as.numeric(fs::file_size(file_out))
# }else{
#   f.size=0
# }

#Scenario differentiation
if (sce=="Histo") code="h"
if (sce=="SCF") code="scf"
if (sce=="WStat") code="wcf"
if (sce=="RWStat") code="rwcf"

#Loop on all pixels within squarq
#Load the file
#loading the files as netcdf (needs to be checked offline)
if (code=="h"){
  filename=paste0("dis_",Nsq,"_1951_2020_",code,"_RNetwork")
}
if (code=="scf"){
  filename=paste0("dis_",Nsq,"_1951_2020_",code)
}
if (code=="wcf"){
  filename=paste0("dis_",Nsq,"_1951_2020_",code)
}
if (code=="rwcf"){
  filename=paste0("dis_",Nsq,"_1951_2020_",code)
}
dists=disNcopenloc(filename,hydroDir,outhybas,1)
df.dis=dists 
print(paste0("hazard: ",haz," opening square ", Nsq, " /88"))
timeStamps=(as.Date(df.dis$time,origin="1979-01-01"))
timeStamps=as.POSIXct(timeStamps-1/24)
txx=timeStamps
df.dis$timeStamps=txx

names(df.dis)[c(1,2)]=c("dis","outlets")

#loading the frost file for drought
if (haz=="drought"){
  load(file=paste0(hydroDir,"/Drought/catchment_frost.Rdata"))
  rmv=which(year(frostcat$time)==1950)
  frostcat=frostcat[-rmv,]
  #remove first day
  frostcat=frostcat[-1,]
  Catchmentrivers7=read.csv(paste0(hydroDir,"/GeoData/HYBAS07/from_hybas_eu_onlyid.csv"),encoding = "UTF-8", header = T, stringsAsFactors = F)
  outletname="GeoData/HYBAS07/outletsv8_hybas07_01min"
  outhyb07=outletopen(hydroDir,outletname,nrspace)
  catmatch=match(outhyb07$outlets,Catchmentrivers7$pointid)
  mycat=Catchmentrivers7[catmatch,]
  
  hybas07 <- read_sf(dsn = paste0(hydroDir,"/GeoData/HYBAS07/hybas_eu_lev07_v1c.shp"))
  hybasf7=fortify(hybas07) 
  Catamere07=inner_join(hybasf7,Catchmentrivers7,by= "HYBAS_ID")
  Catamere07$llcoord=paste(round(Catamere07$POINT_X,4),round(Catamere07$POINT_Y,4),sep=" ") 
  Catf7=inner_join(Catamere07,outhybas,by= c("llcoord"="latlong"))
  st_geometry(Catf7)=NULL	
  tail="low"
  
}

ThDir<-paste0(hydroDir,"/Thresholds")

TH=read.csv(paste0(ThDir,"/trenTH_SCF_",tail,"_",Nsq,".csv"))

#retain threshold from SCF run
thresh_vec=data.frame(TH$cid,TH$Th_new)
names(thresh_vec)=c("cid","th")

thresh_vec$cid=as.numeric(thresh_vec$cid)
Nsq=as.numeric(Nsq)
thresh_vec$cid=thresh_vec$cid-Nsq*10000
thresh_vec$cid=thresh_vec$cid+Nsq*100000

#load directly extremes
xtrempoints<-read.csv(paste0(ThDir,"/xtrempoints_SCF_",tail,"_",Nsq,".csv"))
xtrempoints$cid=xtrempoints$cid-Nsq*10000
xtrempoints$cid=xtrempoints$cid+Nsq*100000

startid=1
endid=length(unikout)
endid=20

TrendSave=c()
SdSave=c()
catlist=c()
epy=-1
for (idfix in startid:endid){
  #idfix=2
  start_time <- Sys.time()
  print(paste0("hazard:",haz," square: ", Nsq, " pixel: ",idfix,"/",endid))
  catch=as.numeric(unikout[idfix])
  
  
  timeStamps=txx
  xid=xtrempoints$id[which(xtrempoints$cid==catch)]
  thresh=thresh_vec[which(thresh_vec$cid==catch),]
  thresh=thresh$th
  frosttime=NA
  df.disX=disNcopenloc(filename,hydroDir,outhybas,idfix)
  timeAndSeries=data.frame(txx,df.disX$outlets)
  names(timeAndSeries)=c("date","Qs")
  rmv=which(as.integer(format(timeAndSeries$date, "%Y"))==1950)
  if (length(rmv)>0){
    timeAndSeries=timeAndSeries[-rmv,]
  }
  #remove first day to avoid errors
  timeAndSeries=timeAndSeries[-1,]
  if (haz=="drought"){
    trans="rev"
    #seasonal split
    catmat=Catf7[which(Catf7$outlets==catch),]
    Tcatmat=mycat[which(mycat$HYBAS_ID==catmat$HYBAS_ID),]
    Tcatchment=which(colnames(frostcat)==Tcatmat$pointid)
    
    intermit=interid(timeAndSeries,trans,WindowSize=7)
    interflag=intermit$flags[2]
    timeAndSeries=data.frame(timeAndSeries$date,intermit$trdis$Q7)
    #remove frost timesteps, this can be modified to do the anlysis only on frost moments
    if (length(Tcatchment)>0){
      frostserie=data.frame(frostcat[,1],frostcat[,Tcatchment])
      frosttime=which(frostserie[,2]<0)
    }else{
      frosttime=NA
    }
    minPeakDistanceInDays=30
    tail="low"
  }else if (haz=="flood"){
    
    minPeakDistanceInDays=7
    interflag=0
    timeAndSeries <- max_daily_value(timeAndSeries)
    tail="high"
    trans="ori"
  
  }
  
  names(timeAndSeries)=c("timestamp","dis")
  dt1=min(diff(timeAndSeries$timestamp),na.rm=T)
  dt=as.numeric(dt1)
  tdim=attributes(dt1)$units
  if (tdim=="hours") dt=dt/24
  if (dt==1){
    timeDays=timeAndSeries$timestamp
  }else{
    timeDays=unique(as.Date(timeAndSeries$timestamp))
  }
  
  bounds=c(year(timeDays[1]),year(timeDays[length(timeDays)]))
  realbound=bounds
  tbound=c(as.Date(paste0(realbound[1],"-12-31")),as.Date(paste0(realbound[2],"-12-31")))
  Impdates=seq(tbound[1],tbound[2],by="1 year")
  
  nv=length(unique(timeAndSeries$dis))
  if(length(na.omit(timeAndSeries$dis))>1 & interflag<3 & nv>15){
    
    if (length(which(is.na(timeAndSeries$dis)))>0){
      print("Na alert")
      seriefill=tsEvaFillSeries(timeAndSeries$timestamp,timeAndSeries$dis)
      timeAndSeries$dis=seriefill
    }
    names(timeAndSeries)=c("timestamp","data")
    
    if (haz=="drought" & length(!is.na(frosttime))>1){
      if (season=="nonfrost"){
        print("nonfrost season")
        timeAndSeries$data[frosttime]=NA
      }else if (season=="frost"){
        print("frost season")
        timeAndSeries$data[-frosttime]=NA
      }else if (season=="year"){
        print("no seasonal divide")
      }else {print("season must be frost or nonfrost")}
    }else{
      print("no frost season for this river")
    }
    
    timeWindow = 365.25*30; #time windows in days, the correction is done within the functions
    windowSize=366
    
    timeStamps=timeAndSeries$timestamp
    cat(paste0("\nsquare: ", Nsq, " pixel: ",idfix,"/",endid))
    
    timeStamps=as.POSIXct(timeAndSeries[,1])
    dt1=min(diff(timeStamps),na.rm=T)
    dt=as.numeric(dt1)
    tdim=attributes(dt1)$units
    if (tdim=="hours") dt=dt/24
    if (tdim=="seconds") dt=dt/3600
    series=timeAndSeries[,2]
    
    if (epy==-1){
      if (tail=="high") epy=3
      if (tail=="low") epy=2
    }
    
    #If the biggest event is more than 100 times greater than the second biggest event
    sanitycheck = computeAnnualMaxima(timeAndSeries);
    anmax=sanitycheck$annualMax[order(sanitycheck$annualMax,decreasing = T)]
    aloc=sanitycheck$annualMaxIndx[order(sanitycheck$annualMax,decreasing = T)][1]
    yrs=year(sanitycheck$annualMaxDate[order(sanitycheck$annualMax,decreasing = T)][1])
    x = anmax[1]/anmax[2]
    if (x > 100)
    {
      message(paste0("biggest event ", x, " times bigger than second biggest"))
      message ("removing this event from timeserie, reruning first steps")
      series[(aloc[1]-50):(aloc[1]+50)]=mean(series)
    }
    
    if (minPeakDistanceInDays == -1) stop('label parameter minPeakDistanceInDays must be set')
    
    if (tail=="low"){
      #default 7-day flow for low flow, can be modified by user
      lowdt=7
      start_index=1
      indices_to_extract <- seq(from = start_index, to = length(series), by = lowdt/dt)
      series=series[indices_to_extract]
      timeStamps=timeStamps[indices_to_extract]
      shape_bnd=c(-2,0)
      if (trans=="rev"){
        series=-1*series
      }else if(trans=="inv"){
        series=1/series
      }else if (trans=="lninv"){
        series=-log(series)
      }
    }
    #plot(series)
    if (sce =="SCF"){
      TrendTh=thresh
      message(paste0('\nevaluating long term variations of the peaks'))
      if (is.na(TrendTh)){
        TrendTh=try(tsEvaFindTrendThreshold(series, timeStamps, timeWindow),T)
        if (inherits(TrendTh, "try-error") || length(TrendTh) == 0 || all(is.na(TrendTh))) {
          TrendTh <- 0.1
        }
        print(TrendTh)
        trasfData = tsEvaTransformSeriesToStationaryPeakTrend( timeStamps, series, timeWindow, TrendTh);
      }else{
        if (TrendTh=="MMX"){
          trasfData = tsEvaTransformSeriesToStationaryMMXTrend( timeStamps, series, timeWindow);
          message("using MMX trend")
        }else{
          trasfData = tsEvaTransformSeriesToStationaryPeakTrend( timeStamps, series, timeWindow, TrendTh);
          c=0
          while ((length(unique(trasfData$trendSeries))<2)) {
            trasfData = tsEvaTransformSeriesToStationaryPeakTrend( timeStamps, series, timeWindow, TrendTh=TrendTh-c);
            c=c+0.1
          }
        }
      }
    
    }else{
      print("taking extremes from CF1")
      trasfData = tsEvaTransformSeriesToStationaryPeakTrendX( timeStamps, series, timeWindow, xid);
    }
   
    trendS=trasfData$trendSeries
    if (tail=="low"){
      trendS=-trasfData$trendSeries
    }
    sdS=trasfData$stdDevSeries
    
    #extract trend and variance once per year
    
    
  #Saving main outputs
  catlist=c(catlist,catch)
  TrendSave=cbind(TrendSave,trendS)
  SdSave=cbind(SdSave,sdS)
  end_time <- Sys.time()
  cat(paste0("\nloop duration: ",round(end_time-start_time,2)," seconds\n"))
  
  
  }else{print("no analysis")}
  
}
TrendSave=data.frame(TrendSave)
SdSave=data.frame(SdSave)
colnames(TrendSave)=colnames(SdSave)=catlist
TrendSave=cbind(timeStamps,TrendSave)
SdSave=cbind(timeStamps,SdSave)

Results=list(Trend=TrendSave,Variability=SdSave)
save(Results, file=paste0(hydroDir,"TrendVar/TrendVarX_",outlets,haz,"_",Nsq,"_1951_2020.Rdata"))