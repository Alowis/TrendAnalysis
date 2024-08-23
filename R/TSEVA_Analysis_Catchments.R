##########################################################################################
###################   DEMO of TSEVA for different locations   ############################
##########################################################################################
#Import functions from TSEVA
source("D:/tilloal/Documents/LFRuns_utils/TSEVA_demo/demo_functions.R")
library(xts)
library(RtsEva)
library(ggplot2)


#Experimental TSEVA
TsEvaNs<- function(timeAndSeries, timeWindow, transfType='trendPeaks',minPeakDistanceInDays=10,
                   seasonalityVar=NA,minEventsPerYear=-1, gevMaxima='annual',
                   ciPercentile=90, gevType = 'GEV', evdType = c('GEV', 'GPD'),
                   tail="high", epy=-1, lowdt=7, trans=NULL, TrendTh=NA){
  
  
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
  
  if ( transfType != 'trend' & transfType != 'seasonal' & transfType != 'trendCIPercentile'
       & transfType != 'seasonalCIPercentile' & transfType != 'trendPeaks'){
    stop('\nnonStationaryEvaJRCApproach: transfType can be in (trend, seasonal, trendCIPercentile, trendPeaks)')}
  
  if (minPeakDistanceInDays == -1) stop('label parameter minPeakDistanceInDays must be set')
  
  nonStationaryEvaParams = c()
  stationaryTransformData = c()
  
  # default shape parameter bounds
  shape_bnd=c(-0.5,1)
  
  if (tail=="low"){
    #default 7-day flow for low flow, can be modified by user
    start_index=1
    indices_to_extract <- seq(from = start_index, to = length(series), by = lowdt/dt)
    series=series[indices_to_extract]
    timeStamps=timeStamps[indices_to_extract]
    shape_bnd=c(-1,0)
    if (trans=="rev"){
      series=-1*series
    }else if(trans=="inv"){
      series=1/series
    }else if (trans=="lninv"){
      series=-log(series)
    }
    
  }
  if (transfType == 'trend'){
    message('\nevaluating long term variations of extremes')
    trasfData = tsEvaTransformSeriesToStationaryTrendOnly(timeStamps, series, timeWindow)
    gevMaxima = 'annual';
    potEventsPerYear = epy;
    minEventsPerYear = 1
    
    
  }else  if (transfType == 'trendChange'){
    message('\nevaluating long term variations of extremes and change point detection')
    trasfData = tsEvaTransformSeriesToStationaryTrendAndChangepts(timeStamps, series, timeWindow)
    gevMaxima = 'annual';
    potEventsPerYear = epy;
    minEventsPerYear = 1
    
  }else if (transfType == 'seasonal'){
    message('\nevaluating long term an seasonal variations of extremes')
    trasfData = tsEvaTransformSeriesToStationaryMultiplicativeSeasonality(timeStamps, series, timeWindow, seasonalityVar=seasonalityVar)
    gevMaxima = 'monthly'
    potEventsPerYear = 12
    minEventsPerYear = 12
    
  } else if (transfType == 'trendCIPercentile') {
    if (is.na(ciPercentile)){
      stop('For trendCIPercentile transformation the label parameter cipercentile is mandatory')
    }
    message(paste0('\nevaluating long term variations of extremes using the ', ciPercentile, 'th percentile'))
    trasfData = tsEvaTransformSeriesToStationaryTrendOnly_ciPercentile( timeStamps, series, timeWindow, ciPercentile);
    gevMaxima = 'annual'
    potEventsPerYear = epy
    minEventsPerYear = 1
    
  }else if (transfType == 'trendPeaks') {
    message(paste0('\nevaluating long term variations of the peaks'))
    if (is.na(TrendTh)){
      TrendTh=try(tsEvaFindTrendThreshold(series, timeStamps, timeWindow),T)
      if(length(TrendTh)==0){
        TrendTh=NA
      }
    }
    trasfData = tsEvaTransformSeriesToStationaryPeakTrend( timeStamps, series, timeWindow, TrendTh);
    gevMaxima = 'annual'
    potEventsPerYear = epy
    minEventsPerYear = 1
    
  } else  if (transfType == 'trendChangeCIPercentile'){
    if (is.na(ciPercentile)){
      stop('For trendCIPercentile transformation the label parameter cipercentile is mandatory')
    }
    message('\n evaluating long term variations of extremes using the ', ciPercentile, 'th percentile and change point detection')
    trasfData = tsEvaTransformSeriesToStationaryTrendAndChangepts_ciPercentile(timeStamps, series, timeWindow,ciPercentile)
    gevMaxima = 'annual';
    potEventsPerYear = epy;
    minEventsPerYear = 0
    
  } else if (transfType == 'seasonalCIPercentile') {
    if (is.na(ciPercentile)) stop('For seasonalCIPercentile transformation the label parameter cipercentile is mandatory')
    message(paste0('\nevaluating long term variations of extremes using the ', ciPercentile, 'th percentile\n'))
    trasfData = tsEvaTransformSeriesToStatSeasonal_ciPercentile( timeStamps, series, timeWindow, ciPercentile)
    gevMaxima = 'monthly'
    potEventsPerYear = 12
    minEventsPerYear = 6
  }
  
  
  dtn=min(diff(trasfData$timeStamps),na.rm=T)
  dtn=as.numeric(dtn)
  tdim=attributes(dtn)$units
  if (dtn<1) {
    pace=1/dtn
    tsDaily=seq(1,length(trasfData$timeStamps),by=pace)
    trasfData$stdDevSeriesOr=trasfData$stdDevSeries
    trasfData$trendSeriesOr=trasfData$trendSeries
    trasfData$stdDevErrorOr=trasfData$stdDevError
    trasfData$stdDevSeries=trasfData$stdDevSeries[tsDaily]
    trasfData$trendSeries=trasfData$trendSeries[tsDaily]
    trasfData$stdDevError=trasfData$stdDevError[tsDaily]
  }
  
  ms = data.frame(trasfData$timeStamps, trasfData$stationarySeries)
  minPeakDistance = minPeakDistanceInDays/dtn;
  
  #estimating the non stationary EVA parameters
  message('\nExecuting stationary eva')
  pointData = tsEvaSampleData(ms, potEventsPerYear, minEventsPerYear, minPeakDistanceInDays,tail);
  evaAlphaCI = .68; # in a gaussian approximation alphaCI~68% corresponds to 1 sigma confidence
  eva = tsEVstatistics(pointData, evaAlphaCI, gevMaxima, gevType, evdType,shape_bnd);
  
  if (eva$isValid==FALSE) {
    message("problem in the computation of EVA statistics")
  }
  
  eva[[2]]$GPDstat$thresholdError <- pointData$POT$thresholdError
  
  # !!! Assuming a Gaussian approximation to compute the standard errors for
  # the GEV parameters
  if (eva[[2]][[1]]$method[1]!="No fit") {
    epsilonGevX <- eva[[2]][[1]]$parameters[3]
    errEpsilonX <- epsilonGevX - eva[[2]][[1]]$paramCIs[1,3]
    sigmaGevX <- eva[[2]][[1]]$parameters[2]
    errSigmaGevX <- sigmaGevX - eva[[2]][[1]]$paramCIs[1, 2]
    muGevX <- eva[[2]][[1]]$parameters[1]
    errMuGevX <- muGevX - eva[[2]][[1]]$paramCIs[1, 1]
    
    message('\nTransforming to non stationary eva ...\n')
    epsilonGevNS = epsilonGevX;
    errEpsilonGevNS = errEpsilonX;
    sigmaGevNS = trasfData$stdDevSeries*sigmaGevX;
    
    #propagating the errors on stdDevSeries and sigmaGevX to sigmaGevNs.
    # err(sigmaNs) = sqrt{ [sigmaX*err(stdDev)]^2 + [stdDev*err(sigmaX)]^2 }
    # the error on sigmaGevNs is time dependant.
    errSigmaGevFit = trasfData$stdDevSeries*errSigmaGevX;
    errSigmaGevTransf = sigmaGevX*trasfData$stdDevError;
    errSigmaGevNS = (  errSigmaGevTransf^2   +  errSigmaGevFit^2  )^.5;
    muGevNS = trasfData$stdDevSeries*muGevX + trasfData$trendSeries;
    
    # propagating the errors on stdDevSeries, trendSeries and sigmaGevX to muGevNS.
    # err(muNs) = sqrt{ [muX*err(stdDev)]^2 + [stdDev*err(muX)]^2 + err(trend)^2 }
    # the error on muGevNS is time dependant.
    errMuGevFit = trasfData$stdDevSeries*errMuGevX;
    errMuGevTransf = (  (muGevX*trasfData$stdDevError)^2 + trasfData$trendError^2  )^.5;
    errMuGevNS = (  errMuGevTransf^2   +  errMuGevFit^2  )^.5;
    gevParams=c()
    gevParams$epsilon = epsilonGevNS;
    gevParams$sigma = sigmaGevNS;
    gevParams$mu = muGevNS;
    gevParams$annualMax=trasfData$nonStatSeries[pointData$annualMaxIndx]
    gevParams$monthlyMax=trasfData$nonStatSeries[pointData$monthlyMaxIndx]
    gevParams$annualMaxIndx=pointData$annualMaxIndx
    gevParams$monthlyMaxIndx=pointData$monthlyMaxIndx
    
    
    if(tolower(gevMaxima) == "annual") {
      gevParams$timeDelta <- 365.25
      gevParams$timeDeltaYears <- 1
    } else if(tolower(gevMaxima) == "monthly") {
      gevParams$timeDelta <- 365.25/12
      gevParams$timeDeltaYears <- 1/12
    }
    gevParamStdErr=c()
    gevParamStdErr$epsilonErr <- errEpsilonGevNS
    
    gevParamStdErr$sigmaErrFit <- errSigmaGevFit
    gevParamStdErr$sigmaErrTransf <- errSigmaGevTransf
    gevParamStdErr$sigmaErr <- errSigmaGevNS
    
    gevParamStdErr$muErrFit <- errMuGevFit
    gevParamStdErr$muErrTransf <- errMuGevTransf
    gevParamStdErr$muErr <- errMuGevNS
    
    gevObj=list()
    gevObj$method <- eva[[2]][[1]]$method
    gevObj$parameters <- gevParams
    gevObj$paramErr <- gevParamStdErr
    gevObj$stationaryParams <- eva[[2]][[1]]
    gevObj$objs$monthlyMaxIndexes <- pointData$monthlyMaxIndexes
  }else{
    
    epsilonGevX <- eva[[2]][[1]]$parameters[3]
    sigmaGevX <- eva[[2]][[1]]$parameters[2]
    muGevX <- eva[[2]][[1]]$parameters[1]
    message('\nTransforming to non stationary eva ...\n')
    epsilonGevNS = epsilonGevX;
    sigmaGevNS = trasfData$stdDevSeries*sigmaGevX;
    muGevNS = trasfData$stdDevSeries*muGevX + trasfData$trendSeries;
    
    gevParams=c()
    gevParams$epsilon = epsilonGevNS;
    gevParams$sigma = sigmaGevNS;
    gevParams$mu = muGevNS;
    gevParams$annualMax=trasfData$nonStatSeries[pointData$annualMaxIndx]
    gevParams$monthlyMax=trasfData$nonStatSeries[pointData$monthlyMaxIndx]
    gevParams$annualMaxIndx=pointData$annualMaxIndx
    gevParams$monthlyMaxIndx=pointData$monthlyMaxIndx
    
    if(tolower(gevMaxima) == "annual") {
      gevParams$timeDelta <- 365.25
      gevParams$timeDeltaYears <- 1
    } else if(tolower(gevMaxima) == "monthly") {
      gevParams$timeDelta <- 365.25/12
      gevParams$timeDeltaYears <- 1/12
    }
    
    gevObj=list()
    gevObj$method = "No fit";
    gevObj$parameters = gevParams;
    gevObj$paramErr = NULL;
    gevObj$stationaryParams = NULL;
    gevObj$objs.monthlyMaxIndexes = NULL;
  }
  
  # estimating the non stationary GPD parameters
  # !!! Assuming a Gaussian approximation to compute the standard errors for
  # the GPD parameters
  if (eva[[2]][[2]]$method!="No fit") {
    epsilonPotX <- eva[[2]][[2]]$parameters[2]
    errEpsilonPotX <- epsilonPotX - eva[[2]][[2]]$paramCIs[1,1]
    sigmaPotX <- eva[[2]][[2]]$parameters[1]
    errSigmaPotX <- sigmaPotX - eva[[2]][[2]]$paramCIs[1, 2]
    thresholdPotX = eva[[2]][[2]]$parameters[3];
    errThresholdPotX = eva[[2]][[2]]$thresholdError;
    nPotPeaks = eva[[2]][[2]]$parameters[5];
    percentilePotX = eva[[2]][[2]]$parameters[6];
    
    dtPeaks = minPeakDistance;
    timeStamps=as.Date(timeStamps)
    dtPotX = as.numeric(timeStamps[length(timeStamps)] - timeStamps[1])/length(series)*dtPeaks;
    epsilonPotNS = epsilonPotX;
    errEpsilonPotNS = errEpsilonPotX;
    sigmaPotNS = sigmaPotX*trasfData$stdDevSeries;
    
    # propagating the errors on stdDevSeries and sigmaPotX to sigmaPotNs.
    # err(sigmaNs) = sqrt{ [sigmaX*err(stdDev)]^2 + [stdDev*err(sigmaX)]^2 }
    # the error on sigmaGevNs is time dependant.
    errSigmaPotFit = trasfData$stdDevSeries*errSigmaPotX;
    errSigmaPotTransf = sigmaPotX*trasfData$stdDevError;
    errSigmaPotNS = (  errSigmaPotTransf^2   +  errSigmaPotFit^2  )^.5;
    thresholdPotNS = thresholdPotX*trasfData$stdDevSeries + trasfData$trendSeries;
    # propagating the errors on stdDevSeries and trendSeries to thresholdPotNs.
    # err(thresholdPotNs) = sqrt{ [thresholdPotX*err(stdDev)]^2 + err(trend)^2 }
    # the error on thresholdPotNs is constant.
    thresholdErrFit = 0;
    
    thresholdErrTransf = ((trasfData$stdDevSeries*errThresholdPotX)^2 + (thresholdPotX*trasfData$stdDevError)^2  +  trasfData$trendError^2)^.5;
    thresholdErr = thresholdErrTransf;
    
    potParams=c()
    potParams$epsilon = epsilonPotNS;
    potParams$sigma = sigmaPotNS;
    potParams$threshold = thresholdPotNS;
    potParams$percentile = percentilePotX;
    potParams$timeDelta = dtPotX;
    potParams$timeDeltaYears = dtPotX/365.25;
    potParams$timeHorizonStart = min(trasfData$timeStamps);
    potParams$timeHorizonEnd = max(trasfData$timeStamps);
    potParams$peaks=trasfData$nonStatSeries[pointData$POT$ipeaks]
    potParams$peakID=pointData$POT$ipeaks
    potParams$peakST=pointData$POT$stpeaks
    potParams$peakEN=pointData$POT$endpeaks
    potParams$nPeaks = nPotPeaks;
    
    
    potParamStdErr=c()
    potParamStdErr$epsilonErr = errEpsilonPotNS;
    potParamStdErr$sigmaErrFit = errSigmaPotFit;
    potParamStdErr$sigmaErrTransf = errSigmaPotTransf;
    potParamStdErr$sigmaErr = errSigmaPotNS;
    potParamStdErr$thresholdErrFit = thresholdErrFit;
    potParamStdErr$thresholdErrTransf = thresholdErrTransf;
    potParamStdErr$thresholdErr = thresholdErr;
    
    potObj=list()
    potObj$method = eva[[2]][[2]]$method;
    potObj$parameters = potParams;
    potObj$paramErr = potParamStdErr;
    potObj$stationaryParams = eva[[2]][[2]];
    potObj$objs = NULL;
  }else{
    
    dtPeaks = minPeakDistance;
    timeStamps=as.Date(timeStamps)
    dtPotX = as.numeric(timeStamps[length(timeStamps)] - timeStamps[1])/length(series)*dtPeaks;
    thresholdPotX = pointData$POT$threshold
    thresholdPotNS = thresholdPotX*trasfData$stdDevSeries + trasfData$trendSeries;
    
    epsilonPotX <- pointData$POT$pars[2]
    sigmaPotX <- pointData$POT$pars[1]
    epsilonPotNS = epsilonPotX;
    sigmaPotNS = sigmaPotX*trasfData$stdDevSeries;
    thresholdPotNS = thresholdPotX*trasfData$stdDevSeries + trasfData$trendSeries;
    
    potParams=c()
    potParams$epsilon = epsilonPotNS;
    potParams$sigma = sigmaPotNS;
    potParams$threshold = thresholdPotNS;
    potParams$percentile = pointData$POT$percentile;
    potParams$timeDelta = dtPotX;
    potParams$timeDeltaYears = dtPotX/365.2425;
    potParams$timeHorizonStart = min(trasfData$timeStamps);
    potParams$timeHorizonEnd = max(trasfData$timeStamps);
    potParams$peaks=trasfData$nonStatSeries[pointData$POT$ipeaks]
    potParams$peakID=pointData$POT$ipeaks
    potParams$peakST=pointData$POT$stpeaks
    potParams$peakEN=pointData$POT$endpeaks
    potParams$nPeaks = length(pointData$POT$peaks);
    
    potObj=list()
    potObj$method = "No fit";
    potObj$parameters = potParams;
    potObj$paramErr = NULL;
    potObj$stationaryParams = NULL;
    potObj$objs = NULL;
  }
  
  # setting output objects
  nonStationaryEvaParams <- list(gevObj=gevObj, potObj=potObj)
  stationaryTransformData <- trasfData
  return(list(nonStationaryEvaParams=nonStationaryEvaParams,stationaryTransformData=stationaryTransformData))
  
}
##########################################################################################






main_path = "D:/tilloal/Documents/06_Floodrivers/"
valid_path = paste0(main_path)
dis_path<-paste0(main_path,'HERA_CFS/')
dis_path<-paste0(main_path,'HERA/')
dis_path<-paste0(main_path,'HERA_Rstat/')
setwd(valid_path)


# # Data generation -----------------------------
# #Load all HERA simulated discharge and create a single file
# #create one file with all years
# y=1951
# Q_d <- read.csv(paste0(dis_path,'HERA_CFR6h_H07_', y, '.csv'), header = F)  # CSVs with observations
# Station_data_IDs <- as.vector(t(Q_d[1, ]))[-1]
# Q_sim <- Q_d
# yrlist=c(1952:2020)
# for (y in yrlist){
#   print(y)
#   Q_d <- read.csv(paste0(dis_path,'HERA_CFR6h_H07_', y, '.csv'), header = F)  # CSVs with observations
#   Q_s <- Q_d[-1, ]
#   Q_sim=rbind(Q_sim,Q_s)
# }
# 
# save(Q_sim,file=paste0(dis_path,"HERA_CFR6h_H07_19502020.Rdata"))



#load(file=paste0(dis_path,"HERA_CFS6h_H07_19502020.Rdata"))
# load(file=paste0(dis_path,"HERA6h_H07_19502020.Rdata"))
# #load(file=paste0(dis_path,"HERACFR6h_H07_19502020.Rdata"))
# 
# 
# Station_data_IDs <- as.vector(t(Q_sim[1, ]))[-1]
# 
# tqx=as.POSIXct((Q_sim[,1]*3600*24), origin="1979-01-01 00:00:00")
# txx=tqx[-1]
#######################  Arguments importation #############################
haz = "drought"
var = "dis"
outlets="Hybas07"
outletname = "outletsv8_hybas07_01min"
season="year"
hydroDir<-("D:/tilloal/Documents/LFRuns_utils/data")
transftypes=c("ori","rev","inv","lninv")
trendtypes=c("trend","trendPeaks")
trendtrans=expand.grid(transftypes,trendtypes)
if (haz=="drought") tt=6
if (haz=="flood") tt=5
print(trendtrans[tt,])
trans=trendtrans[tt,1]
######################################################################################

#Load ancilliary inputs
#load frost days
frostcat=readRDS(file="D:/tilloal/Documents/LFRuns_utils/data/Drought/catchment_frost.RDS")

#load reservoir influence
reservoirdog=outletopen(paste0(hydroDir,"/reservoirs"),"res_ratio_diff_2020-1951")
reservoirdog$llcoord=paste(round(reservoirdog$Var1,4),round(reservoirdog$Var2,4),sep=" ") 
#Hybas07
Catchmentrivers7=read.csv(paste0(hydroDir,"/Catchments/from_hybas_eu_onlyid.csv"),encoding = "UTF-8", header = T, stringsAsFactors = F)
outhyb07=outletopen(hydroDir,outletname)
outhyb07$llcoord=paste(round(outhyb07$Var1,4),round(outhyb07$Var2,4),sep=" ") 
catmatch=match(outhyb07$outlets,Catchmentrivers7$pointid)
mycat=Catchmentrivers7[catmatch,]

hybas07 <- read_sf(dsn = paste0(hydroDir,"/Catchments/hydrosheds/hybas_eu_lev07_v1c.shp"))
hybasf7=fortify(hybas07) 
Catamere07=inner_join(hybasf7,mycat,by= "HYBAS_ID")
Catamere07$llcoord=paste(round(Catamere07$POINT_X,4),round(Catamere07$POINT_Y,4),sep=" ") 
Catf7=Catamere07

st_geometry(Catf7)=NULL


#Ok retry the plot with these IDs

reservoir_bg=inner_join(reservoirdog,outhyb07,by="llcoord")

#load the threshold table

load("D:/tilloal/Documents/LFRuns_utils/TrendAnalysis/thresholds_catchments.rdata")
load("D:/tilloal/Documents/LFRuns_utils/TrendAnalysis/thresholds_low_catchments.rdata")
thtable=dx
#scenari=c("Histo","RWstat","SocCF")
scenari=c("RWstat","SocCF")
for (sce in scenari){
  print(sce)
  if (sce=="Histo"){
    dis_path<-paste0(main_path,'HERA/')
    load(file=paste0(dis_path,"HERA6h_H07_19502020.Rdata"))
    Station_data_IDs <- as.vector(t(Q_sim[1, ]))[-1]
    tqx=as.POSIXct((Q_sim[,1]*3600*24), origin="1979-01-01 00:00:00")
    txx=tqx[-1]
    
    #remove locations without discharge
    wolala=which(is.na(Q_sim[2,]))-1
    Q_sim2=Q_sim[,-which(is.na(Q_sim[2,]))]
    Station_data_IDs=Station_data_IDs[-wolala]
    
  }
  if(sce=="SocCF"){
    dis_path<-paste0(main_path,'HERA_CFS/')
    load(file=paste0(dis_path,"HERA_CFS6h_H07_19502020.Rdata"))
    Station_data_IDs <- as.vector(t(Q_sim[1, ]))[-1]
    tqx=as.POSIXct((Q_sim[,1]*3600*24), origin="1979-01-01 00:00:00")
    txx=tqx[-1]
    
    #remove locations without discharge
    wolala=which(is.na(Q_sim[2,]))-1
    Q_sim2=Q_sim[,-which(is.na(Q_sim[2,]))]
    Station_data_IDs=Station_data_IDs[-wolala]
  }
  if(sce=="RWstat"){
    dis_path<-paste0(main_path,'HERA_RWstat/')
    load(file=paste0(dis_path,"HERA_CFR6h_H07_19502020.Rdata"))
    Station_data_IDs <- as.vector(t(Q_sim[1, ]))[-1]
    tqx=as.POSIXct((Q_sim[,1]*3600*24), origin="1979-01-01 00:00:00")
    txx=tqx[-1]
    
    #remove locations without discharge
    wolala=which(is.na(Q_sim[2,]))-1
    Q_sim2=Q_sim[,-which(is.na(Q_sim[2,]))]
    Station_data_IDs=Station_data_IDs[-wolala]
  }
  #start the loop
  RetPerGPD=c()
  RetPerGEV=c()
  RetLevGEV=c()
  RetLevGPD=c()
  parlist=c()
  peaklist=c()
  catlist=c()
  IRES=c()
  Reservoir_i=c()
  for (id in 1:length(Station_data_IDs)){
    #id=1993
    id=which(Station_data_IDs==23610)
    Reservoir_alteration=0
    outlet=Q_sim2[1,id+1]
    df.dis=Q_sim2[,id+1][-1]
    stationID=Station_data_IDs[id]
    
    #Use this until i finish to generate csv of hybas outlets
    #remove 1950
    rmv=which(year(txx)==1950)
    df.dis=data.frame(txx,df.dis)
    names(df.dis)[c(1,2)]=c("time","dis")
    #remove 1950 which is not reliable
    if (length(rmv)>0){
      df.dis=df.dis[-rmv,]
    }
    start_time <- Sys.time()
    print(paste0(id,"/",length(Station_data_IDs)))
    catch=Station_data_IDs[id]
    timeStamps=txx
    dt = difftime(timeStamps[2],timeStamps[1],units="days")
    dt= as.numeric(dt)
    
    if (haz=="drought"){
      #seasonality divide: frost vs non frost
      minPeakDistanceInDays=30
      percentile=95
      catmat=Catf7[which(Catf7$pointid==catch),]
      Tcatmat=mycat[which(mycat$HYBAS_ID==catmat$HYBAS_ID),]
      Tcatchment=which(colnames(frostcat)==Tcatmat$pointid)
      data=df.dis
      names(data)=c("date","Qs")
      intermit=interid(data, WindowSize=7)
      #intermit$trdis$Q7[which(intermit$trdis$Q7<1e-4)]=NA
      interflag=intermit$flags[2]
      if (!exists("trans")){trans="rev"}
      print(paste0(trans," transformation used for low flows"))
      if (length(rmv)>0){
        series=data.frame(txx[-rmv],intermit$trdis$Q7)
      } else {
      series=data.frame(txx,intermit$trdis$Q7)
      }
      TrendTh=thtable$common[which(thtable$catch==stationID)]
      #remove frost timesteps, this can be modified to do the anlysis only on frost moments
      if (length(Tcatchment)>0){
        frostserie=data.frame(frostcat[,1],frostcat[,Tcatchment])
        names(frostserie)=c("time","Ta")
        frostserie=frostserie[-rmv,]
        frosttime=which(frostserie[,2]<0)
        
      }else{
        frosttime=NA
      }
      trans="rev"
      tail="low"
      
    }else if (haz=="flood"){
      minPeakDistanceInDays=7
      data=df.dis
      catmat=Catf7[which(Catf7$pointid==catch),]
      names(data)=c("date","Qs")
      # Extract the maximum daily value
      series <- max_daily_value(data)
      # series=data.frame(txx,df.disX$outlets)
      interflag=0
      frosttime=NA
      trans="ori"
      tail="high"
      TrendTh=thtable$common[which(thtable$catch==stationID)]
      print(paste0("trend threshold is: ",TrendTh))
    }
    
    names(series)=c("timestamp","dis")
    # plot(series)
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
    
    if(length(na.omit(series$dis))>1 & interflag<1 & nv>15){
      
      if (haz=="drought" & length(!is.na(frosttime))>1){
        if (season=="nonfrost"){
          print("nonfrost season")
          timeAndSeries$data[frosttime]=NA
        }else if (season=="frost"){
          print("frost season")
          timeAndSeries$data[-frosttime]=NA
        }else {print("season must be frost or nonfrost")}
      }
      tsm=1/dt
      series=timeAndSeries[,2]
      timeWindow = 365.25*30; #time windows in days, the correction is done within the functions
      windowSize=366
      timeStamps=timeAndSeries$timestamp
  
      reservoir_loc=reservoir_bg[match(catch,reservoir_bg$outlets.y),]
      Reservoir_alteration=reservoir_loc$outlets.x
  
    Nonstat<-TsEvaNs(timeAndSeries, timeWindow, transfType=trendtrans[tt,2],
                     minPeakDistanceInDays = minPeakDistanceInDays, tail=tail, 
                     trans=trendtrans[tt,1], TrendTh = TrendTh)
    nonStationaryEvaParams=Nonstat[[1]]
    
    print(paste0("GPD shape: ",round(nonStationaryEvaParams$potObj$parameters$epsilon,2)))
    print(paste0("GEV shape: ",round(nonStationaryEvaParams$gevObj$parameters$epsilon,2)))
    #   
    stationaryTransformData=Nonstat[[2]]
    plot(stationaryTransformData$stdDevSeries)
    which(is.na(stationaryTransformData$stdDevSeries))
    #   
    ExRange= c(min(nonStationaryEvaParams$potObj$parameters$peaks),max(nonStationaryEvaParams$potObj$parameters$peaks))
  
      if (haz=="flood") wr2 <- c(seq(min(ExRange),max(ExRange),length.out=700))
      if (haz=="drought") wr2 <- c(seq(1.1*min(ExRange),0.1*max(ExRange),length.out=700))
  
      # Plot1= tsEvaPlotGPDImageScFromAnalysisObj(wr2, nonStationaryEvaParams, stationaryTransformData, minYear = '1950',trans=trendtrans[tt,1])
      # timeIndex=1
      # Plot2 = tsEvaPlotReturnLevelsGPDFromAnalysisObj(nonStationaryEvaParams, stationaryTransformData, timeIndex, trans=trendtrans[tt,1],ylabel="Discharge (m3/s)",ope=T)
      # print(Plot2)
      #ggsave(paste0(hydroDir,"/TSEVA_hybas/",haz,"/GPDBeam_",trans,"_Cat",stationID,".jpg"),Plot2, width=24, height=20, units=c("cm"),dpi=300)
      #Plot3 = tsEvaPlotReturnLevelsGEVFromAnalysisObj(nonStationaryEvaParams, stationaryTransformData, timeIndex, trans="rev",ylabel="Discharge (m3/s)")
  
      stationaryTransformData$timeStampsDay=unique(as.Date(stationaryTransformData$timeStamps))
      pikos=data.frame(nonStationaryEvaParams$potObj$parameters$peaks,nonStationaryEvaParams$potObj$parameters$peakID,nonStationaryEvaParams$potObj$parameters$peakST, nonStationaryEvaParams$potObj$parameters$peakEN)
      names(pikos)=c("value","timeID","tIDstart","tIDend")
      pikos$time=stationaryTransformData$timeStamps[pikos$timeID]
      pikos$catch=rep(catch,length(pikos[,1]))
      #Here I need to convert the timeStamp to a daily one if dt is not 1
      dt1=min(diff(timeStamps),na.rm=T)
      dt=as.numeric(dt1)
      tdim=attributes(dt1)$units
      if (tdim=="hours") dt=dt/24
      if (dt==1){
        timeDays=stationaryTransformData$timeStamps
      }else{
        timeDays=stationaryTransformData$timeStampsDay
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
      RLevs100=ComputeReturnLevels(nonStationaryEvaParams, RPgoal, timeIndex,trans=trendtrans[tt,1])
  
      if (RLevs100$Fit=="No fit"){
  
        RLgpd=nonStationaryEvaParams$gevObj$parameters$annualMax
        names(RLgpd)=year(Impdates)
  
        RLgev=nonStationaryEvaParams$gevObj$parameters$annualMax
        names(RLgev)=year(Impdates)
  
        nRPgev=rep(NA, length(Impdates))
        names(nRPgev)=year(Impdates)
  
        nRPgpd=rep(NA, length(Impdates))
        names(nRPgpd)=year(Impdates)
  
        params=data.frame(matrix(ncol=17,nrow=(length(Impdates)-1)))
        params[,1]=rep(catch,7)
        params[,2]=year(Impdates)[-1]
        params[,4]=rep(interflag,7)
        if (is.null(colnames(parlist))){
          colnames(params)=rep("nom",17)
        }else{
          colnames(params)=colnames(parlist)
        }
      }else{
        #####
        RLgev=RLevs100$ReturnLevels[2]
        RLgpd=RLevs100$ReturnLevels[3]
        ERgev=RLevs100$ReturnLevels[4]
        ERgpd=RLevs100$ReturnLevels[5]
        nRPgev=nRPgpd=10
        params=c()
        for (t in 2:length(Impdates)){
          timeIndex=tindexes[t]
          RLevs100i=ComputeReturnLevels(nonStationaryEvaParams, RPgoal, timeIndex,trans=trendtrans[tt,1])
          params=c(catch,year(Impdates[t]),timeIndex,RLevs100i$Params)
          names(params)[1:3]=c("catchment","Year","timeIndex")
  
          Rper=RPcalc(params,RPiGEV=RLevs100$ReturnLevels[2],RPiGPD=RLevs100$ReturnLevels[3])
          nRPgpd=c(nRPgpd,Rper[2])
          nRPgev=c(nRPgev,Rper[1])
          RLgev=cbind(RLgev,RLevs100i$ReturnLevels[2])
          RLgpd=cbind(RLgpd,RLevs100i$ReturnLevels[3])
          ERgev=cbind(ERgev,RLevs100i$ReturnLevels[4])
          ERgpd=cbind(ERgpd,RLevs100i$ReturnLevels[5])
          if (length(parlist)>1) colnames(parlist)=names(params)
          parlist=rbind(parlist,params)
        }
  
        RLgev=as.data.frame(RLgev)
        names(RLgev)=year(Impdates)
        rownames(RLgev)=RPgoal
  
        RLgpd=as.data.frame(RLgpd)
        names(RLgpd)=year(Impdates)
        rownames(RLgpd)=RPgoal
  
        nRPgev=as.data.frame(t(nRPgev))
        names(nRPgev)=year(Impdates)
  
        nRPgpd=as.data.frame(t(nRPgpd))
        names(nRPgpd)=year(Impdates)
  
        peaklist=rbind(peaklist,pikos)
      }
    }else{
      cat(paste0("\n No values in this pixel \n or intermittent river (flag = ",interflag,")"))
      if (is.na(interflag)) interflag=-9999
      if (interflag>0){
        filling=intermit$flags[3]
        datex=yday(intermit$DaysBlow$time)
        dtect=c(diff(datex),-1)
        last_days <- intermit$DaysBlow$time[which(dtect<0)]
        tindexes=match(last_days,intermit$DaysBlow$time)
        oops=intermit$DaysBlow[tindexes,]
        # verif=series$timestamp[which(!is.na(match(year(series$timestamp),year(Impdates))))]
        # dfrep=data.frame(year(verif),oops)
        # yagg=aggregate(list(m=dfrep$oops),
        #                by = list(year=dfrep$year.verif.),
        #                FUN = function(x) c(max=max(x)))
  
      }else{
        filling=NA
      }
      RLgpd=oops$RP
      names(RLgpd)=year(Impdates)
  
      RLgev=oops$RP
      names(RLgev)=year(Impdates)
  
      nRPgev=rep(NA, length(Impdates))
      names(nRPgev)=year(Impdates)
  
      nRPgpd=rep(NA, length(Impdates))
      names(nRPgpd)=year(Impdates)
  
      params=data.frame(matrix(ncol=17,nrow=(length(Impdates)-1)))
      params[,1]=rep(catch,length(Impdates)-1)
      params[,2]=year(Impdates)[-1]
      params[,4]=rep(interflag,length(Impdates)-1)
      if (is.null(colnames(parlist))){
        colnames(params)=rep("nom",17)
      }else{
        colnames(params)=colnames(parlist)
      }
      parlist=as.data.frame((rbind(parlist,params)))
  
      pikos=data.frame(matrix(ncol=6,nrow=1))
      pikos[,1]=NA
      pikos[,6]=catch
      if (is.null(colnames(peaklist))){
        colnames(pikos)=c("value","timeID","tIDstart","tIDend","time","catch")
      }else{
        colnames(pikos)=colnames(peaklist)
      }
      #parlist=as.data.frame(mapply(c,parlist,params))
      peaklist=as.data.frame((rbind(peaklist,pikos)))
  
    }
  
    #Saving main outputs
    catlist=c(catlist,catch)
    IRES=c(IRES,interflag)
    Reservoir_i=c(Reservoir_i,Reservoir_alteration)
    RetLevGEV=rbind(RetLevGEV,RLgev)
    RetLevGPD=rbind(RetLevGPD,RLgpd)
    RetPerGEV=rbind(RetPerGEV,nRPgev)
    RetPerGPD=rbind(RetPerGPD,nRPgpd)
    end_time <- Sys.time()
    cat(paste0("\nloop duration: ",round(end_time-start_time,2)," seconds\n"))
  }
  
  parlist=as.data.frame(parlist)
  Results=list(parameters=parlist,RetLevGEV=RetLevGEV,RetLevGPD=RetLevGPD,RetPerGEV=RetPerGEV,RetPerGPD=RetPerGPD,Peaks=peaklist,catrest=data.frame(catlist,Reservoir_i,IRES))
  Results$Peaks
  
  #save(parlist,file=paste0(hydroDir,"/TSEVA_hybas/outputs/parCat07_",sce,"_",haz,"_1950_2020.Rdata"))
  save(Results, file=paste0(hydroDir,"/TSEVA_hybas/outputs/ResCat07_",sce,"_",haz,"v2_1951_2020.Rdata"))

# save(parlist,file=paste0(hydroDir,"/TSEVA_hybas/outputs/parCat07_ResCF_",haz,"_1950_2020.Rdata"))
# save(Results, file=paste0(hydroDir,"/TSEVA_hybas/outputs/ResCat07_ResCF_",haz,"_1950_2020.Rdata"))

}


