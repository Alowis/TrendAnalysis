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

hybas07 <- read_sf(dsn = paste0(hydroDir,"/Catchments/hydrosheds/hybas_eu_lev07_v1c.shp"))
hybasf7=fortify(hybas07) 
#  Function declaration ---------------------------------------------------


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








# Function that open netcdf outlet files
outletopen=function(dir,outletname,nrspace=rep(NA,5)){
  ncbassin=paste0(dir,"/",outletname,".nc")
  ncb=nc_open(ncbassin)
  name.vb=names(ncb[['var']])
  namev=name.vb[1]
  if ("Band1"%in% name.vb)namev="Band1"
  name.lon="lon"
  name.lat="lat"
  if (!is.na(nrspace[1])){
    start=as.numeric(nrspace[c(2,4)])
    count=as.numeric(nrspace[c(3,5)])-start+1
  }else{
    londat = ncvar_get(ncb,name.lon) 
    llo=length(londat)
    latdat = ncvar_get(ncb,name.lat)
    lla=length(latdat)
    start=c(1,1)
    count=c(llo,lla)
  }
  
  londat = ncvar_get(ncb,name.lon,start=start[1],count=count[1]) 
  llo=length(londat)
  latdat = ncvar_get(ncb,name.lat,start=start[2],count=count[2])
  lla=length(latdat)
  outlets = ncvar_get(ncb,namev,start = start, count= count) 
  outlets=as.vector(outlets)
  outll=expand.grid(londat,latdat)
  lonlatloop=expand.grid(c(1:llo),c(1:lla))
  outll$idlo=lonlatloop$Var1
  outll$idla=lonlatloop$Var2
  
  outll=outll[which(!is.na(outlets)),]
  outlets=outlets[which(!is.na(outlets))]
  outll=data.frame(outlets,outll)
  return (outll)
}

#Function to compute return levels associated to a RP (not used in this script)
ComputeReturnLevels<-function(nonStationaryEvaParams, RPgoal, timeIndex){
  #GEV
  epsilonGEV <- nonStationaryEvaParams[[1]]$parameters$epsilon
  sigmaGEV <- mean(nonStationaryEvaParams[[1]]$parameters$sigma[timeIndex])
  muGEV <- mean(nonStationaryEvaParams[[1]]$parameters$mu[timeIndex])
  dtSampleYears <- nonStationaryEvaParams[[1]]$parameters$timeDeltaYears
  epsilonStdErrGEV <- nonStationaryEvaParams[[1]]$paramErr$epsilonErr
  sigmaStdErrGEV <- mean(nonStationaryEvaParams[[1]]$paramErr$sigmaErr[timeIndex])
  muStdErrGEV <- mean(nonStationaryEvaParams[[1]]$paramErr$muErr[timeIndex])
  
  #GPD
  epsilonGPD <- nonStationaryEvaParams[[2]]$parameters$epsilon
  sigmaGPD <- mean(nonStationaryEvaParams[[2]]$parameters$sigma[timeIndex])
  thresholdGPD <- mean(nonStationaryEvaParams[[2]]$parameters$threshold[timeIndex])
  nPeaks <- nonStationaryEvaParams[[2]]$parameters$nPeaks
  epsilonStdErrGPD <- nonStationaryEvaParams[[2]]$paramErr$epsilonErr
  sigmaStdErrGPD <- mean(nonStationaryEvaParams[[2]]$paramErr$sigmaErr[timeIndex])
  thresholdStdErrGPD <- mean(nonStationaryEvaParams[[2]]$paramErr$thresholdErr[timeIndex])
  thStart <- nonStationaryEvaParams[[2]]$parameters$timeHorizonStart
  thEnd <- nonStationaryEvaParams[[2]]$parameters$timeHorizonEnd
  sampleTimeHorizon <- as.numeric((thEnd - thStart)/365.2425)
  
  returnLevelsGEV <- tsEvaComputeReturnLevelsGEV(epsilonGEV, sigmaGEV, muGEV, epsilonStdErrGEV, sigmaStdErrGEV, muStdErrGEV, RPgoal)
  
  returnLevelsGPD <- tsEvaComputeReturnLevelsGPD(epsilonGPD, sigmaGPD, thresholdGPD, epsilonStdErrGPD, sigmaStdErrGPD, thresholdStdErrGPD,
                                                 nPeaks, sampleTimeHorizon, RPgoal)
  rlevGEV=returnLevelsGEV$returnLevels
  rlevGPD=returnLevelsGPD$returnLevels
  
  ParamGEV=c(epsilonGEV,sigmaGEV,muGEV,epsilonStdErrGEV, sigmaStdErrGEV, muStdErrGEV)
  names(ParamGEV)=c("epsilonGEV","sigmaGEV","muGEV","epsilonStdErrGEV","sigmaStdErrGEV","muStdErrGEV")
  
  ParamGPD=c(epsilonGPD,sigmaGPD,thresholdGPD,epsilonStdErrGPD,sigmaStdErrGPD, thresholdStdErrGPD,nPeaks,sampleTimeHorizon)
  names(ParamGPD)=c("epsilonGPD","sigmaGPD","thresholdGPD","epsilonStdErrGPD","sigmaStdErrGPD","thresholdStdErrGPD","nPeaks","SampleTimeHorizon")
  
  return(list(ReturnLevels=c(ReturnPeriod=RPgoal, GEV=as.numeric(rlevGEV),GPD=as.numeric(rlevGPD)),Params=c(ParamGEV,ParamGPD)))
  
}

#function to obtain RL for amy return period
GPDLargeRLs<-function(Paramsfl, RPgoal){
  #GPD
  epsilonGPD <- Paramsfl$epsilonGPD
  sigmaGPD <-Paramsfl$sigmaGPD
  thresholdGPD <- Paramsfl$thresholdGPD
  nPeaks <- Paramsfl$nPeaks
  epsilonStdErrGPD <- Paramsfl$epsilonStdErrGPD
  sigmaStdErrGPD <- Paramsfl$sigmaStdErrGPD
  thresholdStdErrGPD <- Paramsfl$thresholdStdErrGPD
  sampleTimeHorizon<-Paramsfl$SampleTimeHorizon
  
  X0 <- nPeaks/SampleTimeHorizon
  RPiGPD=RPgoal
  qXGPD=(X0*RPiGPD)
  
  # if (epsilonGPD != 0) {
  # estimating the return levels
  returnLevels <- thresholdGPD + sigmaGPD/epsilonGPD*((qXGPD)^epsilonGPD - 1)
  
  # suicide=tsEvaComputeReturnLevelsGPD(epsilonGPD[1], sigmaGPD[1],thresholdGPD[1],epsilonStdErrGPD[1], sigmaStdErrGPD[1],
  #                                     thresholdStdErrGPD[1], nPeaks[1], sampleTimeHorizon[1], 100)
  # estimating the error
  # estimating the differential of returnLevels to the parameters
  # !! ASSUMING NON ZERO ERROR ON THE THRESHOLD AND 0 ERROR ON THE PERCENTILE.
  # THIS IS NOT COMPLETELY CORRECT BECAUSE THE PERCENTILE DEPENDS ON THE
  # THRESHOLD AND HAS THEREFORE AN ERROR RELATED TO THAT OF THE
  # THRESHOLD
  dxm_u <- 1
  dxm_sigma <- 1/epsilonGPD*(qXGPD^epsilonGPD - 1)
  dxm_epsilon <- -sigmaGPD/epsilonGPD^2*((qXGPD)^epsilonGPD - 1) + sigmaGPD/epsilonGPD*log(qXGPD)*qXGPD^epsilonGPD
  
  returnLevelsErr <- sqrt((dxm_u*thresholdStdErrGPD)^2 + (dxm_sigma*sigmaStdErrGPD)^2 + (dxm_epsilon*epsilonStdErrGPD)^2)
  #
  # } else {
  #   returnLevels <- thresholdGPD + sigmaGPD*log(qXGPD)
  #   # estimating the error
  #   # estimating the differential of returnLevels to the parameters
  #   # !! ASSUMING NON ZERO ERROR ON THE THRESHOLD, 0 ERROR ON THE
  #   # PERCENTILE AND 0 ERROR ON EPSILON.
  #   # THIS IS NOT COMPLETELY CORRECT BECAUSE THE PERCENTILE DEPENDS ON THE
  #   # THRESHOLD AND HAS THEREFORE AN ERROR RELATED TO THAT OF THE
  #   # THRESHOLD
  #   dxm_u <- 1
  #   dxm_sigma <- log(qXGPD)
  #   
  #   returnLevelsErr <-(  (dxm_u*thresholdStdErrGPD)^2  +  (dxm_sigma*sigmaStdErrGPD)^2  )^.5;
  # }
  # 
  
  return(list(ReturnLevels=c(ReturnPeriod=RPgoal),GPD=data.frame(returnLevels,returnLevelsErr)))
  
}

#Function to compute return periods associated to a return level. Not used
RPcalc<- function(params,RPiGEV,RPiGPD){
  paramx=data.frame(t(params))
  qxV=1-exp(-(1+paramx$epsilonGEV*(RPiGEV-paramx$muGEV)/paramx$sigmaGEV)^(-1/paramx$epsilonGEV))
  returnPeriodGEV=1/qxV
  X0 <- paramx$nPeaks/paramx$SampleTimeHorizon
  qxD=(((1+paramx$epsilonGPD*(RPiGPD-paramx$thresholdGPD)/paramx$sigmaGPD)^(-1/paramx$epsilonGPD)))
  returnPeriodGPD=1/(X0*qxD)
  
  return(c(GEV=returnPeriodGEV,GPD=returnPeriodGPD))
}

#Initial functions to plot changes and trends at catchment level
plotchangemaps=function(basemap,catmap,datar,law="GPD",type,period=c(1951,1952),parlist,valuenames,haz){
  data=right_join(catmap,datar,by = c("outlets"="unikout"))
  data$outlets
  data$pointid
  datacol=names(data)
  valcol=match(valuenames,datacol)
  datatwin=data
  valcol2=valcol-1
  st_geometry(datatwin) <- NULL
  datatwin=as.data.frame(datatwin)
  class(datatwin)
  dtc=names(datatwin)
  mcor=unique(match(datar$unikout,parlist$catchment))
  
  if (length(which(is.na(mcor)))>0) datar=datar[-which(is.na(mcor)),]
  
  cref=paste0("X",period[1])
  crefloc=match(cref,dtc)
  finalperiod=paste0("X",period[2])
  colsel=match(finalperiod,datacol)
  
  if (type=="RLchange"){
    palet=c(hcl.colors(11, palette = "RdYlBu", alpha = NULL, rev = F, fixup = TRUE))
    title=paste0("Relative change in 10-years ",haz," Return Level between ", period[1], " and ", period[2])
    legend="Relative Change (%)"
    datathin=datatwin[,c(valcol2)]
    tmpval=(datatwin[,valcol2])/(datatwin[,crefloc])*100-100
    if (haz=="drought"){
      for (it in 1:length(tmpval[1,])){
        tmpval[,it][which(data$IRES==1)]=NA
        #tmpval[,it][which(data$IRES==1)]=datathin[which(data$IRES==1),it]/datathin[which(data$IRES==1),crefloc]
        # tmpval[which(data$IRES==1),it][which(is.nan(tmpval[which(data$IRES==1),it]))]=0 #drying did not happen in both periods
        # tmpval[which(data$IRES==1),it][which(is.infinite(tmpval[which(data$IRES==1),it]))]=-100 #drying only happens in recent period
        # tmpval[which(data$IRES==1),it][which(tmpval[which(data$IRES==1),it]==0)]=100 #drying only happens in initial period
        # tmpval[,it][which(is.infinite(tmpval[,it]))]=100
        # tmpval[,it][which(is.na(tmpval[,it]))]=0
        
        tmpval[,it][which(is.infinite(tmpval[,it]))]=100
        #tmpval[,it][which(is.na(tmpval[,it]))]=0
        #more stuff
        #tmpval[,it][which(is.infinite(datatwin[,it]))]=100
      }
      oyea=as.data.frame(datatwin[which(data$IRES==1),c(1,valcol2)])
      oyeb=as.data.frame(datatwin[which(datatwin[,crefloc]==0),c(1,valcol2)])
      oyval=rbind(oyea,oyeb)
    }
    #tmpval[which(data$IRES==1),]=NA
    data[,valcol]=tmpval
    br=c(-100,-80,-60,-40,-20,0,20,40,60,80,100)
    limi=c(-100,100)
    trans=scales::modulus_trans(.8)
    colNA="green"
  }
  if (type=="RPchange"){
    palet=c(hcl.colors(9, palette = "BrBG", alpha = NULL, rev = F, fixup = TRUE))
    title=paste0(period[2]," return period of a 100-years Return Level ",haz," in ", period[1])
    legend="Return Period (Years)"
    br=c(0.1,1,10,100,1000)
    limi=c(0.1,1000)
    trans="log"
    newRP=RPchangeCal(parlist, yi=period[1], yf=period[2], RetLev=datar,law, valuenames)
    ow=match(data$outlets,newRP$catchment)
    oyea=data[which(data$IRES==1),colsel]
    st_geometry(oyea) <- NULL
    
    data[,colsel]=newRP$newRP[ow]
    data[which(data$IRES==1),colsel]=oyea
    colNA="red"
    
  }
  
  names(data)[colsel]="fillplot"
  pl1=ggplot(basemap) +
    geom_sf(fill="white")+
    geom_sf(data=data,aes(fill=fillplot,geometry=geometry),color="transparent")+ 
    geom_sf(fill=NA, color="grey") +
    scale_fill_gradientn(
      colors=palet,
      breaks=br,limits=limi,trans=trans,
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
  
  
  #Second plot showing IRES, dam influenced rivers and zero flow rivers
  
  if (haz=="drought"){
    dt2=data[which(data$Reservoir_i>0.1),]
    dt3=oyea
    dt3=data[which(data$IRES==1),]
    dt3$flag=1
    dt4=data[which(datatwin[,crefloc]==0),]
    dt4$flag=2
    dt5=data[which(datatwin[,colsel-1]==0),]
    dt5$flag=3
    dtf=rbind(dt3,dt4,dt5)
    dtff=aggregate(list(fl=dtf$flag),
                   by = list(ID=dtf$HYBAS_ID),
                   FUN = function(x) c(sum=sum(x)))
    length(unique(dtf$HYBAS_ID))
    
    dtf=right_join(dtf,dtff,by=c("HYBAS_ID"="ID"))
    dtf=dtf[,c(1,18,19,97:99)]
    und=match(unique(dtf$HYBAS_ID),dtf$HYBAS_ID)
    dtf=dtf[und,]
    factor(dtf$fl)
    dtf$flv=factor(dtf$fl)
    
    pl2=ggplot(basemap) +
      geom_sf(fill="white")+
      geom_sf(data=data,aes(geometry=geometry),fill="transparent",color="grey")+ 
      geom_sf(fill=NA, color="grey") +
      geom_sf(data=dt2,aes(geometry=geometry,color="darkblue"),fill="transparent")+ 
      geom_sf(data=dtf,aes(geometry=geometry,fill=flv),color="transparent",show.legend = FALSE)+ 
      geom_point(data=dtf,aes(x=POINT_X,y=POINT_Y,fill=flv), shape = 21, size = 0, show.legend = TRUE) +
      scale_fill_manual(values=c("orange","cyan","orangered","grey30"),
                        labels = c("IRES","10yRL[1]=0","10yRL[2]=0","10yRL[1]=10yRL[2]=0","tg"), 
                        name="Low flows")+
      coord_sf(xlim = c(min(nco[,1]),max(nco[,1])), ylim = c(min(nco[,2]),max(nco[,2])))+
      labs(x="Longitude", y = "Latitude")+
      scale_color_identity(name = "Reservoir construction\n1951-2020",
                           guide = "legend",
                           labels= " ")+
      guides(fill = guide_legend(override.aes = list(size=5)))+
      theme(axis.title=element_text(size=tsize),
            panel.background = element_rect(fill = "aliceblue", colour = "grey1"),
            panel.border = element_rect(linetype = "solid", fill = NA, colour="black"),
            legend.title = element_text(size=tsize),
            legend.text = element_text(size=osize),
            legend.position = "right",
            panel.grid.major = element_line(colour = "grey70"),
            panel.grid.minor = element_line(colour = "grey90"),
            legend.key = element_rect(fill = "transparent", colour = "transparent"),
            legend.key.size = unit(.8, "cm"))+
      ggtitle("Reservoirs and zero flow")
  }else{
    pl2=NA
  }
  #save this as a shapefile for visu:
  
  # savecrap=inner_join(data,Catamere07,by= c("HYBAS_ID"))
  # st_write(data, paste0(hydroDir,"/disasterdrought3.shp"))
  # 
  return(list(pl1,pl2))
}



#Initial functions to plot changes and trends at catchment level
haz="drought"
datar=datarD
plotchangemaps_qsp=function(basemap,catmap,datar,upArea, GHR_riv, HydroRsf, law="GPD",type,period=c(1951,2020),parlist,valuenames,haz){
  
  data=right_join(GHR_riv,datar,by = c("outlets"="unikout"))
  data=right_join(catmap,data,by = c("outlets"="outlets"))
  names(data)[28]="HydroR"
  upag=match(data$outlets,UpArea$outlets )
  data$uparea=UpArea$upa[upag]
  data$outlets
  data$pointid
  datacol=names(data)
  valcol=match(valuenames,datacol)
  # data=data[which(data$IRES==1),]
  datatwin=data
  st_geometry(datatwin) <- NULL
  valcol2=valcol-1
  
  datatwin=as.data.frame(datatwin)
  class(datatwin)
  dtc=names(datatwin)
  mcor=unique(match(datar$unikout,parlist$catchment))
  
  if (length(which(is.na(mcor)))>0) datar=datar[-which(is.na(mcor)),]
  
  cref=paste0("X",period[1])
  crefloc=match(cref,dtc)
  finalperiod=paste0("X",period[2])
  colsel=match(finalperiod,datacol)
  
  
  palet=c(hcl.colors(11, palette = "RdYlBu", alpha = NULL, rev = F, fixup = TRUE))
  title=paste0("Change in 10-years ",haz," Return Level between ", period[1], " and ", period[2])
  legend="Change in specific discharge (l/s/km2)"
  datathin=datatwin[,c(valcol2)]
  #value_L_s_m2 <- (depth_mm_day / 1000) / 86400 * 1000
  tmpval=(datatwin[,valcol2]-datatwin[,crefloc])
  tmpval=tmpval*1000/(data$uparea)
  if (haz=="drought"){
    data[,valcol]=tmpval*100
    br=c(-150,-100,-50,0,50,100,150)
    labels=br/100
    limi=c(-150,150)
  }
  if (haz=="flood"){
    data[,valcol]=tmpval
    br=c(-30,-20,-10,0,10,20,30)
    labels=br
    limi=c(-30,30)
  }
  trans=scales::modulus_trans(.8)
  colNA="gray10"
  
  names(data)[colsel]="fillplot"
  data[which(data$IRES==1),colsel]=NA
  pl1=ggplot(basemap) +
    geom_sf(fill="white")+
    geom_sf(data=data,aes(fill=fillplot,geometry=geometry),color="transparent")+ 
    geom_sf(fill=NA, color="grey") +
    scale_fill_gradientn(
      colors=palet,
      breaks=br,labels=labels, limits=limi,trans=trans,
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
  
  # savecrap=inner_join(data,Catamere07,by= c("HYBAS_ID"))
  # st_write(data, paste0(hydroDir,"/disasterdrought3.shp"))
  # 
  
  
  pointagg=aggregate(list(Rchange=data$fillplot),
                     by = list(HydroR=data$HydroR),
                     FUN = function(x) c(mean=median(x,na.rm=T),dev=sd(x,na.rm=T),len=length(x),q1=quantile(x, 0.25, na.rm=T),q3=quantile(x, 0.75, na.rm=T)))
  pointagg <- do.call(data.frame, pointagg)
  
  
  #natch pointagg with hybasf
  pointplot=inner_join(HydroRsf,pointagg,by= c("Id"="HydroR"))
  
  pl2=ggplot(basemap) +
    geom_sf(fill="white")+
    geom_sf(data=pointplot,aes(fill=Rchange.mean,geometry=geometry),color="transparent")+ 
    geom_sf(fill=NA, color="grey") +
    scale_fill_gradientn(
      colors=palet,
      breaks=br,labels=labels, limits=limi,trans=trans,
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
  
  
  
  #plot with kriging interpolation
  
  #reduced dtset
  colnames(data)
  dts=data[,c(1,21,22,23,colsel,106)]
  st_geometry(dts) <- NULL
  
  #generate the grid:
  ##############################################################################
  #interpolating the sen-slopes via kriging
  ##############################################################################
  
  library("sp")
  library("rworldmap")
  newmap <- getMap(resolution = "high")
  europe_frame <- data.frame(country=NA,long=NA,lat=NA,plot_group=NA)
  for(j in c(1:nrow(newmap@data))) {
    print(j)
    if(newmap@data$REGION[j]!="Europe"|is.na(newmap@data$REGION[j])) {next}
    else {
      n_polys = length(newmap@polygons[[j]]@Polygons)
      for(k in c(1:n_polys)) {
        tmp_frame <- data.frame(country=rep(newmap@polygons[[j]]@ID,nrow(newmap@polygons[[j]]@Polygons[[k]]@coords)),
                                long=newmap@polygons[[j]]@Polygons[[k]]@coords[,1],
                                lat=newmap@polygons[[j]]@Polygons[[k]]@coords[,2],
                                plot_group = rep(paste(newmap@polygons[[j]]@ID,k,sep="-"),nrow(newmap@polygons[[j]]@Polygons[[k]]@coords)))
        europe_frame <- rbind(europe_frame,tmp_frame)
      }
    }
  }
  europe_frame <- europe_frame[-1,]
  
  
  #in order to perform kriging we need a different coordinate system than long-lat and a grid
  
  laea_proj="+proj=laea +lat_0=52 +lon_0=10 +x_0=4321000 +y_0=3210000 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs +type=crs"
  lambert_proj="+proj=laea +lat_0=55 +lon_0=20 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=km +no_defs"
  
  spatial_krig = SpatialPointsDataFrame(coords=dts[,c("Var1.x", "Var2.x")], data=dts[,c("fillplot","uparea")],proj4string=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"))
  
  spatial_krig_trafo = spTransform(spatial_krig,CRS(laea_proj))
  
  dts$laea_x = spatial_krig_trafo@coords[,1]
  dts$laea_y = spatial_krig_trafo@coords[,2]
  
  europe_krig = SpatialPointsDataFrame(coords=europe_frame[,c("long", "lat")], data=europe_frame[,c(1:4)],proj4string=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"))
  europe_krig_trafo = spTransform(europe_krig,CRS(laea_proj))
  
  europe_frame$laea_x = europe_krig_trafo@coords[,1]
  europe_frame$laea_y = europe_krig_trafo@coords[,2]
  
  
  #cover of europe for spatial grid
  europe_cover = data.frame(x=c(-2700,-2700,2000,2000),y=c(-2500,2000,2000,-2500),data=rep("data",4))
  #another europe cover
  lims=c(xmin=2665274, xmax=6528705, ymin=1464363, ymax=5360738)
  europe_cover = data.frame(x=c(2000000,2000000,7000000,7000000),y=c(1000000,6000000,6000000,1000000),data=rep("data",4))
  
  coordinates(europe_cover) =~x+y
  my_big_grid = makegrid(europe_cover,cellsize=20000)
  keep_point = rep(0,nrow(my_big_grid))
  relevant_countries = unique(europe_frame$plot_group)
  relevant_countries <- relevant_countries[-grep("Rus|Aze|Turk|Ukr|Isra", relevant_countries, ignore.case = TRUE)]
  
  #for the grid, we check which points are contained in any country
  for(country in relevant_countries)  {
    tmp_country = europe_frame[which(europe_frame$plot_group==country),]
    pt_in_country_tmp = point.in.polygon(point.x=my_big_grid[,1], point.y=my_big_grid[,2], pol.x=tmp_country[,"laea_x"], pol.y=tmp_country[,"laea_y"], mode.checked=FALSE)
    if(length(which(pt_in_country_tmp==1))>0) {
      keep_point[which(pt_in_country_tmp==1)] = 1
    }
    print(which(country==relevant_countries)/length(relevant_countries))
  }
  hist(keep_point)
  my_grid = my_big_grid[which(keep_point==1),]
  names(my_grid) = c("lambert_x","lambert_y")
  
  grid_krig = SpatialPointsDataFrame(coords=my_grid[,c(1,2)], data=my_grid[,c(1,2)],proj4string=CRS(laea_proj))
  #kriging
  library("automap")
  kriging_result <- autoKrige(fillplot~1, input_data=spatial_krig_trafo, new_data = grid_krig)
  
  
  # Convert kriging results to a data frame
  kriging_data <- as.data.frame(kriging_result$krige_output)
  
  # Extract the coordinates from the SpatialPointsDataFrame
  coords <- coordinates(kriging_result$krige_output)
  
  # Combine the coordinates with the kriging data
  kriging_data$lon <- coords[, 1]
  kriging_data$lat <- coords[, 2]
  
  
  library(dplyr)
  # Transform back to longlat projection if necessary
  kriging_data_out = SpatialPointsDataFrame(coords=kriging_data[,c("lon", "lat")], data=kriging_data[,c(3:5)],proj4string=CRS(laea_proj))
  kriging_data_longlat <- spTransform(kriging_data_out, CRS("+proj=longlat +datum=WGS84"))
  
  kriging_data_longlat=as.data.frame(kriging_data_out)
  
  
  # Create a ggplot map of the kriging results
  m3<-kriging_data_longlat
  m4 <- m3 %>%
    # create a new variable from count
    mutate(countv1=cut(var1.pred,breaks=c(-100,-24,-12,-5,-2, 0,2, 5, 12,100),
                       labels=c("< -24", "-24 - -12", "-12 - -5", "-5 - -2","-2 - 0","0 - 2", "0 - 5", "5 - 12","> 12"))) %>%
    # change level order
    mutate(countv1=factor(as.character(countv1), levels=rev(levels(countv1))))
  
  outletname="upArea_European_01min"
  outll=outletopen(hydroDir,outletname)
  cord.dec=outll[,c(2,3)]
  cord.dec = SpatialPoints(cord.dec, proj4string=CRS("+proj=longlat"))
  cord.UTM <- spTransform(cord.dec, CRS("+init=epsg:3035"))
  nco=cord.UTM@coords
  coord.lamb=spTransform(cord.dec, CRS(laea_proj))
  nci=coord.lamb@coords
  
  coord.laea=spTransform(cord.dec, CRS(laea_proj))
  nci=coord.laea@coords
  
  world <- ne_countries(scale = "medium", returnclass = "sf")
  Europe <- world[which(world$continent == "Europe"),]
  e2=st_transform(Europe,  crs=3035)
  w2=st_transform(world,  crs=3035)
  basemap=w2
  basemap=st_transform(basemap, CRS(laea_proj))
  
  Recapoint <- st_as_sf(dts, coords = c("Var1.x", "Var2.x"), crs = 4326)
  Recapoint <- st_transform(Recapoint, crs = 3035)
  mp <- Recapoint %>%
    # create a new variable from count
    mutate(countv1=cut(fillplot, breaks=c(-100,-24,-12,-5,-2, 0,2, 5, 12,100),
                       labels=c("< -24", "-24 - -12", "-12 - -5", "-5 - -2","-2 - 0","0 - 2", "0 - 5", "5 - 12","> 12"))) %>%
    # change level order
    mutate(countv1=factor(as.character(countv1), levels=rev(levels(countv1))))
  
  
  tsize=12
  osize=10
  
  paletc=c(hcl.colors(9, palette = "RdYlBu", alpha = NULL, rev = T, fixup = TRUE))
  paletf=c(hcl.colors(7, palette = "RdYlBu", alpha = NULL, rev = T, fixup = TRUE))
  pl3<-ggplot(basemap) +
    geom_sf(fill="white")+
    geom_tile(data=m4, aes(x = lon, y = lat, fill = countv1)) +
    #metR::geom_contour_fill(data=kriging_data_longlat, aes(x = lon, y = lat, z = var1.pred),binwidth=2,alpha=1) +
    geom_sf(fill=NA, color="grey") +
    geom_sf(data=mp,aes(geometry=geometry,color=countv1),size=1.5,alpha=.5,shape=16,stroke=0.5)+ 
    geom_sf(data=mp,aes(geometry=geometry),color="black",size=1.5,alpha=.5,shape=21,stroke=0.5)+ 
    #scale_fill_distiller(palette = "RdYlBu",limits=c(-12,12),oob = scales::squish,guide="coloursteps",direction=1)+
    #scale_fill_gradientn(colors=paletx,limits=c(-12,12),oob = scales::squish) +
    scale_fill_manual(values=paletc, na.value = "grey90")+
    scale_color_manual(values=paletc, na.value = "grey90")+
    coord_sf(xlim = c(min(nci[,1]),max(nci[,1])), ylim = c(min(nci[,2]),max(nci[,2])))+
    labs(x="Longitude", y = "Latitude")+
    # guides(fill = guide_coloursteps(barwidth = 1, barheight = 10))+
    labs(fill = "Change [%] ") +
    guides( color= "none",
            fill = guide_legend(theme = theme(
              legend.title.position = "right")))+
    theme(axis.title=element_text(size=tsize),
          panel.background = element_rect(fill = "aliceblue", colour = "grey1"),
          panel.border = element_rect(linetype = "solid", fill = NA, colour="black"),
          legend.title = element_text(angle = 90 , size=tsize),
          legend.text = element_text(size=osize),
          legend.key.height= unit(1.6, 'cm'),
          legend.key.width= unit(.2, 'cm'),
          legend.position = "right",
          legend.text.position = "left",
          panel.grid.major = element_line(colour = "grey70"),
          panel.grid.minor = element_line(colour = "grey90"),
          legend.key = element_rect(fill = "transparent", colour = "transparent"),
          legend.key.size = unit(.8, "cm"))+
    ggtitle(paste0("Change in 10y flood return level (",period[1],"-",period[2],")"))
  pl3

  
  return(list(pl1, pl2, pl3))
}

plotTrendSig=function(basemap,catmap,datar,type,period=c(1950,2020),valuenames,Impdates){
  
  
  data=right_join(catmap,datar,by = c("outlets"="unikout"))
  datacol=names(data)
  valcol=match(valuenames,datacol)
  datatwin=data
  valcol2=valcol-1
  st_geometry(datatwin) <- NULL
  datatwin=as.data.frame(datatwin)
  class(datatwin)
  dtc=names(datatwin)
  dtreg=datatwin[,c(41:48)]
  # miniTS=ts(data=(as.numeric(dtreg[2,])),start=year(as.Date(paste0(Impdates[1],"-06-01"))),
  #           end=year(as.Date(paste0(Impdates[8],"-06-01"))),frequency = .1)
  # minReg=as.data.frame(miniTS)
  # minReg$time=as.numeric((Impdates))
  # 
  # mk=MannKendall(miniTS)
  # mkt=mk$tau
  # mks=mk$sl
  # plot(miniTS)
  # optn=lm(x~time, data=minReg)
  # abline(optn)
  mcor=unique(match(datar$unikout,parlist$catchment))
  
  if (length(which(is.na(mcor)))>0) datar=datar[-which(is.na(mcor)),]
  
  cref=paste0("Y",period[1])
  crefloc=match(cref,dtc)
  finalperiod=paste0("Y",period[2])
  colsel=match(finalperiod,datacol)
  tmpval=(datatwin[,valcol2])/(datatwin[,crefloc])*100-100
  data[,valcol]=tmpval
  palet=c(hcl.colors(11, palette = "RdBu", alpha = NULL, rev = F, fixup = TRUE))
  title=paste0("Relative change in 100-years flood Return Level (", period[1], "-", period[2],")")
  #legend="Kendall tau"
  legend="Relative change (%)"
  tmpval=(datatwin[,valcol2])
  mkta=c()
  mksa=c()
  for (it in 1:length(tmpval[,1])){
    mks=NA
    mkt=NA
    miniTS=as.numeric(tmpval[it,])
    if (!is.infinite(miniTS[1])){
      #mk=MannKendall(miniTS)
      mk=mmkh(miniTS,ci=0.95)
      
      #compute trend as well with sen.slope
      mkt=mk$tau
      mks=mk$sl
    }else{
      mkt=NA
      mks=NA
    }
    mkta=c(mkta,mkt)
    mksa=c(mksa,mks)
  }
  # br=c(-1,-.80,-.60,-.40,-.20,0,.20,.40,.60,.80,1)
  # limi=c(-1,1)
  trans=scales::modulus_trans(.6)
  colNA="darkgrey"
  #title="Trends in 100 years Re flood magnitude"
  
  data$mkta=mkta
  data$sl=mksa
  names(data)[colsel]="fillplot"
  
  datasl=data[which(data$sl<=0.1),]
  datasl$centroid <- st_centroid(datasl$geometry)
  datasl$siglvl=1
  datasl$siglvl[which(datasl$sl<=0.05)]=2
  st_geometry(datasl) <- "centroid"
  datasig_f=st_as_sf(datasl,coords=c("POINT_X","POINT_Y"),crs=3035)
  datasig_f$siglvl=as.factor(datasig_f$siglvl)
  datasig_f$sign="red"
  datasig_f$sign[which(datasig_f$fillplot>=0)]="blue"
  
  br=c(-100,-80,-60,-40,-20,0,20,40,60,80,100)
  limi=c(-100,100)
  tsize=14
  osize=12
  ocrap<-ggplot(basemap) +
    geom_sf(fill="white")+
    geom_sf(data=data,aes(fill=fillplot,geometry=geometry),color="transparent")+ 
    geom_sf(fill=NA, color="grey") +
    # geom_sf(data=datasig_f,aes(size=siglvl),fill="grey",color="transparent",shape=21,alpha=0.1)+ 
    geom_sf(data=datasig_f,aes(size=siglvl,colour=sign),shape=8, alpha=0.4,lwd=.5)+ 
    scale_color_manual(values=c("blue" ="darkblue","red"="darkred"))+
    scale_size_manual(values=c(.1,.3),labels=c("p<=0.1","p<=0.05"), name="Trend significance level")+
    coord_sf(xlim = c(min(nco[,1]),max(nco[,1])), ylim = c(min(nco[,2]),max(nco[,2])))+
    # coord_sf(xlim = c(min(nco[,1]),max(nco[,1])), ylim = c(min(nco[,2]),max(nco[,2])))+
    # geom_point(data=event,aes(x=lon,y=lat,col=dr))+
    # scale_color_distiller(palette = "Spectral",
    #                       direction = -1, limit=c(0,max(eventp$dr)),oob = scales::squish)+
    scale_fill_gradientn(
      colors=palet,
      breaks=br,limits=limi,trans=trans,
      oob = scales::squish,na.value=colNA, name=legend)   +
    labs(x="Longitude", y = "Latitude")+
    guides(fill = guide_colourbar(barwidth = 1, barheight = 10),colour="none")+
    theme(axis.title=element_text(size=tsize),
          panel.background = element_rect(fill = "aliceblue", colour = "grey1"),
          panel.border = element_rect(linetype = "solid", fill = NA, colour="black"),
          legend.title = element_text(size=tsize),
          legend.text = element_text(size=osize),
          legend.position = "right",
          panel.grid.major = element_line(colour = "grey70"),
          panel.grid.minor = element_line(colour = "grey90"),
          legend.key = element_rect(fill = "transparent", colour = "transparent"),
          legend.key.size = unit(.8, "cm"))+
    ggtitle(title)
  
  ocrap
  
}

#Function displaying changes in both drought and floods
plotbichangemaps=function(basemap,catmap,datarD,datarF,law="GPD",type,period=c(1951,2020),valuenames){
  
  
  catmap$outlets=catmap$pointid
  data1=right_join(catmap,datarD,by = c("outlets"="unikout"))
  data2=right_join(catmap,datarF,by = c("outlets"="unikout"))
  datacol=names(data1)
  valcol=match(valuenames,datacol)
  datatwin1=data1
  datatwin2=data2
  valcol2=valcol-1
  st_geometry(datatwin2) <- NULL
  st_geometry(datatwin1) <- NULL
  datatwin1=as.data.frame(datatwin1)
  datatwin2=as.data.frame(datatwin2)
  class(datatwin1)
  dtc=names(datatwin1)
  mcor=unique(match(datarD$unikout,parlist$catchment))
  if (length(which(is.na(mcor)))>0) datar1=datar1[-which(is.na(mcor)),]
  
  cref=paste0("X",period[1])
  crefloc=match(cref,dtc)
  finalperiod=paste0("X",period[2])
  colsel=match(finalperiod,datacol)
  
  
  palet=c(hcl.colors(11, palette = "BrBG", alpha = NULL, rev = F, fixup = TRUE))
  title=paste0("Relative change in 100-years flood Return Level between ", period[1], " and ", period[2])
  legend="Relative Change (%)"
  datathin=datatwin1[valcol2]
  tmpval1=(datatwin1[,valcol2])/(datatwin1[,crefloc])*100-100
  tmpval2=(datatwin2[,valcol2])/(datatwin2[,crefloc])*100-100
  for (it in 1:length(tmpval1[1,])){
    tmpval1[,it][which(data1$IRES==1)]=datathin[which(data1$IRES==1),it]/datathin[which(data1$IRES==1),crefloc]
    tmpval1[which(data1$IRES==1),it][which(is.nan(tmpval1[which(data1$IRES==1),it]))]=0 #drying did not happen in both periods
    tmpval1[which(data1$IRES==1),it][which(is.infinite(tmpval1[which(data1$IRES==1),it]))]=-100 #drying only happens in recent period
    tmpval1[which(data1$IRES==1),it][which(tmpval1[which(data1$IRES==1),it]==0)]=100 #drying only happens in initial period
    tmpval1[,it][which(is.infinite(tmpval1[,it]))]=100
    tmpval1[,it][which(is.na(tmpval1[,it]))]=0
    
    tmpval1[,it][which(tmpval1[,it]>100)]=100
    tmpval1[,it][which(tmpval1[,it]<(-100))]=-100
    
    tmpval2[,it][which(tmpval2[,it]>100)]=100
    tmpval2[,it][which(tmpval2[,it]<(-100))]=-100
  }
  data1[,valcol]=tmpval1
  data2[,valcol]=tmpval2
  br=c(-100,-80,-60,-40,-20,0,20,40,60,80,100)
  limi=c(-100,100)
  trans=scales::modulus_trans(.8)
  colNA="transparent"
  
  st_geometry(data2) <- NULL
  data2=as.data.frame(data2)
  datamatch=match(data1$outlets,data2$outlets)
  names(data1)[colsel]="fillplot"
  
  data1$fillplot2=data2[datamatch,colsel-1]
  
  data1$fillplot=-data1$fillplot
  plot(data1$fillplot,data1$fillplot2)
  
  quantile(data1$fillplot2,.5,na.rm=T)
  databi <- bi_class(data1, x = fillplot, y = fillplot2, style = "quantile", dim = 3, keep_factors = TRUE, dig_lab=2)
  
  
  alterclass=data.frame(data1$fillplot)
  alterclass$class=NA
  alterclass$class[which(alterclass[,1]<(-10))]=1
  alterclass$class[which(alterclass[,1]>=-10 & alterclass[,1]<10)]=2
  alterclass$class[which(alterclass[,1]>10)]=3
  
  c1=alterclass$class
  alterclass=data.frame(data1$fillplot2)
  alterclass$class=NA
  alterclass$class[which(alterclass[,1]<(-10))]=1
  alterclass$class[which(alterclass[,1]>=-10 & alterclass[,1]<10)]=2
  alterclass$class[which(alterclass[,1]>10)]=3
  c2=alterclass$class
  
  cx=paste(c1,c2,sep="-")
  
  databi$bi_class=cx
  # create map
  map <- ggplot(basemap) +
    geom_sf(fill="white")+
    geom_sf(data = databi, mapping = aes(fill = bi_class), color = "transparent", size = 0.01, show.legend = FALSE) +
    coord_sf(xlim = c(min(nco[,1]),max(nco[,1])), ylim = c(min(nco[,2]),max(nco[,2])))+
    bi_scale_fill(pal = "BlueOr", dim = 3, na.value=colNA ) +
    labs()+
    theme(axis.title=element_text(size=tsize),
          panel.background = element_rect(fill = "aliceblue", colour = "grey1"),
          panel.border = element_rect(linetype = "solid", fill = NA, colour="black"),
          legend.title = element_text(size=tsize),
          legend.text = element_text(size=osize),
          legend.position = "bottom",
          panel.grid.major = element_line(colour = "grey70"),
          panel.grid.minor = element_line(colour = "grey90"),
          legend.key = element_rect(fill = "transparent", colour = "transparent"),
          legend.key.size = unit(.8, "cm"))
  
  map
  
  legend <- bi_legend(pal = "BlueOr",
                      dim = 3,
                      xlab = "  -  Drought changes  +  ",
                      ylab = "  -  Flood changes  +  ",
                      size = 16,
                      arrows = FALSE)
  
  
  
  pl=ggarrange(map, legend, 
               labels = c("Map", "Key"),
               ncol = 2, nrow = 1,widths = c(2,1), heights=c(1,1), vjust=-1)
  return(pl)
  
  
}

#Function displaying changes in both drought and floods
plotbichangemaps_qsp=function(basemap,catmap,datarD,datarF,UpArea,law="GPD",type,period=c(1951,2020),valuenames){
  
  
  catmap$outlets=catmap$pointid
  data1=right_join(catmap,datarD,by = c("outlets"="unikout"))
  data2=right_join(catmap,datarF,by = c("outlets"="unikout"))
  
  upag=match(data1$outlets,UpArea$outlets )
  data1$uparea=UpArea$upa[upag]
  upag=match(data2$outlets,UpArea$outlets )
  data2$uparea=UpArea$upa[upag]
  
  datacol=names(data1)
  valcol=match(valuenames,datacol)
  datatwin1=data1
  datatwin2=data2
  valcol2=valcol-1
  st_geometry(datatwin2) <- NULL
  st_geometry(datatwin1) <- NULL
  datatwin1=as.data.frame(datatwin1)
  datatwin2=as.data.frame(datatwin2)
  class(datatwin1)
  dtc=names(datatwin1)
  mcor=unique(match(datarD$unikout,parlist$catchment))
  if (length(which(is.na(mcor)))>0) datar1=datar1[-which(is.na(mcor)),]
  
  cref=paste0("X",period[1])
  crefloc=match(cref,dtc)
  finalperiod=paste0("X",period[2])
  colsel=match(finalperiod,datacol)
  
  
  palet=c(hcl.colors(11, palette = "BrBG", alpha = NULL, rev = F, fixup = TRUE))
  title=paste0("Relative change in 100-years flood Return Level between ", period[1], " and ", period[2])
  legend="Relative Change (%)"
  datathin=datatwin1[valcol2]
  
  tmpval1=(datatwin1[,valcol2]-datatwin1[,crefloc])
  tmpval2=(datatwin2[,valcol2]-datatwin2[,crefloc])
  
  tmpval1=tmpval1*1000/(data1$uparea)
  tmpval2=tmpval2*1000/(data2$uparea)
  # tmpval1=(datatwin1[,valcol2])/(datatwin1[,crefloc])*100-100
  # tmpval2=(datatwin2[,valcol2])/(datatwin2[,crefloc])*100-100
  # for (it in 1:length(tmpval1[1,])){
  #   tmpval1[,it][which(data1$IRES==1)]=datathin[which(data1$IRES==1),it]/datathin[which(data1$IRES==1),crefloc]
  #   tmpval1[which(data1$IRES==1),it][which(is.nan(tmpval1[which(data1$IRES==1),it]))]=0 #drying did not happen in both periods
  #   tmpval1[which(data1$IRES==1),it][which(is.infinite(tmpval1[which(data1$IRES==1),it]))]=-100 #drying only happens in recent period
  #   tmpval1[which(data1$IRES==1),it][which(tmpval1[which(data1$IRES==1),it]==0)]=100 #drying only happens in initial period
  #   tmpval1[,it][which(is.infinite(tmpval1[,it]))]=100
  #   tmpval1[,it][which(is.na(tmpval1[,it]))]=0
  #   
  #   tmpval1[,it][which(tmpval1[,it]>100)]=100
  #   tmpval1[,it][which(tmpval1[,it]<(-100))]=-100
  #   
  #   tmpval2[,it][which(tmpval2[,it]>100)]=100
  #   tmpval2[,it][which(tmpval2[,it]<(-100))]=-100
  # }
  data1[,valcol]=tmpval1
  data2[,valcol]=tmpval2
  br=c(-100,-80,-60,-40,-20,0,20,40,60,80,100)
  limi=c(-100,100)
  trans=scales::modulus_trans(.8)
  colNA="transparent"
  
  st_geometry(data2) <- NULL
  data2=as.data.frame(data2)
  datamatch=match(data1$outlets,data2$outlets)
  names(data1)[colsel]="fillplot"
  
  data1$fillplot2=data2[datamatch,colsel-1]
  
  data1$fillplot=-data1$fillplot
  data1$fillplot[which(data1$IRES==1)]=NA
  plot(data1$fillplot,data1$fillplot2, xlim=c(-2,2), ylim=c(-50,50))
  
  quantile(data1$fillplot2,.5,na.rm=T)
  databi <- bi_class(data1, x = fillplot, y = fillplot2, style = "quantile", dim = 3, keep_factors = TRUE, dig_lab=2)
  
  plot(data1$fillplot)
  quantile(data1$fillplot,c(0.75),na.rm=T)-quantile(data1$fillplot,c(0.25),na.rm=T)/5
  alterclass=data.frame(data1$fillplot)
  alterclass$class=NA
  alterclass$class[which(alterclass[,1]<(-0.10))]=1
  alterclass$class[which(alterclass[,1]>=-0.10 & alterclass[,1]<0.10)]=2
  alterclass$class[which(alterclass[,1]>0.10)]=3
  c1=alterclass$class
  
  
  hist(data1$fillplot2)
  quantile(data1$fillplot2,c(0.75),na.rm=T)-quantile(data1$fillplot2,c(0.25),na.rm=T)/5
  alterclass=data.frame(data1$fillplot2)
  alterclass$class=NA
  alterclass$class[which(alterclass[,1]<(-5))]=1
  alterclass$class[which(alterclass[,1]>=-5 & alterclass[,1]<5)]=2
  alterclass$class[which(alterclass[,1]>5)]=3
  c2=alterclass$class
  
  c1[which(is.na(c1))]=2
  cx=paste(c1,c2,sep="-")
  
  databi$bi_class=cx
  # create map
  map <- ggplot(basemap) +
    geom_sf(fill="white")+
    geom_sf(data = databi, mapping = aes(fill = bi_class), alpha=0.9, color = "transparent", size = 0.01, show.legend = FALSE) +
    geom_sf(fill=NA, color="gray42") +
    coord_sf(xlim = c(min(nco[,1]),max(nco[,1])), ylim = c(min(nco[,2]),max(nco[,2])))+
    bi_scale_fill(pal = "BlueOr", dim = 3, na.value=colNA ) +
    labs()+
    theme(axis.title=element_text(size=tsize),
          panel.background = element_rect(fill = "aliceblue", colour = "grey1"),
          panel.border = element_rect(linetype = "solid", fill = NA, colour="black"),
          legend.title = element_text(size=tsize),
          legend.text = element_text(size=osize),
          legend.position = "bottom",
          panel.grid.major = element_line(colour = "grey70"),
          panel.grid.minor = element_line(colour = "grey90"),
          legend.key = element_rect(fill = "transparent", colour = "transparent"),
          legend.key.size = unit(.8, "cm"))
  
  map
  
  legend <- bi_legend(pal = "BlueOr",
                      dim = 3,
                      xlab = "  -  Drought changes  +  ",
                      ylab = "  -  Flood changes  +  ",
                      size = 16,
                      arrows = FALSE)
  
  #reproduce my color palette: 
  loscolors=c("#174f28","#845e29","#dd6a29","#167984","#819185","#d8a386","#169dd0","#7ebbd2","#d3d3d3")
  al=0.8
  dataX=data1[-which(is.na(data1$fillplot)),]
  dataX$density <- get_density(dataX$fillplot, dataX$fillplot2, n = 100)
  monstre_legend<-ggplot() +
    annotate("rect", xmin=-Inf, xmax=-0.1, ymin=-Inf, ymax=-5, fill=loscolors[9], alpha=al) +
    annotate("rect", xmin=-Inf, xmax=-0.1, ymin=-5, ymax=5, fill=loscolors[8], alpha=al) +
    annotate("rect", xmin=-Inf, xmax=-0.1, ymin=5, ymax=Inf, fill=loscolors[7], alpha=al) +
    
    annotate("rect", xmin=-0.1, xmax=0.1, ymin=-Inf, ymax=-5, fill=loscolors[6], alpha=al) +
    annotate("rect", xmin=-0.1, xmax=0.1, ymin=-5, ymax=5, fill=loscolors[5], alpha=al) +
    annotate("rect", xmin=-0.1, xmax=0.1, ymin=5, ymax=Inf, fill=loscolors[4], alpha=al) +
    
    annotate("rect", xmin=0.1, xmax=Inf, ymin=-Inf, ymax=-5, fill=loscolors[3], alpha=al) +
    annotate("rect", xmin=0.1, xmax=Inf, ymin=-5, ymax=5, fill=loscolors[2], alpha=al) +
    annotate("rect", xmin=0.1, xmax=Inf, ymin=5, ymax=Inf, fill=loscolors[1], alpha=al) +
    
    geom_point(data=dataX, aes(x=fillplot, y=fillplot2, col=density),alpha=0.6,size=2, shape=16) +
    #geom_point(data= Rsig, aes(x=ObsChange, y=SimChange),fill="transparent", color="gray25",shape=21,size=2)+
    #geom_abline(slope=1,intercept=0,col="gray25",lwd=1.5,alpha=.5)+
    # 
    # annotate("label", x=-90, y=100, label= paste0("N = ",ln2),size=5)+
    # annotate("label", x=90, y=100, label= paste0("N = ",lp1),size=5)+
    # annotate("label", x=90, y=-100, label= paste0("N = ",lp2),size=5)+
    # annotate("label", x=-90, y=-100, label= paste0("N = ",ln1),size=5)+
    scale_color_viridis(option="F")+
    scale_x_continuous(name="10-y Drought RL change (l/s/km2)",breaks = seq(-5,5,by=1), labels=seq(5,-5,by=-1),limits = c(-2.5,2.5))+
    scale_y_continuous(name="10-y Flood RL change (l/s/km2)",breaks = seq(-100,100,by=50), limits = c(-100,100))+
    # scale_color_gradientn(
    #   colors=palet, limits=c(-1,1),oob = scales::squish,
    #   name="temporal correlation",breaks=seq(-1,1, by=0.2))+
    # scale_alpha_continuous(name="trend significance",range = c(0.5, 1))+
    theme(axis.title=element_text(size=tsize),
          panel.background = element_rect(fill = "white", colour = "grey1"),
          panel.border = element_rect(linetype = "solid", fill = NA, colour="black"),
          legend.title = element_text(size=tsize),
          legend.text = element_text(size=osize),
          legend.position = "none",
          panel.grid.major = element_line(colour = "grey70"),
          panel.grid.minor = element_line(colour = "grey90"),
          legend.key = element_rect(fill = "transparent", colour = "transparent"),
          legend.key.size = unit(.8, "cm"))
  
  
  
  
  pl=ggarrange(map, 
               ggarrange(legend, monstre_legend, 
                         ncol = 1, nrow = 2),
               ncol=2, nrow=1, widths = c(3,2), heights=c(1,1),vjust=-1)
  
  pl<-annotate_figure(pl, top = text_grob(paste0("Joint changes in floods and droughts (",period[1],"-",period[2],")"), 
                                          color = "black", face = "bold", size = 14))
  pl
  return(pl)
  
  
}

#Main function. Diplay changes in river flood at pixel and catchment level
plotchangemapix=function(basemap,catmap,datar,law="GPD",type,period=c(1950,2020),parlist,hybasf,haz="flood"){
  
  
  #data=right_join(outloc,datar,by = c("outlets"="unikout"))
  data=datar
  
  data$unikout=data$outlets
  #st_geometry(data) <- NULL
  datacol=names(data)
  years=c(period[1]:period[2])
  valuenames=paste0("Y",years)
  
  valcol=match(valuenames,datacol)
  datatwin=data
  valcol2=valcol
  datacl=data[c(valcol,which(names(data)=="outlets"))]
  
  datatwin=as.data.frame(datatwin)
  class(datatwin)
  dtc=names(datatwin)
  mcor=unique(match(datar$unikout,parlist$catchment))
  
  if (length(which(is.na(mcor)))>0) datar=datar[-which(is.na(mcor)),]
  
  cref=paste0("Y",period[1])
  crefloc=match(cref,dtc)
  finalperiod=paste0("Y",period[2])
  colsel=match(finalperiod,datacol)
  
  if (type=="RLchange"){
    palet=c(hcl.colors(11, palette = "BrBG", alpha = NULL, rev = F, fixup = TRUE))
    title=paste0("Relative change in 100-years flood Return Level between ", period[1], " and ", period[2])
    legend="Relative Change (%)"
    tmpval=(datatwin[,valcol2])/(datatwin[,crefloc])*100-100
    for (it in 1:length(tmpval[1,])){
      tmpval[,it][which(tmpval[,it]>800)]=NA
    }
    data[,valcol]=tmpval
    br=c(-100,-80,-50,-20,0,20,50,80,100)
    if(haz=="drought"){
      limi=c(-100,100)
    }else if (haz=="flood"){
      limi=c(-50,50)
    }
    trans=scales::modulus_trans(.8)
    colNA="gray"
  }
  if (type=="RPchange"){
    palet=c(hcl.colors(9, palette = "BrBG", alpha = NULL, rev = T, fixup = TRUE))
    title=paste0(period[2]," return period of a 100-years \n Return Level ",haz," in ", period[1])
    legend="Return Period (Years)"
    br=c(10,20,50,100,200,500,1000)
    limi=c(10,1000)
    trans="log"
    if(haz=="drought"){
      vc= match(valuenames,names(datacl))
      datacl[,vc]=1/datacl[,vc]
    }
    newRP=RPchangeCal(parlist, yi=period[1], yf=period[2], RetLev=datacl,law, valuenames)
    ow=match(data$outlets,newRP$catchment)
    data[,colsel]=newRP$newRP[ow]
    colNA="black"
    
  }
  
  
  #there is something to be fixed here
  names(data)[colsel]="fillplot"
  
  #match here with catchments
  st_geometry(catmap) <- NULL
  ohoho=catmap[,c(1:21)]
  #need to check what happens here (clean catcmentproj)
  
  data2=inner_join(data,catmap[,c(1:20)],by=c("latlong"="llcoord"))
  points <- st_as_sf(data2[,c(1:80)], coords = c("Var1", "Var2"), crs = 4326)
  points <- st_transform(points, crs = 3035)
  
  rev=FALSE
  if(type=="RPchange") rev=TRUE
  #palet=c(hcl.colors(11, palette = "RdYlBu", alpha = NULL, rev = rev, fixup = TRUE))
  pl1=ggplot(basemap) +
    geom_sf(fill="gray85")+
    geom_sf(fill=NA, color="grey") +
    geom_sf(data=points,aes(col=fillplot,geometry=geometry),alpha=1,size=0.25,stroke=0,shape=15)+ 
    coord_sf(xlim = c(min(nco[,1]),max(nco[,1])), ylim = c(min(nco[,2]),max(nco[,2])))+
    # coord_sf(xlim = c(min(nco[,1]),max(nco[,1])), ylim = c(min(nco[,2]),max(nco[,2])))+
    # geom_point(data=event,aes(x=lon,y=lat,col=dr))+
    # scale_color_distiller(palette = "Spectral",
    #                       direction = -1, limit=c(0,max(eventp$dr)),oob = scales::squish)+
    scale_color_gradientn(
      colors=palet,
      breaks=br,limits=limi,trans=trans,
      oob = scales::squish,na.value=colNA, name=legend)   +
    labs(x="Longitude", y = "Latitude")+
    guides(color = guide_colourbar(barwidth = 12, barheight = 1))+
    theme(axis.title=element_text(size=tsize),
          panel.background = element_rect(fill = "aliceblue", colour = "grey1"),
          panel.border = element_rect(linetype = "solid", fill = NA, colour="black"),
          legend.title = element_text(size=tsize),
          legend.text = element_text(size=osize),
          legend.position = "bottom",
          panel.grid.major = element_line(colour = "grey70"),
          panel.grid.minor = element_line(colour = "grey90"),
          legend.key = element_rect(fill = "transparent", colour = "transparent"),
          legend.key.size = unit(.8, "cm"))+
    ggtitle(title)
  
  #aggregate by hybas07 catchment
  
  if (type=="RLchange"){
    pointagg=aggregate(list(Rchange=points$fillplot),
                       by = list(HYBAS_ID=points$HYBAS_ID),
                       FUN = function(x) c(mean=median(x,na.rm=T),dev=sd(x,na.rm=T),len=length(x)))
    pointagg <- do.call(data.frame, pointagg)
  }
  if (type=="RPchange"){
    pointagg=aggregate(list(Rchange=points$fillplot),
                       by = list(HYBAS_ID=points$HYBAS_ID),
                       FUN = function(x) c(mean=median(x,na.rm=T),dev=sd(x,na.rm=T),len=length(x),q1=quantile(x, 0.25, na.rm=T),q3=quantile(x, 0.75, na.rm=T)))
    pointagg <- do.call(data.frame, pointagg)
  }
  
  #natch pointagg with hybasf
  pointplot=inner_join(hybasf,pointagg,by= "HYBAS_ID")
  
  pl2=ggplot(basemap) +
    geom_sf(fill="gray95")+
    geom_sf(data=pointplot,aes(fill=Rchange.mean,geometry=geometry),color="transparent")+ 
    geom_sf(fill=NA, color="grey") +
    coord_sf(xlim = c(min(nco[,1]),max(nco[,1])), ylim = c(min(nco[,2]),max(nco[,2])))+
    # coord_sf(xlim = c(min(nco[,1]),max(nco[,1])), ylim = c(min(nco[,2]),max(nco[,2])))+
    # geom_point(data=event,aes(x=lon,y=lat,col=dr))+
    # scale_color_distiller(palette = "Spectral",
    #                       direction = -1, limit=c(0,max(eventp$dr)),oob = scales::squish)+
    scale_fill_gradientn(
      colors=palet,
      breaks=br,limits=limi,trans=trans,
      oob = scales::squish,na.value=colNA, name=legend)   +
    labs(x="Longitude", y = "Latitude")+
    guides(fill = guide_colourbar(barwidth = 15, barheight = .8))+
    theme(axis.title=element_text(size=tsize),
          panel.background = element_rect(fill = "aliceblue", colour = "grey1"),
          panel.border = element_rect(linetype = "solid", fill = NA, colour="black"),
          legend.title = element_text(size=tsize),
          legend.text = element_text(size=osize),
          legend.position = "bottom",
          panel.grid.major = element_line(colour = "grey70"),
          panel.grid.minor = element_line(colour = "grey90"),
          legend.key = element_rect(fill = "transparent", colour = "transparent"),
          legend.key.size = unit(.8, "cm"))+
    ggtitle(title)
  
  return(list(pl1,pl2,points,pointplot))
  
  
}

#Function to plot changes in RL along with statistical significance of the trend. do significance at catchment level by averaging.
plotTrendSipix=function(basemap,datar,period=c(1950,2020),hybasf,valuenames,nco){
  
  
  #data=right_join(outloc,datar,by = c("outlets"="unikout"))
  data=datar
  data$unikout=data$outlets
  data$Y2020=data$Y2019
  st_geometry(data) <- NULL
  
  
  #mcor=unique(match(data$unikout,parlist$catchment))
  
  years=c(period[1]:period[2])
  valuenames=paste0("Y",years)
  datacol=names(data)
  
  
  #this need to be changed to make advantage of the whole data
  valcol=match(valuenames,datacol)
  datatwin=data
  valcol2=valcol
  #st_geometry(datatwin) <- NULL
  datatwin=as.data.frame(datatwin)
  class(datatwin)
  dtc=names(datatwin)
  #dtreg=datatwin[,c(41:48)]
  # miniTS=ts(data=(as.numeric(dtreg[2,])),start=year(as.Date(paste0(Impdates[1],"-06-01"))),
  #           end=year(as.Date(paste0(Impdates[8],"-06-01"))),frequency = .1)
  # minReg=as.data.frame(miniTS)
  # minReg$time=as.numeric((Impdates))
  # 
  # mk=MannKendall(miniTS)
  # mkt=mk$tau
  # mks=mk$sl
  # plot(miniTS)
  # optn=lm(x~time, data=minReg)
  # abline(optn)
  # mcor=unique(match(data$unikout,parlist$catchment))
  # 
  # if (length(which(is.na(mcor)))>0) datar=datar[-which(is.na(mcor)),]
  cref=paste0("Y",period[1])
  crefloc=match(cref,dtc)
  finalperiod=paste0("Y",period[2])
  colsel=match(finalperiod,datacol)
  tmpval=(datatwin[,valcol2])/(datatwin[,crefloc])*100-100
  data[,valcol]=tmpval
  palet=c(hcl.colors(11, palette = "RdBu", alpha = NULL, rev = F, fixup = TRUE))
  title=paste0("Relative change in 100-years flood Return Level (", period[1], "-", period[2],")")
  
  #aggregate by hybas07 catchment
  
  pointagg=aggregate(list(Rchange=data[,valcol]),
                     by = list(HYBAS_ID=data$HYBAS_ID),
                     FUN = function(x) c(mean=mean(x)))
  pointagg <- do.call(data.frame, pointagg)
  
  
  #legend="Kendall tau"
  legend="Relative change (%)"
  tmpval=(datatwin[,valcol2])
  mkta=c()
  mksa=c()
  for (it in 1:length(pointagg[,1])){
    if (it%%1000==0) print(it)
    mks=NA
    mkt=NA
    miniTS=as.numeric(pointagg[it,])
    print(it)
    if (!is.na(miniTS[2])){
      #mk=MannKendall(miniTS[-1])
      mk2=mmkh(miniTS[-1],ci=0.95)
      mk=data.frame(tau=mk2[6],sl=mk2[2])
      
      #compute trend as well with sen.slope
      mkt=mk$tau
      mks=mk$sl
    }else{
      mkt=NA
      mks=NA
    }
    mkta=c(mkta,mkt)
    mksa=c(mksa,mks)
  }
  # br=c(-1,-.80,-.60,-.40,-.20,0,.20,.40,.60,.80,1)
  # limi=c(-1,1)
  trans=scales::modulus_trans(.6)
  colNA="darkgrey"
  #title="Trends in 100 years Re flood magnitude"
  
  pointagg$mkta=mkta
  pointagg$sl=mksa
  
  
  
  #give a projection to pixels
  names(data)[colsel]="fillplot"
  points <- st_as_sf(data, coords = c("Var1", "Var2"), crs = 4326)
  points <- st_transform(points, crs = 3035)
  
  
  #natch pointagg with hybasf
  pointsag=inner_join(hybasf,pointagg,by= "HYBAS_ID")
  
  pointsag$siglvl=0
  pointsag$siglvl[which(pointsag$sl<=0.1)]=1
  
  catsig=pointsag[which(pointsag$siglvl>0),]
  catsig <- st_transform(catsig, crs = 3035)
  ## create a grid of points
  grdpts <- sf::st_make_grid(catsig, what = "centers",cellsize = 40000)
  my.points <- sf::st_sf(grdpts)
  
  sf_use_s2(FALSE)
  pointsInside <- sf::st_join(x = my.points, y = catsig, left = FALSE)
  pointsInside$sign="positive"
  pointsInside$sign[which(pointsInside$mkta<0)]="negative"
  pointsag$siglvl=factor(pointsag$siglvl)
  
  # #TBC with the river pixels and better colors
  #   ggplot(basemap) +
  #   geom_sf(data=pointsInside,aes(geometry=geometry, fill=sign),alpha=.7,size=2,stroke=0,shape=21,color="black")+
  #   coord_sf(xlim = c(min(nco[,1]),max(nco[,1])), ylim = c(min(nco[,2]),max(nco[,2])))+
  #   scale_fill_manual(values=c("tomato4","steelblue4"), name="Significant trends")
  #   
  # ggplot(regiod)+
  # geom_sf(mapping=aes(geometry=geometry,group=name,fill=name))
  br=c(-100,-80,-50,-25,-10,0,10,25,50,80,100)
  limi=c(-55,55)
  tsize=16
  osize=12
  ocrap<-ggplot(basemap) +
    geom_sf(fill="gray85",color="darkgrey",size=0.5)+
    geom_sf(data=points,aes(col=fillplot,geometry=geometry),alpha=1,size=0.25,stroke=0,shape=15)+ 
    geom_sf(data=pointsInside,aes(geometry=geometry, fill=sign),alpha=.7,size=.9,stroke=0,shape=21,color="black")+
    #geom_sf(data=pointsag,aes(fill=siglvl,geometry=geometry),alpha=0.4,color="transparent")+ 
    geom_sf(fill=NA, color="grey") +
    #geom_sf(regiod,mapping=aes(geometry=geometry,group=name),fill=NA,color="orange")+
    coord_sf(xlim = c(min(nco[,1]),max(nco[,1])), ylim = c(min(nco[,2]),max(nco[,2])))+
    # geom_sf(data=datasl1,aes(col=fillplot,geometry=geometry),alpha=0.8,size=0.05,stroke=0,shape=1)+ 
    # geom_sf(data=datasl2,aes(col=fillplot,geometry=geometry),alpha=0.9,size=0.1,stroke=0,shape=1)+ 
    #geom_sf(data=data,aes(fill=fillplot,geometry=geometry),color="transparent")+ 
    #geom_sf(data=datasig_f,aes(size=siglvl),fill="grey",color="transparent",shape=21,alpha=0.1)+ 
    #scale_color_manual(values=c("blue" ="darkblue","red"="darkred"))+
    #geom_col_pattern(data=pointsag,aes(fill=siglvl,geometry=geometry),colour='black', pattern = 'circle') +
    scale_fill_manual(values=c("tomato4","steelblue4"), name="Significant trends")+
    ggtitle(title)+
    # coord_sf(xlim = c(min(nco[,1]),max(nco[,1])), ylim = c(min(nco[,2]),max(nco[,2])))+
    # geom_point(data=event,aes(x=lon,y=lat,col=dr))+
    # scale_color_distiller(palette = "Spectral",
    #                       direction = -1, limit=c(0,max(eventp$dr)),oob = scales::squish)+
    scale_color_gradientn(
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
  
}


#Ridgeplots of changes over a given area
plotHistoDates=function(outloc,datar,law="GPD",type,period=c(1951,2020),parlist,valuenames){
  
  
  #data=right_join(outloc,datar,by = c("outlets"="unikout"))
  data=datar
  data$unikout=data$outlets
  datacol=names(data)
  valcol=match(valuenames,datacol)
  datatwin=data
  valcol2=valcol
  
  datatwin=as.data.frame(datatwin)
  class(datatwin)
  dtc=names(datatwin)
  mcor=unique(match(datar$unikout,parlist$catchment))
  
  if (length(which(is.na(mcor)))>0) datar=datar[-which(is.na(mcor)),]
  
  cref=paste0("Y",period[1])
  crefloc=match(cref,dtc)
  # finalperiod=paste0("Y",period[2])
  # colsel=match(finalperiod,datacol)
  
  palet=c(hcl.colors(11, palette = "BrBG", alpha = NULL, rev = F, fixup = TRUE))
  title=paste0("Relative change in 100-years flood Return Level between ", period[1], " and ", period[2])
  legend="Relative Change (%)"
  tmpval=(datatwin[,valcol2])/(datatwin[,crefloc])*100-100
  for (it in 1:length(tmpval[1,])){
    tmpval[,it][which(tmpval[,it]>800)]=NA
  }
  data[,valcol]=tmpval
  br=c(-100,-80,-60,-40,-20,0,20,40,60,80,100)
  limi=c(-100,100)
  trans=scales::modulus_trans(.8)
  colNA="transparent"
  
  ret=list()
  
  subdata=inner_join(outloc,data,by=c("outlets"="unikout"))
  
  Regname=unique(subdata$name)
  if (is.null(Regname)) Regname = "Europe"
  valcol3=match(valuenames,colnames(subdata))
  crefloc2=match(cref,colnames(subdata))
  subdata[,crefloc2]=NA
  nref=colnames(subdata)[crefloc2]
  ptref=data.frame(x=0,nref)
  dfd=subdata[,valcol3]
  
  dfplot= dfd %>% gather("Decade","Val")
  
  stats=aggregate(list(val=dfplot$Val),
                  by = list(Decades=dfplot$Decade),
                  FUN = function(x) c(mean=mean(x,na.rm=T),dev=sd(x,na.rm=T),q25=quantile(x,0.25,na.rm=T),q75=quantile(x,0.75,na.rm=T),median=quantile(x,0.5,na.rm=T)))
  stats <- do.call(data.frame, stats)
  stats=data.frame(rep(Regname,length(stats$Decades)),stats)
  names(stats)[1]="Region"
  
  ret[[1]]=stats
  # Plot
  ret[[2]]<-ggplot(dfplot, aes(x = Val, y = Decade, fill = stat(x)),alpha=.1) +
    geom_density_ridges_gradient(alpha=0.4,scale = 2, rel_min_height = 0.002, quantile_lines = TRUE,
                                 size=1,bandwidth=2,  gradient_lwd = 1.) +
    scale_fill_gradientn( name= "100 years flood RL \n change compared to 1950  ", colors=palet,limits=c(-30,30),
                          oob = scales::squish) +
    geom_point(data=ptref,aes(x=x, y = nref), col="black",size=7,pch=21)+
    scale_x_continuous(expand = c(0, 0),limits=c(-100,100),name= "Change in Return Level",breaks = seq(-100,100,by=20)) +
    scale_y_discrete(expand = c(0.1, 0.1),labels = c("1950","1960","1970","1980","1990","2000","2010","2020"), name="Years") +
    theme_ridges() +
    theme(
      panel.spacing = unit(0.3, "lines"),
      strip.text.x = element_text(size = 8),
      legend.position = "bottom",
      legend.key.width= unit(2.5, 'cm'),
      legend.text =element_text(size=14),
      panel.border = element_rect(linetype = "solid", fill = NA, colour="black"),
      legend.title = element_text(size=14),
      plot.title = element_text(size=20)
    )+
    ggtitle(Regname)
  
  return(ret)
}

#Function to compute changes in return periods associated to a return level.
RPchangeCal=function(parlist, yi, yf, RetLev,law,valuenames){
  #Use TSEVA function for return level computation
  #1 selection of the year I want
  #check if the computation can be done on all the points at the same time
  #ay ay ay
  cname=names(RetLev)
  cna=match(valuenames,cname)
  parlistinit=parlist[which(parlist$Year==yi),]
  parlistf=parlist[which(parlist$Year==yf),]
  cref=paste0("Y",yi)
  crefloc=match(cref,cname)
  RLi=RetLev[,crefloc]
  
  
  if (law=="GEV"){
    paramx=data.frame((parlistf))
    qxV=1-exp(-(1+paramx$epsilonGEV*(RLi-paramx$muGEV)/paramx$sigmaGEV)^(-1/paramx$epsilonGEV))
    returnPeriods=1/qxV
  } else if (law=="GPD"){
    paramx=data.frame((parlistf))
    X0 <- paramx$nPeaks/paramx$SampleTimeHorizon
    qxD=(((1+paramx$epsilonGPD*(RLi-paramx$thresholdGPD)/paramx$sigmaGPD)^(-1/paramx$epsilonGPD)))
    returnPeriods=1/(X0*qxD)
  }
  min(returnPeriods)
  
  #jpp=inner_join(paramx$catchment,returnPeriods)
  return(data.frame(catchment=paramx$catchment,newRP=returnPeriods))
}

#Function to return a plot of the temporal distribution of changes in RL
#fix this

FascPlotChange<- function(datar,outloc,UpArea,period=c(1951,2020),lims=c(-1, 1),q1,q2){
  
  seqy=c(period[1]:period[2])
  data=datar
  subdata=inner_join(outloc,data,by=c("outlets"="unikout"))
  upag=match(subdata$outlets,UpArea$outlets )
  subdata$uparea=UpArea$upa[upag]
  
  datacol=names(subdata)
  valuenames=paste0("X",seq(period[1],period[2]))
  valcol=match(valuenames,datacol)
  oldata=subdata
  subdata=subdata[,valcol]
  
  #datrd=subdata/subdata[,1]*100-100
  datrd=subdata-subdata[,1]
  datrd=datrd*1000/oldata$uparea
  #smooth
  datrp=data.frame(unikout=datar$unikout,datrd)
  lwtf=c()
  yvf=c()
  for (y in 1:length(seqy)){
    print(y)
    yx=seqy[y]
    daty=datrp[,which(names(datrp)==paste0("X",yx))]
    blss=which(daty>1000)
    datrd[blss,which(names(datrd)==paste0("",yx))]=NA
    rmy=datrp$unikout[blss]
    
    #compare with list already existing
    mex=na.omit(match(lwtf,rmy))
    if (length(mex)>0) rmy=rmy[-mex]
    yv=rep(yx,length(rmy))
    yvf=c(yvf,yv)
    lwtf=c(lwtf,rmy)
  }
  
  outrm=data.frame(yvf,lwtf)
  
  dfplot= datrd %>% gather("Year","pixel")
  dfplot$year=as.numeric(substr(dfplot$Year,2,5))
  if (length(which(is.na(dfplot$pixel)))>0) dfplot=dfplot[-which(is.na(dfplot$pixel)),]
  if (length(which(is.infinite(dfplot$pixel)))>0) dfplot=dfplot[-which(is.infinite(dfplot$pixel)),]
  
  countfunction<-function(values){
    l1=length(which(values<0.3 & values>0))
    l2=length(which(values>0.3))
    l3=length(which(values<(-0.3)))
    l4=length(which(values>-0.3 & values<0))
    return(c(l2,l1,-l4,-l3))
  }
  dfq=aggregate(list(change=dfplot$pixel),
                by = list(Year=dfplot$year),
                FUN = function(x) c(bg=countfunction(x)))
  dfq <- do.call(data.frame, dfq)
  
  data_long <- gather(dfq, group, value, change.bg1:change.bg4) %>%
    arrange(factor(Year)) %>% 
    mutate(x=factor(Year, levels=unique(Year)))
  
  
  colorz = c('change.bg4'='#d73027','change.bg3'='#fee090','change.bg2'='lightblue','change.bg1'='royalblue')
  labs=rev(c("< -0.1","-0.1 - 0","0 - 0.1", ">0.1"))
  legend_order=c('change.bg1', 'change.bg2', 'change.bg3','change.bg4')
  data_long$p2 <- relevel(data_long$position, 'change.bg1')
  data_long$position <- factor(data_long$group, levels=c('change.bg1', 'change.bg2', 'change.bg4','change.bg3'))
  
  ggplot(data_long, aes(fill=position, y=value, x=Year)) + 
    geom_bar(position="stack", stat="identity")+
    scale_fill_manual(values = colorz, breaks = legend_order,name="change in qsp",labels=labs)  + 
    scale_y_continuous(name="Number of catchments")+
    scale_x_continuous(name="Years") +
    # scale_x_discrete(labels= c("< -0.41","-0.41-0.2", "0.2-0.5","0.5-0.7","0.7-0.8", ">0.8"), name="KGE")+
    theme(axis.title=element_text(size=16, face="bold"),
          axis.text = element_text(size=13),
          panel.background = element_rect(fill = "white", colour = "white"),
          panel.grid = element_blank(),
          panel.border = element_rect(linetype = "solid", fill = NA, colour="black"),
          legend.title = element_text(size=14),
          legend.text = element_text(size=12),
          legend.position = "bottom",
          panel.grid.major = element_line(colour = "grey80"),
          panel.grid.minor.x = element_line(colour = "grey90",linetype="dashed"),
          legend.key = element_rect(fill = "transparent", colour = "transparent"),
          legend.key.size = unit(.8, "cm"))
  
  
  dfquantiles=aggregate(list(change=dfplot$pixel),
                        by = list(Year=dfplot$year),
                        FUN = function(x) c(mean=mean(x,na.rm=T),q1=quantile(x, q1, na.rm=T),q3=quantile(x, q2, na.rm=T)))
  
  dfquantiles <- do.call(data.frame, dfquantiles)
  names(dfquantiles)=c("Year", "mean","qlow","qhigh")
  
  
  
  dtfinal=c()
  xmm=c(-10,10)
  for( y in 1:length(seqy)){
    print(y)
    yx=seqy[y]
    dfp1=dfplot[which(dfplot$year==yx),]
    dfp1=dfp1[which(!is.na(dfp1$pixel)),]
    merde=kde(dfp1$pixel,bgridsize=1500,xmin=xmm[1],xmax=xmm[2],h=0.1)
    
    dfqy=dfquantiles[which(dfquantiles$Year==yx),]
    ziz=data.frame(ev.points=merde$eval.points,estimates=merde$estimate)
    ziz$estimates=ziz$estimates/sum(ziz$estimates)
    sum(ziz$estimates)
    dtz=cbind(rep(yx,1500),ziz)
    dtz[,3][which(dtz[,2]<dfqy$qlow)]=NA
    dtz[,3][which(dtz[,2]>dfqy$qhigh)]=NA
    dtfinal=rbind(dtfinal,dtz)
  }
  
  dfqplus=seq(period[1],period[2],by=0.1)
  dfql=interp1(dfquantiles$Year,dfquantiles$qlow,dfqplus)
  dfqh=interp1(dfquantiles$Year,dfquantiles$qhigh,dfqplus)
  dfqplus=data.frame(dfqplus,dfql,dfqh)
  names(dfqplus)=c("Year","qlow","qhigh")
  
  #alpha=1+(log(values)))
  
  names(dtfinal)=c("Year","Change","values")
  ptref=data.frame(x=1950,y=0)
  
  palet=c(hcl.colors(11, palette = "YlGnBu", alpha = NULL, rev = T, fixup = TRUE))
  lims=c(-10,10)
  plo=ggplot(dtfinal, aes(x = Year, y = Change))+
    geom_raster(aes(fill=values,alpha=(log(values))),interpolate = T)+
    scale_fill_gradientn(colours = palet,
                         n.breaks=10,guide='coloursteps',na.value = "transparent",oob = scales::squish,limits=c(0,0.005))+
    scale_y_continuous(
      n.breaks = 10,name= "Change in Return Level",
      limit=lims,expand = c(0, 0))+
    scale_x_continuous(
      n.breaks = 10,name= "Years",
      limit=c(1949, 2020),expand = c(0, 0))+
    geom_line(data=dfquantiles,aes(x=Year,y=mean),size=2,col="purple",alpha=0.9)+
    geom_line(data=dfquantiles,aes(x=Year,y=qlow),size=1,col="grey",linetype="longdash")+
    geom_line(data=dfquantiles,aes(x=Year,y=qhigh),size=1,col="grey",linetype="longdash")+
    geom_hline(yintercept=0)+
    geom_point(data=ptref,aes(x=x, y = y), col="black",fill="white",size=7,pch=21)+
    theme_bw() +
    theme(
      panel.spacing = unit(0.3, "lines"),
      strip.text.x = element_text(size = 8),
      legend.position = "none",
      legend.key.width= unit(2.5, 'cm'),
      legend.text =element_text(size=14),
      panel.border = element_rect(linetype = "solid", fill = NA, colour="black"),
      legend.title = element_text(size=14),
      plot.title = element_text(size=20)
    )+
    ggtitle(unique(outloc$name))
  
  plo
  return(plo)
}

FascPlotBiChange<- function(datarD, datarF,outloc,UpArea,period=c(1951,2020),lims=c(-1, 1),q1,q2){
  
  catmap$outlets=catmap$pointid
  data1=right_join(outloc,datarD,by = c("outlets"="unikout"))
  data2=right_join(outloc,datarF,by = c("outlets"="unikout"))
  
  upag=match(data1$outlets,UpArea$outlets )
  data1$uparea=UpArea$upa[upag]
  upag=match(data2$outlets,UpArea$outlets )
  data2$uparea=UpArea$upa[upag]
  
  datacol=names(data1)
  valuenames=paste0("X",seq(period[1],period[2]))
  valcol=match(valuenames,datacol)
  
  datatwin1=data1
  datatwin2=data2
  valcol2=valcol-1
  st_geometry(datatwin2) <- NULL
  st_geometry(datatwin1) <- NULL
  datatwin1=as.data.frame(datatwin1)
  datatwin2=as.data.frame(datatwin2)
  class(datatwin1)
  dtc=names(datatwin1)
  mcor=unique(match(datarD$unikout,parlist$catchment))
  if (length(which(is.na(mcor)))>0) datar1=datar1[-which(is.na(mcor)),]
  
  cref=paste0("X",period[1])
  crefloc=match(cref,dtc)
  finalperiod=paste0("X",period[2])
  colsel=match(finalperiod,datacol)
  
  seqy=c(period[1]:period[2])
  
  oldata=data1
  subdata=data1[,valcol]
  subdata2=data2[,valcol]
  #datrd=subdata/subdata[,1]*100-100
  datrd=subdata-subdata[,1]
  #datrd=subdata
  datrd=datrd*1000/oldata$uparea
  
  datrf=subdata2-subdata2[,1]
  #datrf=subdata2
  datrf=datrf*1000/data2$uparea
  #smooth
  datrp=data.frame(unikout=data1$outlets,datrd)
  datrq=data.frame(unikout=data2$outlets,datrf)
  lwtf=c()
  yvf=c()
  for (y in 1:length(seqy)){
    print(y)
    yx=seqy[y]
    daty=datrp[,which(names(datrp)==paste0("X",yx))]
    datz=datrq[,which(names(datrq)==paste0("X",yx))]
    blss=which(daty>1000)
    datrd[blss,which(names(datrd)==paste0("",yx))]=NA
    rmy=datrp$unikout[blss]
    
    #compare with list already existing
    mex=na.omit(match(lwtf,rmy))
    if (length(mex)>0) rmy=rmy[-mex]
    yv=rep(yx,length(rmy))
    yvf=c(yvf,yv)
    lwtf=c(lwtf,rmy)
  }
  
  outrm=data.frame(yvf,lwtf)
  
  dfplot= datrd %>% gather("Year","pixel")
  dfplot$year=as.numeric(substr(dfplot$Year,2,5))
  
  dfplot2=datrf %>% gather("Year","pixel")
  dfplot2$year=as.numeric(substr(dfplot2$Year,2,5))
  
  dfplot$haz="low"
  dfplot2$haz="high"
  dfplot2$pixel=dfplot2$pixel/10
  dfplotH=rbind(dfplot,dfplot2)
  Dcad=seq(1955,2015,by =10)
  mb=which(!is.na(match(dfplotH$year,Dcad)))
  
  dfplotX=dfplotH[mb,]
  
  q1=0.25
  q2=0.75
  #plot with curve and confidence interval shaded colors
  dfquantiles=aggregate(list(change=dfplotX$pixel),
                        by = list(Year=dfplotX$year,haz=dfplotX$haz),
                        FUN = function(x) c(median=median(x,na.rm=T),q1=quantile(x, q1, na.rm=T),q3=quantile(x, q2, na.rm=T)))
  
  dfquantiles <- do.call(data.frame, dfquantiles)
  names(dfquantiles)=c("Year","haz", "median","qlow","qhigh")
  
  ggplot(dfquantiles, aes(x=Year,y=median, color=haz))+
    geom_line()+
    geom_ribbon(aes(ymin=dfquantiles$qlow, ymax=dfquantiles$qhigh,fill=dfquantiles$haz), linetype=2, alpha=0.1)
  
  #Change of strategy, I can do violin plots for both with change
  #first I do the violin for one variable
  ggplot(dfplotX, aes(x=factor(year), y=pixel, fill=factor(haz))) +
    geom_boxplot(notch=F,position=position_dodge(.9),alpha=.8,draw_quantiles = c(0.25, 0.5, 0.75),linewidth=0.8,scale="width",outlier.alpha = 0.4,trim=T)+
    
    
    
    
    scale_y_continuous(limits = c(-2,2),name="qsp change drought",
                       sec.axis = sec_axis( trans=~.*10, name="Second Axis"))
  scale_x_discrete(labels=c("1" = "100-200", "2" = "200-500",
                            "3" = "500-1000","4" = "1000-10 000",
                            "5" = "10 000-100 000","6" = ">100 000"),name="Catchment Area (km2)")+
    scale_fill_gradientn(
      colors=palet, n.breaks=6,limits=c(0.4,0.8)) +
    geom_text(data=data.frame(), aes(x=names(meds), y=meds-0.05, label=agUpA$upav), col='black', size=4,fontface="bold")+
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
  
  
  
  
  if (length(which(is.na(dfplot$pixel)))>0) dfplot=dfplot[-which(is.na(dfplot$pixel)),]
  if (length(which(is.infinite(dfplot$pixel)))>0) dfplot=dfplot[-which(is.infinite(dfplot$pixel)),]
  
  countfunction<-function(values, bounds){
    l1=length(which(values<bounds & values>0))
    l2=length(which(values>bounds))
    l3=length(which(values<(-bounds)))
    l4=length(which(values>(-bounds) & values<0))
    return(c(l2,l1,l4,l3))
  }
  dfq=aggregate(list(change=dfplot$pixel),
                by = list(Year=dfplot$year),
                FUN = function(x,...) c(bg=countfunction(x,0.1)))
  dfq <- do.call(data.frame, dfq)
  
  dfq2=aggregate(list(change=dfplot2$pixel),
                 by = list(Year=dfplot2$year),
                 FUN = function(x,...) c(bg=countfunction(x,5)))
  dfq2 <- do.call(data.frame, dfq2)
  
  data_long <- gather(dfq, group, value, change.bg1:change.bg4) %>%
    arrange(factor(Year)) %>% 
    mutate(x=factor(Year, levels=unique(Year)))
  
  
  colorz = c('change.bg4'='#d73027','change.bg3'='#fee090','change.bg2'='lightblue','change.bg1'='royalblue')
  labs=rev(c("< -0.1","-0.1 - 0","0 - 0.1", ">0.1"))
  legend_order=c('change.bg1', 'change.bg2', 'change.bg3','change.bg4')
  data_long$p2 <- relevel(data_long$position, 'change.bg1')
  data_long$position <- factor(data_long$group, levels=c('change.bg1', 'change.bg2', 'change.bg4','change.bg3'))
  
  ggplot(data_long, aes(fill=position, y=value, x=Year)) + 
    geom_bar(position="stack", stat="identity")+
    scale_fill_manual(values = colorz, breaks = legend_order,name="change in qsp",labels=labs)  + 
    scale_y_continuous(name="Number of catchments")+
    scale_x_continuous(name="Years") +
    # scale_x_discrete(labels= c("< -0.41","-0.41-0.2", "0.2-0.5","0.5-0.7","0.7-0.8", ">0.8"), name="KGE")+
    theme(axis.title=element_text(size=16, face="bold"),
          axis.text = element_text(size=13),
          panel.background = element_rect(fill = "white", colour = "white"),
          panel.grid = element_blank(),
          panel.border = element_rect(linetype = "solid", fill = NA, colour="black"),
          legend.title = element_text(size=14),
          legend.text = element_text(size=12),
          legend.position = "bottom",
          panel.grid.major = element_line(colour = "grey80"),
          panel.grid.minor.x = element_line(colour = "grey90",linetype="dashed"),
          legend.key = element_rect(fill = "transparent", colour = "transparent"),
          legend.key.size = unit(.8, "cm"))
  
  
  dfquantiles=aggregate(list(change=dfplot$pixel),
                        by = list(Year=dfplot$year),
                        FUN = function(x) c(mean=mean(x,na.rm=T),q1=quantile(x, q1, na.rm=T),q3=quantile(x, q2, na.rm=T)))
  
  dfquantiles <- do.call(data.frame, dfquantiles)
  names(dfquantiles)=c("Year", "mean","qlow","qhigh")
  
  
  
  dtfinal=c()
  xmm=c(-1,1)
  for( y in 1:length(seqy)){
    print(y)
    yx=seqy[y]
    dfp1=dfplot[which(dfplot$year==yx),]
    dfp1=dfp1[which(!is.na(dfp1$pixel)),]
    merde=kde(dfp1$pixel,bgridsize=1000,xmin=xmm[1],xmax=xmm[2],h=0.02)
    
    dfqy=dfquantiles[which(dfquantiles$Year==yx),]
    ziz=data.frame(ev.points=merde$eval.points,estimates=merde$estimate)
    ziz$estimates=ziz$estimates/sum(ziz$estimates)
    sum(ziz$estimates)
    dtz=cbind(rep(yx,1000),ziz)
    dtz[,3][which(dtz[,2]<dfqy$qlow)]=NA
    dtz[,3][which(dtz[,2]>dfqy$qhigh)]=NA
    dtfinal=rbind(dtfinal,dtz)
  }
  
  dfqplus=seq(period[1],period[2],by=0.1)
  dfql=interp1(dfquantiles$Year,dfquantiles$qlow,dfqplus)
  dfqh=interp1(dfquantiles$Year,dfquantiles$qhigh,dfqplus)
  dfqplus=data.frame(dfqplus,dfql,dfqh)
  names(dfqplus)=c("Year","qlow","qhigh")
  
  #alpha=1+(log(values)))
  
  names(dtfinal)=c("Year","Change","values")
  ptref=data.frame(x=1950,y=0)
  
  palet=c(hcl.colors(11, palette = "YlGnBu", alpha = NULL, rev = T, fixup = TRUE))
  
  plo=ggplot(dtfinal, aes(x = Year, y = Change))+
    geom_raster(aes(fill=values),interpolate = T)+
    scale_fill_gradientn(colours = palet,
                         n.breaks=10,guide='coloursteps',na.value = "transparent",oob = scales::squish,limits=c(0,0.005))+
    scale_y_continuous(
      n.breaks = 10,name= "Change in Return Level",
      limit=lims,expand = c(0, 0))+
    scale_x_continuous(
      n.breaks = 10,name= "Years",
      limit=c(1949, 2020),expand = c(0, 0))+
    geom_line(data=dfquantiles,aes(x=Year,y=mean),size=2,col="red",alpha=0.9)+
    geom_line(data=dfquantiles,aes(x=Year,y=qlow),size=1,col="grey",linetype="longdash")+
    geom_line(data=dfquantiles,aes(x=Year,y=qhigh),size=1,col="grey",linetype="longdash")+
    geom_hline(yintercept=0)+
    geom_point(data=ptref,aes(x=x, y = y), col="black",fill="white",size=7,pch=21)+
    theme_bw() +
    theme(
      panel.spacing = unit(0.3, "lines"),
      strip.text.x = element_text(size = 8),
      legend.position = "none",
      legend.key.width= unit(2.5, 'cm'),
      legend.text =element_text(size=14),
      panel.border = element_rect(linetype = "solid", fill = NA, colour="black"),
      legend.title = element_text(size=14),
      plot.title = element_text(size=20)
    )+
    ggtitle(unique(outloc$name))
  
  plo
  return(plo)
}

#Another function that fits sen slope on annual max and annual min (It will be with the peaks but almost the same)
#find inspiration in other functions


#Season computation functions
seasony=function(x){
  theta=x*(2*pi/365.25)
  #plot(cos(theta),sin(theta),xlim=c(-1,1),ylim=c(-1,1))
  xi=1/(length(theta))*sum(cos(na.omit(theta)))
  yi=1/(length(theta))*sum(sin(na.omit(theta)))
  if (xi<=0){
    Di=(atan(yi/xi)+pi)*(365.25/(2*pi))
  }else if(xi>0 & yi>=0){
    Di=(atan(yi/xi))*(365.25/(2*pi))
  }else if(xi>0 & yi<0){
    Di=(atan(yi/xi)+2*pi)*(365.25/(2*pi))
  }
  R=sqrt(xi^2+yi^2)
  return(c(Di,R))
}
season1=function(x){
  l1=length(which(!is.na(x)))
  if(l1>0){
    x=x[which(!is.na(x))]
    theta=x*(2*pi/365.25)
    # plot(theta)
    
    xi=1/(length(theta))*sum(cos(theta))
    yi=1/(length(theta))*sum(sin(theta))
    if (xi<=0){
      Di=(atan(yi/xi)+pi)*(365.25/(2*pi))
    }else if(xi>0 & yi>=0){
      Di=(atan(yi/xi))*(365.25/(2*pi))
    }else if(xi>0 & yi<0){
      Di=(atan(yi/xi)+2*pi)*(365.25/(2*pi))
    }
    R=sqrt(xi^2+yi^2)
  }else{Di=NA}
  return(Di)
}
season2=function(x){
  l1=length(which(!is.na(x)))
  if(l1>0){
    x=x[which(!is.na(x))]
    theta=x*(2*pi/365.25)
    # plot(theta)
    
    xi=1/(length(theta))*sum(cos(theta))
    yi=1/(length(theta))*sum(sin(theta))
    if (xi<=0){
      Di=(atan(yi/xi)+pi)*(365.25/(2*pi))
    }else if(xi>0 & yi>=0){
      Di=(atan(yi/xi))*(365.25/(2*pi))
    }else if(xi>0 & yi<0){
      Di=(atan(yi/xi)+2*pi)*(365.25/(2*pi))
    }
    R=sqrt(xi^2+yi^2)
  }else{Di=NA
  R=NA}
  return(R)
}

#Concentration of floods around peaks
freqsum=function(piks,seriag){
  bo1=piks-60
  bo2=piks+60
  ss1=0
  ss2=0
  ss=sum(seriag$Freq[max((piks-60),1):min((piks+60),366)])/sum(seriag$Freq)
  
  if (bo1<1) ss1=sum(seriag$Freq[(366+bo1):366])/sum(seriag$Freq)
  if (bo2>366)ss2=sum(seriag$Freq[1:(bo2-366)])/sum(seriag$Freq)
  
  ss=ss+ss1+ss2
  return(ss)
}

#Seasonality of extremes, 30 years moving window
RunningSeason = function(series, timestamps, windowSize, nyears,yplot){
  rnseas = matrix(nrow=nyears,ncol= 1)
  rncon = matrix(nrow=nyears,ncol= 1)
  rnnp = matrix(nrow=nyears,ncol= 1)
  ssave = matrix(nrow=nyears,ncol= 1)
  rnpx = matrix(nrow=nyears,ncol= 4)
  saglist=list()
  
  dx = floor(windowSize/2)
  l = length(series)
  dayvec=data.frame(day=seq(1,366))
  id=0
  for (ii in seq(1,l, by=366)){
    id=id+1
    minindx = max(ii - dx, 1);
    maxindx = min(ii + dx, l);
    subsrs = series[minindx:maxindx]
    subsrs = subsrs[which(!is.na(subsrs))]
    if (length(!is.na(subsrs))>1){
      seasonsub=seasony(subsrs)
      
      #new seasonality method can detect several peaks
      seriag=table(subsrs)
      days=as.numeric(names(seriag))
      rdays=data.frame(d=seq(1:366))
      seriag=data.frame(days,seriag)
      ziz=rdays$d[match(rdays$d,days)]
      seriag=seriag[ziz,]
      seriag$days=rdays$d
      
      seriag=seriag[order(seriag$days),]
      seriag$Freq[which(is.na(seriag$Freq))]=0
      seriag$fravg=movingFun(seriag$Freq, n= 30, fun = mean, circular = TRUE)
      
      seriag$fravg=seriag$fravg/max(seriag$fravg)
      th=0.1
      pks <- findpeaks(seriag$fravg,minpeakdistance = 90, minpeakheight = th, zero = "-", npeaks = 1)
      #if no peaks found change data order
      if (is.null(pks)){
        x=seriag$fravg
        lng=length(x)
        hn <- 90
        x <- c(x[(lng-hn+1):lng], x[1:(lng-hn+1)])
        pks <- findpeaks(x,minpeakdistance = 90, minpeakheight = th, zero = "-")
        #deshifting
        pks[,c(2:4)]=pks[,c(2:4)]-90
        pks[,c(2:4)][which(pks[,c(2:4)]<0)]=pks[,c(2:4)][which(pks[,c(2:4)]<0)]+366
      }
      np=length(pks[,1])
      rnpx[id,1] = pks[1,2]
      
      rnpx[id,2]=NA
      rnpx[id,3]=NA
      rnpx[id,4]=NA
      ss=freqsum(pks[1,2],seriag)
      npx=1
      while (ss<0.75){
        npx=npx+1
        seriag$fravg=seriag$fravg/max(seriag$fravg)
        #can 2 peaks include more than 75% of flood days?
        pks <- findpeaks(seriag$fravg,minpeakdistance = 90, minpeakheight = 0.1, zero = "-", npeaks = npx)
        npv=length(pks[,1])
        if(npv==2){
          rnpx[id,2] = pks[2,2]
          pk1 = pks[1,2]
          pk2 = pks[2,2]
          dp=diff(pks[,2])
          dp[which(dp<=-182)] = 365.25+dp[which(dp<=-182)]
          dp[which(dp>=182)] = 365.25-dp[which(dp>=182)]
          wd=which(abs(dp)<60)
          if (length(wd)>0){
            rnpx[id,wd]=seasony(c(rnpx[id,wd],rnpx[id,wd+1]))[1]
            rnpx[id,wd+1]=NA
            ss=freqsum(rnpx[id,wd],seriag)
          }else{
            ss=ss+freqsum(pk2,seriag)
          }
        }
        if(npv==3){
          pk1 = pks[1,2]
          pk2 = pks[2,2]
          pk3 = pks[3,2]
          dp=diff(pks[,2])
          dp[which(dp<=-182)] = 365.25+dp[which(dp<=-182)]
          dp[which(dp>=182)] = 365.25-dp[which(dp>=182)]
          wd=which(abs(dp)<90)
          if (length(wd)>0){
            rnpx[id,wd]=seasony(c(rnpx[id,wd],rnpx[id,wd+1]))[1]
            rnpx[id,wd+1]=NA
            ss=freqsum(rnpx[id,wd])
            
          }else{
            ss=ss+freqsum(pk2,seriag)+freqsum(pk3,seriag)
          }
        }
        if(npx>=4) break
      }
      np=length(which(!is.na(rnpx[id,])))
      ssave[id]=ss
      rnnp[id] = np
      rnseas[id] = seasonsub[1]
      rncon[id] = seasonsub[2]
      saglist=c(saglist,list(seriag))
      
    }else{
      rnseas[id] = NA
      rncon[id] = NA
      ssave[id]= NA
      rnnp[id] = NA
    }
  }
  
  
  timeD=timestamps[seq(1,l, by=366)]
  rnpx=as.data.frame(rnpx)
  names(rnpx)=c("pk1","pk2","pk3","pk4")
  
  
  # Plot section
  yplot=70
  seriag=saglist[[yplot]]
  th=max(seriag$fravg)/2
  rnpx$y=rep(1,71)
  plot=ggplot(seriag) +
    geom_segment(aes(x=days,xend=days,y=0,yend=fravg),col="darkblue") +
    # Make custom panel grid
    geom_hline(
      aes(yintercept = y), 
      data.frame(y = c(0.25,0.50,0.75,1)),
      color = "lightgrey",lwd=1
    ) + 
    geom_vline(xintercept = seasonsub[1], col="red", lwd=2)+
    #geom_hline(yintercept = th, col="blue", lwd=2)+
    geom_point(data=rnpx[yplot,],aes(x=pk1,y=y),pch=16,size=4, col="royalblue")+
    geom_point(data=rnpx[yplot,],aes(x=pk2,y=y),pch=16,size=4, col="chartreuse3")+
    coord_polar() +
    scale_x_continuous(limits = c(0,366), breaks=seq(366,30.5,-30.5),
                       label=direction_labeller)+
    scale_y_continuous(lim=c(0,1))+
    theme(axis.title=element_text(size=16, face="bold"),
          axis.text = element_text(size=14),
          panel.background = element_rect(fill = "white", colour = "white"),
          #panel.grid = element_blank(),
          panel.grid.major.x = element_blank(),
          panel.border = element_rect(linetype = "solid", fill = NA, colour="black"),
          legend.title = element_text(size=14),
          legend.text = element_text(size=12),
          # panel.grid.major = element_line(colour = "grey70"),
          panel.grid.minor = element_line(colour = "grey90"),
          legend.key = element_rect(fill = "transparent", colour = "transparent"),
          legend.key.size = unit(.8, "cm"))+
    labs( y = "Frequency", x= "Days")
  plot
  return(list(data.frame(timeD,rnseas,rncon, rnnp, ssave,rnpx),plot))
}

#for seasonal plots
direction_labeller <- function(x){
  ifelse(x %% 30.5 == 0, c('Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct',"Nov",'Dec')[1+(as.integer(x/30.5) %% 12)], '')
}

# Function to plot the distribution of pixels for the four seasons
plotSeasonFlood=function(seasonP,nco,tsize,osize){
  
  out=list()
  
  out[[1]]=ggplot(w2) +
    geom_sf(fill="white")+
    geom_sf(data=seasonP,aes(color=factor(summer),geometry=geometry),alpha=1,size=0.1,stroke=0,shape=15)+ 
    geom_sf(fill=NA, color="grey") +
    coord_sf(xlim = c(min(nco[,1]),max(nco[,1])), ylim = c(min(nco[,2]),max(nco[,2])))+
    scale_color_manual(
      values=c("1" ="darkred"),
      na.value="grey")+
    guides(colour = guide_legend(override.aes = list(size = 10)))+
    # guides(color = guide_coloursteps(barwidth = 1, barheight = 10))+
    # scale_fill_manual(values=paletx,
    #                   na.value="grey", name="Date") +
    labs(x="Longitude", y = "Latitude")+
    theme(axis.title=element_text(size=tsize),
          panel.background = element_rect(fill = "aliceblue", colour = "grey1"),
          panel.border = element_rect(linetype = "solid", fill = NA, colour="black"),
          legend.title = element_text(size=tsize),
          legend.text = element_text(size=osize),
          panel.grid.major = element_line(colour = "grey70"),
          panel.grid.minor = element_line(colour = "grey90"),
          legend.position = "none",
          legend.key = element_rect(fill = "transparent", colour = "transparent"),
          legend.key.size = unit(.8, "cm"))+
    theme(plot.margin = margin(0.2,0.2,0.2,0.2, "cm"))+
    ggtitle("Summer")
  
  out[[2]]=ggplot(w2) +
    geom_sf(fill="white")+
    geom_sf(data=seasonP,aes(color=factor(autumn),geometry=geometry),alpha=1,size=0.1,stroke=0,shape=15)+ 
    geom_sf(fill=NA, color="grey") +
    coord_sf(xlim = c(min(nco[,1]),max(nco[,1])), ylim = c(min(nco[,2]),max(nco[,2])))+
    scale_color_manual(
      values=c("1" ="darkgreen"),
      na.value="grey")+
    guides(colour = guide_legend(override.aes = list(size = 10)))+
    # guides(color = guide_coloursteps(barwidth = 1, barheight = 10))+
    # scale_fill_manual(values=paletx,
    #                   na.value="grey", name="Date") +
    labs(x="Longitude", y = "Latitude")+
    theme(axis.title=element_text(size=tsize),
          panel.background = element_rect(fill = "aliceblue", colour = "grey1"),
          panel.border = element_rect(linetype = "solid", fill = NA, colour="black"),
          legend.title = element_text(size=tsize),
          legend.text = element_text(size=osize),
          panel.grid.major = element_line(colour = "grey70"),
          panel.grid.minor = element_line(colour = "grey90"),
          legend.position = "none",
          legend.key = element_rect(fill = "transparent", colour = "transparent"),
          legend.key.size = unit(.8, "cm"))+
    theme(plot.margin = margin(0.2,0.2,0.2,0.2, "cm"))+
    ggtitle("Autumn")
  
  out[[3]]=ggplot(w2) +
    geom_sf(fill="white")+
    geom_sf(data=seasonP,aes(color=factor(winter),geometry=geometry),alpha=1,size=0.1,stroke=0,shape=15)+ 
    geom_sf(fill=NA, color="grey") +
    coord_sf(xlim = c(min(nco[,1]),max(nco[,1])), ylim = c(min(nco[,2]),max(nco[,2])))+
    scale_color_manual(
      values=c("1" ="blue"),
      na.value="grey")+
    guides(colour = guide_legend(override.aes = list(size = 10)))+
    # guides(color = guide_coloursteps(barwidth = 1, barheight = 10))+
    # scale_fill_manual(values=paletx,
    #                   na.value="grey", name="Date") +
    labs(x="Longitude", y = "Latitude")+
    theme(axis.title=element_text(size=tsize),
          panel.background = element_rect(fill = "aliceblue", colour = "grey1"),
          panel.border = element_rect(linetype = "solid", fill = NA, colour="black"),
          legend.title = element_text(size=tsize),
          legend.text = element_text(size=osize),
          panel.grid.major = element_line(colour = "grey70"),
          panel.grid.minor = element_line(colour = "grey90"),
          legend.position = "none",
          legend.key = element_rect(fill = "transparent", colour = "transparent"),
          legend.key.size = unit(.8, "cm"))+
    theme(plot.margin = margin(0.2,0.2,0.2,0.2, "cm"))+
    ggtitle("Winter")
  
  out[[4]]=ggplot(w2) +
    geom_sf(fill="white")+
    geom_sf(data=seasonP,aes(color=factor(spring),geometry=geometry),alpha=1,size=0.1,stroke=0,shape=15)+ 
    geom_sf(fill=NA, color="grey") +
    coord_sf(xlim = c(min(nco[,1]),max(nco[,1])), ylim = c(min(nco[,2]),max(nco[,2])))+
    scale_color_manual(
      values=c("1" ="magenta"),
      na.value="grey")+
    guides(colour = guide_legend(override.aes = list(size = 10)))+
    # guides(color = guide_coloursteps(barwidth = 1, barheight = 10))+
    # scale_fill_manual(values=paletx,
    #                   na.value="grey", name="Date") +
    labs(x="Longitude", y = "Latitude")+
    theme(axis.title=element_text(size=tsize),
          panel.background = element_rect(fill = "aliceblue", colour = "grey1"),
          panel.border = element_rect(linetype = "solid", fill = NA, colour="black"),
          legend.title = element_text(size=tsize),
          legend.text = element_text(size=osize),
          panel.grid.major = element_line(colour = "grey70"),
          panel.grid.minor = element_line(colour = "grey90"),
          legend.position = "none",
          legend.key = element_rect(fill = "transparent", colour = "transparent"),
          legend.key.size = unit(.8, "cm"))+
    theme(plot.margin = margin(0.2,0.2,0.2,0.2, "cm"))+
    ggtitle("Spring")
  
  return(out)
  
}

get_density <- function(x, y, ...) {
  dens <- MASS::kde2d(x, y, ...)
  ix <- findInterval(x, dens$x)
  iy <- findInterval(y, dens$y)
  ii <- cbind(ix, iy)
  return(dens$z[ii])
}

UpAopen=function(dir,outletname,Sloc_final){
  ncbassin=paste0(dir,outletname)
  ncb=nc_open(ncbassin)
  name.vb=names(ncb[['var']])
  namev=name.vb[2]
  #time <- ncvar_get(ncb,"time")
  
  #timestamp corretion
  name.lon="lon"
  name.lat="lat"
  londat = ncvar_get(ncb,name.lon) 
  llo=length(londat)
  latdat = ncvar_get(ncb,name.lat)
  lla=length(latdat)
  start=c(1,1)
  count=c(llo,lla)
  
  
  londat = ncvar_get(ncb,name.lon,start=start[1],count=count[1]) 
  llo=length(londat)
  latdat = ncvar_get(ncb,name.lat,start=start[2],count=count[2])
  lla=length(latdat)
  outlets = ncvar_get(ncb,namev,start = start, count= count) 
  outlets=as.vector(outlets)/1000000
  outll=expand.grid(londat,latdat)
  lonlatloop=expand.grid(c(1:llo),c(1:lla))
  outll$upa=outlets
  outll$idlo=lonlatloop$Var1
  outll$idla=lonlatloop$Var2
  
  outll$idlalo=paste(outll$idlo,outll$idla,sep=" ")
  outfinal=inner_join(outll, Sloc_final, by="idlalo")
  return (outfinal)
}

# Load inputs from HPC computation ----------------------------------------




# Pre-loaded results ------------------------------------------------------



## Spatial data for catchments ----
#Set data directory
hydroDir<-("D:/tilloal/Documents/LFRuns_utils/data")


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
# st_write(Catfchier, "catcheck.shp")
unikout=Results$catrest$catlist

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

#Unsure running every scenario separately is the best approach

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

dis_path<-paste0(main_path,'HERA_Rstat/')
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

stid=303662
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



#TSEVA for each scenario
catmat=Catf7[which(Catf7$pointid==stid),]
minPeakDistanceInDays=7
interflag=0
frosttime=NA
trans="ori"
tail="high"



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
#####
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

## 1. selection of data to be plotted ----
#Loading saved results in .Rdata

load(file=paste0(hydroDir,"/TSEVA_hybas/outputs/parCat07_Histo_flood_1950_2020.Rdata"))
load(file=paste0(hydroDir,"/TSEVA_hybas/outputs/ResCat07_Histo_flood_1950_2020.Rdata"))
Flood_histo=Results

load(file=paste0(hydroDir,"/TSEVA_hybas/outputs/parCat07_SocCF_flood_1950_2020.Rdata"))
load(file=paste0(hydroDir,"/TSEVA_hybas/outputs/ResCat07_SocCF_flood_1950_2020.Rdata"))
Flood_SocCF=Results
length(which(Results$catrest[,3]>0))


load(file=paste0(hydroDir,"/TSEVA_hybas/outputs/parCat07_ResCF_flood_1950_2020.Rdata"))
load(file=paste0(hydroDir,"/TSEVA_hybas/outputs/ResCat07_ResCF_flood_1950_2020.Rdata"))
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



bigfish=Flood_histo$RetPerGPD
bigfish=data.frame(bigfish[,70],unikout=Flood_ResCF$catrest$catlist)
#reproduce what was done before for all catchments


#look at a pixel with high negative change

#create a time vector with years
dates <- seq.Date(ymd("1951-06-01"), ymd("2020-06-01"), by = "year")
dates= as.POSIXct(dates, format = "%Y-%m-%d")
dfall=c()
for (id in 1:length(Station_data_IDs)){
  
  print(id)
  stid=Station_data_IDs[id]
  # stid=362439
  catch=stid
  
  RLGPDH=as.vector(t(RLGPD_H[which(RLGPD_H$unikout==stid),-71]))
  RLGPDR=as.vector(t(RLGPD_RCF[which(RLGPD_RCF$unikout==stid),-71]))
  RLGPDS=as.vector(t(RLGPD_SCF[which(RLGPD_SCF$unikout==stid),-71]))
  
  
  pid=Peaks_H[which(Peaks_H$catch==stid),]
  psid=Peaks_SCF[which(Peaks_SCF$catch==stid),]
  prid=Peaks_RCF[which(Peaks_RCF$catch==stid),]
  
  
  # plot(pid$time,pid$value,pch=16,type="p")
  # points(psid$time,psid$value,col="red",type="p",pch=16)
  # points(prid$time,prid$value,col="blue",type="p",pch=16)
  # lines(dates,RLGPDH,col=1,lwd=2)
  # lines(dates,RLGPDR,col="blue",lwd=2)
  # lines(dates,RLGPDS,col="red",lwd=2)
  
  
  #1 I remove the climate trend
  RLGPDSdt=RLGPDS-(RLGPDS-RLGPDS[1])
  RLGPDRdt=RLGPDR-(RLGPDS-RLGPDS[1])
  RLGPDHdt=RLGPDH-(RLGPDS-RLGPDS[1])
  
  #Now I remove the socioeconomic trend
  RLGPDRdtp=RLGPDRdt-(RLGPDRdt-RLGPDSdt)
  RLGPDHdtp=RLGPDHdt-(RLGPDRdt-RLGPDSdt)

  #And the reservoir
  RLGPDHdtpx=RLGPDHdtp-(RLGPDHdtp-RLGPDRdtp)

  
  #change from climate
  Climtrend=(RLGPDS-RLGPDS[1])

  #change from socioeco
  #Soctrend=(RLGPDRdt-RLGPDRdt[1])+(RLGPDRdt-RLGPDSdt)
  Soctrend=(RLGPDRdt-RLGPDSdt)

  
  
  #change from reservoirs
  Restrend=(RLGPDHdtp-RLGPDRdtp)
  #Restrend=(RLGPDHdtp-RLGPDHdtp[1])+(RLGPDHdtp-RLGPDRdtp)
  points(Restrend,col=2)
  
  #  compute change in specific discharge instead of relative change
  upag=match(stid ,UpArea$outlets)
  upareaL=UpArea$upa[upag]
  
  df <- data.frame(
    time = c(year(timeDays[1]):year(timeDays[length(timeDays)])),
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

ggplot(df_l2020, aes(x=factor(driver), y=value)) +
  geom_boxplot(notch=F,position=position_dodge(.9),alpha=.8,aes(fill=driver),linewidth=0.8,outlier.alpha = 0.4)+
  scale_y_continuous(limits = c(-200,200),name="Contribution to change in 10y Flood RL(l/s/km2) (1951-2020)",breaks = seq(-100,100,by=10),minor_breaks = seq(-100,100,5))+
  scale_x_discrete(labels=c("Dres" = "Reservoirs", "Dsoc" = "Land use + water use",
                            "Dclim" = "Climate"),name="Driver")+
  coord_cartesian(ylim = c(-25,25))+
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





#load upstream area
main_path = 'D:/tilloal/Documents/06_Floodrivers/'
valid_path = paste0(main_path,'DataPaper/')
outletname="/GIS/upArea_European_01min.nc"
dir=valid_path
outf$idlalo=paste(outf$idlo, outf$idla, sep=" ")
UpArea=UpAopen(valid_path,outletname,outf)
head(UpArea)


#Load rasterized hydroregions

GridHR=raster( paste0(hydroDir,"/HydroRegions_raster_WGS84.tif"))
GHR=as.data.frame(GridHR,xy=T)
GHR=GHR[which(!is.na(GHR[,3])),]
GHR$llcoord=paste(round(GHR$x,4),round(GHR$y,4),sep=" ") 
GHR_riv=inner_join(GHR,outf,by= c("llcoord"="latlong"))
GHshpp <- read_sf(dsn ="Z:/ClimateRun4/nahaUsers/tilloal/HydroRegions/her_all_adjusted.shp")
HydroRsf=fortify(GHshpp) 
### Changes in RLs or RPs at pixel and catchment levels ----



datarF=RLGPD_H
datarF$IRES=0
Plot.change=plotchangemaps_qsp(basemap,catmap=catmap,datarF, UpArea,GHR_riv,HydroRsf, law="GPD",type="RLchange",period=c(1951,2020),parlist,valuenames,haz=haz)
Plot.change[[3]]

#Deep comparison between 1951 for all catchments

data_fy=df_2020=dfall[which(dfall$time==2020), ]


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

#load reservoir influence
reservoirdog=outletopen(paste0(hydroDir,"/reservoirs"),"res_ratio_diff_2020-1951")
reservoirdog$llcoord=paste(round(reservoirdog$Var1,4),round(reservoirdog$Var2,4),sep=" ") 
reservoir_bg=inner_join(reservoirdog,outf,by=c("llcoord"="latlong"))


hybmatch=match(data_fy$catch,cst7$outlets)
data_fy$hybasID=cst7$HYBAS_ID[hybmatch]

driver="socioeconomic changes"
driver="climate"
#catchment plot with difference
data=right_join(catmap,data_fy,by = c("outlets"="catch"))
title=paste0("10y flood RL change due to ", driver ,"(1951-2020)")
legend="change(l/s/km2)"
palet=c(hcl.colors(11, palette = "RdYlBu", alpha = NULL, rev = F, fixup = TRUE))
br=c(-100,-80,-50,-30,-20,-10,-8,-6,-4,-2,0,2,4,6,8,10,20,30,50,80,100)
br=c(-100,-80,-50,-20,-10,0,10,20,30,50,80,100)
limi=c(-50,50)
trans=scales::modulus_trans(.8)
colNA="darkgrey"

ggplot(basemap) +
  geom_sf(fill="white")+
  geom_sf(data=data,aes(fill=Dsoc,geometry=geometry),color="transparent")+ 
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

#maybe I should compute this change is qsp


#compare this with the chance from the SocCF scenario

Plot.change=plotchangemaps(basemap,catmap=catmap,RLGPD_SCF, law="GPD",type="RLchange",period=c(1951,2020),parlist,valuenames,haz="flood")
Plot.change[[3]]

#success

#check individual catchments with unexpected changes to see if thresold is problematic (E.g murcia)


#Look by hydroregions as well could be cool


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






