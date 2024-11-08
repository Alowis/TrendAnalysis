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
library(RtsEva)
#  Function declaration ---------------------------------------------------


###########################  FUNCTIONS   #################################################
timeAndSeries=timeAndSeriesSCF
TrendTh=thresh$th
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
    # if (is.na(TrendTh)){
    #   TrendTh=try(tsEvaFindTrendThreshold(series, timeStamps, timeWindow),T)
    #   if(length(TrendTh)==0){
    #     TrendTh=0.1
    #   }
    # }
    trasfData = tsEvaTransformSeriesToStationaryPeakTrend( timeStamps, series, timeWindow, TrendTh);
    #trasfData = tsEvaTransformSeriesToStationaryMMXTrend( timeStamps, series, timeWindow);
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
  print(ms[c(1:100),])
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


tsEvaTransformSeriesToStationaryMMXTrend<-function (timeStamps, series, timeWindow) 
{
  tserie=data.frame(timeStamps,series)
  monthly_max=tserie %>% 
    mutate(month = floor_date(timeStamps, "month")) %>% 
    group_by(month) %>% 
    summarise(max_value = max(series))
  
  
  serieb <- series
  tm=na.omit(match(as.Date(as.character(monthly_max$month)),as.Date(timeStamps)))
  print(serieb[tm])
  serieb[-tm] <- NA
  rs <- tsEvaDetrendTimeSeries(timeStamps, serieb, timeWindow)
  detrendSeries <- series - rs@trendSeries
  detrendSerie1 <- serieb - rs@trendSeries
  qd2 <- min(detrendSerie1, na.rm = T)
  nRunMn <- rs@nRunMn
  varianceSeries <- tsEvaNanRunningVariance(detrendSerie1, 
                                            nRunMn)
  varianceSeries <- tsEvaNanRunningMean(varianceSeries, ceiling(nRunMn/2))
  stdDevSeries1 <- varianceSeries^0.5
  stdDevSeries <- stdDevSeries1
  avgStdDev <- mean(stdDevSeries)
  S <- 2
  N <- timeWindow * 4
  stdDevError <- avgStdDev * (2 * S^2/N^3)^(1/4)
  statSeries <- detrendSeries/stdDevSeries
  xtremS <- statSeries
  xtremS <- na.omit(xtremS)
  allS <- na.omit(serieb)
  statSer3Mom <- tsEvaNanRunningStatistics(statSeries, nRunMn)$rn3mom
  statSer4Mom <- tsEvaNanRunningStatistics(statSeries, nRunMn)$rn4mom
  statSer3Mom <- tsEvaNanRunningMean(statSer3Mom, ceiling(nRunMn))
  statSer4Mom <- tsEvaNanRunningMean(statSer4Mom, ceiling(nRunMn))
  trendError <- mean(stdDevSeries)/N^0.5
  trasfData <- list(runningStatsMulteplicity = nRunMn, stationarySeries = statSeries, 
                    trendSeries = rs@trendSeries, trendSeriesNonSeasonal = NULL, 
                    trendError = trendError, stdDevSeries = stdDevSeries, 
                    stdDevSeriesNonSeasonal = NULL, stdDevError = stdDevError * 
                      rep(1, length(stdDevSeries)), timeStamps = timeStamps, 
                    nonStatSeries = series, statSer3Mom = statSer3Mom, statSer4Mom = statSer4Mom)
  return(trasfData)
}



dir=hydroDir
# outletopen=function(dir,outletname,nrspace=rep(NA,5)){
# 
#   #opening one
#   outletname="efas_rnet_100km_01min"
#   ncbassin=paste0(dir,"/",outletname,".nc")
#   ncb=nc_open(ncbassin)
#   name.vb=names(ncb[['var']])
#   namev=name.vb[1]
#   #time <- ncvar_get(ncb,"time")
#   
#   #timestamp corretion
#   name.lon="lon"
#   name.lat="lat"
#   if (!is.na(nrspace[1])){
#     start=as.numeric(nrspace[c(2,4)])
#     count=as.numeric(nrspace[c(3,5)])-start+1
#   }else{
#     start=c(1,1)
#     count=c(llo,lla)
#   }
#   
#   
#   
#   londat = ncvar_get(ncb,name.lon,start=start[1],count=count[1]) 
#   llo=length(londat)
#   latdat = ncvar_get(ncb,name.lat,start=start[2],count=count[2])
#   lla=length(latdat)
#   outlets = ncvar_get(ncb,namev,start = start, count= count) 
#   image(outlets)
#   
#   #opening 2
#   outletname=filename
#   ncbassin=paste0(dir,"/",outletname,".nc")
#   nc2=nc_open(ncbassin)
#   name.vb=names(nc2[['var']])
#   namev=name.vb[1]
#   #time <- ncvar_get(ncb,"time")
#   
#   #timestamp corretion
#   name.lon="lon"
#   name.lat="lat"
#   
#   londat = ncvar_get(ncb,name.lon,start=start[1],count=count[1]) 
#   llo=length(londat)
#   latdat = ncvar_get(ncb,name.lat,start=start[2],count=count[2])
#   lla=length(latdat)
#   
#   start = c(1,1,1)
#   count= c(llo,lla,1)
#   discharge = ncvar_get(nc2,namev,start = start, count= count) 
#   image(discharge)
#   
#   
#   
#   
#   rownames(outlets) <- round(londat,4)
#   colnames(outlets) <- round(latdat,4)
#   
#   outX=outlets
#   for (i in 1:length(londat)){
#     for (j in 1:length(latdat)){
#       print(i)
#       print(j)
#       outX[i,j]=paste(rownames(outX)[i],colnames(outX)[j],sep=" ")
#     }
#   }
# 
#   vec <- as.vector(outlets)
#   
#   
#   veX=as.vector(outX)
#   veX2=veX[-which(is.na(vec))]
#   llon=c()
#   llat=c()
#   for (k in 1:length(veX2)){
#     print(k)
#     lon <- as.numeric(strsplit(veX2[k], " ")[[1]][1])
#     lat <- as.numeric(strsplit(veX2[k], " ")[[1]][2])
#     llon=c(llon,lon)
#     llat=c(llat,lat)
#   }
# 
#   outll=expand.grid(londat,latdat)
#   outll$bg=veX
#   lonlatloop=expand.grid(c(1:llo),c(lla:1))
#   
#   
#   
#   # verif=data.frame(dspr$Var1,lonlatloop$Var1)
#   # verif$f=verif$lonlatloop.Var1-verif$dspr.Var1
#   outll$idlo=lonlatloop$Var1
#   outll$idla=lonlatloop$Var2
#   outlx=outll[which(!is.na(vec)),]
#   
#   vec <- as.vector(outlets)
#  
#   
# 
#   dis2=discharge[which(!is.na(vec))]
#   veX2=data.frame(dis2,llon,llat)
#   
#   outll=data.frame(dis2,outlx)
#   
#   # 
#   # plot(outll$Var2,veX2$llat)
#   # 
#   # ggplot(data=outll,aes(x=idlo,y=idla,col=dis2))+
#   #   scale_color_gradient(limits=c(0,100),na.value="darkred")+
#   #   geom_point()
#   # 
#   # ggplot(data=veX2,aes(x=llon,y=llat,col=dis2))+
#   #   scale_color_gradient(limits=c(0,100),na.value="darkred")+
#   #   geom_point()
#   # 
#   # df=df[which(!is.na(df$Value)),]
#   # 
#   # plot(outll$Var1)
#   # plot(df$Row_Index,df$Col_Index)
#   # 
#   # 
#   # plot(df$Col)
#   
#   return (outll)
# }

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
disNcopen=function(fname,dir,outloc){
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
  
  outllplus=matrix(-9999, nrow = lt*length(outloc[,1]), ncol = 5)
  outllplus=as.data.frame(outllplus)
  for ( idc in 1:length(outloc[,1])){
    print(idc)
    idm=1+(idc-1)*lt
    outloc[idc,]
    start=c(outloc$idlo[idc],outloc$idla[idc],1)
    count=c(1,1,lt)
    outlets = ncvar_get(ncd,namev,start = start, count= count)
    outlets=as.vector(outlets)
    outid=rep(outloc[idc,1],length(time))
    #lonlatloop=expand.grid(c(1:lla),c(1:lt))
    lon=rep(londat[idc],length(time))
    lat=rep(latdat[idc],length(time))
    outll=data.frame(outlets,outid,lon,lat,time)
    names(outllplus)=names(outll)
    fck=length(c(idm:(idm+lt-1)))
    outllplus[c(idm:(idm+lt-1)),]= outll
  }
  return (outllplus)
}
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
  lon=rep(outloc$Var1[idc],length(time))
  lat=rep(outloc$Var2[idc],length(time))
  outll=data.frame(outlets,outid,lon,lat,time)
  
  return (outll)
}
ComputeReturnLevels<-function(nonStationaryEvaParams, RPgoal, timeIndex){
  
  #GEV
  epsilonGEV <- nonStationaryEvaParams[[1]]$parameters$epsilon
  sigmaGEV <- mean(nonStationaryEvaParams[[1]]$parameters$sigma[timeIndex])
  muGEV <- mean(nonStationaryEvaParams[[1]]$parameters$mu[timeIndex])
  dtSampleYears <- nonStationaryEvaParams[[1]]$parameters$timeDeltaYears
  
  #GPD
  epsilonGPD <- nonStationaryEvaParams[[2]]$parameters$epsilon
  sigmaGPD <- mean(nonStationaryEvaParams[[2]]$parameters$sigma[timeIndex])
  thresholdGPD <- mean(nonStationaryEvaParams[[2]]$parameters$threshold[timeIndex])
  nPeaks <- nonStationaryEvaParams[[2]]$parameters$nPeaks
  thStart <- nonStationaryEvaParams[[2]]$parameters$timeHorizonStart
  thEnd <- nonStationaryEvaParams[[2]]$parameters$timeHorizonEnd
  sampleTimeHorizon <- as.numeric((thEnd - thStart)/365.2425)
  
  if (nonStationaryEvaParams[[1]]$method=="No fit"){
    print("could not fit EVD to this pixel")
    ParamGEV=c(epsilonGEV,sigmaGEV,muGEV,NA, NA, NA)
    names(ParamGEV)=c("epsilonGEV","sigmaGEV","muGEV","epsilonStdErrGEV","sigmaStdErrGEV","muStdErrGEV")
    
    ParamGPD=c(epsilonGPD,sigmaGPD,thresholdGPD,NA,NA, NA,nPeaks,sampleTimeHorizon)
    names(ParamGPD)=c("epsilonGPD","sigmaGPD","thresholdGPD","epsilonStdErrGPD","sigmaStdErrGPD","thresholdStdErrGPD","nPeaks","SampleTimeHorizon")
    return(list(Fit="No fit",Params=c(ParamGEV,ParamGPD)))
  }else{
    #GEV
    # epsilonGEV <- nonStationaryEvaParams[[1]]$parameters$epsilon
    # sigmaGEV <- mean(nonStationaryEvaParams[[1]]$parameters$sigma[timeIndex])
    # muGEV <- mean(nonStationaryEvaParams[[1]]$parameters$mu[timeIndex])
    # dtSampleYears <- nonStationaryEvaParams[[1]]$parameters$timeDeltaYears
    epsilonStdErrGEV <- nonStationaryEvaParams[[1]]$paramErr$epsilonErr
    sigmaStdErrGEV <- mean(nonStationaryEvaParams[[1]]$paramErr$sigmaErr[timeIndex])
    muStdErrGEV <- mean(nonStationaryEvaParams[[1]]$paramErr$muErr[timeIndex])
    
    #GPD
    # epsilonGPD <- nonStationaryEvaParams[[2]]$parameters$epsilon
    # sigmaGPD <- mean(nonStationaryEvaParams[[2]]$parameters$sigma[timeIndex])
    # thresholdGPD <- mean(nonStationaryEvaParams[[2]]$parameters$threshold[timeIndex])
    # nPeaks <- nonStationaryEvaParams[[2]]$parameters$nPeaks
    epsilonStdErrGPD <- nonStationaryEvaParams[[2]]$paramErr$epsilonErr
    sigmaStdErrGPD <- mean(nonStationaryEvaParams[[2]]$paramErr$sigmaErr[timeIndex])
    thresholdStdErrGPD <- mean(nonStationaryEvaParams[[2]]$paramErr$thresholdErr[timeIndex])
    # thStart <- nonStationaryEvaParams[[2]]$parameters$timeHorizonStart
    # thEnd <- nonStationaryEvaParams[[2]]$parameters$timeHorizonEnd
    # sampleTimeHorizon <- as.numeric((thEnd - thStart)/365.2425)
    
    returnLevelsGEV <- tsEvaComputeReturnLevelsGEV(epsilonGEV, sigmaGEV, muGEV, epsilonStdErrGEV, sigmaStdErrGEV, muStdErrGEV, RPgoal)
    
    returnLevelsGPD <- tsEvaComputeReturnLevelsGPD(epsilonGPD, sigmaGPD, thresholdGPD, epsilonStdErrGPD, sigmaStdErrGPD, thresholdStdErrGPD,
                                                   nPeaks, sampleTimeHorizon, RPgoal)
    rlevGEV=returnLevelsGEV$returnLevels
    rlevGPD=returnLevelsGPD$returnLevels
    
    errGEV=returnLevelsGEV$returnLevelsErr
    errGPD=returnLevelsGPD$returnLevelsErr
    
    ParamGEV=c(epsilonGEV,sigmaGEV,muGEV,epsilonStdErrGEV, sigmaStdErrGEV, muStdErrGEV)
    names(ParamGEV)=c("epsilonGEV","sigmaGEV","muGEV","epsilonStdErrGEV","sigmaStdErrGEV","muStdErrGEV")
    
    ParamGPD=c(epsilonGPD,sigmaGPD,thresholdGPD,epsilonStdErrGPD,sigmaStdErrGPD, thresholdStdErrGPD,nPeaks,sampleTimeHorizon)
    names(ParamGPD)=c("epsilonGPD","sigmaGPD","thresholdGPD","epsilonStdErrGPD","sigmaStdErrGPD","thresholdStdErrGPD","nPeaks","SampleTimeHorizon")
    return(list(Fit="Fitted",ReturnLevels=c(ReturnPeriod=RPgoal, GEV=as.numeric(rlevGEV),GPD=as.numeric(rlevGPD),errGEV=as.numeric(errGEV),errGPD=as.numeric(errGPD)),Params=c(ParamGEV,ParamGPD)))
  }  
  
  
}
RPcalc<- function(params,RPiGEV,RPiGPD){
  paramx=data.frame(t(params))
  qxV=1-exp(-(1+paramx$epsilonGEV*(RPiGEV-paramx$muGEV)/paramx$sigmaGEV)^(-1/paramx$epsilonGEV))
  if (is.na(qxV)){
    returnPeriodGEV=9999
  }else{
    returnPeriodGEV=1/qxV
  }
  X0 <- paramx$nPeaks/paramx$SampleTimeHorizon
  qxD=(((1+paramx$epsilonGPD*(RPiGPD-paramx$thresholdGPD)/paramx$sigmaGPD)^(-1/paramx$epsilonGPD)))
  if (is.na(qxD)){
    returnPeriodGPD=9999
  }else{
    returnPeriodGPD=1/(X0*qxD)
  }
  return(c(GEV=returnPeriodGEV,GPD=returnPeriodGPD))
}
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
  
  return(data.frame(catchment=paramx$catchment,newRP=returnPeriods))
}

interid<- function(data,trans,WindowSize) {
  
  dt1=min(diff(data$date),na.rm=T)
  dt=as.numeric(dt1)
  tdim=attributes(dt1)$units
  if (tdim=="hours") dt=dt/24
  nRunMn = ceiling(WindowSize/dt);
  colnames(data)
  data$Q7=rollmean(data$Qs,nRunMn, align = "right", fill=NA)
  
  if (trans=="rev"){
    data$Qtrans=-(data$Q7)
  }else if(trans=="inv"){
    data$Qtrans=1/data$Q7
  }else if(trans=="lninv"){
    data$Qtrans=-ln(data$Q7)
  }
  if (length(which(is.na(data$Q7)))==length(data$Q7)){
    print("no data in this pixel")
    list0=NA
    dis07=data
    l0=NA
    fl=NA
    mindis=NA
    
  }else{
    if (length(which(is.na(data$Qs)))>0){
      print("Na alert")
      #seriefill=tsEvaFillSeries(data$date,data$Qs)
      #data$Qs=seriefill
    }
    #tinversing the discharge for EVA
    #Identification of the intermittent river
    mindis=min(data$Q7,na.rm=T)
    m0=length(which(data$Qs7==mindis))
    l0=length(which(data$Q7<=1e-4))
    #30y running average on this
    dayysbelow=tsEvaNanRunnigBlowTh(data$Q7,1e-4,4*365*30)
    dayysbelow$time=as.Date(data$date[dayysbelow$time])
    yrtot=length(unique(year(data$date)))
    list0=NA
    fl=0 #flag for intermitent river
    if (l0>=1){
      print("intermittent river")
      fl=1
      data$Qinv=NA
      data$Qlninv
      #keep days with 0 discharge
      list0=data$date[which(data$Qs==0)]
    }else if(mindis>0 & m0>=yrtot){ 
      print(paste0("river with floor low flow ",mindis))
      fl=3
      list0=data$date[which(data$Q7==mindis)]
    }
    dis07=data[,c(2,3,4)]
  }
  #The objective here is to return also the dicharge in a better format for next step of the analysis
  return(list(zerodate=list0,trdis=dis07,DaysBlow=dayysbelow,flags=c(n0d=l0,intertype=fl,mindischarge=mindis)))
}
# Load inputs from HPC computation ----------------------------------------




var = "dis"

Nsq=41
outlets="RNetwork"
#tbound=c(as.Date(paste0("1950-06-01")),as.Date(paste0("2020-06-01")))
#Impdates=seq(tbound[1],tbound[2],by="10 years")
workDir = "D:/tilloal/Documents/LFRuns_utils/data/"
setwd(workDir)

rspace= read.csv(paste0(workDir,"subspace_efas.csv"))
rspace=rspace[,-1]
nrspace=rspace[Nsq,]
print(nrspace)
if (outlets=="Hybas07"){
  outletname="outletsv8_hybas07_01min"
  nameout="UCH07"
  hydroDir<-"/BGFS/CLIMEX/tilloal/HydroMeteo/Timeseries/dis6_UnCalOut7"
  outhybas=outletopen(workDir,outletname,nrspace)
}else if (outlets=="Hybas09"){
  outletname="outlets_hybas09_01min"
  nameout="UCH09"
  hydroDir<-"/BGFS/CLIMEX/tilloal/HydroMeteo/Timeseries/dis6_UnCalOut9"
  outhybas=outletopen(workDir,outletname,nrspace)
}else if (outlets=="RNetwork"){
  outletname="efas_rnet_100km_01min"
  nameout="UCRnet"
  outhybas=outletopen(workDir,outletname,nrspace)
  #hydroDir<-paste0("/BGFS/CLIMEX/tilloal/HydroMeteo/Timeseries/dis6_",foldin)
  Idstart=as.numeric(Nsq)*100000
  if (length(outhybas$outlets)>0){
    outhybas$outlets=seq((Idstart+1),(Idstart+length(outhybas$outlets)))
  }
}
unikout=outhybas$outlets
outhybas$latlong=paste(round(outhybas$Var1,4),round(outhybas$Var2,4),sep=" ")

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

outletname="efas_rnet_100km_01min"
nameout="RNetwork"
### Hybas07 ----
Catchmentrivers7=read.csv(paste0(hydroDir,"/Catchments/from_hybas_eu_onlyid.csv"),encoding = "UTF-8", header = T, stringsAsFactors = F)
hybas07 <- read_sf(dsn = paste0(hydroDir,"/Catchments/hydrosheds/hybas_eu_lev07_v1c.shp"))
hybasf7=fortify(hybas07) 
Catamere07=inner_join(hybasf7,Catchmentrivers7,by= "HYBAS_ID")



outletname="outletsv8_hybas07_01min"
outhyb07=outletopen(hydroDir,outletname,nrspace)
catmatch=match(outhyb07$outlets,Catchmentrivers7$pointid)
mycat=Catchmentrivers7[catmatch,]

Catamere07$llcoord=paste(round(Catamere07$POINT_X,4),round(Catamere07$POINT_Y,4),sep=" ")
Catf7=inner_join(Catamere07,outhybas,by= c("llcoord"="latlong"))

st_geometry(Catf7)=NULL

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



#Load the file
#loading the files as netcdf (needs to be checked offline)
filename2=paste0("timeseries/RWCF/dis_",Nsq,"_1951_2020_rcf")
filename1=paste0("timeseries/SocCF/dis_",Nsq,"_1950_2020_cf")
df.dis$timeStamps=timeStamps
#unikout=unique(df.dis$outid)

names(df.dis)[c(1,2)]=c("dis","outlets")
#unikout=unique(df.dis$outid)


load(file=paste0(workDir,"/Drought/catchment_frost.Rdata"))
#remove first year
rmv=which(year(frostcat$time)==1950)
frostcat=frostcat[-rmv,]

tail="low"
TH1=read.csv(paste0(hydroDir,"/TrendAnalysis/trenTH_Histo_",tail,"_",Nsq,".csv"))
TH2=read.csv(paste0(hydroDir,"/TrendAnalysis/trenTH_SCF_",tail,"_",Nsq,".csv"))
TH3=inner_join(TH1,TH2,by="cid")
thresh_vec=data.frame(TH3$cid, TH3$Th_new.x)
thresh_vec$TH3.Th_new.x[which(is.na(thresh_vec$TH3.Th_new.x))]=TH3$Th_new.y[which(is.na(thresh_vec$TH3.Th_new.x))]
names(thresh_vec)=c("cid","th")

#extract discharge for a given pixel for SocCF and RWCF runs

#catchment1
catch1=4104686
c1=414686
#catchment2
catch2=4104741
c2=414741

#Run1
dists=disNcopenloc(filename1,workDir,outhybas,1)
df.dis=dists 
print(paste0("opening square ", Nsq, " /88"))
timeStamps=unique(as.Date(df.dis$time,origin="1979-01-01"))
timeStamps=as.POSIXct(timeStamps-1/24)
txx=timeStamps
print(tail(txx))
length(txx)
idfix=which(catch1==outhybas$outlets)
timeStamps=txx
dt = difftime(timeStamps[2],timeStamps[1],units="days")
dt= as.numeric(dt)

#catch=472446
thresh=thresh_vec[which(thresh_vec$cid==c2),]
trans="rev"
#seasonal split
catmat=Catf7[which(Catf7$outlets==catch2),]
Tcatmat=mycat[which(mycat$HYBAS_ID==catmat$HYBAS_ID),]
Tcatchment=which(colnames(frostcat)==Tcatmat$pointid)

df.disX=disNcopenloc(filename1,hydroDir,outhybas,idfix)
series=df.disX$outlets
print(mean(series))
#df.q7=tsEvaNanRunningMean(series, windowSize =7)
data=data.frame(txx,series)
names(data)=c("date","Qs")
rmv=which(as.integer(format(data$date, "%Y"))==1950)
if (length(rmv)>0){
  data=data[-rmv,]
}
intermit=interid(data,trans,WindowSize=7)
interflag=intermit$flags[2]
serie1=data.frame(data$date,intermit$trdis$Q7)
#series=data.frame(txx,intermit$invdis$Qinv)
#remove frost timesteps, this can be modified to do the anlysis only on frost moments
if (length(Tcatchment)>0){
  frostserie=data.frame(frostcat[,1],frostcat[,Tcatchment])
  frosttime=which(frostserie[,2]<0)
}else{
  frosttime=NA
}
ciPercentile=80
minPeakDistanceInDays=10


names(serie1)=c("timestamp","dis")
dt1=min(diff(serie1$timestamp),na.rm=T)
dt=as.numeric(dt1)
tdim=attributes(dt1)$units
if (tdim=="hours") dt=dt/24
if (dt==1){
  timeDays=serie1$timestamp
}else{
  timeDays=unique(as.Date(serie1$timestamp))
}


bounds=c(year(timeDays[1]),year(timeDays[length(timeDays)]))
realbound=bounds
tbound=c(as.Date(paste0(realbound[1],"-12-31")),as.Date(paste0(realbound[2],"-12-31")))
Impdates=seq(tbound[1],tbound[2],by="1 years")
datex=yday(timeDays)
dtect=c(diff(datex),-1)
last_days <- timeDays[which(dtect<0)]
tindexes=match(last_days,timeDays)


# Compute return periods and levels
nv=length(unique(serie1$dis))

if (length(which(is.na(serie1$dis)))>0){
  print("Na alert")
  seriefill=tsEvaFillSeries(serie1$timestamp,serie1$dis)
  serie1$dis=seriefill
}
timeAndSeriesSCF=serie1
names(timeAndSeriesSCF)=c("timestamp","data")
#timeAndSeries$data=jitter(timeAndSeries$data)

if (haz=="drought" & length(!is.na(frosttime))>1){
  if (season=="nonfrost"){
    print("nonfrost season")
    timeAndSeriesSCF$data[frosttime]=NA
  }else if (season=="frost"){
    print("frost season")
    timeAndSeriesSCF$data[-frosttime]=NA
  }else if (season=="year"){
    print("no seasonal divide")
  }else {print("season must be frost or nonfrost")}
}else{
  print("no frost season for this river")
}
#I remove the first two months because they sometimes have weird values
tsm=1/dt

rmv=which(as.integer(format(timeAndSeriesSCF$timestamp, "%Y"))==1950)
if (length(rmv)>0){
  timeAndSeriesSCF=timeAndSeriesSCF[-rmv,]
}
# rmv2=which(is.na(timeAndSeriesSCF$data))
# timeAndSeriesSCF=timeAndSeriesSCF[-rmv2,]
# series=timeAndSeriesSCF[,2]

timeWindow = 365.25*30; #time windows in days, the correction is done within the functions
windowSize=366
minPeakDistanceInDays=30
trans="rev"
timeStamps=timeAndSeriesSCF$timestamp
cat(paste0("\nsquare: ", Nsq, " pixel: ",idfix,"/",endid))
Nonstat<-TsEvaNs(timeAndSeriesSCF, timeWindow, transfType='trendPeaks',
                 ciPercentile = ciPercentile, minPeakDistanceInDays = minPeakDistanceInDays,
                 lowdt=7,trans=trans,tail = tail, TrendTh = thresh$th)

nonStationaryEvaParamsSCF=Nonstat[[1]]
stationaryTransformDataSCF=Nonstat[[2]]





#Run2

dists=disNcopenloc(filename2,workDir,outhybas,1)
df.dis=dists 
print(paste0("opening square ", Nsq, " /88"))
timeStamps=unique(as.Date(df.dis$time,origin="1979-01-01"))
timeStamps=as.POSIXct(timeStamps-1/24)
txx=timeStamps
print(tail(txx))

idfix=which(catch1==outhybas$outlets)
timeStamps=txx
dt = difftime(timeStamps[2],timeStamps[1],units="days")
dt= as.numeric(dt)

#catch=472446
thresh=thresh_vec[which(thresh_vec$cid==c2),]
trans="rev"
#seasonal split
catmat=Catf7[which(Catf7$outlets==catch2),]
Tcatmat=mycat[which(mycat$HYBAS_ID==catmat$HYBAS_ID),]
Tcatchment=which(colnames(frostcat)==Tcatmat$pointid)

df.disX=disNcopenloc(filename2,hydroDir,outhybas,idfix)
series=df.disX$outlets
print(mean(series))
#df.q7=tsEvaNanRunningMean(series, windowSize =7)

data=data.frame(txx,series)
names(data)=c("date","Qs")
intermit=interid(data,trans,WindowSize=7)
interflag=intermit$flags[2]
series=data.frame(txx,intermit$trdis$Q7)
#series=data.frame(txx,intermit$invdis$Qinv)
#remove frost timesteps, this can be modified to do the anlysis only on frost moments
if (length(Tcatchment)>0){
  frostserie=data.frame(frostcat[,1],frostcat[,Tcatchment])
  frosttime=which(frostserie[,2]<0)
}else{
  frosttime=NA
}
ciPercentile=80
minPeakDistanceInDays=30




names(series)=c("timestamp","dis")
dt1=min(diff(series$timestamp),na.rm=T)
dt=as.numeric(dt1)
tdim=attributes(dt1)$units
if (tdim=="hours") dt=dt/24
if (dt==1){
  timeDays=serie1$timestamp
}else{
  timeDays=unique(as.Date(series$timestamp))
}


bounds=c(year(timeDays[1]),year(timeDays[length(timeDays)]))
realbound=bounds
tbound=c(as.Date(paste0(realbound[1],"-12-31")),as.Date(paste0(realbound[2],"-12-31")))
Impdates=seq(tbound[1],tbound[2],by="1 years")
datex=yday(timeDays)
dtect=c(diff(datex),-1)
last_days <- timeDays[which(dtect<0)]
tindexes=match(last_days,timeDays)


# Compute return periods and levels
plot(series)
nv=length(unique(series$dis))

if (length(which(is.na(series$dis)))>0){
  print("Na alert")
  seriefill=tsEvaFillSeries(series$timestamp,series$dis)
  series$dis=seriefill
}
timeAndSeriesRWCF=series
names(timeAndSeriesRWCF)=c("timestamp","data")
#timeAndSeries$data=jitter(timeAndSeries$data)

if (haz=="drought" & length(!is.na(frosttime))>1){
  if (season=="nonfrost"){
    print("nonfrost season")
    timeAndSeriesRWCF$data[frosttime]=NA
  }else if (season=="frost"){
    print("frost season")
    timeAndSeriesRWCF$data[-frosttime]=NA
  }else if (season=="year"){
    print("no seasonal divide")
  }else {print("season must be frost or nonfrost")}
}else{
  print("no frost season for this river")
}
#I remove the first two months because they sometimes have weird values
tsm=1/dt

rmv=which(as.integer(format(timeAndSeriesRWCF$timestamp, "%Y"))==1950)
rmv=which(year(timeStamps)==1950)
if (length(rmv)>0){
  timeAndSeriesRWCF=timeAndSeriesSCF[-rmv,]
}
# rmv2=which(is.na(timeAndSeriesRWCF$data))
# timeAndSeriesRWCF=timeAndSeriesRWCF[-rmv2,]
#series=timeAndSeriesRWCF[,2]



rm(tsEvaSampleData)
tsEvaSampleData <- function(ms, meanEventsPerYear,minEventsPerYear, minPeakDistanceInDays,tail=NA) {
  
  pctsDesired = c(90, 95, 99, 99.9)
  args <- list(meanEventsPerYear = meanEventsPerYear,
               minEventsPerYear = minEventsPerYear,
               potPercentiles = c(seq(70,90,by=1), seq(91,99.5,by=0.5)))
  meanEventsPerYear = args$meanEventsPerYear
  minEventsPerYear = args$minEventsPerYear
  potPercentiles = args$potPercentiles
  if(is.na(tail)) stop("tail for POT selection needs to be 'high' or 'low'")
  
  POTData <- tsGetPOT(ms, potPercentiles, meanEventsPerYear,minEventsPerYear,minPeakDistanceInDays, tail)
  
  vals <- quantile(ms[,2], pctsDesired/100,na.rm=T)
  percentiles <- list(precentiles = pctsDesired, values = vals)
  
  pointData <- list()
  pointData$completeSeries <- ms
  pointData$POT <- POTData
  pointDataA <- computeAnnualMaxima(ms)
  pointDataM <- computeMonthlyMaxima(ms)
  
  yrs <- unique(as.numeric(format(as.Date(ms[,1]+3600), "%Y")))
  yrs <- yrs - min(yrs)
  pointData$years <- seq(min(yrs),max(yrs),1)
  
  pointData$Percentiles <- percentiles
  pointData$annualMax=pointDataA$annualMax
  pointData$annualMaxDate=pointDataA$annualMaxDate
  pointData$annualMaxIndx=pointDataA$annualMaxIndx
  pointData$monthlyMax=pointDataM$monthlyMax
  pointData$monthlyMaxDate=pointDataM$monthlyMaxDate
  pointData$monthlyMaxIndx=pointDataM$monthlyMaxIndx
  
  return(pointData)
}



tsGetPOT <- function(ms, pcts, desiredEventsPerYear,minEventsPerYear, minPeakDistanceInDays, tail) {
  
  if (minPeakDistanceInDays == -1) {
    stop("label parameter 'minPeakDistanceInDays' must be set")
  }
  dt1=min(diff(ms[,1]),na.rm=T)
  dt=as.numeric(dt1)
  tdim=attributes(dt1)$units
  if (tdim=="hours") dt=dt/24
  if (tdim=="seconds") dt=dt/3600
  minPeakDistance <- minPeakDistanceInDays/dt
  minRunDistance <- minPeakDistance
  nyears <- round(as.numeric((max(ms[,1]) - min(ms[,1]))/365.25))
  if (length(pcts) == 1) {
    pcts = c(pcts - 3, pcts)
    desiredEventsPerYear = -1
  }
  
  numperyear <- rep(NA, length(pcts))
  minnumperyear <- rep(NA, length(pcts))
  thrsdts <- rep(NA, length(pcts))
  gpp=rep(NA, length(pcts))
  devpp=rep(NA, length(pcts))
  dej=0
  skip=0
  trip=NA
  perfpen=0
  for (ipp in 1:length(pcts)) {
    #Skip is used to prevent finding peaks for unappropriate thresholds
    if (skip>0) {
      skip=skip-1
    }else{
      if(dej==0){
        thrsdt <- quantile(ms[,2],pcts[ipp]/100,na.rm=T)
        thrsdts[ipp] <- thrsdt
        ms[,2][which(is.na(ms[,2]))]=-9999
        minEventsPerYear=1
        
        if(tail=="high") {
          #boundaries of shape parameter
          shape_bnd=c(-0.5,1)
          pks <- pracma::findpeaks(ms[,2],minpeakdistance = minPeakDistance, minpeakheight = thrsdt)
        }
        if(tail=="low") {
          pks <- declustpeaks(data = ms[,2] ,minpeakdistance = minPeakDistance ,minrundistance = minRunDistance, qt=thrsdt)
          shape_bnd=c(-2,0)
        }
        numperyear[ipp] <- length(pks[,1])/nyears
        #print(numperyear[ipp])
        if(numperyear[ipp]>=3*desiredEventsPerYear & ipp<(length(pcts)-5)) skip = floor(length(pcts)/8)
        if(numperyear[ipp]<0.9*minEventsPerYear) {
          perfpen=(pcts[ipp])*100
        }
        if(numperyear[ipp]<(0.7*minEventsPerYear)) {
          perfpen=(pcts[ipp])*1000
        }
        if(numperyear[ipp]<=desiredEventsPerYear+1 & dej==0){
          fgpd=suppressWarnings(try(POT::fitgpd(pks[,1], threshold = thrsdt, est = "mle",method="BFGS",std.err.type = "expected")))
          if(inherits(fgpd, "try-error")){
            gpdpar=9999
            deviance=9999
            devpp[ipp]=1e9
            gpp[ipp]=9999
          }else {
            gpdpar=fgpd$fitted.values
            deviance=fgpd$deviance
            devpp[ipp]=AIC(fgpd)+perfpen
            gpp[ipp]=gpdpar[2]
          }
          nperYear <- tsGetNumberPerYear(ms, pks[,2])
          minnumperyear[ipp] <- min(nperYear$Freq, na.rm = TRUE)
        }
      }
    }
  }
  md= abs(pks[1,1]-pks[2,1])
  devpp[1]=NA
  if(is.na(trip)){
    isok=F
    devpx=devpp
    count=length(devpx)
    while(isok==F){
      #safety measure for stability of parameter
      dshap=c(0,diff(gpp))
      #Penalizing fits with positive shape parameters for low tail
      if(tail=="low") {
        #for very bounded distributions
        if (md<0.1){
          devpp[which(gpp>=-0.5)]=devpp[which(gpp>=-0.5)]+9999
        }else{
          devpp[which(gpp>=0)]=devpp[which(gpp>=0)]+9999
        }

      }
      devpp[which(abs(dshap)>0.5)]=devpp[which(abs(dshap)>0.5)]+99999
      trip=which.min(devpp)
      #message(paste0("shape outside boudaries: ",round(gpp[trip],2)))
      #isok=T
      #trip=which.min(devpx)
      isok=dplyr::between(round(gpp[trip],1), shape_bnd[1], shape_bnd[2])
      count=count+1
      if(isok==F)devpx[trip]=devpx[trip]+9999
      # if(count>(length(devpx)-1)){
      #   #safety measure for stability of parameter
      #   dshap=c(0,diff(gpp))
      #   plot(pcts,dshap)
      #   #Penalizing fits with positive shape parameters for low tail
      #   if(tail=="low") devpp[which(gpp>=0)]=devpp[which(gpp>=0)]+99999
      #   devpp[which(abs(dshap)>0.5)]=devpp[which(abs(dshap)>0.5)]+9999
      #   trip=which.min(devpp)
      #   message(paste0("shape outside boudaries: ",round(gpp[trip],2)))
      #   isok=T
      # }
    }
  }
  # plot(pcts,devpp,ylim=c(0,1e5))
  # plot(pcts,gpp)
  # print(devpp)
  # print(pcts)
  message(paste0("\nmax threshold is: ", pcts[trip],"%"))
  message(paste0("\naverage number of events per year = ",round(numperyear[trip],1) ))
  
  diffNPerYear <- mean(diff(na.omit(rev(numperyear)), na.rm = TRUE))
  if (diffNPerYear == 0) diffNPerYear <- 1
  diffNPerYear <- 1
  thresholdError <- -mean(diff(na.omit(thrsdts))/diffNPerYear)/2
  indexp <- trip
  if (!is.na(indexp)) {
    thrsd <- quantile(ms[,2],pcts[indexp]/100)
    pct <- pcts[indexp]
  } else {
    thrsd <- 0
    pct
  }
  # Find peaks in the second column of the matrix 'ms'
  if(tail=="high") pks_and_locs <- pracma::findpeaks(ms[,2],minpeakdistance = minPeakDistance, minpeakheight = thrsd)
  if(tail=="low") pks_and_locs <- declustpeaks(data = ms[,2] ,minpeakdistance = minPeakDistance ,minrundistance = minRunDistance, qt=thrsd)
  
  # Assign peaks and peak locations to separate variables
  pks <- pks_and_locs[,1]
  locs <- pks_and_locs[,2]
  st<-pks_and_locs[,3]
  end=pks_and_locs[,4]
  # Create a list to store results
  POTdata <- list()
  # Assign values to the fields of the list
  POTdata[['threshold']] <- thrsd
  POTdata[['thresholdError']] <- thresholdError
  POTdata[['percentile']] <- pct
  POTdata[['peaks']] <- pks
  POTdata[['stpeaks']] <- st
  POTdata[['endpeaks']] <- end
  POTdata[['ipeaks']] <- locs
  POTdata[['time']] <- ms[locs, 1]
  POTdata[['pars']] <- gpdpar
  
  
  return(POTdata)
}


timeWindow = 365.25*30; #time windows in days, the correction is done within the functions
windowSize=366
trans="rev"
minPeakDistanceInDays=30
timeStamps=timeAndSeriesRWCF$timestamp
cat(paste0("\nsquare: ", Nsq, " pixel: ",idfix,"/",endid))
Nonstat<-TsEvaNs(timeAndSeriesRWCF, timeWindow, transfType='trendPeaks',
                 ciPercentile = ciPercentile, minPeakDistanceInDays = minPeakDistanceInDays,
                 lowdt=7,trans=trans,tail = tail, TrendTh = thresh$th)

nonStationaryEvaParamsRWCF=Nonstat[[1]]
stationaryTransformDataRWCF=Nonstat[[2]]



ExRange= c(min(nonStationaryEvaParamsRWCF$potObj$parameters$peaks),max(nonStationaryEvaParamsRWCF$potObj$parameters$peaks))
haz="drought"
if (haz=="flood") wr2 <- c(seq(0.8*min(ExRange),1.2*max(ExRange),length.out=300))
if (haz=="drought") wr2 <- c(seq(1.2*min(ExRange),0.2*max(ExRange),length.out=300))

Plot1= tsEvaPlotGPDImageScFromAnalysisObj(wr2, nonStationaryEvaParamsRWCF, stationaryTransformDataRWCF,trans=trans)

Plot1

timeIndex=10
Plot2 = tsEvaPlotReturnLevelsGPDFromAnalysisObj(nonStationaryEvaParamsRWCF, stationaryTransformDataRWCF, timeIndex, trans=trans,ylabel="Discharge (m3/s)")
Plot2

Plot2 = tsEvaPlotReturnLevelsGEVFromAnalysisObj(nonStationaryEvaParamsSCF, stationaryTransformDataSCF, timeIndex, trans=trans,ylabel="Discharge (m3/s)")
Plot2

stationaryTransformDataSCF$timeStampsDay=unique(as.Date(stationaryTransformDataSCF$timeStamps))
pikosS = data.frame(nonStationaryEvaParamsSCF$potObj$parameters$peaks,
                    nonStationaryEvaParamsSCF$potObj$parameters$peakID,
                    nonStationaryEvaParamsSCF$potObj$parameters$peakST,
                    nonStationaryEvaParamsSCF$potObj$parameters$peakEN)
names(pikosS) = c("value", "timeID", "tIDstart", "tIDend")
pikosS$time = stationaryTransformDataSCF$timeStamps[pikosS$timeID]
pikosS$catch = rep(catch1, length(pikosS[,1]))

stationaryTransformDataRWCF$timeStampsDay = unique(as.Date(stationaryTransformDataRWCF$timeStamps))
pikosR = data.frame(nonStationaryEvaParamsRWCF$potObj$parameters$peaks,
                    nonStationaryEvaParamsRWCF$potObj$parameters$peakID,
                    nonStationaryEvaParamsRWCF$potObj$parameters$peakST,
                    nonStationaryEvaParamsRWCF$potObj$parameters$peakEN)
names(pikosR) = c("value", "timeID", "tIDstart", "tIDend")
pikosR$time = stationaryTransformDataRWCF$timeStamps[pikosR$timeID]
pikosR$catch = rep(catch1, length(pikosR[,1]))


merdax=(timeAndSeriesSCF$data-timeAndSeriesRWCF$data)
merdax[c(20:100)]
plot(pikosS$timeID,pikosS$value,pch=15)
points(pikosR$timeID,pikosR$value,col=2)




#Here I need to convert the timeStamp to a daily one if dt is not 1
dt1=min(diff(timeStamps),na.rm=T)
dt=as.numeric(dt1)
tdim=attributes(dt1)$units
if (tdim=="hours") dt=dt/24
if (dt==1){
  timeDays=stationaryTransformDataSCF$timeStamps
}else{
  timeDays=stationaryTransformDataSCF$timeStampsDay
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
timeIndex=tindexes[7]
RLevs100R=ComputeReturnLevels(nonStationaryEvaParamsRWCF, RPgoal, timeIndex)
RLevs100S=ComputeReturnLevels(nonStationaryEvaParamsSCF, RPgoal, timeIndex)

RLevs100R
RLevs100S

paramsS=c(catch1,year(Impdates[1]),timeIndex,RLevs100S$Params)
names(paramsS)[1:3]=c("catchment","Year","timeIndex")

RLgevS=RLevs100S$ReturnLevels[2]
RLgpdS=RLevs100S$ReturnLevels[3]
ERgevS=RLevs100S$ReturnLevels[4]
ERgpdS=RLevs100S$ReturnLevels[5]
nRPgevS=nRPgpdS=10
paramsS=c()
for (t in 2:length(Impdates)){
  timeIndex=tindexes[t]
  RLevs100i=ComputeReturnLevels(nonStationaryEvaParamsSCF, RPgoal, timeIndex)
  params=c(catch1,year(Impdates[t]),timeIndex,RLevs100i$Params)
  names(params)[1:3]=c("catchment","Year","timeIndex")
  
  Rper=RPcalc(params,RPiGEV=RLevs100S$ReturnLevels[2],RPiGPD=RLevs100S$ReturnLevels[3])
  nRPgpdS=c(nRPgpdS,Rper[2])
  nRPgevS=c(nRPgevS,Rper[1])
  RLgevS=cbind(RLgevS,RLevs100i$ReturnLevels[2])
  RLgpdS=cbind(RLgpdS,RLevs100i$ReturnLevels[3])
  ERgevS=cbind(ERgevS,RLevs100i$ReturnLevels[4])
  ERgpdS=cbind(ERgpdS,RLevs100i$ReturnLevels[5])
  # if (length(parlist)>1) colnames(parlist)=names(params)
  # parlist=rbind(parlist,params)
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




RLgevR=RLevs100R$ReturnLevels[2]
RLgpdR=RLevs100R$ReturnLevels[3]
ERgevR=RLevs100R$ReturnLevels[4]
ERgpdR=RLevs100R$ReturnLevels[5]
nRPgevR=nRPgpdR=10
paramsR=c()
for (t in 2:length(Impdates)){
  timeIndex=tindexes[t]
  RLevs100i=ComputeReturnLevels(nonStationaryEvaParamsRWCF, RPgoal, timeIndex)
  params=c(catch1,year(Impdates[t]),timeIndex,RLevs100i$Params)
  names(params)[1:3]=c("catchment","Year","timeIndex")
  
  Rper=RPcalc(params,RPiGEV=RLevs100R$ReturnLevels[2],RPiGPD=RLevs100R$ReturnLevels[3])
  nRPgpdR=c(nRPgpdR,Rper[2])
  nRPgevR=c(nRPgevR,Rper[1])
  RLgevR=cbind(RLgevR,RLevs100i$ReturnLevels[2])
  RLgpdR=cbind(RLgpdR,RLevs100i$ReturnLevels[3])
  ERgevR=cbind(ERgevR,RLevs100i$ReturnLevels[4])
  ERgpdR=cbind(ERgpdR,RLevs100i$ReturnLevels[5])
  # if (length(parlist)>1) colnames(parlist)=names(params)
  # parlist=rbind(parlist,params)
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


plot(nonStationaryEvaParamsSCF$gevObj$parameters$annualMaxIndx,nonStationaryEvaParamsSCF$gevObj$parameters$annualMax,type="p",pch=16)
lines(nonStationaryEvaParamsRWCF$gevObj$parameters$annualMaxIndx,nonStationaryEvaParamsRWCF$gevObj$parameters$annualMax,type="p",col="blue",pch=16)
#lines(nonStationaryEvaParamsS$gevObj$parameters$annualMaxIndx,nonStationaryEvaParamsS$gevObj$parameters$annualMax,type="p",col="red",pch=16)

lines(parlistS$timeIndex , t(RLgpdS))
lines(parlistS$timeIndex,t(RLgpdR),col="blue")
lines(parlistS$timeIndex,t(RLgpdS),col="red")


RLGPDR=as.vector(-t(RLgpdR))
RLGPDS=as.vector(-t(RLgpdS))

RLGPDR[which(RLGPDR<0)]=0
RLGPDS[which(RLGPDS<0)]=0
plot(as.numeric(RLGPDS))
points(as.numeric(RLGPDR),col=4)

#Socioeconomic changes
DSoc=(RLGPDR-RLGPDS)/mean(RLGPDS)*100
plot(DSoc)

#difference attributed to reservoir konstraktion
#Reservoir changes
DRes=(RLGPDH-RLGPDR)/RLGPDH[1]*100
plot(DRes)



#Climate change
DClim=(RLGPDS-RLGPDS[1])


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






