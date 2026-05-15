##########################################################################################
############   THIS SCRPT IS FOR RUNNING 1 SQUARE OF THE HERA DOMAIN ON A HPC ############
##########################################################################################

setwd(dirname(rstudioapi::getSourceEditorContext()$path))
getwd()
source("functions_trends.R")
suppressWarnings(suppressMessages(library(ncdf4)))
suppressWarnings(suppressMessages(library(sf)))
suppressWarnings(suppressMessages(library(rnaturalearth)))
suppressWarnings(suppressMessages(library(rnaturalearthdata)))
suppressWarnings(suppressMessages(library(rgeos)))
suppressWarnings(suppressMessages(library(dplyr)))
suppressWarnings(suppressMessages(library(fs)))
suppressWarnings(suppressMessages(library(tsibble)))
suppressWarnings(suppressMessages(library(ggplot2)))
suppressWarnings(suppressMessages(library(scales)))
suppressWarnings(suppressMessages(library(pracma)))
suppressWarnings(suppressMessages(library(lubridate)))
suppressWarnings(suppressMessages(library(xts)))
# suppressWarnings(suppressMessages(library(texmex)))
# suppressWarnings(suppressMessages(library(ggplot2)))
# suppressWarnings(suppressMessages(library(moments)))
# suppressWarnings(suppressMessages(library(dplyr)))
# suppressWarnings(suppressMessages(library(changepoint)))
# suppressWarnings(suppressMessages(library(iemisc)))
suppressWarnings(suppressMessages(library(evd)))
suppressWarnings(suppressMessages(library(POT)))
suppressWarnings(suppressMessages(library(RtsEva)))


###########################  FUNCTIONS   #################################################


#New peak tremd function

tsEvaTransformSeriesToStationaryMMXTrend<-function (timeStamps, series, timeWindow) 
{
  tserie=data.frame(timeStamps,series)
  monthly_max=tserie %>% 
    mutate(month = floor_date(timeStamps, "month")) %>% 
    group_by(month) %>% 
    summarise(max_value = max(series))
  
  
  serieb <- series
  tm=na.omit(match(as.Date(as.character(monthly_max$month)),as.Date(timeStamps)))
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

#Experimental TSEVA function

TsEvaNs<- function(timeAndSeries, timeWindow, transfType='trendPeaks',minPeakDistanceInDays=10,
                   seasonalityVar=NA,minEventsPerYear=-1, gevMaxima='annual',
                   ciPercentile=90, gevType = 'GEV', evdType = c('GEV', 'GPD'),
                   tail="high", epy=-1, lowdt=7, trans=NULL, TrendTh=NA,shape_bnd=NA){
  
  
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
  
  
  minEventsPerYear = 2
  if (tail=="low"){
    #default 7-day flow for low flow, can be modified by user
    start_index=1
    indices_to_extract <- seq(from = start_index, to = length(series), by = lowdt/dt)
    series=series[indices_to_extract]
    timeStamps=timeStamps[indices_to_extract]
    if (is.na(shape_bnd)){
      shape_bnd=c(-1,0)
    }
    minEventsPerYear = 1
    if (trans=="rev"){
      series=-1*series
    }else if(trans=="inv"){
      series=1/series
    }else if (trans=="lninv"){
      series=-log(series)
    }
  }else{  # default shape parameter bounds
    if (is.na(shape_bnd)){
      shape_bnd=c(-0.5,1)
    }
  }
  
  if (transfType == 'trend'){
    message('\nevaluating long term variations of extremes')
    trasfData = tsEvaTransformSeriesToStationaryTrendOnly(timeStamps, series, timeWindow)
    gevMaxima = 'annual';
    potEventsPerYear = epy;
    
    
  }else  if (transfType == 'trendChange'){
    message('\nevaluating long term variations of extremes and change point detection')
    trasfData = tsEvaTransformSeriesToStationaryTrendAndChangepts(timeStamps, series, timeWindow)
    gevMaxima = 'annual';
    potEventsPerYear = epy;
    
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
    
  }else if (transfType == 'trendPeaks') {
    print(TrendTh)
    message(paste0('\nevaluating long term variations of the peaks'))
    if (is.na(TrendTh)){
      TrendTh=try(tsEvaFindTrendThreshold(series, timeStamps, timeWindow),T)
      if (inherits(TrendTh, "try-error") || length(TrendTh) == 0 || all(is.na(TrendTh))) {
        TrendTh <- 0.1
      }
      trasfData = tsEvaTransformSeriesToStationaryPeakTrend( timeStamps, series, timeWindow, TrendTh);
    }else{
      if (TrendTh=="MMX"){
        trasfData = tsEvaTransformSeriesToStationaryMMXTrend( timeStamps, series, timeWindow);
        print("using MMX trend")
      }else{
        trasfData = tsEvaTransformSeriesToStationaryPeakTrend( timeStamps, series, timeWindow, TrendTh);
        c=0
        while ((length(unique(trasfData$trendSeries))<2)) {
          trasfData = tsEvaTransformSeriesToStationaryPeakTrend( timeStamps, series, timeWindow, TrendTh=TrendTh-c);
          c=c+0.1
        }
      }
    }
    gevMaxima = 'annual'
    potEventsPerYear = epy
    
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
  pointData = tsEvaSampleDataX(ms, potEventsPerYear, minEventsPerYear, minPeakDistanceInDays,tail, shape_bnd);
  evaAlphaCI = .68; # in a gaussian approximation alphaCI~68% corresponds to 1 sigma confidence
  eva = tsEVstatisticsX(pointData, evaAlphaCI, gevMaxima, gevType, evdType,shape_bnd);
  
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

tsEVstatisticsX <- function(pointData, alphaCI = 0.95, gevMaxima = 'annual', gevType = 'GEV', evdType = c('GEV', 'GPD'),shape_bnd=c(-0.5,1)) {
  # Create empty data structures
  EVmeta <- list()
  EVdata <- list()
  isValid <- TRUE
  
  minGEVSample <- 7
  if(is.null(alphaCI)){alphaCI <- .95}
  if(is.null(gevMaxima)){gevMaxima <- 'annual'}
  if(is.null(gevType)){gevType <- 'GEV'}
  if(is.null(evdType)){evdType <- c('GEV', 'GPD')}
  
  # Define Tr vector
  Tr <- c(5,10,20,50,100,200,500,1000)
  EVmeta$Tr <- Tr
  nyears <- length(pointData$annualMax)[1]
  
  imethod <- 1
  methodname <- 'GEVstat'
  paramEsts <- numeric(3)
  paramCIs <- matrix(NA, nrow = 2, ncol = 3)
  rlvls <- numeric(length(Tr))
  
  # GEV statistics
  if (('GEV' %in% evdType) && !is.null(pointData$annualMax)) {
    if (gevMaxima == 'annual') {
      tmpmat <- pointData$annualMax
    } else if (gevMaxima == 'monthly') {
      tmpmat <- pointData$monthlyMax
    } else {
      stop(paste0('Invalid gevMaxima type: ', gevMaxima))
    }
    iIN <- length(tmpmat)
    if (sum(iIN) >= minGEVSample) {
      tmp <- data.frame(yr=year(pointData$annualMaxDate),dt=tmpmat)
      # Perform GEV/Gumbel fitting and computation of return levels
      if (gevType == "GEV"){
        stdfit=TRUE
        #try to fit GEV with bounded shape parameters and stderr, reduces the constrains if no fit.
        fit <- suppressWarnings(try(evd::fgev(x=tmp$dt,method="L-BFGS-B",lower=c(-Inf,-Inf,shape_bnd[1]),upper=c(Inf,Inf,shape_bnd[2]),std.err = T),TRUE))
        if(inherits(fit, "try-error")){
          stdfit=FALSE
          message("Not able to fit GEV stderr")
          fit <- suppressWarnings(try(evd::fgev(x=tmp$dt,method="L-BFGS-B",lower=c(-Inf,-Inf,shape_bnd[1]),upper=c(Inf,Inf,shape_bnd[2]),std.err = F),TRUE))
        }
        if(inherits(fit, "try-error")){
          stdfit=FALSE
          message("Not able to fit GEV with constrained parameters")
          print(tmp)
          fit <- suppressWarnings(try(evd::fgev(x=tmp$dt,std.err = F),TRUE))
          print(fit)
        }
        paramEsts <- c(mu=fit$par[1],sigma=fit$par[2],xi=fit$par[3])
        alphaCIx=1-alphaCI
        if (stdfit==TRUE){
          probs <- c(alphaCIx/2, 1-alphaCIx/2)
          # Compute the CI for k using a normal distribution for khat.
          kci <- try(stats::qnorm(probs, paramEsts[3], fit$std.err[3]),TRUE)
          kci[kci < -1] <- -1
          # Compute the CI for sigma using a normal approximation for log(sigmahat)
          # and transform back to the original scale.
          lnsigci <- try(stats::qnorm(probs, log(paramEsts[2]), fit$std.err[2]/paramEsts[2]),silent=T)
          muci <- stats::qnorm(probs, paramEsts[1], fit$std.err[1])
          paramCIs <- cbind(muci=muci, sigci=exp(lnsigci), kci=kci)
        }else{
          #not able to generate parameters CI
          paramCIs <- cbind(muci=NA, sigci=NA, kci=NA)
        }
      } else if (gevType == "gumbel"){
        fit <- texmex::evm(y = .data$dt, data = tmp, family = gumbel)
        paramEsts <- c(fit$par[1], exp(fit$par[2]),0)
        
        alphaCIx=1-alphaCI
        probs <- c(alphaCIx/2, 1-alphaCIx/2)
        # Compute the CI for k using a normal distribution for khat.
        if(is.character(fit$se[1])){
          message("Gumbel fitted")
          #methodname <- 'No fit'
          gevType = "gumbel"
        }else{
          message("Gumbel fitted")
          #gevType = "gumbel"
          kci <- c(NA,NA)
          # Compute the CI for sigma using a normal approximation for log(sigmahat)
          # and transform back to the original scale.
          lnsigci <- try(stats::qnorm(probs, log(paramEsts[2]), fit$se[2]))
          
          muci <- stats::qnorm(probs, paramEsts[1], fit$se[1])
          paramCIs <- cbind(muci=muci, sigci=exp(lnsigci), kci=kci)
        }
      }
      else {
        stop("tsEVstatistics: invalid gevType: ", gevType, ". Can be only GEV, or Gumbel")
      }
      rlvls <- qgev(1-1/Tr, paramEsts[1], paramEsts[2], paramEsts[3])
      
    }else{
      message("Could not fit GEV")
      methodname <- 'No fit'
      rlvls=NA
      paramEsts=NA
      paramCIs=NA
      
    }
  }
  
  EVdata$GEVstat <- list(method=methodname,
                         values=rlvls,
                         parameters=paramEsts,
                         paramCIs = paramCIs)
  
  # Create output structures for GEV statistics
  # GPD statistics
  imethod <- 2
  methodname <- 'GPDstat'
  paramEsts <- numeric(6)
  paramCIs <- matrix(NA, nrow = 2, ncol = 3)
  rlvls <- numeric(length(Tr))
  
  if (('GPD' %in% evdType) && !is.null(pointData$annualMax)) {
    # Perform GPD fitting and computation of return levels
    message("Fitted GPD")
    ik <- 1
    th=pointData$POT$threshold
    d1 <- pointData$POT$peaks
    
    fit=pointData$POT$pars
    # fit <- suppressWarnings(try(POT::fitgpd(pks[,1], threshold = thrsdt, est = "mle",
    #                                         method = "L-BFGS-B", scale= pointData$POT$pars[1], shape= pointData$POT$pars[2]),TRUE))
    if(!inherits(fit, "try-error")){
      ksi <- fit$param[2]
      sgm <- fit$param[1]
      fit$std.err
      alphaCIx=1-alphaCI
      probs <- c(alphaCIx/2, 1-alphaCIx/2)
      kci <- try(stats::qnorm(probs, ksi, fit$std.err[2]),silent=T)
      kci[kci < -1] <- -1
      # Compute the CI for sigma using a normal approximation for log(sigmahat)
      # and transform back to the original scale.
      lnsigci <- try(stats::qnorm(probs, log(sgm), fit$std.err[1]/sgm))
      paramCIs <- cbind(kci, sigci=exp(lnsigci))
      
      # Create output structures for GPD statistics
      
      # Assign values to paramEstsall
      paramEstsall <- c(sgm, ksi, pointData$POT$threshold, length(d1), length(pointData$POT$peaks), pointData$POT$percentile)
      # Assign values to rlvls
      rlvls <- pointData$POT$threshold + (sgm/ksi) * ((((length(d1)/length(pointData$POT$peaks))*(1/Tr))^(-ksi))-1)
      
      EVdata$GPDstat <- list(method=methodname,
                             values=rlvls,
                             parameters=paramEstsall,
                             paramCIs = paramCIs)
    }else {
      methodname <- 'No fit'
      ik <- 1
      th=pointData$POT$threshold
      d1 <- pointData$POT$peaks
      paramEstsall <- c(pointData$POT$pars[1], pointData$POT$pars[2],
                        pointData$POT$threshold, length(d1),
                        length(pointData$POT$peaks), pointData$POT$percentile)
      EVdata$GPDstat <-  list(method=methodname,
                              values=NA,
                              parameters=paramEstsall,
                              paramCIs = NA)
      message("could not estimate GPD: bounded distribution")
    }
  } else {
    methodname <- 'No fit'
    ik <- 1
    th=pointData$POT$threshold
    d1 <- pointData$POT$peaks
    paramEstsall <- c(pointData$POT$pars[1], pointData$POT$pars[2],
                      pointData$POT$threshold, length(d1),
                      length(pointData$POT$peaks), pointData$POT$percentile)
    EVdata$GPDstat <-  list(method=methodname,
                            values=NA,
                            parameters=paramEstsall,
                            paramCIs = NA)
    message("could not estimate GPD: bounded distribution")
  }
  
  # Return outputs
  return(list(EVmeta = EVmeta, EVdata = EVdata, isValid = isValid))
}

tsEvaSampleData1 <- function(ms, meanEventsPerYear,minEventsPerYear, minPeakDistanceInDays,tail=NA) {
  
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
tsEvaSampleDataX <- function(ms, meanEventsPerYear,minEventsPerYear, minPeakDistanceInDays,tail=NA,shape_bnd) {
  
  pctsDesired = c(90, 95, 99, 99.9)
  args <- list(meanEventsPerYear = meanEventsPerYear,
               minEventsPerYear = minEventsPerYear,
               potPercentiles = c(seq(70,90,by=1), seq(91,95,by=0.5),seq(95.1,99.5,by=0.1)))
  meanEventsPerYear = args$meanEventsPerYear
  minEventsPerYear = args$minEventsPerYear
  potPercentiles = args$potPercentiles
  if(is.na(tail)) stop("tail for POT selection needs to be 'high' or 'low'")
  
  POTData <- tsGetPOTX(ms, potPercentiles, meanEventsPerYear,minEventsPerYear,minPeakDistanceInDays, tail, shape_bnd)
  
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
tsGetPOTX <- function(ms, pcts, desiredEventsPerYear,minEventsPerYear, minPeakDistanceInDays, tail, shape_bnd) {
  
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
  
  # --- Pre-allocation ---
  n_pcts        <- length(pcts)
  numperyear    <- rep(NA, n_pcts)
  minnumperyear <- rep(NA, n_pcts)
  thrsdts       <- rep(NA, n_pcts)
  gpp           <- rep(NA, n_pcts)
  devpp         <- rep(NA, n_pcts)
  fitlist       <-  vector(mode = "list", length = n_pcts)
  
  # Pre-clean data: Replace NAs once outside the loop for speed
  ms_clean <- ms[,2]
  ms_clean2=ms_clean
  if (length(which(is.na(ms_clean))>0)){
    ms_clean2<-ms_clean[-which(is.na(ms_clean))]
    ms_clean[which(is.na(ms_clean))]=-9999
  }
  if (length(which(is.infinite(ms_clean))>0)){
    print("bg")
    ms_clean2<-ms_clean[-which(is.infinite(ms_clean))]
    ms_clean[which(is.infinite(ms_clean))]=-9999
  }
  #ms_clean[is.na(ms_clean)] <- -9999
  #adf.test(ms_clean)
  skip <- 0
  #minEventsPerYear <- 1
  
  for (ipp in 1:n_pcts) {
    
    # 1. Manage the skip logic
    if (skip > 0) {
      skip <- skip - 1
      next # Move to the next iteration immediately
    }
    
    # 2. Threshold Calculation
    thrsdt <- stats::quantile(ms_clean2, pcts[ipp]/100, na.rm = TRUE)
    thrsdts[ipp] <- thrsdt
    
    # 3. Peak Finding (High vs Low Tail)
    pks <- if(tail == "high") {
      pracma::findpeaks(ms_clean, minpeakdistance = minPeakDistance, minpeakheight = thrsdt)
    } else {
      declustpeaks(data = ms_clean, minpeakdistance = minPeakDistance, 
                   minrundistance = minRunDistance, qt = thrsdt)
    }
    
    #MannKendall(pks[order(pks[,2]),1])
    # Check if any peaks were found
    if (is.null(pks) || length(pks) == 0) next
    
    # 4. Frequency Calculations
    numperyear[ipp] <- nrow(pks) / nyears
    
    # 5. Continuous Multiplicative Penalty
    # Define how 'aggressive' the penalty should be (Strength)
    # A strength of 10 means at 0 events, AIC is multiplied by 11.
    penalty_strength <- 1000
    penalty_factor <- 1
    
    if (numperyear[ipp] < minEventsPerYear) {
      # Calculate deficit (0 at the limit, 1 at zero events)
      deficit <- 1 - (numperyear[ipp] / minEventsPerYear)
      
      # Progressive multiplier: starts at 1 and grows quadratically
      penalty_factor <- 1 + (deficit^2 * penalty_strength)
    }
    
    # 6. GPD Fitting
    if (numperyear[ipp] <= (desiredEventsPerYear + 1)) {
      
      # print("hello")
      # fgpd <- suppressWarnings(try(
      #   POT::fitgpd(pks[,1], threshold = thrsdt, est = "mle",
      #               method = "BFGS", shape = gpdshape_bnd[1] ,std.err.type = "observed"),silent=TRUE))   # Upper bounds for [scale, shape]
      
      fgpd <- suppressWarnings(try(POT::fitgpd(pks[,1], threshold = thrsdt, est = "mle",
                                               method = "L-BFGS-B", lower = c(1e-6, shape_bnd[1]), # Lower bounds for [scale, shape]
                                               upper = c(Inf, shape_bnd[2]),    # Upper bounds for [scale, shape]
                                               std.err.type = "observed"),
                                   silent = TRUE))
      fitlist[[ipp]]=fgpd
      
      if (inherits(fgpd, "try-error")) {
        devpp[ipp] <- 1e9 # Keep a massive constant for hard failures
        gpp[ipp]   <- NA  
      } else {
        # Apply the Multiplicative Penalty to the AIC
        # Note: AIC can be negative. We use abs() or a shift if necessary, 
        # but usually, for GPD, we focus on the magnitude of the deviance.
        current_aic <- fgpd$deviance
        
        # Logic: If AIC is positive, multiply to make it larger (worse).
        # If AIC is negative, we must be careful: multiplying by >1 makes it 'better'.
        # We use a conditional check to ensure the penalty always hurts the score.
        if (current_aic > 0) {
          devpp[ipp] <- current_aic + penalty_factor
        } else {
          # If AIC is negative, we divide by the factor to move it closer to zero (worse)
          devpp[ipp] <- current_aic + penalty_factor
        }
        
        gpp[ipp] <- fgpd$fitted.values[2]
      }
    }
  }
  #peaks with lowest threshold (retrieving the two largest peaks)
  pkx <- declustpeaks(data = ms_clean ,minpeakdistance = minPeakDistance ,minrundistance = minRunDistance, qt=stats::quantile(ms_clean,pcts[1]/100,na.rm=T))
  #md= abs(pkx[1,1]-pkx[2,1])
  devpp[1]=NA
  
  # --- 1. Pre-calculate Distances and Constraints ---
  # Distance to nearest boundary
  dist_to_bnd <- pmin(abs(gpp - shape_bnd[1]), abs(gpp - shape_bnd[2]))
  # Calculate step changes in shape (stability check)
  dshap <- c(0, diff(gpp))
  
  # Determine the threshold for "low tail" logic based on peak distance
  low_tail_threshold <- shape_bnd[2]
  
  # --- 2. Define the "Search Space" using a Mask ---
  # Start with your original deviance (AIC)
  devpx <- devpp
  
  # Apply Hard Constraints (Eliminate impossible or unstable fits)
  # Instead of +9999, we use NA to strictly remove them from consideration
  if (tail == "low") {
    devpx[gpp >= low_tail_threshold] <- NA
  }
  
  # Stability check: Remove shapes with massive jumps
  devpx[abs(dshap) > 0.5] <- NA
  print(devpx)
  # --- 3. Selection Logic ---
  
  # Check if we have any valid candidates left within the boundaries
  trip_in_bounds <- which(!is.na(devpx) & gpp >= shape_bnd[1] & gpp <= shape_bnd[2])
  
  if (length(trip_in_bounds) > 0) {
    # BEST CASE: Best AIC within your safe boundaries
    trip <- trip_in_bounds[which.min(devpx[trip_in_bounds])]
    message(paste0("Optimal shape found: ", round(gpp[trip], 3)))
    
  } else {
    # FALLBACK: No valid candidates inside boundaries. 
    # We balance AIC performance with proximity to the boundary.
    
    # Remove NA values for the calculation
    
    valid_indices <- which(!is.na(devpx))
    print(valid_indices)
    if (length(valid_indices) == 0) {
      # Extreme Safety: if EVERY fit was rejected by constraints, revert to raw devpp
      valid_indices <- 1:length(devpp)
      devpx <- devpp 
    }
    
    # Normalize remaining AIC and Distances to 0-1 scale for fair comparison
    v_dev   <- devpx[valid_indices]
    v_dist  <- dist_to_bnd[valid_indices]
    
    norm_dev  <- (v_dev - min(v_dev, na.rm=T)) / (diff(range(v_dev, na.rm=T)) + 1e-10)
    norm_dist <- (v_dist - min(v_dist, na.rm=T)) / (diff(range(v_dist, na.rm=T)) + 1e-10)
    
    # Calculate combined score (lower is better)
    combined_score <- (0.5 * norm_dev) + (0.5 * norm_dist)
    
    # Map back to the original index
    trip <- valid_indices[which.min(combined_score)]
    
    message(paste0("Shape outside boundaries. Fallback selected: ", round(gpp[trip], 3)))
  }
  
  isok <- TRUE
  print(numperyear)
  message(paste0("\nmax threshold is: ", pcts[trip],"%"))
  message(paste0("\nshape parameter is: ", round(gpp[trip],2)))
  message(paste0("\naverage number of events per year = ",round(numperyear[trip],1) ))
  
  diffNPerYear <- mean(diff(stats::na.omit(rev(numperyear)), na.rm = TRUE))
  if (diffNPerYear == 0) diffNPerYear <- 1
  diffNPerYear <- 1
  thresholdError <- -mean(diff(stats::na.omit(thrsdts))/diffNPerYear)/2
  indexp <- trip
  print(indexp)
  if (length(indexp)>0) {
    thrsd <- stats::quantile(ms_clean2,pcts[indexp]/100)
    pct <- pcts[indexp]
    fit=fitlist[[indexp]]
  } else {
    thrsd <- 0
    pct=0
  }
  # Find peaks in the second column of the matrix 'ms'
  if(tail=="high") pks_and_locs <- pracma::findpeaks(ms_clean,minpeakdistance = minPeakDistance, minpeakheight = thrsd)
  if(tail=="low") pks_and_locs <- declustpeaks(data = ms_clean ,minpeakdistance = minPeakDistance ,minrundistance = minRunDistance, qt=thrsd)
  
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
  POTdata[['pars']] <- fit
  
  
  return(POTdata)
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
        #print(shape_bnd)
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
  
  #peaks with lowest threshold (retrieving the two largest peaks)
  pkx <- declustpeaks(data = ms[,2] ,minpeakdistance = minPeakDistance ,minrundistance = minRunDistance, qt=quantile(ms[,2],pcts[1]/100,na.rm=T))
  md= abs(pkx[1,1]-pkx[2,1])
  devpp[1]=NA
  if(is.na(trip)){
    isok=F
    devpx=devpp
    count=1
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
      if(count>(length(devpx)-1)){
        #safety measure for stability of parameter
        trip=which.min(devpp)
        message(paste0("shape outside boudaries: ",round(gpp[trip],2)))
        isok=T
      }
    }
  }
  # plot(pcts,devpp,ylim=c(0,1e5))
  # plot(pcts,gpp)
  # print(devpp)
  # print(pcts)
  message(paste0("\nmax threshold is: ", pcts[trip],"%"))
  message(paste0("\nshape parameter is: ", round(gpp[trip],2)))
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


outletopen=function(dir,outletname,nrspace=rep(NA,5)){
  ncbassin=paste0(dir,"/",outletname,".nc")
  ncb=nc_open(ncbassin)
  name.vb=names(ncb[['var']])
  namev=name.vb[1]
  #time <- ncvar_get(ncb,"time")
  
  #timestamp corretion
  name.lon="lon"
  name.lat="lat"
  if (!is.na(nrspace[1])){
    start=as.numeric(nrspace[c(2,4)])
    count=as.numeric(nrspace[c(3,5)])-start+1
  }else{
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
  lon=rep(londat[idc],length(time))
  lat=rep(latdat[idc],length(time))
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
plotchangemaps=function(basemap,catmap,datar,law="GPD",type,period=c(1950,1990),parlist,valuenames){
  
  
  data=full_join(catmap,datar,by = c("outlets"="unikout"))
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
      tmpval[,it][which(tmpval[,it]>80)]=NA
    }
    data[,valcol]=tmpval
    br=c(-100,-80,-60,-40,-20,0,20,40,60,80,100)
    limi=c(-100,100)
    trans=scales::modulus_trans(.8)
    colNA="darkgrey"
  }
  if (type=="RPchange"){
    palet=c(hcl.colors(9, palette = "BrBG", alpha = NULL, rev = T, fixup = TRUE))
    title=paste0(period[2]," return period of a 100-years Return Level flood in ", period[1])
    legend="Return Period (Years)"
    br=c(1,10,100,1000,10000)
    limi=c(1,10000)
    trans="log"
    newRP=RPchangeCal(parlist, yi=period[1], yf=period[2], RetLev=datar,law, valuenames)
    ow=match(data$outlets,newRP$catchment)
    data[,colsel]=newRP$newRP[ow]
    colNA="black"
    
  }
  
  
  names(data)[colsel]="fillplot"
  
  ggplot(basemap) +
    geom_sf(fill="white")+
    geom_sf(data=data,aes(fill=fillplot,geometry=geometry),color="transparent")+ 
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
plotRegime <- function(data, catch, seuil=NULL){
  #~ data <- as.simu(data)
  names(data)=c("date","Q")
  if( is.null(seuil) ) seuil <- -Inf
  ## Suppression des lignes sans couples Qobs/Qsim
  data <- subset( data, !is.na(Q) & Q > seuil)
  ## Preparation de la sequence de mois
  deb=data$date[1]
  mois.deb <-  seq(as.Date(deb), by="month", length=12)
  
  period=paste0(format(range(data$date)[1],"%b %Y"),"-",format(range(data$date)[2],"%b %Y"))
  main = bquote(.(catch)~ " Regime | Period: "~ .(period))
  jours <- as.numeric(format(data$date,"%j"))
  ## Iddinces des lignes de meme categorie
  ind.j <- tapply(seq(length(jours)), jours, c)
  ind.j <- ind.j[-366]
  ## Calcul
  Qc <- data.frame(date=as.numeric(names(ind.j)),
                   mean=sapply(ind.j, function(x) mean(data$Q[x], na.rm=TRUE)),
                   q10=sapply(ind.j, function(x) quantile(data$Q[x],0.1, na.rm=TRUE)),
                   q90=sapply(ind.j, function(x) quantile(data$Q[x],.9, na.rm=TRUE)))
  
  qlim=c(0,1.2*max(Qc[,4]))
  plot(Qc$date, Qc$mean, type="n", axes=FALSE, ylim=qlim,
       xlab = NA, ylab = expression(paste("Debit (",m^3/s,")")))
  mtext(main,3,font = 2,line = 0.5,cex = .75)
  lines(Qc$date, Qc$mean, col="darkblue",lwd=2)
  lines(Qc$date, Qc$q10, col="blue",lwd=1,lty=2)
  lines(Qc$date, Qc$q90, col="blue",lwd=1,lty=2)
  abline(v = format(mois.deb,"%j"), col="lightgrey", lty=3)
  abline(h = quantile(data$Q,.95),col="red",lwd=2)
  polygon(c(Qc$date,rev(Qc$date)),c(Qc$q10,rev(Qc$q90)),
          col=alpha("lightblue",.4),border="transparent")
  axis(2)
  axis(1, format(mois.deb,"%j"), label=format(mois.deb,"%b"),cex.axis=.8)
  box()
  legend("topleft", leg=c("mean","Q90","Q10"),
         lwd=c(2,2,1,2), col=c("darkblue","blue","blue"),
         cex=.8, lty=c(1,2,2),bty="n")
  return(Qc)
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
    start_index=1
    indices_to_extract <- seq(from = start_index, to = length(data$Q7), by = WindowSize/dt)
    datat=data$Q7[indices_to_extract]
    l0=length(which(datat<=1e-4))
    f0=1/(l0/(length(datat)/52))
    dayysbelow=tsEvaNanRunnigBlowTh(series=data$Q7,threshold=1e-4,windowSize=4*365*30)
    
    dayysbelow$time=as.Date(data$date[dayysbelow$time])
    yrtot=length(unique(year(data$date)))
    list0=NA
    fl=0 #flag for intermitent river
    if (l0>=7){
      print("intermittent river class2")
      fl=3
      data$Qinv=NA
      data$Qlninv
      #keep days with 0 discharge
      list0=data$date[which(data$Qs==0)]
    }else if(l0>1){ 
      print("intermittent river class1 ")
      fl=1
      list0=data$date[which(data$Q7==mindis)]
    }else if(mindis>0 & m0>=yrtot){ 
      print(paste0("river with floor low flow ",mindis))
      fl=2
      list0=data$date[which(data$Q7==mindis)]
    }
    dis07=data[,c(2,3,4)]
  }
  #The objective here is to return also the dicharge in a better format for next step of the analysis
  return(list(zerodate=list0,trdis=dis07,DaysBlow=dayysbelow,flags=c(n0d=l0,intertype=fl,mindischarge=mindis)))
}

#Open upstream area file for selected pixels
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
  
  #outll$idlalo=paste(outll$idlo,outll$idla,sep=" ")
  outll$latlong=paste(round(outll$Var1,4),round(outll$Var2,4),sep=" ")
  outfinal=inner_join(outll, Sloc_final, by="latlong")
  return (outfinal)
}




#Set data directory
hydroDir<-("D:/tilloal/Documents/LFRuns_utils/ChangingHydroExtremes/data")

# Arguments importation #############################

#default: tail is high
tail="high"
haz = "drought"
var = "dis"
outlets="RNetwork"
outletname <- "/GeoData/efas_rnet_100km_01min"
season="year"
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


#LOAD THE RESULTS FROM THE SOCIOCF RUN
haz="drought"
if (haz == "drought") namefile="Drought.nonfrost.SocCF"
if (haz == "flood") namefile="flood.year.SocCF"
load(file=paste0(hydroDir,"/",haz,"/peaks.",namefile,".Rdata"))
load(file=paste0(hydroDir,"/",haz,"/params.",namefile,".Rdata"))
#keep only points in the desired square
Paramsfl$square=round(Paramsfl$catchment/100000)
Paramsfl=Paramsfl[which(Paramsfl$square==Nsq),]

Peaksave$square=round(Peaksave$catch/100000)
Peaksave=Peaksave[which(Peaksave$square==Nsq),]
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

filenameSCF=paste0("dis_",Nsq,"_1951_2020_scf_RNetwork")

dists=disNcopenloc(filename,hydroDir,outhybas,1)
df.dis=dists 
print(paste0("hazard: ",haz," opening square ", Nsq, " /88"))
timeStamps=(as.Date(df.dis$time,origin="1979-01-01"))
timeStamps=as.POSIXct(timeStamps-1/24)
txx=timeStamps
df.dis$timeStamps=txx

names(df.dis)[c(1,2)]=c("dis","outlets")
haz="drought"
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
TH1=read.csv(paste0(ThDir,"/trenTH_Histo_",tail,"_",Nsq,".csv"))
TH2=read.csv(paste0(ThDir,"/trenTH_SCF_",tail,"_",Nsq,".csv"))
TH3=inner_join(TH2,TH1,by="cid")

#retain thresholds fro historical run unless it is NA
thresh_vec=data.frame(TH3$cid, TH3$Th_new.x)
if(length(which(is.na(thresh_vec$TH3.Th_new.x)))>0){
  print("corr")
  thresh_vec$TH3.Th_new.x[which(is.na(thresh_vec$TH3.Th_new.x))]=TH3$Th_new.y[which(is.na(thresh_vec$TH3.Th_new.x))]
}
names(thresh_vec)=c("cid","th")
thresh_vec$cid=as.numeric(thresh_vec$cid)
Nsq=as.numeric(Nsq)
thresh_vec$cid=thresh_vec$cid-Nsq*10000
thresh_vec$cid=thresh_vec$cid+Nsq*100000

#Bigloop ----------
startid=1
endid=length(unikout)
endid=7
haz="drought"
RetPerGPD=c()
RetPerGEV=c()
RetLevGEV=c()
RetLevGPD=c()
parlist=c()
peaklist=c()
catlist=c()
IRES=c()
for (idfix in startid:endid){
  start_time <- Sys.time()
  print(paste0("hazard:",haz," square: ", Nsq, " pixel: ",idfix,"/",endid))
  catch=as.numeric(unikout[idfix])
  
  
  timeStamps=txx
  thresh=thresh_vec[which(thresh_vec$cid==catch),]
  thresh=thresh$th
  #  shape_bnd=rbind(shape_bnd_gpd[idfix,],shape_bnd_gev[idfix,])
  frosttime=NA
  df.disX=disNcopenloc(filename,hydroDir,outhybas,idfix)
  series=data.frame(txx,df.disX$outlets)
  names(series)=c("date","Qs")
  rmv=which(as.integer(format(series$date, "%Y"))==1950)
  if (length(rmv)>0){
    series=series[-rmv,]
  }
  #remove first day to avoid errors
  series=series[-1,]
  if (haz=="drought"){
    trans="rev"
    #seasonal split
    catmat=Catf7[which(Catf7$outlets==catch),]
    Tcatmat=mycat[which(mycat$HYBAS_ID==catmat$HYBAS_ID),]
    Tcatchment=which(colnames(frostcat)==Tcatmat$pointid)
    
    intermit=interid(series,trans,WindowSize=7)
    interflag=intermit$flags[2]
    series=data.frame(series$date,intermit$trdis$Q7)
    #remove frost timesteps, this can be modified to do the anlysis only on frost moments
    if (length(Tcatchment)>0){
      frostserie=data.frame(frostcat[,1],frostcat[,Tcatchment])
      frosttime=which(frostserie[,2]<0)
    }else{
      frosttime=NA
    }
    ciPercentile=80
    minPeakDistanceInDays=30
    tail="low"
  }else if (haz=="flood"){
    
    ciPercentile=95
    minPeakDistanceInDays=7
    interflag=0
    series <- max_daily_value(series)
    tail="high"
    trans="ori"
  }
  
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
  
  nv=length(unique(series$dis))
  if(length(na.omit(series$dis))>1 & interflag<3 & nv>15){
    
    if (length(which(is.na(series$dis)))>0){
      print("Na alert")
      seriefill=tsEvaFillSeries(series$timestamp,series$dis)
      series$dis=seriefill
    }
    timeAndSeries=series
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
    series=timeAndSeries[,2]
    
    timeWindow = 365.25*30; #time windows in days, the correction is done within the functions
    windowSize=366
    
    timeStamps=timeAndSeries$timestamp
    cat(paste0("\nsquare: ", Nsq, " pixel: ",idfix,"/",endid))
    #TSEVA-----
    Nonstat<-TsEvaNs(timeAndSeries, timeWindow, transfType='trendPeaks',ciPercentile = ciPercentile, 
                      minPeakDistanceInDays = minPeakDistanceInDays,lowdt=7,trans=trans,tail = tail,TrendTh = thresh,shape_bnd=NA)
    nonStationaryEvaParams=Nonstat[[1]]
    stationaryTransformData=Nonstat[[2]]
    oo=nonStationaryEvaParams$gevObj$parameters$epsilon[1]
    print(oo)
    
    stationaryTransformData$timeStampsDay=unique(as.Date(stationaryTransformData$timeStamps))
    pikos=data.frame(nonStationaryEvaParams$potObj$parameters$peaks,nonStationaryEvaParams$potObj$parameters$peakID,nonStationaryEvaParams$potObj$parameters$peakST, nonStationaryEvaParams$potObj$parameters$peakEN)
    names(pikos)=c("value","timeID","tIDstart","tIDend")
    pikos$time=timeStamps[pikos$timeID]
    pikos$catch=rep(catch,length(pikos[,1]))
    
    dt1=min(diff(timeStamps),na.rm=T)
    dt=as.numeric(dt1)
    tdim=attributes(dt1)$units
    if (tdim=="hours") dt=dt/24
    
    #change this part
    if (dt==1){
      timeDays=stationaryTransformData$timeStamps
    }else{
      timeDays=stationaryTransformData$timeStampsDay
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
    RPgoal=10
    timeIndex=tindexes[1]
    RLevs100=ComputeReturnLevels(nonStationaryEvaParams, RPgoal, timeIndex)
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
      RLgev=RLevs100$ReturnLevels[2]
      RLgpd=RLevs100$ReturnLevels[3]
      ERgev=RLevs100$ReturnLevels[4]
      ERgpd=RLevs100$ReturnLevels[5]
      
      nRPgev=nRPgpd=10
      params=c()
      for (t in 1:length(Impdates)){
        timeIndex=tindexes[t]
        RLevs100i=ComputeReturnLevels(nonStationaryEvaParams, RPgoal, timeIndex)
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
    cat(paste0("\n No values in this pixel ",idfix," \n or intermittent river (flag = ",interflag,")"))
    if (is.na(interflag)) interflag=-9999
    if (interflag>0){
      filling=intermit$flags[3]
      datex=yday(intermit$DaysBlow$time)
      dtect=c(diff(datex),-1)
      last_days <- intermit$DaysBlow$time[which(dtect<0)]
      tindexes=match(last_days,intermit$DaysBlow$time)
      oops=intermit$DaysBlow[tindexes,]
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
    peaklist=as.data.frame((rbind(peaklist,pikos)))
    
  }
  #Saving main outputs
  catlist=c(catlist,catch)
  IRES=c(IRES,interflag)
  
  RetLevGEV=rbind(RetLevGEV,RLgev)
  
  RetLevGPD=rbind(RetLevGPD,RLgpd)
  
  RetPerGEV=rbind(RetPerGEV,nRPgev)
  
  RetPerGPD=rbind(RetPerGPD,nRPgpd)
  end_time <- Sys.time()
  cat(paste0("\nloop duration: ",round(end_time-start_time,2)," seconds\n"))
  
}
parlist=as.data.frame(parlist)
epsilons=parlist[which(parlist$Year==1952),]

shape_bnd
mate=match(epsilons$catchment,EpsilonTab$catchment)
epsilonI=EpsilonTab[mate,]

matp=match(epsilons$catchment,ParamsflH$catchment)
ParamHI=ParamsflH[matp,]
ParamSI=ParamsflSCF[matp,]

plot(ParamHI$sigmaGPD,epsilons$sigmaGPD)
plot(ParamHI$thresholdGPD,ParamSI$thresholdGPD)
hist((epsilonI$epsilonH-epsilonI$epsilonS))
hist((epsilonI$epsilonH-epsilons$epsilonGPD))

Results=list(parameters=parlist,RetLevGEV=RetLevGEV,RetLevGPD=RetLevGPD,RetPerGEV=RetPerGEV,RetPerGPD=RetPerGPD,Peaks=peaklist,catrest=data.frame(catlist,IRES))

plot(peaklist$value[which(peaklist$catch==catlist[1])][order(peaklist$timeID)],ylim=c(0,1.4*max(t(RetLevGPD[1,]))))
plot(t(RetLevGPD[1,]),type="l")








#load results from previous method and check
id=4
haz="flood"
if (haz == "drought") namefile="Drought.nonfrost.Histo"
if (haz == "flood") namefile="flood.year.Histo"
load(file=paste0(hydroDir,"/",haz,"/RL100.",namefile,".Rdata"))
#keep only points in the desired square
RLGPDfl_1=RLGPDfl[which(RLGPDfl$unikout==catlist[id])]

if (haz == "flood") namefile="flood.year.SocCF"
if (haz == "drought") namefile="Drought.nonfrost.SocCF"
load(file=paste0(hydroDir,"/",haz,"/RL100.",namefile,".Rdata"))
load(file=paste0(hydroDir,"/",haz,"/RL100.",namefile,".Rdata"))
RLGPDfl_2=RLGPDfl[which(RLGPDfl$unikout==catlist[id])]

#plot(peaklist$value[which(peaklist$catch==catlist[1])][order(peaklist$timeID)],ylim=c(0,1.4*max(t(RetLevGPD[2,]))))
corRL=RetLevGPD[id,1]-t(RLGPDfl_2)[1]
plot(t(RetLevGPD[id,]),type="l",ylim=c(0.6*min(RetLevGPD[id,]),1.4*max(RetLevGPD[id,])))
lines(t(RLGPDfl_1)[-71],col=2)
lines(t(RLGPDfl_2)[-71],col=3)

plot(t(RLGPDfl_1)[-71]-t(RLGPDfl_2)[-71])
plot(t(RetLevGPD[id,])-t(RLGPDfl_2)[-71])






