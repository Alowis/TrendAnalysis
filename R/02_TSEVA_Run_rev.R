##########################################################################################
############   THIS SCRPT IS FOR RUNNING 1 SQUARE OF THE HERA DOMAIN ON A HPC ############
##########################################################################################

setwd(dirname(rstudioapi::getSourceEditorContext()$path))
getwd()
source("functions_trends.R")

tsGetPOT_2 <- function(ms, pcts, desiredEventsPerYear,minEventsPerYear, minPeakDistanceInDays, tail) {
  
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
        thrsdt <- stats::quantile(ms[,2],pcts[ipp]/100,na.rm=T)
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
            devpp[ipp]=1e9
            gpp[ipp]=9999
          } else {
            gpdpar = fgpd$fitted.values # [scale, shape]
            
            # 1. Calculate the Cumulative Distribution Function (CDF) values for the peaks
            # Using the GPD formula: F(x) = 1 - (1 + shape * (x-thresh)/scale)^(-1/shape)
            scaled_peaks = (pks[,1] - thrsdt) / gpdpar[1]
            if(abs(gpdpar[2]) < 1e-10) { # Handle case where shape is nearly 0 (Exponential)
              z = 1 - exp(-scaled_peaks)
            } else {
              z = 1 - (1 + gpdpar[2] * scaled_peaks)^(-1/gpdpar[2])
            }
            
            # 2. Calculate Right-Tail Weighted Anderson-Darling (ADR)
            # This version emphasizes the upper tail of the distribution
            z = sort(z)
            n = length(z)
            i = 1:n
            # ADR formula: -n/2 - 2*sum(z) - sum((2*i-1)*log(1-z))/n (simplified version)
            # Note: We use a version that targets the right tail specifically
            adr_stat = -n - (1/n) * sum((2*i - 1) * log(z) + (2*n + 1 - 2*i) * log(1 - z))
            
            # 3. Store the statistic plus your existing performance penalties
            devpp[ipp] = adr_stat + perfpen 
            gpp[ipp] = gpdpar[2]
          }
          
          nperYear <- tsGetNumberPerYear(ms, pks[,2])
          minnumperyear[ipp] <- min(nperYear$Freq, na.rm = TRUE)
          
          # if(numperyear[ipp]<=desiredEventsPerYear+1 & dej==0){
          # fgpd=suppressWarnings(try(POT::fitgpd(pks[,1], threshold = thrsdt, est = "mle",method="BFGS",std.err.type = "expected")))
          # if(inherits(fgpd, "try-error")){
          # gpdpar=9999
          # deviance=9999
          # devpp[ipp]=1e9
          # gpp[ipp]=9999
          # }else {
          # gpdpar=fgpd$fitted.values
          # deviance=fgpd$deviance
          # devpp[ipp]=stats::AIC(fgpd)+perfpen
          # gpp[ipp]=gpdpar[2]
          # }
          # nperYear <- tsGetNumberPerYear(ms, pks[,2])
          # minnumperyear[ipp] <- min(nperYear$Freq, na.rm = TRUE)
        }
      }
    }
  }
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
        thrsdt <- stats::quantile(ms[,2],pcts[ipp]/100,na.rm=T)
        thrsdts[ipp] <- thrsdt
        ms[,2][which(is.na(ms[,2]))]=-9999
        minEventsPerYear=1
        
        if(tail=="high") {
          # #boundaries of shape parameter
          # shape_bnd=c(-0.5,1)
          pks <- pracma::findpeaks(ms[,2],minpeakdistance = minPeakDistance, minpeakheight = thrsdt)
        }
        if(tail=="low") {
          pks <- declustpeaks(data = ms[,2] ,minpeakdistance = minPeakDistance ,minrundistance = minRunDistance, qt=thrsdt)
          # shape_bnd=c(-2,0)
        }
        numperyear[ipp] <- length(pks[,1])/nyears
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
            devpp[ipp]=stats::AIC(fgpd)+perfpen
            gpp[ipp]=gpdpar[2]
          }
          nperYear <- tsGetNumberPerYear(ms, pks[,2])
          minnumperyear[ipp] <- min(nperYear$Freq, na.rm = TRUE)
        }
      }
    }
  }
  
  #peaks with lowest threshold (retrieving the two largest peaks)
  pkx <- declustpeaks(data = ms[,2] ,minpeakdistance = minPeakDistance ,minrundistance = minRunDistance, qt=stats::quantile(ms[,2],pcts[1]/100,na.rm=T))
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
      devpp[which(abs(dshap)>0.5)]=devpp[which(abs(dshap)>0.5)]+9999
      trip=which.min(devpp)
      #message(paste0("shape outside boudaries: ",round(gpp[trip],2)))
      #isok=T
      #trip=which.min(devpx)
      isok=dplyr::between(round(gpp[trip],1), shape_bnd[1], shape_bnd[2])
      count=count+1
      if(isok==F)devpx[trip]=NA
      if(count>(length(devpx)-1)){
        #safety measure for stability of parameter
        trip=which.min(devpp)
        message(paste0("shape outside boudaries: ",round(gpp[trip],2)))
        isok=T
      }
    }
  }
  message(paste0("\nmax threshold is: ", pcts[trip],"%"))
  message(paste0("\nshape parameter is: ", round(gpp[trip],2)))
  message(paste0("\naverage number of events per year = ",round(numperyear[trip],1) ))
  
  diffNPerYear <- mean(diff(stats::na.omit(rev(numperyear)), na.rm = TRUE))
  if (diffNPerYear == 0) diffNPerYear <- 1
  diffNPerYear <- 1
  thresholdError <- -mean(diff(stats::na.omit(thrsdts))/diffNPerYear)/2
  indexp <- trip
  if (!is.na(indexp)) {
    thrsd <- stats::quantile(ms[,2],pcts[indexp]/100)
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
#New function test
RPcalcGPD<-function(params,RPiGPD){
  paramx=data.frame(t(params))
  X0 <- paramx$nPeaks/paramx$SampleTimeHorizon
  qxD=(((1+paramx$epsilonGPD*(RPiGPD-paramx$thresholdGPD)/paramx$sigmaGPD)^(-1/paramx$epsilonGPD)))
  if (is.na(qxD)){
    returnPeriodGPD=9999
  }else{
    returnPeriodGPD=1/(X0*qxD)
  }
  return(c(GPD=returnPeriodGPD))
}
ComputeReturnLevelsGPD<-function(nonStationaryEvaParams, RPgoal, timeIndex){
  
  
  #GPD
  epsilonGPD <- nonStationaryEvaParams[[1]]$parameters$epsilon
  sigmaGPD <- mean(nonStationaryEvaParams[[1]]$parameters$sigma[timeIndex])
  thresholdGPD <- mean(nonStationaryEvaParams[[1]]$parameters$threshold[timeIndex])
  nPeaks <- nonStationaryEvaParams[[1]]$parameters$nPeaks
  thStart <- nonStationaryEvaParams[[1]]$parameters$timeHorizonStart
  thEnd <- nonStationaryEvaParams[[1]]$parameters$timeHorizonEnd
  sampleTimeHorizon <- as.numeric((thEnd - thStart)/365.2425)
  
  if (nonStationaryEvaParams[[1]]$method=="No fit"){
    
    ParamGPD=c(epsilonGPD,sigmaGPD,thresholdGPD,NA,NA, NA,nPeaks,sampleTimeHorizon)
    names(ParamGPD)=c("epsilonGPD","sigmaGPD","thresholdGPD","epsilonStdErrGPD","sigmaStdErrGPD","thresholdStdErrGPD","nPeaks","SampleTimeHorizon")
    return(list(Fit="No fit",Params=c(ParamGEV,ParamGPD)))
  }else{
    
    
    #GPD
    # epsilonGPD <- nonStationaryEvaParams[[2]]$parameters$epsilon
    # sigmaGPD <- mean(nonStationaryEvaParams[[2]]$parameters$sigma[timeIndex])
    # thresholdGPD <- mean(nonStationaryEvaParams[[2]]$parameters$threshold[timeIndex])
    # nPeaks <- nonStationaryEvaParams[[2]]$parameters$nPeaks
    epsilonStdErrGPD <- nonStationaryEvaParams[[1]]$paramErr$epsilonErr
    sigmaStdErrGPD <- mean(nonStationaryEvaParams[[1]]$paramErr$sigmaErr[timeIndex])
    thresholdStdErrGPD <- mean(nonStationaryEvaParams[[1]]$paramErr$thresholdErr[timeIndex])
    # thStart <- nonStationaryEvaParams[[2]]$parameters$timeHorizonStart
    # thEnd <- nonStationaryEvaParams[[2]]$parameters$timeHorizonEnd
    # sampleTimeHorizon <- as.numeric((thEnd - thStart)/365.2425)
    
    #returnLevelsGEV <- tsEvaComputeReturnLevelsGEV(epsilonGEV, sigmaGEV, muGEV, epsilonStdErrGEV, sigmaStdErrGEV, muStdErrGEV, RPgoal)
    
    returnLevelsGPD <- tsEvaComputeReturnLevelsGPD(epsilonGPD, sigmaGPD, thresholdGPD, epsilonStdErrGPD, sigmaStdErrGPD, thresholdStdErrGPD,
                                                   nPeaks, sampleTimeHorizon, RPgoal)
    
    rlevGPD=returnLevelsGPD$returnLevels
    
    
    errGPD=returnLevelsGPD$returnLevelsErr
    
    ParamGPD=c(epsilonGPD,sigmaGPD,thresholdGPD,epsilonStdErrGPD,sigmaStdErrGPD, thresholdStdErrGPD,nPeaks,sampleTimeHorizon)
    names(ParamGPD)=c("epsilonGPD","sigmaGPD","thresholdGPD","epsilonStdErrGPD","sigmaStdErrGPD","thresholdStdErrGPD","nPeaks","SampleTimeHorizon")
    return(list(Fit="Fitted",ReturnLevels=c(ReturnPeriod=RPgoal,GPD=as.numeric(rlevGPD),errGPD=as.numeric(errGPD)),Params=c(ParamGPD)))
  }  
  
  
}


TsEvaNsX<- function(timeAndSeries, timeWindow, transfType='trendPeaks',minPeakDistanceInDays=10,
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
  
  
  minEventsPerYear = 1.6
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
  pointData = tsEvaSampleData1(ms, potEventsPerYear, minEventsPerYear, minPeakDistanceInDays,tail);
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




TsEvaNs_fromExtremes <- function(
    timeStamps,
    extremePoints,        # data.frame: col1=timestamps, col2=extreme values (pre-filtered)
    originalPoints,       # data.frame: col1=timestamps, col2=original series values at same indices
    series,
    Peaks_0,
    Peaks_1,
    Peaks_tid,
    Peaks_loc,
    trasfData,            # output$stationaryTransformData from a previous TsEvaNs call
    timeWindow,
    minPeakDistanceInDays = 10,
    gevMaxima = 'annual',
    gevType = 'GEV',
    evdType = c('GEV', 'GPD'),
    tail = "high",
    epy = 3,
    shape_bnd = c(-0.5, 1)
) {
  
  # --- 1. Parse inputs -------------------------------------------------------
  
  newSeries   <- extremePoints[, 2]
  origSeries  <- originalPoints[, 2]   # original values at the same time indices
  
  dt1  <- min(diff(timeStamps), na.rm = TRUE)
  dtn  <- as.numeric(dt1)
  tdim <- attributes(dt1)$units
  if (!is.null(tdim)) {
    if (tdim == "hours")   dtn <- dtn / 24
    if (tdim == "seconds") dtn <- dtn / 3600
  }

  nRunMn <- ceiling(timeWindow / dtn)
  
  # --- 2. Adjust trend -------------------------------------------------------
  # The caller supplies trasfData (from a parent TsEvaNs run).
  # We extend its trendSeries by the point-wise shift: orig - new.
  # This absorbs the difference between the full-series trend and the
  # subset-specific level into the non-stationary location parameter.
  
  
  pointShift <-  newSeries - origSeries
  N <- timeWindow * 4
  trend_error <- sd(pointShift) / sqrt(N)

  lsm=length(pointShift)/10
  plot(pointShift)
  smoothShift=tsEvaNanRunningMean(pointShift,lsm)
  plot(smoothShift,col=2)
  
  # trasfData$trendSeries is defined over the full original time axis;
  # we subset it to the extreme-point positions and add the shift.
  trendSerieX=trasfData$trendSeries
  trendSerieX[-ids_above_threshold]=NA
  adjustedTrend <- trasfData$trendSeries[ids_above_threshold] + smoothShift
  trendSerieX[ids_above_threshold]=adjustedTrend
  # plot(trendSerieX)
  # points(trasfData$trendSeries)
  trendSerieX=approx(trendSerieX,n=length(trendSerieX))$y

  series_above_threshold<-series
  series_above_threshold[-ids_above_threshold]=NA

  
  detrend_result <- tsEvaDetrendTimeSeries(timeStamps, series_above_threshold, 
                                           timeWindow)
  
  n_run_mean <- detrend_result@nRunMn
  trendSerieX=tsEvaNanRunningMean(trendSerieX,n_run_mean/60)

  lli=c(min(min(trasfData$trendSeries,min(trendSerieX))),max(max(trasfData$trendSeries,max(trendSerieX))))
  #plot(detrend_result@trendSeries,type="l",col=3)
  plot(trendSerieX,type="l")
  lines(trasfData$trendSeries,col=2)
  lines(detrend_result@trendSeries,type="l",col=3)
  trendSerieX=detrend_result@trendSeries

  detrended_series_above_threshold <- series_above_threshold - 
    trendSerieX

  variance_series <- tsEvaNanRunningVariance(detrended_series_above_threshold, 
                                             n_run_mean)
  variance_series <- tsEvaNanRunningMean(variance_series, ceiling(n_run_mean/2))
  std_dev_series <- sqrt(variance_series)

  plot(std_dev_series)
  lines(trasfData$stdDevSeries)
  avg_std_dev <- mean(std_dev_series, na.rm = TRUE)
  adjustedStdDev<-std_dev_series[ids_above_threshold]
  S <- 2
  N <- timeWindow * 4
  adjustedStdDevError <- avg_std_dev * (2 * S^2/N^3)^(1/4)* rep(1, length(std_dev_series))
  adjustedTrendError  <- trend_error
  
  # --- 3. Build a minimal stationarised series for EVA ----------------------
  # Stationarise: remove trend and rescale by stdDev, exactly as TsEvaNs does
  stationarySeries <- (series - trendSerieX) / std_dev_series
  # mean(abs(stationarySeries-trasfData$stationarySeries),na.rm=T)
  # stationarySeries <- (series - trasfData$trendSeries) / trasfData$stdDevSeries
  # mean(abs(stationarySeries-trasfData$stationarySeries),na.rm=T)
  # plot(stationarySeries-trasfData$stationarySeries)
  # points(trasfData$stationarySeries, col=2)
  adjusted_tpeaks=trendSerieX[Peaks_tid][order(Peaks_tid)]
  adjusted_stdpeaks=std_dev_series[Peaks_tid][order(Peaks_tid)]
  stationary_peaks= (Peaks_1-adjusted_tpeaks)/adjusted_stdpeaks
  # adjusted_tpeaks=trasfData$trendSeries[Peaks_tid][order(Peaks_tid)]
  # adjusted_stdpeaks=trasfData$stdDevSeries[Peaks_tid][order(Peaks_tid)]
  # stationary_peaks_0= (Peaks_0-adjusted_tpeaks)/adjusted_stdpeaks
  # plot(stationary_peaks)
  # points(stationary_peaks_0,col=2)
  # sd(stationary_peaks)
  # adjusted_tpeaks=trasfData$trendSeries[Peaks_loc$timeID][order(Peaks_loc$timeID)]
  # adjusted_stdpeaks=trasfData$stdDevSeries[Peaks_loc$timeID][order(Peaks_loc$timeID)]
  # stationary_peaks= (Peaks_loc$value[order(Peaks_loc$timeID)]-adjusted_tpeaks)/adjusted_stdpeaks
  # plot(stationary_peaks)
  # sd(stationary_peaks)
  
  #plot(Peaks_loc$timeID,Peaks_loc$value)

  #ms <- data.frame(timeStamps, stationarySeries)
  
  # --- 4. Set POT/GEV sampling parameters -----------------------------------
  potEventsPerYear <- epy
  
  
  # --- 5. GPD threshold = min of the (stationarised) subset -----------------
  # In the stationary space the threshold corresponds to the minimum peak,
  # i.e. every point in the subset is "above threshold" by construction.
  gpdThreshold <- min(stationary_peaks, na.rm = TRUE)
  
  percentile <- ecdf(stationarySeries)(gpdThreshold)
  # --- 6. Run stationary EVA ------------------------------------------------
  message('\nExecuting stationary EVA on extreme-point subset')
  alphaCI=0.68
  Tr <- c(5, 10, 20, 50, 100, 200, 500, 1000)
  # Override the POT threshold with the subset minimum
  EVdata=list()
  ik <- 1
  th = gpdThreshold
  d1 <- stationary_peaks
  fit <- suppressWarnings(try(POT::fitgpd(d1, threshold = th, 
                                          est = "mle", method = "BFGS", std.err.type = "expected"), 
                              TRUE))
  if (!inherits(fit, "try-error")) {
    ksi <- fit$par[2]
    sgm <- fit$par[1]
    alphaCIx = 1 - alphaCI
    probs <- c(alphaCIx/2, 1 - alphaCIx/2)
    kci <- try(stats::qnorm(probs, ksi, fit$std.err[2]), 
               silent = T)
    kci[kci < -1] <- -1
    lnsigci <- try(stats::qnorm(probs, log(sgm), fit$std.err[1]/sgm))
    paramCIs <- cbind(kci, sigci = exp(lnsigci))
    paramEstsall <- c(sgm, ksi, gpdThreshold, 
                      length(d1), length(stationary_peaks), percentile)
    rlvls <- gpdThreshold + (sgm/ksi) * ((((length(d1)/length(stationary_peaks)) * 
                                             (1/Tr))^(-ksi)) - 1)
    EVdata$GPDstat <- list(method = "GPD", values = rlvls, 
                           parameters = paramEstsall, paramCIs = paramCIs)
  }else {
    methodname <- "No fit"
    ik <- 1
    th = gpdThreshold
    d1 <- stationary_peaks
    paramEstsall <- c(NA, NA, 
                      gpdThreshold, length(d1), length(stationary_peaks), 
                      percentile)
    EVdata$GPDstat <- list(method = methodname, values = NA, 
                           parameters = paramEstsall, paramCIs = NA)
    message("could not estimate GPD: bounded distribution")
  }
  eva <- EVdata

    if (isFALSE(eva$isValid)) {
      message("Problem in the computation of EVA statistics")
    }
  
  # --- 8. Back-transform GPD parameters to non-stationary space -------------
  if (eva[[1]]$method != "No fit") {
    
    epsilonPotX    <- eva[[1]]$parameters[2]
    errEpsilonPotX <- epsilonPotX - eva[[1]]$paramCIs[1, 1]
    sigmaPotX      <- eva[[1]]$parameters[1]
    errSigmaPotX   <- sigmaPotX   - eva[[1]]$paramCIs[1, 2]
    thresholdPotX  <- eva[[1]]$parameters[3]   # == gpdThreshold in stationary space
    errThresholdPotX <- eva[[1]]$thresholdError
    nPotPeaks      <- eva[[1]]$parameters[5]
    percentilePotX <- eva[[1]]$parameters[6]
    
    minPeakDistance <- minPeakDistanceInDays / dtn
    tsDate  <- as.Date(timeStamps)
    dtPotX  <- as.numeric(tsDate[length(tsDate)] - tsDate[1]) /
      length(newSeries) * minPeakDistance
    
    epsilonPotNS <- epsilonPotX
    errEpsilonPotNS <- errEpsilonPotX
    
    sigmaPotNS <- sigmaPotX * std_dev_series
    errSigmaPotFit   <- std_dev_series * errSigmaPotX
    errSigmaPotTransf <- sigmaPotX * adjustedStdDevError
    errSigmaPotNS    <- sqrt(errSigmaPotTransf^2 + errSigmaPotFit^2)
    
    # Threshold back-transform: same rule as mu (trend + stdDev scaling)
    thresholdPotNS <- thresholdPotX * std_dev_series + trendSerieX
    errThresholdPotX=0
    thresholdErrTransf <- sqrt(
      (std_dev_series  * errThresholdPotX)^2 +
        (thresholdPotX   * adjustedStdDevError)^2 +
        adjustedTrendError^2
    )
    thresholdErr <- thresholdErrTransf
    
    potParams <- list(
      epsilon          = epsilonPotNS,
      sigma            = sigmaPotNS,
      threshold        = thresholdPotNS,
      percentile       = percentilePotX,
      timeDelta        = dtPotX,
      timeDeltaYears   = dtPotX / 365.25,
      timeHorizonStart = min(timeStamps),
      timeHorizonEnd   = max(timeStamps),
      peaks            = Peaks_1,
      peakID           = Peaks_tid,
      peakST           = Peaks_loc$tIDstart,
      peakEN           = Peaks_loc$tIDend,
      nPeaks           = nPotPeaks
    )
    
    potParamStdErr <- list(
      epsilonErr        = errEpsilonPotNS,
      sigmaErrFit       = errSigmaPotFit,
      sigmaErrTransf    = errSigmaPotTransf,
      sigmaErr          = errSigmaPotNS,
      thresholdErrFit   = 0,
      thresholdErrTransf = thresholdErrTransf,
      thresholdErr      = thresholdErr
    )
    
    potObj <- list(
      method           = eva[[1]]$method,
      parameters       = potParams,
      paramErr         = potParamStdErr,
      stationaryParams = eva[[1]],
      objs             = NULL
    )
    
  } else {
    
    # No-fit fallback for GPD
    minPeakDistance <- minPeakDistanceInDays / dtn
    tsDate  <- as.Date(timeStamps)
    dtPotX  <- as.numeric(tsDate[length(tsDate)] - tsDate[1]) /
      length(newSeries) * minPeakDistance
    
    thresholdPotX  <- gpdThreshold   # == gpdThreshold
    epsilonPotX    <- NA
    sigmaPotX      <- NA
    
    potParams <- list(
      epsilon          = epsilonPotX,
      sigma            = sigmaPotX  * std_dev_series,
      threshold        = thresholdPotX * std_dev_series + trendSerieX,
      percentile       = percentile,
      timeDelta        = dtPotX,
      timeDeltaYears   = dtPotX / 365.2425,
      timeHorizonStart = min(timeStamps),
      timeHorizonEnd   = max(timeStamps),
      peaks            = Peaks_1,
      peakID           = Peaks_tid,
      peakST           = Peaks_loc$tIDstart,
      peakEN           = Peaks_loc$tIDend,
      nPeaks           = nPotPeaks
    )
    
    potObj <- list(
      method           = "No fit",
      parameters       = potParams,
      paramErr         = NULL,
      stationaryParams = NULL,
      objs             = NULL
    )
  }
  
  # --- 9. Return (same structure as TsEvaNs) ---------------------------------
  nonStationaryEvaParams <- list(potObj = potObj)
  
  return(list(
    nonStationaryEvaParams = nonStationaryEvaParams,
    stationaryTransformData = list(
      timeStamps       = timeStamps,
      stationarySeries = stationarySeries,
      trendSeries      = trendSerieX,
      stdDevSeries     = std_dev_series,
      stdDevError      = adjustedStdDevError,
      trendError       = adjustedTrendError,
      nonStatSeries    = series
    )
  ))
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
  idfix=2
  start_time <- Sys.time()
  print(paste0("hazard:",haz," square: ", Nsq, " pixel: ",idfix,"/",endid))
  catch=as.numeric(unikout[idfix])
  Peaks_loc=Peaksave[which(Peaksave$catch==catch),]
  Params_loc=Paramsfl[which(Paramsfl$catchment==catch),]
  timeStamps=txx
  thresh=thresh_vec[which(thresh_vec$cid==catch),]
  thresh=thresh$th
  frosttime=NA
  interflag=NA
  #do the analysis of trend and transfo for the SCF
  
  df.dis0=disNcopenloc(filenameSCF,hydroDir,outhybas,idfix)
  series0=data.frame(txx,df.dis0$outlets)
  names(series0)=c("date","Qs")
  rmv=which(as.integer(format(series0$date, "%Y"))==1950)
  if (length(rmv)>0){
    series0=series0[-rmv,]
  }
  #remove first day to avoid errors
  series0=series0[-1,]
  if (haz=="drought"){
    trans="rev"
    #seasonal split
    catmat=Catf7[which(Catf7$outlets==catch),]
    Tcatmat=mycat[which(mycat$HYBAS_ID==catmat$HYBAS_ID),]
    Tcatchment=which(colnames(frostcat)==Tcatmat$pointid)
    
    intermit=interid(series0,trans,WindowSize=7)
    interflag=intermit$flags[2]
    series0=data.frame(series0$date,intermit$trdis$Q7)
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
    series0 <- max_daily_value(series0)
    tail="high"
    trans="ori"
  }
  
  names(series0)=c("timestamp","dis")
  dt1=min(diff(series0$timestamp),na.rm=T)
  dt=as.numeric(dt1)
  tdim=attributes(dt1)$units
  if (tdim=="hours") dt=dt/24
  if (dt==1){
    timeDays=series0$timestamp
  }else{
    timeDays=unique(as.Date(series0$timestamp))
  }
  
  bounds=c(year(timeDays[1]),year(timeDays[length(timeDays)]))
  realbound=bounds
  tbound=c(as.Date(paste0(realbound[1],"-12-31")),as.Date(paste0(realbound[2],"-12-31")))
  Impdates=seq(tbound[1],tbound[2],by="1 year")
  
  nv=length(unique(series0$dis))
  if(length(na.omit(series0$dis))>1 & interflag<3 & nv>15){
    
    if (length(which(is.na(series0$dis)))>0){
      print("Na alert")
      seriefill=tsEvaFillSeries(series0$timestamp,series0$dis)
      series0$dis=seriefill
    }
    timeAndSeries0=series0
    names(timeAndSeries0)=c("timestamp","data")
    
    if (haz=="drought" & length(!is.na(frosttime))>1){
      if (season=="nonfrost"){
        print("nonfrost season")
        timeAndSeries0$data[frosttime]=NA
      }else if (season=="frost"){
        print("frost season")
        timeAndSeries0$data[-frosttime]=NA
      }else if (season=="year"){
        print("no seasonal divide")
      }else {print("season must be frost or nonfrost")}
    }else{
      print("no frost season for this river")
    }
    series0=timeAndSeries0[,2]
    
    timeWindow = 365.25*30; #time windows in days, the correction is done within the functions
    windowSize=366
    
    cat(paste0("\nsquare: ", Nsq, " pixel: ",idfix,"/",endid))
    #just compute trend with trendpeaks
    
    timeStamps=as.POSIXct(timeAndSeries0[,1])
    dt1=min(diff(timeStamps),na.rm=T)
    dt=as.numeric(dt1)
    tdim=attributes(dt1)$units
    if (tdim=="hours") dt=dt/24
    if (tdim=="seconds") dt=dt/3600

    if (tail=="high") epy=3
    if (tail=="low") epy=2
    
    #If the biggest event is more than 100 times greater than the second biggest event
    sanitycheck = computeAnnualMaxima(timeAndSeries0);
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
    
    stationaryTransformData = c()
    
    minEventsPerYear = 1
    if (tail=="low"){
      #default 7-day flow for low flow, can be modified by user
      start_index=1
      lowdt=7
      indices_to_extract <- seq(from = start_index, to = length(series0), by = lowdt/dt)
      series0=series0[indices_to_extract]
      timeStamps=timeStamps[indices_to_extract]
      shape_bnd=c(-1,0)
      minEventsPerYear = 1
      if (trans=="rev"){
        series0=-1*series0
      }else if(trans=="inv"){
        series0=1/series0
      }else if (trans=="lninv"){
        series0=-log(series0)
      }
    }
      Peaks_0=series0[Peaks_loc$timeID][order(Peaks_loc$timeID)]
      Peaks_0t=timeStamps[Peaks_loc$timeID][order(Peaks_loc$timeID)]
      #plot(Peaks_0t,Peaks_0)
      TrendTh=thresh
      print(TrendTh)
      message(paste0('\nevaluating long term variations of the peaks'))
      if (is.na(TrendTh)){
        TrendTh=try(tsEvaFindTrendThreshold(series0, timeStamps, timeWindow),T)
        if (inherits(TrendTh, "try-error") || length(TrendTh) == 0 || all(is.na(TrendTh))) {
          TrendTh <- 0.1
        }
        trasfData = tsEvaTransformSeriesToStationaryPeakTrend( timeStamps, series0, timeWindow, TrendTh);
      }else{
          trasfData = tsEvaTransformSeriesToStationaryPeakTrend( timeStamps, series0, timeWindow, TrendTh);
          c=0
          while ((length(unique(trasfData$trendSeries))<2)) {
            trasfData = tsEvaTransformSeriesToStationaryPeakTrend( timeStamps, series0, timeWindow, TrendTh=TrendTh-c);
            c=c+0.1
          }
          TrendTh=TrendTh+c
        }
      TrendTh
      #plot(series0)
      qd <- stats::quantile(series0, TrendTh, na.rm = T)
      series_above_threshold <- series0[which(series0 >= qd)]
      ids_above_threshold=which(series0 >= qd)
      
    gevMaxima = 'annual'
    potEventsPerYear = epy

    
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
    
  ##new serie --------
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
  # if(length(na.omit(series$dis))>1 & interflag<3 & nv>15){
    
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
    if (tail=="low"){
      #default 7-day flow for low flow, can be modified by user
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
    Peaks_1=series[Peaks_loc$timeID][order(Peaks_loc$timeID)]
    Peaks_tid=Peaks_loc$timeID
    Peaks_t=timeStamps[Peaks_loc$timeID][order(Peaks_loc$timeID)]
    plot(Peaks_t,Peaks_1)
    points(Peaks_t,Peaks_0,col=2, pch=3)
    plot(Peaks_1-Peaks_0)
    cat(paste0("\nsquare: ", Nsq, " pixel: ",idfix,"/",endid))
    extremePoints=data.frame(timeStamps[ids_above_threshold],series[ids_above_threshold])
    series_above_threshold=series
    series_above_threshold[-ids_above_threshold]=NA
    detrend_result <- tsEvaDetrendTimeSeries(timeStamps, series_above_threshold, 
                                             timeWindow)
    plot(detrend_result@trendSeries)
    detrended_series <- series - detrend_result@trendSeries
    detrended_series_above_threshold <- series_above_threshold - 
      detrend_result@trendSeries
    n_run_mean <- detrend_result@nRunMn
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

    
    originalPoints=data.frame(timeStamps[ids_above_threshold],series0[ids_above_threshold])
    
    
    plot(extremePoints$series.ids_above_threshold.)
    points(originalPoints$series0.ids_above_threshold.,col=2)
    points(originalPoints,col=2)
    lines(trasfData$timeStamps,trasfData$trendSeries)
    plot(extremePoints,ylim=c(-0.1,-0.01))
    points(originalPoints, col=2)
    plot(trasfData$timeStamps, trasfData$trendSeries,col=2)
    lines(timeStamps,trendS)
    lines(transformed_data$timeStamps, transformed_data$trendSeries,col=1)
    
   
    plot(transformed_data$timeStamps, transformed_data$stdDevSeries,col=1)
    lines(trasfData$timeStamps, trasfData$stdDevSeries,col=2)
    
    #plot(extremePoints,originalPoints)
    Nonstat<-TsEvaNs_fromExtremes(    timeStamps,
                                      extremePoints,        # data.frame: col1=timestamps, col2=extreme values (pre-filtered)
                                      originalPoints,       # data.frame: col1=timestamps, col2=original series values at same indices
                                      series,
                                      Peaks_0,
                                      Peaks_1,
                                      Peaks_tid,
                                      Peaks_loc,
                                      trasfData,            # output$stationaryTransformData from a previous TsEvaNs call
                                      timeWindow,
                                      minPeakDistanceInDays = 7,
                                      gevMaxima = 'annual',
                                      gevType = 'GEV',
                                      evdType = c('GEV', 'GPD'),
                                      tail = "high",
                                      epy = 3,
                                      shape_bnd = c(-0.5, 1))
      #TsEvaNs(timeAndSeries, timeWindow, transfType='trendPeaks',ciPercentile = ciPercentile, minPeakDistanceInDays = minPeakDistanceInDays,lowdt=7,trans=trans,tail = tail,TrendTh = thresh)
    nonStationaryEvaParams=Nonstat[[1]]
    stationaryTransformData=Nonstat[[2]]
    
    stationaryTransformData$timeStampsDay=unique(as.Date(stationaryTransformData$timeStamps))
    pikos=data.frame(nonStationaryEvaParams$potObj$parameters$peaks,
                     nonStationaryEvaParams$potObj$parameters$peakID,
                     nonStationaryEvaParams$potObj$parameters$peakST,
                     nonStationaryEvaParams$potObj$parameters$peakEN)
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
    RLevs100=ComputeReturnLevelsGPD(nonStationaryEvaParams, RPgoal, timeIndex)
    if (RLevs100$Fit=="No fit"){
      
      RLgpd=computeAnnualMaxima(timeAndSeries)[[1]]
      names(RLgpd)=year(Impdates)
      
      RLgev=computeAnnualMaxima(timeAndSeries)[[1]]
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
      RLgpd=RLevs100$ReturnLevels[2]
      ERgpd=RLevs100$ReturnLevels[3]
      
      nRPgev=nRPgpd=10
      params=c()
      for (t in 2:length(Impdates)){
        timeIndex=tindexes[t]
        RLevs100i=ComputeReturnLevelsGPD(nonStationaryEvaParams, RPgoal, timeIndex)
        params=c(catch,year(Impdates[t]),timeIndex,RLevs100i$Params)
        names(params)[1:3]=c("catchment","Year","timeIndex")
        
        Rper=RPcalcGPD(params,RPiGPD=RLevs100$ReturnLevels[2])
        nRPgpd=c(nRPgpd,Rper[2])
   
        RLgpd=cbind(RLgpd,RLevs100i$ReturnLevels[2])
        ERgpd=cbind(ERgpd,RLevs100i$ReturnLevels[3])
        if (length(parlist)>1) colnames(parlist)=names(params)
        parlist=rbind(parlist,params)
      }
      
      # RLgev=as.data.frame(RLgev)
      # names(RLgev)=year(Impdates)
      # rownames(RLgev)=RPgoal
      
      RLgpd=as.data.frame(RLgpd)
      names(RLgpd)=year(Impdates)
      rownames(RLgpd)=RPgoal
      
      # nRPgev=as.data.frame(t(nRPgev))
      # names(nRPgev)=year(Impdates)
      # 
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
    
    params=data.frame(matrix(ncol=11,nrow=(length(Impdates)-1)))
    params[,1]=rep(catch,length(Impdates)-1)
    params[,2]=year(Impdates)[-1]
    params[,4]=rep(interflag,length(Impdates)-1)
    if (is.null(colnames(parlist))){
      colnames(params)=rep("nom",11)
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
  
  #RetLevGEV=rbind(RetLevGEV,RLgev)
  
  RetLevGPD=rbind(RetLevGPD,RLgpd)
  
  #RetPerGEV=rbind(RetPerGEV,nRPgev)
  
  RetPerGPD=rbind(RetPerGPD,nRPgpd)
  end_time <- Sys.time()
  cat(paste0("\nloop duration: ",round(end_time-start_time,2)," seconds\n"))
  
}

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







Results=list(parameters=parlist,RetLevGEV=RetLevGEV,RetLevGPD=RetLevGPD,RetPerGEV=RetPerGEV,RetPerGPD=RetPerGPD,Peaks=peaklist,catrest=data.frame(catlist,IRES))


tsGetPOT_2 <- function(ms, pcts, desiredEventsPerYear,minEventsPerYear, minPeakDistanceInDays, tail) {
  
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
        thrsdt <- stats::quantile(ms[,2],pcts[ipp]/100,na.rm=T)
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
            devpp[ipp]=1e9
            gpp[ipp]=9999
          } else {
            gpdpar = fgpd$fitted.values # [scale, shape]
            
            # 1. Calculate the Cumulative Distribution Function (CDF) values for the peaks
            # Using the GPD formula: F(x) = 1 - (1 + shape * (x-thresh)/scale)^(-1/shape)
            scaled_peaks = (pks[,1] - thrsdt) / gpdpar[1]
            if(abs(gpdpar[2]) < 1e-10) { # Handle case where shape is nearly 0 (Exponential)
              z = 1 - exp(-scaled_peaks)
            } else {
              z = 1 - (1 + gpdpar[2] * scaled_peaks)^(-1/gpdpar[2])
            }
            
            # 2. Calculate Right-Tail Weighted Anderson-Darling (ADR)
            # This version emphasizes the upper tail of the distribution
            z = sort(z)
            n = length(z)
            i = 1:n
            # ADR formula: -n/2 - 2*sum(z) - sum((2*i-1)*log(1-z))/n (simplified version)
            # Note: We use a version that targets the right tail specifically
            adr_stat = -n - (1/n) * sum((2*i - 1) * log(z) + (2*n + 1 - 2*i) * log(1 - z))
            
            # 3. Store the statistic plus your existing performance penalties
            devpp[ipp] = adr_stat + perfpen 
            gpp[ipp] = gpdpar[2]
          }
          
          nperYear <- tsGetNumberPerYear(ms, pks[,2])
          minnumperyear[ipp] <- min(nperYear$Freq, na.rm = TRUE)
          
          # if(numperyear[ipp]<=desiredEventsPerYear+1 & dej==0){
          # fgpd=suppressWarnings(try(POT::fitgpd(pks[,1], threshold = thrsdt, est = "mle",method="BFGS",std.err.type = "expected")))
          # if(inherits(fgpd, "try-error")){
          # gpdpar=9999
          # deviance=9999
          # devpp[ipp]=1e9
          # gpp[ipp]=9999
          # }else {
          # gpdpar=fgpd$fitted.values
          # deviance=fgpd$deviance
          # devpp[ipp]=stats::AIC(fgpd)+perfpen
          # gpp[ipp]=gpdpar[2]
          # }
          # nperYear <- tsGetNumberPerYear(ms, pks[,2])
          # minnumperyear[ipp] <- min(nperYear$Freq, na.rm = TRUE)
        }
      }
    }
  }
}







