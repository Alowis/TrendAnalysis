
#Library calling ----
required_packages <- c(
  "ncdf4", 
  "sf", 
  "rnaturalearth", 
  "rnaturalearthdata", 
  "rgeos", 
  "dplyr", 
  "RtsEva", 
  "evd",
  "POT",
  "lubridate", 
  "zoo",
  "fs",
  "Kendall", 
  "biscale", 
  "cowplot", 
  "ggpubr", 
  "ggridges", 
  "ggplot2", 
  "viridis", 
  "hrbrthemes", 
  "tidyverse", 
  "raster", 
  "modifiedmk", 
  "ks", 
  "xts",
  "pracma", 
  "data.table", 
  "matrixStats", 
  "tsibble",
  "scales",
  "foreach", 
  "doParallel",
  "biscale",
  "ggnewscale",
  "pracma",
  "data.table",
  "matrixStats"
)

#Install packages if missing and load them
for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    suppressMessages(install.packages(pkg))
  }
  if (!(pkg %in% .packages())) {
    suppressPackageStartupMessages(library(pkg, character.only = TRUE))
  }
}


#1 Miscellaneous Functions ---------------------------------------------------

#1.1 Functions for trend threshold selection

# Checks if there are not gaps bigger than two years in extreme value timeseries
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

# Finds the optimal extreme threshold for trend estimation
tsEvaFindTrendThreshold2<-function(series, timeStamps, timeWindow){
  ptn = timeStamps[which(!is.na(series))]
  bounds = unique(lubridate::year(ptn))
  nr <- rep(1, length(series))
  #initialization of the normalized trend (no trend)
  nr = nr + rnorm(length(series), 0, 1e-05)
  sts <- c()
  lnegs = c()
  pctd = c()
  pcts <- seq(0.4, 0.95, by = 0.05)
  #iteration over threshold values
  for (iter in 1:length(pcts)) {
    thrsdt <- quantile(series, pcts[iter], na.rm = TRUE)
    series_no_na <- series
    series_no_na[which(is.na(series_no_na))] <- -9999
    serieb <- series_no_na
    timeb = timeStamps
    timeb = timeb[-which(serieb < thrsdt)]
    serieb[which(serieb < thrsdt)] <- NA
    checkY=check_timeserie2(timeb,bounds)
    #if a gap in extreme is larger than two years, stop the loop
    if (checkY == FALSE) {
      print(paste0("not all years - q= ",pcts[iter]))
      break
    }
    #compute the trend on extremes
    rs <- tsEvaDetrendTimeSeries(timeStamps, serieb, timeWindow, 
                                 fast = T)
    #normalize the trend
    norm_trend <- rs@trendSeries/mean(rs@trendSeries, na.rm = TRUE)
    
    #identification of trends leading to negative flow values (for low flows)
    dtr1 = serieb - rs@trendSeries
    lneg = length(which(dtr1 < 0))
    #stability of the trend between different threshold values
    stab <- cor(nr, norm_trend, use = "pairwise.complete.obs")
    if (iter == 1) nr <- norm_trend

    lnegs = c(lnegs, lneg)
    sts <- c(sts, stab)
    pctd = c(pctd, pcts[iter])
  }
  dow=abs(diff(sts))[-1]
  
  #by default, the selected threshold is the highest possible value
  rval = pctd[length(pctd)]
  
  #if trend abruptly changes, select threshold before the abrupt change
  if(max(dow)>0.2){
    print("breaking point")
    rval = pctd[which.max(dow)+1]
  }
  #if trend creates negative flow value, select the threshold with the least negative values
  if (sum(lnegs) > 1) {
    rval = pctd[which.min(lnegs)]
  }
  
  return(rval)
}

# open netcdf outlet files
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

#compute return levels associated to a RP (not used in this script)
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

#obtain RL for any return period
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

#Open reservoir location file
ReservoirOpen=function(dir,outletname,Sloc_final){
  ncbassin=paste0(dir,outletname)
  ncb=nc_open(ncbassin)
  name.vb=names(ncb[['var']])
  namev=name.vb[1]
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

#Extract discharge timeseries for selected locations
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

#Extract discharge timeserie for one location
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

#Detection of intermittent rivers
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

#2  Univariate analysis  ---------------------------------------------------

calculate_return_levels <- function(ParamSpecial, XX = NULL,ci=1) {
  # If XX is not provided, calculate it based on nPeaks
  if (is.null(XX)) {
    X0 <- ParamSpecial$nPeaks / 70
    XX <- X0 * 10
  }
  
  # Extract parameters from ParamSpecial
  npars <- length(ParamSpecial$sigmaGPD)
  nt <- 1
  sigma_ <- ParamSpecial$sigmaGPD
  sigmaStdErr_ <- ParamSpecial$sigmaStdErrGPD
  threshold_ <- ParamSpecial$thresholdGPD
  thresholdStdErr_ <- ParamSpecial$thresholdStdErrGPD
  epsilon <- ParamSpecial$epsilonGPD
  epsilonStdErr <- ParamSpecial$epsilonStdErrGPD
  
  # Calculate return levels
  returnLevels <- threshold_ + sigma_ / epsilon * ((XX)^epsilon - 1)
  
  # Calculate derivatives
  dxm_u <- 1
  dxm_sigma <- 1 / epsilon * (XX^epsilon - 1)
  dxm_epsilon <- -sigma_ / epsilon^2 * ((XX)^epsilon - 1) + 
    sigma_ / epsilon * log(XX) * XX^epsilon
  
  # Calculate return level errors
  returnLevelsErr <- sqrt((dxm_u * ci*thresholdStdErr_)^2 + 
                            (dxm_sigma * ci*sigmaStdErr_)^2 + (dxm_epsilon * ci*epsilonStdErr)^2)
  
  # Replace large errors with NA
  returnLevelsErr[which(returnLevelsErr > 10000)] <- NA
  
  # Return a data frame
  data.frame(returnLevelErr = returnLevelsErr, returnLevels = returnLevels)
}

weighted_average <- function(r,point, points, rivermask, max_distance) {
  # Calculate distances between the point of interest and all other points
  distance3 <- spDists(points, point, longlat=T)
  # Calculate weights as the inverse of the distance
  weights <- 1 / (10 + (distance3))
  # Set weights to 0 for points beyond the maximum distance
  weights[distance3 > max_distance] <- 0
  
  # Normalize weights
  weights <- weights / sum(weights,na.rm=T)
  # Calculate the weighted average
  weighted_avg <- sum(weights * points$z, na.rm=T)
  return(weighted_avg)
}

neighbour_finder <- function(point, points, max_distance) {
  
  distance3 <- spDists(points, point, longlat=T)
  neighbours=points$z[which(distance3 < max_distance)]
  return(neighbours)
}

ComputeChange<- function(Drivertrend, unikout, DataI, outhybas07, 
                         parameters, rmpixels, UpAvec, GNF, Regio, yrname, change, eps=0.1) {
  
  # Data preparation
  data <- data.frame(Drivertrend, unikout = unikout)
  rmp2 <- na.omit(unique(parameters$catchment[rmpixels]))
  data <- data[-match(rmp2, data$unikout), ]
  check=rowMeans(data[,-71])
  pb=which(abs(check)>10000)
  ouy=data[c(pb[1:length(pb)]),]
  data[pb,]=NA
  
  
  data <- right_join(data, UpAvec, by = c("unikout" = "outl2"))
  DataI <- right_join(DataI, UpAvec, by = c("unikout" = "outl2"))
  # Aggregation choice
  DataC <- right_join(GHR_riv, data, by = c("outl2" = "unikout"))
  DataC=data.frame(DataC)
  length(which(is.na(DataC$Var2)))
  # Matching NUTS3 and NUTS2 IDs
  HRM <- match(DataC$HydroRegions_raster_WGS84, Regio$Id)
  BRM <- match(DataC$outl2, biogeo_rivers$outl2)
  # N2M <- match(DataC$NUTS3_Raster2ID, NUTS3$N2ID)
  
  # Adding NUTS3 and NUTS2 IDs to the data
  DataC$Regio_id <- Regio$Id[HRM]
  DataC$Biogeo_id <- biogeo_rivers$code[BRM]
  # DataC$NUTS2_id <- NUTS3$NUTS2_ID[N2M]
  yrange=match(yrname,colnames(DataC))
  
  
  # Normalize the data by area (upa)
  DataX=DataC
  DataC[, yrange] <- DataC[, yrange] * 1000 / DataC$upa
  
  #this line to go to relative differences
  irange=match(yrname,colnames(DataI))
  #important line for control of relative difference
  eps=rep(eps,length(DataI$unikout))
  
  
  #matching the two tables:
  matx=match(DataX$outl2,DataI$unikout)
  DataI=DataI[matx,]
  
  DataX[,yrange]=DataX[,yrange]/(DataI[,irange]+eps)*100
  # min(DataX$Y2015,na.rm=T)
  return(data=DataX)
}

calculatePoints <- function(trendPlot, yrlist, pointagg, Regio, GHshpp, datap) {
  
  # Initialize empty lists for storing results
  mkta <- c()
  mksa <- c()
  chlist <- c()
  
  tmpval <- trendPlot[,-1]
  
  for (it in 1:length(pointagg[,1])) {
    # print(it)
    mks <- NA
    mkt <- NA
    miniTS <- as.numeric(tmpval[it,])
    
    # Compute the differences
    dmt <- diff(miniTS)
    
    # Find points where the trend changes direction
    sign_change <- which(diff(sign(dmt)) != 0)
    sign_change <- sign_change + 2
    
    # Remove changes that occur over periods shorter than 3 years
    dsc <- diff(sign_change)
    rms <- which(dsc < 3)
    sign_change <- sign_change[-c(rms, (rms + 1))]
    
    # Determine the direction of change and corresponding years
    dirchange <- sign(c(dmt[1], dmt[sign_change - 1]))
    yrchange <- c(yrlist[1], yrlist[sign_change])
    
    # Perform Mann-Kendall test if data exists and positive trend is detected
    if (!is.na(miniTS[2]) & max(abs(diff(miniTS[-1])), na.rm = TRUE) > 0) {
      mk2 <- mmkh(miniTS[-1], ci = 0.95)
      mk <- data.frame(tau = mk2[7], sl = mk2[2])
      
      mkt <- mk$tau
      mks <- mk$sl
      #print(mkt)
    } else {
      mkt <- NA
      mks <- NA
    }
    
    mkta <- c(mkta, mkt)
    mksa <- c(mksa, mks)
    
    # Store change information
    changes <- data.frame(rep(it, length(dirchange)), dirchange, yrchange)
    chlist <- rbind(chlist, changes)
  }
  
  # Add Mann-Kendall test results to pointagg
  pointagg$mkta <- mkta
  pointagg$sl <- mksa
  
  # Merge with NUTS3 data to get spatial points
  pointsag <- inner_join(Regio, pointagg, by = c("CODEB" = "HydroR"))
  
  # Assign significance levels
  pointsag$siglvl <- 0
  pointsag$siglvl[which(pointsag$sl <= 0.05)] <- 1
  
  # Define change categories
  pointsag$change <- 0
  pointsag$change[which(pointsag$Rchange_rel.mean > 0)] <- 1
  pointsag$change[which(pointsag$Rchange_rel.mean < 0)] <- -1
  pointsag$change[which(pointsag$Rchange_rel.mean > 0 & pointsag$siglvl > 0)] <- 2
  pointsag$change[which(pointsag$Rchange_rel.mean < 0 & pointsag$siglvl > 0)] <- -2
  
  
  
  # Assign slope values
  pointsag$tslop <- pointsag$mkta
  pointsag$tslop[which(pointsag$siglvl < 1)] <- NA
  
  # Filter significant points and convert to spatial data
  catsig <- pointsag[which(pointsag$siglvl > 0), ]
  catsig <- st_transform(catsig, crs = 3035)
  
  # Create grid points
  pointsInside=NA
  if (length(catsig$Id)>0){
    grdpts <- sf::st_make_grid(catsig, what = "centers", cellsize = 30000)
    my.points <- sf::st_sf(grdpts)
    sf_use_s2(FALSE)
    
    # Find points inside the significant categories
    pointsInside <- sf::st_join(x = my.points, y = catsig, left = FALSE)
    pointsInside$sign <- "increase"
    pointsInside$sign[which(pointsInside$mkta <= 0)] <- "decrease"
  }
  # Create pagg (significant points data)
  pagg <- catsig
  pagg$sign <- "increase"
  pagg$sign[which(pagg$tslop < 0)] <- "decrease"
  
  # Create final points
  if (length(is.na(datap$Var1))>0) datap=datap[-which(is.na(datap$Var1)),]
  points <- st_as_sf(datap, coords = c("Var1", "Var2"), crs = 4326)
  points <- st_transform(points, crs = 3035)
  
  return(list(PagD=pointsag,points=points,psp=pointsInside))
}

processTrendData <- function(trendData, DataTr, id_var = "HydroR") {
  trtest <- suppressWarnings(melt(trendData, id.vars = id_var, variable.name = "variable", value.name = "value"))
  
  # Convert variable to character and extract the year
  trtest$variable <- as.character(trtest$variable)
  craplife <- data.frame(strsplit(trtest$variable, ".Y"))
  trtest$yr <- as.numeric(craplife[2,])
  
  
  decades=c(1955,1965,1975,1985,1995,2005,2015)
  
  valcol=which(!is.na(match(trtest$yr,decades)))
  # Calculate the decade
  trtest=trtest[valcol,]
  trtest$decad=trtest$yr
  
  # Aggregate by decade and location
  trtime <- aggregate(list(value = trtest$value),
                      by = list(yr = trtest$decad, loc = trtest[[id_var]]),
                      FUN = function(x) c(mean = mean(x, na.rm = TRUE),
                                          l = length(x),
                                          ql = quantile(x, 0.035, na.rm = TRUE),
                                          qh = quantile(x, 0.975, na.rm = TRUE)))
  
  # Convert to a data frame
  tData <- do.call(data.frame, trtime)
  names(tData)[c(3, 4, 5, 6)] <- c("changeC", "length", "cq1", "cq2")
  
  #Aggregation for biogeoregions
  
  decades=c("Y1955","Y1965","Y1975","Y1985","Y1995","Y2005","Y2015")
  valcol=match(decades,colnames(DataTr))
  wc=match("Biogeo_id",colnames(DataTr))
  valcol=c(wc,valcol)
  dfd=DataTr[,valcol]
  
  dftest <- suppressWarnings(melt(dfd, id.vars = "Biogeo_id", variable.name = "variable", value.name = "value"))
  
  dftime <- aggregate(list(value = dftest$value),
                      by = list(yr = dftest$variable, loc = dftest[["Biogeo_id"]]),
                      FUN = function(x) c(mean = mean(x, na.rm = TRUE),
                                          l = length(x),
                                          med= median(x, na.rm=T),
                                          ql = quantile(x, 0.25, na.rm = TRUE),
                                          qh = quantile(x, 0.75, na.rm = TRUE),
                                          w1 = quantile(x, 0.025, na.rm = TRUE),
                                          w2 = quantile(x, 0.975, na.rm = TRUE)))
  BgData <- do.call(data.frame, dftime)
  BgData$decad=seq(1950,2010,by=10)
  
  
  # Aggregate only by decade (without location)
  trtime_global <- aggregate(list(value = trtest$value),
                             by = list(yr = trtest$decad),
                             FUN = function(x) c(mean = mean(x, na.rm = TRUE),
                                                 l = length(x),
                                                 med= median(x, na.rm=T),
                                                 ql = quantile(x, 0.25, na.rm = TRUE),
                                                 qh = quantile(x, 0.75, na.rm = TRUE),
                                                 w1 = quantile(x, 0.025, na.rm = TRUE),
                                                 w2 = quantile(x, 0.975, na.rm = TRUE)))
  
  # Convert to a data frame
  tGlobal <- do.call(data.frame, trtime_global)
  # Return both the processed data with location (tData) and the global trend (tGlobal)
  return(list(tData = tData, tGlobal = tGlobal, BgData=BgData))
}

UpATrendData <- function(DataTr, id_var = "upagroup") {
  decades=c("Y2015")
  
  DataTr$upagroup=1
  DataTr$upagroup[which(DataTr$upa>200 & DataTr$upa<=500)]=2
  DataTr$upagroup[which(DataTr$upa>500 & DataTr$upa<=1000)]=3
  DataTr$upagroup[which(DataTr$upa>1000 & DataTr$upa<=10000)]=4
  DataTr$upagroup[which(DataTr$upa>10000)]=5
  valcol=match(decades,colnames(DataTr))
  vi=match(id_var,colnames(DataTr))
  valcol=c(vi,valcol)
  dfd=DataTr[,valcol]
  
  dftest <- suppressWarnings(melt(dfd, id.vars = id_var, variable.name = "variable", value.name = "value"))
  
  dftime <- aggregate(list(value = dftest$value),
                      by = list(yr = dftest$variable, loc = dftest[["upagroup"]]),
                      FUN = function(x) c(mean = mean(x, na.rm = TRUE),
                                          l = length(x),
                                          med= median(x, na.rm=T),
                                          ql = quantile(x, 0.25, na.rm = TRUE),
                                          qh = quantile(x, 0.75, na.rm = TRUE),
                                          w1 = quantile(x, 0.025, na.rm = TRUE),
                                          w2 = quantile(x, 0.975, na.rm = TRUE)))
  BgData <- do.call(data.frame, dftime)
  
  return(Cdata=BgData)
}


#3 Bivariate analysis -----

calculateTrendSig <- function(trendPlot, pointagg) {
  
  # Initialize empty lists for storing results
  mkta <- c()
  mksa <- c()
  
  tmpval <- trendPlot[,-1]
  
  for (it in 1:length(pointagg[,1])) {
    # print(it)
    mks <- NA
    mkt <- NA
    miniTS <- as.numeric((tmpval[it,]))
    miniTS=miniTS[-which(is.na(miniTS))]
    plot(miniTS)
    # Perform Mann-Kendall test if data exists and positive trend is detected
    if (!is.na(miniTS[2]) & max(abs(diff(miniTS[-1])), na.rm = TRUE) > 0) {
      mk2 <- mmkh(miniTS[-1], ci = 0.95)
      mk <- data.frame(tau = mk2[7], sl = mk2[2])
      
      mkt <- mk$tau
      mks <- mk$sl
      #print(mkt)
    } else {
      mkt <- NA
      mks <- NA
    }
    
    mkta <- c(mkta, mkt)
    mksa <- c(mksa, mks)
  }
  
  # Add Mann-Kendall test results to pointagg
  pointagg$mkta <- mkta
  pointagg$sl <- mksa
  
  pointagg$siglvl <- 0
  pointagg$siglvl[which(pointagg$sl <= 0.05)] <- 1
  
  # Define change categories
  pointagg$change <- 0
  pointagg$change[which(pointagg$Rchange.Y2020> 0)] <- 1
  pointagg$change[which(pointagg$Rchange.Y2020< 0)] <- -1
  pointagg$change[which(pointagg$Rchange.Y2020> 0 & pointagg$siglvl > 0)] <- 2
  pointagg$change[which(pointagg$Rchange.Y2020< 0 & pointagg$siglvl > 0)] <- -2
  
  # Assign slope values
  pointagg$tslop <- pointagg$mkta
  pointagg$tslop[which(pointagg$siglvl < 1)] <- NA
  
  return(pointagg)
}

# Define the get_density function
get_density <- function(x, y, n = 200) {
  dens <- MASS::kde2d(x = x, y = y, n = n)
  ix <- findInterval(x, dens$x)
  iy <- findInterval(y, dens$y)
  ii <- cbind(ix, iy)
  return(dens$z[ii])
}

