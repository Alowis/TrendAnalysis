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
library(data.table)


# THIS SCRIPT WILL BE A CLEAN AND SHORT VERSION OF THE INITIAL SCRIPT
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







#  Function declaration ---------------------------------------------------


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
plotchangemaps=function(basemap,catmap,datar,law="GPD",type,period=c(1950,2020),parlist,valuenames){
  
  
  data=right_join(catmap,datar,by = c("outlets"="unikout"))
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
      tmpval[,it][which(tmpval[,it]>800)]=NA
    }
    data[,valcol]=tmpval
    br=c(-100,-80,-60,-40,-20,0,20,40,60,80,100)
    limi=c(-100,100)
    trans=scales::modulus_trans(.8)
    colNA="transparent"
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
plotbichangemaps=function(basemap,catmap,datar1,datar2,law="GPD",type,period=c(1950,1990),parlist,valuenames){
  
  
  catmap$outlets=catmap$pointid
  data1=right_join(catmap,datar1,by = c("outlets"="unikout"))
  data2=right_join(catmap,datar2,by = c("outlets"="unikout"))
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
  mcor=unique(match(datar1$unikout,parlist$catchment))
  if (length(which(is.na(mcor)))>0) datar1=datar1[-which(is.na(mcor)),]
  
  cref=paste0("Y",period[1])
  crefloc=match(cref,dtc)
  finalperiod=paste0("Y",period[2])
  colsel=match(finalperiod,datacol)
  
  
  palet=c(hcl.colors(11, palette = "BrBG", alpha = NULL, rev = F, fixup = TRUE))
  title=paste0("Relative change in 100-years flood Return Level between ", period[1], " and ", period[2])
  legend="Relative Change (%)"
  tmpval1=(datatwin1[,valcol2])/(datatwin1[,crefloc])*100-100
  tmpval2=(datatwin2[,valcol2])/(datatwin2[,crefloc])*100-100
  for (it in 1:length(tmpval1[1,])){
    tmpval1[,it][which(tmpval1[,it]>200)]=NA
    tmpval2[,it][which(tmpval2[,it]>200)]=NA
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
  
  # combine map with legend
  finalPlot <- ggdraw() +
    draw_plot(map, 0, 0, 1, 1) +
    draw_plot(legend, 0.64, 0.7, 0.2, 0.2)
  
  finalPlot
  
  ggarrange(map, legend, 
            labels = c("Map", "Key"),
            ncol = 2, nrow = 1,widths = c(2,1), heights=c(1,1), vjust=-1)
  
  
}

#Main function. Diplay changes in river flood at pixel and catchment level
plotchangemapix=function(basemap,catmap,datar,law="GPD",type,period=c(1950,2019),parlist,hybasf,haz="flood"){
  
  
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
  #st_geometry(data) <- NULL
  
  
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
                     by = list(HYBAS_ID=data$outlets),
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

#Function to compare two scenarios
plotSceComp=function(basemap,data1,data2,date=c(2020),hybasf,valuenames,nco){
  
  
  #data=right_join(outloc,datar,by = c("outlets"="unikout"))
  data1=datariH
  data1$unikout=data1$outlets
  data1$Y2020=data1$Y2019
  
  data2=datariSCF
  data2$unikout=data2$outlets
  data2$Y2020=data2$Y2019
  
  datacomb=inner_join(data1,data2,by="outlets")
  
  
  #mcor=unique(match(data$unikout,parlist$catchment))
  date=1950
  seq=c(1950:2020)
  datacol=names(datacomb)
  aystat=c()
  for (it in 1:length(seq)){
    yr=seq[it]
    dloc1=paste0("Y",yr,".x")
    dloc2=paste0("Y",yr,".y")
    #this need to be changed to make advantage of the whole data
    vc1=match(dloc1,datacol)
    vc2=match(dloc2,datacol)
    tmpval=(datacomb[,vc1])/(datacomb[,vc2])*100-100
    stats=c(yr,mean(tmpval),median(abs(tmpval)),median(tmpval),sd(tmpval))
    aystat=rbind(aystat,stats)
  }
  plot(aystat[,2])
  dateloc1=paste0("Y",date,".x")
  dateloc2=paste0("Y",date,".y")
  datacol=names(datacomb)
  
  vcol1=match(dateloc1,datacol)
  vcol2=match(dateloc2,datacol)
  tmpval=(datacomb[,vcol1])/(datacomb[,vcol2])*100-100
  datacomb$plot=tmpval
  
  #difference in changes
  ploc1x=paste0("Y",period[1],".x")
  ploc1y=paste0("Y",period[1],".y")
  
  ploc2x=paste0("Y",period[2],".x")
  ploc2y=paste0("Y",period[2],".y")
  
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
    geom_sf(data=points,aes(col=plot,geometry=geometry),alpha=1,size=0.25,stroke=0,shape=15)+ 
    geom_sf(fill=NA, color="grey") +
    #geom_sf(regiod,mapping=aes(geometry=geometry,group=name),fill=NA,color="orange")+
    coord_sf(xlim = c(min(nco[,1]),max(nco[,1])), ylim = c(min(nco[,2]),max(nco[,2])))+
    ggtitle(title)+
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
  
  points <- st_as_sf(datacomb[,c(1,2,3,5,157)], coords = c("Var1.x", "Var2.x"), crs = 4326)
  points <- st_transform(points, crs = 3035)
  title=paste0("Contribution of Socioeconomic changes to total changes in 100-y RL flood")
  legend="Contribution (%)"
  palet=c(hcl.colors(11, palette = "RdBu", alpha = NULL, rev = F, fixup = TRUE))
  br=c(-100,-80,-50,-30,-20,-10,-8,-6,-4,-2,0,2,4,6,8,10,20,30,50,80,100)
  limi=c(-10,10)
  ocrap2<-ggplot(basemap) +
    geom_sf(fill="gray85",color="darkgrey",size=0.5)+
    geom_sf(data=points,aes(col=changediff,geometry=geometry),alpha=1,size=0.25,stroke=0,shape=15)+ 
    geom_sf(fill=NA, color="grey") +
    #geom_sf(regiod,mapping=aes(geometry=geometry,group=name),fill=NA,color="orange")+
    coord_sf(xlim = c(min(nco[,1]),max(nco[,1])), ylim = c(min(nco[,2]),max(nco[,2])))+
    ggtitle(title)+
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


#Ridgeplots of changes over a given area
plotHistoDates=function(outloc,datar,law="GPD",type,period=c(1950,2020),parlist,valuenames){
  
  
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
FascPlotChange<- function(datar,outloc,period=c(1952,2020),lims=c(-60, 60),q1,q2){
  
  seqy=c(period[1]:period[2])
  data=datar
  subdata=inner_join(outloc,data,by=c("outlets"="unikout"))
  
  datacol=names(subdata)
  valuenames=paste0("Y",seq(period[1],period[2]))
  valcol=match(valuenames,datacol)
  
  subdata=subdata[,valcol]
  
  datrd=subdata/subdata[,1]*100-100
  #smooth
  datrp=data.frame(unikout=datar$unikout,datrd)
  lwtf=c()
  yvf=c()
  for (y in 1:length(seqy)){
    print(y)
    yx=seqy[y]
    daty=datrp[,which(names(datrp)==paste0("Y",yx))]
    blss=which(daty>1000)
    datrd[blss,which(names(datrd)==paste0("Y",yx))]=NA
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
  
  dfquantiles=aggregate(list(change=dfplot$pixel),
                        by = list(Year=dfplot$year),
                        FUN = function(x) c(mean=mean(x,na.rm=T),q1=quantile(x, q1, na.rm=T),q3=quantile(x, q2, na.rm=T)))
  
  dfquantiles <- do.call(data.frame, dfquantiles)
  names(dfquantiles)=c("Year", "mean","qlow","qhigh")
  
  
  
  dtfinal=c()
  for( y in 1:length(seqy)){
    print(y)
    yx=seqy[y]
    dfp1=dfplot[which(dfplot$year==yx),]
    dfp1=dfp1[which(!is.na(dfp1$pixel)),]
    merde=kde(dfp1$pixel,bgridsize=1000,xmin=-100,xmax=100,h=1)
    
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
    geom_raster(aes(fill=values,interpolate = T))+
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



# Pre-loaded results -----------
#Set data directory
hydroDir<-("D:/tilloal/Documents/LFRuns_utils/data")

#still need to create the outlets file outf
if (!exists("outf")){
  outf=c()
  for( Nsq in 1:88){
    print(Nsq)
    rspace= read.csv(paste0(hydroDir,"/subspace_efas.csv"))
    rspace=rspace[,-1]
    nrspace=rspace[Nsq,]
    #outletname="outletsv8_hybas07_01min"
    #outletname="outlets_hybas09_01min"
    outletname="efas_rnet_100km_01min"
    
    outhybas=outletopen(hydroDir,outletname,nrspace)
    Idstart=as.numeric(Nsq)*10000
    Idstart2=as.numeric(Nsq)*100000
    if (length(outhybas$outlets)>0){
      outhybas$outlets=seq((Idstart+1),(Idstart+length(outhybas$outlets)))
      outhybas$outl2=seq((Idstart2+1),(Idstart2+length(outhybas$outlets)))
      outhybas$latlong=paste(round(outhybas$Var1,4),round(outhybas$Var2,4),sep=" ")
      #outcut=which(!is.na(match(outhybas$outlets,parlist$catchment)))
      # zebi=seq(parlist$catchment[1],parlist$catchment[length(parlist$catchment)])
      # outcut=which(!is.na(match(outhybas$outlets,zebi)))
      outhloc=outhybas
      outf=rbind(outf,outhloc)
    }
  }
}
## Spatial data for catchments ----

### Hybas07 ----
Catchmentrivers7=read.csv(paste0(hydroDir,"/Catchments/from_hybas_eu_onlyid.csv"),encoding = "UTF-8", header = T, stringsAsFactors = F)
hybas07 <- read_sf(dsn = paste0(hydroDir,"/Catchments/hydrosheds/hybas_eu_lev07_v1c.shp"))
hybasf7=fortify(hybas07)
Catamere07=inner_join(hybasf7,Catchmentrivers7,by= "HYBAS_ID")
Catamere07$llcoord=paste(round(Catamere07$POINT_X,4),round(Catamere07$POINT_Y,4),sep=" ")
cst7=right_join(Catamere07,outf,by= c("llcoord"="latlong"))
GNF=cst7
st_geometry(GNF)=NULL
rm(Catamere07)
#cst7=st_transform(cst7,  crs=3035)

outlethybas07="outletsv8_hybas07_01min"

outhybas07=outletopen(hydroDir,outlethybas07)


#matching outlets with pixel Ids

outhybas07$latlong=paste(round(outhybas07$Var1,4),round(outhybas07$Var2,4),sep=" ")

mhy=match(outhybas07$latlong,outf$latlong)
outhybas07$outID=outf$outl2[mhy]


### Plot parameters ----
palet2=c(hcl.colors(9, palette = "Blues", alpha = NULL, rev = TRUE, fixup = TRUE))
outletname="efas_rnet_100km_01min"
outll=outletopen(hydroDir,outletname)
cord.dec=outll[,c(2,3)]
cord.dec = SpatialPoints(cord.dec, proj4string=CRS("+proj=longlat"))
cord.UTM <- spTransform(cord.dec, CRS("+init=epsg:3035"))
nco=cord.UTM@coords
world <- ne_countries(scale = "medium", returnclass = "sf")
Europe <- world[which(world$continent == "Europe"),]
e2=st_transform(Europe,  crs=3035)
w2=st_transform(world,  crs=3035)
tsize=12
osize=12
Impdates=seq(1950,2020,by=10)
valuenames=paste0("Y",Impdates)
catmap=cst7
rm(cst7)
basemap=w2


#load UpArea
#load upstream area
main_path = 'D:/tilloal/Documents/06_Floodrivers/'
valid_path = paste0(main_path,'DataPaper/')
outletname="/GIS/upArea_European_01min.nc"
dir=valid_path
outf$idlalo=paste(outf$idlo, outf$idla, sep=" ")
UpArea=UpAopen(valid_path,outletname,outf)
head(UpArea)



##Loading saved results in .Rdata ---------------------------

#load historical run
haz="Drought"
namefile="Drought.nonfrost.Histo3"

#haz="Flood"
#namefile="flood.Histo4"

load(file=paste0(hydroDir,"/",haz,"/params.",namefile,".Rdata"))
load(file=paste0(hydroDir,"/",haz,"/RL100.",namefile,".Rdata"))
load(file=paste0(hydroDir,"/",haz,"/peaks.",namefile,".Rdata"))
#Peaksave=data.table(Peaksave)
gc()

Paramsfl=data.table(Paramsfl[,-c(4:9,17)])
ParamsflH=Paramsfl
PeakH=Peaksave
RLGPDflH=RLGPDfl
rm(Paramsfl,RLGPDfl)
gc()

#load results from Socio-CF run
namefile="Drought.nonfrost.SocCF3"
#namefile="flood.socCF4"
load(file=paste0(hydroDir,"/",haz,"/params.",namefile,".Rdata"))
load(file=paste0(hydroDir,"/",haz,"/RL100.",namefile,".Rdata"))
load(file=paste0(hydroDir,"/",haz,"/peaks.",namefile,".Rdata"))

RLGPDflSCF=RLGPDfl
Paramsfl=Paramsfl[,-c(4:9,17)]
ParamsflSCF=data.table(Paramsfl)
PeakSCF=Peaksave
rm(Paramsfl,RLGPDfl)
gc()

#load results from Res+WU CF run
namefile="Drought.nonfrost.RWCF3"
#namefile="flood.RWCF4"
load(file=paste0(hydroDir,"/",haz,"/params.",namefile,".Rdata"))
load(file=paste0(hydroDir,"/",haz,"/RL100.",namefile,".Rdata"))
load(file=paste0(hydroDir,"/",haz,"/peaks.",namefile,".Rdata"))
RLGPDflRWCF=RLGPDfl
Paramsfl=Paramsfl[,-c(4:9,17)]
ParamsflRWCF=data.table(Paramsfl)
PeakRWCF=Peaksave
rm(Paramsfl,RLGPDfl)
gc()

#load results from Res CF run
namefile="Drought.nonfrost.ResCF3"
#namefile="flood.RCF4"
load(file=paste0(hydroDir,"/",haz,"/params.",namefile,".Rdata"))
load(file=paste0(hydroDir,"/",haz,"/RL100.",namefile,".Rdata"))
load(file=paste0(hydroDir,"/",haz,"/peaks.",namefile,".Rdata"))

RLGPDflRCF=RLGPDfl
Paramsfl=Paramsfl[,-c(4:9,17)]
ParamsflRCF=data.table(Paramsfl)
PeakRCF=Peaksave
rm(Paramsfl,RLGPDfl)
gc()


###### SINGLE CATCHMENT PLOTS #########


#Load Arno low flow
RPGPDfl=c()
RLGPDfl=c()
naallfl=c()
Paramsfl=c()
parlist=c()
Peaks=c()
Peaksave=c()
lf=list.files(path ="D:/tilloal/Documents/06_Floodrivers/WD_runs_Arno/TSEVA/out", full.names = TRUE, recursive = TRUE)
lf
file=lf[4]
load(file)
fils=sub(".*/", "", file)
tt=unlist(strsplit(fils, "[_]"))
Nsq=53
if (cal==F){
  if (Nsq>88) Nsq=floor(Nsq/10)
}else if (cal==T){
  Nsq=floor(Nsq/10)
}

out.type=sub("\\_.*", "", fils)
print(Nsq)
parlist=Results$parameters
parlist=data.table(parlist)

# unikpar=unique(parlist$catchment)
# dparmerde=diff(unikpar)
# plot(dparmerde)

rspace= read.csv(paste0(hydroDir,"/subspace_efas.csv"))
rspace=rspace[,-1]
nrspace=rspace[Nsq,]

#specify with which outlets I am working with
#outletname="outlets_hybas09_01min"
outletname="efas_rnet_100km_01min"

outhybas=outletopen(hydroDir,outletname,nrspace)
Idstart=as.numeric(Nsq)*10000
Idstart2=as.numeric(Nsq)*100000
outf=c()
if (length(outhybas$outlets)>0){
  outhybas$outlets=seq((Idstart+1),(Idstart+length(outhybas$outlets)))
  outhybas$outl2=seq((Idstart2+1),(Idstart2+length(outhybas$outlets)))
  outhybas$latlong=paste(round(outhybas$Var1,4),round(outhybas$Var2,4),sep=" ")
  #outcut=which(!is.na(match(outhybas$outlets,parlist$catchment)))
  zebi=seq(parlist$catchment[1],parlist$catchment[length(parlist$catchment)])
  outcut=which(!is.na(match(outhybas$outl2,zebi)))
  outhloc=outhybas[outcut,]
  outf=rbind(outf,outhloc)
}
unikout=unique(outhloc$outl2)

pc=c()


RetLevGPD=Results$RetLevGPD
RetPerGPD=Results$RetPerGPD
RetLevGEV=Results$RetLevGEV
RetPerGEV=Results$RetPerGEV

rm50=which(colnames(RetLevGPD)=="1950")

if (length(rm50)>0){
  print("remove")
  RetLevGPD=RetLevGPD[,-rm50]
  RetPerGPD=RetPerGPD[,-rm50]
  RetLevGEV=RetLevGEV[,-rm50]
  RetPerGEV=RetPerGEV[,-rm50]
}



Impdates=names(RetLevGPD)

names(RetLevGPD)=paste0("Y",Impdates)
names(RetPerGPD)=paste0("Y",Impdates)
RetLevGPD=data.table(RetLevGPD,unikout)
RetPerGPD=data.table(RetPerGPD,unikout)

names(RetLevGEV)=paste0("Y",Impdates)
names(RetPerGEV)=paste0("Y",Impdates)
RetLevGEV=data.table(RetLevGEV,unikout)
RetPerGEV=data.table(RetPerGEV,unikout)

RPGPDfl=rbind(RPGPDfl,RetPerGPD)
RLGPDfl=rbind(RLGPDfl,RetLevGPD)
Paramsfl=rbind(Paramsfl,parlist)
Peakfoir=Results$Peaks
myIDs=paste(Peakfoir$timeID,Peakfoir$catch,sep=" ")
unid=unique(myIDs)
locu=match(unid,myIDs)
Peaks=Peakfoir[locu,]
Peaks$time=as.POSIXct(Peaks$time,origin="1970-01-01 00:00:00")
Peaks$d=yday(Peaks$time)
Peaks$year=year(Peaks$time)
Peaksave=rbind(Peaksave,Peaks)


# save(RLGPDfl,file=paste0("D:/tilloal/Documents/06_Floodrivers/WD_runs_Arno/TSEVA/out/RL100.drought.nonfrost.Arno.NGWR.Rdata"))
# save(Paramsfl,file=paste0("D:/tilloal/Documents/06_Floodrivers/WD_runs_Arno/TSEVA/out/params.drought.nonfrost.Arno.NGWR.Rdata"))
# save(Peaksave,file=paste0("D:/tilloal/Documents/06_Floodrivers/WD_runs_Arno/TSEVA/out/Peaks.drought.nonfrost.Arno.NGWR.Rdata"))

RLGPDArno=RLGPDfl
ParamsArno=Paramsfl
PeakArno=Peaksave

load(file=paste0("D:/tilloal/Documents/06_Floodrivers/WD_runs_Arno/TSEVA/out/RL100.drought.nonfrost.Arno.NGWR.Rdata"))
load(file=paste0("D:/tilloal/Documents/06_Floodrivers/WD_runs_Arno/TSEVA/out/params.drought.nonfrost.Arno.NGWR.Rdata"))
load(file=paste0("D:/tilloal/Documents/06_Floodrivers/WD_runs_Arno/TSEVA/out/Peaks.drought.nonfrost.Arno.NGWR.Rdata"))

RLGPDArnold=RLGPDfl
ParamsArnold=Paramsfl
PeakArnold=Peaksave

#create a time vector with years
dates <- seq.Date(ymd("1951-06-01"), ymd("2020-06-01"), by = "year")
dates= as.Date(dates, format = "%Y-%m-%d")

time=seq.Date(ymd("1951-01-01"), ymd("2020-12-31"), by = "day")
start_index=1
indices_to_extract <- seq(from = start_index, to = length(time), by = 7)

timeStamps=time[indices_to_extract]


#I pick a pixel 

plot_riverchange<-function(plot_inputs,main,dates){
  
  yr.deb <-  seq(as.Date("1950-01-15"), by="5 years", length=14)
  Climtrend=plot_inputs$Climtrend
  Soctrend=plot_inputs$Soctrend
  Restrend=plot_inputs$Restrend
  Wutrend=plot_inputs$Wutrend
  RLGPDH=plot_inputs$RLGPDH
  RLGPDS=plot_inputs$RLGPDS
  RLGPDR=plot_inputs$RLGPDR
  RLGPDW=plot_inputs$RLGPDW
  RLGPDNRF=plot_inputs$RLGPDNRF
  
  pid=plot_inputs$pid
  pwid=plot_inputs$pwid
  prid=plot_inputs$prid
  psid=plot_inputs$psid
  pnrid=plot_inputs$pnrfid
  
  ctrend=c(Climtrend)
  lctrend=c(Climtrend+Soctrend)
  wlctrend=c(lctrend+Wutrend)
  rwlctrend=c(wlctrend+Restrend)
  
  lutrend=c(lctrend,rev(ctrend))
  wutrend=c(wlctrend,rev(lctrend))
  retrend=c(rwlctrend,rev(wlctrend))
  pdates=c(dates,rev(dates))
  
  plot(dates,RLGPDH)
  qlim=c(.8*min(c(pid$value,prid$value,pnrid$value, psid$value)), 1.2*max(c(pid$value,prid$value,psid$value)))
  
  plot(pid$time2, pid$value, col=alpha("grey",.8) ,pch=16,axes=FALSE,xaxs="i",yaxs="i",ylim=qlim,
       xlab = NA, ylab="")
  
  points(psid$time2,psid$value,col=alpha("royalblue",.8),type="p",pch=1)
  points(prid$time2,prid$value,col=alpha("darkgreen",.8),type="p",pch=3)
  points(pwid$time2,pwid$value,col=alpha("orange",.8),type="p",pch=4)
  points(pnrid$time2,pnrid$value,col=alpha("darkred",.8),type="p",pch=4)
  
  lines(dates,RLGPDH,col="lightgrey",lwd=2,lty=1)
  lines(dates,RLGPDS,col="royalblue",lwd=2,lty=2)
  lines(dates,RLGPDR,col="darkgreen",lwd=2,lty=3)
  lines(dates,RLGPDW,col="orange",lwd=2,lty=4)
  lines(dates,RLGPDNRF,col="darkred",lwd=2,lty=4)
  mtext(main,3,font = 2,line = 0.5,cex = 1.5)
  abline(v = yr.deb, col="lightgrey", lty=2)
  abline(h = 0, col=1)
  axis(2,cex.axis=1)
  title(ylab = expression(paste("Q (",m^3/s,")")),cex.lab=1.5,line=2, xlab="years")
  axis(1, yr.deb, label=format(yr.deb,"%Y"),cex.axis=1)
  box()
  ## Trace des erreurs absolues
  # polygon(c(dates[1],dates,dates[length(dates)]),c(0,abs(Climtrend),0),
  #         col=alpha("royalblue",.5),border="royalblue")
  # polygon(pdates,abs(lutrend),
  #         col=alpha("orange",.5),border="orange")
  # polygon(pdates,abs(wutrend),
  #         col=alpha("grey",.5),border="darkgreen")
  # polygon(pdates,abs(retrend),
  #         col=alpha("grey",.5),border="darkred")
  ## Trace de la grille mensuelle
  
  ## Definition de la legende
  legend("topleft", leg=c("Historical 10y RL","Res+WU static 10y RL","Res static 10y RL","Socio static 10y RL","Historical - No Return flow from GW"),
         lwd=c(2,2,2,2), col=c("grey","orange","darkgreen","royalblue","darkred"),
         cex=1, lty=c(1,2,3,4,5), bg=alpha("white",.6))
} 


#who is the Arno here
pix=4103620

pix=5300097

#Arno
pix=5302991

#Arno gauge
pix=5303259

haz="Drought"

ylplot=seq(as.POSIXct("1950-06-01"),as.POSIXct("2020-06-01"),by="year")
Peak1=PeakSCF[which(PeakSCF$catch==pix),]

Peak2=PeakRWCF[which(PeakRWCF$catch==pix),]

Peak3=PeakRCF[which(PeakRCF$catch==pix),]

Peak4=PeakH[which(PeakH$catch==pix),]

Peak5=PeakArno[which(PeakArno$catch==pix),]

if (haz=="Drought"){
  
  Peak1$time2=timeStamps[Peak1$timeID]
  Peak2$time2=timeStamps[Peak2$timeID]
  Peak3$time2=timeStamps[Peak3$timeID]
  Peak4$time2=timeStamps[Peak4$timeID]
  Peak5$time2=timeStamps[Peak5$timeID]
  
  Peak1$value=-Peak1$value
  Peak2$value=-Peak2$value
  Peak3$value=-Peak3$value
  Peak4$value=-Peak4$value
  Peak5$value=-Peak5$value
  
  RL1=as.numeric(-RLGPDflSCF[which(RLGPDflSCF$unikout==pix),-71])
  RL2=as.numeric(-RLGPDflRWCF[which(RLGPDflRWCF$unikout==pix),-71])
  RL3=as.numeric(-RLGPDflRCF[which(RLGPDflRCF$unikout==pix),-71])
  RL4=as.numeric(-RLGPDflH[which(RLGPDflH$unikout==pix),-71])
  RL5=as.numeric(-RLGPDArno[which(RLGPDArno$unikout==pix),-71])
}else{
  
  Peak1$time2=as.Date(Peak1$time)
  Peak2$time2=as.Date(Peak2$time)
  Peak3$time2=as.Date(Peak3$time)
  Peak4$time2=as.Date(Peak4$time)
  Peak5$time2=as.Date(Peak5$time)
  
  RL1=as.numeric(RLGPDflSCF[which(RLGPDflSCF$unikout==pix),-71])
  RL2=as.numeric(RLGPDflRWCF[which(RLGPDflRWCF$unikout==pix),-71])
  RL3=as.numeric(RLGPDflRCF[which(RLGPDflRCF$unikout==pix),-71])
  RL4=as.numeric(RLGPDflH[which(RLGPDflH$unikout==pix),-71])
  RL5=as.numeric(RLGPDArno[which(RLGPDArno$unikout==pix),-71])
  
}

plot(Peak1$time2,Peak1$value,pch=16,col="blue")
points(Peak2$time2,Peak2$value,pch=16,col="orange")
points(Peak3$time2,Peak3$value,pch=16,col="darkgreen")
points(Peak5$time2,Peak5$value,pch=16,col="darkred")
points(Peak4$time2,Peak4$value,pch=16, col=1)


RL1[which(RL1<0)]=0
RL2[which(RL2<0)]=0
RL3[which(RL3<0)]=0
RL4[which(RL4<0)]=0
RL5[which(RL5<0)]=0

#change from climate
Climtrend=RL1-RL1[1]
Soctrend=(RL2-RL1)
Wutrend=RL3-RL2
Restrend=RL4-RL3


# Define the size in centimeters
width_cm <- 20
height_cm <- 15

# Define the resolution in DPI (300 for high quality)
dpi <- 300

# Convert centimeters to pixels
width_px <- width_cm * (dpi / 2.54)
height_px <- height_cm * (dpi / 2.54)
# Open a JPEG device
jpeg(paste("D:/tilloal/Documents/LFRuns_utils/TrendAnalysis/plots/",pix, "_drought_scenarii_ArnoVico.jpg", sep = ""),
     width = width_px, height = height_px, quality = 100, res=dpi)

# create a function from this  
main="Reno @ Bologna"
main="Arno @ Firenze"
main="Arno @ Vicopisano"


plot_inputs=list(Climtrend=Climtrend,Soctrend=Soctrend,Restrend=Restrend,Wutrend=Wutrend,
                 RLGPDH=RL4,RLGPDR=RL3,RLGPDS=RL1,RLGPDW=RL2,RLGPDNRF=RL5,
                 pid=Peak4,psid=Peak1,prid=Peak3,pwid=Peak2,pnrfid=Peak5)

test=plot_riverchange(plot_inputs,main,dates)

dev.off()



#play with different scenarii

cor_st=RL1-RL5

cor_hist=RL5 - (RL4-cor_st[1])
cor_res= - (RL3-RL4)
plot(cor_res)
plot(RL4-cor_st[1]+cor_hist)
dc=RL4[1]-RL3[1]

plot(RL5,ylim=c(0,5))
points(RL3-cor_st[1]+cor_hist+cor_res[1],col=2)
points(RL1-cor_st[1],col=3)
points(RL2-cor_st[1],col=4)
points(RL4-cor_st[1]+cor_hist,col=5)


RLt1=RL1-cor_st[1]
RLt2=RL2-cor_st[1]
RLt3=RL3-cor_st[1]+cor_hist+cor_res[1]



# Define the size in centimeters
width_cm <- 20
height_cm <- 15

# Define the resolution in DPI (300 for high quality)
dpi <- 300

# Convert centimeters to pixels
width_px <- width_cm * (dpi / 2.54)
height_px <- height_cm * (dpi / 2.54)
# Open a JPEG device
jpeg(paste("D:/tilloal/Documents/LFRuns_utils/TrendAnalysis/plots/",pix, "_drought_correction_ArnoVico.jpg", sep = ""),
     width = width_px, height = height_px, quality = 100, res=dpi)

# create a function from this  
main="Reno @ Bologna"
main="Arno @ Firenze"
main="Arno @ Vicopisano"

plot(Peak5$time2, Peak5$value, col=alpha("darkred",.8) ,pch=16,axes=FALSE,xaxs="i",yaxs="i",ylim=c(0,8),
     xlab = NA, ylab="")

lines(dates,RL5,col="darkred",lwd=2,lty=1)
lines(dates,RLt1,col="royalblue",lwd=2,lty=2)
lines(dates,RLt3,col="darkgreen",lwd=2,lty=3)
lines(dates,RLt2,col="orange",lwd=2,lty=4)
mtext(main,3,font = 2,line = 0.5,cex = 1.5)
abline(v = yr.deb, col="lightgrey", lty=2)
abline(h = 0, col=1)
axis(2,cex.axis=1)
title(ylab = expression(paste("Q (",m^3/s,")")),cex.lab=1.5,line=2, xlab="years")
axis(1, yr.deb, label=format(yr.deb,"%Y"),cex.axis=1)
box()
## Trace des erreurs absolues
# polygon(c(dates[1],dates,dates[length(dates)]),c(0,abs(Climtrend),0),
#         col=alpha("royalblue",.5),border="royalblue")
# polygon(pdates,abs(lutrend),
#         col=alpha("orange",.5),border="orange")
# polygon(pdates,abs(wutrend),
#         col=alpha("grey",.5),border="darkgreen")
# polygon(pdates,abs(retrend),
#         col=alpha("grey",.5),border="darkred")
## Trace de la grille mensuelle

## Definition de la legende
legend("bottomleft", leg=c("Historical","Res+WU static (cor)","Res static (cor)","Socio static (cor)"),
       lwd=c(2,2,2,2), col=c("darkred","orange","darkgreen","royalblue"),
       cex=1, lty=c(1,2,3,4,5), bg=alpha("white",.6))


dev.off()



#impact of LZ on low flows
pix=5303259
PeakAo=PeakArno[which(PeakArno$catch==pix),]
PeakAod=PeakArnold[which(PeakArnold$catch==pix),]

PeakAo$time2=timeStamps[PeakAo$timeID]
PeakAod$time2=timeStamps[PeakAod$timeID]

PeakAo$value=PeakAo$value
PeakAod$value=PeakAod$value


RLArno=as.numeric(-RLGPDArno[which(RLGPDArno$unikout==pix),-71])
RLArnold=as.numeric(-RLGPDArnold[which(RLGPDArnold$unikout==pix),-71])


plot(PeakAo$time2,PeakAo$value,pch=16,col="blue")
points(PeakAod$time2,PeakAod$value,pch=16,col="orange")


RLArno[which(RLArno<0)]=0
RLArnold[which(RLArnold<0)]=0

yr.deb <-  seq(as.Date("1950-01-15"), by="5 years", length=14)


# Open a JPEG device
jpeg(paste("D:/tilloal/Documents/LFRuns_utils/TrendAnalysis/plots/",pix, "_LZImpact.jpg", sep = ""),
     width = width_px, height = height_px, quality = 100, res=dpi)

# create a function from this  
main="LZ impact Arno @ Vicopisano"

plot(PeakAo$time2, PeakAo$value, col=alpha("darkred",.8) ,pch=16,axes=FALSE,xaxs="i",yaxs="i",ylim=c(0,8),
     xlab = NA, ylab="")
points(PeakAod$time2, PeakAod$value, col=alpha("royalblue",.8) ,pch=10)
lines(dates,RLArno,col="darkred",lwd=2,lty=1)
lines(dates,RLArnold,col="royalblue",lwd=2,lty=2)
abline(v = yr.deb, col="lightgrey", lty=2)
abline(h = 0, col=1)
axis(2,cex.axis=1)
title(ylab = expression(paste("Q (",m^3/s,")")),cex.lab=1.5,line=2, xlab="years")
axis(1, yr.deb, label=format(yr.deb,"%Y"),cex.axis=1)
box()
## Trace des erreurs absolues
# polygon(c(dates[1],dates,dates[length(dates)]),c(0,abs(Climtrend),0),
#         col=alpha("royalblue",.5),border="royalblue")
# polygon(pdates,abs(lutrend),
#         col=alpha("orange",.5),border="orange")
# polygon(pdates,abs(wutrend),
#         col=alpha("grey",.5),border="darkgreen")
# polygon(pdates,abs(retrend),
#         col=alpha("grey",.5),border="darkred")
## Trace de la grille mensuelle

## Definition de la legende
legend("bottomleft", leg=c("LZ_calibrated","LZ_default"),
       lwd=c(2,2), col=c("darkred","royalblue"),
       cex=1, lty=c(1,2), bg=alpha("white",.6))

dev.off()
#compare with observations:

#Load my observations

dis_obs=read.csv(file="D:/tilloal/Documents/06_Floodrivers/WD_runs_Arno/observations/port_TOS01005191.csv",sep=",",header = F)

dis_meta= dis_obs[c(1:16),]
dis_data= dis_obs[-c(1:16),-c(4,5)]
colnames(dis_data)=c("date","Q","data_type")
dis_data=dis_data[-c(1,2,3),]
dis_data$Q=as.numeric(dis_data$Q)
dis_data$time=as.Date(dis_data$date,format = "%d/%m/%Y")
plot(dis_data$time,dis_data$Q)
dis_data$time[which.max(dis_data$Q)]


#extract location of the river gauge
meta_coord=data.frame(lon=dis_meta$V5[6],lat=dis_meta$V3[6])

#Load the netcdf file
#' Extract Discharge Data from NetCDF File at Specific Location
#'
#' @description
#' The `disNcopenloc` function opens a NetCDF file and extracts discharge data
#' for a specified location based on its identifier. It retrieves a time series
#' for the given location and constructs a data frame with discharge values along
#' with corresponding longitude, latitude, and time.
#'
#' @param fname A character string specifying the file name (without the extension)
#'   of the NetCDF file containing the discharge data.
#' @param dir A character string specifying the directory path where the NetCDF
#'   file is located.
#' @param outloc A data frame containing location information, including IDs and
#'   corresponding grid indices for latitude ('idla') and longitude ('idlo').
#' @param idc An integer representing the location identifier index for which the
#'   discharge data will be extracted.
#'
#' @return A data frame with columns for discharge values ('outlets'), location
#'   identifier ('outid'), longitude ('lon'), latitude ('lat'), and time ('time').
#'
#' @examples
#' # Assuming the NetCDF file 'discharge_data.nc' is in the directory 'data/',
#' # 'locations' is a data frame with grid indices, and 'location_id' is the
#' # identifier index:
#' discharge_series <- disNcopenloc("discharge_data", "data/", locations, location_id)
#' # The result is a data frame with the discharge time series for the specified location.
#'
#' @export
#' @importFrom ncdf4 nc_open ncvar_get
#'
#' @seealso
#' `ncdf4::nc_open` and `ncdf4::ncvar_get` for functions used to access and
#' extract data from NetCDF files.
#'
disNcopenloc=function(fname,dir,outloc,idc){
  ncdis=paste0(dir,"/",fname,".nc")
  ncd=nc_open(ncdis)
  name.vb=names(ncd[['var']])
  namev="dis"
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


#match pixel 


### Distance between "official gauges" and EFAS points -----------------------
points=outf[,c(1,2,3)]
Vsfloc=st_as_sf(points, coords = c("Var1", "Var2"), crs = 4326)
Statloc=meta_coord
Stati=st_as_sf(Statloc, coords = c("lon", "lat"), crs = 4326)
dist=c()
for (r in Vsfloc$outlets){
  cat(paste0(r,"\n"))
  v1=Vsfloc[which(Vsfloc$outlets==r),]
  v2=Stati
  oula=st_distance(v1,v2)
  oula=oula/1000
  dist=c(dist,oula)
}

plot(as.numeric(dist[which(dist<10)]))


#match with upstream area
upR=inner_join(outf,UpArea,by="outl2")
upR$distance=dist

#take the one with smallest distance.... validated
matchpix=upR[which.min(upR$distance),]

#Plots of discharge from NGWRF and classic run (timeseries/validation)
fname="dis_53_1950_2020_cf"
dir="D:/tilloal/Documents/LFRuns_utils/data/timeseries/Validations"
  
dis_mod=disNcopenloc(fname,dir,matchpix,1)
txx=as.POSIXct((dis_mod$time*3600*24 -3600), origin="1979-01-01 00:00:00")

#remove 1950
rmv=which(year(txx)==1950)
#remove 1950 which is not reliable
if (length(rmv)>0){
  dis_mod=dis_mod[-rmv,]
}
plot(dis_mod$outlets[1:3000])

data=data.frame(time=txx[-rmv],Q=dis_mod$outlets)
names(data)=c("date","Qs")

data$date2=data$date-(3600)
data$day=as.Date(data$date2)

dth1=data
daily_data<-aggregate(list(Q=data$Qs),
 by = list(days=data$day),
 FUN = function(x) mean = mean(x,na.rm=T))


plot(daily_data)


#load default HERA
valid_path = paste0(main_path,'DataPaper/')
fname="dis_1951_2020_arno"
dir="D:/tilloal/Documents/06_Floodrivers/WD_Runs_Arno/out"

dis_mod2=disNcopenloc(fname,dir,matchpix,1)
txx=as.POSIXct((dis_mod2$time*3600*24 -3600), origin="1979-01-01 00:00:00")

data=data.frame(time=txx,Q=dis_mod2$outlets)
names(data)=c("date","Qs")

data$date2=data$date-(3600)
data$day=as.Date(data$date2)

data[which(yday(data$day)<6),]=dth1[which(yday(data$day)<6),]
data[which(yday(data$day)>=365),]=dth1[which(yday(data$day)>=365),]

daily_data2<-aggregate(list(Q=data$Qs),
                      by = list(days=data$day),
                      FUN = function(x) mean = mean(x,na.rm=T))

#load LZ cal 
valid_path = paste0(main_path,'DataPaper/')
fname="dis_1951_2020_arnoZ"
dir="D:/tilloal/Documents/06_Floodrivers/WD_Runs_Arno/out"

dis_mod3=disNcopenloc(fname,dir,matchpix,1)
txx=as.POSIXct((dis_mod2$time*3600*24 -3600), origin="1979-01-01 00:00:00")

data=data.frame(time=txx,Q=dis_mod3$outlets)
names(data)=c("date","Qs")

data$date2=data$date-(3600)
data$day=as.Date(data$date2)

data[which(yday(data$day)<6),]=dth1[which(yday(data$day)<6),]
data[which(yday(data$day)>=365),]=dth1[which(yday(data$day)>=365),]

daily_data3<-aggregate(list(Q=data$Qs),
                       by = list(days=data$day),
                       FUN = function(x) mean = mean(x,na.rm=T))

#daily_data2[which(yday(daily_data2$days)==1),]=daily_data2[which(yday(daily_data2$days)==1)-1,]+rnorm(1,0,10)
plot(daily_data2)
plot(daily_data3)
plot(dis_mod$outlets[1000:1600])
points(data$Qs[1000:1600],col=2)
library(hydroGOF)
kge_2runs=KGE(daily_data2$Q,daily_data3$Q, na.rm=TRUE, method="2012",out.type="full")



#Now KGE against observations

md=which(dis_data$time>=daily_data$days[1])
dis_datac=dis_data[md,]

dis_data_kge=left_join(daily_data,dis_datac,by=c("days"="time"))

kge_h1=KGE(dis_data_kge$Q.x,dis_data_kge$Q.y, na.rm=TRUE, method="2012",out.type="full")

dis_data_kge2=left_join(daily_data2,dis_datac,by=c("days"="time"))

kge_h2=KGE(dis_data_kge2$Q.x,dis_data_kge2$Q.y, na.rm=TRUE, method="2012",out.type="full")

dis_data_kge3=left_join(daily_data3,dis_datac,by=c("days"="time"))

kge_h3=KGE(dis_data_kge3$Q.x,dis_data_kge3$Q.y, na.rm=TRUE, method="2012",out.type="full")

dis_data_kge4=left_join(daily_data3,daily_data2,by=c("days"="days"))
#Now regime plots

source("~/06_Floodrivers/DataPaper/Code/HERA/functions.R")

#put data into right format


  
plotQj <- function(data,catch,UpAloc,run,labels){
  #~ data <- as.simu(data)
  names(data)=c("date1","Q1","date2","Q2")
  ## Suppression des lignes sans couples Qobs/Qsim
  data <- subset( data, !is.na(Q1) & !is.na(Q2))
  ## Preparation de la sequence de mois
  mois.deb <-  seq(as.Date("1950-01-15"), by="month", length=12)
  ## Vecteur de catÃ©gories
  period1=paste0(format(range(data$date1)[1],"%Y"),"-",format(range(data$date1)[2],"%Y"))
  period2=paste0(format(range(data$date2)[1],"%Y"),"-",format(range(data$date2)[2],"%Y"))
  print(period2)
  n.area=UpAloc
  name.riv=catch
  #main = bquote(.(name.riv)~ " Regime (Area="~.(n.area)~ km^2~") | Periods: "~ .(period1) ~" vs" ~ .(period2))
  main = bquote(.(name.riv)~ "Regime ( UpArea="~.(n.area)~ km^2~") | Period: "~ .(period2))
  ## Suppression des lignes sans couples Qobs/Qsim
  ## Preparation de la sequence de mois
  jours <- as.numeric(format(data$date1,"%j"))
  ## Iddinces des lignes de mÃ©me catÃ©gorie
  ind.j <- tapply(seq(length(jours)), jours, c)
  ind.j <- ind.j[-366]
  
  Qc <- data.frame(date=as.numeric(names(ind.j)),
                   obs=sapply(ind.j, function(x) mean(data$Q1[x], na.rm=TRUE)),
                   sim=sapply(ind.j, function(x) mean(data$Q2[x], na.rm=TRUE)),
                   q25o=sapply(ind.j, function(x) quantile(data$Q1[x],0.25, na.rm=TRUE)),
                   q75o=sapply(ind.j, function(x) quantile(data$Q1[x],.75, na.rm=TRUE)),
                   q25s=sapply(ind.j, function(x) quantile(data$Q2[x],0.25, na.rm=TRUE)),
                   q75s=sapply(ind.j, function(x) quantile(data$Q2[x],.75, na.rm=TRUE)))
  
  md=mean(data$Q1)
  qlim=c(0, max(unlist(Qc[,c(5,7)]),2*md))
  
  
  ## ParamÃ©tres graphiques
  plot(Qc$date, Qc$obs, type="n", axes=FALSE, ylim=qlim,xaxs="i",yaxs="i",
       xlab = NA, ylab="")
  mtext(main,3,font = 2,line = 0.5,cex = 1.2)
  abline(v = format(mois.deb,"%j"), col="lightgrey", lty=2)
  lines(Qc$date, Qc$obs, col="blue",lwd=2)
  lines(Qc$date, Qc$sim, col="red",lwd=2)
  axis(2,cex.axis=1)
  title(ylab = expression(paste("Q (",m^3/s,")")),cex.lab=1.5,line=2)
  axis(1, format(mois.deb,"%j"), label=format(mois.deb,"%b"),cex.axis=1)
  box()
  polygon(c(Qc$date,rev(Qc$date)),c(Qc$q25o,rev(Qc$q75o)),
          col=alpha("lightblue",.3),border="transparent")
  
  polygon(c(Qc$date,rev(Qc$date)),c(Qc$q25s,rev(Qc$q75s)),
          col=alpha("indianred",.3),border="transparent")
  ## Trace des erreurs absolues
  polygon(c(0,Qc$date,365),c(0,abs(Qc$sim-Qc$obs),0),
          col=alpha("lightgrey",.4),border="grey")
  ## Trace de la grille mensuelle
  
  text(x = max(Qc$date)-90, y = max(qlim) * 0.85, labels = labels, pos = 4, cex = .9, col = "black")
  
  
  ## Trace des courbes de dÃ©bits journaliers simules et observes
  
  ## Definition de la legende
  legend("topleft", leg=c(run,"observed","deviation"),
         lwd=c(2,2,1,2), col=c("red","blue","grey"),
         cex=1, lty=c(1,1,2,1), bg=alpha("white",.6))
  #return(Qc)
} 


#Define the size in centimeters
width_cm <- 20
height_cm <- 15

# Define the resolution in DPI (300 for high quality)
dpi <- 300

# Convert centimeters to pixels
width_px <- width_cm * (dpi / 2.54)
height_px <- height_cm * (dpi / 2.54)


dataR1=data.frame(date1=dis_data_kge$days,Q1=dis_data_kge$Q.y,date2=dis_data_kge$days,Q2=dis_data_kge$Q.x)
dataR2=data.frame(date1=dis_data_kge2$days,Q1=dis_data_kge2$Q.y,date2=dis_data_kge2$days,Q2=dis_data_kge2$Q.x)
dataR3=data.frame(date1=dis_data_kge3$days,Q1=dis_data_kge3$Q.y,date2=dis_data_kge3$days,Q2=dis_data_kge3$Q.x)

dataRX=data.frame(date1=dis_data_kge4$days,Q1=dis_data_kge4$Q.y,date2=dis_data_kge4$days,Q2=dis_data_kge4$Q.x)


statX=c(round(kge_2runs$KGE.value,2),round(kge_2runs$KGE.elements,2))

labels <- paste(
  "KGE':", statX[1], "\n",
  "r:", statX[2], "\n",
  "Beta:", statX[3], "\n",
  "Gamma:", statX[4]
)
plotQj(dataRX,"Arno @ Vicopisano",8900,"HERA_ori",labels)


stat1=c(round(kge_h3$KGE.value,2),round(kge_h3$KGE.elements,2))

labels <- paste(
  "KGE':", stat1[1], "\n",
  "r:", stat1[2], "\n",
  "Beta:", stat1[3], "\n",
  "Gamma:", stat1[4]
)

catch="Arno @ Vicopisano"
jpeg(paste("D:/tilloal/Documents/LFRuns_utils/TrendAnalysis/plots/Arno_vico_Hera_LZcal.jpg", sep = ""),
     width = width_px, height = height_px, quality = 100, res=dpi)

plotQj(dataR3,"Arno @ Vicopisano",8900,"HERA_ori",labels)
dev.off()


stat2=c(round(kge_h2$KGE.value,2),round(kge_h2$KGE.elements,2))
labels <- paste(
  "KGE':", stat2[1], "\n",
  "r:", stat2[2], "\n",
  "Beta:", stat2[3], "\n",
  "Gamma:", stat2[4]
)

jpeg(paste("D:/tilloal/Documents/LFRuns_utils/TrendAnalysis/plots/Arno_vico_Hera_modified.jpg", sep = ""),
     width = width_px, height = height_px, quality = 100, res=dpi)

plotQj(dataR2,"Arno @ Vicopisano",8900,"HERA_mod",labels)
dev.off()

#Specialized performance indicators for low flow

dis_data_kge=left_join(daily_data,dis_datac,by=c("days"="time"))

kge1sq_h1=KGE(1/dis_data_kge$Q.x,1/dis_data_kge$Q.y, na.rm=TRUE, method="2012",out.type="full")

dis_data_kge2=left_join(daily_data2,dis_datac,by=c("days"="time"))

kge1sq_h2=KGE(1/dis_data_kge2$Q.x,1/dis_data_kge2$Q.y, na.rm=TRUE, method="2012",out.type="full")




#################################################################################################
#TSEVA on observed data, for comparisaon with modelled values
#################################################################################################

#Experimental TSEVA function

library(RtsEva)

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
    print(TrendTh)
    message(paste0('\nevaluating long term variations of the peaks'))
    if (is.na(TrendTh)){
      TrendTh=try(tsEvaFindTrendThreshold(series, timeStamps, timeWindow),T)
      if(length(TrendTh)==0){
        TrendTh=0.1
      }
      trasfData = tsEvaTransformSeriesToStationaryPeakTrend( timeStamps, series, timeWindow, TrendTh);
    }else{
      if (TrendTh=="MMX"){
        trasfData = tsEvaTransformSeriesToStationaryMMXTrend( timeStamps, series, timeWindow);
        print("using MMX trend")
      }else{
        trasfData = tsEvaTransformSeriesToStationaryPeakTrend( timeStamps, series, timeWindow, TrendTh);
      }
    }
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
        if(numperyear[ipp]>=3*desiredEventsPerYear & ipp<(length(pcts)-5)) skip = floor(length(pcts)/12)
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

#Now fitting TSEVA on observations between 1950 and 2020


workDir = "D:/tilloal/Documents/LFRuns_utils/data/"
haz="drought"
tail="low"
ThDir<-paste0(workDir,"TrendAnalysis/trenTH")
THX=c()
Nsq=53
if (file.exists(paste0(ThDir,"/trenTH_Histo_",tail,"_",Nsq,".csv"))){
  print(Nsq)
  TH1=read.csv(paste0(ThDir,"/trenTH_Histo_",tail,"_",Nsq,".csv"))
  TH2=read.csv(paste0(ThDir,"/trenTH_SCF_",tail,"_",Nsq,".csv"))
  TH3=inner_join(TH1,TH2,by="cid")
  TH3$cid=TH3$cid-Nsq*10000
  TH3$cid=TH3$cid+Nsq*100000
  TH3=inner_join(TH3,matchpix,by=c("cid"="outl2"))
  THX=rbind(THX,TH3)
}

#retain thresholds fro historical run unless it is NA
thresh_vec=data.frame(THX$cid, THX$Th_new.x)
#print(thresh_vec[c(1:10),])
if(length(which(is.na(thresh_vec$TH3.Th_new.x)))>0){
  print("corr")
  thresh_vec$TH3.Th_new.x[which(is.na(thresh_vec$TH3.Th_new.x))]=THX$Th_new.y[which(is.na(thresh_vec$TH3.Th_new.x))]
}
names(thresh_vec)=c("cid","th")
thresh_vec$cid=as.numeric(thresh_vec$cid)


#loading the frost file for drought
if (haz=="drought"){
  load(file=paste0(workDir,"Drought/catchment_frost.Rdata"))
  rmv=which(year(frostcat$time)==1950)
  frostcat=frostcat[-rmv,]
  #remove first day
  frostcat=frostcat[-1,]
  Catchmentrivers7=read.csv(paste0(workDir,"Catchments/from_hybas_eu_onlyid.csv"),encoding = "UTF-8", header = T, stringsAsFactors = F)
  outletname="outletsv8_hybas07_01min"
  outhyb07=outletopen(workDir,outletname,nrspace)
  catmatch=match(outhyb07$outlets,Catchmentrivers7$pointid)
  mycat=Catchmentrivers7[catmatch,]
  
  hybas07 <- read_sf(dsn = paste0(workDir,"Catchments/hydrosheds/hybas_eu_lev07_v1c.shp"))
  hybasf7=fortify(hybas07) 
  Catamere07=inner_join(hybasf7,Catchmentrivers7,by= "HYBAS_ID")
  Catamere07$llcoord=paste(round(Catamere07$POINT_X,4),round(Catamere07$POINT_Y,4),sep=" ") 
  Catf7=inner_join(Catamere07,outhybas,by= c("llcoord"="latlong"))
  st_geometry(Catf7)=NULL	
  tail="low"
  
}

parlist=c()

uparea=matchpix$upa
timeAndSeries=dis_data_kge[,c(1,4)]

txx=timeAndSeries$days
start_time <- Sys.time()
catch=matchpix$outl2
catch2=matchpix$outlets.x
timeStamps=txx
thresh=thresh_vec[which(thresh_vec$cid==catch),]
thresh=thresh$th
frosttime=NA
series=timeAndSeries

names(series)=c("date","Qs")
rmv=which(as.integer(format(series$date, "%Y"))==1950)
if (length(rmv)>0){
  series=series[-rmv,]
}
#remove double days
doubday=which(diff(series$date)==0)+1

if (length(doubday)>0){
  series=series[-doubday,]
}
dt1=min(diff(series$date),na.rm=T)
if (haz=="drought"){
  trans="rev"
  #seasonal split
  catmat=Catf7[which(Catf7$outlets==catch2),]
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
  #series = df.dis[which(df.dis$outlets==catch),c(6,1)] 
  # df.disX=disNcopenloc(filename,hydroDir,outhybas,idfix)
  # series=data.frame(txx,df.disX$outlets)
  # names(series)=c("date","Qs")
  ciPercentile=95
  minPeakDistanceInDays=7
  interflag=0
  # series$date=as.POSIXlt(series$date)
  #series <- max_daily_value(series)
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

if (length(which(is.na(series$dis)))>0){
  print("Na alert")
  seriefill=tsEvaFillSeries(series$timestamp,series$dis)
  series$dis=seriefill
}
timeAndSeries=series
names(timeAndSeries)=c("timestamp","data")
#timeAndSeries$data=jitter(timeAndSeries$data)

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
#I remove the first year (1950)
# rmv=which(as.integer(format(timeAndSeries$timestamp, "%Y"))==1950)
# if (length(rmv)>0){
# timeAndSeries=timeAndSeries[-rmv,]
# }
series=timeAndSeries[,2]
timeWindow = 365.25*30; #time windows in days, the correction is done within the functions
windowSize=366

timeStamps=timeAndSeries$timestamp
#cat(paste0("\nsquare: ", Nsq, " pixel: ",idfix,"/",endid))
Nonstat<-TsEvaNs(timeAndSeries, timeWindow, transfType='trendPeaks',ciPercentile = ciPercentile, minPeakDistanceInDays = minPeakDistanceInDays,lowdt=7,trans=trans,tail = tail,TrendTh = thresh)
nonStationaryEvaParams=Nonstat[[1]]
stationaryTransformData=Nonstat[[2]]

stationaryTransformData$timeStampsDay=unique(as.Date(stationaryTransformData$timeStamps))
pikos=data.frame(nonStationaryEvaParams$potObj$parameters$peaks,nonStationaryEvaParams$potObj$parameters$peakID,nonStationaryEvaParams$potObj$parameters$peakST, nonStationaryEvaParams$potObj$parameters$peakEN)
names(pikos)=c("value","timeID","tIDstart","tIDend")
pikos$time=timeStamps[pikos$timeID]
pikos$catch=rep(catch,length(pikos[,1]))
#Here I need to convert the timeStamp to a daily one if dt is not 1

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
  names(RLgpd)=year(Impdates)[-c(1,2)]
  
  RLgev=nonStationaryEvaParams$gevObj$parameters$annualMax
  names(RLgev)=year(Impdates)[-c(1,2)]
  
  nRPgev=rep(NA, length(Impdates))
  names(nRPgev)=year(Impdates)[-c(1,2)]
  
  nRPgpd=rep(NA, length(Impdates))
  names(nRPgpd)=year(Impdates)[-c(1,2)]
  
  params=data.frame(matrix(ncol=17,nrow=(length(Impdates)-3)))
  params[,1]=rep(catch,7)
  params[,2]=year(Impdates)[-c(1:2)]
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
  for (t in 4:length(Impdates)){
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

#Saving main outputs
catlist=c(catlist,catch)
IRES=c(IRES,interflag)
#Reservoir_i=c(Reservoir_i,Reservoir_alteration)

RetLevGEV=rbind(RetLevGEV,RLgev)

RetLevGPD=rbind(RetLevGPD,RLgpd)

RetPerGEV=rbind(RetPerGEV,nRPgev)

RetPerGPD=rbind(RetPerGPD,nRPgpd)
end_time <- Sys.time()
cat(paste0("\nloop duration: ",round(end_time-start_time,2)," seconds\n"))


catch="Arno @ Vicopisano"
jpeg(paste("D:/tilloal/Documents/LFRuns_utils/TrendAnalysis/plots/Arno_vico_Hera_obs.jpg", sep = ""),
     width = width_px, height = height_px, quality = 100, res=dpi)



yr.deb <-  seq(as.Date("1950-01-15"), by="5 years", length=14)
RL6=as.numeric(-RLgpd)[-c(71,72,73)]

pikos$time2=timeStamps[pikos$timeID]

plot(Peak5$time2, Peak5$value, col=alpha("darkred",.8) ,pch=16,axes=FALSE,xaxs="i",yaxs="i",ylim=c(0,20),
     xlab = NA, ylab="")
points(Peak1$time2, Peak1$value, col=alpha("black",.8) ,pch=4)
points(pikos$time2, -pikos$value, col=alpha("violet",.8) ,pch=3)
lines(dates,RL5,col="darkred",lwd=2,lty=1)
lines(dates,RL6,col="violet",lwd=2,lty=1)
lines(dates,RL1,col="royalblue",lwd=2,lty=2)
lines(dates,RL3,col="darkgreen",lwd=2,lty=3)
lines(dates,RL2,col="orange",lwd=2,lty=4)
mtext(main,3,font = 2,line = 0.5,cex = 1.5)
abline(v = yr.deb, col="lightgrey", lty=2)
abline(h = 0, col=1)
axis(2,cex.axis=1)
title(ylab = expression(paste("Q (",m^3/s,")")),cex.lab=1.5,line=2, xlab="years")
axis(1, yr.deb, label=format(yr.deb,"%Y"),cex.axis=1)
box()
## Trace des erreurs absolues
# polygon(c(dates[1],dates,dates[length(dates)]),c(0,abs(Climtrend),0),
#         col=alpha("royalblue",.5),border="royalblue")
# polygon(pdates,abs(lutrend),
#         col=alpha("orange",.5),border="orange")
# polygon(pdates,abs(wutrend),
#         col=alpha("grey",.5),border="darkgreen")
# polygon(pdates,abs(retrend),
#         col=alpha("grey",.5),border="darkred")
## Trace de la grille mensuelle

## Definition de la legende
legend("topleft", leg=c("Historical noRF","Res+WU static","Res static","Socio static","Observed"),
       lwd=c(2,2,2,2), col=c("darkred","orange","darkgreen","royalblue","violet"),
       cex=1, lty=c(1,2,3,4,1), bg=alpha("white",.6))
dev.off()

