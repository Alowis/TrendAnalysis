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
                     by = list(Id=data$outlets),
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
#Set data directory
hydroDir<-("D:/tilloal/Documents/LFRuns_utils/data")


# Define the get_density function
get_density <- function(x, y, n = 200) {
  dens <- MASS::kde2d(x = x, y = y, n = n)
  ix <- findInterval(x, dens$x)
  iy <- findInterval(y, dens$y)
  ii <- cbind(ix, iy)
  return(dens$z[ii])
}

# Pre-loaded results -----------


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

#load results from previous script 

#floods
load(file=paste0(hydroDir,"/TSEVA/output_plots/outputs_flood_year_relxHR_new.Rdata"))
#droughts
load(file=paste0(hydroDir,"/TSEVA/output_plots/outputs_drought_nonfrost_relxHR_new2.Rdata"))
#I extract the trend at MUTS3 level fist

FloodTrends=Output_fl_year$TrendRegio

#Output_dr_nonfrost=Output_dr_year
DroughtTrends=Output_dr_nonfrost$TrendRegio

driver=unique(FloodTrends$driver)
driver=unique(DroughtTrends$driver)
#I need to aggregate all the changes
# TotalFloodTrend=FloodTrends[which(FloodTrends$driver==driver[1]),c(2:71)]+
#   FloodTrends[which(FloodTrends$driver==driver[2]),c(2:71)]+
#   FloodTrends[which(FloodTrends$driver==driver[3]),c(2:71)]+
#   FloodTrends[which(FloodTrends$driver==driver[4]),c(2:71)]
#   
# TotalDroughtTrend=DroughtTrends[which(DroughtTrends$driver==driver[1]),c(2:71)]+
#   DroughtTrends[which(DroughtTrends$driver==driver[2]),c(2:71)]+
#   DroughtTrends[which(DroughtTrends$driver==driver[3]),c(2:71)]+
#   DroughtTrends[which(DroughtTrends$driver==driver[4]),c(2:71)]

#I only look at the climate trend


### Hybas07 ----
Catchmentrivers7=read.csv(paste0(hydroDir,"/Catchments/from_hybas_eu_onlyid.csv"),encoding = "UTF-8", header = T, stringsAsFactors = F)
hybas07 <- read_sf(dsn = paste0(hydroDir,"/Catchments/hydrosheds/hybas_eu_lev07_v1c.shp"))
hybasf7=fortify(hybas07)
Catamere07=inner_join(hybasf7,Catchmentrivers7,by= "HYBAS_ID")
Catamere07$llcoord=paste(round(Catamere07$POINT_X,4),round(Catamere07$POINT_Y,4),sep=" ")
cst7=right_join(Catamere07,outf,by= c("llcoord"="latlong"))
GNF=cst7
st_geometry(GNF)=NULL

UnHY=unique(GNF$HYBAS_ID)


### HydroRegions ----

GridHR=raster( paste0(hydroDir,"/HydroRegions_raster_WGS84.tif"))
GHR=as.data.frame(GridHR,xy=T)
GHR=GHR[which(!is.na(GHR[,3])),]
a1=(unique(GHR$HydroRegions_raster_WGS84))
GHR$llcoord=paste(round(GHR$x,4),round(GHR$y,4),sep=" ")
GHR_riv=inner_join(GHR,outf,by= c("llcoord"="latlong"))
a2=(unique(GHR_riv$HydroRegions_raster_WGS84))

ups=which(is.na(match(a1,a2)))
mierda=a1[ups]
GHshpp <- read_sf(dsn ="Z:/ClimateRun4/nahaUsers/tilloal/HydroRegions/her_all_adjusted_cor.shp")
HydroRsf=fortify(GHshpp)
length(unique(GHshpp$IRST_NAMEB))

#correction of the Baltic shield
# GHshppCor=GHshpp[which(GHshpp$IRST_NAMEB=="BALTIC SHIELD"),]
# GHshppCor$IRST_NAMEB=c("FINLAND BALTIC COAST SOUTH","FINLAND BALTIC COAST SOUTH","FINLAND NORTH","FINLAND BALTIC COAST NORTH","FINLAND BALTIC COAST CENTRAL","FINLAND CENTRAL")
# GHshpp$IRST_NAMEB[which(GHshpp$IRST_NAMEB=="BALTIC SHIELD")]=c("FINLAND BALTIC COAST SOUTH","FINLAND BALTIC COAST SOUTH","FINLAND NORTH","FINLAND BALTIC COAST NORTH","FINLAND BALTIC COAST CENTRAL","FINLAND CENTRAL")
# st_write(GHshpp,"Z:/ClimateRun4/nahaUsers/tilloal/HydroRegions/her_all_adjusted_cor.shp")

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
basemap=w2

# Load required libraries
library(ggplot2)


#1. TOTAL Trend ------------

## catchment level ---------------
TotalFloodTrend=FloodTrends[which(FloodTrends$driver==driver[5]),c(2:71)]
TotalDroughtTrend=DroughtTrends[which(DroughtTrends$driver==driver[5]),c(2:71)]
TotalFloodTrend$HydroR=FloodTrends[which(FloodTrends$driver==driver[1]),1]
TotalDroughtTrend$HydroR=DroughtTrends[which(DroughtTrends$driver==driver[1]),1]

FlPoint=inner_join(GHshpp,TotalFloodTrend,by=c("CODEB"="HydroR"))
st_geometry(FlPoint)<-NULL
Flplot=inner_join(HydroRsf,FlPoint,by= c("Id"))
period=c(1951,2020)
haz="flood"

Flplot <- st_transform(Flplot, crs = 3035)
DrPoint=inner_join(GHshpp,TotalDroughtTrend,by=c("CODEB"="HydroR"))

st_geometry(DrPoint)<-NULL
Drplot=inner_join(HydroRsf,DrPoint,by= c("Id"))
period=c(1951,2020)
haz="flood"
Drplot <- st_transform(Drplot, crs = 3035)

Flplot$d2020=Drplot$Rchange.Y2020
Flplot$f2020=Flplot$Rchange.Y2020
FlplotTot=Flplot

#Univariate plot for verification
colNA="transparent"
# legend="Relative change (%)"
tsize=16
osize=12
legend="Change in Qsp \n(l/s/km2)"
palet=c(hcl.colors(11, palette = "RdYlBu", alpha = NULL, rev = F, fixup = TRUE))
# 
# br=seq(-100,100,.1)
# labels=br
# limi=c(-50,50)
# ggplot(basemap) +
#   geom_sf(fill="white")+
#   geom_sf(data = Flplot, mapping = aes(fill = f2020), alpha=0.9, color = "transparent", size = 0.01, show.legend = TRUE) +
#   geom_sf(fill=NA, color="gray42") +
#   coord_sf(xlim = c(min(nco[,1]),max(nco[,1])), ylim = c(min(nco[,2]),max(nco[,2])))+
#   scale_fill_gradientn(
#     colors=palet,
#     breaks=br,limits=limi,
#     oob = scales::squish,na.value=colNA, name=legend)  +
#   labs()+
#   theme(axis.title=element_text(size=tsize),
#         panel.background = element_rect(fill = "aliceblue", colour = "grey1"),
#         panel.border = element_rect(linetype = "solid", fill = NA, colour="black"),
#         legend.title = element_text(size=tsize),
#         legend.text = element_text(size=osize),
#         legend.position = "bottom",
#         panel.grid.major = element_line(colour = "grey70"),
#         panel.grid.minor = element_line(colour = "grey90"),
#         legend.key = element_rect(colour = "transparent"),
#         legend.key.size = unit(.8, "cm"))
# 

#create a new variable which make flood and drought on the same scale
#normalize flood change
a=sd(FlplotTot$f2020,na.rm=T)
b=sd(FlplotTot$d2020,na.rm=T)
a=b=1
databitot=FlplotTot
databitot$x=databitot$f2020/a
databitot$y=databitot$d2020/b

breaker1=0
breaker2=5

### flood -----------------
alterclass=data.frame(databitot$x)
alterclass$class=NA
alterclass$class[which(alterclass[,1]<=(-breaker2))]=1
alterclass$class[which(alterclass[,1]>(-breaker2) & alterclass[,1]<(breaker1))]=2
alterclass$class[which(alterclass[,1]>=(breaker1) & alterclass[,1]<(breaker2))]=3
alterclass$class[which(alterclass[,1]>=breaker2)]=4
c1=alterclass$class

### drought ----------------
alterclass=data.frame(databitot$y)
alterclass$class=NA
alterclass$class[which(alterclass[,1]<=(-breaker2))]=1
alterclass$class[which(alterclass[,1]>(-breaker2) & alterclass[,1]<(breaker1))]=2
alterclass$class[which(alterclass[,1]>=(breaker1) & alterclass[,1]<(breaker2))]=3
alterclass$class[which(alterclass[,1]>=breaker2)]=4
c2=alterclass$class

cx=paste(c2,c1,sep="-")
databitot$bi_class=cx
databitot$combined_category <-databitot$bi_class

## Pixel level ------------------------

FloodTrendsP=Output_fl_year$TrendPix
DroughtTrendsP=Output_dr_nonfrost$TrendPix

TotalFloodTrendPix=FloodTrendsP[which(FloodTrendsP$driver==driver[5]),]
TotalDroughtTrendPix=DroughtTrendsP[which(DroughtTrendsP$driver==driver[5]),]

ClimFloodTrendPix=FloodTrendsP[which(FloodTrendsP$driver==driver[1]),]
ClimDroughtTrendPix=DroughtTrendsP[which(DroughtTrendsP$driver==driver[1]),]

LuFloodTrendPix=FloodTrendsP[which(FloodTrendsP$driver==driver[2]),]
LuDroughtTrendPix=DroughtTrendsP[which(DroughtTrendsP$driver==driver[2]),]

ResFloodTrendPix=FloodTrendsP[which(FloodTrendsP$driver==driver[3]),]
ResDroughtTrendPix=DroughtTrendsP[which(DroughtTrendsP$driver==driver[3]),]

WuFloodTrendPix=FloodTrendsP[which(FloodTrendsP$driver==driver[4]),]
WuDroughtTrendPix=DroughtTrendsP[which(DroughtTrendsP$driver==driver[4]),]

# Remove NA rows 
TotalFloodTrendPix = TotalFloodTrendPix[which(!is.na(TotalFloodTrendPix$Var1)),]
TotalDroughtTrendPix = TotalDroughtTrendPix[which(!is.na(TotalDroughtTrendPix$Var1)),]

ClimFloodTrendPix = ClimFloodTrendPix[which(!is.na(ClimFloodTrendPix$Var1)),]
ClimDroughtTrendPix = ClimDroughtTrendPix[which(!is.na(ClimDroughtTrendPix$Var1)),]

LuFloodTrendPix = LuFloodTrendPix[which(!is.na(LuFloodTrendPix$Var1)),]
LuDroughtTrendPix = LuDroughtTrendPix[which(!is.na(LuDroughtTrendPix$Var1)),]

ResFloodTrendPix = ResFloodTrendPix[which(!is.na(ResFloodTrendPix$Var1)),]
ResDroughtTrendPix = ResDroughtTrendPix[which(!is.na(ResDroughtTrendPix$Var1)),]

WuFloodTrendPix = WuFloodTrendPix[which(!is.na(WuFloodTrendPix$Var1)),]
WuDroughtTrendPix = WuDroughtTrendPix[which(!is.na(WuDroughtTrendPix$Var1)),]


points <- st_as_sf(TotalFloodTrendPix, coords = c("Var1", "Var2"), crs = 4326)
points <- st_transform(points, crs = 3035)
Flpixplot=points
Flpixplot$d2020=TotalDroughtTrendPix$Y2020

#normalize changes
a=sd(Flpixplot$Y2020,na.rm=T)
b=sd(Flpixplot$d2020,na.rm=T)
a=b=1
databipi=Flpixplot
databipi$x=databipi$Y2020/a
databipi$y=databipi$d2020/b


breaker1=0
breaker2=5

### floods --------------------
alterclass=data.frame(databipi$x)
alterclass$class=NA
alterclass$class[which(alterclass[,1]<=(-breaker2))]=1
alterclass$class[which(alterclass[,1]>(-breaker2) & alterclass[,1]<(breaker1))]=2
alterclass$class[which(alterclass[,1]>=(breaker1) & alterclass[,1]<(breaker2))]=3
alterclass$class[which(alterclass[,1]>=breaker2)]=4
c1=alterclass$class

### droughts --------------------
alterclass=data.frame(databipi$y)
alterclass$class=NA
alterclass$class[which(alterclass[,1]<=(-breaker2))]=1
alterclass$class[which(alterclass[,1]>(-breaker2) & alterclass[,1]<(breaker1))]=2
alterclass$class[which(alterclass[,1]>=(breaker1) & alterclass[,1]<(breaker2))]=3
alterclass$class[which(alterclass[,1]>=breaker2)]=4
c2=alterclass$class

cx=paste(c2,c1,sep="-")
databipi$bi_class=cx

## plot section --------------

colors <- c(
  "1-1" = "#dd6a29",  # high x, low y
  "2-1" = "#d9926a",  # medium-high x, low y
  "3-1" = "#d6b3a0",  # medium-low x, low y
  "4-1" = "#d3d3d3",  # low x, low y

  "1-2" = "#a36229",  # high x, medium-low y
  "2-2" = "#a08769",  # medium-high x, medium-low y
  "3-2" = "#9ea69f",  # medium-low x, medium-low y
  "4-2" = "#9cc4d2",  # low x, medium-low y
  
  "1-3" = "#635929",  # high x, medium-high y
  "2-3" = "#617b69",  # medium-high x, medium-high y
  "3-3" = "#60979f",  # medium-low x, medium-high y
  "4-3" = "#5fb2d1",  # low x, medium-high y
  
  "1-4" = "#174f28",  # high x, high y
  "2-4" = "#166d68",  # medium-high x, high y
  "3-4" = "#16869e",  # medium-low x, high y
  "4-4" = "#169dd0"  # low x, high y

)


colorp <- c(
  "1-1" = "#dd6a40",  # high x, low y
  "2-1" = "#d9926a",  # medium-high x, low y
  "3-1" = "#DEB887",  # medium-low x, low y
  "4-1" = "#FFD39B",  # low x, low y
  
  "1-2" = "#a36229",  # high x, medium-low y
  "2-2" = "#999999",  # medium-high x, medium-low y
  "3-2" = "#999999",  # medium-low x, medium-low y
  "4-2" = "#9cc4d2",  # low x, medium-low y
  
  "1-3" = "#635929",  # high x, medium-high y
  "2-3" = "#999999",  # medium-high x, medium-high y
  "3-3" = "#999999",  # medium-low x, medium-high y
  "4-3" = "#5fb2d1",  # low x, medium-high y
  
  "1-4" = "#174f28",  # high x, high y
  "2-4" = "#166d68",  # medium-high x, high y
  "3-4" = "#16869e",  # medium-low x, high y
  "4-4" = "#169dd0"  # low x, high y
  
)

colors2 <- c(
  "1-1" = "#dd6a29",  # high x, low y
  "2-1" = "#dd6a29",  # medium-high x, low y
  "3-1" = "#ffe1ff",  # medium-low x, low y
  "4-1" = "#ffe1ff",  # low x, low y
  
  "1-2" = "#dd6a29",  # high x, medium-low y
  "2-2" = "#9ea69f",  # medium-high x, medium-low y
  "3-2" = "#9ea69f",  # medium-low x, medium-low y
  "4-2" = "#ffe1ff",  # low x, medium-low y
  
  "1-3" = "#174f28",  # high x, medium-high y
  "2-3" = "#9ea69f",  # medium-high x, medium-high y
  "3-3" = "#9ea69f",  # medium-low x, medium-high y
  "4-3" = "#169dd0",  # low x, medium-high y
  
  "1-4" = "#174f28",  # high x, high y
  "2-4" = "#174f28",  # medium-high x, high y
  "3-4" = "#169dd0",  # medium-low x, high y
  "4-4" = "#169dd0"  # low x, high y
  
)

loscolors=c("Accelerating" = "#174f28","Drying" = "#dd6a29","Stable"="gray60","Wetting" = "#169dd0","Decelerating" = "burlywood")

bi_pal(pal = colors, dim = 4)
#colors=(c("#dd6a29","#169dd0","#7ebbd2","#d3d3d3","#167984","#174f28","#845e29","#819185","#d8a386"))


ggplot(databitot, aes(x = d2020, y = f2020, fill = combined_category)) +
  geom_point(shape = 21, color = "black", size = 3,show.legend = FALSE) +
  #bi_scale_fill(pal = "BlueOr", dim = 3, na.value=colNA )+
  scale_fill_manual(values = colorp) +
  coord_cartesian(xlim=c(-50,50),ylim=c(-50,50)) +

  # Add main quadrant labels
  annotate("text", x = 18, y = 40, label = "Wetting", color = "gray12", size = 5, fontface = "bold") +
  annotate("text", x = -18, y = 40, label = "Accelerating", color = "gray12", size = 5, fontface = "bold") +
  annotate("text", x = -18, y = -40, label = "Drying", color = "gray12", size = 5, fontface = "bold") +
  annotate("text", x = 18, y = -40, label = "Decelerating", color = "gray12", size = 5, fontface = "bold") +

  # Set axis labels
  labs(x = "Change in drought flows (l/s/km2)",
       y = "Change in flood flows (l/s/km2)") +

  # Customize the theme
  theme(axis.title=element_text(size=tsize),
        panel.background = element_rect(fill = "transparent", colour = "grey1"),
        panel.border = element_rect(linetype = "solid", fill = NA, colour="black"),
        legend.title = element_text(size=tsize),
        legend.text = element_text(size=osize),
        legend.position = "bottom",
        panel.grid.major = element_line(colour = "grey70"),
        panel.grid.minor = element_line(colour = "grey90"),
        legend.key = element_rect(fill = "transparent", colour = "transparent"),
        legend.key.size = unit(.8, "cm"))


map <- ggplot(basemap) +
  geom_sf(fill="white")+
  geom_sf(data = databitot, mapping = aes(fill = combined_category), alpha=0.7, color = "transparent", size = 0.01,show.legend = F) +
  geom_sf(data = databipi, mapping = aes(col = bi_class,geometry=geometry,size=upa), alpha=1,stroke=0,shape=15, show.legend = FALSE) +
  geom_sf(fill=NA, color="gray42") +
  # bi_scale_fill(pal = "BlueOr", dim = 3, na.value=colNA ) +
  scale_fill_manual(values = colorp) +
  coord_sf(xlim = c(min(nco[,1]),max(nco[,1])), ylim = c(min(nco[,2]),max(nco[,2])))+
  scale_size(range = c(0.08, 0.4), trans="sqrt",name= expression(paste("Upstream area ", (km^2),
                                                                       sep = " ")),
             breaks=c(101,1000,10000,100000,500000), labels=c("100","1000", "10 000", "100 000", "500 000"),
             guide = "none")+
  #bi_scale_color(pal = "BlueOr", dim = 3, na.value=colNA ) +
  scale_color_manual(values = colorp) +
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



legend <- bi_legend(pal = colorp,
                    dim = 4,
                    xlab = "  +  Drought intensity  -  ",
                    ylab = "  -  Flood intensity  +  ",
                    size = 16,
                    arrows = FALSE)

pl=ggarrange(map, legend,
             labels = c("Map", "Key"),
             ncol = 2, nrow = 1,widths = c(2,1), heights=c(1,1), vjust=-1)


pl
ggsave(paste0("D:/tilloal/Documents/LFRuns_utils/TrendAnalysis/plots/totchange_bvPIX_HR_n3.jpg"), pl, width=20, height=20, units=c("cm"),dpi=500)



#2. CLIMATE TREND ------------------


ClimateFloodTrend=FloodTrends[which(FloodTrends$driver==driver[1]),]
ClimateDroughtTrend=DroughtTrends[which(DroughtTrends$driver==driver[1]),]
ClimateFloodTrend$HydroR=FloodTrends[which(FloodTrends$driver==driver[1]),1]
ClimateDroughtTrend$HydroR=DroughtTrends[which(DroughtTrends$driver==driver[1]),1]


FlPoint=full_join(GHshpp,ClimateFloodTrend,by=c("CODEB"="HydroR"))
st_geometry(FlPoint)<-NULL
Flplot=inner_join(HydroRsf,FlPoint,by= c("Id"))
period=c(1951,2020)
haz="flood"

Flplot <- st_transform(Flplot, crs = 3035)

DrPoint=full_join(GHshpp,ClimateDroughtTrend,by=c("CODEB"="HydroR"))
st_geometry(DrPoint)<-NULL
period=c(1951,2020)
haz="flood"
#Drplot <- st_transform(DrPoint, crs = 3035)
Drplot=inner_join(HydroRsf,DrPoint,by= c("Id"))

mv=match(Drplot$Id,Flplot$Id)
plot(diff(mv))
Flplot$d2020=Drplot$Rchange.Y2020
Flplot$f2020=Flplot$Rchange.Y2020

#Univariate plot for verification
colNA="transparent"
# legend="Relative change (%)"
tsize=16
osize=12
legend="Change in Q (%)"
palet=c(hcl.colors(11, palette = "RdYlBu", alpha = NULL, rev = F, fixup = TRUE))

uplot=F

if (uplot==T){
  br=seq(-100,100,2)
  labels=br
  limi=c(-20,20)
  ggplot(basemap) +
    geom_sf(fill="white")+
    geom_sf(data = Flplot, mapping = aes(fill = d2020), alpha=0.9, color = "transparent", size = 0.01, show.legend = TRUE) +
    geom_sf(fill=NA, color="gray42") +
    coord_sf(xlim = c(min(nco[,1]),max(nco[,1])), ylim = c(min(nco[,2]),max(nco[,2])))+
    scale_fill_gradientn(
      colors=palet,
      breaks=br,limits=limi,
      oob = scales::squish,na.value=colNA, name=legend)  +
    labs()+
    theme(axis.title=element_text(size=tsize),
          panel.background = element_rect(fill = "aliceblue", colour = "grey1"),
          panel.border = element_rect(linetype = "solid", fill = NA, colour="black"),
          legend.title = element_text(size=tsize),
          legend.text = element_text(size=osize),
          legend.position = "bottom",
          panel.grid.major = element_line(colour = "grey70"),
          panel.grid.minor = element_line(colour = "grey90"),
          legend.key = element_rect(colour = "transparent"),
          legend.key.size = unit(.8, "cm"))
}

#normalize change
FlplotClim=Flplot
a=sd(FlplotClim$f2020,na.rm=T)
b=sd(FlplotClim$d2020,na.rm=T)
a=b=1
databiclim=FlplotClim
databiclim$x=databiclim$f2020/a
databiclim$y=databiclim$d2020/b

hist(databiclim$x,breaks=100,xlim=c(-5,5))
hist(databiclim$y,breaks=100,xlim=c(-5,5))
sd(databiclim$x)

breaker1=0
breaker2=5

### floods --------------------
alterclass=data.frame(databiclim$x)
alterclass$class=NA
alterclass$class[which(alterclass[,1]<=(-breaker2))]=1
alterclass$class[which(alterclass[,1]>(-breaker2) & alterclass[,1]<(breaker1))]=2
alterclass$class[which(alterclass[,1]>=(breaker1) & alterclass[,1]<(breaker2))]=3
alterclass$class[which(alterclass[,1]>=breaker2)]=4
c1=alterclass$class

### droughts --------------------
alterclass=data.frame(databiclim$y)
alterclass$class=NA
alterclass$class[which(alterclass[,1]<=(-breaker2))]=1
alterclass$class[which(alterclass[,1]>(-breaker2) & alterclass[,1]<(breaker1))]=2
alterclass$class[which(alterclass[,1]>=(breaker1) & alterclass[,1]<(breaker2))]=3
alterclass$class[which(alterclass[,1]>=breaker2)]=4
c2=alterclass$class

cx=paste(c2,c1,sep="-")
databiclim$bi_class=cx
# Combine the main category and subcategory to make a label
databiclim$combined_category <-databiclim$bi_class

merdas=ClimateDroughtTrend$Rchange.Y2020[which(ClimateDroughtTrend$HydroR==2070016940)]
merdav=databiclim$d2020[which(databiclim$Id==2070016940)]
## Pixel level --------------

a=sd(ClimFloodTrendPix$Y2020,na.rm=T)
b=sd(ClimDroughtTrendPix$Y2020,na.rm=T)
a=b=1
databipic=Flpixplot
databipic$x=ClimFloodTrendPix$Y2020/a
databipic$y=ClimDroughtTrendPix$Y2020/b


sd(databipic$x,na.rm=T)
breaker1=0
breaker2=5

alterclass=data.frame(databipic$x)
alterclass$class=NA
alterclass$class[which(alterclass[,1]<=(-breaker2))]=1
alterclass$class[which(alterclass[,1]>(-breaker2) & alterclass[,1]<(breaker1))]=2
alterclass$class[which(alterclass[,1]>=(breaker1) & alterclass[,1]<(breaker2))]=3
alterclass$class[which(alterclass[,1]>=breaker2)]=4

c1=alterclass$class
alterclass=data.frame(databipic$y)
alterclass$class=NA
alterclass$class[which(alterclass[,1]<=(-breaker2))]=1
alterclass$class[which(alterclass[,1]>(-breaker2) & alterclass[,1]<(breaker1))]=2
alterclass$class[which(alterclass[,1]>=(breaker1) & alterclass[,1]<(breaker2))]=3
alterclass$class[which(alterclass[,1]>=breaker2)]=4

c2=alterclass$class

cx=paste(c2,c1,sep="-")

databipic$bi_class=cx

mbicli=data.frame(x=mean(databiclim$d2020,na.rm=T),y=mean(databiclim$f2020,na.rm=T),xsd=sd(databiclim$d2020,na.rm=T),ysd=sd(databiclim$f2020,na.rm=T))

mbicli=data.frame(x=mean(databiclim$d2020,na.rm=T),y=mean(databiclim$f2020,na.rm=T),
                 xq1=quantile(databiclim$d2020,0.05,na.rm=T),yq1=quantile(databiclim$f2020,0.05,na.rm=T),
                 xq2=quantile(databiclim$d2020,0.95,na.rm=T),yq2=quantile(databiclim$f2020,0.95,na.rm=T))



ggplot(databiclim, aes(x = d2020, y = f2020, fill = combined_category)) +
  geom_point(shape = 21, color = "black", size = 3,show.legend = FALSE) +
  #bi_scale_fill(pal = "BlueOr", dim = 3, na.value=colNA )+
  scale_fill_manual(values = colors) +
  coord_cartesian(xlim=c(-100,100),ylim=c(-100,100)) +

  # Add main quadrant labels
  annotate("text", x = 1.8, y = 40, label = "Wetting", color = "gray12", size = 5, fontface = "bold") +
  annotate("text", x = -1.8, y = 40, label = "Accelerating", color = "gray12", size = 5, fontface = "bold") +
  annotate("text", x = -1.8, y = -40, label = "Drying", color = "gray12", size = 5, fontface = "bold") +
  annotate("text", x = 1.8, y = -40, label = "Decelerating", color = "gray12", size = 5, fontface = "bold") +

  # Set axis labels
  labs(x = "Change in drought flows (l/s/km2)",
       y = "Change in flood flows (l/s/km2)") +

  # Customize the theme
  theme(axis.title=element_text(size=tsize),
        panel.background = element_rect(fill = "transparent", colour = "grey1"),
        panel.border = element_rect(linetype = "solid", fill = NA, colour="black"),
        legend.title = element_text(size=tsize),
        legend.text = element_text(size=osize),
        legend.position = "bottom",
        panel.grid.major = element_line(colour = "grey70"),
        panel.grid.minor = element_line(colour = "grey90"),
        legend.key = element_rect(fill = "transparent", colour = "transparent"),
        legend.key.size = unit(.8, "cm"))




map <- ggplot(basemap) +
  geom_sf(fill="white")+
  geom_sf(data = databiclim, mapping = aes(fill = combined_category), alpha=0.7, color = "transparent", size = 0.01,show.legend = F) +
  geom_sf(data = databipic, mapping = aes(col = bi_class,geometry=geometry,size=upa), alpha=1,stroke=0,shape=15, show.legend = FALSE) +
  geom_sf(fill=NA, color="gray42") +
  # bi_scale_fill(pal = "BlueOr", dim = 3, na.value=colNA ) +
  scale_fill_manual(values = colorp) +
  coord_sf(xlim = c(min(nco[,1]),max(nco[,1])), ylim = c(min(nco[,2]),max(nco[,2])))+
  scale_size(range = c(0.08, 0.4), trans="sqrt",name= expression(paste("Upstream area ", (km^2),
                                                                       sep = " ")),
             breaks=c(101,1000,10000,100000,500000), labels=c("100","1000", "10 000", "100 000", "500 000"),
             guide = "none")+
  #bi_scale_color(pal = "BlueOr", dim = 3, na.value=colNA ) +
  scale_color_manual(values = colorp, na.value=NA) +
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



legend <- bi_legend(pal = colorp,
                    dim = 4,
                    xlab = "  -  Drought flows  +  ",
                    ylab = "  -  Flood flows  +  ",
                    size = 16,
                    arrows = FALSE)

pl=ggarrange(map, legend,
             labels = c("Map", "Key"),
             ncol = 2, nrow = 1,widths = c(2,1), heights=c(1,1), vjust=-1)



ggsave(paste0("D:/tilloal/Documents/LFRuns_utils/TrendAnalysis/plots/climchange_HR_new3.jpg"), pl, width=20, height=20, units=c("cm"),dpi=800)


#3. LAND USE TREND --------------

LanduseFloodTrend=FloodTrends[which(FloodTrends$driver==driver[2]),]
LanduseDroughtTrend=DroughtTrends[which(DroughtTrends$driver==driver[2]),]
LanduseFloodTrend$HydroR=FloodTrends[which(FloodTrends$driver==driver[1]),1]
LanduseDroughtTrend$HydroR=DroughtTrends[which(DroughtTrends$driver==driver[1]),1]

LanduseFloodTrend$f2020=LanduseFloodTrend$Rchange.Y2020
LanduseDroughtTrend$f2020=LanduseDroughtTrend$Rchange.Y2020

FlPoint=full_join(GHshpp,LanduseFloodTrend,by=c("CODEB"="HydroR"))
st_geometry(FlPoint)<-NULL
Flplot=inner_join(HydroRsf,FlPoint,by= c("Id"))
period=c(1951,2020)
haz="flood"

Flplot <- st_transform(Flplot, crs = 3035)
DrPoint=full_join(GHshpp,LanduseDroughtTrend,by=c("CODEB"="HydroR"))
st_geometry(DrPoint)<-NULL
Drplot=inner_join(HydroRsf,DrPoint,by= c("Id"))
period=c(1951,2020)
haz="flood"
Drplot <- st_transform(Drplot, crs = 3035)
Flplot$d2020=Drplot$Rchange.Y2020
Flplot$f2020=Flplot$Rchange.Y2020

FlplotLu=Flplot
a=sd(FlplotLu$f2020,na.rm=T)
b=sd(FlplotLu$d2020,na.rm=T)
a=b=1
databilu=FlplotLu
databilu$x=databilu$f2020/a
databilu$y=databilu$d2020/b

breaker1=0
breaker2=5


### floods --------------------
alterclass=data.frame(databilu$x)
alterclass$class=NA
alterclass$class[which(alterclass[,1]<=(-breaker2))]=1
alterclass$class[which(alterclass[,1]>(-breaker2) & alterclass[,1]<(breaker1))]=2
alterclass$class[which(alterclass[,1]>=(breaker1) & alterclass[,1]<(breaker2))]=3
alterclass$class[which(alterclass[,1]>=breaker2)]=4
c1=alterclass$class

### droughts --------------------
alterclass=data.frame(databilu$y)
alterclass$class=NA
alterclass$class[which(alterclass[,1]<=(-breaker2))]=1
alterclass$class[which(alterclass[,1]>(-breaker2) & alterclass[,1]<(breaker1))]=2
alterclass$class[which(alterclass[,1]>=(breaker1) & alterclass[,1]<(breaker2))]=3
alterclass$class[which(alterclass[,1]>=breaker2)]=4
c2=alterclass$class

cx=paste(c2,c1,sep="-")
databilu$bi_class=cx
# Combine the main category and subcategory to make a label
databilu$combined_category <-databilu$bi_class

## Pixel level

luflood=LuFloodTrendPix$Y2020
ludrought=LuDroughtTrendPix$Y2020
a=sd(luflood,na.rm=T)
b=sd(ludrought,na.rm=T)
a=b=1
databipilu=Flpixplot
databipilu$x=luflood/a
databipilu$y=ludrought/b


sd(databipilu$x,na.rm=T)

# databiclim <- bi_class(Flplot, x = x, y = y, style = "fisher", dim = 3, keep_factors = TRUE, dig_lab=2)

# breaker1=0
# breaker2=0.25
# breaker2=10

alterclass=data.frame(databipilu$x)
alterclass$class=NA
alterclass$class[which(alterclass[,1]<=(-breaker2))]=1
alterclass$class[which(alterclass[,1]>(-breaker2) & alterclass[,1]<(breaker1))]=2
alterclass$class[which(alterclass[,1]>=(breaker1) & alterclass[,1]<(breaker2))]=3
alterclass$class[which(alterclass[,1]>=breaker2)]=4

c1=alterclass$class
alterclass=data.frame(databipilu$y)
alterclass$class=NA
alterclass$class[which(alterclass[,1]<=(-breaker2))]=1
alterclass$class[which(alterclass[,1]>(-breaker2) & alterclass[,1]<(breaker1))]=2
alterclass$class[which(alterclass[,1]>=(breaker1) & alterclass[,1]<(breaker2))]=3
alterclass$class[which(alterclass[,1]>=breaker2)]=4

c2=alterclass$class

cx=paste(c2,c1,sep="-")

databipilu$bi_class=cx

mbilu=data.frame(x=mean(databilu$d2020,na.rm=T),y=mean(databilu$f2020,na.rm=T),xsd=sd(databilu$d2020,na.rm=T),ysd=sd(databilu$f2020,na.rm=T))


mbilu=data.frame(x=mean(databilu$d2020,na.rm=T),y=mean(databilu$f2020,na.rm=T),
                  xq1=quantile(databilu$d2020,0.05,na.rm=T),yq1=quantile(databilu$f2020,0.05,na.rm=T),
                  xq2=quantile(databilu$d2020,0.95,na.rm=T),yq2=quantile(databilu$f2020,0.95,na.rm=T))

# 
# 
# ggplot(databilu, aes(x = d2020, y = f2020, fill = combined_category)) +
#   geom_point(shape = 21, color = "black", size = 3,show.legend = FALSE) +
#   #bi_scale_fill(pal = "BlueOr", dim = 3, na.value=colNA )+
#   scale_fill_manual(values = colors) +
#   coord_cartesian(xlim=c(-100,100),ylim=c(-100,100)) +
#   geom_point(data=mbilu, aes(x = x,y = y), fill="red",color="red",size=4)+
#   
#   # Add main quadrant labels
#   annotate("text", x = 1.8, y = 40, label = "Wetting", color = "gray12", size = 5, fontface = "bold") +
#   annotate("text", x = -1.8, y = 40, label = "Accelerating", color = "gray12", size = 5, fontface = "bold") +
#   annotate("text", x = -1.8, y = -40, label = "Drying", color = "gray12", size = 5, fontface = "bold") +
#   annotate("text", x = 1.8, y = -40, label = "Decelerating", color = "gray12", size = 5, fontface = "bold") +
#   
#   # Set axis labels
#   labs(x = "Change in drought flows (l/s/km2)", 
#        y = "Change in flood flows (l/s/km2)") +
#   
#   # Customize the theme
#   theme(axis.title=element_text(size=tsize),
#         panel.background = element_rect(fill = "transparent", colour = "grey1"),
#         panel.border = element_rect(linetype = "solid", fill = NA, colour="black"),
#         legend.title = element_text(size=tsize),
#         legend.text = element_text(size=osize),
#         legend.position = "bottom",
#         panel.grid.major = element_line(colour = "grey70"),
#         panel.grid.minor = element_line(colour = "grey90"),
#         legend.key = element_rect(fill = "transparent", colour = "transparent"),
#         legend.key.size = unit(.8, "cm"))
# 
# 
# 
# 
# map <- ggplot(basemap) +
#   geom_sf(fill="white")+
#   geom_sf(data = databilu, mapping = aes(fill = combined_category), alpha=0.7, color = "transparent", size = 0.01,show.legend = F) +
#   geom_sf(data = databipilu, mapping = aes(col = bi_class,geometry=geometry,size=upa), alpha=1,stroke=0,shape=15, show.legend = FALSE) +
#   geom_sf(fill=NA, color="gray42") +
#   # bi_scale_fill(pal = "BlueOr", dim = 3, na.value=colNA ) +
#   scale_fill_manual(values = colors) +
#   coord_sf(xlim = c(min(nco[,1]),max(nco[,1])), ylim = c(min(nco[,2]),max(nco[,2])))+
#   scale_size(range = c(0.08, 0.4), trans="sqrt",name= expression(paste("Upstream area ", (km^2),
#                                                                        sep = " ")),
#              breaks=c(101,1000,10000,100000,500000), labels=c("100","1000", "10 000", "100 000", "500 000"),
#              guide = "none")+
#   #bi_scale_color(pal = "BlueOr", dim = 3, na.value=colNA ) +
#   scale_color_manual(values = colors) +
#   labs()+
#   theme(axis.title=element_text(size=tsize),
#         panel.background = element_rect(fill = "aliceblue", colour = "grey1"),
#         panel.border = element_rect(linetype = "solid", fill = NA, colour="black"),
#         legend.title = element_text(size=tsize),
#         legend.text = element_text(size=osize),
#         legend.position = "bottom",
#         panel.grid.major = element_line(colour = "grey70"),
#         panel.grid.minor = element_line(colour = "grey90"),
#         legend.key = element_rect(fill = "transparent", colour = "transparent"),
#         legend.key.size = unit(.8, "cm"))
# 
# 
# 
# legend <- bi_legend(pal = colors,
#                     dim = 4,
#                     xlab = "  +  Drought changes  -  ",
#                     ylab = "  -  Flood changes  +  ",
#                     size = 16,
#                     arrows = FALSE)
# 
# pl=ggarrange(map, legend, 
#              labels = c("Map", "Key"),
#              ncol = 2, nrow = 1,widths = c(2,1), heights=c(1,1), vjust=-1)
# 
# pl
# 
# ggsave(paste0("D:/tilloal/Documents/LFRuns_utils/TrendAnalysis/plots/luchange_bvrel10p.jpg"), pl, width=20, height=20, units=c("cm"),dpi=800) 
# 


#4. RESERVOIR trend --------------

# Set up the Reservoir
ReservoirFloodTrend=FloodTrends[which(FloodTrends$driver==driver[3]),c(2:71)]
ReservoirDroughtTrend=DroughtTrends[which(DroughtTrends$driver==driver[3]),c(2:71)]
ReservoirFloodTrend$HydroR=FloodTrends[which(FloodTrends$driver==driver[1]),1]
ReservoirDroughtTrend$HydroR=DroughtTrends[which(DroughtTrends$driver==driver[1]),1]

FlPoint=full_join(GHshpp,ReservoirFloodTrend,by=c("CODEB"="HydroR"))
st_geometry(FlPoint)<-NULL
Flplot=inner_join(HydroRsf,FlPoint,by= c("Id"))
period=c(1951,2020)
haz="flood"

Flplot <- st_transform(Flplot, crs = 3035)

DrPoint=full_join(GHshpp,ReservoirDroughtTrend,by=c("CODEB"="HydroR"))
st_geometry(DrPoint)<-NULL
Drplot=inner_join(HydroRsf,DrPoint,by= c("Id"))
period=c(1951,2020)
haz="flood"
Drplot <- st_transform(Drplot, crs = 3035)

Flplot$d2020=Drplot$Rchange.Y2020
Flplot$f2020=Flplot$Rchange.Y2020

FlplotRes=Flplot
a=sd(FlplotRes$f2020,na.rm=T)
b=sd(FlplotRes$d2020,na.rm=T)
a=b=1
databire=FlplotRes
databire$x=databire$f2020/a
databire$y=databire$d2020/b

# breaker1=0
# breaker2=0.25
# breaker2=10

alterclass=data.frame(databire$x)
alterclass$class=NA
alterclass$class[which(alterclass[,1]<=(-breaker2))]=1
alterclass$class[which(alterclass[,1]>(-breaker2) & alterclass[,1]<(breaker1))]=2
alterclass$class[which(alterclass[,1]>=(breaker1) & alterclass[,1]<(breaker2))]=3
alterclass$class[which(alterclass[,1]>=breaker2)]=4

c1=alterclass$class
alterclass=data.frame(databire$y)
alterclass$class=NA
alterclass$class[which(alterclass[,1]<=(-breaker2))]=1
alterclass$class[which(alterclass[,1]>(-breaker2) & alterclass[,1]<(breaker1))]=2
alterclass$class[which(alterclass[,1]>=(breaker1) & alterclass[,1]<(breaker2))]=3
alterclass$class[which(alterclass[,1]>=breaker2)]=4

c2=alterclass$class

cx=paste(c2,c1,sep="-")

databire$bi_class=cx
# Combine the main category and subcategory to make a label
databire$combined_category <-databire$bi_class


#Reservoir signal
crsiflood=ResFloodTrendPix$Y2020
crsidrought=ResDroughtTrendPix$Y2020
a=sd(crsiflood,na.rm=T)
b=sd(crsidrought,na.rm=T)
a=b=1
databipire=Flpixplot
databipire$x=crsiflood/a
databipire$y=crsidrought/b


sd(databipire$x,na.rm=T)

# databiclim <- bi_class(Flplot, x = x, y = y, style = "fisher", dim = 3, keep_factors = TRUE, dig_lab=2)

# breaker1=0
# breaker2=0.25
# breaker2=10

alterclass=data.frame(databipire$x)
alterclass$class=NA
alterclass$class[which(alterclass[,1]<=(-breaker2))]=1
alterclass$class[which(alterclass[,1]>(-breaker2) & alterclass[,1]<(breaker1))]=2
alterclass$class[which(alterclass[,1]>=(breaker1) & alterclass[,1]<(breaker2))]=3
alterclass$class[which(alterclass[,1]>=breaker2)]=4

c1=alterclass$class
alterclass=data.frame(databipire$y)
alterclass$class=NA
alterclass$class[which(alterclass[,1]<=(-breaker2))]=1
alterclass$class[which(alterclass[,1]>(-breaker2) & alterclass[,1]<(breaker1))]=2
alterclass$class[which(alterclass[,1]>=(breaker1) & alterclass[,1]<(breaker2))]=3
alterclass$class[which(alterclass[,1]>=breaker2)]=4

c2=alterclass$class

cx=paste(c2,c1,sep="-")

databipire$bi_class=cx

mbires=data.frame(x=mean(databire$d2020,na.rm=T),y=mean(databire$f2020,na.rm=T),
                  xq1=quantile(databire$d2020,0.05,na.rm=T),yq1=quantile(databire$f2020,0.05,na.rm=T),
                  xq2=quantile(databire$d2020,0.95,na.rm=T),yq2=quantile(databire$f2020,0.95,na.rm=T))


#5. Water Demand signal ------------------


WaterDemandFloodTrend=FloodTrends[which(FloodTrends$driver==driver[4]),c(2:71)]
WaterDemandDroughtTrend=DroughtTrends[which(DroughtTrends$driver==driver[4]),c(2:71)]
WaterDemandFloodTrend$HydroR=FloodTrends[which(FloodTrends$driver==driver[4]),1]
WaterDemandDroughtTrend$HydroR=DroughtTrends[which(DroughtTrends$driver==driver[4]),1]

FlPoint=full_join(GHshpp,WaterDemandFloodTrend,by=c("CODEB"="HydroR"))
st_geometry(FlPoint)<-NULL
Flplot=inner_join(HydroRsf,FlPoint,by= c("Id"))
period=c(1951,2020)
haz="flood"

Flplot <- st_transform(Flplot, crs = 3035)

DrPoint=full_join(GHshpp,WaterDemandDroughtTrend,by=c("CODEB"="HydroR"))
st_geometry(DrPoint)<-NULL
Drplot=inner_join(HydroRsf,DrPoint,by= c("Id"))
period=c(1951,2020)
haz="flood"
Drplot <- st_transform(Drplot, crs = 3035)

Flplot$d2020=Drplot$Rchange.Y2020
Flplot$f2020=Flplot$Rchange.Y2020

FlplotWD=Flplot
a=sd(FlplotWD$f2020,na.rm=T)
b=sd(FlplotWD$d2020,na.rm=T)
a=b=1
databiwd=FlplotWD
databiwd$x=databiwd$f2020/a
databiwd$y=databiwd$d2020/b

# breaker1=0
# breaker2=0.25
# breaker2=10

alterclass=data.frame(databiwd$x)
alterclass$class=NA
alterclass$class[which(alterclass[,1]<=(-breaker2))]=1
alterclass$class[which(alterclass[,1]>(-breaker2) & alterclass[,1]<(breaker1))]=2
alterclass$class[which(alterclass[,1]>=(breaker1) & alterclass[,1]<(breaker2))]=3
alterclass$class[which(alterclass[,1]>=breaker2)]=4

c1=alterclass$class
alterclass=data.frame(databiwd$y)
alterclass$class=NA
alterclass$class[which(alterclass[,1]<=(-breaker2))]=1
alterclass$class[which(alterclass[,1]>(-breaker2) & alterclass[,1]<(breaker1))]=2
alterclass$class[which(alterclass[,1]>=(breaker1) & alterclass[,1]<(breaker2))]=3
alterclass$class[which(alterclass[,1]>=breaker2)]=4

c2=alterclass$class

cx=paste(c2,c1,sep="-")

databiwd$bi_class=cx
# Combine the main category and subcategory to make a label
databiwd$combined_category <-databiwd$bi_class


## Pixel level ---------------
crsiflood=WuFloodTrendPix$Y2020
crsidrought=WuDroughtTrendPix$Y2020
a=sd(crsiflood,na.rm=T)
b=sd(crsidrought,na.rm=T)
a=b=1
databipiwd=Flpixplot
databipiwd$x=crsiflood/a
databipiwd$y=crsidrought/b


sd(databipiwd$x,na.rm=T)

# databiclim <- bi_class(Flplot, x = x, y = y, style = "fisher", dim = 3, keep_factors = TRUE, dig_lab=2)

# breaker1=0
# breaker2=0.25
# breaker2=10

alterclass=data.frame(databipiwd$x)
alterclass$class=NA
alterclass$class[which(alterclass[,1]<=(-breaker2))]=1
alterclass$class[which(alterclass[,1]>(-breaker2) & alterclass[,1]<(breaker1))]=2
alterclass$class[which(alterclass[,1]>=(breaker1) & alterclass[,1]<(breaker2))]=3
alterclass$class[which(alterclass[,1]>=breaker2)]=4

c1=alterclass$class
alterclass=data.frame(databipiwd$y)
alterclass$class=NA
alterclass$class[which(alterclass[,1]<=(-breaker2))]=1
alterclass$class[which(alterclass[,1]>(-breaker2) & alterclass[,1]<(breaker1))]=2
alterclass$class[which(alterclass[,1]>=(breaker1) & alterclass[,1]<(breaker2))]=3
alterclass$class[which(alterclass[,1]>=breaker2)]=4

c2=alterclass$class

cx=paste(c2,c1,sep="-")

databipiwd$bi_class=cx

#mbiwd=data.frame(x=mean(databiwd$d2020,na.rm=T),y=mean(databiwd$f2020,na.rm=T),xsd=sd(databiwd$d2020,na.rm=T),ysd=sd(databiwd$f2020,na.rm=T))


mbiwd=data.frame(x=mean(databiwd$d2020,na.rm=T),y=mean(databiwd$f2020,na.rm=T),
                  xq1=quantile(databiwd$d2020,0.05,na.rm=T),yq1=quantile(databiwd$f2020,0.05,na.rm=T),
                  xq2=quantile(databiwd$d2020,0.95,na.rm=T),yq2=quantile(databiwd$f2020,0.95,na.rm=T))


# Bivariate plots of mean contribution ----------------------



biogeo <- read_sf(dsn = paste0(hydroDir,"/eea_3035_biogeo-regions_2016/BiogeoRegions2016_wag84.shp"))
biogeof=fortify(biogeo)
st_geometry(biogeof)<-NULL
biogeoregions=raster( paste0(hydroDir,"/eea_3035_biogeo-regions_2016/Biogeo_rasterized_wsg84.tif"))
Gbiogeoregions=as.data.frame(biogeoregions,xy=T)
biogeomatch=inner_join(biogeof,Gbiogeoregions,by= c("PK_UID"="Biogeo_rasterized_wsg84"))
biogeomatch$latlong=paste(round(biogeomatch$x,4),round(biogeomatch$y,4),sep=" ")
biogeo_rivers=right_join(biogeomatch,outf, by="latlong")



#match biogeoregions with HRs
databitotxHR=databipi
RegioRLi=aggregate(list(val=databitotxHR$y),
                   by = list(HydroR=databitotxHR$HydroRegions_raster_WGS84),
                   FUN = function(x) c(mean=mean(x,na.rm=T),dev=sd(x,na.rm=T),len=length(x),med=median(x,na.rm=T),q1=quantile(x, 0.05, na.rm=T),q3=quantile(x, 0.95, na.rm=T)))
RegioRLi <- do.call(data.frame, RegioRLi)

HRM=na.omit((match(RegioRLi$HydroR,HydroRsf$Id)))
HydroRsf_dom=HydroRsf[HRM,]


layerMatch=match(GHR_riv$outl2,biogeo_rivers$outl2)

GHR_riv$biogeoR=biogeo_rivers$code[layerMatch]
unique(GHR_riv$HydroRegions_raster_WGS84)
verif=(match(GHR_riv$HydroRegions_raster_WGS84,HydroRsf_dom$Id))
unique(verif)
unique(HydroRsf_dom$Id[verif])
GHR_riv$HydrRName=HydroRsf_dom$IRST_NAMEB[verif]


bhp_m=na.omit(unique(match(HydroRsf_dom$Id,GHR_riv$HydroRegions_raster_WGS84)))
HydroRsf_dom$biogeoreg=GHR_riv$biogeoR[bhp_m]


sample=GHR_riv[bhp_m,]


RegioAg = aggregate(list(oc=GHR_riv$HydrRName),
                    by = list(BR=GHR_riv$biogeoR,HR=GHR_riv$HydroRegions_raster_WGS84),
                    FUN = function(x) c(len=length(x)))
RegioAg <- do.call(data.frame, RegioAg)

ur=unique(RegioAg$HR)

#matching biogeoregions and hybas07
databiwdxBG=inner_join(biogeo_rivers,databipiwd,by=c("outl2"))
databiluxBG=inner_join(biogeo_rivers,databipilu,by=c("outl2"))
databiclixBG=inner_join(biogeo_rivers,databipic,by=c("outl2"))
databiresxBG=inner_join(biogeo_rivers,databipire,by=c("outl2"))
databitotxBG=inner_join(biogeo_rivers,databipi,by=c("outl2"))


#matching hydroregions and hybas07
databiwdxHR=inner_join(GHR_riv,databipiwd,by=c("outl2"))
databiluxHR=inner_join(GHR_riv,databipilu,by=c("outl2"))
databiclixHR=inner_join(GHR_riv,databipic,by=c("outl2"))
databiresxHR=inner_join(GHR_riv,databipire,by=c("outl2"))
databitotxHR=inner_join(GHR_riv,databipi,by=c("outl2"))
# 
# databiwdxHR=databipiwd
# databiluxHR=databipilu
# databiclixHR=databipic
# databiresxHR=databipire
# databitotxHR=databipi

mbf=rbind(mbicli,mbires,mbilu,mbiwd)
mbf$names=c('Climate',
            'Reservoirs',
            'Landuse',
            'WaterDemand')
# mbf$lbx=mbf$x-mbf$xsd
# mbf$hbx=mbf$x+mbf$xsd
# mbf$lby=mbf$y-mbf$ysd
# mbf$hby=mbf$x+mbf$ysd

loscolors=c("Accelerating" = "#174f28","Drying" = "#dd6a29","Stable"="gray60","Wetting" = "#169dd0","Decelerating" = "burlywood")
clabels=c("Climate","Land use","Reservoirs", "Water demand")
colorn = c("WaterDemand" ='limegreen',"Reservoirs" ='tomato4',"Landuse" ='orange',"Climate" ='royalblue')
tsize=osize=20
bpc=ggplot(mbf, aes(x = x, y = y, fill = names)) +
  geom_linerange(data=mbf,aes(x=x, ymin=yq1,ymax=yq2,color = names,group=factor(names)),lwd=1,alpha=0.8) +
  geom_linerange(data=mbf,aes(y=y, xmin=xq1,xmax=xq2,color = names,group=factor(names)),lwd=1,alpha=0.8) +
  geom_point(shape = 21, color = "black", size = 5, show.legend = TRUE) +
  scale_fill_manual(values = colorn, name = "Drivers", labels = clabels) +
  scale_color_manual(values = colorn, name = "Drivers", labels = clabels) +
  coord_cartesian(xlim=c(-100,100),ylim=c(-50,50)) +
  annotate("text", x = 50, y = 20, label = "Wetting", color = "#169dd0", size = 7, fontface = "bold") +
  annotate("text", x = -50, y = 20, label = "Accelerating", color = "#174f28", size = 7, fontface = "bold") +
  annotate("text", x = -50, y = -20, label = "Drying", color = "#dd6a29", size = 7, fontface = "bold") +
  annotate("text", x = 50, y = -20, label = "Decelerating", color = "burlywood", size = 7, fontface = "bold") +
  labs(x = "Change in drought (%)", 
       y = "Change in flood (%)") +
  scale_x_continuous(trans=scales::modulus_trans(.5),
                     breaks=seq(-100,100,50),
                     minor_breaks = seq(-100,100,10)) +
  scale_y_continuous(trans=scales::modulus_trans(.5),
                     breaks=seq(-100,100,25),
                     minor_breaks = seq(-100,100,10)) +
  guides( fill = guide_legend(override.aes = list(size = 6)))+
  theme(
    axis.title = element_text(size = 18, face = "bold"),
    title = element_text(size = 22, face = "bold"),
    axis.text = element_text(size = 16, face = "bold"),
    panel.background = element_rect(fill = "white", colour = "white"),
    panel.grid = element_blank(),
    panel.border = element_rect(linetype = "solid", fill = NA, colour = "black"),
    legend.title = element_text(size = 20, face = "bold"),
    legend.text = element_text(size = 16),
    #legend.key = element_rect(fill = "lightgray", colour = "lightgray"),
    legend.position = "right",
    panel.grid.major = element_line(colour = "grey60"),
    panel.grid.minor= element_line(colour = "grey70", linetype = "dashed"),
    #legend.key = element_rect(fill = "white", colour = "white"),
    #legend.key.size = unit(1, "cm")
  )

bpc
ggsave(paste0("D:/tilloal/Documents/LFRuns_utils/TrendAnalysis/plots/Contribution_bv.jpg"),bpc, width=25, height=20, units=c("cm"),dpi=300) 


length(which(is.na(databiresxHR$x.y)))

# databiwdxHR$y.y[which(abs(databiwdxHR$y.y)>100)]=100
mbiwd1=aggregate(list(val=databiwdxHR$y.y),
                 by = list(HR=databiwdxHR$HydrRName),
                 FUN = function(x) c(dmean=mean(x,na.rm=T),dmed=median(x,na.rm=T),dq1=quantile(x,0.05,na.rm=T),dq2=quantile(x,0.95,na.rm=T),l=length(x),sd=sd(x,na.rm=T)))
mbiwd1 <- do.call(data.frame, mbiwd1)

mbicli1=aggregate(list(val=databiclixHR$y.y),
                  by = list(HR=databiclixHR$HydrRName),
                  FUN = function(x) c(dmean=mean(x,na.rm=T),dmed=median(x,na.rm=T),dq1=quantile(x,0.05,na.rm=T),dq2=quantile(x,0.95,na.rm=T),l=length(x),sd=sd(x,na.rm=T)))
mbicli1 <- do.call(data.frame, mbicli1)

#databiresxHR$y[which(abs(databiresxHR$y)==0)]=NA
mbires1=aggregate(list(val=databiresxHR$y.y),
                  by = list(HR=databiresxHR$HydrRName),
                  FUN = function(x) c(dmean=mean(x,na.rm=T),dmed=median(x,na.rm=T),dq1=quantile(x,0.05,na.rm=T),dq2=quantile(x,0.95,na.rm=T),l=length(x),sd=sd(x,na.rm=T)))
mbires1 <- do.call(data.frame, mbires1)

# mbires1$se=mbires1$val.sd/sqrt(mbires1$val.l)
# databiluxHR$y.y[which(abs(databiluxHR$y.y)>100)]=100
mbilu1=aggregate(list(val=databiluxHR$y.y),
                 by = list(HR=databiluxHR$HydrRName),
                 FUN = function(x) c(dmean=mean(x,na.rm=T),dmed=median(x,na.rm=T),dq1=quantile(x,0.05,na.rm=T),dq2=quantile(x,0.95,na.rm=T),l=length(x),sd=sd(x,na.rm=T)))
mbilu1 <- do.call(data.frame, mbilu1)


#databitotxHR$y.y[which(abs(databitotxHR$y.y)>100)] =100
mball1=aggregate(list(val=databitotxHR$y.y),
                 by = list(HR=databitotxHR$HydrRName),
                 FUN = function(x) c(dmean=mean(x,na.rm=T),dmed=median(x,na.rm=T),dq1=quantile(x,0.05,na.rm=T),dq2=quantile(x,0.95,na.rm=T),l=length(x),sd=sd(x,na.rm=T)))
mball1 <- do.call(data.frame, mball1)
#mball1$se=mball1$val.sd/sqrt(mball1$val.l)  


#databiwdxHR$x.y[which(abs(databiwdxHR$x.y)<=1e-3)]=NA
mbiwd2=aggregate(list(val=databiwdxHR$x.y),
                 by = list(HR=databiwdxHR$HydrRName),
                 FUN = function(x) c(fmean=mean(x,na.rm=T),fq1=quantile(x,0.05,na.rm=T),fq2=quantile(x,0.95,na.rm=T),l=length(x),sd=sd(x,na.rm=T)))
mbiwd2 <- do.call(data.frame, mbiwd2)

mbicli2=aggregate(list(val=databiclixHR$x.y),
                  by = list(HR=databiclixHR$HydrRName),
                  FUN = function(x) c(fmean=mean(x,na.rm=T),fq1=quantile(x,0.05,na.rm=T),fq2=quantile(x,0.95,na.rm=T),l=length(x),sd=sd(x,na.rm=T)))
mbicli2 <- do.call(data.frame, mbicli2)

# databiresxHR$x.y[which(abs(databiresxHR$x.y)>1e4)]=1e3
mbires2=aggregate(list(val=databiresxHR$x.y),
                  by = list(HR=databiresxHR$HydrRName),
                  FUN = function(x) c(fmean=mean(x,na.rm=T),fq1=quantile(x,0.05,na.rm=T),fq2=quantile(x,0.95,na.rm=T),l=length(x),sd=sd(x,na.rm=T)))
mbires2 <- do.call(data.frame, mbires2)

# databiluxHR$x.y[which(abs(databiluxHR$x.y)<=1e-3)]=NA
mbilu2=aggregate(list(val=databiluxHR$x.y),
                 by = list(HR=databiluxHR$HydrRName),
                 FUN = function(x) c(fmean=mean(x,na.rm=T),fq1=quantile(x,0.05,na.rm=T),fq2=quantile(x,0.95,na.rm=T),l=length(x),sd=sd(x,na.rm=T)))
mbilu2 <- do.call(data.frame, mbilu2)

mball2=aggregate(list(val=databitotxHR$x.y),
                 by = list(HR=databitotxHR$HydrRName),
                 FUN = function(x) c(fmean=mean(x,na.rm=T),fq1=quantile(x,0.05,na.rm=T),fq2=quantile(x,0.95,na.rm=T),l=length(x),sd=sd(x,na.rm=T)))
mball2 <- do.call(data.frame, mball2)

mbiwd3=inner_join(mbiwd1,mbiwd2,by="HR")
mbires3=inner_join(mbires1,mbires2,by="HR")
mbilu3=inner_join(mbilu1,mbilu2,by="HR")
mbicli3=inner_join(mbicli1,mbicli2,by="HR")
mball3=inner_join(mball1,mball2,by="HR")

RegionName=inner_join(mbicli2,HydroRsf,by=c("HR"="IRST_NAMEB"))

l1=length(mbicli2$HR)

mbfX=rbind(mbicli3,mbires3,mbilu3,mbiwd3,mball3)
mbfX$names=c(rep('Climate',l1),
            rep('Reservoirs',l1),
            rep('Landuse',l1),
            rep('WaterDemand',l1),
            rep('Historical',l1))

zeb=match(mbfX$HR,RegionName$HR)
mbfX$regioname=RegionName$IRST_NAMEB[zeb]


#correlation between change in drought and flood of several drivers

cor.test(mbicli3$val.dmean, mbicli3$val.fmean)
cor.test(mbires3$val.dmean, mbires3$val.fmean)
cor.test(mbilu3$val.dmean, mbilu3$val.fmean)
# mbf$lbx=mbf$x-mbf$xsd
# mbf$hbx=mbf$x+mbf$xsd
# mbf$lby=mbf$y-mbf$ysd
# mbf$hby=mbf$x+mbf$ysd

# breg=unique(mbfX$HR)
# breg=breg[c(1,3,5,6,7)]
# mbk=which(!is.na(match(mbfX$bg,breg)))
# mbfX=mbfX[mbk,]


colnames(mbfX)[c(2,4,5,8,9,10)]=c("x","xq1","xq2","y","yq1","yq2")
mbfX$xl=((mbfX$xq2-mbfX$xq1)+(mbfX$yq2-mbfX$yq1))/2

mbfH=mbfX[which(mbfX$names=="Historical"),]
mbfX=mbfX[-which(mbfX$names=="Historical"),]

loscolors=c("Accelerating" = "#174f28","Drying" = "#dd6a29","Stable"="gray60","Wetting" = "#169dd0","Decelerating" = "burlywood")
clabels=c("Climate","Land use","Reservoirs", "Water demand")
colorn = c("WaterDemand" ='limegreen',"Reservoirs" ='tomato4',"Landuse" ='orange',"Climate" ='royalblue')
shapex=c("Alpine"=0,"Atlantic"=1,"Boreal"=2, "Continental"=3,"Mediterranean"=4)
tsize=osize=20
lalim=60
bpc=ggplot() +
  geom_vline(xintercept = 0, 
  color = "black", size = 1.5) +
  geom_hline(yintercept = 0, 
             color = "black", size = 1.5) +
  #geom_linerange(data=mbfH,aes(x=x, ymin=yq1,ymax=yq2,color = names,group=factor(names)),lwd=2,alpha=0.5) +
  #geom_linerange(data=mbfH,aes(y=y, xmin=xq1,xmax=xq2,color = names,group=factor(names)),lwd=2,alpha=0.5) +
  # geom_point(data=mbf, aes(x = x, y = y, fill = names), shape = 21, color = "black", size = 5, show.legend = FALSE) +
  geom_point(data=mbfX, aes(x = x, y = y, color = names, size=val.l.x), show.legend = TRUE, stroke=1, alpha=.5) +
  #scale_fill_manual(values = colorn, name = "Drivers", labels = clabels) +
  scale_color_manual(values = colorn, name = "Drivers", labels = clabels) +
  #scale_shape_manual(values = shapex, name = "Regions") +
  scale_size(range = c(2, 6))+
  coord_cartesian(xlim=c(-lalim,lalim),ylim=c(-lalim,lalim)) +
  annotate("text", x = 35, y = 45, label = "Wetting", color = "#169dd0", size = 6, fontface = "bold") +
  annotate("text", x = -35, y = 45, label = "Accelerating", color = "#174f28", size = 6, fontface = "bold") +
  annotate("text", x = -35, y = -45, label = "Drying", color = "#dd6a29", size = 6, fontface = "bold") +
  annotate("text", x = 35, y = -45, label = "Decelerating", color = "burlywood", size = 6, fontface = "bold") +
  labs(x = "Change in drought intensity (%)", 
       y = "Change in flood intensity (%)") +
  scale_x_continuous(trans=scales::modulus_trans(.5),
                     breaks=c(-200,-100,-50,-20,-10,-5,0,5,10,20,50,100,200),
                     labels=rev(c(-200,-100,-50,-20,-10,-5,0,5,10,20,50,100,200)),
                     minor_breaks = seq(-100,100,10)) +
  scale_y_continuous(trans=scales::modulus_trans(.5),
                     breaks=c(-200,-100,-50,-20,-10,-5,0,5,10,20,50,100,200),
                     minor_breaks = seq(-100,100,10), expand = c(0,0)) +
  guides( color = guide_legend(override.aes = list(size = 6)),size="none")+
  theme(
    axis.title = element_text(size = 18, face = "bold"),
    title = element_text(size = 22, face = "bold"),
    axis.text = element_text(size = 16, face = "bold"),
    panel.background = element_rect(fill = "white", colour = "white"),
    panel.grid = element_blank(),
    panel.border = element_rect(linetype = "solid", fill = NA, colour = "black"),
    legend.title = element_text(size = 20, face = "bold"),
    legend.text = element_text(size = 16),
    #legend.key = element_rect(fill = "lightgray", colour = "lightgray"),
    legend.position = "right",
    panel.grid.major = element_line(colour = "grey60"),
    panel.grid.minor= element_line(colour = "grey70", linetype = "dashed"),
    #legend.key = element_rect(fill = "white", colour = "white"),
    #legend.key.size = unit(1, "cm")
  )

bpc
ggsave(paste0("D:/tilloal/Documents/LFRuns_utils/TrendAnalysis/plots/Contribution_HR_new3.jpg"),bpc, width=25, height=20, units=c("cm"),dpi=300) 



#I can also plot all catchment values and color them by region, although it may be a mess

mb1=mbfX[which(mbfX$bg=="Mediterranean"),]
mb1=mbfX[which(mbfX$bg=="Continental"),]

bpc=ggplot(mb1, aes(x = x, y = y, fill = names)) +
  geom_point(data=mb1, aes(x = x, y = y, color = names, size=val.l.x), show.legend = TRUE, stroke=1, alpha=.5) +
  #scale_fill_manual(values = colorn, name = "Drivers", labels = clabels) +
  scale_color_manual(values = colorn, name = "Drivers", labels = clabels) +
  #scale_shape_manual(values = shapex, name = "Regions") +
  coord_cartesian(xlim=c(-100,100),ylim=c(-50,50)) +
  annotate("text", x = 50, y = 20, label = "Wetting", color = "#169dd0", size = 7, fontface = "bold") +
  annotate("text", x = -50, y = 20, label = "Accelerating", color = "#174f28", size = 7, fontface = "bold") +
  annotate("text", x = -50, y = -20, label = "Drying", color = "#dd6a29", size = 7, fontface = "bold") +
  annotate("text", x = 50, y = -20, label = "Decelerating", color = "burlywood", size = 7, fontface = "bold") +
  labs(x = "Change in drought (%)", 
       y = "Change in flood (%)") +
  scale_x_continuous(trans=scales::modulus_trans(.5),
                     breaks=seq(-100,100,50),
                     minor_breaks = seq(-100,100,10)) +
  scale_y_continuous(trans=scales::modulus_trans(.5),
                     breaks=seq(-100,100,25),
                     minor_breaks = seq(-100,100,10)) +
  guides( fill = guide_legend(override.aes = list(size = 6)))+
  theme(
    axis.title = element_text(size = 18, face = "bold"),
    title = element_text(size = 22, face = "bold"),
    axis.text = element_text(size = 16, face = "bold"),
    panel.background = element_rect(fill = "white", colour = "white"),
    panel.grid = element_blank(),
    panel.border = element_rect(linetype = "solid", fill = NA, colour = "black"),
    legend.title = element_text(size = 20, face = "bold"),
    legend.text = element_text(size = 16),
    #legend.key = element_rect(fill = "lightgray", colour = "lightgray"),
    legend.position = "right",
    panel.grid.major = element_line(colour = "grey60"),
    panel.grid.minor= element_line(colour = "grey70", linetype = "dashed"),
    #legend.key = element_rect(fill = "white", colour = "white"),
    #legend.key.size = unit(1, "cm")
  )

bpc
ggsave(paste0("D:/tilloal/Documents/LFRuns_utils/TrendAnalysis/plots/Contribution_bv.jpg"),bpc, width=25, height=20, units=c("cm"),dpi=300) 

length(unique(GHshpp$IRST_NAMEB))

#plot of historical change with colors by biogeoregions

#fin fix or just use something else as colors
Hydroplot=c()
for (id in 1:length(ur)){
  myr=ur[id]
  maty=which(!is.na((match(RegioAg$HR,myr))))
  lar=RegioAg[maty,]
  if(length(lar$oc)>1){
    sumi=sum(lar$oc)
    lar=lar[which.max(lar$oc),]
    lar$oc=lar$oc/sumi*100
  }else{  lar$oc=100}
  Hydroplot=rbind(Hydroplot,lar)
} 

bioplot <- inner_join(HydroRsf_dom,Hydroplot, by=c("Id"="HR"))
bioplot$BR[which(bioplot$BR=="Pannonian")]="Continental"
bioplot$BR[which(bioplot$BR=="Steppic")]="Continental"
#nutplot=nutplot[-which(is.na(nutplot$maxcol)),]
br=c(-50,-20,-10,0,10,20,50)
labels=br
limi=c(-60,60)
colx=c("Alpine"="deeppink","Atlantic"="forestgreen","Boreal"="darkviolet", "Continental"="orange3","Mediterranean"="gold")

ggplot(basemap) +
  geom_sf(fill="gray95")+
  geom_sf(fill=NA, color="grey") +
  geom_sf(data=bioplot,aes(fill=factor(BR),geometry=geometry),color="black",alpha=.7,size=0.0001)+ 
  coord_sf(xlim = c(min(nco[,1]),max(nco[,1])), ylim = c(min(nco[,2]),max(nco[,2])))+
  scale_fill_manual(values = colx, name="BioGeoRegions") +
  labs(x="Longitude", y = "Latitude")+
  guides(color = guide_colourbar(barwidth = 12, barheight = 1))+
  theme(axis.title=element_text(size=tsize),
        panel.background = element_rect(fill = "aliceblue", colour = "grey1"),
        panel.border = element_rect(linetype = "solid", fill = NA, colour="black"),
        legend.title = element_text(size=tsize),
        legend.text = element_text(size=osize),
        legend.position = "right",
        panel.grid.major = element_line(colour = "grey70"),
        panel.grid.minor = element_line(colour = "grey90"),
        legend.key = element_rect(fill = "transparent", colour = "transparent"),
        legend.key.size = unit(.8, "cm"))

ggsave(paste0("D:/tilloal/Documents/LFRuns_utils/TrendAnalysis/plots/Map_HRxBG2.jpg"), width=25, height=20, units=c("cm"),dpi=300) 

h2p=match(mbfH$HR,bioplot$IRST_NAMEB)
mbfH$biogeo=bioplot$BR[h2p]
colx=c("Alpine"="deeppink","Atlantic"="forestgreen","Boreal"="darkviolet", "Continental"="orange3","Mediterranean"="gold")
bpc=ggplot() +
  geom_vline(xintercept = 0, 
             color = "black", size = 1.5) +
  geom_hline(yintercept = 0, 
             color = "black", size = 1.5) +
  #geom_linerange(data=mbfH,aes(x=x, ymin=yq1,ymax=yq2,color = names,group=factor(names)),lwd=2,alpha=0.5) +
  #geom_linerange(data=mbfH,aes(y=y, xmin=xq1,xmax=xq2,color = names,group=factor(names)),lwd=2,alpha=0.5) +
  # geom_point(data=mbf, aes(x = x, y = y, fill = names), shape = 21, color = "black", size = 5, show.legend = FALSE) +
  geom_point(data=mbfH, aes(x = x, y = y, color=biogeo,size=val.l.x), show.legend = TRUE, stroke=1, alpha=.6) +
  #scale_fill_manual(values = colorn, name = "Drivers", labels = clabels) +
  #scale_color_manual(values = colorn, name = "Drivers", labels = clabels) +
  scale_color_manual(values = colx, name = "Regions") +
  scale_size(range = c(3, 8))+
  coord_cartesian(xlim=c(-lalim,lalim+10),ylim=c(-lalim,lalim)) +
  annotate("text", x = 40, y = 45, label = "Wetting", color = "#169dd0", size = 6, fontface = "bold") +
  annotate("text", x = -35, y = 45, label = "Accelerating", color = "#174f28", size = 6, fontface = "bold") +
  annotate("text", x = -35, y = -45, label = "Drying", color = "#dd6a29", size = 6, fontface = "bold") +
  annotate("text", x = 40, y = -45, label = "Decelerating", color = "burlywood", size = 6, fontface = "bold") +
  labs(x = "Change in drought intensity (%)", 
       y = "Change in flood intensity (%)") +
  scale_x_continuous(trans=scales::modulus_trans(.5),
                     breaks=c(-200,-100,-50,-20,-10,-5,0,5,10,20,50,100,200),
                     labels=rev(c(-200,-100,-50,-20,-10,-5,0,5,10,20,50,100,200)),
                     minor_breaks = seq(-100,100,10)) +
  scale_y_continuous(trans=scales::modulus_trans(.5),
                     breaks=c(-200,-100,-50,-10,-5,0,5,10,50,100,200),
                     minor_breaks = seq(-100,100,10), expand = c(0,0)) +
  guides( color = guide_legend(override.aes = list(size = 6)),size="none")+
  theme(
    axis.title = element_text(size = 18, face = "bold"),
    title = element_text(size = 22, face = "bold"),
    axis.text = element_text(size = 16, face = "bold"),
    panel.background = element_rect(fill = "white", colour = "white"),
    panel.grid = element_blank(),
    panel.border = element_rect(linetype = "solid", fill = NA, colour = "black"),
    legend.title = element_text(size = 20, face = "bold"),
    legend.text = element_text(size = 16),
    #legend.key = element_rect(fill = "lightgray", colour = "lightgray"),
    legend.position = "right",
    panel.grid.major = element_line(colour = "grey60"),
    panel.grid.minor= element_line(colour = "grey70", linetype = "dashed"),
    #legend.key = element_rect(fill = "white", colour = "white"),
    #legend.key.size = unit(1, "cm")
  )

bpc
ggsave(paste0("D:/tilloal/Documents/LFRuns_utils/TrendAnalysis/plots/Contribution_HRxBG_new3.jpg"),bpc, width=25, height=20, units=c("cm"),dpi=300) 

v=0

#wetting
length(which(mbfH$x>v & mbfH$y>v))/length(mbfH$HR)*100
#accelerating
length(which(mbfH$x<(-v) & mbfH$y>v))/length(mbfH$HR)*100
#decelerating
length(which(mbfH$x>v & mbfH$y<(-v)))/length(mbfH$HR)*100
#drying
length(which(mbfH$x<(-v) & mbfH$y<(-v)))/length(mbfH$HR)*100


mbfC=mbfX[which(mbfX$names=="Climate"),]

#wetting
length(which(mbfC$x>v & mbfC$y>v))/length(mbfC$HR)*100
#accelerating
length(which(mbfC$x<(-v) & mbfC$y>v))/length(mbfC$HR)*100
#decelerating
length(which(mbfC$x>v & mbfC$y<(-v)))/length(mbfC$HR)*100
#drying
length(which(mbfC$x<(-v) & mbfC$y<(-v)))/length(mbfC$HR)*100

mbfL=mbfX[which(mbfX$names=="Landuse"),]
length(which(mbfL$x>v & mbfL$y>v))/length(mbfL$HR)*100


mbfR=mbfX[which(mbfX$names=="Reservoirs"),]
length(which(mbfL$x>v & mbfL$y>v))/length(mbfL$HR)*100

mbfL$xcl=mbfL$x+mbfC$x
mbfL$ycl=mbfL$y+mbfC$y

mbfL$xcr=mbfR$x+mbfC$x
mbfL$ycr=mbfR$y+mbfC$y

mbfL$xcrl=mbfR$x+mbfC$x+mbfL$x
mbfL$ycrl=mbfR$y+mbfC$y+mbfL$y

mbfL$xrl=mbfR$x+mbfL$x
mbfL$yrl=mbfR$y+mbfL$y

abline(0,1)
abline(0,0)
abline(v=0)
mbfH$cat="wetting"
mbfH$cat[(which(mbfH$x>v & mbfH$y<(-v)))]="decelerating"
mbfH$cat[(which(mbfH$x<(-v) & mbfH$y<(-v)))]="drying"
mbfH$cat[(which(mbfH$x<(-v) & mbfH$y>v))]="accelerating"

mbfH$cat2="wetting"
mbfH$cat2[(which(mbfL$xcr>v & mbfL$ycr<(-v)))]="decelerating"
mbfH$cat2[(which(mbfL$xcr<(-v) & mbfL$ycr<(-v)))]="drying"
mbfH$cat2[(which(mbfL$xcr<(-v) & mbfL$ycr>v))]="accelerating"

mbfH$cat3="wetting"
mbfH$cat3[(which(mbfC$x>v & mbfC$y<(-v)))]="decelerating"
mbfH$cat3[(which(mbfC$x<(-v) & mbfC$y<(-v)))]="drying"
mbfH$cat3[(which(mbfC$x<(-v) & mbfC$y>v))]="accelerating"


mbfMAAX=inner_join(mbfH,mbfL,by="HR")
mbana=mbfH[which(mbfH$cat2=="decelerating"),]

length(which(mbana$cat3=="decelerating"))
length(which(mbana$cat3=="accelerating"))
length(which(mbana$cat3=="drying"))
length(which(mbana$cat3=="wetting"))

16/(24)
write.csv(mbfH,file=paste0("D:/tilloal/Documents/LFRuns_utils/TrendAnalysis/Histo_res_bvHR_new2.csv"))

write.csv(mbfX,file=paste0("D:/tilloal/Documents/LFRuns_utils/TrendAnalysis/Drivers_res_bvHR_new2.csv"))




#Standard error plot
st_geometry(databiclim) <- NULL
databic2=inner_join(mbfH,GHshpp,by=c("HR"="IRST_NAMEB"))
#databic2=full_join(databic2,mbfH,by=c("IRST_NAMEB.y"="HR"))
databic2$xdr=2*(databic2$val.sd.x/sqrt(databic2$val.l.x))*1.96
databic2$xfl=2*(databic2$val.sd.y/sqrt(databic2$val.l.x))*1.96
plot(databic2$xdr)
palet=c(hcl.colors(11, palette = "OrRd", alpha = NULL, rev = T, fixup = TRUE))

map <- ggplot(basemap) +
  geom_sf(fill="white")+
  geom_sf(data = databic2, mapping = aes(fill = xdr, geometry=geometry), alpha=1) +
  geom_sf(fill=NA, color="gray42") +
  # bi_scale_fill(pal = "BlueOr", dim = 3, na.value=colNA ) +
  # scale_fill_manual(values = colorp) +
  scale_fill_gradientn(
    colors=palet, name="95% CI of the mean", limits=c(0,5), oob = scales::squish)   +
  coord_sf(xlim = c(min(nco[,1]),max(nco[,1])), ylim = c(min(nco[,2]),max(nco[,2])))+
  #bi_scale_color(pal = "BlueOr", dim = 3, na.value=colNA ) +
  # scale_color_manual(values = colorp) +
  ggtitle("Changes in drought flows")+
  labs()+
  theme(axis.title=element_text(size=tsize),
        panel.background = element_rect(fill = "aliceblue", colour = "grey1"),
        panel.border = element_rect(linetype = "solid", fill = NA, colour="black"),
        legend.title = element_text(size=tsize),
        legend.text = element_text(size=osize),
        legend.position = "right",
        panel.grid.major = element_line(colour = "grey70"),
        panel.grid.minor = element_line(colour = "grey90"),
        legend.key = element_rect(fill = "transparent", colour = "transparent"),
        legend.key.size = unit(.8, "cm"))


map

map <- ggplot(basemap) +
  geom_sf(fill="white")+
  geom_sf(data = databic2, mapping = aes(fill = xfl,geometry=geometry ), alpha=0.7, color = "transparent", size = 0.01) +
  geom_sf(fill=NA, color="gray42") +
  # bi_scale_fill(pal = "BlueOr", dim = 3, na.value=colNA ) +
  # scale_fill_manual(values = colorp) +
  scale_fill_gradientn(
    colors=palet, name="95% CI of mean", limits=c(0,5), oob = scales::squish)   +
  coord_sf(xlim = c(min(nco[,1]),max(nco[,1])), ylim = c(min(nco[,2]),max(nco[,2])))+
  #bi_scale_color(pal = "BlueOr", dim = 3, na.value=colNA ) +
  # scale_color_manual(values = colorp) +
  labs()+
  ggtitle("Changes in flood flows")+
  theme(axis.title=element_text(size=tsize),
        panel.background = element_rect(fill = "aliceblue", colour = "grey1"),
        panel.border = element_rect(linetype = "solid", fill = NA, colour="black"),
        legend.title = element_text(size=tsize),
        legend.text = element_text(size=osize),
        legend.position = "right",
        panel.grid.major = element_line(colour = "grey70"),
        panel.grid.minor = element_line(colour = "grey90"),
        legend.key = element_rect(fill = "transparent", colour = "transparent"),
        legend.key.size = unit(.8, "cm"))


map

#\hydroregions plot
color = grDevices::colors()[grep('gr(a|e)y|white|black', grDevices::colors(), invert = T)]
coHR=sample(color, length(unique(databic2$HR)))
map <- ggplot(basemap) +
  geom_sf(fill="white")+
  geom_sf(data = databic2, mapping = aes(fill = factor(HR),geometry=geometry), alpha=0.7, color = "black", size = 0.001,show.legend = FALSE) +
  geom_sf(fill=NA, color="gray42") +
  # bi_scale_fill(pal = "BlueOr", dim = 3, na.value=colNA ) +
  # scale_fill_gradientn(
  #   colors=palet, name="95% CI of mean", limits=c(0,10), oob = scales::squish)   +
  coord_sf(xlim = c(min(nco[,1]),max(nco[,1])), ylim = c(min(nco[,2]),max(nco[,2])))+
  #bi_scale_color(pal = "BlueOr", dim = 3, na.value=colNA ) +
  scale_fill_manual(values = coHR) +
  labs()+
  theme(axis.title=element_text(size=tsize),
        panel.background = element_rect(fill = "aliceblue", colour = "grey1"),
        panel.border = element_rect(linetype = "solid", fill = NA, colour="black"),
        legend.title = element_text(size=tsize),
        legend.text = element_text(size=osize),
        legend.position = "right",
        panel.grid.major = element_line(colour = "grey70"),
        panel.grid.minor = element_line(colour = "grey90"),
        legend.key = element_rect(fill = "transparent", colour = "transparent"),
        legend.key.size = unit(.8, "cm"))


map

ggsave(paste0("D:/tilloal/Documents/LFRuns_utils/TrendAnalysis/plots/mapreg_HR.jpg"), map, width=30, height=20, units=c("cm"),dpi=300) 
# Combination of drivers ------------------------
databicr=databiclim
databicr$x= databire$x + databiclim$x
databicr$y= databire$y + databiclim$y


hist(databire$x,breaks=100,xlim=c(-50,50))
hist(databicr$y,breaks=100,xlim=c(-5,5))

breaker1=0
breaker2=10

alterclass=data.frame(databicr$x)
alterclass$class=NA
alterclass$class[which(alterclass[,1]<=(-breaker2))]=1
alterclass$class[which(alterclass[,1]>(-breaker2) & alterclass[,1]<(breaker1))]=2
alterclass$class[which(alterclass[,1]>=(breaker1) & alterclass[,1]<(breaker2))]=3
alterclass$class[which(alterclass[,1]>=breaker2)]=4
c1=alterclass$class

alterclass=data.frame(databicr$y)
alterclass$class=NA
alterclass$class[which(alterclass[,1]<=(-breaker2))]=1
alterclass$class[which(alterclass[,1]>(-breaker2) & alterclass[,1]<(breaker1))]=2
alterclass$class[which(alterclass[,1]>=(breaker1) & alterclass[,1]<(breaker2))]=3
alterclass$class[which(alterclass[,1]>=breaker2)]=4
c2=alterclass$class

cx=paste(c2,c1,sep="-")
databicr$bi_class=cx
# Combine the main category and subcategory to make a label
databicr$combined_category <-databicr$bi_class

# 
# ggplot(databicr, aes(x = d2020, y = f2020, fill = combined_category)) +
#   geom_point(shape = 21, color = "black", size = 3,show.legend = FALSE) +
#   #bi_scale_fill(pal = "BlueOr", dim = 3, na.value=colNA )+
#   scale_fill_manual(values = colors) +
#   coord_cartesian(xlim=c(-2,2),ylim=c(-50,50)) +
#   
#   # Add main quadrant labels
#   annotate("text", x = 1.8, y = 40, label = "Wetting", color = "gray12", size = 5, fontface = "bold") +
#   annotate("text", x = -1.8, y = 40, label = "Accelerating", color = "gray12", size = 5, fontface = "bold") +
#   annotate("text", x = -1.8, y = -40, label = "Drying", color = "gray12", size = 5, fontface = "bold") +
#   annotate("text", x = 1.8, y = -40, label = "Decelerating", color = "gray12", size = 5, fontface = "bold") +
#   
#   # Set axis labels
#   labs(x = "Change in drought flows (l/s/km2)", 
#        y = "Change in flood flows (l/s/km2)") +
#   
#   # Customize the theme
#   theme(axis.title=element_text(size=tsize),
#         panel.background = element_rect(fill = "transparent", colour = "grey1"),
#         panel.border = element_rect(linetype = "solid", fill = NA, colour="black"),
#         legend.title = element_text(size=tsize),
#         legend.text = element_text(size=osize),
#         legend.position = "bottom",
#         panel.grid.major = element_line(colour = "grey70"),
#         panel.grid.minor = element_line(colour = "grey90"),
#         legend.key = element_rect(fill = "transparent", colour = "transparent"),
#         legend.key.size = unit(.8, "cm"))
# 
# 


databicrl=databiclim

databicrl$x=databilu$x + databire$x + databiclim$x
databicrl$y=databilu$y + databire$y + databiclim$y


hist(databicrl$x,breaks=100,xlim=c(-5,5))
hist(databicrl$y,breaks=100,xlim=c(-5,5))

# breaker1=0
# breaker2=10

alterclass=data.frame(databicrl$x)
alterclass$class=NA
alterclass$class[which(alterclass[,1]<=(-breaker2))]=1
alterclass$class[which(alterclass[,1]>(-breaker2) & alterclass[,1]<(breaker1))]=2
alterclass$class[which(alterclass[,1]>=(breaker1) & alterclass[,1]<(breaker2))]=3
alterclass$class[which(alterclass[,1]>=breaker2)]=4
c1=alterclass$class

alterclass=data.frame(databicrl$y)
alterclass$class=NA
alterclass$class[which(alterclass[,1]<=(-breaker2))]=1
alterclass$class[which(alterclass[,1]>(-breaker2) & alterclass[,1]<(breaker1))]=2
alterclass$class[which(alterclass[,1]>=(breaker1) & alterclass[,1]<(breaker2))]=3
alterclass$class[which(alterclass[,1]>=breaker2)]=4
c2=alterclass$class
cx=paste(c2,c1,sep="-")
databicrl$bi_class=cx
# Combine the main category and subcategory to make a label
databicrl$combined_category <-databicrl$bi_class


## Pixel level -----------
databipicr=databipic
databipicr$x= databipire$x + databipic$x
databipicr$y= databipire$y + databipic$y


# breaker1=0
# breaker2=0.25

alterclass=data.frame(databipicr$x)
alterclass$class=NA
alterclass$class[which(alterclass[,1]<=(-breaker2))]=1
alterclass$class[which(alterclass[,1]>(-breaker2) & alterclass[,1]<(breaker1))]=2
alterclass$class[which(alterclass[,1]>=(breaker1) & alterclass[,1]<(breaker2))]=3
alterclass$class[which(alterclass[,1]>=breaker2)]=4
c1=alterclass$class

alterclass=data.frame(databipicr$y)
alterclass$class=NA
alterclass$class[which(alterclass[,1]<=(-breaker2))]=1
alterclass$class[which(alterclass[,1]>(-breaker2) & alterclass[,1]<(breaker1))]=2
alterclass$class[which(alterclass[,1]>=(breaker1) & alterclass[,1]<(breaker2))]=3
alterclass$class[which(alterclass[,1]>=breaker2)]=4
c2=alterclass$class

cx=paste(c2,c1,sep="-")
databipicr$bi_class=cx
# Combine the main category and subcategory to make a label
databipicr$combined_category <-databipicr$bi_class


databipicrl=databipic

databipicrl$x=databipilu$x + databipire$x + databipic$x
databipicrl$y=databipilu$y + databipire$y + databipic$y

# breaker1=0
# breaker2=0.25

alterclass=data.frame(databipicrl$x)
alterclass$class=NA
alterclass$class[which(alterclass[,1]<=(-breaker2))]=1
alterclass$class[which(alterclass[,1]>(-breaker2) & alterclass[,1]<(breaker1))]=2
alterclass$class[which(alterclass[,1]>=(breaker1) & alterclass[,1]<(breaker2))]=3
alterclass$class[which(alterclass[,1]>=breaker2)]=4
c1=alterclass$class

alterclass=data.frame(databipicrl$y)
alterclass$class=NA
alterclass$class[which(alterclass[,1]<=(-breaker2))]=1
alterclass$class[which(alterclass[,1]>(-breaker2) & alterclass[,1]<(breaker1))]=2
alterclass$class[which(alterclass[,1]>=(breaker1) & alterclass[,1]<(breaker2))]=3
alterclass$class[which(alterclass[,1]>=breaker2)]=4
c2=alterclass$class
cx=paste(c2,c1,sep="-")
databipicrl$bi_class=cx
# Combine the main category and subcategory to make a label
databipicrl$combined_category <-databipicrl$bi_class


# Some plots here --------------
# 
# 
# ggplot(databicrl, aes(x = d2020, y = f2020, fill = combined_category)) +
#   geom_point(shape = 21, color = "black", size = 3,show.legend = FALSE) +
#   #bi_scale_fill(pal = "BlueOr", dim = 3, na.value=colNA )+
#   scale_fill_manual(values = colors) +
#   coord_cartesian(xlim=c(-2,2),ylim=c(-50,50)) +
#   
#   # Add main quadrant labels
#   annotate("text", x = 1.8, y = 40, label = "Wetting", color = "gray12", size = 5, fontface = "bold") +
#   annotate("text", x = -1.8, y = 40, label = "Accelerating", color = "gray12", size = 5, fontface = "bold") +
#   annotate("text", x = -1.8, y = -40, label = "Drying", color = "gray12", size = 5, fontface = "bold") +
#   annotate("text", x = 1.8, y = -40, label = "Decelerating", color = "gray12", size = 5, fontface = "bold") +
#   
#   # Set axis labels
#   labs(x = "Change in drought flows (l/s/km2)", 
#        y = "Change in flood flows (l/s/km2)") +
#   
#   # Customize the theme
#   theme(axis.title=element_text(size=tsize),
#         panel.background = element_rect(fill = "transparent", colour = "grey1"),
#         panel.border = element_rect(linetype = "solid", fill = NA, colour="black"),
#         legend.title = element_text(size=tsize),
#         legend.text = element_text(size=osize),
#         legend.position = "bottom",
#         panel.grid.major = element_line(colour = "grey70"),
#         panel.grid.minor = element_line(colour = "grey90"),
#         legend.key = element_rect(fill = "transparent", colour = "transparent"),
#         legend.key.size = unit(.8, "cm"))
# 
# 
# 
# 
# 
# colNA="grey"
# map <- ggplot(basemap) +
#   geom_sf(fill="white")+
#   geom_sf(data = databicrl, mapping = aes(fill = combined_category), alpha=0.9, color = "transparent", size = 0.01,show.legend = F) +
#   geom_sf(fill=NA, color="gray42") +
#   coord_sf(xlim = c(min(nco[,1]),max(nco[,1])), ylim = c(min(nco[,2]),max(nco[,2])))+
#   # bi_scale_fill(pal = "BlueOr", dim = 3, na.value=colNA ) +
#   scale_fill_manual(values = colors) +
#   labs()+
#   theme(axis.title=element_text(size=tsize),
#         panel.background = element_rect(fill = "aliceblue", colour = "grey1"),
#         panel.border = element_rect(linetype = "solid", fill = NA, colour="black"),
#         legend.title = element_text(size=tsize),
#         legend.text = element_text(size=osize),
#         legend.position = "bottom",
#         panel.grid.major = element_line(colour = "grey70"),
#         panel.grid.minor = element_line(colour = "grey90"),
#         legend.key = element_rect(fill = "transparent", colour = "transparent"),
#         legend.key.size = unit(.8, "cm"))
# 
# map
# 
# 
# 
# 
# bi_pal(pal = colors, dim = 4)
# 
# legend <- bi_legend(pal = colors,
#                     dim = 4,
#                     xlab = "  +  Drought changes  -  ",
#                     ylab = "  -  Flood changes  +  ",
#                     size = 16,
#                     arrows = FALSE)
# 
# legend
# 
# # combine map with legend
# # finalPlot <- ggdraw() +
# #   draw_plot(map, 0, 0, 1, 1) +
# #   draw_plot(legend, 0.7, .65, 0.2, 0.2)
# # 
# # finalPlot
# 
# pl=ggarrange(map, legend, 
#              labels = c("Map", "Key"),
#              ncol = 2, nrow = 1,widths = c(2,1), heights=c(2,1), vjust=-1)
# pl


# SANKEY AT CATCHMENT LEVEL ------------

#Sankey diagram of transers between classes
databiclim=databiclim[-which(is.na(databiclim$d2020)),]
databicrl=databicrl[-which(is.na(databicrl$d2020)),]
databicr=databicr[-which(is.na(databicr$d2020)),]
databitot=databitot[-which(is.na(databitot$d2020)),]


databiclim$maxicat="Wetting"
databiclim$maxicat[which(databiclim$combined_category == "1-1" |
                         databiclim$combined_category == "1-2" |
                         databiclim$combined_category == "2-1" )]="Drying"

databiclim$maxicat[which(databiclim$combined_category == "3-1" |
                           databiclim$combined_category == "4-2" |
                           databiclim$combined_category == "4-1" )]="Decelerating"

databiclim$maxicat[which(databiclim$combined_category == "1-3" |
                           databiclim$combined_category == "1-4" |
                           databiclim$combined_category == "2-4" )]="Accelerating"

databiclim$maxicat[which(databiclim$combined_category == "3-3" |
                          databiclim$combined_category == "3-2" |
                          databiclim$combined_category == "2-3" |
                          databiclim$combined_category == "2-2" )]="Stable"

databicr$maxicat="Wetting"
databicr$maxicat[which(databicr$combined_category == "1-1" |
                           databicr$combined_category == "1-2" |
                           databicr$combined_category == "2-1")]="Drying"

databicr$maxicat[which(databicr$combined_category == "3-1" |
                           databicr$combined_category == "4-2" |
                           databicr$combined_category == "4-1" )]="Decelerating"

databicr$maxicat[which(databicr$combined_category == "1-3" |
                           databicr$combined_category == "1-4" |
                           databicr$combined_category == "2-4" )]="Accelerating"

databicr$maxicat[which(databicr$combined_category == "3-3" |
                           databicr$combined_category == "3-2" |
                           databicr$combined_category == "2-3" |
                           databicr$combined_category == "2-2" )]="Stable"

databicrl$maxicat="Wetting"
databicrl$maxicat[which(databicrl$combined_category == "1-1" |
                           databicrl$combined_category == "1-2" |
                           databicrl$combined_category == "2-1" )]="Drying"

databicrl$maxicat[which(databicrl$combined_category == "3-1" |
                           databicrl$combined_category == "4-2" |
                           databicrl$combined_category == "4-1" )]="Decelerating"

databicrl$maxicat[which(databicrl$combined_category == "1-3" |
                           databicrl$combined_category == "1-4" |
                           databicrl$combined_category == "2-4" )]="Accelerating"

databicrl$maxicat[which(databicrl$combined_category == "3-3" |
                           databicrl$combined_category == "3-2" |
                           databicrl$combined_category == "2-3" |
                           databicrl$combined_category == "2-2" )]="Stable"


databitot$maxicat="Wetting"
databitot$maxicat[which(databitot$combined_category == "1-1" |
                           databitot$combined_category == "1-2" |
                           databitot$combined_category == "2-1" )]="Drying"

databitot$maxicat[which(databitot$combined_category == "3-1" |
                           databitot$combined_category == "4-2" |
                           databitot$combined_category == "4-1")]="Decelerating"

databitot$maxicat[which(databitot$combined_category == "1-3" |
                           databitot$combined_category == "1-4" |
                           databitot$combined_category == "2-4" )]="Accelerating"

databitot$maxicat[which(databitot$combined_category == "3-3" |
                           databitot$combined_category == "3-2" |
                           databitot$combined_category == "2-3" |
                           databitot$combined_category == "2-2" )]="Stable"

class1=databiclim$maxicat
class2=databicr$maxicat
class3=databicrl$maxicat
class4=databitot$maxicat
# class1=databiclim$combined_category
# class2=databicrl$combined_category
ltot=length(databiclim$bi_class)
loscolors=c("Accelerating" = "#174f28","Drying" = "#dd6a29","Stable"="gray60","Wetting" = "#169dd0","Decelerating" = "burlywood")
#8 categories 
# Define the nodes
brk=unique(class1)
library(ggsankey)


d <- data.frame(cbind(class1,class2,class3,class4))
names(d) <- c('OnlyClimate',
              'ClimateReservoirs',
              'ClimateReservoirLanduse',
              'AllDrivers')

df <- d%>%
  make_long(OnlyClimate, 
            ClimateReservoirs,
            ClimateReservoirLanduse,
            AllDrivers)


df$node <- factor(df$node,levels = brk)
df$next_node <- factor(df$next_node,levels = brk)

df$idnext=c(1:length(df$next_node))
reagg <- df%>%
  dplyr::group_by(node,x)%>%  # Here we are grouping the data by node and then we are taking the frequency of it 
  tally()

df2 <- merge(df, 
             reagg, 
             by.x = 'node', 
             by.y = 'node', 
             all.x = FALSE)
df2=df2[which(df2$x.x==df2$x.y),]
df2$np=round(df2$n/ltot*100)
pl <- ggplot(df2, aes(x = x.x,                        
                      next_x = next_x,                                     
                      node = node,
                      next_node = next_node,        
                      fill = factor(node),
                      label= paste0(node,"\n", np,"%")))



pl <- pl +geom_sankey(flow.alpha = 0.5,          #This Creates the transparency of your node 
                      node.color = "black",     # This is your node color        
                      show.legend = TRUE)        # This determines if you want your legend to show

pl <- pl + geom_sankey_label(size = 3, 
                             color = "black", 
                             fill = "white")

pl 

pl <- pl + theme_sankey(base_size = 18) 
pl <- pl + theme(legend.position = 'none')
pl <- pl + theme(axis.title = element_blank(),
                 axis.text.y = element_blank(),
                 axis.ticks = element_blank(),
                 panel.grid = element_blank())

pl=pl + scale_fill_manual(values = loscolors) +
  ggtitle("Change in the Europe's hydrological cycle between 1951 and 2020")

pl


ggsave(paste0("D:/tilloal/Documents/LFRuns_utils/TrendAnalysis/plots/Sankey_changes_catNew.jpg"), pl, width=30, height=20, units=c("cm"),dpi=300) 

mean_lf=c(mean(databitot$d2020,na.rm=T),mean(databiclim$d2020,na.rm=T),mean(databicr$d2020,na.rm=T),mean(databicrl$d2020,na.rm=T))
mean_scen=data.frame()

ggplot() +
  geom_point(data=databitot, aes(x = d2020, y = f2020, fill = combined_category),
             shape = 21, color = "transparent", fill="blue",size = 3,alpha=0.1,show.legend = FALSE) +
  geom_point(data=databiclim, aes(x = d2020, y = f2020, fill = combined_category),
             shape = 21, color = "transparent", fill= "red",alpha=0.1, size = 3,show.legend = FALSE) +
  geom_point(data=databicr, aes(x = d2020, y = f2020, fill = combined_category),
             shape = 21, color = "transparent", fill= "darkgreen",alpha=0.1,size = 3,show.legend = FALSE) +
  geom_point(data=databicrl, aes(x = d2020, y = f2020, fill = combined_category),
             shape = 21, color = "transparent", fill= "purple",alpha=0.1,size = 3,show.legend = FALSE) +
  #gg_bagplot(data=databicrl, d2020, f2020, color = "#00659e", scatterplot = FALSE)
  #bi_scale_fill(pal = "BlueOr", dim = 3, na.value=colNA )+
  #scale_fill_manual(values = colors) +
  coord_cartesian(xlim=c(-100,100),ylim=c(-50,50)) +
  
  # Add main quadrant labels
  annotate("text", x = 1.8, y = 40, label = "Wetting", color = "gray12", size = 5, fontface = "bold") +
  annotate("text", x = -1.8, y = 40, label = "Accelerating", color = "gray12", size = 5, fontface = "bold") +
  annotate("text", x = -1.8, y = -40, label = "Drying", color = "gray12", size = 5, fontface = "bold") +
  annotate("text", x = 1.8, y = -40, label = "Decelerating", color = "gray12", size = 5, fontface = "bold") +
  
  # Set axis labels
  labs(x = "Change in drought flows (l/s/km2)", 
       y = "Change in flood flows (l/s/km2)") +
  
  # Customize the theme
  theme(axis.title=element_text(size=tsize),
        panel.background = element_rect(fill = "transparent", colour = "grey1"),
        panel.border = element_rect(linetype = "solid", fill = NA, colour="black"),
        legend.title = element_text(size=tsize),
        legend.text = element_text(size=osize),
        legend.position = "bottom",
        panel.grid.major = element_line(colour = "grey70"),
        panel.grid.minor = element_line(colour = "grey90"),
        legend.key = element_rect(fill = "transparent", colour = "transparent"),
        legend.key.size = unit(.8, "cm"))



#reproduce my color palette: 
loscolors=c("#174f28","#845e29","#dd6a29","#167984","#819185","#d8a386","#169dd0","#7ebbd2","#d3d3d3")
al=0.8
dataX=databitot[-which(is.na(databitot$d2020)),]
#dataX=dataX[-which(is.na(dataX$f2020)),]
library(MASS)
dataX$density <- get_density(dataX$d2020, dataX$f2020, n = 300)
monstre_legend<-ggplot() +
  
  # annotate("rect", xmin=-Inf, xmax=-0.1, ymin=-Inf, ymax=-5, fill=loscolors[9], alpha=al) +
  # annotate("rect", xmin=-Inf, xmax=-0.1, ymin=-5, ymax=5, fill=loscolors[8], alpha=al) +
  # annotate("rect", xmin=-Inf, xmax=-0.1, ymin=5, ymax=Inf, fill=loscolors[7], alpha=al) +
  # 
  # annotate("rect", xmin=-0.1, xmax=0.1, ymin=-Inf, ymax=-5, fill=loscolors[6], alpha=al) +
  # annotate("rect", xmin=-0.1, xmax=0.1, ymin=-5, ymax=5, fill=loscolors[5], alpha=al) +
  # annotate("rect", xmin=-0.1, xmax=0.1, ymin=5, ymax=Inf, fill=loscolors[4], alpha=al) +
  # 
  # annotate("rect", xmin=0.1, xmax=Inf, ymin=-Inf, ymax=-5, fill=loscolors[3], alpha=al) +
  # annotate("rect", xmin=0.1, xmax=Inf, ymin=-5, ymax=5, fill=loscolors[2], alpha=al) +
  # annotate("rect", xmin=0.1, xmax=Inf, ymin=5, ymax=Inf, fill=loscolors[1], alpha=al) +
  geom_tile(data=data, aes(x = x*b, y = y*a, fill = combined_category),alpha=0.5) +
  scale_fill_manual(values = colors) +
  coord_cartesian(xlim=c(-2.5,2.5),ylim=c(-100,100), expand = FALSE) +
  
  # Add main quadrant labels
  annotate("text", x = 4, y = 4, label = "Wetting", color = "white", size = 5, fontface = "bold") +
  annotate("text", x = -4, y = 4, label = "Accelerating", color = "white", size = 5, fontface = "bold") +
  annotate("text", x = -4, y = -4, label = "Drying", color = "white", size = 5, fontface = "bold") +
  annotate("text", x = 4, y = -4, label = "Decelerating", color = "white", size = 5, fontface = "bold") +
  
  # Set axis labels
  labs(x = "Change in drought flows (% per decade)", 
       y = "Change in flood flows (% per decade)") +
  geom_point(data=dataX, aes(x=d2020, y=f2020, col=density),alpha=0.5,size=2, shape=16) +
  #geom_point(data= Rsig, aes(x=ObsChange, y=SimChange),fill="transparent", color="gray25",shape=21,size=2)+
  #geom_abline(slope=1,intercept=0,col="gray25",lwd=1.5,alpha=.5)+
  # 
  # annotate("label", x=-90, y=100, label= paste0("N = ",ln2),size=5)+
  # annotate("label", x=90, y=100, label= paste0("N = ",lp1),size=5)+
  # annotate("label", x=90, y=-100, label= paste0("N = ",lp2),size=5)+
  # annotate("label", x=-90, y=-100, label= paste0("N = ",ln1),size=5)+
  scale_color_viridis(option="F")+
  scale_x_continuous(name="10-y Drought RL change (l/s/km2)",breaks = seq(-5,5,by=1), labels=seq(-5,5,by=1),limits = c(-2.5,2.5))+
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

monstre_legend

monstre_legend <- monstre_legend + 
  theme(aspect.ratio = 1) 

pl=ggarrange(map, 
             ggarrange(legend, monstre_legend, 
                       ncol = 1, nrow = 2),
             ncol=2, nrow=1, widths = c(3,3), heights=c(1,1),vjust=-1)

pl=ggarrange(map, monstre_legend, 
             ncol = 2, nrow = 1,widths = c(1,1), vjust=-1)

pl<-annotate_figure(pl, top = text_grob(paste0("Joint changes in floods and droughts (",period[1],"-",period[2],")"), 
                                        color = "black", face = "bold", size = 14))

ggsave(paste0("D:/tilloal/Documents/LFRuns_utils/TrendAnalysis/plots/bvHydro_totalchange.jpg"), pl, width=20, height=20, units=c("cm"),dpi=1000) 



#Continue this with improved plot to check what happens



############### PIXEL LEVEL TREND ######################


####### SANKEY for pixels #################

mean(databipi$x,na.rm=T)
mean(databipic$y,na.rm=T)
mean(databipire$y,na.rm=T)
mean(databipicrl$x,na.rm=T)
#Sankey diagram of transers between classes
databipi=databipi[-which(is.na(databipicrl$x)),]
databipic=databipic[-which(is.na(databipicrl$x)),]
databipicr=databipicr[-which(is.na(databipicrl$x)),]
databipicrl=databipicrl[-which(is.na(databipicrl$x)),]


databipic$maxicat="Wetting"
databipic$maxicat[which(databipic$bi_class == "1-1" |
                           databipic$bi_class == "1-2" |
                           databipic$bi_class == "2-1" )]="Drying"

databipic$maxicat[which(databipic$bi_class == "3-1" |
                           databipic$bi_class == "4-2" |
                           databipic$bi_class == "4-1" )]="Decelerating"

databipic$maxicat[which(databipic$bi_class == "1-3" |
                           databipic$bi_class == "1-4" |
                           databipic$bi_class == "2-4" )]="Accelerating"

databipic$maxicat[which(databipic$bi_class == "3-3" |
                           databipic$bi_class == "3-2" |
                           databipic$bi_class == "2-3" |
                           databipic$bi_class == "2-2" )]="Stable"

databipicr$maxicat="Wetting"
databipicr$maxicat[which(databipicr$bi_class == "1-1" |
                          databipicr$bi_class == "1-2" |
                          databipicr$bi_class == "2-1")]="Drying"

databipicr$maxicat[which(databipicr$bi_class == "3-1" |
                          databipicr$bi_class == "4-2" |
                          databipicr$bi_class == "4-1" )]="Decelerating"

databipicr$maxicat[which(databipicr$bi_class == "1-3" |
                          databipicr$bi_class == "1-4" |
                          databipicr$bi_class == "2-4" )]="Accelerating"

databipicr$maxicat[which(databipicr$bi_class == "3-3" |
                          databipicr$bi_class == "3-2" |
                          databipicr$bi_class == "2-3" |
                          databipicr$bi_class == "2-2" )]="Stable"

databipicrl$maxicat="Wetting"
databipicrl$maxicat[which(databipicrl$bi_class == "1-1" |
                           databipicrl$bi_class == "1-2" |
                           databipicrl$bi_class == "2-1" )]="Drying"

databipicrl$maxicat[which(databipicrl$bi_class == "3-1" |
                           databipicrl$bi_class == "4-2" |
                           databipicrl$bi_class == "4-1" )]="Decelerating"

databipicrl$maxicat[which(databipicrl$bi_class == "1-3" |
                           databipicrl$bi_class == "1-4" |
                           databipicrl$bi_class == "2-4" )]="Accelerating"

databipicrl$maxicat[which(databipicrl$bi_class == "3-3" |
                           databipicrl$bi_class == "3-2" |
                           databipicrl$bi_class == "2-3" |
                           databipicrl$bi_class == "2-2" )]="Stable"


databipi$maxicat="Wetting"
databipi$maxicat[which(databipi$bi_class == "1-1" |
                          databipi$bi_class == "1-2" |
                          databipi$bi_class == "2-1" )]="Drying"

databipi$maxicat[which(databipi$bi_class == "3-1" |
                          databipi$bi_class == "4-2" |
                          databipi$bi_class == "4-1")]="Decelerating"

databipi$maxicat[which(databipi$bi_class == "1-3" |
                          databipi$bi_class == "1-4" |
                          databipi$bi_class == "2-4" )]="Accelerating"

databipi$maxicat[which(databipi$bi_class == "3-3" |
                          databipi$bi_class == "3-2" |
                          databipi$bi_class == "2-3" |
                          databipi$bi_class == "2-2" )]="Stable"




#how many pixel are accelerating, deceleration, etch in each biogeoregion

bgname=c("Alpine","Atlantic","Boreal", "Continental","Mediterranean")
categorix=c("Accelerating","Drying","Stable","Wetting" ,"Decelerating")

#only climate driver
# Initialize an empty data frame with 'bgname' as row names and 'categorix' as column names
results_df <- data.frame(matrix(0, nrow = length(bgname), ncol = length(categorix)))
rownames(results_df) <- bgname
colnames(results_df) <- categorix

# Loop through each 'neu' in 'bgname'
for (neu in bgname) {
  databipicBG <- databipic[which(databipic$Biogeo_id == neu), ]
  
  # Loop through each 'cat' in 'categorix'
  for (cat in categorix) {
    lol <- length(which(databipicBG$maxicat == cat)) / length(databipicBG$maxicat) * 100
    
    # Fill the data frame with the percentage for the current 'neu' and 'cat'
    results_df[neu, cat] <- round(lol)
  }
}

# Print the results as a table
print(results_df)

#All drivers
# Initialize an empty data frame with 'bgname' as row names and 'categorix' as column names
results_df_all <- data.frame(matrix(0, nrow = length(bgname), ncol = length(categorix)))
rownames(results_df_all) <- bgname
colnames(results_df_all) <- categorix

# Loop through each 'neu' in 'bgname'
for (neu in bgname) {
  databipicBG <- databipi[which(databipi$Biogeo_id == neu), ]
  
  # Loop through each 'cat' in 'categorix'
  for (cat in categorix) {
    lol <- length(which(databipicBG$maxicat == cat)) / length(databipicBG$maxicat) * 100
    
    # Fill the data frame with the percentage for the current 'neu' and 'cat'
    results_df_all[neu, cat] <- round(lol)
  }
}

# Print the results as a table
print(results_df_all)


# Loop through each 'cat' in 'categorix'
results_eu=c()
for (cat in categorix) {
  lol <- length(which(databipic$maxicat == cat)) / length(databipic$maxicat) * 100
  
  # Fill the data frame with the percentage for the current 'neu' and 'cat'
  results_eu <- c(results_eu,round(lol))
}

# Print the results as a table
print(results_df)


class1=databipic$maxicat
class2=databipicr$maxicat
class3=databipicrl$maxicat
class4=databipi$maxicat

ltot=length(databipi$bi_class)
# class1=databiclim$combined_category
# class2=databicrl$combined_category

loscolors=c("Accelerating" = "#174f28","Drying" = "#dd6a29","Stable"="gray60","Wetting" = "#169dd0","Decelerating" = "burlywood")
#8 categories 
# Define the nodes
brk=unique(class1)
library(ggsankey)


d <- data.frame(cbind(class1,class2,class3,class4))
names(d) <- c('OnlyClimate',
              'ClimateReservoirs',
              'ClimateReservoirLanduse',
              'AllDrivers')

df <- d%>%
  make_long(OnlyClimate, 
            ClimateReservoirs,
            ClimateReservoirLanduse,
            AllDrivers)


df$node <- factor(df$node,levels = brk)
df$next_node <- factor(df$next_node,levels = brk)

df$idnext=c(1:length(df$next_node))
reagg <- df%>%
  dplyr::group_by(node,x)%>%  # Here we are grouping the data by node and then we are taking the frequency of it 
  tally()

df2 <- merge(df, 
             reagg, 
             by.x = 'node', 
             by.y = 'node', 
             all.x = FALSE)
df2=df2[which(df2$x.x==df2$x.y),]

df2$np=round(df2$n/ltot*100)
pl <- ggplot(df2, aes(x = x.x,                        
                      next_x = next_x,                                     
                      node = node,
                      next_node = next_node,        
                      fill = factor(node),
                      label= paste0(node,"\n", np,"%")))



pl <- pl +geom_sankey(flow.alpha = 0.5,          #This Creates the transparency of your node 
                      node.color = "black",     # This is your node color        
                      show.legend = TRUE)        # This determines if you want your legend to show

pl <- pl + geom_sankey_label(size = 3, 
                             color = "black", 
                             fill = "white")

pl 

pl <- pl + theme_sankey(base_size = 18) 
pl <- pl + theme(legend.position = 'none')
pl <- pl + theme(axis.title = element_blank(),
                 axis.text.y = element_blank(),
                 axis.ticks = element_blank(),
                 panel.grid = element_blank())

pl=pl + scale_fill_manual(values = loscolors) +
  ggtitle("Change in the Europe's hydrological cycle between 1951 and 2020")

pl


ggsave(paste0("D:/tilloal/Documents/LFRuns_utils/TrendAnalysis/plots/Sankey_changes_pixelsNew3.jpg"), pl, width=30, height=10, units=c("cm"),dpi=300) 


# part with barplots ------------------
#barplot of contribution by driver for each category
ptnnn=cbind(ClimFloodTrendPix,ResFloodTrendPix$Y2020,LuFloodTrendPix$Y2020,WuFloodTrendPix$Y2020)
ptnnn=FloodTrends
mx=(match(FloodTrendsP$outl2,ptnnn$outl2))
dt_acc=data.frame(d[mx,],ptnnn)

fagg = aggregate(list(val=dt_acc$Y2020),
                       by = list(Group=dt_acc$AllDrivers,driver=dt_acc$driver),
                       FUN = function(x) c(mean2=mean(x,na.rm=T),dev2=sd(x,na.rm=T),len2=length(x),median2=quantile(x,0.5,na.rm=T)))
fagg <- do.call(data.frame, fagg)


ptnnn=DroughtTrends
mx=(match(DroughtTrends$outl2,ptnnn$outl2))
dt_acc=data.frame(d[mx,],ptnnn)

dagg =aggregate(list(val=dt_acc$Y2020),
                by = list(Group=dt_acc$AllDrivers,driver=dt_acc$driver),
                FUN = function(x) c(mean=mean(x,na.rm=T),dev=sd(x,na.rm=T),len=length(x),median=quantile(x,0.5,na.rm=T)))
dagg <- do.call(data.frame, dagg)

fdagg=cbind(fagg,dagg)
grp=unique(fagg$Group)
sel_agg=c()
for (f in 1:5){
  sel=fdagg[which(fagg$Group==grp[f]),]
  sel=sel[-which(sel$driver=="Total"),]
  tdchange=sum(abs(sel$val.mean))
  tfchange=sum(abs(sel$val.mean2))
  sel$rvd=(abs(sel$val.mean))/tdchange*100
  sel$rvf=(abs(sel$val.mean2))/tfchange*100
  sel_agg=rbind(sel_agg,sel)
}

sel_agg$signlow="positive"
sel_agg$signlow[which(sel_agg$val.mean<0)]="negative"

sel_agg$signhi="positive"
sel_agg$signhi[which(sel_agg$val.mean2<0)]="negative"
sel_agg=sel_agg[,-c(1,2)]

colorb=rev(c("royalblue","red2"))
colorz = c("Clim" ='dodgerblue4',"Landuse" ='gold4',"Reservoirs" ='firebrick4',"Wateruse"="olivedrab")
colorn = c("wateruse" ='limegreen',"Reservoirs" ='tomato4',"Landuse" ='orange',"Clim" ='royalblue')

saveplot=list()
for (idg in 1:5){
  sel_aggp=sel_agg[which(sel_agg$Group==grp[idg]),]
  f_sel=round(fdagg$val.mean2[which(fdagg$Group==grp[idg] & fdagg$driver=="Total")],2)
  if (f_sel>0) f_sel=paste0("+",f_sel)
  d_sel=round(fdagg$val.mean[which(fdagg$Group==grp[idg] & fdagg$driver=="Total")],2)
  if (d_sel>0) d_sel=paste0("+",d_sel)
  # Customize the barplot
  saveplot[[idg]]<-ggplot() +
    geom_bar(data=sel_aggp, aes(x = driver, y = rvf, fill = factor(signhi)),
             position="dodge", stat="identity",color="black") +
    geom_bar(data=sel_aggp, aes(x = driver, y = -rvd, fill = factor(signlow)),
             position="dodge", stat="identity",color="black")+
    scale_fill_manual(values = colorb,name="change \ndirection",labels=c("-","+")) +
    annotate("text", x = 3.6, y = 80, label = paste0(f_sel," l/s/km2"), color = "black", size = 5, fontface = "bold") +
    annotate("text", x = 3.6, y = -80, label = paste0(d_sel," l/s/km2"), color = "black", size = 5, fontface = "bold") +
    scale_y_continuous(
      expand = c(0,0),
      name = "share of change contribution (%)",
      breaks=seq(-100,100,100),
      labels=c(100,0,100),
      minor_breaks = seq(-2000,2000,10),
      limits=c(-100,100),
      sec.axis = sec_axis( transform=~.*1, name="",
                           breaks=seq(-20,80,100), labels=c("Droughts","Floods"))
    )+
    # geom_point(data=agbar, aes(x=drivers, y = fac*(mean), color = factor(chsig),group=factor(chsig)),
    #            position=position_dodge(width=.9),size=3) +
    # guides(color = guide_legend(override.aes = list(color = colorb)))+
    # scale_x_discrete(labels=xlabels)+
    # geom_linerange(data=agbar,aes(x=drivers, ymin=fac*(q1),ymax=fac*(q2),color = factor(chsig),group=factor(chsig)),
    #                position = position_dodge2(width = .9),lwd=1) +
    scale_color_manual(values = colorb, name="Change group",labels=ncol,guide = "none") +
    theme(axis.title=element_text(size=16, face="bold"),
          axis.text.y = element_text( face="bold",size=12),
          axis.text.y.right =element_text(angle=270, face="bold",size=16, color="black"),
          axis.text.x = element_text(size=12,face="bold"),
          panel.background = element_rect(fill = "white", colour = "white"),
          panel.grid = element_blank(),
          panel.border = element_rect(linetype = "solid", fill = NA, colour="black"),
          legend.title = element_text(size=12, face="bold"),
          legend.text = element_text(size=12),
          axis.ticks.y = element_blank(),  # Remove x-axis tickmarks
          panel.grid.major.y =element_line(color = "grey20"),
          panel.grid.major.x =element_line(color = "transparent"),
          legend.position = "right",
          title = element_text(size=16, face="bold"),
          panel.grid.major = element_line(colour = "grey0"),
          panel.grid.minor.y = element_line(colour = "grey85",linetype="dashed"),
          legend.key = element_rect(fill = "transparent", colour = "transparent"),
          legend.key.size = unit(.8, "cm"))+
    ggtitle(grp[idg])
  
  
  
}

saveplot[[2]]

plf=ggarrange(saveplot[[1]], saveplot[[2]], saveplot[[3]],
              saveplot[[4]],saveplot[[5]],
             ncol = 2, nrow = 3,widths = c(1,1), vjust=-1)

plf

ggsave(paste0("D:/tilloal/Documents/LFRuns_utils/TrendAnalysis/plots/change_by_class.svg"), plf, width=30, height=40, units=c("cm"),dpi=300) 


#redo this plot by Biogeoregions

mean(FloodTrendsP$Y2020[which(driver==driver[3])])

biogeo <- read_sf(dsn = paste0(hydroDir,"/eea_3035_biogeo-regions_2016/BiogeoRegions2016_wag84.shp"))
biogeof=fortify(biogeo)
st_geometry(biogeof)<-NULL
biogeoregions=raster( paste0(hydroDir,"/eea_3035_biogeo-regions_2016/Biogeo_rasterized_wsg84.tif"))
Gbiogeoregions=as.data.frame(biogeoregions,xy=T)
biogeomatch=inner_join(biogeof,Gbiogeoregions,by= c("PK_UID"="Biogeo_rasterized_wsg84"))
biogeomatch$latlong=paste(round(biogeomatch$x,4),round(biogeomatch$y,4),sep=" ")
biogeo_rivers=right_join(biogeomatch,outf, by="latlong")


#matching biogeoregions and hybas07
bio_acc=inner_join(biogeo_rivers,dt_acc,by=c("outl2"))


#barplot of contribution by driver for each category
ptnnn=cbind(ClimFloodTrendPix,ResFloodTrendPix$Y2020,LuFloodTrendPix$Y2020,WuFloodTrendPix$Y2020)
ptnnn=FloodTrendsP
mx=(match(FloodTrendsP$outl2,ptnnn$outl2))
dt_acc=data.frame(d[mx,],ptnnn)
bio_acc=inner_join(biogeo_rivers,dt_acc,by=c("outl2"))
# bio_acc$Y2020[which(bio_acc$Y2020==0)]=NA
fagg =aggregate(list(val=bio_acc$Y2020),
                by = list(Group=bio_acc$code,driver=bio_acc$driver),
                FUN = function(x) c(mean2=mean(x,na.rm=T),dev2=sd(x,na.rm=T),len2=length(x),median2=quantile(x,0.5,na.rm=T)))
fagg <- do.call(data.frame, fagg)


ptnnn=DroughtTrendsP
mx=(match(DroughtTrendsP$outl2,ptnnn$outl2))
dt_acc=data.frame(d[mx,],ptnnn)
bio_acc=inner_join(biogeo_rivers,dt_acc,by=c("outl2"))
#put a limit to max change: 1000
bio_acc$Y2020[which(bio_acc$Y2020>200)]=200
#or 
#bio_acc$Y2020[which(bio_acc$Y2020==0)]=NA
#bio_acc$Y2020[which(bio_acc$Y2020<(-200))]=200
dagg =aggregate(list(val=bio_acc$Y2020),
                by = list(Group=bio_acc$code,driver=bio_acc$driver),
                FUN = function(x) c(mean=mean(x,na.rm=T),dev=sd(x,na.rm=T),len=length(x),median=quantile(x,0.5,na.rm=T)))
dagg <- do.call(data.frame, dagg)

fdagg=cbind(fagg,dagg)
grp=unique(fagg$Group)[-c(4,8,9)]
sel_agg=c()
tfsav=c()
tdsav=c()
for (f in 1:6){
  sel=fdagg[which(fagg$Group==grp[f]),]
  sel=sel[-which(sel$driver=="Total"),]
  tdchange=sum(abs(sel$val.mean))
  tfchange=sum(abs(sel$val.mean2))
  sel$rvd=(abs(sel$val.mean))/tdchange*100
  sel$rvf=(abs(sel$val.mean2))/tfchange*100
  # sel$rvd=(abs(sel$val.mean))
  # sel$rvf=(abs(sel$val.mean2))
  sel_agg=rbind(sel_agg,sel)
  tfsav=c(tfsav,tfchange)
  tdsav=c(tdsav,tdchange)
}

sel_agg$signlow=sel_agg$val.mean/mean(tdsav)*100
#sel_agg$signlow[which(sel_agg$val.mean<0)]="negative"
sel_agg$signhi=sel_agg$val.mean2/mean(tfsav)*100
#sel_agg$signhi[which(sel_agg$val.mean2<0)]="negative"
sel_agg=sel_agg[,-c(1,2)]

colorb=rev(c("royalblue","red2"))
colorz = c("Clim" ='dodgerblue4',"Landuse" ='gold4',"Reservoirs" ='firebrick4',"Wateruse"="olivedrab")
colorn = c("wateruse" ='limegreen',"Reservoirs" ='tomato4',"Landuse" ='orange',"Clim" ='royalblue')
palet=c(hcl.colors(11, palette = "RdYlBu", alpha = NULL, rev = F, fixup = TRUE))
paletf=c(hcl.colors(11, palette = "RdBu", alpha = NULL, rev = F, fixup = TRUE))
saveplot=list()
for (idg in 1:6){
  #idg=6
  sel_aggp=sel_agg[which(sel_agg$Group==grp[idg]),]
  f_sel=round(fdagg$val.mean2[which(fdagg$Group==grp[idg] & fdagg$driver=="Total")],1)
  if (f_sel>0) f_sel=paste0("+",f_sel)
  d_sel=round(fdagg$val.mean[which(fdagg$Group==grp[idg] & fdagg$driver=="Total")],1)
  if (d_sel>0) d_sel=paste0("+",d_sel)
  # Customize the barplot
  saveplot[[idg]]<-ggplot() +
    geom_bar(data=sel_aggp, aes(x = driver, y = rvf, fill = signhi),
             position="dodge", stat="identity",color="black") +
    geom_bar(data=sel_aggp, aes(x = driver, y = -rvd, fill = signlow),
             position="dodge", stat="identity",color="black")+
    # scale_fill_manual(values = colorb,name="change \ndirection",labels=c("-","+")) +
    scale_fill_gradientn(
      colors=paletf,
      breaks=c(-100,100),labels=c(" - "," + "),limits=c(-100,100),trans=scales::modulus_trans(.2),
      oob = scales::squish,na.value=colNA, name="Change \nMagnitude") +
    annotate("text", x = 3.6, y = 80, label = paste0(f_sel," %"), color = "black", size = 5, fontface = "bold") +
    annotate("text", x = 3.6, y = -80, label = paste0(d_sel," %"), color = "black", size = 5, fontface = "bold") +
    scale_y_continuous(
      expand = c(0,0),
      name = "change contribution (%)",
      breaks=seq(-100,100,100),
      labels=c(100,0,100),
      minor_breaks = seq(-2000,2000,10),
      limits=c(-100,100),
      sec.axis = sec_axis( transform=~.*1, name="",
                           breaks=seq(-20,80,100), labels=c("Droughts","Floods"))
    )+
    guides(fill = guide_colourbar(barwidth = 1, barheight = 6))+
    # geom_point(data=agbar, aes(x=drivers, y = fac*(mean), color = factor(chsig),group=factor(chsig)),
    #            position=position_dodge(width=.9),size=3) +
    # guides(color = guide_legend(override.aes = list(color = colorb)))+
    # scale_x_discrete(labels=xlabels)+
    # geom_linerange(data=agbar,aes(x=drivers, ymin=fac*(q1),ymax=fac*(q2),color = factor(chsig),group=factor(chsig)),
    #                position = position_dodge2(width = .9),lwd=1) +
    scale_color_manual(values = colorb, name="Change group",labels=ncol,guide = "none") +
    theme(axis.title=element_text(size=16, face="bold"),
          axis.text.y = element_text( face="bold",size=12),
          axis.text.y.right =element_text(angle=270, face="bold",size=16, color="black"),
          axis.text.x = element_text(size=12,face="bold"),
          panel.background = element_rect(fill = "white", colour = "white"),
          panel.grid = element_blank(),
          panel.border = element_rect(linetype = "solid", fill = NA, colour="black"),
          legend.title = element_text(size=12, face="bold"),
          legend.text = element_text(size=12),
          axis.ticks.y = element_blank(),  # Remove x-axis tickmarks
          panel.grid.major.y =element_line(color = "grey20"),
          panel.grid.major.x =element_line(color = "transparent"),
          legend.position = "right",
          title = element_text(size=16, face="bold"),
          panel.grid.major = element_line(colour = "grey0"),
          panel.grid.minor.y = element_line(colour = "grey50",linetype="dashed"),
          legend.key = element_rect(fill = "transparent", colour = "transparent"),
          legend.key.size = unit(.8, "cm"))+
    ggtitle(grp[idg])
  
  
  
}

saveplot[[6]]

plf=ggarrange(saveplot[[1]], saveplot[[2]], saveplot[[3]],
              saveplot[[4]],saveplot[[5]],saveplot[[6]],
              ncol = 2, nrow = 3,widths = c(1,1), vjust=-1)

plf

ggsave(paste0("D:/tilloal/Documents/LFRuns_utils/TrendAnalysis/plots/change_by_biogeoregion_200.jpg"), plf, width=30, height=40, units=c("cm"),dpi=300) 





HRM=na.omit((match(RegioRLi$HydroR,HydroRsf$Id)))
HydroRsf_dom=HydroRsf[HRM,]

bhp_m=na.omit(unique(match(HydroRsf_dom$Id,bio_hybas$Id)))
HydroRsf_dom$biogeoreg=bio_hybas$code[bhp_m]



# DELIRIUMS

r=1
h=0.8
l=0.2

circle_centers <- list(
  Wetting = c(r, r),
  Accelerating = c(-r, r),
  Drying = c(-r, -r),
  Decelerating = c(r, -r)
)



# Create data frame for the grid
data <- expand.grid(x = seq(-8, 8, by = .03), y = seq(-8, 8, by = .03))

# Assign categories based on quadrants
data$category <- with(data,
                      ifelse(x >= 0 & y >= 0, "Wetting",
                             ifelse(x >= 0 & y < 0, "Decelerating",
                                    ifelse(x < 0 & y < 0, "Drying", "Accelerating"))))


# Calculate the radius (distance) for each category's quarter-circle center
data$radius <- with(data, ifelse(category == "Wetting", sqrt((x - circle_centers$Wetting[1])^2 + (y - circle_centers$Wetting[2])^2),
                                 ifelse(category == "Accelerating", sqrt((x - circle_centers$Accelerating[1])^2 + (y - circle_centers$Accelerating[2])^2),
                                        ifelse(category == "Drying", sqrt((x - circle_centers$Drying[1])^2 + (y - circle_centers$Drying[2])^2),
                                               sqrt((x - circle_centers$Decelerating[1])^2 + (y - circle_centers$Decelerating[2])^2)))))

# Assign subcategories based on:
# - Low: Within the quarter circle (radius <= 3)
# - Medium: Within 3 < radius <= 6, in the straight-line continuation
# - High: Beyond radius 6, in the straight-line continuation


data$subcategory <- with(data, 
                         ifelse(radius > h, "Low",
                                ifelse(radius > l, "Medium", "High")))



data$subcategory <- with(data, 
                         # Drying category
                         ifelse(category == "Drying" & x < -r & y <= -(h), "High",
                                
                                ifelse(category == "Drying" & x <= -(h) &  y < -r, "High",
                                       ifelse(category == "Drying" & x < -r &  y < -l, "Medium",
                                              ifelse(category == "Drying" & x < -l &  y < -r, "Medium",
                                                     
                                                     # Wetting category
                                                     ifelse(category == "Wetting" & x > r & y >= h, "High",
                                                            ifelse(category == "Wetting" & x >= h &  y > r, "High",
                                                                   ifelse(category == "Wetting" & x > r &  y > l, "Medium",
                                                                          ifelse(category == "Wetting" & x > l &  y > r, "Medium",
                                                                                 
                                                                                 # Decelerating category
                                                                                 ifelse(category == "Decelerating" & x > r & y <= -(h), "High",
                                                                                        ifelse(category == "Decelerating" & x >= h & y < -r, "High",
                                                                                               ifelse(category == "Decelerating" & x > r & y < -l, "Medium",
                                                                                                      ifelse(category == "Decelerating" & x > l & y < -r, "Medium",
                                                                                                             
                                                                                                             # Accelerating category
                                                                                                             ifelse(category == "Accelerating" & x < -r & y >= h, "High",
                                                                                                                    ifelse(category == "Accelerating" & x <= -(h) & y > r, "High",
                                                                                                                           ifelse(category == "Accelerating" & x < -r & y > l, "Medium",
                                                                                                                                  ifelse(category == "Accelerating" & x < -l & y > r, "Medium",
                                                                                                                                         # If no condition is met, keep the current subcategory value
                                                                                                                                         subcategory)))))))))))))))))




# Combine the main category and subcategory to make a label
data$combined_category <- paste(data$subcategory, data$category)




#Same shit at pixel level





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




