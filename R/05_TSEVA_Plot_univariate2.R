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



### European Biogeo regions ----

biogeo <- read_sf(dsn = paste0(hydroDir,"/eea_3035_biogeo-regions_2016/BiogeoRegions2016_wag84.shp"))
biogeof=fortify(biogeo)
st_geometry(biogeof)<-NULL
biogeoregions=raster( paste0(hydroDir,"/eea_3035_biogeo-regions_2016/Biogeo_rasterized_wsg84.tif"))
Gbiogeoregions=as.data.frame(biogeoregions,xy=T)
biogeomatch=inner_join(biogeof,Gbiogeoregions,by= c("PK_UID"="Biogeo_rasterized_wsg84"))
biogeomatch$latlong=paste(round(biogeomatch$x,4),round(biogeomatch$y,4),sep=" ")
biogeo_rivers=right_join(biogeomatch,outf, by="latlong")


### HydroRegions ----

# GridHR=raster( paste0(hydroDir,"/HydroRegions_raster_WGS84.tif"))
# GHR=as.data.frame(GridHR,xy=T)
# GHR=GHR[which(!is.na(GHR[,3])),]
# GHR$llcoord=paste(round(GHR$x,4),round(GHR$y,4),sep=" ") 
# GHR_riv=inner_join(GHR,outf,by= c("llcoord"="latlong"))
# GHshpp <- read_sf(dsn ="Z:/ClimateRun4/nahaUsers/tilloal/HydroRegions/her_all_adjusted.shp")
# HydroRsf=fortify(GHshpp) 


### NUTS3 ----


# NUTS3 <- read_sf(dsn = paste0(hydroDir,"/Countries/NUTS3/NUTS3_Extended_domain.shp"))
# NUTS3$N3ID=c(1:length(NUTS3$NUTS_ID))
# N2ID=unique(NUTS3$NUTS2_ID)
# N2IDn=c(1:length(N2ID))
# mati=match(NUTS3$NUTS2_ID,N2ID)
# NUTS3$N2ID=N2IDn[mati]
# st_write(NUTS3, paste0(hydroDir,"/Countries/NUTS3/NUTS3_modified.shp"), driver = "ESRI Shapefile")

# NUTS3 <- read_sf(dsn = paste0(hydroDir,"/Countries/NUTS3/NUTS3_modified.shp"))
# GridNUTS3=raster( paste0(hydroDir,"/Countries/NUTS3/NUTS3_Raster3ID.tif"))
# GN3=as.data.frame(GridNUTS3,xy=T)
# GN3=GN3[which(!is.na(GN3[,3])),]
# GN3$llcoord=paste(round(GN3$x,4),round(GN3$y,4),sep=" ") 
# GN3_riv=right_join(GN3,outf,by= c("llcoord"="latlong"))
# 
# GridNUTS2=raster( paste0(hydroDir,"/Countries/NUTS3/NUTS3_Raster2ID.tif"))
# GN2=as.data.frame(GridNUTS2,xy=T)
# GN2=GN2[which(!is.na(GN2[,3])),]
# GN2$llcoord=paste(round(GN2$x,4),round(GN2$y,4),sep=" ") 
# GN2_riv=right_join(GN2,outf,by= c("llcoord"="latlong"))
# 
# GNF=right_join(GN3,GN2_riv,by="llcoord")
# 
# GNUTS3sf=fortify(NUTS3) 
# 
# GNFx=GNF[which(is.na(GNF$NUTS3_Raster3ID)),]
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
# haz="Drought"
# namefile="Drought.nonfrost.Histo3"

haz="Flood"
namefile="flood.Histo4"

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
namefile="flood.socCF4"
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
namefile="flood.RWCF4"
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
namefile="flood.RCF4"
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
  
  pid=plot_inputs$pid
  pwid=plot_inputs$pwid
  prid=plot_inputs$prid
  psid=plot_inputs$psid
  
  ctrend=c(Climtrend)
  lctrend=c(Climtrend+Soctrend)
  wlctrend=c(lctrend+Wutrend)
  rwlctrend=c(wlctrend+Restrend)
  
  lutrend=c(lctrend,rev(ctrend))
  wutrend=c(wlctrend,rev(lctrend))
  retrend=c(rwlctrend,rev(wlctrend))
  pdates=c(dates,rev(dates))
  
  plot(dates,RLGPDH)
  qlim=c(.8*min(c(pid$value,prid$value,psid$value)), 1.2*max(c(pid$value,prid$value,psid$value)))
  
  plot(pid$time2, pid$value, col=alpha("grey",.8) ,pch=16,axes=FALSE,xaxs="i",yaxs="i",ylim=qlim,
       xlab = NA, ylab="")
  
  points(psid$time2,psid$value,col=alpha("royalblue",.8),type="p",pch=1)
  points(prid$time2,prid$value,col=alpha("darkgreen",.8),type="p",pch=3)
  points(pwid$time2,pwid$value,col=alpha("orange",.8),type="p",pch=4)
  
  lines(dates,RLGPDH,col="lightgrey",lwd=2,lty=1)
  lines(dates,RLGPDS,col="royalblue",lwd=2,lty=2)
  lines(dates,RLGPDR,col="darkgreen",lwd=2,lty=3)
  lines(dates,RLGPDW,col="orange",lwd=2,lty=4)
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
  legend("topleft", leg=c("Historical 10y RL","Res+WU static 10y RL","Res static 10y RL","Socio static 10y RL"),
         lwd=c(2,2,2,2), col=c("grey","orange","darkgreen","royalblue"),
         cex=1, lty=c(1,2,3,4), bg=alpha("white",.6))
} 


#who is the Arno here
pix=4103620
pix=5300097
pix=5302991

ylplot=seq(as.POSIXct("1950-06-01"),as.POSIXct("2020-06-01"),by="year")
Peak1=PeakSCF[which(PeakSCF$catch==pix),]

Peak2=PeakRWCF[which(PeakRWCF$catch==pix),]

Peak3=PeakRCF[which(PeakRCF$catch==pix),]

Peak4=PeakH[which(PeakH$catch==pix),]

if (haz=="Drought"){
  
  Peak1$time2=timeStamps[Peak1$timeID]
  Peak2$time2=timeStamps[Peak2$timeID]
  Peak3$time2=timeStamps[Peak3$timeID]
  Peak4$time2=timeStamps[Peak4$timeID]
  
  Peak1$value=-Peak1$value
  Peak2$value=-Peak2$value
  Peak3$value=-Peak3$value
  Peak4$value=-Peak4$value
  
  RL1=as.numeric(-RLGPDflSCF[which(RLGPDflSCF$unikout==pix),-71])
  RL2=as.numeric(-RLGPDflRWCF[which(RLGPDflRWCF$unikout==pix),-71])
  RL3=as.numeric(-RLGPDflRCF[which(RLGPDflRCF$unikout==pix),-71])
  RL4=as.numeric(-RLGPDflH[which(RLGPDflH$unikout==pix),-71])
}else{
  
  Peak1$time2=as.Date(Peak1$time)
  Peak2$time2=as.Date(Peak2$time)
  Peak3$time2=as.Date(Peak3$time)
  Peak4$time2=as.Date(Peak4$time)
  
  RL1=as.numeric(RLGPDflSCF[which(RLGPDflSCF$unikout==pix),-71])
  RL2=as.numeric(RLGPDflRWCF[which(RLGPDflRWCF$unikout==pix),-71])
  RL3=as.numeric(RLGPDflRCF[which(RLGPDflRCF$unikout==pix),-71])
  RL4=as.numeric(RLGPDflH[which(RLGPDflH$unikout==pix),-71])
}

plot(Peak1$time2,Peak1$value,pch=16,col="blue")
points(Peak2$time2,Peak2$value,pch=16,col="orange")
points(Peak3$time2,Peak3$value,pch=16,col="darkgreen")
points(Peak4$time2,Peak4$value,pch=16, col=1)


RL1[which(RL1<0)]=0
RL2[which(RL2<0)]=0
RL3[which(RL3<0)]=0
RL4[which(RL4<0)]=0

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
jpeg(paste("D:/tilloal/Documents/LFRuns_utils/TrendAnalysis/plots/",pix, "_drought_scenarii.jpg", sep = ""),
     width = width_px, height = height_px, quality = 100, res=dpi)

# create a function from this  
main="Reno @ Bologna"
main="Arno @ Firenze"



plot_inputs=list(Climtrend=Climtrend,Soctrend=Soctrend,Restrend=Restrend,Wutrend=Wutrend,
                 RLGPDH=RL4,RLGPDR=RL3,RLGPDS=RL1,RLGPDW=RL2,
                 pid=Peak4,psid=Peak1,prid=Peak3,pwid=Peak2)

test=plot_riverchange(plot_inputs,main,dates)

dev.off()
####################################


RLGPDflSCF=as.data.frame(RLGPDflSCF)
RLGPDflH=as.data.frame(RLGPDflH)
RLGPDflRCF=as.data.frame(RLGPDflRCF)
RLGPDflRWCF=as.data.frame(RLGPDflRWCF)

### load IRES status

if (haz=="Drought"){
  hazard="Drought"
  season="nonfrost"
  mmx=""
  
  load(file=paste0(hydroDir,"/",haz,"/IRES.",season,".Histo",mmx,".Rdata"))
  IRES_Histo=IRES_save
  load(file=paste0(hydroDir,"/",hazard,"/IRES.",season,".ResCF",mmx,".Rdata"))
  IRES_ResCF=IRES_save
  load(file=paste0(hydroDir,"/",hazard,"/IRES.",season,".RWCF",mmx,".Rdata"))
  IRES_RWCF=IRES_save
  load(file=paste0(hydroDir,"/",hazard,"/IRES.",season,".SocCF",mmx,".Rdata"))
  IRES_SocCF=IRES_save
  
  
  #who is missing in SocCF
  
  msoc=match(IRES_Histo$catlist,IRES_SocCF$catlist)
  IRES_Histo$catlist[which(is.na(msoc))]
  
  
  IRES_comb=full_join(IRES_Histo,IRES_SocCF,by="catlist")
  IRES_comb=full_join(IRES_comb,IRES_ResCF,by="catlist")
  IRES_comb=full_join(IRES_comb,IRES_RWCF,by="catlist")
  IRES_comb$IRES.y[which(is.na(IRES_comb$IRES.y))]=0
  IRES_comb$IRES.x.x[which(is.na(IRES_comb$IRES.x.x))]=0
  IRES_comb$IRES.y.y[which(is.na(IRES_comb$IRES.y.y))]=0
  IRES_comb$gen_IR=ceiling((IRES_comb$IRES.x+IRES_comb$IRES.y+IRES_comb$IRES.x.x+IRES_comb$IRES.y.y)/4)
  IRES_comb$dtectblss=(IRES_comb$IRES.x+IRES_comb$IRES.y+IRES_comb$IRES.x.x+IRES_comb$IRES.y.y)
  names(IRES_comb)[c(2,3,4,5)]=c("Histo","SCF","RCF","RWCF")
  IR_locs=which(IRES_comb$gen_IR==1)
  length(which(IRES_comb$gen_IR==1))


}

#Identify shitty locations and find parameters
for (run in 1:4){
  print(run)
  if (run==1) RLcheck=RLGPDflSCF
  if (run==2) RLcheck=RLGPDflRWCF
  if (run==3) RLcheck=RLGPDflRCF
  if (run==4) RLcheck=RLGPDflH
  
  pbshit=RLcheck$unikout[which(is.nan((RLcheck$Y2020)))]
  if (length(pbshit)>1){
    for (f in 1:length(pbshit)){
      fi=pbshit[f]
      print(fi)
      c=RLcheck[which(RLcheck$unikout==fi),]
      last=which(is.nan(as.numeric(c)))[1]-1
      c[which(is.nan(as.numeric(c)))]=c[last]
      RLcheck[which(RLcheck$unikout==fi),]=c
    }
    if (run==1) RLGPDflSCF=RLcheck
    if (run==2) RLGPDflRWCF=RLcheck
    if (run==3) RLGPDflRCF=RLcheck
    if (run==4) RLGPDflH=RLcheck
  }
}


### For drought ----------------

if (haz=="Drought"){
#Histo
  RLIRH<- RLGPDflH[IR_locs,]
  RLGPDflH[IR_locs,c(1:70)]<-NA
  
  RLGPDflH[,c(1:70)]=-RLGPDflH[,c(1:70)]
  for (id in c(1:70)){
    RLGPDflH[which(RLGPDflH[,id]<0),id]=0
    RLGPDflH[which(is.infinite(RLGPDflH[,id])),id]=NA
  }
  
  RLIRRCF<- RLGPDflRCF[IR_locs,]
  RLGPDflRCF[IR_locs,c(1:70)]<-NA
  
  RLGPDflRCF[,c(1:70)]=-RLGPDflRCF[,c(1:70)]
  for (id in c(1:70)){
    RLGPDflRCF[which(RLGPDflRCF[,id]<0),id]=0
    RLGPDflRCF[which(is.infinite(RLGPDflRCF[,id])),id]=NA
  }
  
  RLIRSCF<- RLGPDflSCF[IR_locs,]
  RLGPDflSCF[IR_locs,c(1:70)]<-NA
  
  RLGPDflSCF[,c(1:70)]=-RLGPDflSCF[,c(1:70)]
  for (id in c(1:70)){
    RLGPDflSCF[which(RLGPDflSCF[,id]<0),id]=0
    RLGPDflSCF[which(is.infinite(RLGPDflSCF[,id])),id]=NA
  }
  
  
  RLIRRWCF<- RLGPDflRWCF[IR_locs,]
  RLGPDflRWCF[IR_locs,c(1:70)]<-NA
  
  RLGPDflRWCF[,c(1:70)]=-RLGPDflRWCF[,c(1:70)]
  for (id in c(1:70)){
    RLGPDflRWCF[which(RLGPDflRWCF[,id]<0),id]=0
    RLGPDflRWCF[which(is.infinite(RLGPDflRWCF[,id])),id]=NA
  }
  
  
  
  
  IRpoints=inner_join(IRES_comb,UpArea,by=c("catlist"="outl2"))
  
  #Plot where IRs are
  
  colIR=c("0"="royalblue","1"="orange","2"="orangered","3"="tomato","4"="purple")
  points <- st_as_sf(IRpoints, coords = c("Var1.x", "Var2.x"), crs = 4326)
  points <- st_transform(points, crs = 3035)
  
IRmap=ggplot(basemap) +
    geom_sf(fill="gray95",color="gray10",size=0.5)+
    geom_sf(data=points,aes(col=factor(dtectblss),geometry=geometry,size=upa),alpha=.9,stroke=0,shape=15)+
    coord_sf(xlim = c(min(nco[,1]),max(nco[,1])), ylim = c(min(nco[,2]),max(nco[,2])))+
    scale_colour_manual(values = colIR, name="IR", labels=c("0"="perennial","1"="IRES")) +
    scale_size(range = c(0.08, 0.4), trans="sqrt",name= expression(paste("Upstream area ", (km^2),
                                                                         sep = " ")),
               breaks=c(101,1000,10000,100000,500000), labels=c("100","1000", "10 000", "100 000", "500 000"),
               guide = "none")+
    #ggtitle(title)+
    labs(x="Longitude", y = "Latitude")+
    guides(colour = guide_legend(override.aes = list(size = 10)))+
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
  
  ggsave(paste0("D:/tilloal/Documents/LFRuns_utils/TrendAnalysis/plots/IRES_allruns3.jpg"), IRmap, width=20, height=20, units=c("cm"),dpi=1000) 



}





## Change attribution-----------------------
#Before change attribution, I compute total change to compare with final results

#mato=match(RLGPDflSCF[,71],RLGPDflH[,71])


### change from climate--------------------
Climtrend=RLGPDflSCF[,c(1:70)]-(RLGPDflSCF[,1])

#1 I remove the climate trend
RLGPDflSt=RLGPDflSCF[,c(1:70)]-Climtrend
RLGPDflRWt=RLGPDflRWCF[,c(1:70)]-Climtrend
RLGPDflRt=RLGPDflRCF[,c(1:70)]-Climtrend
RLGPDflHt=RLGPDflH[,c(1:70)]-Climtrend

### change from socoeconomy -----------



## socioeconomic trend (land use) -----------
RLGPDflRWtp=RLGPDflRWt-(RLGPDflRWt-RLGPDflSt)
RLGPDflRtp=RLGPDflRt-(RLGPDflRWt-RLGPDflSt)
RLGPDflHtp=RLGPDflHt-(RLGPDflRWt-RLGPDflSt)

Soctrend=(RLGPDflRWt-RLGPDflSt)


#try change between two runs without cchange removal
Soctrend=(RLGPDflRWCF[,-71]-RLGPDflSCF[,-71])


plot(as.numeric(RLGPDflRWCF[6786,-71]))
lines(as.numeric(RLGPDflSCF[6786,-71]))
RLGPDflSt[6786,]

#plot(as.numeric(Soctrend[6786,]))

## the water use trend ----------
RLGPDflRtpx=RLGPDflRtp-(RLGPDflRtp-RLGPDflRWtp)
RLGPDflHtpx=RLGPDflHtp-(RLGPDflRtp-RLGPDflRWtp)
WUtrend=(RLGPDflRtp-RLGPDflRWtp)
# 
# plot(as.numeric(RLGPDflRCF[6786,-71]),type="l",col=3)
# lines(as.numeric(RLGPDflH[6786,-71]),col=1)
# lines(as.numeric(RLGPDflRWCF[6786,-71]),col=2)
# lines(as.numeric(RLGPDflSCF[6786,-71]),col=4)


#plot(as.numeric(RLGPDflRCF[6786,-71])-as.numeric(RLGPDflRWCF[6786,-71]))
                                   
WUtrend=(RLGPDflRCF[,-71]-RLGPDflRWCF[,-71])

## the reservoir trend ----------
RLGPDflHtpxf=RLGPDflHtpx-(RLGPDflHtpx-RLGPDflRtpx)
Restrend=(RLGPDflHtpx-RLGPDflRtpx)


Restrend=RLGPDflH[,-71]-RLGPDflRCF[,-71]

#2 ways to buid total trend

Totaltrend=Climtrend+Soctrend+WUtrend+Restrend
#Totaltrend=RLGPDflH[,c(1:70)]-(RLGPDflH[,1])

#plot(as.numeric(Restrend[6786,]))

#Restrend[-rmat,]<-0
### Upstream area vector --------
UpAvec=UpArea[,c(12,3)]




### shape parameter vector to remove pixels with crazy shape parameter --------
Shapepar1=ParamsflH[,c(1,2,4)]
Shapepar2=ParamsflSCF[,c(1,2,4)]
Shapepar3=ParamsflRWCF[,c(1,2,4)]
Shapepar4=ParamsflRCF[,c(1,2,4)]
rmfuckers1=which(Shapepar1$epsilonGPD>1.5)
rmfuckers2=which(Shapepar2$epsilonGPD>1.5)
rmfuckers3=which(Shapepar3$epsilonGPD>1.5)
rmfuckers4=which(Shapepar4$epsilonGPD>1.5)

rmfuckers=unique(c(rmfuckers1,rmfuckers2,rmfuckers3,rmfuckers4))


## Aggregation by Regions of initial RL ----------------
data = data.frame(RLGPDflSCF)
data2 = data.frame(RLGPDflRWCF)
rmfck2=unique(ParamsflSCF$catchment[rmfuckers])

dataR=data[match(rmfck2,data$unikout),]
data=data[-match(rmfck2,data$unikout),]
data2=data2[-match(rmfck2,data2$unikout),]

#mean 10y RL over the period for each pixel

mdat=rowMeans(data[,c(1:70)],na.rm=T)
mdat2=rowMeans(data2[,c(1:70)],na.rm=T)
mdat3=(mdat+mdat2)/2

mdat2=rowMeans(data[,c(1:70)],na.rm=T)

mdat=data[,1]
mdat[(which(mdat==0))]=1
hist(mdat,xlim=c(0,1000),breaks=1000)
data$mdat=mdat
data=right_join(data,UpAvec,by = c("unikout"="outl2"))
data=data[,c(71,72,73)]


### Spatial aggregation to desired regions: Hybas07 ----------

DataI=right_join(GNF,data,by = c("outl2"="unikout"))


# Dis <- st_as_sf(DataI, coords = c("Var1", "Var2"), crs = 4326)
# 
# Dis <- st_transform(Dis, crs = 3035)
# Dis=Dis[,c(1,8,94,95,96)]
# st_write(obj=Dis,dsn=paste0(hydroDir,"/TSEVA/Init_RL_droughtRWCF.shp"))

HRM=match(DataI$HYBAS_ID,hybasf7$HYBAS_ID)

#HRM=match(DataI$NUTS3_Raster3ID,NUTS3$N3ID)
#N2M=match(DataI$NUTS3_Raster2ID,NUTS3$N2ID)

### RL in 1951 aggregated to region ----------

DataI$Hid=hybasf7$HYBAS_ID[HRM]
#DataI$NUTS3_id=NUTS3$NUTS_ID[HRM]
# DataI$NUTS2_id=NUTS3$NUTS2_ID[N2M]

#DataC$upaHR=GHshpp$SURF_KM2[HRM]
# DataCL=DataC[which(DataC$upa>DataC$upaHR),]
# DataC=DataC[-which(DataC$upa>DataC$upaHR)]

DataI$mdat=DataI$mdat * 1000 / (DataI$upa)

mean(DataI$mdat,na.rm=T)
RegioRLi=aggregate(list(RL10=DataI$mdat),
                   by = list(HydroR=DataI$Hid),
                   FUN = function(x) c(mean=mean(x,na.rm=T),dev=sd(x,na.rm=T),len=length(x),med=median(x,na.rm=T),q1=quantile(x, 0.05, na.rm=T),q3=quantile(x, 0.95, na.rm=T)))
RegioRLi <- do.call(data.frame, RegioRLi)

# ziz=match(NUTS3$NUTS_ID,NUTSRLi$HydroR)
# NUTMiss=NUTS3[which(is.na(ziz)),]
# 
# NUTS2RLi=aggregate(list(RL10=DataI$mdat),
#                    by = list(HydroR=DataI$NUTS2_id),
#                    FUN = function(x) c(mean=mean(x,na.rm=T),dev=sd(x,na.rm=T),len=length(x),med=median(x,na.rm=T),q1=quantile(x, 0.05, na.rm=T),q3=quantile(x, 0.95, na.rm=T)))
# NUTS2RLi <- do.call(data.frame, NUTS2RLi)
# 
# 
# replace= data.frame(N3=NUTMiss$NUTS_ID,NUTS2RLi[match(NUTMiss$NUTS2_ID,NUTS2RLi$HydroR),])
# replace=replace[-which(is.na(replace$RL10.mean)),]
# replace$HydroR=replace$N3
# replace=replace[,-1]
# 
# NUTSRLi=rbind(NUTSRLi,replace)




#matching biogeoregions and hybas07
bio_hybas=inner_join(biogeo_rivers,GNF,by=c("outl2"))

HRM=na.omit((match(RegioRLi$HydroR,hybasf7$HYBAS_ID)))
hybasf7_dom=hybasf7[HRM,]

bhp_m=na.omit(unique(match(hybasf7_dom$HYBAS_ID,bio_hybas$HYBAS_ID)))
hybasf7_dom$biogeoreg=bio_hybas$code[bhp_m]


#plot verification

#pointplot=inner_join(HydroRsf,pointap,by= c("Id"="Id"))


bioplot <- st_transform(hybasf7_dom, crs = 3035)
#nutplot=nutplot[-which(is.na(nutplot$maxcol)),]
br=c(-50,-20,-10,0,10,20,50)
labels=br
limi=c(-50,50)

ggplot(basemap) +
  geom_sf(fill="gray95")+
  geom_sf(fill=NA, color="grey") +
  geom_sf(data=bioplot,aes(fill=factor(biogeoreg),geometry=geometry),color="transparent",alpha=1,size=0.25,stroke=0,shape=15)+ 
  coord_sf(xlim = c(min(nco[,1]),max(nco[,1])), ylim = c(min(nco[,2]),max(nco[,2])))+
  # scale_fill_manual(values = colorz, name=" ") +
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
        legend.key.size = unit(.8, "cm"))





## Aggregation by Regions of all trends ----------------

### parameter definition
yrange=c(14:83)
yrange=c(26:95)
unikout=RLGPDflSCF[, 71]
parameters=data.frame(ParamsflSCF)
rmpixels=rmfuckers
Regio=hybasf7
RegionAggregate<- function(Drivertrend, unikout, DataI, outhybas07, 
                           parameters, rmpixels, UpAvec, GNF, Regio, 
                           GHshpp, yrange) {
  
 # Drivertrend=Restrend
  # Data preparation
  eps=0.1
  data <- data.frame(Drivertrend, unikout = unikout)
  rmfck2 <- unique(parameters$catchment[rmpixels])
  data <- data[-match(rmfck2, data$unikout), ]
  check=rowMeans(data[,-71])
  pb=which(abs(check)>10000)
  ouy=data[c(pb[1:length(pb)]),]
  data[pb,]=NA

  
  data <- right_join(data, UpAvec, by = c("unikout" = "outl2"))

  # Aggregation choice
  DataC <- right_join(GNF, data, by = c("outl2" = "unikout"))
  DataC=data.frame(DataC)
  length(which(is.na(DataC$Var2)))
  # Matching NUTS3 and NUTS2 IDs
  HRM <- match(DataC$HYBAS_ID, Regio$HYBAS_ID)
  BRM <- match(DataC$outl2, biogeo_rivers$outl2)
  # N2M <- match(DataC$NUTS3_Raster2ID, NUTS3$N2ID)
  
  # Adding NUTS3 and NUTS2 IDs to the data
  DataC$Regio_id <- Regio$HYBAS_ID[HRM]
  DataC$Biogeo_id <- biogeo_rivers$code[BRM]
  # DataC$NUTS2_id <- NUTS3$NUTS2_ID[N2M]
  
  # Normalize the data by area (upa)
  DataC[, yrange] <- DataC[, yrange] * 1000 / DataC$upa
  DataX=DataC
  
  #this line to go to relative differences
  DataX[,yrange]=DataX[,yrange]/(DataI$mdat+eps)*100
  
  # Point aggregation by NUTS3
  # Here I aggregate relative changes, while at pixel level I show absolute changes
  pointagg <- aggregate(list(Rchange_rel = DataX$Y2020,Rchange_qsp=DataC$Y2020),
                        by = list(HydroR = DataX$HYBAS_ID),
                        FUN = function(x) c(mean = mean(x, na.rm = TRUE), 
                                            dev = sd(x, na.rm = TRUE), 
                                            len = length(x), 
                                            med = median(x, na.rm = TRUE), 
                                            q1 = quantile(x, 0.05, na.rm = TRUE), 
                                            q3 = quantile(x, 0.95, na.rm = TRUE)))
  pointD <- do.call(data.frame, pointagg)
  
  
  
  # Extraction of outlets
  mdx=match(outhybas07$outID,DataC$outl2)
  pointOut <- DataC[mdx,]
  
  #pointD$Rchange.rel=pointD$Rchange.mean/RegioRLi$RL10.mean *100
  #pointD$RLi=RegioRLi$RL10.mean
  
 # ziz <- match(Regio$HYBAS_ID, pointD$HydroR)
  #NUTMiss <- NUTS3[which(is.na(ziz)), ]
  
  # # Point aggregation by NUTS2
  # pointagg2 <- aggregate(list(Rchange = DataC$Y2020),
  #                        by = list(HydroR = DataC$NUTS2_id),
  #                        FUN = function(x) c(mean = mean(x, na.rm = TRUE), 
  #                                            dev = sd(x, na.rm = TRUE), 
  #                                            len = length(x), 
  #                                            med = median(x, na.rm = TRUE), 
  #                                            q1 = quantile(x, 0.05, na.rm = TRUE), 
  #                                            q3 = quantile(x, 0.95, na.rm = TRUE)))
  # pointD2 <- do.call(data.frame, pointagg2)
  # 
  # # Replace missing NUTS3 data with NUTS2 data
  # replace <- data.frame(N3 = NUTMiss$NUTS_ID, pointD2[match(NUTMiss$NUTS2_ID, pointD2$HydroR), ])
  # replace <- replace[-which(is.na(replace$Rchange.mean)), ]
  # replace$HydroR <- replace$N3
  # replace <- replace[, -1]
  # 
  # pointD <- rbind(pointD, replace)
  
  # Climate trend aggregation by NUTS3
  trendagg <- aggregate(list(Rchange = DataC[, yrange]),
                        by = list(HydroR = DataC$HYBAS_ID),
                        FUN = function(x) c(mean = mean(x, na.rm = TRUE)))
  trendD <- do.call(data.frame, trendagg)
  
  
  # Climate trend aggregation by Biogeoregions
  # trendagg <- aggregate(list(Rchange = DataC[, yrange]),
  #                       by = list(HydroR = DataC$Biogeo_id),
  #                       FUN = function(x) c(mean = mean(x, na.rm = TRUE)))
  # trendD <- do.call(data.frame, trendagg)
  
  # # Climate trend aggregation by NUTS2
  # trendagg2 <- aggregate(list(Rchange = DataC[, yrange]),
  #                        by = list(HydroR = DataC$HYBAS_ID),
  #                        FUN = function(x) c(mean = mean(x, na.rm = TRUE)))
  # trendD2 <- do.call(data.frame, trendagg2)
  # 
  # # Replace missing NUTS3 data in trendClim with NUTS2 data
  # replace <- data.frame(N3 = NUTMiss$NUTS_ID, trendD2[match(NUTMiss$NUTS2_ID, trendD2$HydroR), ])
  # replace <- replace[-which(is.na(replace$Rchange.mean)), ]
  # replace$HydroR <- replace$N3
  # replace <- replace[, -1]
  # 
  # trendD <- rbind(trendD, replace)
  
  # Return the final trendClim as the main output
  return(list(trend=trendD,change=pointD,ChangeOut=pointOut,data=DataC))
}

TotAgg=RegionAggregate(Drivertrend=Totaltrend, unikout, DataI, outhybas07,parameters, rmpixels, UpAvec, GNF, hybasf7, GHshpp, yrange)
ClimAgg=RegionAggregate(Drivertrend=Climtrend, unikout, DataI,outhybas07,parameters, rmpixels, UpAvec, GNF, hybasf7, GHshpp, yrange)
LuAgg=RegionAggregate(Drivertrend=Soctrend, unikout, DataI,outhybas07,parameters, rmpixels, UpAvec, GNF, hybasf7, GHshpp, yrange)
ResAgg=RegionAggregate(Drivertrend=Restrend, unikout, DataI, outhybas07,parameters, rmpixels, UpAvec, GNF, hybasf7, GHshpp, yrange)
WuAgg=RegionAggregate(Drivertrend=WUtrend, unikout, DataI, outhybas07,parameters, rmpixels, UpAvec, GNF, hybasf7, GHshpp, yrange)

#### Plot 1: boxplot of the different contributions... in 2020 ----------

pointTot=TotAgg$change
pointClim=ClimAgg$change
pointSoc=LuAgg$change
pointRes=ResAgg$change
pointWu=WuAgg$change

DataT=TotAgg$data
DataC=ClimAgg$data
DataL=LuAgg$data
DataR=ResAgg$data
DataW=WuAgg$data

OutletsT=TotAgg$ChangeOut
OutletsC=ClimAgg$ChangeOut
OutletsL=LuAgg$ChangeOut
OutletsR=ResAgg$ChangeOut
OutletsW=WuAgg$ChangeOut



#match dataT with Biogeoregions

mbio=match(biogeo_rivers$outl2,DataT$outl2)
DataT$biogeo= biogeo_rivers$code[mbio]


# correction of ResTrend

RestrendCor <- data.frame(Restrend, unikout = RLGPDflH$unikout)
r_path = paste0(hydroDir,'/reservoirs/')
outletname="res_ratio_diff_2020-1951.nc"
new_reservoirs=ReservoirOpen(r_path,outletname,outf)

new_reservoirs=new_reservoirs[which(new_reservoirs$upa>0),]
rmat=match(new_reservoirs$outl2,DataW$outl2)
Rcrap=DataR[-rmat,]
Rcrap=Rcrap[which(abs(Rcrap$Y2020)>1e-1),]
#Look at rivers with reservoir influence
LUtrendX=DataL[rmat,]
LUtrendY=DataL[-rmat,]

length(LUtrendX$Y2020[which(LUtrendX$Y2020<0)])
#extract and identify these locations
LUtrend_o=LUtrendX[which(LUtrendX$Y2020<0),]
#can I relate this to the reservoirs?
pwu <- st_as_sf(Rcrap, coords = c("Var1", "Var2"), crs = 4326)
pwu <- st_transform(pwu, crs = 3035)
blss=ggplot(basemap) +
  geom_sf(fill="gray95",color="gray10",size=0.5)+
  geom_sf(data=pwu,aes(geometry=geometry,size=upa,col=Y2020),alpha=.9,stroke=0,shape=15)+
  coord_sf(xlim = c(min(nco[,1]),max(nco[,1])), ylim = c(min(nco[,2]),max(nco[,2])))+
  # scale_colour_manual(values = colorz, name="Largest change driver", labels=lab1) +
  # scale_fill_manual(values = colorz, name="Largest change driver",labels=lab1) +
  scale_size(range = c(0.08, 0.4), trans="sqrt",name= expression(paste("Upstream area ", (km^2),
                                                                       sep = " ")),
             breaks=c(101,1000,10000,100000,500000), labels=c("100","1000", "10 000", "100 000", "500 000"),
             guide = "none")+
  #ggtitle(title)+
  labs(x="Longitude", y = "Latitude")+
  guides(colour = guide_legend(override.aes = list(size = 10)))+
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

ggsave(paste0("D:/tilloal/Documents/LFRuns_utils/TrendAnalysis/plots/Blss_Rcrap.jpg"), blss,width=20, height=20, units=c("cm"),dpi=1000) 



#match this with fracwater maps


















df_2020=rbind(pointClim[,c(1,2)],pointSoc[,c(1,2)],pointRes[,c(1,2)],pointWu[,c(1,2)])
idvec=c(rep("Clim",length(pointClim$HydroR)),rep("LUC",length(pointClim$HydroR)),
        rep("Res",length(pointClim$HydroR)),rep("WU",length(pointClim$HydroR)))
names(df_2020)[2]="value"
df_2020$driver=idvec
npl="Contribution to change in 10y Flood RL(% of 1951 10Y flood) (1951-2020)"
br=seq(-100,100,by=10)
mb=seq(-100,100,5)
yl=c(-50,50)
npl="Contribution to change in 10y Flood RL(l/s/km2) (1951-2020)"
br=seq(-1,1,by=.1)
mb=seq(-1,1,.05)
yl=c(-.5,0.5)

br=seq(-100,100,by=5)
mb=seq(-100,100,1)
yl=c(-20,20)
meds <- c(by(df_2020$value, df_2020$driver, median))
q <- c(by(df_2020$value, df_2020$driver, quantile))
merdecol=match(df_2020$driver,names(meds))
df_2020$col=meds[merdecol]
palet=c(hcl.colors(9, palette = "Blues", alpha = NULL, rev = T, fixup = TRUE))
colorz = c("Res" ='tomato4',"LUC" ='orange',"Clim" ='royalblue',"WU"="limegreen")

ggplot(df_2020, aes(x=factor(driver), y=value,fill=driver)) +
  ggdist::stat_halfeye(adjust = 2, width = .6, justification = -.3, .width = .1,scale=0.6,
                       trim=F, point_colour = NA, normalize="groups") + 
  #ggdist::stat_dots(side = "left", dotsize = .1, justification = 1.1, binwidth = .1) +
  #geom_boxplot(width = .1, outlier.shape = NA) +
  geom_boxplot(notch=F,width = .1,position=position_dodge(.9),alpha=.8,aes(fill=driver),linewidth=0.8,outlier.alpha = 0.4)+
  scale_y_continuous(limits = c(-200,200),name=npl,breaks = br,minor_breaks = mb)+
  scale_x_discrete(labels=c("Res" = "Reservoirs", "LUC" = "Land use",
                            "Clim" = "Climate","WU"="Water use"),name="Driver")+
  coord_cartesian(ylim = yl)+
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


br=seq(-100,100,by=5)
labels=br
colNA="transparent"
# Plot of ordered change by region, can be important
palet=c(hcl.colors(11, palette = "RdYlBu", alpha = NULL, rev = F, fixup = TRUE))
title=paste0("Change in 10-years ",haz," Return Level between ", period[1], " and ", period[2])
legend="Change in specific discharge (l/s/km2)"
pointP=pointClim
pointP=pointP[order(pointP$Rchange.mean),]
pointP$id=c(1:length(pointP$Rchange.mean))

limi=c(-5000,5000)
ggplot() +
  geom_point(data=pointP, aes(y=id, x=Rchange.mean, color=Rchange.mean)) + 
  geom_point(data=pointP, aes(y=id, x=Rchange.mean), pch=21, colour="gray3",alpha=1) + 
  geom_vline(xintercept = 0,lwd=1, col="black")+
  scale_x_continuous(limits=limi,breaks=br)+
  scale_color_gradientn(colors=palet,
                        breaks=br,limits=c(-100,100),
                        oob = scales::squish,na.value=colNA, name=legend)+
  coord_cartesian(xlim=c(-100,100))+
  geom_segment(data=pointP,aes(y=id, yend=id, x=Rchange.q1.5.,
                    xend=Rchange.q3.95.,color=Rchange.mean), alpha=.5)+
  guides(color = guide_colourbar(barwidth = 12, barheight = 0.5,reverse=F, title.position = "top")) +
  theme(axis.title=element_text(size=14, face="bold"),
        axis.text = element_text(size=10),
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





PointAI=data.frame(clim=abs(pointClim$Rchange_qsp.mean),lu=abs(pointSoc$Rchange_qsp.mean),
                   res=abs(pointRes$Rchange_qsp.mean),wu=abs(pointWu$Rchange_qsp.mean))

max_col_numbers <- apply(PointAI, 1, function(x) which.max(x))
max_col_numbers<-as.numeric(max_col_numbers)

colorz = c("4"="limegreen","3" ='tomato4',"2" ='orange',"1" ='royalblue')


#pointap=full_join(GHshpp,pointClim,by=c("CODEB"="HydroR"))
pointap=full_join(Regio,pointClim,by=c("HYBAS_ID"="HydroR"))
st_geometry(pointap)<-NULL
#match pointagg with hybasf
cmat=match(Regio$HYBAS_ID,pointClim$HydroR)
pointap=pointap[which(!is.na(cmat)),]
cmat=cmat[which(!is.na(cmat))]
pointap$maxcol=max_col_numbers[cmat]
#Map plot of which driver is the largest
#pointap=pointap[-which(is.na(pointap$maxcol)),]

#pointplot=inner_join(HydroRsf,pointap,by= c("Id"="Id"))
pointplot=inner_join(hybasf7,pointap,by= c("HYBAS_ID"))
period=c(1951,2020)

nutplot <- st_transform(pointplot, crs = 3035)
#nutplot=nutplot[-which(is.na(nutplot$maxcol)),]
br=c(-50,-20,-10,0,10,20,50)
labels=br
limi=c(-50,50)

pl1=ggplot(basemap) +
  geom_sf(fill="gray95")+
  geom_sf(fill=NA, color="grey") +
  geom_sf(data=nutplot,aes(fill=factor(maxcol),geometry=geometry),color="transparent",alpha=1,size=0.25,stroke=0,shape=15)+ 
  coord_sf(xlim = c(min(nco[,1]),max(nco[,1])), ylim = c(min(nco[,2]),max(nco[,2])))+
  scale_fill_manual(values = colorz, name=" ") +
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
        legend.key.size = unit(.8, "cm"))
pl1

#Main driver at pixel level

#pixel level drivers 
idl=c(25,95,96)
DataC1=DataC[,idl]
DataL1=DataL[,idl]
DataW1=DataW[,idl]
DataR1=DataR[,idl]


DataAI=data.frame(clim=abs(DataC1$Y2020),lu=abs(DataL1$Y2020),
                   resw=abs(DataR1$Y2020),wu=abs(DataW1$Y2020))

max_col_numbers <- apply(DataAI, 1, function(x) which.max(x))
max_col_numbers=as.numeric(max_col_numbers)

datap=DataC
datap$maxcol=max_col_numbers
#datap=datap[-which(is.na(datap$maxcol)),]
#datap=datap[-which(is.na(datap$Var1)),]
points <- st_as_sf(datap, coords = c("Var1", "Var2"), crs = 4326)
points <- st_transform(points, crs = 3035)
points=points[-which(is.na(points$maxcol)),]
points[which(points$maxcol==4),]
if (length(which(is.na(nutplot$maxcol)))>0) nutplot=nutplot[-which(is.na(nutplot$maxcol)),]
lab1=c("Climate","Land use","Reservoirs","Water demand")

ok<-ggplot(basemap) +
  geom_sf(fill="gray95",color="gray10",size=0.5)+
  #geom_sf(data=pointsag,aes(fill=mkta,geometry=geometry),alpha=0.3,color="transparent")+
  geom_sf(data=nutplot,aes(fill=factor(maxcol),geometry=geometry),color="transparent",alpha=.6,size=0.25,stroke=0,shape=15)+ 
  geom_sf(data=points,aes(col=factor(maxcol),geometry=geometry,size=upa),alpha=.9,stroke=0,shape=15)+
  coord_sf(xlim = c(min(nco[,1]),max(nco[,1])), ylim = c(min(nco[,2]),max(nco[,2])))+
  scale_colour_manual(values = colorz, name="Largest change driver", labels=lab1) +
  scale_fill_manual(values = colorz, name="Largest change driver",labels=lab1) +
  scale_size(range = c(0.08, 0.4), trans="sqrt",name= expression(paste("Upstream area ", (km^2),
                                                                  sep = " ")),
             breaks=c(101,1000,10000,100000,500000), labels=c("100","1000", "10 000", "100 000", "500 000"),
             guide = "none")+
  #ggtitle(title)+
  labs(x="Longitude", y = "Latitude")+
  guides(colour = guide_legend(override.aes = list(size = 10)))+
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

mmx="th"

ggsave(paste0("D:/tilloal/Documents/LFRuns_utils/TrendAnalysis/plots/largestdriver",haz,mmx,"qsp.jpg"), ok, width=20, height=20, units=c("cm"),dpi=1000) 


trtest <- melt(DataC1, id.vars = "outl2", variable.name = "variable", value.name = "value")





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
    if (!is.na(miniTS[2]) & max(diff(miniTS[-1]), na.rm = TRUE) > 0) {
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
  pointsag <- inner_join(Regio, pointagg, by = c("HYBAS_ID" = "HydroR"))
  
  # Define change categories
  pointsag$change <- 0
  pointsag$change[which(pointsag$Rchange.mean > 0)] <- 1
  pointsag$change[which(pointsag$Rchange.mean < 0)] <- -1
  pointsag$change[which(pointsag$change == 1 & pointsag$Rchange.q1.5. >= 0)] <- 2
  pointsag$change[which(pointsag$change == -1 & pointsag$Rchange.q3.95. <= 0)] <- -2
  pointsag$change[which(pointsag$Rchange.mean * pointsag$Rchange.med < 0)] <- 0
  
  # Assign significance levels
  pointsag$siglvl <- 0
  pointsag$siglvl[which(pointsag$sl <= 0.05)] <- 1
  
  # Assign slope values
  pointsag$tslop <- pointsag$mkta
  pointsag$tslop[which(pointsag$siglvl < 1)] <- NA
  
  # Filter significant points and convert to spatial data
  catsig <- pointsag[which(pointsag$siglvl > 0), ]
  catsig <- st_transform(catsig, crs = 3035)
  
  # Create grid points
  if (length(catsig$HYBAS_ID)>0){
    grdpts <- sf::st_make_grid(catsig, what = "centers", cellsize = 30000)
    my.points <- sf::st_sf(grdpts)
    sf_use_s2(FALSE)
    
    # Find points inside the significant categories
    pointsInside <- sf::st_join(x = my.points, y = catsig, left = FALSE)
    pointsInside$sign <- "positive"
    pointsInside$sign[which(pointsInside$mkta <= 0)] <- "negative"
  }
  # Create pagg (significant points data)
  pagg <- catsig
  pagg$sign <- "positive"
  pagg$sign[which(pagg$tslop < 0)] <- "negative"
  
  # Create final points
  points <- st_as_sf(datap, coords = c("Var1", "Var2"), crs = 4326)
  points <- st_transform(points, crs = 3035)
  
  return(list(PagD=pointsag,points=points,psp=pointsInside))
}


yrlist=c(1951:2020)
driverlist=c("climate","landuse","reservoirs","wateruse","all")
driver=driverlist[1]
Dchangelist=list()
for (driver in driverlist){
  print(driver)
  if (driver=="climate"){
    trendPlot=ClimAgg$trend
    datap=DataC
    pointagg=pointClim
  }
  if (driver=="landuse"){
    trendPlot=LuAgg$trend
    datap=DataR
    pointagg=pointSoc
  }
  if (driver=="reservoirs"){
    trendPlot=ResAgg$trend
    datap=DataR
    pointagg=pointRes
  }
  if (driver=="wateruse"){
    trendPlot=WuAgg$trend
    datap=DataW
    pointagg=pointWu
  }
  if (driver=="all"){
    trendPlot=TotAgg$trend
    datap=DataT
    pointagg=pointTot
  }
  Pplot=calculatePoints(trendPlot, yrlist, pointagg, Regio, GHshpp, datap)
  
  save=F
  
  if (save==T){ 
    colNA="transparent"
    if (haz=="Flood"){
      if (driver=="climate" | driver=="all"){
        br=c(-50,-25,-10,-5,-2,0,2,5,10,25,50)
        labels=br
        limi=c(-50,50)
        tsize=16
        osize=12
        
        #specifically designed plot for climate
        library(ggnewscale)
        titleX=paste0("Change in 10-year ",haz," attributed \nto climatic changes (l/s/km2) - 1951-2020")
        legend="Change \n(l/s/km2)"
        legend2="Change (%)"
        palet=c(hcl.colors(11, palette = "RdYlBu", alpha = NULL, rev = F, fixup = TRUE))
        paletf=c(hcl.colors(11, palette = "RdBu", alpha = NULL, rev = F, fixup = TRUE))
        points=Pplot$points
        pag=Pplot$PagD
        pointsInside=Pplot$psp
        # quantile(points$Y2020,0.95,na.rm=T)
        ocrap<-ggplot(basemap) +
          geom_sf(fill="gray95",color="darkgrey",size=0.5)+
          geom_sf(data=pag,aes(fill=Rchange_qsp.mean,geometry=geometry),alpha=0.2,color="transparent")+
          geom_sf(data=points,aes(col=Y2020,geometry=geometry,size=upa),alpha=.9,stroke=0,shape=15)+ 
          
          scale_size(range = c(0.08, 0.4), trans="sqrt",name= expression(paste("Upstream area ", (km^2),
                                                                               sep = " ")),
                     breaks=c(101,1000,10000,100000,500000), labels=c("100","1000", "10 000", "100 000", "500 000"),
                     guide = "none")+
          scale_fill_gradientn(
            colors=paletf,
            breaks=br,limits=c(-25,25),trans=scales::modulus_trans(.3),
            oob = scales::squish,na.value=colNA, name=legend2)   +
          guides(fill = "none")+
          new_scale_fill()+
          geom_sf(data=pointsInside,aes(geometry=geometry, fill=sign),alpha=.6,size=.5,stroke=0,shape=21,color="black")+
          scale_fill_manual(values=c("tomato4","steelblue4"), name="Significant trends")+
          coord_sf(xlim = c(min(nco[,1]),max(nco[,1])), ylim = c(min(nco[,2]),max(nco[,2])))+
          scale_color_gradientn(
            colors=palet,
            breaks=br,limits=limi,trans=scales::modulus_trans(.3),
            oob = scales::squish,na.value=colNA, name=legend)   +
          labs(x="Longitude", y = "Latitude")+
          guides(colour = guide_colourbar(barwidth = 1.5, barheight = 14),
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
                legend.key.size = unit(1, "cm"))+
          ggtitle(titleX)
        
        ggsave(paste0("D:/tilloal/Documents/LFRuns_utils/TrendAnalysis/plots/mapUltimeS_climate_",haz,mmx,"qsp.jpg"), ocrap, width=22, height=20, units=c("cm"),dpi=1000) 
        
      }else{
          br=c(-20,-10,-5,-2,0,2,5,10,20)
          labels=br
          limi=c(-10,10)
          tsize=16
          osize=12
          titleX=paste0("Change in 10-year ",haz," attributed \nto Socioeconomic changes (l/s/km2) - 1951-2020")
    #custom plot combining the effect of dams at pixe level and land+water use at catchment level
          legend="Change \n(l/s/km2)"
          legend2="Change (%)"
          palet=c(hcl.colors(11, palette = "BrBg", alpha = NULL, rev = F, fixup = TRUE))
          paletf=c(hcl.colors(11, palette = "RdBu", alpha = NULL, rev = F, fixup = TRUE))
          points=Pplot$points
          pag=Pplot$PagD
          pointsInside=Pplot$psp
          # quantile(points$Y2020,0.95,na.rm=T)
          ocrap<-ggplot(basemap) +
            geom_sf(fill="gray95",color="darkgrey",size=0.5)+
            geom_sf(data=pag,aes(fill=Rchange_qsp.mean,geometry=geometry),alpha=0.6,color="transparent")+
            geom_sf(data=points,aes(col=Y2020,geometry=geometry,size=upa),alpha=.9,stroke=0,shape=15)+ 
            
            scale_size(range = c(0.1, 0.5), trans="sqrt",name= expression(paste("Upstream area ", (km^2),
                                                                                 sep = " ")),
                       breaks=c(101,1000,10000,100000,500000), labels=c("100","1000", "10 000", "100 000", "500 000"),
                       guide = "none")+
            scale_fill_gradientn(
              colors=palet,
              breaks=br,limits=c(-10,10),trans=scales::modulus_trans(.3),
              oob = scales::squish,na.value=colNA, name="Land + \nWater use ")   +
            #guides(fill = "none")+
            # new_scale_fill()+
            coord_sf(xlim = c(min(nco[,1]),max(nco[,1])), ylim = c(min(nco[,2]),max(nco[,2])))+
            scale_color_gradientn(
              colors=paletf,
              breaks=br,limits=limi,trans=scales::modulus_trans(.3),
              oob = scales::squish,na.value=colNA, name="Reservoirs \n")   +
            labs(x="Longitude", y = "Latitude")+
            guides(colour = guide_colourbar(barheight = 16, barwidth = .6),
                   fill = guide_colourbar(barheight = 16, barwidth = .6))+
            theme(axis.title=element_text(size=tsize),
                  title = element_text(size=osize),
                  axis.text=element_text(size=osize),
                  panel.background = element_rect(fill = "aliceblue", colour = "grey1"),
                  panel.border = element_rect(linetype = "solid", fill = NA, colour="black"),
                  legend.text = element_text(size=8),
                  legend.title = element_text(size = osize, margin = margin(t = 2, r = 2, b = 6, l = 0)),
                  legend.spacing.x = unit(0.2, "cm"),
                  legend.position = "right",
                  legend.box = "horizontal",  # Stack legends vertically
                  panel.grid.major = element_line(colour = "grey70"),
                  panel.grid.minor = element_line(colour = "grey90"),
                  legend.key = element_rect(fill = "transparent", colour = "transparent"),
                  legend.key.size = unit(1, "cm"))+
            ggtitle(titleX)
          
          ggsave(paste0("D:/tilloal/Documents/LFRuns_utils/TrendAnalysis/plots/mapUltime_Socio_",haz,mmx,"qsp.jpg"), ocrap, width=23, height=20, units=c("cm"),dpi=1000) 
          
          
          
      }
      }else if(haz=="Drought"){
        br=c(-50,-20,-10,0,10,20,50)
        labels=br/100
        limi=c(-50,50)
        
        tsize=16
        osize=12
        
        #specifically designed plot for climate
        library(ggnewscale)
        titleX=paste0("Change in 10-year ",haz," attributed \nto climatic changes (l/s/km2) - 1951-2020")
        legend="Change \n(l/s/km2)"
        legend2="Change (%)"
        palet=c(hcl.colors(11, palette = "RdYlBu", alpha = NULL, rev = F, fixup = TRUE))
        paletf=c(hcl.colors(11, palette = "RdBu", alpha = NULL, rev = F, fixup = TRUE))
        points=Pplot$points
        pag=Pplot$PagD
        pointsInside=Pplot$psp
        # quantile(points$Y2020,0.95,na.rm=T)
        ocrap<-ggplot(basemap) +
          geom_sf(fill="gray95",color="darkgrey",size=0.5)+
          geom_sf(data=pag,aes(fill=Rchange_qsp.mean*100,geometry=geometry),alpha=0.2,color="transparent")+
          geom_sf(data=points,aes(col=Y2020*100,geometry=geometry,size=upa),alpha=.9,stroke=0,shape=15)+ 
          
          scale_size(range = c(0.08, 0.4), trans="sqrt",name= expression(paste("Upstream area ", (km^2),
                                                                               sep = " ")),
                     breaks=c(101,1000,10000,100000,500000), labels=c("100","1000", "10 000", "100 000", "500 000"),
                     guide = "none")+
          scale_fill_gradientn(
            colors=paletf,
            breaks=br,limits=limi,labels = labels, trans=scales::modulus_trans(.5),
            oob = scales::squish,na.value=colNA, name=legend2)   +
          guides(fill = "none")+
          new_scale_fill()+
          geom_sf(data=pointsInside,aes(geometry=geometry, fill=sign),alpha=.6,size=.5,stroke=0,shape=21,color="black")+
          scale_fill_manual(values=c("tomato4","steelblue4"), name="Significant trends")+
          coord_sf(xlim = c(min(nco[,1]),max(nco[,1])), ylim = c(min(nco[,2]),max(nco[,2])))+
          scale_color_gradientn(
            colors=palet,
            breaks=br,limits=limi,labels = labels,trans=scales::modulus_trans(.5),
            oob = scales::squish,na.value=colNA, name=legend)   +
          labs(x="Longitude", y = "Latitude")+
          guides(colour = guide_colourbar(barwidth = 1.5, barheight = 14),
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
                legend.key.size = unit(1, "cm"))+
          ggtitle(titleX)
        
        ggsave(paste0("D:/tilloal/Documents/LFRuns_utils/TrendAnalysis/plots/mapUltimeS_climate_",haz,mmx,"qsp.jpg"), ocrap, width=22, height=20, units=c("cm"),dpi=1000) 
        
        
        br=c(-50,-20,-10,0,10,20,50)
        labels=br/100
        limi=c(-50,50)
        tsize=16
        osize=12
        titleX=paste0("Change in 10-year ",haz," attributed \nto Socioeconomic changes (l/s/km2) - 1951-2020")
        #custom plot combining the effect of dams at pixe level and land+water use at catchment level
        legend="Change \n(l/s/km2)"
        legend2="Change (%)"
        palet=c(hcl.colors(11, palette = "BrBg", alpha = NULL, rev = F, fixup = TRUE))
        paletf=c(hcl.colors(11, palette = "RdBu", alpha = NULL, rev = F, fixup = TRUE))
        points=Pplot$points
        pag=Pplot$PagD
        pointsInside=Pplot$psp
        # quantile(points$Y2020,0.95,na.rm=T)
        ocrap<-ggplot(basemap) +
          geom_sf(fill="gray95",color="darkgrey",size=0.5)+
          geom_sf(data=pag,aes(fill=Rchange_qsp.mean*100,geometry=geometry),alpha=0.6,color="transparent")+
          geom_sf(data=points,aes(col=Y2020*100,geometry=geometry,size=upa),alpha=.9,stroke=0,shape=15)+ 
          
          scale_size(range = c(0.1, 0.5), trans="sqrt",name= expression(paste("Upstream area ", (km^2),
                                                                              sep = " ")),
                     breaks=c(101,1000,10000,100000,500000), labels=c("100","1000", "10 000", "100 000", "500 000"),
                     guide = "none")+
          scale_fill_gradientn(
            colors=palet,
            breaks=br,limits=limi,labels = labels,trans=scales::modulus_trans(.5),
            oob = scales::squish,na.value=colNA, name="Land + \nWater use ")   +
          #guides(fill = "none")+
          # new_scale_fill()+
          coord_sf(xlim = c(min(nco[,1]),max(nco[,1])), ylim = c(min(nco[,2]),max(nco[,2])))+
          scale_color_gradientn(
            colors=paletf,
            breaks=br,limits=limi,labels = labels, trans=scales::modulus_trans(.5),
            oob = scales::squish,na.value=colNA, name="Reservoirs \n")   +
          labs(x="Longitude", y = "Latitude")+
          guides(colour = guide_colourbar(barheight = 16, barwidth = .6),
                 fill = guide_colourbar(barheight = 16, barwidth = .6))+
          theme(axis.title=element_text(size=tsize),
                title = element_text(size=osize),
                axis.text=element_text(size=osize),
                panel.background = element_rect(fill = "aliceblue", colour = "grey1"),
                panel.border = element_rect(linetype = "solid", fill = NA, colour="black"),
                legend.text = element_text(size=8),
                legend.title = element_text(size = osize, margin = margin(t = 2, r = 2, b = 6, l = 0)),
                legend.spacing.x = unit(0.2, "cm"),
                legend.position = "right",
                legend.box = "horizontal",  # Stack legends vertically
                panel.grid.major = element_line(colour = "grey70"),
                panel.grid.minor = element_line(colour = "grey90"),
                legend.key = element_rect(fill = "transparent", colour = "transparent"),
                legend.key.size = unit(1, "cm"))+
          ggtitle(titleX)
        
        ggsave(paste0("D:/tilloal/Documents/LFRuns_utils/TrendAnalysis/plots/mapUltime_Socio_",haz,mmx,"qsp.jpg"), ocrap, width=23, height=20, units=c("cm"),dpi=1000) 
        
        
        
        
      }

  #Plot parameters
  colNA="white"
  #legend="Relative change (%)"
  tsize=16
  osize=12
  
  legend="Change \n(l/s/km2)"
  legend2="Change (%)"
  palet=c(hcl.colors(11, palette = "RdYlBu", alpha = NULL, rev = F, fixup = TRUE))
  paletf=c(hcl.colors(11, palette = "RdBu", alpha = NULL, rev = F, fixup = TRUE))
  points=Pplot$points
  pag=Pplot$PagD
  # quantile(points$Y2020,0.95,na.rm=T)
  ocrap<-ggplot(basemap) +
    geom_sf(fill="gray95",color="darkgrey",size=0.5)+
    geom_sf(data=pag,aes(fill=Rchange.mean,geometry=geometry),alpha=0.4,color="transparent")+
   #geom_sf(data=pointsInside,aes(geometry=geometry, fill=sign),alpha=.7,size=.9,stroke=0,shape=21,color="black")+
    geom_sf(data=points,aes(col=Y2020,geometry=geometry,size=upa),alpha=.9,stroke=0,shape=15)+ 
  
    scale_size(range = c(0.08, 0.4), trans="sqrt",name= expression(paste("Upstream area ", (km^2),
                                                                         sep = " ")),
               breaks=c(101,1000,10000,100000,500000), labels=c("100","1000", "10 000", "100 000", "500 000"),
               guide = "none")+
    #scale_fill_manual(values=c("tomato4","steelblue4"), name="Significant trends")+
    scale_fill_gradientn(
      colors=palet,
      breaks=br,limits=limi,
      oob = scales::squish,na.value=colNA, name=legend2)   +
    coord_sf(xlim = c(min(nco[,1]),max(nco[,1])), ylim = c(min(nco[,2]),max(nco[,2])))+
    scale_color_gradientn(
      colors=palet,
      breaks=br,limits=limi,
      oob = scales::squish,na.value=colNA, name=legend)   +
    labs(x="Longitude", y = "Latitude")+
    guides(colour = guide_colourbar(barwidth = 1.5, barheight = 14),
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
          legend.key.size = unit(1, "cm"))+
    ggtitle(titleX)
  
    ggsave(paste0("D:/tilloal/Documents/LFRuns_utils/TrendAnalysis/plots/maptext_",driver,"_",haz,mmx,"qsp.jpg"), ocrap, width=20, height=20, units=c("cm"),dpi=1000) 
  }
  
  Dchangelist=c(Dchangelist,list(Pplot))
}
#create a barplot with 1 category per driver and number of locations with increasing|decreasing trend
lm=match(pointsag$HYBAS_ID,RegioRLi$HydroR)
RegioRLi$HydroR[lm]

paggC=Dchangelist[[1]]$PagD
paggL=Dchangelist[[2]]$PagD
paggR=Dchangelist[[3]]$PagD
paggW=Dchangelist[[4]]$PagD
paggC$driver="Clim"
paggL$driver="Lu"
paggR$driver="Res"
paggW$driver="Wu"
colnames(paggC)
colnames(paggR)
pointsAD=rbind(paggC,paggL,paggR,paggW)



abc=abs(pointsAD$Rchange_qsp.mean)
thcor=quantile(abc[which(abc>0)],0.1,na.rm=T)
pointsAD$change[which(abs(pointsAD$Rchange_qsp.mean)<thcor)]=0
pointsAD$Rchange_qsp.mean[which(abs(pointsAD$Rchange_qsp.mean)>1e6)]=pointsAD$Rchange_qsp.med[which(abs(pointsAD$Rchange_qsp.mean)>1e6)]
# plot(pointsag$Rchange.mean)
# sd(pointsag$Rchange.mean)
# pointsag$change=0
# pointsag$change[which(pointsag$Rchange.mean>0)]=1
# pointsag$change[which(pointsag$Rchange.mean<(-0))]=-1
# pointsag$change[which(pointsag$change==1 & pointsag$Rchange.q1.5.>0)]=2
# pointsag$change[which(pointsag$change==-1 & pointsag$Rchange.q3.95.<0)]=-2
# pointsag$change[which(pointsag$Rchange.mean*pointsag$Rchange.med<0)]=0
# hist(pointsag$Rchange.mean)


colorb=rev(c("royalblue","lightblue","lightyellow","rosybrown1","red2"))

colorbp=rev(c("darkblue","dodgerblue","orange","indianred1","brown"))

#ncol=c("negative significant","negative","none","positive","positive significant")
ncol=c("- -","-","=","+","+ +")

xlabels=c("Climate","Land use","Reservoirs","Water demand")

agbar=aggregate(list(pointsAD$Rchange_qsp.mean),
                by = list(drivers=pointsAD$driver, chsig=pointsAD$change),
                FUN = function(x) c(mean=mean(x,na.rm=T),l=length(x),ql=quantile(x,0.025,na.rm=T),qh=quantile(x,0.975,na.rm=T)))
agbar <- do.call(data.frame, agbar)
names(agbar)[c(3,4,5,6)]=c("mean","length","q1","q2")
haz="Drought"
fac=1000


# Customize the barplot
ggplot(agbar, aes(x = drivers, y = length, fill = factor(chsig))) +
  geom_bar(position="dodge", stat="identity") +
  scale_fill_manual(values = colorb, name="Change group",labels=ncol) +
  scale_y_continuous(
    name = "Cound of NUTS3",
    breaks=seq(0,2000,100),
    minor_breaks = seq(-2000,2000,100),
    sec.axis = sec_axis( transform=~.*1/fac, name="Mean change (% of 1951 10Y RL)",
                         breaks=seq(-500,500,25))
  )+
  geom_point(data=agbar, aes(x=drivers, y = fac*(mean), color = factor(chsig),group=factor(chsig)),
             position=position_dodge(width=.9),size=3) +
  guides(color = guide_legend(override.aes = list(color = colorb)))+
  scale_x_discrete(labels=xlabels)+
  geom_linerange(data=agbar,aes(x=drivers, ymin=fac*(q1),ymax=fac*(q2),color = factor(chsig),group=factor(chsig)),
                 position = position_dodge2(width = .9),lwd=1) +
  scale_color_manual(values = colorbp, name="Change group",labels=ncol,guide = "none") +
  theme(axis.title=element_text(size=14, face="bold"),
        axis.text = element_text(size=10),
        axis.text.x = element_text(size=12,face="bold"),
        panel.background = element_rect(fill = "white", colour = "white"),
        panel.grid = element_blank(),
        panel.border = element_rect(linetype = "solid", fill = NA, colour="black"),
        legend.title = element_text(size=14, face="bold"),
        legend.text = element_text(size=12),
        axis.ticks.y = element_blank(),  # Remove x-axis tickmarks
        panel.grid.major.y =element_line(color = "lightgrey"),
        panel.grid.minor.y = element_line(color = "lightgray"),
        legend.position = "right",
        panel.grid.major = element_line(colour = "grey80"),
        panel.grid.minor.x = element_line(colour = "grey90",linetype="dashed"),
        legend.key = element_rect(fill = "transparent", colour = "transparent"),
        legend.key.size = unit(.8, "cm"))

# dev.off()
ggsave("D:/tilloal/Documents/LFRuns_utils/TrendAnalysis/plots/chang_by_drivers_droughtth11.jpg", width=30, height=20, units=c("cm"),dpi=1000) 


trendClim


#redo that plot here, remove bars and show change in time
trendData=ClimAgg$trend
processTrendData <- function(trendData, DataTr, id_var = "HydroR") {
  # Melt the trend data
  # trendData=ClimAgg$trend
  trtest <- suppressWarnings(melt(trendData, id.vars = id_var, variable.name = "variable", value.name = "value"))
  
  # Convert variable to character and extract the year
  trtest$variable <- as.character(trtest$variable)
  craplife <- data.frame(strsplit(trtest$variable, ".Y"))
  trtest$yr <- as.numeric(craplife[2,])
  
  # Calculate the decade
  trtest$decad <- 10 * round(trtest$yr / 10)
  
  # Aggregate by decade and location
  trtime <- aggregate(list(value = trtest$value),
                      by = list(yr = trtest$decad, loc = trtest[[id_var]]),
                      FUN = function(x) c(mean = mean(x, na.rm = TRUE),
                                          l = length(x),
                                          ql = quantile(x, 0.05, na.rm = TRUE),
                                          qh = quantile(x, 0.95, na.rm = TRUE)))
  
  # Convert to a data frame
  tData <- do.call(data.frame, trtime)
  names(tData)[c(3, 4, 5, 6)] <- c("changeC", "length", "cq1", "cq2")
  
  #Aggregation for biogeoregions
  # I selection only decade years
  
  decades=c("Y1955","Y1965","Y1975","Y1985","Y1995","Y2005","Y2015")
  #Insipration for aggregation by Hydroregions
  
  #ideally i would compare with the approach by catchment and then by bioregion
  valcol=match(decades,colnames(DataTr))
  valcol=c(98,valcol)
  dfd=DataTr[,valcol]
  
  dftest <- suppressWarnings(melt(dfd, id.vars = "Biogeo_id", variable.name = "variable", value.name = "value"))
  
  dftime <- aggregate(list(value = dftest$value),
                      by = list(yr = dftest$variable, loc = dftest[["Biogeo_id"]]),
                      FUN = function(x) c(mean = mean(x, na.rm = TRUE),
                                          l = length(x),
                                          med= median(x, na.rm=T),
                                          ql = quantile(x, 0.25, na.rm = TRUE),
                                          qh = quantile(x, 0.75, na.rm = TRUE),
                                          w1 = quantile(x, 0.1, na.rm = TRUE),
                                          w2 = quantile(x, 0.9, na.rm = TRUE)))
  BgData <- do.call(data.frame, dftime)
  # iqr=BgData$value.qh.75.-BgData$value.ql.25.
  # BgData$w1=BgData$value.ql.25.-1.5*iqr
  # BgData$w2=BgData$value.qh.75.+1.5*iqr
  
  BgData$decad=seq(1950,2010,by=10)
  
  
  # Aggregate only by decade (without location)
  trtime_global <- aggregate(list(value = trtest$value),
                             by = list(yr = trtest$decad),
                             FUN = function(x) c(mean = mean(x, na.rm = TRUE),
                                                 l = length(x),
                                                 med= median(x, na.rm=T),
                                                 ql = quantile(x, 0.25, na.rm = TRUE),
                                                 qh = quantile(x, 0.75, na.rm = TRUE),
                                                 w1 = quantile(x, 0.1, na.rm = TRUE),
                                                 w2 = quantile(x, 0.9, na.rm = TRUE)))
  
  # Convert to a data frame
  tGlobal <- do.call(data.frame, trtime_global)
 # names(tGlobal)[c(2, 3, 4, 5)] <- c("changeC", "length", "cq1", "cq2")
  
  # Return both the processed data with location (tData) and the global trend (tGlobal)
  return(list(tData = tData, tGlobal = tGlobal, BgData=BgData))
}


TdataClim=processTrendData(ClimAgg$trend,ClimAgg$data, id_var = "HydroR")
TdataLuse=processTrendData(LuAgg$trend,LuAgg$data,id_var = "HydroR")
TdataRes=processTrendData(ResAgg$trend,ResAgg$data,id_var = "HydroR")
TdataWuse=processTrendData(WuAgg$trend,WuAgg$data,id_var = "HydroR")



bio_names=unique(biogeo$code)

#I keep only the bioregion that I like
bio_names=bio_names[c(1,3,4,6,7,9,11)]

trtF=rbind(TdataClim$tGlobal,TdataLuse$tGlobal,TdataRes$tGlobal,TdataWuse$tGlobal)

trtF$year=rep(seq(1950,2020,10),4)
trtF$driver=c(rep("Clim",8),rep("LUC",8),rep("Res",8),rep("WU",8))

names(trtF)[c(2,4,5,6,7,8)]=c("changeC","med","cq1","cq2","w1","w2")

colorz = c("Clim" ='dodgerblue4',"LUC" ='gold4',"Res" ='firebrick4',"WU"="olivedrab")
colorn = c("WU" ='limegreen',"Res" ='tomato4',"LUC" ='orange',"Clim" ='royalblue')

trtF=trtF[-which(trtF$yr==2020),]
#use pointap for bars
library(ggnewscale)
fac=1

xlabs=seq(1950,2010,10)
clabels=c("Climate","Land use","Reservoirs", "Water demand")

nplot="Mean change (% of 1951 10Y RL)"
br=seq(-50,100,10)
nplot="Mean change (l/s/km2)"
#br=seq(-.5,1,.1)
br=seq(-50,100,10)
ggplot() +
  # IQR represented as rectangles
  geom_linerange(data=trtF,aes(x=year, ymin=fac*(w1),ymax=fac*w2,color = factor(driver),group=factor(driver)),
                 position = position_dodge2(width = 9),lwd=1,alpha=0.8) +
  scale_color_manual(values = colorn, name = "Drivers", labels = clabels) +
  new_scale_color()+
  
  geom_rect(data=trtF, aes(xmin = year - 4.5, xmax = year + 4.5, 
                           ymin = cq1, ymax = cq2, fill = factor(driver), group = factor(driver)), 
            alpha = 0.5, position = position_dodge(width = 9)) +
  
  # Median as points over the IQR
  geom_point(data=trtF, aes(x = year, y = changeC, color = factor(driver), group = factor(driver)),
             position = position_dodge(width = 9), size = 3) +
  
  # Median as a horizontal line within the IQR rectangle
  geom_rect(data=trtF, aes(xmin = year - 4.5, xmax = year + 4.5, 
                           ymin = med-1e-3, ymax = med+1e-3, fill = factor(driver), group = factor(driver)), 
            alpha = 1, position = position_dodge(width = 9)) +
  
  # Y-axis settings
  scale_y_continuous(name = nplot, breaks = br) +
  
  # X-axis settings
  scale_x_continuous(breaks = xlabs, labels = xlabs, name = "Decades") +
  
  # Manual fill and color scales
  scale_fill_manual(values = colorn, name = "Drivers", labels = clabels) +
  scale_color_manual(values = colorz, name = "Drivers", labels = clabels) +
  
  # Customize legend
  guides(color = guide_legend(override.aes = list(color = colorz))) +
  
  # Theme customization
  theme(
    axis.title = element_text(size = 18, face = "bold"),
    title = element_text(size = 22, face = "bold"),
    axis.text = element_text(size = 16),
    axis.text.x = element_text(size = 16, face = "bold"),
    panel.background = element_rect(fill = "white", colour = "white"),
    panel.grid = element_blank(),
    panel.border = element_rect(linetype = "solid", fill = NA, colour = "black"),
    legend.title = element_text(size = 20, face = "bold"),
    legend.text = element_text(size = 16),
    axis.ticks.y = element_blank(),
    panel.grid.major.y = element_line(color = "lightgray"),
    panel.grid.minor.y = element_line(color = "lightgray"),
    legend.position = "right",
    panel.grid.major = element_line(colour = "grey80"),
    panel.grid.minor.x = element_line(colour = "grey90", linetype = "dashed"),
    legend.key = element_rect(fill = "transparent", colour = "transparent"),
    legend.key.size = unit(0.8, "cm")
  ) +
  
  # Add title
  ggtitle("Europe")


ggsave(paste0("D:/tilloal/Documents/LFRuns_utils/TrendAnalysis/plots/bxp_in_time_EU_",haz,"_qsp.jpg"), width=30, height=20, units=c("cm"),dpi=1000) 



library(ggnewscale)
for (bn in bio_names){
  # bn=bio_names[1]
  print(bn)
  trtF=rbind(TdataClim$BgData[which(TdataClim$BgData$loc==bn),],TdataLuse$BgData[which(TdataLuse$BgData$loc==bn),],
             TdataRes$BgData[which(TdataRes$BgData$loc==bn),],TdataWuse$BgData[which(TdataWuse$BgData$loc==bn),])
  
  
  trtF$year=rep(seq(1950,2010,10),4)
  trtF$driver=c(rep("Clim",7),rep("LUC",7),rep("Res",7),rep("WU",7))
  
  names(trtF)[c(3,5,6,7,8,9)]=c("changeC","med","cq1","cq2","w1","w2")
  
  # 
  # 
  # trtAI=data.frame(clim=abs(TdataClim$tData$changeC),lu=abs(TdataLuse$tData$changeC),
  #                    resw=abs(TdataRes$tData$changeC),wu=abs(TdataWuse$tData$changeC))
  # 
  # max_col_numbers <- apply(trtAI, 1, function(x) which.max(x))
  # 
  # total_abschange <- aggregate(list(trtF$changeC),
  #                             by = list(yr=trtF$yr),
  #                             FUN = function(x) c(sum=sum(abs(x),na.rm=T),l=length(x)))
  # total_abschange <- do.call(data.frame, total_abschange)
  # names(total_abschange)[c(2,3)]=c("sum","length")
  # 
  # v_tabc=rep(total_abschange$sum,4)
  # 
  # TdataClim$tData$maxcol=max_col_numbers
  # 
  # 
  # 
  # 
  # tCbarbis=aggregate(list(trtF$changeC),
  #                 by = list(yr=trtF$yr,mcol=trtF$driver),
  #                 FUN = function(x) c(mean=abs(x),l=length(x)))
  # tCbarbis <- do.call(data.frame, tCbarbis)
  # names(tCbarbis)[c(3,4)]=c("mean","length")
  # 
  # tCbarbis$contrib=tCbarbis$mean/v_tabc*100
  # #now aggregate tcl
  # tCbar=aggregate(list(TdataClim$tData$changeC),
  #                  by = list(yr=TdataClim$tData$yr,mcol=TdataClim$tData$maxcol),
  #                  FUN = function(x) c(mean=mean(x,na.rm=T),l=length(x)))
  # tCbar <- do.call(data.frame, tCbar)
  # names(tCbar)[c(3,4)]=c("mean","length")
  # 
  # #no influence of Water demand, I add it
  # tcp=tCbar[c(1:7),]
  # tcp$mcol=4
  # tcp$mean=0
  # tcp$length=0
  # 
  # tCbar=rbind(tCbar,tcp)
  
  colorz = c("Clim" ='dodgerblue4',"LUC" ='gold4',"Res" ='firebrick4',"WU"="olivedrab")
  colorn = c("WU" ='limegreen',"Res" ='tomato4',"LUC" ='orange',"Clim" ='royalblue')
  
  
  #use pointap for bars
  
  fac=1
  
  xlabs=seq(1950,2010,10)
  clabels=c("Climate","Land use","Reservoirs", "Water demand")
  # Customize the barplot
  
  # tCbar=tCbar[-which(tCbar$yr==2020),]
  # tCbarbis=tCbarbis[-which(tCbarbis$yr==2020),]
  # trtF=trtF[-which(trtF$yr==2020),]
  
  nplot="Mean change (% of 1951 10Y RL)"
  br=seq(-50,100,10)
  nplot="Mean change (l/s/km2)"
  br=seq(-.5,1,.1)
  #br=seq(-50,100,10)
  # ggplot() +
  #   geom_bar(data=trtF, aes(x = year, y = changeC, fill = factor(driver)),
  #                        position="dodge", stat="identity",alpha=0.6) +
  #   scale_fill_manual(values = colorn, name="Drivers",labels=clabels) +
  #   geom_point(data=trtF, aes(x=year, y = fac*(changeC), color = factor(driver),group=factor(driver)),
  #              position=position_dodge(width=9),size=3) +
  #   scale_y_continuous(name=nplot,
  #                     breaks=br )+
  #   guides(color = guide_legend(override.aes = list(color = colorz)))+
  #   scale_x_continuous(breaks=xlabs,labels=xlabs,name="Decades")+
  #   geom_linerange(data=trtF,aes(x=year, ymin=fac*(cq1-1e-5),ymax=fac*cq2,color = factor(driver),group=factor(driver)),
  #                  position = position_dodge2(width = 9),lwd=1) +
  #   scale_color_manual(values = colorz, name="Drivers",labels=clabels) +
  #   theme(axis.title=element_text(size=18, face="bold"),
  #         title = element_text(size=22, face="bold"),
  #         axis.text = element_text(size=16),
  #         axis.text.x = element_text(size=16,face="bold"),
  #         panel.background = element_rect(fill = "white", colour = "white"),
  #         panel.grid = element_blank(),
  #         panel.border = element_rect(linetype = "solid", fill = NA, colour="black"),
  #         legend.title = element_text(size=20, face="bold"),
  #         legend.text = element_text(size=16),
  #         axis.ticks.y = element_blank(),  # Remove x-axis tickmarks
  #         panel.grid.major.y = element_line(color = "lightgray"),
  #         panel.grid.minor.y = element_line(color = "lightgray"),
  #         legend.position = "right",
  #         panel.grid.major = element_line(colour = "grey80"),
  #         panel.grid.minor.x = element_line(colour = "grey90",linetype="dashed"),
  #         legend.key = element_rect(fill = "transparent", colour = "transparent"),
  #         legend.key.size = unit(.8, "cm"))+
  #   ggtitle(bn)
  
 
  
  
  ggplot() +
    # IQR represented as rectangles
    geom_linerange(data=trtF,aes(x=year, ymin=fac*(w1),ymax=fac*w2,color = factor(driver),group=factor(driver)),
                   position = position_dodge2(width = 9),lwd=1,alpha=0.8) +
    scale_color_manual(values = colorn, name = "Drivers", labels = clabels) +
    new_scale_color()+
    
    geom_rect(data=trtF, aes(xmin = year - 4.5, xmax = year + 4.5, 
                             ymin = cq1, ymax = cq2, fill = factor(driver), group = factor(driver)), 
              alpha = 0.5, position = position_dodge(width = 9)) +
    
    # Median as points over the IQR
    geom_point(data=trtF, aes(x = year, y = changeC, color = factor(driver), group = factor(driver)),
               position = position_dodge(width = 9), size = 3) +
    
    # Median as a horizontal line within the IQR rectangle
    geom_rect(data=trtF, aes(xmin = year - 4.5, xmax = year + 4.5, 
                             ymin = med-1e-3, ymax = med+1e-3, fill = factor(driver), group = factor(driver)), 
              alpha = 1, position = position_dodge(width = 9)) +
    
    # Y-axis settings
    scale_y_continuous(name = nplot, breaks = br) +
    
    # X-axis settings
    scale_x_continuous(breaks = xlabs, labels = xlabs, name = "Decades") +
    
    # Manual fill and color scales
    scale_fill_manual(values = colorn, name = "Drivers", labels = clabels) +
    scale_color_manual(values = colorz, name = "Drivers", labels = clabels) +
    
    # Customize legend
    guides(color = guide_legend(override.aes = list(color = colorz))) +
    
    # Theme customization
    theme(
      axis.title = element_text(size = 18, face = "bold"),
      title = element_text(size = 22, face = "bold"),
      axis.text = element_text(size = 16),
      axis.text.x = element_text(size = 16, face = "bold"),
      panel.background = element_rect(fill = "white", colour = "white"),
      panel.grid = element_blank(),
      panel.border = element_rect(linetype = "solid", fill = NA, colour = "black"),
      legend.title = element_text(size = 20, face = "bold"),
      legend.text = element_text(size = 16),
      axis.ticks.y = element_blank(),
      panel.grid.major.y = element_line(color = "lightgray"),
      panel.grid.minor.y = element_line(color = "lightgray"),
      legend.position = "right",
      panel.grid.major = element_line(colour = "grey80"),
      panel.grid.minor.x = element_line(colour = "grey90", linetype = "dashed"),
      legend.key = element_rect(fill = "transparent", colour = "transparent"),
      legend.key.size = unit(0.8, "cm")
    ) +
    
    # Add title
    ggtitle(bn)
  
  
  
  
  
  
  
  
   
  ggsave(paste0("D:/tilloal/Documents/LFRuns_utils/TrendAnalysis/plots/bxp_in_time",bn,"_",haz,"_qsp.jpg"), width=30, height=20, units=c("cm"),dpi=1000) 
  
}

# ggplot() +
#   geom_bar(data=tCbarbis, aes(x = yr, y = contrib, fill = factor(mcol)),
#            position="dodge", stat="identity",alpha=0.6) +
#   scale_fill_manual(values = colorn, name="Drivers",labels=clabels) +
#   scale_y_continuous(
#     name = "Percentage of total absolute change",
#     breaks=seq(0,100,20),
#     minor_breaks = seq(-100,100,10),
#     sec.axis = sec_axis( transform=~.*1, name="Mean change (% of mean 10Y RL)",
#                          breaks=seq(-50,100,5))
#   )+
#   geom_point(data=trtF, aes(x=year, y = fac*(changeC), color = factor(driver),group=factor(driver)),
#              position=position_dodge(width=8),size=3) +
#   guides(color = guide_legend(override.aes = list(color = colorz)))+
#   scale_x_continuous(breaks=xlabs,labels=xlabs,name="Decades")+
#   geom_linerange(data=trtF,aes(x=year, ymin=fac*(cq1-1e-5),ymax=fac*cq2,color = factor(driver),group=factor(driver)),
#                  position = position_dodge2(width = 8),lwd=1) +
#   scale_color_manual(values = colorz, name="Drivers",labels=clabels) +
#   theme(axis.title=element_text(size=14, face="bold"),
#         axis.text = element_text(size=12),
#         axis.text.x = element_text(size=12,face="bold"),
#         panel.background = element_rect(fill = "white", colour = "white"),
#         panel.grid = element_blank(),
#         panel.border = element_rect(linetype = "solid", fill = NA, colour="black"),
#         legend.title = element_text(size=14, face="bold"),
#         legend.text = element_text(size=12),
#         axis.ticks.y = element_blank(),  # Remove x-axis tickmarks
#         panel.grid.major.y = element_line(color = "lightgray"),
#         panel.grid.minor.y = element_line(color = "lightgray"),
#         legend.position = "right",
#         panel.grid.major = element_line(colour = "grey80"),
#         panel.grid.minor.x = element_line(colour = "grey90",linetype="dashed"),
#         legend.key = element_rect(fill = "transparent", colour = "transparent"),
#         legend.key.size = unit(.8, "cm"))
# 


#Focus on as given location Landuse trend

DataL41=DataL[floor(DataL$outl2/100000)==41,]


Nsq=41
print(Nsq)
rspace= read.csv(paste0(hydroDir,"/subspace_efas.csv"))
rspace=rspace[,-1]
nrspace=rspace[Nsq,]
#outletname="outletsv8_hybas07_01min"
#outletname="outlets_hybas09_01min"
outletname="efas_rnet_100km_01min"
outll=outletopen(hydroDir,outletname,nrspace)
cord.dec=outll[,c(2,3)]
cord.dec = SpatialPoints(cord.dec, proj4string=CRS("+proj=longlat"))
cord.UTM <- spTransform(cord.dec, CRS("+init=epsg:3035"))
nc41=cord.UTM@coords


#plot it and identify shitty points

#load trend threshold:
trenth=read.csv(file="D:/tilloal/Documents/LFRuns_utils/data/TrendAnalysis/trenTH_Histo_low_41.csv")

#plot some parameters
ParamsflSCF=data.frame(ParamsflSCF)
ParamsflSC41=ParamsflSCF[which(ParamsflSCF$Year==2020),]
mp=match(DataL41$outl2,ParamsflSC41$catchment)
ParamsflSC41=ParamsflSC41[mp,]
DataL41=inner_join(DataL41,ParamsflSC41,by=c("outl2"="catchment"))
p41 <- st_as_sf(DataL41, coords = c("Var1", "Var2"), crs = 4326)
p41 <- st_transform(p41, crs = 3035)


# ParamsflRWCF=data.frame(ParamsflSCF)
ParamsflRW41=ParamsflRWCF[which(ParamsflRWCF$Year==2020),]
mp=match(DataL41$outl2,ParamsflRW41$catchment)
ParamsflRW41=ParamsflRW41[mp,]

Dpar=data.frame(outlet=ParamsflRW41$catchment,
           epsilon=ParamsflRW41$epsilonGPD - ParamsflSC41$epsilonGPD,
           th=(ParamsflRW41$thresholdGPD-ParamsflSC41$thresholdGPD),
           sigma=(ParamsflRW41$sigmaGPD- ParamsflSC41$sigmaGPD))


DataL41=inner_join(Dpar,DataL41,by=c("outlet"="outl2"))

ValuScf=inner_join(DataL41,RLGPDflSCF,by=c("outlet"="unikout"))
plot(ValuScf$Y2020.y,ValuScf$Y2020.x,log="x")

thcon=inner_join(DataL41,trenth,by=c("outlets"="cid"))

plot(thcon$Th_new,ValuScf$Y2020.x,log="y")

p41 <- st_as_sf(thcon, coords = c("Var1", "Var2"), crs = 4326)
p41 <- st_transform(p41, crs = 3035)

pal=c(hcl.colors(11, palette = "RdYlBu", alpha = NULL, rev = T, fixup = TRUE))
ggplot(basemap) +
  geom_sf(fill="gray95",color="darkgrey",size=0.5)+
  geom_sf(data=p41,aes(col=Th_new,geometry=geometry,size=upa),alpha=.9,stroke=0,shape=15)+ 
  
  scale_size(range = c(0.8, 2), trans="sqrt",name= expression(paste("Upstream area ", (km^2),
                                                                       sep = " ")),
             breaks=c(101,1000,10000,100000,500000), labels=c("100","1000", "10 000", "100 000", "500 000"),
             guide = "none")+
  coord_sf(xlim = c(min(nc41[,1]),max(nc41[,1])), ylim = c(min(nc41[,2]),max(nc41[,2])))+
  scale_color_gradientn(
    colors=pal,
    limits=c(0.5,1),
    oob = scales::squish,na.value=colNA)   +
  labs(x="Longitude", y = "Latitude")+
  guides(colour = guide_colourbar(barwidth = 1.5, barheight = 14),
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


#extract pixels with highest diff

Hdiff=DataL41


#save some important outputs
TotAgg$data$driver="Total"
ClimAgg$data$driver="Clim"
ResAgg$data$driver="Reservoirs"
LuAgg$data$driver="Landuse"
WuAgg$data$driver="Wateruse"
Alltrend=rbind(ClimAgg$data,ResAgg$data,LuAgg$data,WuAgg$data,TotAgg$data)

TotAgg$ChangeOut$driver="Total"
ClimAgg$ChangeOut$driver="Clim"
ResAgg$ChangeOut$driver="Reservoirs"
LuAgg$ChangeOut$driver="Landuse"
WuAgg$ChangeOut$driver="Wateruse"

TotAgg$trend$driver="Total"
ClimAgg$trend$driver="Clim"
ResAgg$trend$driver="Reservoirs"
LuAgg$trend$driver="Landuse"
WuAgg$trend$driver="Wateruse"

trendRegio=rbind(ClimAgg$trend,LuAgg$trend,ResAgg$trend,WuAgg$trend,TotAgg$trend)

trendOutlets=rbind(ClimAgg$ChangeOut,LuAgg$ChangeOut,ResAgg$ChangeOut,WuAgg$ChangeOut,TotAgg$ChangeOut)



Output_fl_year=list(TrendPix=Alltrend,TrendRegio=trendRegio,TrendOutlets=trendOutlets,Out2020=pointsAD,DataI=DataI)
save(Output_fl_year,file=paste0(hydroDir,"/TSEVA/output_plots/outputs_flood_year_qsp3.Rdata"))

#Output_dr_nonfrost=list(TrendPix=Alltrend,TrendRegio=trendRegio,TrendOutlets=trendOutlets,Out2020=pointsAD,DataI=DataI)
#save(Output_dr_nonfrost,file=paste0(hydroDir,"/TSEVA/output_plots/outputs_drought_nonfrost_qsp3.Rdata"))
# I need to compare this with the change in reservoir influence

#Ok keep going and clean the script


# datathin=datatwin[,c(valcol2)]
#value_L_s_m2 <- (depth_mm_day / 1000) / 86400 * 1000
#tmpval=(datatwin[,valcol2]-datatwin[,crefloc])

#not change but raw values
# if (haz=="drought"){
#   data[,valcol]=tmpval*100
#   br=c(-150,-100,-50,0,50,100,150)
#   labels=br/100
#   limi=c(-150,150)
# }
# if (haz=="flood"){
#   # data[,valcol]=tmpval
#   br=c(-50,-20,-10,0,10,20,50)
#   labels=br
#   limi=c(-50,50)
# }
trans=scales::modulus_trans(.8)
colNA="gray10"
