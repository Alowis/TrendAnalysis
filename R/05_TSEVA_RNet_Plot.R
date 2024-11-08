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
#Set data directory
hydroDir<-("D:/tilloal/Documents/LFRuns_utils/data")


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
## Spatial data for catchments ----

### Hybas07 ----
Catchmentrivers7=read.csv(paste0(hydroDir,"/Catchments/from_hybas_eu_onlyid.csv"),encoding = "UTF-8", header = T, stringsAsFactors = F)
hybas07 <- read_sf(dsn = paste0(hydroDir,"/Catchments/hydrosheds/hybas_eu_lev07_v1c.shp"))
hybasf7=fortify(hybas07) 
Catamere07=inner_join(hybasf7,Catchmentrivers7,by= "HYBAS_ID")
Catamere07$llcoord=paste(round(Catamere07$POINT_X,4),round(Catamere07$POINT_Y,4),sep=" ") 
cst7=inner_join(Catamere07,outf,by= c("llcoord"="latlong"))
cst7=st_transform(cst7,  crs=3035)


### HydroRegions ----

GridHR=raster( paste0(hydroDir,"/HydroRegions_raster_WGS84.tif"))
GHR=as.data.frame(GridHR,xy=T)
GHR=GHR[which(!is.na(GHR[,3])),]
GHR$llcoord=paste(round(GHR$x,4),round(GHR$y,4),sep=" ") 
GHR_riv=inner_join(GHR,outf,by= c("llcoord"="latlong"))
GHshpp <- read_sf(dsn ="Z:/ClimateRun4/nahaUsers/tilloal/HydroRegions/her_all_adjusted.shp")
HydroRsf=fortify(GHshpp) 


### NUTS3 ----


# NUTS3 <- read_sf(dsn = paste0(hydroDir,"/Countries/NUTS3/NUTS3_Extended_domain.shp"))
# NUTS3$N3ID=c(1:length(NUTS3$NUTS_ID))
# N2ID=unique(NUTS3$NUTS2_ID)
# N2IDn=c(1:length(N2ID))
# mati=match(NUTS3$NUTS2_ID,N2ID)
# NUTS3$N2ID=N2IDn[mati]
# st_write(NUTS3, paste0(hydroDir,"/Countries/NUTS3/NUTS3_modified.shp"), driver = "ESRI Shapefile")

NUTS3 <- read_sf(dsn = paste0(hydroDir,"/Countries/NUTS3/NUTS3_modified.shp"))
GridNUTS3=raster( paste0(hydroDir,"/Countries/NUTS3/NUTS3_Raster3ID.tif"))
GN3=as.data.frame(GridNUTS3,xy=T)
GN3=GN3[which(!is.na(GN3[,3])),]
GN3$llcoord=paste(round(GN3$x,4),round(GN3$y,4),sep=" ") 
GN3_riv=right_join(GN3,outf,by= c("llcoord"="latlong"))

GridNUTS2=raster( paste0(hydroDir,"/Countries/NUTS3/NUTS3_Raster2ID.tif"))
GN2=as.data.frame(GridNUTS2,xy=T)
GN2=GN2[which(!is.na(GN2[,3])),]
GN2$llcoord=paste(round(GN2$x,4),round(GN2$y,4),sep=" ") 
GN2_riv=right_join(GN2,outf,by= c("llcoord"="latlong"))

GNF=right_join(GN3,GN2_riv,by="llcoord")

GNUTS3sf=fortify(NUTS3) 

GNFx=GNF[which(is.na(GNF$NUTS3_Raster3ID)),]
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
haz="Flood"
namefile="flood.Histo4"
load(file=paste0(hydroDir,"/",haz,"/params.",namefile,".Rdata"))
load(file=paste0(hydroDir,"/",haz,"/RL100.",namefile,".Rdata"))
#Peaksave=data.table(Peaksave)
gc()

Paramsfl=data.table(Paramsfl[,-c(4:9,17)])
ParamsflH=Paramsfl
RLGPDflH=RLGPDfl
rm(Paramsfl,RLGPDfl)
gc()

#load results from Socio-CF run
namefile="flood.SocCF4"
load(file=paste0(hydroDir,"/",haz,"/params.",namefile,".Rdata"))
load(file=paste0(hydroDir,"/",haz,"/RL100.",namefile,".Rdata"))

RLGPDflSCF=RLGPDfl
Paramsfl=Paramsfl[,-c(4:9,17)]
ParamsflSCF=data.table(Paramsfl)
rm(Paramsfl,RLGPDfl)
gc()

#load results from Res+WU CF run
namefile="flood.RWCF4"
load(file=paste0(hydroDir,"/",haz,"/params.",namefile,".Rdata"))
load(file=paste0(hydroDir,"/",haz,"/RL100.",namefile,".Rdata"))

RLGPDflRWCF=RLGPDfl
Paramsfl=Paramsfl[,-c(4:9,17)]
ParamsflRWCF=data.table(Paramsfl)
rm(Paramsfl,RLGPDfl)
gc()

#load results from Res CF run
namefile="flood.RCF4"
load(file=paste0(hydroDir,"/",haz,"/params.",namefile,".Rdata"))
load(file=paste0(hydroDir,"/",haz,"/RL100.",namefile,".Rdata"))

RLGPDflRCF=RLGPDfl
Paramsfl=Paramsfl[,-c(4:9,17)]
ParamsflRCF=data.table(Paramsfl)
rm(Paramsfl,RLGPDfl)
gc()



RLGPDflSCF=as.data.frame(RLGPDflSCF)
RLGPDflH=as.data.frame(RLGPDflH)
RLGPDflRCF=as.data.frame(RLGPDflRCF)
RLGPDflRWCF=as.data.frame(RLGPDflRWCF)


### For drought ----------------

if (haz=="Drought"){
#Histo
  RLGPDflH[,c(1:70)]=-RLGPDflH[,c(1:70)]
  for (id in c(1:70)){
    RLGPDflH[which(RLGPDflH[,id]<0),id]=0
    RLGPDflH[which(is.infinite(RLGPDflH[,id])),id]=NA
  }
  
  
  RLGPDflRCF[,c(1:70)]=-RLGPDflRCF[,c(1:70)]
  for (id in c(1:70)){
    RLGPDflRCF[which(RLGPDflRCF[,id]<0),id]=0
    RLGPDflRCF[which(is.infinite(RLGPDflRCF[,id])),id]=NA
  }
  
  
  RLGPDflSCF[,c(1:70)]=-RLGPDflSCF[,c(1:70)]
  for (id in c(1:70)){
    RLGPDflSCF[which(RLGPDflSCF[,id]<0),id]=0
    RLGPDflSCF[which(is.infinite(RLGPDflSCF[,id])),id]=NA
  }
  
  
  RLGPDflRWCF[,c(1:70)]=-RLGPDflRWCF[,c(1:70)]
  for (id in c(1:70)){
    RLGPDflRWCF[which(RLGPDflRWCF[,id]<0),id]=0
    RLGPDflRWCF[which(is.infinite(RLGPDflRWCF[,id])),id]=NA
  }
}
## Change attribution-----------------------


### change from climate--------------------
Climtrend=RLGPDflSCF[,c(1:70)]-(RLGPDflSCF[,1])

mato=match(RLGPDflSCF[,71],RLGPDflH[,71])
#1 I remove the climate trend
RLGPDflSt=RLGPDflSCF[,c(1:70)]-Climtrend
RLGPDflRWt=RLGPDflRWCF[mato,c(1:70)]-Climtrend
RLGPDflRt=RLGPDflRCF[mato,c(1:70)]-Climtrend
RLGPDflHt=RLGPDflH[mato,c(1:70)]-Climtrend

### change from socoeconomy -----------


## socioeconomic trend (land use) -----------
RLGPDflRWtp=RLGPDflRWt-(RLGPDflRWt-RLGPDflSt)
RLGPDflRtp=RLGPDflRt-(RLGPDflRt-RLGPDflSt)
RLGPDflHtp=RLGPDflHt-(RLGPDflRWt-RLGPDflSt)

Soctrend=(RLGPDflRWt-RLGPDflSt)

plot(as.numeric(Soctrend[6786,]))

## the water use trend ----------
RLGPDflRtpx=RLGPDflRtp-(RLGPDflRtp-RLGPDflRWtp)
RLGPDflHtpx=RLGPDflHtp-(RLGPDflRtp-RLGPDflRWtp)
WUtrend=(RLGPDflRtp-RLGPDflRWtp)

plot(as.numeric(Restrend[68845,]))

## the reservoir trend ----------
RLGPDflHtpxf=RLGPDflHtpx-(RLGPDflHtpx-RLGPDflRtpx)
Restrend=(RLGPDflHtpx-RLGPDflRtpx)

plot(as.numeric(Restrend[68640,]))
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

rmfuckers=unique(c(rmfuckers1,rmfuckers2,rmfuckers3))


## Aggregation by Regions of initial RL ----------------
data = data.frame(RLGPDflSCF$Y1951,unikout=RLGPDflSCF[,71])
rmfck2=unique(ParamsflSCF$catchment[rmfuckers])
data=data[-match(rmfck2,data$unikout),]
data=right_join(data,UpAvec,by = c("unikout"="outl2"))


### Spatial aggregation to desired regions: NUTS3 ----------

DataI=right_join(GNF,data,by = c("outl2"="unikout"))

HRM=match(DataI$NUTS3_Raster3ID,NUTS3$N3ID)
N2M=match(DataI$NUTS3_Raster2ID,NUTS3$N2ID)

### RL in 1951 aggregated to NUTS3 ----------
DataI$NUTS3_id=NUTS3$NUTS_ID[HRM]
DataI$NUTS2_id=NUTS3$NUTS2_ID[N2M]

#DataC$upaHR=GHshpp$SURF_KM2[HRM]
# DataCL=DataC[which(DataC$upa>DataC$upaHR),]
# DataC=DataC[-which(DataC$upa>DataC$upaHR)]

DataI$RLGPDflSCF.Y1951=DataI$RLGPDflSCF.Y1951*1000/(DataI$upa)

NUTSRLi=aggregate(list(RL10=DataI$RLGPDflSCF.Y1951),
                   by = list(HydroR=DataI$NUTS3_id),
                   FUN = function(x) c(mean=mean(x,na.rm=T),dev=sd(x,na.rm=T),len=length(x),med=median(x,na.rm=T),q1=quantile(x, 0.05, na.rm=T),q3=quantile(x, 0.95, na.rm=T)))
NUTSRLi <- do.call(data.frame, NUTSRLi)

ziz=match(NUTS3$NUTS_ID,NUTSRLi$HydroR)
NUTMiss=NUTS3[which(is.na(ziz)),]

NUTS2RLi=aggregate(list(RL10=DataI$RLGPDflSCF.Y1951),
                   by = list(HydroR=DataI$NUTS2_id),
                   FUN = function(x) c(mean=mean(x,na.rm=T),dev=sd(x,na.rm=T),len=length(x),med=median(x,na.rm=T),q1=quantile(x, 0.05, na.rm=T),q3=quantile(x, 0.95, na.rm=T)))
NUTS2RLi <- do.call(data.frame, NUTS2RLi)


replace= data.frame(N3=NUTMiss$NUTS_ID,NUTS2RLi[match(NUTMiss$NUTS2_ID,NUTS2RLi$HydroR),])
replace=replace[-which(is.na(replace$RL10.mean)),]
replace$HydroR=replace$N3
replace=replace[,-1]

NUTSRLi=rbind(NUTSRLi,replace)



## Aggregation by Regions of all trends ----------------

### parameter definition
yrange=c(14:83)
unikout=RLGPDflSCF[, 71]
parameters=ParamsflSCF
rmpixels=rmfuckers

RegionAggregate<- function(Drivertrend, unikout, parameters, rmpixels, UpAvec, GNF, NUTS3, GHshpp, yrange) {
  
  # Data preparation
  data <- data.frame(Drivertrend, unikout = unikout)
  rmfck2 <- unique(parameters$catchment[rmpixels])
  data <- data[-match(rmfck2, data$unikout), ]
  data <- right_join(data, UpAvec, by = c("unikout" = "outl2"))
  
  # Aggregation choice
  DataC <- right_join(GNF, data, by = c("outl2" = "unikout"))
  
  # Matching NUTS3 and NUTS2 IDs
  HRM <- match(DataC$NUTS3_Raster3ID, NUTS3$N3ID)
  N2M <- match(DataC$NUTS3_Raster2ID, NUTS3$N2ID)
  
  # Adding NUTS3 and NUTS2 IDs to the data
  DataC$NUTS3_id <- NUTS3$NUTS_ID[HRM]
  DataC$NUTS2_id <- NUTS3$NUTS2_ID[N2M]
  
  # Normalize the data by area (upa)
  DataC[, yrange] <- DataC[, yrange] * 1000 / DataC$upa
  
  # Point aggregation by NUTS3
  pointagg <- aggregate(list(Rchange = DataC$Y2020),
                        by = list(HydroR = DataC$NUTS3_id),
                        FUN = function(x) c(mean = mean(x, na.rm = TRUE), 
                                            dev = sd(x, na.rm = TRUE), 
                                            len = length(x), 
                                            med = median(x, na.rm = TRUE), 
                                            q1 = quantile(x, 0.05, na.rm = TRUE), 
                                            q3 = quantile(x, 0.95, na.rm = TRUE)))
  pointD <- do.call(data.frame, pointagg)
  
  # Handling missing NUTS3 IDs
  ziz <- match(NUTS3$NUTS_ID, pointClim$HydroR)
  NUTMiss <- NUTS3[which(is.na(ziz)), ]
  
  # Point aggregation by NUTS2
  pointagg2 <- aggregate(list(Rchange = DataC$Y2020),
                         by = list(HydroR = DataC$NUTS2_id),
                         FUN = function(x) c(mean = mean(x, na.rm = TRUE), 
                                             dev = sd(x, na.rm = TRUE), 
                                             len = length(x), 
                                             med = median(x, na.rm = TRUE), 
                                             q1 = quantile(x, 0.05, na.rm = TRUE), 
                                             q3 = quantile(x, 0.95, na.rm = TRUE)))
  pointD2 <- do.call(data.frame, pointagg2)
  
  # Replace missing NUTS3 data with NUTS2 data
  replace <- data.frame(N3 = NUTMiss$NUTS_ID, pointD2[match(NUTMiss$NUTS2_ID, pointD2$HydroR), ])
  replace <- replace[-which(is.na(replace$Rchange.mean)), ]
  replace$HydroR <- replace$N3
  replace <- replace[, -1]
  
  pointD <- rbind(pointD, replace)
  
  # Climate trend aggregation by NUTS3
  trendagg <- aggregate(list(Rchange = DataC[, yrange]),
                        by = list(HydroR = DataC$NUTS3_id),
                        FUN = function(x) c(mean = mean(x, na.rm = TRUE)))
  trendD <- do.call(data.frame, trendagg)
  
  # Climate trend aggregation by NUTS2
  trendagg2 <- aggregate(list(Rchange = DataC[, yrange]),
                         by = list(HydroR = DataC$NUTS2_id),
                         FUN = function(x) c(mean = mean(x, na.rm = TRUE)))
  trendD2 <- do.call(data.frame, trendagg2)
  
  # Replace missing NUTS3 data in trendClim with NUTS2 data
  replace <- data.frame(N3 = NUTMiss$NUTS_ID, trendD2[match(NUTMiss$NUTS2_ID, trendD2$HydroR), ])
  replace <- replace[-which(is.na(replace$Rchange.mean)), ]
  replace$HydroR <- replace$N3
  replace <- replace[, -1]
  
  trendD <- rbind(trendD, replace)
  
  # Return the final trendClim as the main output
  return(list(trend=trendD,change=pointD))
}

ClimAgg=RegionAggregate(Drivertrend=Climtrend, unikout, parameters, rmpixels, UpAvec, GNF, NUTS3, GHshpp, yrange)
LuAgg=RegionAggregate(Drivertrend=Soctrend, unikout, parameters, rmpixels, UpAvec, GNF, NUTS3, GHshpp, yrange)
ResAgg=RegionAggregate(Drivertrend=Restrend, unikout, parameters, rmpixels, UpAvec, GNF, NUTS3, GHshpp, yrange)
WuAgg=RegionAggregate(Drivertrend=WUtrend, unikout, parameters, rmpixels, UpAvec, GNF, NUTS3, GHshpp, yrange)

#### Plot 1: boxplot of the different contributions... in 2020 ----------

pointClim=ClimAgg$change
pointSoc=LuAgg$change
pointRes=ResAgg$change
pointWu=WuAgg$change
df_2020=rbind(pointClim[,c(1,2)],pointSoc[,c(1,2)],pointRes[,c(1,2)],pointWu[,c(1,2)])
idvec=c(rep("Clim",length(pointClim$HydroR)),rep("LUC",length(pointClim$HydroR)),
        rep("Res",length(pointClim$HydroR)),rep("WU",length(pointClim$HydroR)))
names(df_2020)[2]="value"
df_2020$driver=idvec

meds <- c(by(df_2020$value, df_2020$driver, median))
q <- c(by(df_2020$value, df_2020$driver, quantile))
merdecol=match(df_2020$driver,names(meds))
df_2020$col=meds[merdecol]
palet=c(hcl.colors(9, palette = "Blues", alpha = NULL, rev = T, fixup = TRUE))
colorz = c("Res" ='tomato4',"LUC" ='orange',"Clim" ='royalblue',"WU"="limegreen")

ggplot(df_2020, aes(x=factor(driver), y=value,fill=driver)) +
  ggdist::stat_halfeye(adjust = 1, width = .6, justification = -.3, .width = .1,scale=0.6,
                       trim=F, point_colour = NA, normalize="groups") + 
  #ggdist::stat_dots(side = "left", dotsize = .1, justification = 1.1, binwidth = .1) +
  #geom_boxplot(width = .1, outlier.shape = NA) +
  geom_boxplot(notch=F,width = .1,position=position_dodge(.9),alpha=.8,aes(fill=driver),linewidth=0.8,outlier.alpha = 0.4)+
  scale_y_continuous(limits = c(-200,200),name="Contribution to change in 10y Flood RL(l/s/km2) (1951-2020)",breaks = seq(-100,100,by=10),minor_breaks = seq(-100,100,5))+
  scale_x_discrete(labels=c("Res" = "Reservoirs", "LUC" = "Land use",
                            "Clim" = "Climate","WU"="Water use"),name="Driver")+
  coord_cartesian(ylim = c(-30,30))+
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


br=c(-20,-15,-10,-5,0,5,10,15,20)
labels=br
limi=c(-50,50)
colNA="transparent"
# Plot of ordered change by region, can be important
palet=c(hcl.colors(11, palette = "RdYlBu", alpha = NULL, rev = F, fixup = TRUE))
title=paste0("Change in 10-years ",haz," Return Level between ", period[1], " and ", period[2])
legend="Change in specific discharge (l/s/km2)"
pointP=pointClim
pointP=pointP[order(pointP$Rchange.mean),]
pointP$id=c(1:length(pointP$Rchange.mean))

limi=c(-1000,1000)
ggplot() +
  geom_point(data=pointP, aes(y=id, x=Rchange.mean, color=Rchange.mean)) + 
  geom_point(data=pointP, aes(y=id, x=Rchange.mean), pch=21, colour="gray3",alpha=1) + 
  geom_vline(xintercept = 0,lwd=1, col="black")+
  scale_x_continuous(limits=limi,breaks=br)+
  scale_color_gradientn(colors=palet,
                        breaks=br,limits=c(-10,10),
                        oob = scales::squish,na.value=colNA, name=legend)+
  coord_cartesian(xlim=c(-50,50))+
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





PointAI=data.frame(clim=abs(pointClim$Rchange.mean),lu=abs(pointSoc$Rchange.mean),
                   res=abs(pointRes$Rchange.mean),wu=abs(pointWu$Rchange.mean))

max_col_numbers <- apply(PointAI, 1, function(x) which.max(x))

colorz = c("4"="limegreen","3" ='tomato4',"2" ='orange',"1" ='royalblue')


#pointap=full_join(GHshpp,pointClim,by=c("CODEB"="HydroR"))
pointap=full_join(NUTS3,pointClim,by=c("NUTS_ID"="HydroR"))
st_geometry(pointap)<-NULL
#match pointagg with hybasf
cmat=match(NUTS3$NUTS_ID,pointClim$HydroR)
pointap$maxcol=max_col_numbers[cmat]
#Map plot of which driver is the largest
pointap=pointap[-which(is.na(pointap$maxcol)),]

#pointplot=inner_join(HydroRsf,pointap,by= c("Id"="Id"))
pointplot=inner_join(GNUTS3sf,pointap,by= c("NUTS_ID"))
period=c(1951,2020)
haz="flood"

nutplot <- st_transform(pointplot, crs = 3035)

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

ggsave("D:/tilloal/Documents/LFRuns_utils/TrendAnalysis/plots/map_test_soc1.jpg", pl1, width=20, height=20, units=c("cm"),dpi=1000) 


#Main driver at pixel level

#pixel level drivers 

DataC1=DataC[,c(13,yrange)]
DataL1=DataL[,c(13,83,84)]
DataW1=DataW[,c(13,83,84)]
DataR1=DataR[,c(13,83,84)]


DataAI=data.frame(clim=abs(DataC1$Y2020),lu=abs(DataL1$Y2020),
                   resw=abs(DataR1$Y2020),wu=abs(DataW1$Y2020))

max_col_numbers <- apply(DataAI, 1, function(x) which.max(x))
max_col_numbers=as.numeric(max_col_numbers)

datap=DataC
datap$maxcol=max_col_numbers
datap=datap[-which(is.na(datap$maxcol)),]
points <- st_as_sf(datap, coords = c("Var1", "Var2"), crs = 4326)
points <- st_transform(points, crs = 3035)

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



ggsave("D:/tilloal/Documents/LFRuns_utils/TrendAnalysis/plots/largestdriverDroughtMMX.jpg", ok, width=20, height=20, units=c("cm"),dpi=1000) 


trtest <- melt(DataC1, id.vars = "outl2", variable.name = "variable", value.name = "value")


yrlist=c(1951:2020)
trendPlot=trendRes
datap=DataR
pointagg=pointRes
RLbase=RLGPDflSCF[,c(1,71)]
#significance at HR level
legend="Relative change (%)"
tmpval=(trendPlot[,-1])
mkta=c()
mksa=c()
chlist=c()
for (it in 1:length(pointClim[,1])){
  print(it)
  #it=48
  mks=NA
  mkt=NA
  miniTS=as.numeric(tmpval[it,])
  # plot(miniTS)
  dmt=diff(miniTS)
  #find points where change changes direction
  sign_change <- which(diff(sign(dmt)) != 0)
  sign_change <- sign_change + 2
  #plot(miniTS)
  #remove changes that occur for periods shorter than 3 years
  dsc=diff(sign_change)
  rms=which(dsc<3)
  
  sign_change=sign_change[-c(rms,(rms+1))]
  dirchange=sign(c(dmt[1],dmt[sign_change-1]))
  yrchange <- c(yrlist[1],yrlist[sign_change])
  
  if (!is.na(miniTS[2]) & max(diff(miniTS[-1]),na.rm=T)>0){
    #mk=MannKendall(miniTS[-1])
    mk2=mmkh(miniTS[-1],ci=0.95)
    mk=data.frame(tau=mk2[7],sl=mk2[2])
    
    #compute trend as well with sen.slope
    mkt=mk$tau
    mks=mk$sl
    print(mkt)
  }else{
    mkt=NA
    mks=NA
  }
  mkta=c(mkta,mkt)
  mksa=c(mksa,mks)
  
  changes=data.frame(rep(it,length(dirchange)),dirchange,yrchange)
  chlist=rbind(chlist,changes)
}

pointagg$mkta=mkta
pointagg$sl=mksa
#natch pointagg with hybasf
#pointsag=inner_join(GHshpp,pointagg,by=c("CODEB"="HydroR"))
pointsag=inner_join(NUTS3,pointagg,by=c("NUTS_ID"="HydroR"))
pointsag$change=0
pointsag$change[which(pointsag$Rchange.mean>0)]=1
pointsag$change[which(pointsag$Rchange.mean<(-0))]=-1
pointsag$change[which(pointsag$change==1 & pointsag$Rchange.q1.5.>0)]=2
pointsag$change[which(pointsag$change==-1 & pointsag$Rchange.q3.95.<0)]=-2
pointsag$change[which(pointsag$Rchange.mean*pointsag$Rchange.med<0)]=0

paggR=pointsag




#plot all the changes through all catchments
pa_l <- reshape2::melt(trendPlot, id.vars = "HydroR", variable.name = "Year", value.name = "value")
pa_l$value=as.numeric(pa_l$value)


meds <- c(by(pa_l$value, pa_l$HydroR, mean))

#for climate
br=c(-30,-20,-10,0,10,20,30)
labels=br
limi=c(-30,30)
tsize=16
osize=12
legend="Change in Qsp \n(l/s/km2)"


#for others
br=c(-10,-8,-6,-4,-2,0,2,4,6,8,10)
labels=br
limi=c(-10,10)


colNA="transparent"

palet=c(hcl.colors(11, palette = "RdYlBu", alpha = NULL, rev = F, fixup = TRUE))
paletf=c(hcl.colors(11, palette = "RdBu", alpha = NULL, rev = F, fixup = TRUE))

pointsag$siglvl=0
pointsag$siglvl[which(pointsag$sl<=0.05)]=1
pointsag$tslop=pointsag$mkta
pointsag$tslop[which(pointsag$siglvl<1)]=NA
catsig=pointsag[which(pointsag$siglvl>0),]
catsig <- st_transform(catsig, crs = 3035)
## create a grid of points
grdpts <- sf::st_make_grid(catsig, what = "centers",cellsize = 40000)
my.points <- sf::st_sf(grdpts)

sf_use_s2(FALSE)
pointsInside <- sf::st_join(x = my.points, y = catsig, left = FALSE)
pointsInside$sign="positive"
pointsInside$sign[which(pointsInside$mkta<=0)]="negative"

pag=catsig
hist(pag$change)
pag$sign="positive"
pag$sign[which((pag$tslop<0))]="negative"

pointsag$siglvl=factor(pointsag$siglvl)


#datap=datap[-which(is.na(datap$Var1)),]
points <- st_as_sf(datap, coords = c("Var1", "Var2"), crs = 4326)
points <- st_transform(points, crs = 3035)



#for drought
br=c(-1,-0.5,0,0.5,1)
labels=br
limi=c(-1,1)
tsize=16
osize=12
legend="Change in Qsp \n(l/s/km2)"


# for nonclimate drought

quantile(points$Y2020,0.00000000001,na.rm=T)
br=c(-1,-.5,0,0.5,1)
labels=br
limi=c(-1,1)

ocrap<-ggplot(basemap) +
  geom_sf(fill="gray95",color="darkgrey",size=0.5)+
  geom_sf(data=pag,aes(fill=sign,geometry=geometry),alpha=0.2,color="transparent")+
  geom_sf(data=points,aes(col=Y2020,geometry=geometry,size=upa),alpha=.9,stroke=0,shape=15)+ 
  #geom_sf(data=pointplot,aes(geometry=geometry),color="darkgrey",fill=NA)+ 
  scale_size(range = c(0.08, 0.4), trans="sqrt",name= expression(paste("Upstream area ", (km^2),
                                                                       sep = " ")),
             breaks=c(101,1000,10000,100000,500000), labels=c("100","1000", "10 000", "100 000", "500 000"),
             guide = "none")+
  # scale_fill_gradientn(
  #   colors=paletf,
  #   breaks=seq(-.1,.1,by=0.01), limits=c(-.1,.1),
  #   oob = scales::squish,na.value="transparent", name=legend, guide="none") +
  scale_fill_manual(values=c("tomato4","steelblue4"), name="Significant trends")+
  #new_scale_fill()+
  #geom_sf(data=pointsInside,aes(geometry=geometry, fill=sign),alpha=.8,size=.6,stroke=0,shape=21,color="black")+
  #geom_sf(fill=NA, color="grey") +
  #geom_sf(regiod,mapping=aes(geometry=geometry,group=name),fill=NA,color="orange")+
  coord_sf(xlim = c(min(nco[,1]),max(nco[,1])), ylim = c(min(nco[,2]),max(nco[,2])))+
  # geom_sf(data=datasl1,aes(col=fillplot,geometry=geometry),alpha=0.8,size=0.05,stroke=0,shape=1)+ 
  # geom_sf(data=datasl2,aes(col=fillplot,geometry=geometry),alpha=0.9,size=0.1,stroke=0,shape=1)+ 
  #geom_sf(data=data,aes(fill=fillplot,geometry=geometry),color="transparent")+ 
  #geom_sf(data=datasig_f,aes(size=siglvl),fill="grey",color="transparent",shape=21,alpha=0.1)+ 
  #scale_color_manual(values=c("blue" ="darkblue","red"="darkred"))+
  #geom_col_pattern(data=pointsag,aes(fill=siglvl,geometry=geometry),colour='black', pattern = 'circle') +
 # ggtitle(title)+
  # coord_sf(xlim = c(min(nco[,1]),max(nco[,1])), ylim = c(min(nco[,2]),max(nco[,2])))+
  # geom_point(data=event,aes(x=lon,y=lat,col=dr))+
  # scale_color_distiller(palette = "Spectral",
  #                       direction = -1, limit=c(0,max(eventp$dr)),oob = scales::squish)+
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
        legend.key.size = unit(1, "cm"))



ggsave("D:/tilloal/Documents/LFRuns_utils/TrendAnalysis/plots/mapR_droughtMMX.jpg", ocrap, width=20, height=20, units=c("cm"),dpi=1000) 




#create a barplot with 1 category per driver and number of locations with increasing|decreasing trend
lm=match(pointsag$NUTS_ID,NUTSRLi$HydroR)
NUTSRLi$HydroR[lm]
pointsag$NUTS_ID

paggC$driver="Clim"
paggL$driver="Lu"
paggR$driver="Res"
paggW$driver="Wu"
colnames(paggC)
colnames(paggR)
pointsAD=rbind(paggC,paggL,paggR,paggW)

pointsAD$RChange.o.mean=pointsAD$Rchange.mean/NUTSRLi$RL10.mean[lm]*100


abc=abs(pointsAD$Rchange.mean)
thcor=quantile(abc,0.05,na.rm=T)
pointsAD$change[which(abs(pointsAD$Rchange.mean)<thcor)]=0

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

agbar=aggregate(list(pointsAD$Rchange.mean),
                by = list(drivers=pointsAD$driver, chsig=pointsAD$change),
                FUN = function(x) c(mean=mean(x,na.rm=T),l=length(x),ql=quantile(x,0.025,na.rm=T),qh=quantile(x,0.975,na.rm=T)))
agbar <- do.call(data.frame, agbar)
names(agbar)[c(3,4,5,6)]=c("mean","length","q1","q2")


# Customize the barplot
ggplot(agbar, aes(x = drivers, y = length, fill = factor(chsig))) +
  geom_bar(position="dodge", stat="identity") +
  scale_fill_manual(values = colorb, name="Change group",labels=ncol) +
  scale_y_continuous(
    name = "Cound of NUTS3",
    breaks=seq(0,2000,100),
    minor_breaks = seq(-2000,2000,100),
    sec.axis = sec_axis( transform=~.*0.001, name="Mean change in Qsp (l/s/km2)",
                         breaks=seq(-2,2,.2))
  )+
  geom_point(data=agbar, aes(x=drivers, y = 1000*(mean), color = factor(chsig),group=factor(chsig)),
             position=position_dodge(width=.9),size=3) +
  guides(color = guide_legend(override.aes = list(color = colorb)))+
  scale_x_discrete(labels=xlabels)+
  geom_linerange(data=agbar,aes(x=drivers, ymin=1000*(q1),ymax=1000*(q2),color = factor(chsig),group=factor(chsig)),
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


ggsave("D:/tilloal/Documents/LFRuns_utils/TrendAnalysis/plots/chang_by_drivers_droughtMMX.jpg", width=30, height=20, units=c("cm"),dpi=1000) 


trendClim

trtest <- melt(trendClim, id.vars = "HydroR", variable.name = "variable", value.name = "value")

zob=trtest$variable[1]
trtest$variable=as.character(trtest$variable)
craplife=(strsplit(trtest$variable,".Y"))
craplife=data.frame(craplife)
trtest$yr=as.numeric(craplife[2,])
trtest$decad=10*round(trtest$yr/10)

tCL=trtest
trtime=aggregate(list(tCL$value),
                 by = list(yr=tCL$decad,loc=tCL$HydroR),
                 FUN = function(x) c(mean=mean(x,na.rm=T),l=length(x),ql=quantile(x,0.025,na.rm=T),qh=quantile(x,0.975,na.rm=T)))
tCL <- do.call(data.frame, trtime)
names(tCL)[c(3,4,5,6)]=c("changeC","length","cq1","cq2")

trtime=aggregate(list(trtest$value),
                 by = list(yr=trtest$decad),
                 FUN = function(x) c(mean=mean(x,na.rm=T),l=length(x),ql=quantile(x,0.025,na.rm=T),qh=quantile(x,0.975,na.rm=T)))
trtclim <- do.call(data.frame, trtime)
names(trtclim)[c(2,3,4,5)]=c("changeC","length","cq1","cq2")

trtest <- melt(trendSoc, id.vars = "HydroR", variable.name = "variable", value.name = "value")

zob=trtest$variable[1]

trtest$variable=as.character(trtest$variable)
craplife=(strsplit(trtest$variable,".Y"))
craplife=data.frame(craplife)
trtest$yr=as.numeric(craplife[2,])
trtest$decad=10*round(trtest$yr/10)


tLU=trtest
trtime=aggregate(list(tLU$value),
                 by = list(yr=tLU$decad,loc=tLU$HydroR),
                 FUN = function(x) c(mean=mean(x,na.rm=T),l=length(x),ql=quantile(x,0.025,na.rm=T),qh=quantile(x,0.975,na.rm=T)))
tLU <- do.call(data.frame, trtime)
names(tLU)[c(3,4,5,6)]=c("changeC","length","cq1","cq2")

trtime=aggregate(list(trtest$value),
                 by = list(yr=trtest$decad),
                 FUN = function(x) c(mean=mean(x,na.rm=T),l=length(x),ql=quantile(x,0.025,na.rm=T),qh=quantile(x,0.975,na.rm=T)))
trtLu <- do.call(data.frame, trtime)
names(trtLu)[c(2,3,4,5)]=c("changeC","length","cq1","cq2")



trtest <- melt(trendR, id.vars = "HydroR", variable.name = "variable", value.name = "value")
zob=trtest$variable[1]

trtest$variable=as.character(trtest$variable)
craplife=(strsplit(trtest$variable,".Y"))
craplife=data.frame(craplife)
trtest$yr=as.numeric(craplife[2,])
trtest$decad=10*round(trtest$yr/10)

tR=trtest
trtime=aggregate(list(tR$value),
                 by = list(yr=tR$decad,loc=tR$HydroR),
                 FUN = function(x) c(mean=mean(x,na.rm=T),l=length(x),ql=quantile(x,0.025,na.rm=T),qh=quantile(x,0.975,na.rm=T)))
tR <- do.call(data.frame, trtime)
names(tR)[c(3,4,5,6)]=c("changeC","length","cq1","cq2")


trtime=aggregate(list(trtest$value),
                 by = list(yr=trtest$decad),
                 FUN = function(x) c(mean=mean(x,na.rm=T),l=length(x),ql=quantile(x,0.025,na.rm=T),qh=quantile(x,0.975,na.rm=T)))
trtR <- do.call(data.frame, trtime)
names(trtR)[c(2,3,4,5)]=c("changeC","length","cq1","cq2")




trtest <- melt(trendW, id.vars = "HydroR", variable.name = "variable", value.name = "value")
zob=trtest$variable[1]

trtest$variable=as.character(trtest$variable)
craplife=(strsplit(trtest$variable,".Y"))
craplife=data.frame(craplife)
trtest$yr=as.numeric(craplife[2,])
trtest$decad=10*round(trtest$yr/10)

tW=trtest
trtime=aggregate(list(tW$value),
                 by = list(yr=tW$decad,loc=tW$HydroR),
                 FUN = function(x) c(mean=mean(x,na.rm=T),l=length(x),ql=quantile(x,0.025,na.rm=T),qh=quantile(x,0.975,na.rm=T)))
tW <- do.call(data.frame, trtime)
names(tW)[c(3,4,5,6)]=c("changeC","length","cq1","cq2")


trtime=aggregate(list(trtest$value),
                 by = list(yr=trtest$decad),
                 FUN = function(x) c(mean=mean(x,na.rm=T),l=length(x),ql=quantile(x,0.025,na.rm=T),qh=quantile(x,0.975,na.rm=T)))
trtW <- do.call(data.frame, trtime)
names(trtW)[c(2,3,4,5)]=c("changeC","length","cq1","cq2")




trtF=rbind(trtclim,trtLu,trtR,trtW)



trtF$year=rep(seq(1950,2020,10),4)
trtF$driver=c(rep("Clim",8),rep("LUC",8),rep("Res",8),rep("WU",8))


trtAI=data.frame(clim=abs(tCL$changeC),lu=abs(tLU$changeC),
                   resw=abs(tR$changeC),wu=abs(tW$changeC))

max_col_numbers <- apply(trtAI, 1, function(x) which.max(x))

tCL$maxcol=max_col_numbers

#now aggregate tcl
tCbar=aggregate(list(tCL$changeC),
                 by = list(yr=tCL$yr,mcol=tCL$maxcol),
                 FUN = function(x) c(mean=mean(x,na.rm=T),l=length(x)))
tCbar <- do.call(data.frame, tCbar)
names(tCbar)[c(3,4)]=c("mean","length")

#no influence of Water demand, I add it
tcp=tCbar[c(1:7),]
tcp$mcol=4
tcp$mean=0
tcp$length=0

tCbar=rbind(tCbar,tcp)

colorz = c("Clim" ='dodgerblue4',"LUC" ='gold4',"Res" ='firebrick4',"WU"="olivedrab")
colorn = c("4" ='limegreen',"3" ='tomato4',"2" ='orange',"1" ='royalblue')


#use pointap for bars


xlabs=seq(1950,2020,10)
clabels=c("Climate","Land use","Reservoirs", "Water demand")
# Customize the barplot

tCbar=tCbar[-which(tCbar$yr==2020),]
trtF=trtF[-which(trtF$yr==2020),]
ggplot() +
  geom_bar(data=tCbar, aes(x = yr, y = length, fill = factor(mcol)),
           position="dodge", stat="identity",alpha=0.6) +
  scale_fill_manual(values = colorn, name="Drivers",labels=clabels) +
  scale_y_continuous(
    name = "Cound of NUTS3",
    breaks=seq(0,2000,200),
    minor_breaks = seq(-1000,2000,100),
    sec.axis = sec_axis( transform=~.*0.001, name="Mean change in Qsp (l/s/km2)",
                         breaks=seq(-2,2,.2))
  )+
  geom_point(data=trtF, aes(x=year, y = 1000*(changeC), color = factor(driver),group=factor(driver)),
             position=position_dodge(width=8),size=3) +
  guides(color = guide_legend(override.aes = list(color = colorz)))+
  scale_x_continuous(breaks=xlabs,labels=xlabs,name="Decades")+
  geom_linerange(data=trtF,aes(x=year, ymin=1000*(cq1-1e-5),ymax=1000*cq2,color = factor(driver),group=factor(driver)),
                 position = position_dodge2(width = 8),lwd=1) +
  scale_color_manual(values = colorz, name="Drivers",labels=clabels) +
  theme(axis.title=element_text(size=14, face="bold"),
        axis.text = element_text(size=12),
        axis.text.x = element_text(size=12,face="bold"),
        panel.background = element_rect(fill = "white", colour = "white"),
        panel.grid = element_blank(),
        panel.border = element_rect(linetype = "solid", fill = NA, colour="black"),
        legend.title = element_text(size=14, face="bold"),
        legend.text = element_text(size=12),
        axis.ticks.y = element_blank(),  # Remove x-axis tickmarks
        panel.grid.major.y = element_line(color = "lightgray"),
        panel.grid.minor.y = element_line(color = "lightgray"),
        legend.position = "right",
        panel.grid.major = element_line(colour = "grey80"),
        panel.grid.minor.x = element_line(colour = "grey90",linetype="dashed"),
        legend.key = element_rect(fill = "transparent", colour = "transparent"),
        legend.key.size = unit(.8, "cm"))

ggsave("D:/tilloal/Documents/LFRuns_utils/TrendAnalysis/plots/chang_in_time_droughtMMX.jpg", width=30, height=20, units=c("cm"),dpi=1000) 




#save some important outputs
Climtrend$driver="Clim"
Restrend$driver="Reservoirs"
Soctrend$driver="Landuse"
WUtrend$driver="Wateruse"
Alltrend=rbind(Climtrend,Restrend,Soctrend,WUtrend)


trendClim$driver="Clim"
trendR$driver="Reservoirs"
trendSoc$driver="Landuse"
trendW$driver="Wateruse"

trendNUTS3=rbind(trendClim,trendSoc,trendR,trendW)
pointsAD

trtF
Output_dr_nonfrost_MMX=list(TrendPix=Alltrend,TrendNuts=trendNUTS3,Out2020=pointsAD)

save(Output_dr_nonfrost_MMX,file=paste0(hydroDir,"/TSEVA/output_plots/outputs_drought_nonfrost_MMX.Rdata"))
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



#plot shape parameter


max(Shapepar$epsilonGPD)
hist(Shapepar$epsilonGPD,breaks=100)
mean(Shapepar$epsilonGPD)

#now plot the shape parameter
Shapeplot=inner_join(Shapepar,UpArea,by=c("catchment"= "outl2"))
points <- st_as_sf(Shapeplot, coords = c("Var1.x", "Var2.x"), crs = 4326)
points <- st_transform(points, crs = 3035)

pls=ggplot(basemap) +
  geom_sf(fill="gray85")+
  geom_sf(fill=NA, color="grey") +
  geom_sf(data=points,aes(col=epsilonGPD,geometry=geometry),alpha=1,size=0.25,stroke=0,shape=15)+ 
  coord_sf(xlim = c(min(nco[,1]),max(nco[,1])), ylim = c(min(nco[,2]),max(nco[,2])))+
  scale_color_gradientn(
    colors=palet,
    breaks=c(-0.5,-0.25,0,0.25,0.5,0.75,1),limits=c(-0.5,0.5),
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
        legend.key.size = unit(.8, "cm"))

pls




#identify dominant driver for each region


ziz=full_join(GHshpp,pointagg,by=c("CODEB"="HydroR"))
st_geometry(ziz)<-NULL
#natch pointagg with hybasf
pointplot=inner_join(HydroRsf,ziz,by= c("Id"="Id"))

# legend="10y flood specific discharge (l/s/km2)"
# title=paste0("10-years ",haz," Return Level")
# br=c(1,5,10,25,100,250,1000)
# labels=br
# limi=c(1,1200)

haz="flood"
palet=c(hcl.colors(11, palette = "RdYlBu", alpha = NULL, rev = F, fixup = TRUE))
title=paste0("Change in 10-years ",haz," Return Level between ", period[1], " and ", period[2])
legend="Change in specific discharge (l/s/km2)"
if (haz=="drought"){
  data[,valcol]=tmpval*100
  br=c(-150,-100,-50,0,50,100,150)
  labels=br/100
  limi=c(-150,150)
}
if (haz=="flood"){
  # data[,valcol]=tmpval
  br=c(-20,-15,-10,5,0,5,10,15,20)
  labels=br
  limi=c(-20,20)
}
pl2=ggplot(basemap) +
  geom_sf(fill="white")+
  geom_sf(data=pointplot,aes(fill=Rchange.med,geometry=geometry),color="transparent",alpha=.5)+ 
  geom_sf(data=points,aes(col=Y2020,geometry=geometry),alpha=1,size=0.25,stroke=0,shape=15)+ 
  #geom_sf(fill=NA, color="grey") +
  geom_sf(data=pointplot,aes(geometry=geometry),color="darkgrey",fill=NA)+ 
  scale_fill_gradientn(
    colors=palet,
    breaks=br,labels=labels, limits=limi,
    oob = scales::squish,na.value=colNA, name=legend)   +
  scale_color_gradientn(
    colors=palet,
    breaks=br,limits=limi,
    oob = scales::squish,na.value=colNA, name=legend)   +
  #geom_sf(data=dt3,aes(geometry=geometry,fill=factor(plot)),color="gray12")+ 
  #scale_fill_manual(values=c("1"="purple","-1"="tan2"), name="IRES trend")+
  coord_sf(xlim = c(min(nco[,1]),max(nco[,1])), ylim = c(min(nco[,2]),max(nco[,2])))+
  labs(x="Longitude", y = "Latitude")+
  guides(color = guide_colourbar(barwidth = 20, barheight = .8), fill="none")+
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

min(points$fillplot)
ggsave("D:/tilloal/Documents/LFRuns_utils/TrendAnalysis/plots/mapHR_abs2.jpg", pl2, width=20, height=20, units=c("cm"),dpi=1000) 



pointagg=aggregate(list(Rchange=dataX[,valcol+8]),
                   by = list(HydroR=dataX$HR_id),
                   FUN = function(x) c(mean=mean(x,na.rm=T)))
pointagg <- do.call(data.frame, pointagg)


pointagp=aggregate(list(Rchange=dataX[,valcol+8]),
                   by = list(HydroR=dataX$HR_id),
                   FUN = function(x) c(mean=mean(x,na.rm=T),sd=sd(x,na.rm=T),q10=quantile(x,0.1, na.rm=T),q90=quantile(x,0.9,na.rm=T)))
pointagp <- do.call(data.frame, pointagp)


pap=data.frame(t(pointagp))

tpag=data.frame(t(pointagg))

colnames(tpag)=tpag[1,]
tpag=tpag[-1,]


# Compute the correlation matrix
cor_mat <- cor(tpag[ -c(1,2),])

# Perform hierarchical clustering
hclust_res <- hclust(as.dist(1 - cor_mat))
# Plot the dendrogram
plot(hclust_res, main = "Hierarchical Clustering Dendrogram")


# Extract the cluster assignments


print(dunn_index)
wss <- vector()
dunn_index <- vector()
for (i in 1:10) {
  cluster_assignments <- cutree(hclust_res, k = i)
  wss[i] <- sum(cor_mat[cluster_assignments == cluster_assignments[1]]^2)
  dunn_index[i] <- dunn(as.dist(1 - cor_mat),cluster_assignments)
}

# Plot the WSS as a function of the number of clusters
plot(1:10, wss, type = "b", xlab = "Number of Clusters", ylab = "WSS")

plot(2:10, dunn_index[2:10], type = "b", xlab = "Number of Clusters", ylab = "WSS")

cluster_assignments <- cutree(hclust_res, k = 4)

#legend="Kendall tau"
legend="Relative change (%)"
tmpval=(datatwin[,valcol2])
mkta=c()
mksa=c()
for (it in 1:length(pointagg[,1])){
  print(it)
  mks=NA
  mkt=NA
  miniTS=as.numeric(pointagg[it,])
  if (!is.na(miniTS[2])){
    #mk=MannKendall(miniTS[-1])
    mk2=mmkh(miniTS[-1],ci=0.95)
    mk=data.frame(tau=mk2[6],sl=mk2[2])
    
    #compute trend as well with sen.slope
    mkt=mk$tau
    mks=mk$sl
    print(mkt)
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


#plot all the changes through all catchments
pa_l <- reshape2::melt(pointagg, id.vars = "HydroR", variable.name = "Year", value.name = "value")
pa_l$value=as.numeric(pa_l$value)

tpmean=aggregate(list(rl=pa_l$value),
                 by = list(HydroR=pa_l$HydroR),
                 FUN = function(x) c(mean=mean(x,na.rm=T)))
tpmean <- do.call(data.frame, tpmean)

pa_l$value2=pa_l$value/tp
# Add the cluster assignments to the data frame
pa_l$cluster <- factor(cluster_assignments)

pointagg$mkta=mkta
pointagg$sl=mksa
meds <- c(by(pa_l$value, pa_l$HydroR, mean))

plot(meds)
pa2=pa_l[which(pa_l$cluster==4),]
ggplot(pa2, aes(x = Year, y = value, group = HydroR,col=cluster)) +
  geom_line() +
  scale_y_continuous(limits=c(-10,10))+
  theme_minimal()

plot(as.numeric(pointagg[1,-1]))

pointagg$cluster=factor(cluster_assignments)

#natch pointagg with hybasf
pointsag=inner_join(GHshpp,pointagg,by=c("CODEB"="HydroR"))



br=c(-30,-20,-10,0,10,20,30)
labels=br
limi=c(-30,30)
tsize=16
osize=12
legend="Change in Qsp \n(l/s/km2)"

ggplot(basemap) +
  geom_sf(fill="gray85",color="darkgrey",size=0.5)+
  geom_sf(data=pointsag,aes(fill=cluster,geometry=geometry),alpha=0.8,color="transparent")+
  #geom_sf(data=pointsag,aes(color=mkta,geometry=geometry),alpha=0.4,fill="transparent")+
  scale_color_gradientn(
    colors=palet,
    breaks=seq(-1,1,by=0.2), limits=c(-1,1),
    oob = scales::squish,na.value=colNA, name=legend) +
  coord_sf(xlim = c(min(nco[,1]),max(nco[,1])), ylim = c(min(nco[,2]),max(nco[,2])))+
  # geom_sf(data=datasl1,aes(col=fillplot,geometry=geometry),alpha=0.8,size=0.05,stroke=0,shape=1)+ 
  # geom_sf(data=datasl2,aes(col=fillplot,geometry=geometry),alpha=0.9,size=0.1,stroke=0,shape=1)+ 
  #geom_sf(data=data,aes(fill=fillplot,geometry=geometry),color="transparent")+ 
  #geom_sf(data=datasig_f,aes(size=siglvl),fill="grey",color="transparent",shape=21,alpha=0.1)+ 
  #scale_color_manual(values=c("blue" ="darkblue","red"="darkred"))+
  #geom_col_pattern(data=pointsag,aes(fill=siglvl,geometry=geometry),colour='black', pattern = 'circle') +
  #scale_fill_manual(values=c("tomato4","steelblue4"), name="Significant trends")+
  ggtitle(title)+
  # coord_sf(xlim = c(min(nco[,1]),max(nco[,1])), ylim = c(min(nco[,2]),max(nco[,2])))+
  # geom_point(data=event,aes(x=lon,y=lat,col=dr))+
  # scale_color_distiller(palette = "Spectral",
  #                       direction = -1, limit=c(0,max(eventp$dr)),oob = scales::squish)+
  # scale_color_gradientn(
  #   colors=palet,
  #   breaks=br,limits=limi,trans=trans,
  #   oob = scales::squish,na.value=colNA, name=legend)   +
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


#create a time vector with years
dates <- seq.Date(ymd("1951-06-01"), ymd("2020-06-01"), by = "year")
dates= as.POSIXct(dates, format = "%Y-%m-%d")
dfall=c()
plotOut=FALSE
tail="high"


# Plotting the outputs ----------------------------------------------------

## 1. selection of data to be plotted ----

### drought ----

datar1=RLGPDdr
#problem with 2020, to be solved
datar1$Y2020=datar1$Y2019
length(which(!is.na(datar1$Y2020)))
hist(datar1$Y2020,xlim=c(0,1500), breaks=1000)
#Inverting return period values
datar1[,c(1:71)]=-datar1[,c(1:71)]
datar=data.table(datar1)



### flood ----
RLGPDflSCF=RLGPDfl


datarSCF=RLGPDflSCF
datarSCF$Y2020=datarSCF$Y2019
#add latitude and longitude to input



length((datar$unikout))




datariSCF=inner_join(outf,datarSCF,by = c("outl2"="unikout"))
#datar$Y2020=datar$Y2019
unikout=datar$unikout


# datar=RLGPDfl
# datariH=inner_join(outf,datar,by = c("outlets"="unikout"))
# datar$Y2020=datar$Y2019
# unikout=datar$unikout
# 
# 
# #add latitude and longitude to input
# datari=inner_join(outf,datar,by = c("outlets"="unikout"))
# 
# datariSCF=inner_join(outf,datar,by = c("outlets"="unikout"))
# #join with catchment data
dataricat=inner_join(Catamere07,datariSCF,by = c("llcoord"="latlong"))

#checkpars=inner_join(datari,Paramsfl[which(Paramsfl$Year==2020),],by=c("outlets"="catchment"))

Impdates=seq(1950,2020,by=10)
valuenames=paste0("Y",Impdates)

## 2. different plots of the results ----

### Historical distribution of changes ----
totalch=plotHistoDates(outf,datariSCF,law="GPD",type,period=c(1951,2020),parlist=Paramsfl,valuenames)
totalch[[2]]

#trick to also correct
# v1=Paramsfl[which(Paramsfl$Year==2019),]
# v1=v1[match(unique(v1$catchment),v1$catchment),]
# v1$Year=2020
# Paramsfl[which(Paramsfl$Year==2020),]=v1

### Changes in RLs or RPs at pixel and catchment levels ----
Plot.change=plotchangemapix(basemap,catmap=cst7,datariSCF, law="GPD",type="RLchange",period=c(1950,2020),Paramsfl,hybasf = HydroRsf,haz="flood")

Plot.change[[1]]
ggsave("D:/tilloal/Documents/LFRuns_utils/Figures/RP_change_floodCal3.jpg", width=20, height=15, units=c("cm"),dpi=1500)

Plot.change[[2]]
ggsave("D:/tilloal/Documents/LFRuns_utils/Figures/RL_change_hybas07.jpg", width=20, height=15, units=c("cm"),dpi=1000)

### Changes in RLs with mk test at catchment level ----
Plot.change.sig=plotTrendSipix(basemap,dataricat,period=c(1951,2020),hybasf = hybasf7,valuenames,nco)
Plot.change.sig
ggsave("D:/tilloal/Documents/LFRuns_utils/Figures/RLchange_flood_sig19502020.jpg",Plot.change.sig, width=20, height=15, units=c("cm"),dpi=1500)



#load upstream area
main_path = 'D:/tilloal/Documents/06_Floodrivers/'
valid_path = paste0(main_path,'DataPaper/')
outletname="/GIS/upArea_European_01min.nc"
dir=valid_path
outf$idlalo=paste(outf$idlo, outf$idla, sep=" ")
UpArea=UpAopen(valid_path,outletname,outf)
head(UpArea)

Impdates=seq(1951,2020,by=1)
valuenames=paste0("Y",Impdates)
period=c(1951,2020)
haz="drought"




# TO be improved
plotchangemapix_qspU=function(basemap,catmap,datar,upArea, GHR_riv, HydroRsf, law="GPD",type,period=c(1951,2020),parlist,valuenames,haz){
  
  datar=datariSCF
  data=datar
  # data=right_join(GHR_riv,datar,by = c("outlets"="unikout"))
  # data=right_join(cst7,data,by = c("outlets"="outlets"))
  #names(data)[28]="HydroR"
  upag=match(data$outlets,UpArea$outlets)
  data$uparea=UpArea$upa[upag]
  data$outlets
  data$pointid
  datacol=names(data)
  valcol=match(valuenames,datacol)
  # data=data[which(data$IRES==1),]
  datatwin=data
  st_geometry(datatwin) <- NULL
  valcol2=valcol
  
  datatwin=as.data.frame(datatwin)
  class(datatwin)
  dtc=names(datatwin)
  # mcor=unique(match(datar$unikout,Paramsdr$catchment))
  # 
  # if (length(which(is.na(mcor)))>0) datar=datar[-which(is.na(mcor)),]
  # 
  cref=paste0("Y",period[1])
  crefloc=match(cref,dtc)
  finalperiod=paste0("Y",period[2])
  colsel=match(finalperiod,datacol)
  
  
  palet=c(hcl.colors(11, palette = "RdYlBu", alpha = NULL, rev = F, fixup = TRUE))
  title=paste0("Change in 10-years ",haz," Return Level between ", period[1], " and ", period[2])
  legend="Change in specific discharge (l/s/km2)"
 # datathin=datatwin[,c(valcol2)]
  #value_L_s_m2 <- (depth_mm_day / 1000) / 86400 * 1000
  tmpval=(datatwin[,valcol2]-datatwin[,crefloc])
  
  #not change but raw values
  #tmpval=(datatwin[,valcol2])
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
  

  
  data2=inner_join(data,catmap[,c(1:20)],by=c("latlong"="llcoord"))
  points <- st_as_sf(data2[,c(1:80)], coords = c("Var1", "Var2"), crs = 4326)
  points <- st_transform(points, crs = 3035)
  
  pl1=ggplot(basemap) +
    geom_sf(fill="gray85")+
    geom_sf(fill=NA, color="grey") +
    geom_sf(data=points,aes(col=fillplot,geometry=geometry),alpha=1,size=0.25,stroke=0,shape=15)+ 
    coord_sf(xlim = c(min(nco[,1]),max(nco[,1])), ylim = c(min(nco[,2]),max(nco[,2])))+
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
  
  pl1
  
  ggsave("D:/tilloal/Documents/LFRuns_utils/TrendAnalysis/plots/map_verifx.jpg", pl1, width=20, height=20, units=c("cm"),dpi=1000) 
  # savecrap=inner_join(data,Catamere07,by= c("HYBAS_ID"))
  # st_write(data, paste0(hydroDir,"/disasterdrought3.shp"))
  # 
  
 #  points[68270,]
 #  colnames(Paramsfl)
 #  wtff=Paramsfl[which(Paramsfl$catchment==4300027 ),]
 #  
 # wtf= datatwin[which(datatwin$outlets==4300027),]
 #plot(as.numeric(wtf[1,c(7:77)]))
  #Now aggregate by Hydroregions
  
  rmfuckers=unique(ParamsflSCF$catchment[which(ParamsflSCF$epsilonGPD>1.5)])
  
  data=data[-match(rmfuckers,data$outl2),]
  dataX=right_join(GHR_riv,data,by = c("outl2"="outl2"))
  dataX$uparea
  HRM=match(dataX$HydroRegions_raster_WGS84,GHshpp$Id)
  
  dataX$HR_id=GHshpp$CODEB[HRM]
  dataX$upaHR=GHshpp$SURF_KM2[HRM]
  
  dataXL=dataX[which(dataX$uparea>dataX$upaHR),]
  dataX=dataX[-which(dataX$uparea>dataX$upaHR)]
  pointagg=aggregate(list(Rchange=dataX$fillplot),
                     by = list(HydroR=dataX$HR_id),
                     FUN = function(x) c(mean=median(x,na.rm=T),dev=sd(x,na.rm=T),len=length(x),q1=quantile(x, 0.25, na.rm=T),q3=quantile(x, 0.75, na.rm=T)))
  pointagg <- do.call(data.frame, pointagg)
  
  ziz=full_join(GHshpp,pointagg,by=c("CODEB"="HydroR"))
  st_geometry(ziz)<-NULL
  #natch pointagg with hybasf
  pointplot=inner_join(HydroRsf,ziz,by= c("Id"="Id"))
  
  legend="10y flood specific discharge (l/s/km2)"
  title=paste0("10-years ",haz," Return Level")
  br=c(1,5,10,25,100,250,1000)
  labels=br
  limi=c(1,1200)
  
  pl2=ggplot(basemap) +
    geom_sf(fill="white")+
    geom_sf(data=pointplot,aes(fill=Rchange.mean,geometry=geometry),color="transparent",alpha=.5)+ 
    geom_sf(data=points,aes(col=fillplot,geometry=geometry),alpha=1,size=0.25,stroke=0,shape=15)+ 
    #geom_sf(fill=NA, color="grey") +
    geom_sf(data=pointplot,aes(geometry=geometry),color="darkgrey",fill=NA)+ 
    scale_fill_gradientn(
      colors=palet,
      breaks=br,labels=labels, limits=limi,trans="log",
      oob = scales::squish,na.value=colNA, name=legend)   +
    scale_color_gradientn(
      colors=palet,
      breaks=br,limits=limi,trans="log",
      oob = scales::squish,na.value=colNA, name=legend)   +
    #geom_sf(data=dt3,aes(geometry=geometry,fill=factor(plot)),color="gray12")+ 
    #scale_fill_manual(values=c("1"="purple","-1"="tan2"), name="IRES trend")+
    coord_sf(xlim = c(min(nco[,1]),max(nco[,1])), ylim = c(min(nco[,2]),max(nco[,2])))+
    labs(x="Longitude", y = "Latitude")+
    guides(color = guide_colourbar(barwidth = 20, barheight = .8), fill="none")+
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
  
  min(points$fillplot)
  ggsave("D:/tilloal/Documents/LFRuns_utils/TrendAnalysis/plots/mapHR_abs_flscf.jpg", pl2, width=20, height=20, units=c("cm"),dpi=1000) 

 
  
  pointagg=aggregate(list(Rchange=dataX[,valcol+10]),
                     by = list(HydroR=dataX$HR_id),
                     FUN = function(x) c(mean=mean(x,na.rm=T)))
  pointagg <- do.call(data.frame, pointagg)
  
  
  pointagp=aggregate(list(Rchange=dataX[,valcol+10]),
                     by = list(HydroR=dataX$HR_id),
                     FUN = function(x) c(mean=mean(x,na.rm=T),sd=sd(x,na.rm=T),q10=quantile(x,0.1, na.rm=T),q90=quantile(x,0.9,na.rm=T)))
  pointagp <- do.call(data.frame, pointagp)
  
  
  pap=data.frame(t(pointagp))
  
  tpag=data.frame(t(pointagg))
  
  colnames(tpag)=tpag[1,]
  tpag=tpag[-1,]
  
  
  # Compute the correlation matrix
  cor_mat <- cor(tpag[ -c(1,2),])
  
  # Perform hierarchical clustering
  hclust_res <- hclust(as.dist(1 - cor_mat))
  # Plot the dendrogram
  plot(hclust_res, main = "Hierarchical Clustering Dendrogram")
  
  
  # Extract the cluster assignments

  library(clValid)
  print(dunn_index)
  wss <- vector()
  dunn_index <- vector()
  for (i in 1:10) {
    cluster_assignments <- cutree(hclust_res, k = i)
    wss[i] <- sum(cor_mat[cluster_assignments == cluster_assignments[1]]^2)
    dunn_index[i] <- dunn(as.dist(1 - cor_mat),cluster_assignments)
  }
  
  # Plot the WSS as a function of the number of clusters
  plot(1:10, wss, type = "b", xlab = "Number of Clusters", ylab = "WSS")
  
  plot(2:10, dunn_index[2:10], type = "b", xlab = "Number of Clusters", ylab = "WSS")

  cluster_assignments <- cutree(hclust_res, k = 3)
  
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
  
  
  #plot all the changes through all catchments
  pa_l <- reshape2::melt(pointagg, id.vars = "HydroR", variable.name = "Year", value.name = "value")
  
  
  tpmean=aggregate(list(rl=pa_l$value),
                   by = list(HydroR=pa_l$HydroR),
                   FUN = function(x) c(mean=mean(x,na.rm=T)))
  tpmean <- do.call(data.frame, tpmean)
  
  pa_l$value2=pa_l$value/tpmean$rl
  # Add the cluster assignments to the data frame
  pa_l$cluster <- factor(cluster_assignments)
  
  pointagg$mkta=mkta
  pointagg$sl=mksa
  meds <- c(by(pa_l$value, pa_l$HydroR, mean))
  
  plot(meds)
  pa2=pa_l[which(pa_l$cluster=="4"),]
  ggplot(pa2, aes(x = Year, y = value2, group = HydroR,col=cluster)) +
    geom_line() +
    scale_y_continuous(limits=c(0.5,1.5))+
    theme_minimal()
  
  plot(as.numeric(pointagg[1,-1]))
  
  pointagg$cluster=factor(cluster_assignments)
  
  #natch pointagg with hybasf
  pointsag=inner_join(GHshpp,pointagg,by=c("CODEB"="HydroR"))
  
  
  
  br=c(-30,-20,-10,0,10,20,30)
  labels=br
  limi=c(-30,30)
  tsize=16
  osize=12
  legend="Change in Qsp \n(l/s/km2)"
  
  ggplot(basemap) +
    geom_sf(fill="gray85",color="darkgrey",size=0.5)+
    geom_sf(data=pointsag,aes(fill=cluster,geometry=geometry),alpha=0.8,color="transparent")+
    #geom_sf(data=pointsag,aes(color=mkta,geometry=geometry),alpha=0.4,fill="transparent")+
    scale_color_gradientn(
      colors=palet,
      breaks=seq(-1,1,by=0.2), limits=c(-1,1),
      oob = scales::squish,na.value=colNA, name=legend) +
    coord_sf(xlim = c(min(nco[,1]),max(nco[,1])), ylim = c(min(nco[,2]),max(nco[,2])))+
    # geom_sf(data=datasl1,aes(col=fillplot,geometry=geometry),alpha=0.8,size=0.05,stroke=0,shape=1)+ 
    # geom_sf(data=datasl2,aes(col=fillplot,geometry=geometry),alpha=0.9,size=0.1,stroke=0,shape=1)+ 
    #geom_sf(data=data,aes(fill=fillplot,geometry=geometry),color="transparent")+ 
    #geom_sf(data=datasig_f,aes(size=siglvl),fill="grey",color="transparent",shape=21,alpha=0.1)+ 
    #scale_color_manual(values=c("blue" ="darkblue","red"="darkred"))+
    #geom_col_pattern(data=pointsag,aes(fill=siglvl,geometry=geometry),colour='black', pattern = 'circle') +
    #scale_fill_manual(values=c("tomato4","steelblue4"), name="Significant trends")+
    ggtitle(title)+
    # coord_sf(xlim = c(min(nco[,1]),max(nco[,1])), ylim = c(min(nco[,2]),max(nco[,2])))+
    # geom_point(data=event,aes(x=lon,y=lat,col=dr))+
    # scale_color_distiller(palette = "Spectral",
    #                       direction = -1, limit=c(0,max(eventp$dr)),oob = scales::squish)+
    # scale_color_gradientn(
    #   colors=palet,
    #   breaks=br,limits=limi,trans=trans,
    #   oob = scales::squish,na.value=colNA, name=legend)   +
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
  
  library(ggnewscale)
  br=c(-30,-20,-10,0,10,20,30)
  labels=br
  limi=c(-30,30)
  tsize=16
  osize=12
  legend="Change in Qsp \n(l/s/km2)"
  ocrap<-ggplot(basemap) +
    geom_sf(fill="gray85",color="darkgrey",size=0.5)+
    geom_sf(data=pointsag,aes(fill=mkta,geometry=geometry),alpha=0.4,color="transparent")+
    geom_sf(data=points,aes(col=fillplot,geometry=geometry),alpha=.9,size=0.15,stroke=0,shape=15)+ 
    #geom_sf(data=pointplot,aes(geometry=geometry),color="darkgrey",fill=NA)+ 
    scale_fill_gradientn(
      colors=palet,
      breaks=seq(-1,1,by=0.2), limits=c(-1,1),
      oob = scales::squish,na.value=colNA, name=legend) +
    guides(fill = "none")+
    new_scale_fill()+
    geom_sf(data=pointsInside,aes(geometry=geometry, fill=sign),alpha=.8,size=.6,stroke=0,shape=21,color="black")+
    #geom_sf(fill=NA, color="grey") +
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
  
  ggsave("D:/tilloal/Documents/LFRuns_utils/TrendAnalysis/plots/mapHR_final_flscf8020_sig.jpg", ocrap, width=20, height=20, units=c("cm"),dpi=1000) 
  
  
  #ocrap
  
  return(list(pl1, pl2))
}














### Comparison of two scenarios with mk test at catchment level ----

### Beam plot of change for Europe ----
plotbeam=FascPlotChange(datarSCF,outf,period=c(1950,2020),lims=c(-60, 60),q1=0.25,q2=0.75)

ggsave("D:/tilloal/Documents/LFRuns_utils/Figures/Change_alleurope_faisceau_drought.jpg",plotbeam, width=40, height=15, units=c("cm"),dpi=1500)


### Plot by biogeographic regions ----
biogeo <- read_sf(dsn = paste0(hydroDir,"/eea_3035_biogeo-regions_2016/BiogeoRegions2016_wag84.shp"))
biogeof=fortify(biogeo)
st_geometry(biogeof)<-NULL
biogeoregions=raster( paste0(hydroDir,"/eea_3035_biogeo-regions_2016/Biogeo_rasterized_wsg84.tif"))
Gbiogeoregions=as.data.frame(biogeoregions,xy=T)
biogeomatch=inner_join(biogeof,Gbiogeoregions,by= c("PK_UID"="Biogeo_rasterized_wsg84"))
biogeomatch$latlong=paste(round(biogeomatch$x,4),round(biogeomatch$y,4),sep=" ")
biogeo_rivers=inner_join(biogeomatch,outf, by="latlong")

unikbiog=unique(biogeo_rivers$PK_UID)
i=0
ridgechange=list()
faschange=list()
statregions=c()
for (b in unikbiog){
  i=i+1
  outf9=biogeo_rivers[which(biogeo_rivers$PK_UID==b),]
  rplots=plotHistoDates(outf9,datari,law="GPD",type,period=c(1950),parlist=Paramsfl,valuenames)
  plotv2=FascPlotChange(datar,outf9,period=c(1951,2019))
  faschange[[i]]=plotv2
  ggsave(paste0("D:/tilloal/Documents/LFRuns_utils/Figures/FaisceauChange_",i,".jpg"), plotv2, width=20, height=15, units=c("cm"),dpi=1000)
  ridgechange[[i]]=rplots[[2]]
  statregions=rbind(statregions,rplots[[1]])
  
}

#All results in a single file
plotkeep=faschange[c(1,2,3,5,6,7,8)]

layout_mat<-rbind(c(1,2),
                  c(3,4),
                  c(5,6),
                  c(8,9))
layout_mat

library(gridExtra)
plots<-marrangeGrob(plotkeep, layout_matrix=layout_mat)

ggsave("D:/tilloal/Documents/LFRuns_utils/Figures/FaisceauChange.jpg", plots, width=50, height=60, units=c("cm"),dpi=1000)
write.csv(statregions,file=paste0(hydroDir,"/Flood/stats_by_bgreg.csv"))


### Plot of the shape parameter in every catchment ----
pc=unique(ParamsflH$catchment)
parcat=ParamsflH[match(pc,ParamsflH$catchment),]
parcat$catchment
tsize=16
osize=16
basemap=w2
parplot=inner_join(parcat,UpArea,by=c("catchment"="outl2"))
parplot=full_join(outf,parcat,by = c("outlets"="catchment"))
parpl <- st_as_sf(parplot, coords = c("Var1", "Var2"), crs = 4326)
parpl <- st_transform(parpl, crs = 3035)


palet=c(hcl.colors(9, palette = "BrBG", alpha = NULL, rev = F, fixup = TRUE))
ggplot(basemap) +
  geom_sf(data=parpl,aes(col=epsilonGPD,geometry=geometry),alpha=1,size=0.1,stroke=0,shape=15)+ 
  geom_sf(fill=NA, color="grey") +
  coord_sf(xlim = c(min(nco[,1]),max(nco[,1])), ylim = c(min(nco[,2]),max(nco[,2])))+
  scale_color_gradientn(
    colors=palet, limits=c(-0.5,1),
    oob = scales::squish,na.value="grey", name="shape parameter")   +
  labs(x="Longitude", y = "Latitude")+
  theme(axis.title=element_text(size=tsize),
        panel.background = element_rect(fill = "aliceblue", colour = "grey1"),
        panel.border = element_rect(linetype = "solid", fill = NA, colour="black"),
        legend.title = element_text(size=tsize),
        legend.text = element_text(size=osize),
        panel.grid.major = element_line(colour = "grey70"),
        panel.grid.minor = element_line(colour = "grey90"),
        legend.key = element_rect(fill = "transparent", colour = "transparent"),
        legend.key.size = unit(.8, "cm"))

ggsave("D:/tilloal/Documents/LFRuns_utils/Figures/ShapeParams.jpg", width=30, height=20, units=c("cm"),dpi=1000)

### Plot the 100Y RL at a given time ----
hazard="drought"
chyr=2000
Ychyr=paste0("Y",chyr)

parplot=inner_join(datar,outf,by = c("unikout"="outlets"))
parpl <- st_as_sf(parplot, coords = c("Var1", "Var2"), crs = 4326)
parpl <- st_transform(parpl, crs = 3035)
colplot=match(Ychyr,names(parpl))
names(parpl)[colplot]="Ysel"
palet=c(hcl.colors(9, palette = "YlGnBu", alpha = NULL, rev = T, fixup = TRUE))
max(parpl$Ysel, na.rm=T)

ggplot(basemap) +
  geom_sf(fill="white", color=NA) +
  geom_sf(data=parpl,aes(col=Ysel,geometry=geometry),alpha=1,size=0.25,stroke=0,shape=15)+ 
  geom_sf(fill=NA, color="grey") +
  coord_sf(xlim = c(min(nco[,1]),max(nco[,1])), ylim = c(min(nco[,2]),max(nco[,2])))+
  scale_color_gradientn(
    colors=palet, limits=c(0.1,200),
    breaks=c(0.1,1,10,100),
    trans="log",
    oob = scales::squish,na.value="darkorange", name="Q (m3/s)")   +
  guides(color = guide_colourbar(barwidth = 10, barheight = 1))+
  labs(x="Longitude", y = "Latitude")+
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
  ggtitle(paste0("100 years ", hazard ," Return level - ", chyr))

ggsave("D:/tilloal/Documents/LFRuns_utils/Figures/drought_100yRL_2000.jpg", width=30, height=20, units=c("cm"),dpi=1000)



### Look ar change through different RPs ----
RLsave=c()
for (RP in c(5,10,20,50,100))
{
  print(RP)
  RLs=GPDLargeRLs(Paramsfl,RP)
  
  RLsave=cbind(RLsave,RLs$GPD$returnLevels)
}

RLsave=as.data.frame(RLsave)
Paramplus=cbind(Paramsfl,RLsave)

### Maximum year and work on multi-decadal cycle ----

datar=RLGPDfl
datart=datar[,-c(71,72)]
ymax=apply(datart, 1, which.max)
ymin=apply(datart, 1, which.min)

subdata=datart
datrd=subdata/subdata$Y1950*100-100

yp=c(1950:2019)
maxyears=yp[ymax]
minyears=yp[ymin]
hist(maxyears, breaks=7)
hist(minyears, breaks=7)

datart$maxyear=maxyears
datart$minyear=minyears
datari=inner_join(outf,datar,by = c("outlets"="unikout"))

points <- st_as_sf(datari, coords = c("Var1", "Var2"), crs = 4326)
points <- st_transform(points, crs = 3035)

palet=c(hcl.colors(9, palette = "PRGn", alpha = NULL, rev = T, fixup = TRUE))
ggplot(w2) +
  geom_sf(fill="gray85")+
  geom_sf(fill=NA, color="grey") +
  geom_sf(data=points,aes(col=factor(gt2),geometry=geometry),alpha=1,size=0.25,stroke=0,shape=15)+ 
  coord_sf(xlim = c(min(nco[,1]),max(nco[,1])), ylim = c(min(nco[,2]),max(nco[,2])))+

  scale_color_manual(values=c("orange","green","blue"))  +
  labs(x="Longitude", y = "Latitude")+
  theme(
    panel.background = element_rect(fill = "aliceblue", colour = "grey1"),
    panel.border = element_rect(linetype = "solid", fill = NA, colour="black"),
    legend.position = "right",
    panel.grid.major = element_line(colour = "grey70"),
    panel.grid.minor = element_line(colour = "grey90"),
    legend.key = element_rect(fill = "transparent", colour = "transparent"),
    legend.key.size = unit(.8, "cm"))

ggsave("D:/tilloal/Documents/LFRuns_utils/Figures/maxgroup.jpg", width=30, height=20, units=c("cm"),dpi=1000)


# Link to GMST ------------------------------------------------------------


#extraction of RL100 values per year
subdata=datari[,-c(1:6,77:81)]
datrd=subdata/subdata$Y1950*100-100

#get quantiles of change per year
dfplot= datrd %>% gather("Year","pixel")
dfplot$year=as.numeric(substr(dfplot$Year,2,5))

dfquantiles=aggregate(list(change=dfplot$pixel),
                      by = list(Year=dfplot$year),
                      FUN = function(x) c(mean=mean(x,na.rm=T),q1=quantile(x, 0.025, na.rm=T),q3=quantile(x, 0.975, na.rm=T)))

dfquantiles <- do.call(data.frame, dfquantiles)
names(dfquantiles)=c("Year", "mean","qlow","qhigh")

dfquantiles$Year[which.max(dfquantiles$mean)]

#Read GMST from HadCRUT
GMST=read.csv(file=paste0(hydroDir,"/timeseries/HadCRUT.5.0.1.0.csv"))

#Preindustrial temperature as defined in IPCC
PreIndustrial=mean(GMST$Anomaly..deg.C.[32:61])

#adjusting change to preIndustrial warming
GMST$Anomaly_PreInd=GMST$Anomaly..deg.C.-PreIndustrial

#30 years moving average of temperature
GMST$runningAnomaly=tsEvaNanRunningMean(GMST$Anomaly_PreInd, 30)
plot(GMST$Time,GMST$Anomaly_PreInd,type="o")
abline(h=GMST_sub$runningAnomaly[1])
lines(GMST$Time,GMST$runningAnomaly,col=2, lwd=2)

GMST$Time
tbound=c(1950,2019)
timestamp=seq(tbound[1],tbound[2])
yRL=dfquantiles$Year

GMST_sub=GMST[match(timestamp,GMST$Time),]
GMST_sub$meanchange=dfquantiles$mean
plot(GMST_sub$Time,GMST_sub$Anomaly_PreInd,type="o")
lines(GMST_sub$Time,GMST_sub$runningAnomaly,col=2, lwd=2)


plot(GMST_sub$runningAnomaly,GMST_sub$meanchange,type="o",ylim=c(-10,10))
lines(GMST_sub$runningAnomaly,dfquantiles$qlow)
lines(GMST_sub$runningAnomaly,dfquantiles$qhigh)


#Development part: approximating values to GW times
peakWL=approx(GMST_sub$Time,GMST_sub$runningAnomaly,xout=peakyears)$y
points(jitter(peakWL),pikos$value,pch=16, col="blue")
GWlevels=c(0.5,0.7,1)
WLyear=round(approx(GMST_sub$runningAnomaly,GMST_sub$Time,xout=GWlevels)$y)


# Seasonal analysis at large scale ----------------------------------------

load(file=paste0(hydroDir,"/Flood/params.flood.cal2.Rdata"))
load(file=paste0(hydroDir,"/Flood/peaks.flood.cal2.Rdata"))

#For visualisation purposes, I run the script on one pixel (loop is run on HPC)

#time settings
#max id of the time vector
tidm=25931
timeStamps6h=as.POSIXct(c(0:(tidm*4))*3600*6 -3600, origin="1950-01-03 12:00:00")
tsm=1/0.25
rmv=c(1:(60*tsm))
timeStamps=timeStamps6h[-rmv]
time=unique(as.Date(timeStamps))
time=time[-length(time)]
timeStamps=timeStamps6h[-rmv]


Nsq=42
print(Nsq)
Idstart=as.numeric(Nsq)*100000
nrspace=rspace[Nsq,]
outhybas=outletopen(hydroDir,outletname,nrspace)
outhybas$outlets=seq((Idstart+1),(Idstart+length(outhybas$outlets)))
outhybas$latlong=paste(round(outhybas$Var1,4),round(outhybas$Var2,4),sep=" ")



Nmat=which(!is.na(match(Peaksave$catch,outhybas$outlets)))
Pikday=Peaksave[Nmat,]
unikzero=unique(Pikday$catch)
id=1
cat(paste0("\n",id))
Rpix=Pikday[which(Pikday$catch == unikzero[id]),]

Rpix$d=yday(Rpix$time)
Rpix$year=year(Rpix$time)

#I don't want to look only at the peak but at all drought times

#vector of all drought days
tsEvents=c()
for (ev in 1:length(Rpix$value)){
  event=Rpix[ev,]
  #check that my new bricoled time stamp is valid
  event$time2=timeStamps[event$timeID]
  if (event$time2==event$time){
    seqdates=seq(event$tIDstart,event$tIDend)
    tev=timeStamps[seqdates]
    evn=rep(ev,length(tev))
    tsEv=data.frame(tev,evn)
    tsEvents=rbind(tsEvents,tsEv)
  }else{
    print("times not matching")
  }
}
tsEvents$day=yday(tsEvents$tev)
tsEvents$date=as.Date(tsEvents$tev)
tsEventsD=aggregate(list(ev=tsEvents$evn),
                    by = list(event=tsEvents$evn, date=tsEvents$date),
                    FUN = function(x) c(len=length(x)))
tsEventsD <- do.call(data.frame, tsEventsD)
tsEventsD$day=yday(tsEventsD$date)
tsEventsD$ISev=1
dayvec=data.frame(day=seq(1,366))
nyears=length(unique(year(time)))
tw=30
windowSize=tw*365.25
timestampsD=data.frame(date=time)

tsEventsF=right_join(tsEventsD,timestampsD,by="date")
tsEventsF=tsEventsF[order(tsEventsF$date),]
tsEventsF$ISev[which(is.na(tsEventsF$ISev))]=0

if (length(which(diff(tsEventsF$date)==0))) tsEventsF=tsEventsF[-which(diff(tsEventsF$date)==0),]

series=tsEventsF$day
rseason=RunningSeason(series, timestampsD$date,windowSize, nyears, 50)

rseason$diff=rseason$rnseas-rseason$rnseas[1]
if (length(which(is.na(rseason$diff)))<1){
  rseason$diff[which(rseason$diff>182)]=rseason$diff[which(rseason$diff>182)]-365.25
  rseason$diff[which(rseason$diff<=-182)]=365.25+rseason$diff[which(rseason$diff<=-182)]
}
rseason$catch=rep(unikzero[id],length(rseason$rnseas))

#save example plot
rseason[[2]]
ggsave("flood_seasonality_pix3.jpg", width=30, height=21, units=c("cm"),dpi=1000)


rseasonT=c()
outf=c()
hydroDir<-("D:/tilloal/Documents/LFRuns_utils/data")
outlets="Rnet"
lf=list.files(path = paste0(hydroDir,"/Flood/Peaks/v2"), full.names = TRUE, recursive = TRUE)
lf
for (file in lf){
  load(file)
  fils=sub(".*/", "", file)
  tt=unlist(strsplit(fils, "[_]"))
  Nsq=as.numeric(tt[2])
  out.type=sub("\\_.*", "", fils)
  print(Nsq)
  
  
  
  rspace= read.csv(paste0(hydroDir,"/subspace_efas.csv"))
  rspace=rspace[,-1]
  nrspace=rspace[Nsq,]
  #outletname="outletsv8_hybas07_01min"
  #outletname="outlets_hybas09_01min"
  outletname="efas_rnet_100km_01min"
  
  outhybas=outletopen(hydroDir,outletname,nrspace)
  Idstart=as.numeric(Nsq)*100000
  if (length(outhybas$outlets)>0){
    outhybas$outlets=seq((Idstart+1),(Idstart+length(outhybas$outlets)))
    outhybas$latlong=paste(round(outhybas$Var1,4),round(outhybas$Var2,4),sep=" ")
    # outcut=which(!is.na(match(outhybas$outlets,parlist$catchment)))
    # zebi=unique(parlist$catchment)
    # outhloc=outhybas[outcut,]
    
    outf=rbind(outf,outhybas)
  }
  
  
  rseasonT=rbind(rseasonT,rseasonF)
  
  
}







#comput mean, median, q025 and q975
summarystat=aggregate(list(diff=rseasonT$rnseas),
                      by = list(year=rseasonT$time),
                      FUN = function(x) c(mean=mean(x,na.rm=T),median=median(x,na.rm=T),q025=quantile(x,0.025,na.rm=T),q975=quantile(x,0.975,na.rm=T)))
summarystat=do.call(data.frame, summarystat)

rseasonT$year=year(rseasonT$time)
rseasonP=rseasonT[which(rseasonT$year==2005),]

seasonplot=full_join(rseasonP,outf, by=c("catch"="outlets"))
seasonplot$sq=round(seasonplot$catch/100000)

seasonP <- st_as_sf(seasonplot, coords = c("Var1", "Var2"), crs = 4326)
seasonP <- st_transform(seasonP, crs = 3035)

cord.dec=outf[,c(2,3)]
cord.dec = SpatialPoints(cord.dec, proj4string=CRS("+proj=longlat"))
cord.UTM <- spTransform(cord.dec, CRS("+init=epsg:3035"))
nco=cord.UTM@coords


####################### Seasonal Plots #######################


#Plots for seasonality

## Mean flood date plot
direction_labeller <- function(x){
  ifelse(x %% 30 == 0, rev(c('Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct',"Nov",'Dec'))[(as.integer(x/30) %% 13)], '')
}
my_colors=(c("#2E22EA","#9E3DFB","#F86BE2","#FCCE7B","#C4E416","#4BBA0F","#447D87","#2C24E9"))
x=seq(360,30,-30)
ptn=direction_labeller(seq(360,30,-30))
ptnn=rev(ptn)
paletx=c(hcl.colors(8, palette = "PurpOr", alpha = NULL, rev = TRUE, fixup = TRUE))
tsize=16
osize=16
ggplot(w2) +
  geom_sf(fill="white")+
  geom_sf(data=seasonP,aes(col=rnp1,geometry=geometry),alpha=1,size=0.25,stroke=0,shape=15)+ 
  geom_sf(fill=NA, color="grey") +
  coord_sf(xlim = c(min(nco[,1]),max(nco[,1])), ylim = c(min(nco[,2]),max(nco[,2])))+
  scale_color_gradientn(
    colors=my_colors,
    breaks=seq(15,345,30),
    label=ptn,
    na.value="grey95",
    limits=c(0,366), name="Date")+
  labs(x="Longitude", y = "Latitude")+
  guides(color = guide_colourbar(barwidth = 2, barheight = 15,reverse=T))+
  theme(axis.title=element_text(size=tsize),
        panel.background = element_rect(fill = "aliceblue", colour = "grey1"),
        panel.border = element_rect(linetype = "solid", fill = NA, colour="black"),
        legend.title = element_text(size=tsize),
        legend.text = element_text(size=osize),
        panel.grid.major = element_line(colour = "grey70"),
        panel.grid.minor = element_line(colour = "grey90"),
        legend.key = element_rect(fill = "transparent", colour = "transparent"),
        legend.key.size = unit(.8, "cm"))+
  ggtitle("mean flood date")
ggsave( filename="D:/tilloal/Documents/LFRuns_utils/Figures/mean_flood_date_1965.jpg", width=16.3, height=15, units="cm",dpi=3000)



## Concentration of flood timing plot
paletx=c(hcl.colors(8, palette = "YlOrRd", alpha = NULL, rev = TRUE, fixup = TRUE))

ggplot(w2) +
  geom_sf(fill="white")+
  geom_sf(data=seasonP,aes(col=rncon,geometry=geometry),alpha=1,size=0.25,stroke=0,shape=15)+ 
  geom_sf(fill=NA, color="grey") +
  coord_sf(xlim = c(min(nco[,1]),max(nco[,1])), ylim = c(min(nco[,2]),max(nco[,2])))+
  scale_color_gradientn(
    colors=paletx,
    limits=c(0,1),
    na.value="grey95",oob = scales::squish, name="R")+
  labs(x="Longitude", y = "Latitude")+
  guides(color = guide_colourbar(barwidth = 1, barheight = 10,reverse=F))+
  theme(axis.title=element_text(size=tsize),
        panel.background = element_rect(fill = "aliceblue", colour = "grey1"),
        panel.border = element_rect(linetype = "solid", fill = NA, colour="black"),
        legend.title = element_text(size=tsize),
        legend.text = element_text(size=osize),
        panel.grid.major = element_line(colour = "grey70"),
        panel.grid.minor = element_line(colour = "grey90"),
        legend.key = element_rect(fill = "transparent", colour = "transparent"),
        legend.key.size = unit(.8, "cm"))+
  ggtitle("Concentration around mean flood date")

ggsave( filename="D:/tilloal/Documents/LFRuns_utils/Figures/concentration_flood_season_2005x.jpg", width=16.3, height=15, units="cm",dpi=3000)



## Change in flood timing: 1950-2020
rseason1950=rseasonT[which(rseasonT$year==1965),]
rseason2020=rseasonT[which(rseasonT$year==2005),]

rseason2020$rnseas1950=rseason1950$rnseas
rseason2020$rnnp1950=rseason1950$rnnp
rseason2020$rnp11950=rseason1950$rnp1
rseason2020$rnp21950=rseason1950$rnp2
rseason2020$diff=rseason2020$rnseas-rseason1950$rnseas

rseason2020$diff[which(rseason2020$diff>182)]=rseason2020$diff[which(rseason2020$diff>182)]-365.25
rseason2020$diff[which(rseason2020$diff<=-182)]=365.25+rseason2020$diff[which(rseason2020$diff<=-182)]


rseason2020=full_join(rseason2020,outf, by=c("catch"="outlets"))
rseason2020$sq=round(rseason2020$catch/100000)
rseason2020P <- st_as_sf(rseason2020, coords = c("Var1", "Var2"), crs = 4326)
rseason2020P <- st_transform(rseason2020P, crs = 3035)

#remove differences that are too big for visualisation
rseason2020P$diff2=rseason2020P$diff
rseason2020P$diff2[which(abs(rseason2020P$diff)>60)]=NA
paletx=c(hcl.colors(8, palette = "Spectral", alpha = NULL, rev = TRUE, fixup = TRUE))

ggplot(w2) +
  geom_sf(fill="white")+
  geom_sf(data=rseason2020P,aes(col=diff,geometry=geometry),alpha=1,size=0.25,stroke=0,shape=15)+ 
  geom_sf(fill=NA, color="grey") +
  coord_sf(xlim = c(min(nco[,1]),max(nco[,1])), ylim = c(min(nco[,2]),max(nco[,2])))+
  scale_color_gradientn(
    colors=paletx,
    limits=c(-60,60),
    na.value="grey95",oob = scales::squish, name=" days")+
  labs(x="Longitude", y = "Latitude")+
  theme(axis.title=element_text(size=tsize),
        panel.background = element_rect(fill = "aliceblue", colour = "grey1"),
        panel.border = element_rect(linetype = "solid", fill = NA, colour="black"),
        legend.title = element_text(size=tsize),
        legend.text = element_text(size=osize),
        panel.grid.major = element_line(colour = "grey70"),
        panel.grid.minor = element_line(colour = "grey90"),
        legend.key = element_rect(fill = "transparent", colour = "transparent"),
        legend.key.size = unit(.8, "cm"))+
  ggtitle("change in mean flood timing (2005-1965)")


ggsave( filename="D:/tilloal/Documents/LFRuns_utils/Figures/change_floodTiming_19652005x7.jpg", width=16.3, height=15, units="cm",dpi=3000)


## Aggregate by catchment to decipher a spatial pattern

hist(rseason2020$diff,breaks=365)
mean(rseason2020$diff,na.rm=T)
rseason2020$diff2=rseason2020$diff
rseason2020$diff2[which(abs(rseason2020$diff)>60)]=NA
rseason2020x=inner_join(rseason2020,Catamere07,by=c("latlong"="llcoord"))
rseasonAgg=aggregate(list(Tchange=rseason2020x$diff),
                     by = list(HYBAS_ID=rseason2020x$HYBAS_ID),
                     FUN = function(x) c(med=median(x,na.rm=T),dev=sd(x,na.rm=T),len=length(x)))
rseasonAgg <- do.call(data.frame, rseasonAgg)

rseasonAgg=inner_join(hybasf7,rseasonAgg,by= "HYBAS_ID")

paletx=c(hcl.colors(8, palette = "Spectral", alpha = NULL, rev = TRUE, fixup = TRUE))
ggplot(w2) +
  geom_sf(fill="white")+
  geom_sf(data=rseasonAgg,aes(fill=Tchange.med,geometry=geometry),color="transparent")+ 
  geom_sf(fill=NA, color="gray") +
  coord_sf(xlim = c(min(nco[,1]),max(nco[,1])), ylim = c(min(nco[,2]),max(nco[,2])))+
  scale_fill_gradientn(
    colors=paletx,
    limits=c(-60,60),
    name=" days",
    oob = scales::squish,na.value="gray")   +
  labs(x="Longitude", y = "Latitude")+
  guides(fill = guide_colourbar(barwidth = 2, barheight = 10))+
  theme(axis.title=element_text(size=tsize),
        panel.background = element_rect(fill = "aliceblue", colour = "grey1"),
        panel.border = element_rect(linetype = "solid", fill = NA, colour="black"),
        legend.title = element_text(size=tsize),
        legend.text = element_text(size=osize),
        panel.grid.major = element_line(colour = "grey70"),
        panel.grid.minor = element_line(colour = "grey90"),
        legend.key = element_rect(fill = "transparent", colour = "transparent"),
        legend.key.size = unit(.8, "cm"))

ggsave( filename="D:/tilloal/Documents/LFRuns_utils/Figures/change_floodTiming_19652005_cat.jpg", width=16.3, height=15, units="cm",dpi=3000)


#Change in number of flood season: 1950-2020

rseason2020$chns=rseason2020$rnnp-rseason2020$rnnp1950
rseason2020P <- st_as_sf(rseason2020, coords = c("Var1", "Var2"), crs = 4326)
rseason2020P <- st_transform(rseason2020P, crs = 3035)

paletx=c(hcl.colors(8, palette = "RdBu", alpha = NULL, rev = TRUE, fixup = TRUE))
ggplot(w2) +
  geom_sf(fill="white")+
  geom_sf(data=rseason2020P,aes(col=chns,geometry=geometry),alpha=1,size=0.1,stroke=0,shape=15)+ 
  geom_sf(fill=NA, color="grey") +
  coord_sf(xlim = c(min(nco[,1]),max(nco[,1])), ylim = c(min(nco[,2]),max(nco[,2])))+
  scale_color_gradientn(
    colors=paletx,
    limits=c(-1,1),
    na.value="grey95",oob = scales::squish, name=" change number of flood seasons")+
  labs(x="Longitude", y = "Latitude")+
  theme(axis.title=element_text(size=tsize),
        panel.background = element_rect(fill = "aliceblue", colour = "grey1"),
        panel.border = element_rect(linetype = "solid", fill = NA, colour="black"),
        legend.title = element_text(size=tsize),
        legend.text = element_text(size=osize),
        panel.grid.major = element_line(colour = "grey70"),
        panel.grid.minor = element_line(colour = "grey90"),
        legend.key = element_rect(fill = "transparent", colour = "transparent"),
        legend.key.size = unit(.8, "cm"))+
  ggtitle("change in mean flood seasons (2005-1965)")


my_colors2=rev(c("1"="#8fce00","2"="#4d2aa7","3"="#990000","4"="#4d2aa7"))
ggplot(w2) +
  geom_sf(fill="white")+
  geom_sf(data=seasonP,aes(color=factor(rnnp),geometry=geometry),alpha=1,size=0.3,stroke=0,shape=15)+ 
  geom_sf(fill=NA, color="grey") +
  coord_sf(xlim = c(min(nco[,1]),max(nco[,1])), ylim = c(min(nco[,2]),max(nco[,2])))+
  scale_color_manual(
    values=my_colors2,
    na.value="grey",
    breaks=c("4","3","2","1"),
    name=" # flood \n seasons")+
  guides(colour = guide_legend(override.aes = list(size = 10)))+
  labs(x="Longitude", y = "Latitude")+
  theme(axis.title=element_text(size=tsize),
        panel.background = element_rect(fill = "aliceblue", colour = "grey1"),
        panel.border = element_rect(linetype = "solid", fill = NA, colour="black"),
        legend.title = element_text(size=tsize),
        legend.text = element_text(size=osize),
        panel.grid.major = element_line(colour = "grey70"),
        panel.grid.minor = element_line(colour = "grey90"),
        legend.key = element_rect(fill = "transparent", colour = "transparent"),
        legend.key.size = unit(.8, "cm"))

ggsave( filename="D:/tilloal/Documents/LFRuns_utils/Figures/NfloodSeasons_2005v2.jpg", width=16.3, height=15, units="cm",dpi=1600)

#Peaks in seasons (spring, summer autumn, winter) -------------------------


## different approaches to season clustering ----

seasonP$spring=NA
seasonP$summer=NA
seasonP$autumn=NA
seasonP$winter=NA

seasonP$spring[which(seasonP$rnp1>80 & seasonP$rnp1<172)]=1
seasonP$spring[which(seasonP$rnp2>80 & seasonP$rnp2<172)]=1
seasonP$spring[which(seasonP$rnp3>80 & seasonP$rnp3<172)]=1
seasonP$spring[which(seasonP$rnp4>80 & seasonP$rnp4<172)]=1

seasonP$summer[which(seasonP$rnp1>171 & seasonP$rnp1<266)]=1
seasonP$summer[which(seasonP$rnp2>171 & seasonP$rnp2<266)]=1
seasonP$summer[which(seasonP$rnp3>171 & seasonP$rnp3<266)]=1
seasonP$summer[which(seasonP$rnp4>171 & seasonP$rnp4<266)]=1

seasonP$autumn[which(seasonP$rnp1>265 & seasonP$rnp1<356)]=1
seasonP$autumn[which(seasonP$rnp2>265 & seasonP$rnp2<356)]=1
seasonP$autumn[which(seasonP$rnp3>265 & seasonP$rnp3<356)]=1
seasonP$autumn[which(seasonP$rnp4>265 & seasonP$rnp4<356)]=1

seasonP$winter[which(seasonP$rnp1>355 | seasonP$rnp1<81)]=1
seasonP$winter[which(seasonP$rnp2>355 | seasonP$rnp2<81)]=1
seasonP$winter[which(seasonP$rnp3>355 | seasonP$rnp3<81)]=1
seasonP$winter[which(seasonP$rnp4>355 | seasonP$rnp4<81)]=1


rseason1950=rseasonT[which(rseasonT$year==1965),]

rseason2020=rseasonT[which(rseasonT$year==2005),]
rseason1950$spring=NA
rseason1950$summer=NA
rseason1950$autumn=NA
rseason1950$winter=NA
rseason1950$spring[which(rseason1950$rnp1>80 & rseason1950$rnp1<172)]=rseason1950$rnp1[which(rseason1950$rnp1>80 & rseason1950$rnp1<172)]
rseason1950$spring[which(rseason1950$rnp2>80 & rseason1950$rnp2<172)]=rseason1950$rnp2[which(rseason1950$rnp2>80 & rseason1950$rnp2<172)]

rseason1950$summer[which(rseason1950$rnp1>171 & rseason1950$rnp1<266)]=rseason1950$rnp1[which(rseason1950$rnp1>171 & rseason1950$rnp1<266)]
rseason1950$summer[which(rseason1950$rnp2>171 & rseason1950$rnp2<266)]=rseason1950$rnp2[which(rseason1950$rnp2>171 & rseason1950$rnp2<266)]

rseason1950$autumn[which(rseason1950$rnp1>265 & rseason1950$rnp1<356)]=rseason1950$rnp1[which(rseason1950$rnp1>265 & rseason1950$rnp1<356)]
rseason1950$autumn[which(rseason1950$rnp2>265 & rseason1950$rnp2<356)]=rseason1950$rnp2[which(rseason1950$rnp2>265 & rseason1950$rnp2<356)]

rseason1950$winter[which(rseason1950$rnp1>355 | rseason1950$rnp1<81)]=rseason1950$rnp1[which(rseason1950$rnp1>355 | rseason1950$rnp1<81)]
rseason1950$winter[which(rseason1950$rnp2>355 | rseason1950$rnp2<81)]=rseason1950$rnp2[which(rseason1950$rnp2>355 | rseason1950$rnp2<81)]


rseason2020$spring=NA
rseason2020$summer=NA
rseason2020$autumn=NA
rseason2020$winter=NA
rseason2020$spring[which(rseason2020$rnp1>80 & rseason2020$rnp1<172)]=rseason2020$rnp1[which(rseason2020$rnp1>80 & rseason2020$rnp1<172)]
rseason2020$spring[which(rseason2020$rnp2>80 & rseason2020$rnp2<172)]=rseason2020$rnp2[which(rseason2020$rnp2>80 & rseason2020$rnp2<172)]

rseason2020$summer[which(rseason2020$rnp1>171 & rseason2020$rnp1<266)]=rseason2020$rnp1[which(rseason2020$rnp1>171 & rseason2020$rnp1<266)]
rseason2020$summer[which(rseason2020$rnp2>171 & rseason2020$rnp2<266)]=rseason2020$rnp2[which(rseason2020$rnp2>171 & rseason2020$rnp2<266)]

rseason2020$autumn[which(rseason2020$rnp1>265 & rseason2020$rnp1<356)]=rseason2020$rnp1[which(rseason2020$rnp1>265 & rseason2020$rnp1<356)]
rseason2020$autumn[which(rseason2020$rnp2>265 & rseason2020$rnp2<356)]=rseason2020$rnp2[which(rseason2020$rnp2>265 & rseason2020$rnp2<356)]

rseason2020$winter[which(rseason2020$rnp1>355 | rseason2020$rnp1<81)]=rseason2020$rnp1[which(rseason2020$rnp1>355 | rseason2020$rnp1<81)]
rseason2020$winter[which(rseason2020$rnp2>355 | rseason2020$rnp2<81)]=rseason2020$rnp2[which(rseason2020$rnp2>355 | rseason2020$rnp2<81)]


rseason2020$dwint=rseason2020$winter-rseason1950$winter

rseason2020$dspring[which(rseason2020$dspring>182)]=rseason2020$dspring[which(rseason2020$dspring>182)]-365.25
rseason2020$dspring[which(rseason2020$dspring<=-182)]=365.25+rseason2020$dspring[which(rseason2020$dspring<=-182)]

rseason2020$dspring=rseason2020$spring-rseason1950$spring

rseason2020$dwint[which(rseason2020$dwint>182)]=rseason2020$dwint[which(rseason2020$dwint>182)]-365.25
rseason2020$dwint[which(rseason2020$dwint<=-182)]=365.25+rseason2020$dwint[which(rseason2020$dwint<=-182)]

rseason2020$dsum=rseason2020$summer-rseason1950$summer

rseason2020$dsum[which(rseason2020$dsum>182)]=rseason2020$dsum[which(rseason2020$dsum>182)]-365.25
rseason2020$dsum[which(rseason2020$dsum<=-182)]=365.25+rseason2020$dsum[which(rseason2020$dsum<=-182)]

rseason2020$daut=rseason2020$autumn-rseason1950$autumn

rseason2020$daut[which(rseason2020$daut>182)]=rseason2020$daut[which(rseason2020$daut>182)]-365.25
rseason2020$daut[which(rseason2020$daut<=-182)]=365.25+rseason2020$daut[which(rseason2020$daut<=-182)]


seasonplot=full_join(rseason2020,outf, by=c("catch"="outlets"))
seasonplot$sq=round(seasonplot$catch/100000)
seasonP <- st_as_sf(seasonplot, coords = c("Var1", "Var2"), crs = 4326)
seasonP <- st_transform(seasonP, crs = 3035)

### Change in mean flood timing ----
paletx=c(hcl.colors(8, palette = "Spectral", alpha = NULL, rev = TRUE, fixup = TRUE))
ggplot(w2) +
  geom_sf(fill="white")+
  geom_sf(data=seasonP,aes(col=daut,geometry=geometry),alpha=1,size=0.25,stroke=0,shape=15)+ 
  geom_sf(fill=NA, color="grey") +
  coord_sf(xlim = c(min(nco[,1]),max(nco[,1])), ylim = c(min(nco[,2]),max(nco[,2])))+
  scale_color_gradientn(
    colors=paletx,
    limits=c(-60,60),
    na.value="grey95",oob = scales::squish, name=" days")+
  labs(x="Longitude", y = "Latitude")+
  theme(axis.title=element_text(size=tsize),
        panel.background = element_rect(fill = "aliceblue", colour = "grey1"),
        panel.border = element_rect(linetype = "solid", fill = NA, colour="black"),
        legend.title = element_text(size=tsize),
        legend.text = element_text(size=osize),
        panel.grid.major = element_line(colour = "grey70"),
        panel.grid.minor = element_line(colour = "grey90"),
        legend.key = element_rect(fill = "transparent", colour = "transparent"),
        legend.key.size = unit(.8, "cm"))+
  ggtitle("change in mean flood timing (2005-1965)")


ggsave( filename="D:/tilloal/Documents/LFRuns_utils/Figures/change_floodTiming_19652005x7.jpg", width=16.3, height=15, units="cm",dpi=3000)

### Aggregation to Hybas catchments ----

seasonPx=inner_join(seasonplot,Catamere07,by=c("latlong"="llcoord"))
rseasonAgg=aggregate(list(Tchange=seasonPx$dsum),
                     by = list(HYBAS_ID=seasonPx$HYBAS_ID),
                     FUN = function(x) c(med=mean(x,na.rm=T),dev=sd(x,na.rm=T),len=length(x)))
rseasonAgg <- do.call(data.frame, rseasonAgg)

rseasonAgg=inner_join(hybasf7,rseasonAgg,by= "HYBAS_ID")

paletx=c(hcl.colors(8, palette = "Spectral", alpha = NULL, rev = TRUE, fixup = TRUE))
ggplot(w2) +
  geom_sf(fill="white")+
  geom_sf(data=rseasonAgg,aes(fill=Tchange.med,geometry=geometry),color="transparent")+ 
  geom_sf(fill=NA, color="gray") +
  coord_sf(xlim = c(min(nco[,1]),max(nco[,1])), ylim = c(min(nco[,2]),max(nco[,2])))+
  scale_fill_gradientn(
    colors=paletx,
    limits=c(-60,60),
    name=" days",
    oob = scales::squish,na.value="gray")   +
  labs(x="Longitude", y = "Latitude")+
  guides(fill = guide_colourbar(barwidth = 2, barheight = 10))+
  theme(axis.title=element_text(size=tsize),
        panel.background = element_rect(fill = "aliceblue", colour = "grey1"),
        panel.border = element_rect(linetype = "solid", fill = NA, colour="black"),
        legend.title = element_text(size=tsize),
        legend.text = element_text(size=osize),
        panel.grid.major = element_line(colour = "grey70"),
        panel.grid.minor = element_line(colour = "grey90"),
        legend.key = element_rect(fill = "transparent", colour = "transparent"),
        legend.key.size = unit(.8, "cm"))

ggsave( filename="D:/tilloal/Documents/LFRuns_utils/Figures/change_floodTiming_19652005_cat.jpg", width=16.3, height=15, units="cm",dpi=3000)


mySplot=plotSeasonFlood(seasonP,nco,tsize,osize)
seasons=c("summer","autumn","winter","spring")
library(ggpubr)
saveplot1=ggarrange(mySplot[[1]],NULL, mySplot[[2]],mySplot[[3]],NULL, mySplot[[4]],
                    ncol = 3, nrow = 2,widths = c(1,-0.3,1,1,-0.3,1))

ggsave( filename="D:/tilloal/Documents/LFRuns_utils/Figures/flood_seasonality6.jpg", saveplot1, width=29.3, height=21, units="cm",dpi=2000)

### Histogram of number of pixels experiencing peaks in each season in period 1 and period 2 ----

seqseason=c("spring","summer","autumn","winter")
lseaf=c()
for (i in (1:4)){
  lseas=data.frame(year=c(1965,2005),season=rep(seqseason[i],2),length=c(length(which(!is.na(rseason1950[,11+i])))/length(rseason1950[,11+i]),length(which(!is.na(rseason2020[,11+i])))/length(rseason2020[,11+i])))
  
  lseaf=rbind(lseaf,lseas)
}
ggplot(lseaf, aes(fill=factor(year), x=factor(season), y=length)) +
  # scale_fill_manual(values = cols, name = "Socioeconomic \n Pathways")+
  # labs(title=tit1,x="Warmimg levels (?C)", y = y1)+
  geom_bar(stat = "identity", position="dodge2")+
  geom_vline(xintercept = 1.5,size=1, col="black",lty="dotted")+
  geom_vline(xintercept = 2.5,size=1, col="black",lty="dotted")+
  geom_vline(xintercept = 3.5,size=1, col="black",lty="dotted")+
  # scale_color_manual(values = cols, name = "Socioeconomic \n Pathways",guide="none")+ 
  theme(axis.title=element_text(size=14),
        panel.background = element_rect(fill = "transparent", colour = "grey10"),
        panel.border = element_rect(linetype = "solid", fill = NA, colour="black"),
        legend.title = element_text(size=14),
        legend.text = element_text(size=12),
        panel.grid.major.y = element_line(colour = "grey70"),
        panel.grid.minor.y = element_line(colour = "grey90"),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        legend.key = element_rect(fill = "transparent", colour = "transparent"),
        legend.key.size = unit(.8, "cm"))

