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
plotchangemaps_qsp=function(basemap,catmap,datar,upArea, GHR_riv, HydroRsf, law="GPD",type,period=c(1965,2005),parlist,valuenames,haz){

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
  valcol2=valcol

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
  dts=data[,c(1,21,22,23,88,106)]
  st_geometry(dts) <- NULL
  
  #generate the grid:
  #use UPArea
  
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
  
  europe_cover = data.frame(x=c(-2700,-2700,2000,2000),y=c(-2500,2000,2000,-2500),data=rep("data",4))
  coordinates(europe_cover) =~x+y
  my_big_grid = makegrid(europe_cover,cellsize=20)
  keep_point = rep(0,nrow(my_big_grid))
  relevant_countries = unique(europe_frame$plot_group)
  
  relevant_countries <- relevant_countries[-grep("Rus|Aze|Turk|Ukr|Isra", relevant_countries, ignore.case = TRUE)]

  
  #for the grid, we check which points are contained in any country
  for(country in relevant_countries)  {
    tmp_country = europe_frame[which(europe_frame$plot_group==country),]
    pt_in_country_tmp = point.in.polygon(point.x=my_big_grid[,1], point.y=my_big_grid[,2], pol.x=tmp_country[,"lambert_x"], pol.y=tmp_country[,"lambert_y"], mode.checked=FALSE)
    if(length(which(pt_in_country_tmp==1))>0) {
      keep_point[which(pt_in_country_tmp==1)] = 1
    }
    print(which(country==relevant_countries)/length(relevant_countries))
  }
  
  my_grid = my_big_grid[which(keep_point==1),]
  names(my_grid) = c("lambert_x","lambert_y")
  
  grid_krig = SpatialPointsDataFrame(coords=my_grid[,c(1,2)], data=my_grid[,c(1,2)],proj4string=CRS("+proj=laea +lat_0=55 +lon_0=20 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=km +no_defs"))
  
  spatial_krig = SpatialPointsDataFrame(coords=dts[,c("Var1.x", "Var2.x")], data=dts[,c("fillplot","uparea")],proj4string=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"))
  spatial_krig_trafo = spTransform(spatial_krig,CRS("+proj=laea +lat_0=55 +lon_0=20 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=km +no_defs"))
  
  dts$lambert_x = spatial_krig_trafo@coords[,1]
  dts$lambert_y = spatial_krig_trafo@coords[,2]
  dts$lxy=paste(dts$lambert_x,dts$lambert_y,sep=" ")
  my_big_grid$lxy=paste(my_big_grid$x1,my_big_grid$x2,sep=" ")

  #kriging
  library("automap")
  grid_krig = SpatialPointsDataFrame(coords=my_grid[,c(1,2)], data=my_grid[,c(1,2)],proj4string=CRS("+proj=laea +lat_0=55 +lon_0=20 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=km +no_defs"))
  
  kriging_result <- autoKrige(fillplot~1, input_data=spatial_krig_trafo, new_data = grid_krig)
  # Convert kriging results to a data frame
  kriging_data <- as.data.frame(kriging_result$krige_output)
  
  # Extract the coordinates from the SpatialPointsDataFrame
  coords <- coordinates(kriging_result$krige_output)
  # Combine the coordinates with the kriging data
  kriging_data$lon <- coords[, 1]
  kriging_data$lat <- coords[, 2]
  
  kriging_data$lxy=paste(kriging_data$lambert_x,kriging_data$lambert_y,sep=" ")
  
  # Transform back to longlat projection if necessary
  kriging_data_out = SpatialPointsDataFrame(coords=kriging_data[,c("lon", "lat")], data=kriging_data[,c(3:5)],proj4string=CRS("+proj=laea +lat_0=55 +lon_0=20 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=km +no_defs"))
  kriging_data_longlat <- spTransform(kriging_data_out, CRS("+init=epsg:3035"))
  
  kriging_data_longlat=as.data.frame(kriging_data_out)
  paletx=c(hcl.colors(7, palette = "RdYlBu", alpha = NULL, rev = T, fixup = TRUE))
  # Create a ggplot map of the kriging results
  m3<-kriging_data_longlat
  m4 <- m3 %>%
    # create a new variable from count
    mutate(countv1=cut(var1.pred, breaks=c(-10000,-24,-12,-5, 0, 5, 12, 100000),
                           labels=c("<-24", "-24 - -12", "-12 - -5", "-5 - 0", "0 - 5", "5 - 12", ">12"))) %>%
    # change level order
    mutate(countv1=factor(as.character(countv1), levels=rev(levels(countv1))))
  
  
  basemap=st_transform(basemap, CRS("+proj=laea +lat_0=55 +lon_0=20 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=km +no_defs"))
  pl3<-ggplot(basemap) +
    geom_sf(fill="white")+
    geom_tile(data=m4, aes(x = lon, y = lat, fill = countv1)) +
    #metR::geom_contour_fill(data=kriging_data_longlat, aes(x = lon, y = lat, z = var1.pred),binwidth=2,alpha=1) +
    geom_sf(fill=NA, color="grey") +
    #scale_fill_distiller(palette = "RdYlBu",limits=c(-12,12),oob = scales::squish,guide="coloursteps",direction=1)+
    #scale_fill_gradientn(colors=paletx,limits=c(-12,12),oob = scales::squish) +
    scale_fill_manual(values=paletx, na.value = "grey90")+
    coord_sf(xlim = c(min(nci[,1]),max(nci[,1])), ylim = c(min(nci[,2]),max(nci[,2])))+
    labs(x="Longitude", y = "Latitude")+
    # guides(fill = guide_coloursteps(barwidth = 1, barheight = 10))+
    labs(fill = "Change [%] ") +
    theme(axis.title=element_text(size=tsize),
          panel.background = element_rect(fill = "aliceblue", colour = "grey1"),
          panel.border = element_rect(linetype = "solid", fill = NA, colour="black"),
          legend.title = element_text(size=tsize),
          legend.text = element_text(size=osize),
          legend.position = "right",
          legend.title.position = "top",
          panel.grid.major = element_line(colour = "grey70"),
          panel.grid.minor = element_line(colour = "grey90"),
          legend.key = element_rect(fill = "transparent", colour = "transparent"),
          legend.key.size = unit(.8, "cm"))+
    ggtitle("Change in 10y flood return level (1965-2005)")
  pl3
  
  return(list(pl1, pl2))
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


# Plotting the outputs ----------------------------------------------------

## 1. selection of data to be plotted ----
#Loading saved results in .Rdata

load(file=paste0(hydroDir,"/TSEVA_hybas/outputs/parCat07_Histo_drought_1950_2020.Rdata"))
load(file=paste0(hydroDir,"/TSEVA_hybas/outputs/ResCat07_Histo_drought_1950_2020.Rdata"))

load(file=paste0(hydroDir,"/TSEVA_hybas/outputs/parCat07_SocCF_drought_1950_2020.Rdata"))
load(file=paste0(hydroDir,"/TSEVA_hybas/outputs/ResCat07_SocCF_drought_1950_2020.Rdata"))
length(which(Results$catrest[,3]>0))

### drought ----
RLGPDdr=Results$RetLevGPD
parlistdN=parlist$catchment
RLGPDdr=data.frame(RLGPDdr,unikout=Results$catrest$catlist)
datar=RLGPDdr

#dealing with IR and zero discharge rivers
outltest=Results$catrest
datar=inner_join(datar,outltest,by=c("unikout"="catlist"))
datarD=datar

### flood ----
load(file=paste0(hydroDir,"/TSEVA_hybas/outputs/parCat07_SocCF_flood_1950_2020.Rdata"))
load(file=paste0(hydroDir,"/TSEVA_hybas/outputs/ResCat07_SocCF_flood_1950_2020.Rdata"))
length(which(Results$catrest[,2]>0))
RLGPDfl=Results$RetLevGPD
RLGPDfl=data.frame(RLGPDfl,unikout=Results$catrest$catlist)
datar=RLGPDfl
#dealing with IR and zero discharge rivers
outltest=Results$catrest
datar=inner_join(datar,outltest,by=c("unikout"="catlist"))
datarF=datar


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

#remove negative values (for drought)
for (c in seq(1,(length(datarD[1,])-3))){
  datarD[which(datarD[,c]<0),c]=0
  
}


datarD2=datarD[-which(datarD$IRES==1),]
haz="flood"
datarD=datarD2

Plot.change=plotchangemaps_qsp(basemap,catmap=catmap,datarF, UpArea,GHR_riv,HydroRsf, law="GPD",type="RLchange",period=c(1960,2010),parlist,valuenames,haz=haz)
Plot.change[[2]]



ggsave(paste0("D:/tilloal/Documents/LFRuns_utils/Figures/Hybas07_RLchangeSCF_qsp_",haz,".jpg"),Plot.change[[1]], width=24, height=20, units=c("cm"),dpi=300)

ggsave(paste0("D:/tilloal/Documents/LFRuns_utils/Figures/Hybas07_Reservoir_IRES_.jpg"),Plot.change[[2]], width=24, height=20, units=c("cm"),dpi=300)


haz="flood"
Plot.change=plotchangemaps_qsp(basemap,catmap=catmap,datarF, UpArea,GHR_riv,HydroRsf, law="GPD",type="RLchange",period=c(1965,2005),parlist,valuenames,haz=haz)
Plot.change[[2]]
ggsave(paste0("D:/tilloal/Documents/LFRuns_utils/Figures/Hybas07_RLchange_",haz,".jpg"),Plot.change[[1]], width=24, height=20, units=c("cm"),dpi=300)

library(animation)
yrseq=c(1952:2020)
pllist=list()
  for(i in yrseq) {
    haz="drought"
    print(i)
    Plot.change=plotchangemaps(basemap,catmap=catmap,datarD, law="GPD",type="RLchange",period=c(1951,i),parlist,valuenames,haz=haz)
    pllist=c(pllist,list(Plot.change[[1]]))
  }

saveGIF({
  for (i in 1:(length(pllist))) {
    print(i)
    print(pllist[[i]])
    #Sys.sleep(1)  # pause for 1 second
  }
 }, movie.name = paste0("D:/tilloal/Documents/LFRuns_utils/Figures/Hybas07_RLDroughtchangeF.gif"),ani.width=1000, ami.height=1000, interval = .25)

#Now bivariate change

Plot.cb=plotbichangemaps(basemap,catmap,datarD,datarF,law="GPD",type,period=c(1991,2020),valuenames)

Plot.cb=plotbichangemaps_qsp(basemap,catmap,datarD,datarF,UpArea,law="GPD",type,period=c(1965,2005),valuenames)

Plot.cb
ggsave(paste0("D:/tilloal/Documents/LFRuns_utils/Figures/Hybas07_bichangeSCF_dr-fl_qsp.jpg"),Plot.cb, width=24, height=20, units=c("cm"),dpi=300)



#Now 




RLGPDflSCF=RLGPDfl
datar=RLGPDflSC
#add latitude and longitude to input
datariSCF=inner_join(outf,datarSCF,by = c("outlets"="unikout"))
datar$Y2020=datar$Y2019
unikout=datar$unikout


datar=RLGPDfl
datariH=inner_join(outf,datar,by = c("outlets"="unikout"))
datar$Y2020=datar$Y2019
unikout=datar$unikout


#add latitude and longitude to input
datari=inner_join(outf,datar,by = c("outlets"="unikout"))

datariSCF=inner_join(outf,datar,by = c("outlets"="unikout"))
#join with catchment data
dataricat=inner_join(Catamere07,datari,by = c("llcoord"="latlong"))

#checkpars=inner_join(datari,Paramsfl[which(Paramsfl$Year==2020),],by=c("outlets"="catchment"))

Impdates=seq(1951,2020,by=1)
valuenames=paste0("X",Impdates)

## 2. different plots of the results ----

### Historical distribution of changes ----
totalch=plotHistoDates(outf,datari,law="GPD",type,period=c(1950),parlist=parlist,valuenames)
totalch[[2]]

#trick to also correct
# v1=Paramsfl[which(Paramsfl$Year==2019),]
# v1=v1[match(unique(v1$catchment),v1$catchment),]
# v1$Year=2020
# Paramsfl[which(Paramsfl$Year==2020),]=v1

Plot.change=plotchangemapix(basemap,catmap=cst7,datari, law="GPD",type="RLchange",period=c(1970,2019),Paramsdr,hybasf = hybasf7,haz="drought")

Plot.change[[1]]
ggsave("D:/tilloal/Documents/LFRuns_utils/Figures/RP_change_floodCal3.jpg", width=20, height=15, units=c("cm"),dpi=1500)

Plot.change[[2]]
ggsave("D:/tilloal/Documents/LFRuns_utils/Figures/RL_change_hybas07.jpg", width=20, height=15, units=c("cm"),dpi=1000)

### Changes in RLs with mk test at catchment level ----
Plot.change.sig=plotTrendSipix(basemap,dataricat,period=c(1950,2020),hybasf = hybasf7,valuenames,nco)
Plot.change.sig
ggsave("D:/tilloal/Documents/LFRuns_utils/Figures/RLchange_drought_sig19502020.jpg",Plot.change.sig, width=20, height=15, units=c("cm"),dpi=1500)


### Comparison of two scenarios with mk test at catchment level ----

### Beam plot of change for Europe ----
plotbeam=FascPlotChange(datarD,outf,period=c(1951,2020),lims=c(-100, 100),q1=0.1,q2=0.9)
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
pc=unique(parlist$catchment)
parcat=parlist[match(pc,parlist$catchment),]
tsize=16
osize=16
basemap=w2

parpl=inner_join(catmap,parcat,by = c("outlets"="catchment"))
# parpl <- st_as_sf(parplot, coords = c("Var1", "Var2"), crs = 4326)
# parpl <- st_transform(parpl, crs = 3035)


palet=c(hcl.colors(9, palette = "BrBG", alpha = NULL, rev = F, fixup = TRUE))
ggplot(basemap) +
  geom_sf(data=parpl,aes(fill=epsilonGPD,geometry=geometry))+ 
  geom_sf(fill=NA, color="grey") +
  coord_sf(xlim = c(min(nco[,1]),max(nco[,1])), ylim = c(min(nco[,2]),max(nco[,2])))+
  scale_fill_gradientn(
    colors=palet, limits=c(-1,0),
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
  scale_fill_gradientn(
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
  ggtitle(paste0("100 years ", haz ," Return level - ", chyr))

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

