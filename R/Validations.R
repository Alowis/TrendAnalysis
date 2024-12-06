# ==============================================================================
# Title: Creation of the validation metrics and plots of the HERA dataset, linked to the article:  
# HERA: a high-resolution pan-European hydrological reanalysis (1951-2020)
# Author: Alois Tilloy - Joint Research Centre - Unit C6 
# Date: 2024 - 07 - 25 
# Description:
#   This script allows to generate plots  and assess the skills of the HERA dataset against
#   observed discharge data from 2448 river gauges across Europe
# ==============================================================================

source("~/06_Floodrivers/DataPaper/Code/HERA/functions.R")

main_path = 'D:/tilloal/Documents/06_Floodrivers/'
valid_path = paste0(main_path,'DataPaper/')
dis_path<-paste0(main_path,'dis/calibrated/filtered/Histo/')
setwd(valid_path)

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


#Function for peaks

library(dplyr)
library(tidyr)
library(lubridate)
library(zoo)
library(RtsEva)

#this function will be used to validate how well I reproduce extreme events
#Then Run TSEVA and compare trend sign and magnitude

obs=dis_data_kge2$Q.y
sim=dis_data_kge2$Q.x
quantile=0.1
tail="low"
match_extremes <- function( sim, obs, quantile,tail) {
  # Check if inputs are numeric
  if (!is.numeric(sim) | !is.numeric(obs)) {
    stop("Both sim and obs must be numeric")
  }
  
  # Check if quantile is in (0, 1)
  if (!(0 < quantile & quantile < 1)) {
    stop("Quantile must be in (0, 1)")
  }
  
  # Calculate observed extremes using POT
  thrsdt <- quantile(obs, quantile,na.rm=T)


  obs[which(is.na(obs))]=-9999
  minEventsPerYear=1
  
  if(tail=="high") {
    pks <- pracma::findpeaks(obs,minpeakdistance = 7, minpeakheight = thrsdt)
  }
  if(tail=="low") {
    pks <- declustpeaks(data = obs ,minpeakdistance = 30 ,minrundistance = 30, qt=thrsdt)

  }
  


  # Match simulated values with observed events
  max_in_window <- numeric(length(pks[,1]))
  tmax_in_window <- numeric(length(pks[,1]))
  
  for (i in seq_along(pks[,1])) {
    start_time <- pks[i,2]-1
    end_time <- pks[i,4]+1
    
    snippet <- sim[c(start_time : end_time)]
    
    if (length(snippet) > 0) {
      max_in_window[i] <- max(snippet)
      tmax_in_window[i] <- which.max(snippet)+start_time-1
    } else {
      max_in_window[i] <- NA
      tmax_in_window[i] <- NA
    }
  }
  
  # Add matched values to the data frame
  res_df=as.data.frame(pks)
  names(res_df)=c("observed","peakobs","start","end")
  res_df$model <- max_in_window
  res_df$peakmodel <- tmax_in_window
  
  # Remove rows with NA
  res_df <- res_df[!is.na(res_df$model), ]
  
  
  # Calculate additional metrics
  res_df$diff <- res_df$model - res_df$observed
  res_df$error <- (res_df$diff)
  res_df$error_norm <- (res_df$diff / res_df$observed)
  res_df$tdiff <- (res_df$peakmodel- res_df$peakobs)
  
  return(res_df)
}

# Data generation -----------------------------
#Load all Q_EFAS (HERA simulated discharge) and create a single file

# Part 1: Load required datasets -------------------------------------------

# observed streamflow
Q_data <- read.csv(paste0('out/Q_19502020.csv'), header = F)  # CSVs with observations
max(time)
Q_data=Q_data[-1,-1]
time=as.Date(as.numeric(Q_data$V2)-1,origin="0000-01-01")[-c(1:4)]

### loading stations' initial database------
Q_stations <- st_read(paste0(valid_path,'GIS/Europe_daily_combined_2022_v2_WGS84_NUTS.shp')) 

# simulated streamflow
Q_sim=read.csv(file="out/HERA_Val2_19502020.csv")
HERA_cordata<-read.csv(file="out/HERA_CorStat_19502020.csv")

replax=match(as.numeric(HERA_cordata[1,]),as.numeric(Q_sim[1,]))[-1]
Q_sim[,replax]=HERA_cordata[,-1]

date2=seq(as.Date("1950-01-04"),as.Date("2020-12-31"),by="days")

Station_data_IDs <- as.vector(t(Q_sim[1, -1]))
Station_obs_IDs <- as.numeric(as.vector(t(Q_data[1, -1])))

# Validaton stations
ValidSf=read.csv(file="Stations/Stations_ValidationF.csv")[,-1]
ValidSY=ValidSf[which(ValidSf$removal!="YES"),]
ValidSY=ValidSY[,c(1:16)]

#ValidSY is the final set of stations used
vEFAS=ValidSY[which(ValidSY$csource=="EFAS"),]
length(which(ValidSY$Rlen<(30*365)))

#extract station with more than 30 of record
ValidSYl=ValidSY[which(ValidSY$Rlen>(30*365)),]


#remove simulated and observed streamflow from useless stations

rmv0=which(!is.na(match(Station_data_IDs,ValidSYl$V1)))
rmv1=which(!is.na(match(Station_obs_IDs,ValidSYl$V1)))

Station_data_IDs=Station_data_IDs[rmv0]
Station_obs_IDs=Station_obs_IDs[rmv1]

head(Q_sim)
Q_sim_val=Q_sim[,c(1,rmv0+1)]
Q_obs_val=(Q_data[,c(1,rmv1+1)])
tki=Q_obs_val[c(2:4),]
Q_obs_val=Q_obs_val[-c(2:4),]

Q_obs_val<- data.frame(lapply(Q_obs_val, as.numeric))

rmvt=which(year(time)==1950)


Q_sim_val=Q_sim_val[-(rmvt+1),]
Q_obs_val=Q_obs_val[-(rmvt+1),]

timeo=as.Date(as.numeric(Q_obs_val$V2)-1,origin="0000-01-01")[-1]
max(timeo)
time2=time[-(rmvt)]
## 2.2 Number of years with data ------------

lR=c()
for (id in 1:length(Station_data_IDs)){
  Qs=Q_obs_val[-1,id+1]
  xl=length(Qs)
  l=length(which(!is.na(Qs)))
  lR=c(lR,l)
}
plot(lR[order(lR)])
RecordLen=data.frame(Station_data_IDs,lR)




## First performance check (Thomas' metrics)
library(hydroGOF)
quantile=0.95
kge_h1=c()
kge_h2=c()
nse_h1=c()
nse_h2=c()
extreme_h=c()
extreme_l=c()
cor_h=c()
cor_l=c()
for (s in 1:length(Station_data_IDs)){
  id_HERA=s
  Station=Station_data_IDs[s]
  # Station="6139261"
  # s=which(Station_data_IDs==Station)
  print(s)
  id_obs=match(Station,Station_obs_IDs)
  Station_obs_IDs[id_obs]
  uparea=ValidSYl$upa[which(ValidSYl$V1==Station)]
  woo=ValidSYl[which(ValidSYl$V1==Station),]
  
  HERA_loc=Q_sim_val[,s+1]
  # cat(HERA_loc[1])
  HERA_loc=HERA_loc[-1]
  length(which(!is.na(HERA_loc)))
  Q_loc=as.numeric(Q_obs_val[,id_obs+1])
  # cat(Q_loc[1])
  Q_loc=Q_loc[-1]
  length(which(!is.na(Q_loc)))
  
  # Qmerde=log(Q_loc+eps)
  # if (length(which(is.na(Qmerde)))>0){
  #   break
  # }
  
  qH95=quantile(HERA_loc,.95,na.rm=T)
  qO95=quantile(Q_loc,.95,na.rm=T)
  
  qH10=quantile(HERA_loc,.10,na.rm=T)
  qO10=quantile(Q_loc,.10,na.rm=T)
  
  Qx=Q_loc[which(Q_loc>=qO10)]
  Qxh=HERA_loc[which(Q_loc>=qO10)]
  # plot(HERA_loc,Q_loc)
  eps=0.01*mean(Q_loc,na.rm=T)
  kge_hera=KGE(HERA_loc,Q_loc, na.rm=TRUE, method="2012",out.type="full")
  kge_inv_hera=KGE((1/(HERA_loc+eps)),(1/(Q_loc+eps)), na.rm=TRUE, method="2012",out.type="full")
  nse_hera=NSE(HERA_loc,Q_loc, na.rm=TRUE)
  nse_log_hera=NSE(log(HERA_loc+eps),log(Q_loc+eps), na.rm=TRUE)
  r_hq=cor(HERA_loc[which(Q_loc>qO95)],Q_loc[which(Q_loc>qO95)])
  r_lq=cor(HERA_loc[which(Q_loc<qO10)],Q_loc[which(Q_loc<qO10)])
  ex_ev=match_extremes(HERA_loc,Q_loc,quantile=0.95,tail="high")
  R1_hnorm=ex_ev$error_norm[1]
  R3_hnorm=mean(ex_ev$error_norm[c(1:3)])
  error_hnorm=mean(ex_ev$error_norm)
  
  ex_ev=match_extremes(HERA_loc,Q_loc,quantile=0.1,tail="low")
  R1_lnorm=ex_ev$error_norm[1]
  R3_lnorm=mean(ex_ev$error_norm[c(1:3)])
  error_lnorm=mean(ex_ev$error_norm)
  
  
  #print(c(Station,nse_log_hera))
  if (is.na(nse_log_hera)){
    break
  }
  kge_h1=rbind(kge_h1,c(Station,kge_hera$KGE.value,kge_hera$KGE.elements))
  kge_h2=rbind(kge_h2,c(Station,kge_inv_hera$KGE.value,kge_inv_hera$KGE.elements))
  nse_h1=rbind(nse_h1,c(Station,nse_hera))
  nse_h2=c(nse_h2,nse_log_hera)
  cor_h=c(cor_h,r_hq)
  cor_l=c(cor_l,r_lq)
  extreme_h=rbind(extreme_h,c(R1_hnorm,R3_hnorm,error_hnorm))
  extreme_l=rbind(extreme_l,c(R1_lnorm,R3_lnorm,error_lnorm))
}

median(extreme_h[,3],na.rm=T)
median(kge_h2[,2])


#load past KGE estimates

kgefile="out/EFAS_obs_kgeAY2.csv"

kge=SpatialSkillPlot(ValidSYl,"kge",kgefile)[[1]]
kge_h1=as.data.frame(kge_h1)
kge_h2=as.data.frame(kge_h2)
kge_h1$kgel=kge_h2$V2
nse_h1=as.data.frame(nse_h1)
kge_h1$nse=nse_h1$V2
ptn=inner_join(kge,kge_h1,by="V1")

hist(kge_h2$V2,xlim=c(-1,1),breaks=2000)
median(kge_h2$V2)
#Now some performance plots

cord.dec=ptn[,c(1,2)]
cord.dec = SpatialPoints(cord.dec, proj4string=CRS("+proj=longlat"))
cord.UTM <- spTransform(cord.dec, CRS("+init=epsg:3035"))
nco=cord.UTM@coords
world <- ne_countries(scale = "medium", returnclass = "sf")
Europe <- world[which(world$continent == "Europe"),]
e2=st_transform(Europe,  crs=3035)
w2=st_transform(world,  crs=3035)
tsize=12
osize=12
basemap=w2

points=ptn
points$kgem=(points$V2+points$kgel)/2
metric="nse"

#For correlation
if (metric=="r"){
  palet=c(hcl.colors(9, palette = "YlGnBu", alpha = NULL, rev = T, fixup = TRUE))
  limi=c(0,1)
  var="r"
}
#For Bias
if (metric=="b"){
  palet=c(hcl.colors(9, palette = "RdBu", alpha = NULL, rev = F, fixup = TRUE))
  limi=c(0,2)
  var="b"
}
#For Variability
if (metric=="nse"){
  palet=c(hcl.colors(9, palette = "RdBu", alpha = NULL, rev = F, fixup = TRUE))
  palet=c(hcl.colors(9, palette = "RdYlGn", alpha = NULL, rev = F, fixup = TRUE))
  #palet = c('#f61f0f','#f4cccc','#f6ebeb','#ffffff','#e4edf5','#91bfdb','#1c87d8')
  limi=c(0,1)
  var="nse"
}

if (metric=="kge"){
  colorz = c('#d73027','orange','#fee090','lightblue','royalblue',"darkblue")
  kgelabs=c("< -0.41","-0.41 - 0","0 - 0.2", "0.2 - 0.5","0.5 - 0.75", ">0.75")
  points$skill=as.numeric(points$kgel)
  points$kgecode=0
  points$kgecode[which(is.na(points$skill))]=NA
  points$kgecode[which(points$skill>(-0.41) & points$skill<=0)]=1
  points$kgecode[which(points$skill>0 & points$skill<=0.2)]=2
  points$kgecode[which(points$skill>0.2 & points$skill<=0.5)]=3
  points$kgecode[which(points$skill>0.5 & points$skill<=0.75)]=4
  points$kgecode[which(points$skill>0.75)]=5
  var="kge"
}


parpl <- st_as_sf(points, coords = c("Var1", "Var2"), crs = 4326)
parpl <- st_transform(parpl, crs = 3035)

if (metric =="kge"){
  tsize=size=12
  map=ggplot(basemap) +
    geom_sf(fill="gray95", color=NA) +
    geom_sf(data=parpl,aes(geometry=geometry,colour=factor(kgecode),size=upa),alpha=.7,stroke=0,shape=16)+ 
    geom_sf(data=parpl,aes(geometry=geometry,size=upa),col="grey20",alpha=1,stroke=0.05,shape=1)+ 
    #geom_sf(data=parpef,aes(geometry=geometry),col="black",alpha=1,size=1.4,stroke=0.1,shape=8)+ 
    geom_sf(fill=NA, color="grey20") +
    scale_x_continuous(breaks=seq(-30,40, by=5)) +
    scale_size(range = c(1, 3), trans="sqrt")+
    labs(x="Longitude", y = "Latitude")+
    coord_sf(xlim = c(min(nco[,1]),max(nco[,1])), ylim = c(min(nco[,2]),max(nco[,2])))+
    scale_color_manual(values = colorz, labels= kgelabs, name="KGE'")   +
    labs(x="Longitude", y = "Latitude") +
    guides(colour = guide_legend(override.aes = list(size = 8,shape=15), reverse=T), size = "none")+
    theme(axis.title=element_text(size=tsize),
          panel.background = element_rect(fill = "aliceblue", colour = "grey1"),
          panel.border = element_rect(linetype = "solid", fill = NA, colour="black"),
          legend.title = element_text(size=tsize),
          legend.text = element_text(size=osize),
          legend.position = "right",
          panel.grid.major = element_line(colour = "grey85",linetype="dashed"),
          panel.grid.minor = element_line(colour = "grey90"),
          legend.key = element_rect(fill = "transparent", colour = "transparent"),
          legend.key.size = unit(.8, "cm"))
  
  map
  ggsave("Plots/KGE_low.jpg", map,width=18, height=20, units=c("cm"),dpi=1000)
}else{
  print(metric)
  #plot for correlation even if called skill
  points$skill=points$nse
  
  parpl <- st_as_sf(points, coords = c("Var1", "Var2"), crs = 4326)
  parpl <- st_transform(parpl, crs = 3035)
  map=ggplot(basemap) +
    geom_sf(fill="gray95", color=NA) +
    geom_sf(data=parpl,aes(geometry=geometry,fill=skill,size=upa),color="transparent",alpha=.9,shape=21,stroke=0)+ 
    geom_sf(data=parpl,aes(geometry=geometry,size=upa),col="grey20",alpha=1,stroke=0.05,shape=1)+ 
    geom_sf(fill=NA, color="grey20") +
    scale_x_continuous(breaks=seq(-30,40, by=5)) +
    scale_size(range = c(1, 3), trans="sqrt")+
    scale_fill_gradientn(
      colors=palet,n.breaks=5,oob = scales::squish, limits=limi, name=metric) +
    coord_sf(xlim = c(min(nco[,1]),max(nco[,1])), ylim = c(min(nco[,2]),max(nco[,2])))+
    labs(x="Longitude", y = "Latitude") +
    guides(fill = guide_colourbar(barwidth = 0.5, barheight = 12,reverse=F),size="none")+
    theme(axis.title=element_text(size=tsize),
          panel.background = element_rect(fill = "aliceblue", colour = "grey1"),
          panel.border = element_rect(linetype = "solid", fill = NA, colour="black"),
          legend.title = element_text(size=tsize),
          legend.text = element_text(size=osize),
          legend.position = "right",
          panel.grid.major = element_line(colour = "grey85",linetype="dashed"),
          panel.grid.minor = element_line(colour = "grey90"),
          legend.key = element_rect(fill = "transparent", colour = "transparent"),
          legend.key.size = unit(.8, "cm"))
  
  
  #ggsave("Plots/NSE_v.jpg", map,width=18, height=20, units=c("cm"),dpi=1000)
  
  map
}


#Extreme events metrics

if (metric=="exH"){
  palet1=c(hcl.colors(9, palette = "RdBu", alpha = NULL, rev = F, fixup = TRUE))
  palet2=c(hcl.colors(9, palette = "RdYlGn", alpha = NULL, rev = F, fixup = TRUE))
  #palet = c('#f61f0f','#f4cccc','#f6ebeb','#ffffff','#e4edf5','#91bfdb','#1c87d8')
  limi=c(0,1)
  
  
  #map for bias on peaks other threshold
  parpl$bpeak=extreme_h[,3]
  map1=ggplot(basemap) +
    geom_sf(fill="gray95", color=NA) +
    geom_sf(data=parpl,aes(geometry=geometry,fill=bpeak,size=upa),color="transparent",alpha=.9,shape=21,stroke=0)+ 
    geom_sf(data=parpl,aes(geometry=geometry,size=upa),col="grey20",alpha=1,stroke=0.05,shape=1)+ 
    geom_sf(fill=NA, color="grey20") +
    scale_x_continuous(breaks=seq(-30,40, by=5)) +
    scale_size(range = c(1, 3), trans="sqrt")+
    scale_fill_gradientn(
      colors=palet1,n.breaks=5,oob = scales::squish, limits=limi, name=metric) +
    coord_sf(xlim = c(min(nco[,1]),max(nco[,1])), ylim = c(min(nco[,2]),max(nco[,2])))+
    labs(x="Longitude", y = "Latitude") +
    guides(fill = guide_colourbar(barwidth = 0.5, barheight = 12,reverse=F),size="none")+
    theme(axis.title=element_text(size=tsize),
          panel.background = element_rect(fill = "aliceblue", colour = "grey1"),
          panel.border = element_rect(linetype = "solid", fill = NA, colour="black"),
          legend.title = element_text(size=tsize),
          legend.text = element_text(size=osize),
          legend.position = "right",
          panel.grid.major = element_line(colour = "grey85",linetype="dashed"),
          panel.grid.minor = element_line(colour = "grey90"),
          legend.key = element_rect(fill = "transparent", colour = "transparent"),
          legend.key.size = unit(.8, "cm"))
  
  
  #ggsave("Plots/NSE_v.jpg", map,width=18, height=20, units=c("cm"),dpi=1000)
  
  map1
  
  #map for correlation of extremes
  parpl$cpeak=cor_h
  map2=ggplot(basemap) +
    geom_sf(fill="gray95", color=NA) +
    geom_sf(data=parpl,aes(geometry=geometry,fill=cpeak,size=upa),color="transparent",alpha=.9,shape=21,stroke=0)+ 
    geom_sf(data=parpl,aes(geometry=geometry,size=upa),col="grey20",alpha=1,stroke=0.05,shape=1)+ 
    geom_sf(fill=NA, color="grey20") +
    scale_x_continuous(breaks=seq(-30,40, by=5)) +
    scale_size(range = c(1, 3), trans="sqrt")+
    scale_fill_gradientn(
      colors=palet2,n.breaks=5,oob = scales::squish, limits=limi, name=metric) +
    coord_sf(xlim = c(min(nco[,1]),max(nco[,1])), ylim = c(min(nco[,2]),max(nco[,2])))+
    labs(x="Longitude", y = "Latitude") +
    guides(fill = guide_colourbar(barwidth = 0.5, barheight = 12,reverse=F),size="none")+
    theme(axis.title=element_text(size=tsize),
          panel.background = element_rect(fill = "aliceblue", colour = "grey1"),
          panel.border = element_rect(linetype = "solid", fill = NA, colour="black"),
          legend.title = element_text(size=tsize),
          legend.text = element_text(size=osize),
          legend.position = "right",
          panel.grid.major = element_line(colour = "grey85",linetype="dashed"),
          panel.grid.minor = element_line(colour = "grey90"),
          legend.key = element_rect(fill = "transparent", colour = "transparent"),
          legend.key.size = unit(.8, "cm"))
  
  
  ggsave("Plots/corr_hpeaks.jpg", map2,width=18, height=20, units=c("cm"),dpi=1000)
  
  map2

}
#Now relationship between performances and upstream area


#For now it is KGE but it can be other metrics, such as the lowKGE
plot(points$upa,points$V2,log="x",ylim=c(-1,1))
abline(h=mean(points$V2),col=2,lwd=2)

library(bestNormalize)
library(quantreg)

BN_obj=bestNormalize(points$upa)
normalupa=orderNorm(points$upa)
normalkge <- orderNorm(points$V2)


sqrt(mean(leffect$residuals^2))

dt1=data.frame(lupa=log(points$upa),norkge=normalkge$x.t)


cor(dt1)

ggplot(dt1,aes(x=lupa,y=norkge))+
  geom_point() +
  stat_smooth(method = "lm")



# Generate a sequence of predicted values

leq2 <- coef(rq(norkge~lupa,data=dt1,tau=0.9))
leq1 <- coef(rq(norkge~lupa,data=dt1,tau=0.1))

leffect=lm(norkge~lupa,data=dt1)

x_seq <- seq(min(dt1$lupa), max(dt1$lupa), length.out = 100)

y_pred <- predict(leffect, newdata = data.frame(lupa = x_seq))

q1_pred <- cbind(1,x_seq)%*%leq1
q2_pred <- cbind(1,x_seq)%*%leq2

plot(dt1$lupa,normalkge$x.t)
abline(h=mean(normalkge$x.t),col=2,lwd=2)
lines(x_seq,y_pred)
lines(x_seq,q_pred)
#Create a data frame of the predicted values
df_pred <- data.frame(x = x_seq, y = y_pred)

x_real=exp(x_seq)

y_real=predict(normalkge, newdata = y_pred, inverse = TRUE)

q1_real=predict(normalkge, newdata = q1_pred, inverse = TRUE)
q2_real=predict(normalkge, newdata = q2_pred, inverse = TRUE)


plot(points$upa,points$V2,ylim=c(-5,1))
lines(x_real,y_real,lwd=2,col=2)
lines(x_real,q1_real,lwd=2,col=3)
lines(x_real,q2_real,lwd=2,col=3)
abline(h=median(normalkge$x),col=2,lwd=2)

finalweight=y_real/median(points$V2)-1
finalw1=q1_real/median(points$V2)-1
finalw2=q2_real/median(points$V2)-1


plot(x_real,finalw2,log="x")
plot(BN_obj, leg_loc = "bottomright")




points$nkge=normalkge$x.t

#################################################################################################
#TSEVA on observed data, for comparisaon with modelled values
#################################################################################################

#Experimental TSEVA function

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
  if(is.na(md))md=1
  print(md)
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





# #######################  Arguments importation #############################

haz = "drought"
outlets = "RNetwork"
season= "nonfrost"
tail="low"

workDir = "D:/tilloal/Documents/LFRuns_utils/data/"
setwd(workDir)

rspace= read.csv(paste0(workDir,"/subspace_efas.csv"))
rspace=rspace[,-1]
nrspace=rspace[Nsq,]
print(nrspace)


# Pre-loaded results -----------


station_ido=data.frame(id=Station_obs_IDs)
station_ids=data.frame(id=Station_data_IDs)
station_ido=station_ido$id[order(station_ido$id)]
station_ids=station_ids$id[order(station_ids$id)]

#Set data directory
hydroDir<-("D:/tilloal/Documents/LFRuns_utils/data")
#rm(outf)
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

unikout=outhybas$outlets
ValidSYl$latlong=paste(round(ValidSYl$Var1,4),round(ValidSYl$Var2,4),sep=" ")

ValidSYl_out=inner_join(ValidSYl,outf,by="latlong")

dirout=paste0(workDir,"TSEVA/out/",haz,"/",season)
fex2=file.exists(dirout)
if (fex2==TRUE){
  print("output folder exists")
}else{
  print("output folder does not exist")
}
######################################################################################

txx=time2


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

ThDir<-paste0(workDir,"TrendAnalysis/trenTH")
THX=c()
for (Nsq in 1:88){
  if (file.exists(paste0(ThDir,"/trenTH_Histo_",tail,"_",Nsq,".csv"))){
    print(Nsq)
    TH1=read.csv(paste0(ThDir,"/trenTH_Histo_",tail,"_",Nsq,".csv"))
    TH2=read.csv(paste0(ThDir,"/trenTH_SCF_",tail,"_",Nsq,".csv"))
    TH3=inner_join(TH1,TH2,by="cid")
    TH3$cid=TH3$cid-Nsq*10000
    TH3$cid=TH3$cid+Nsq*100000
    TH3=inner_join(TH3,ValidSYl_out,by=c("cid"="outl2"))
    THX=rbind(THX,TH3)
  }
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

#average aggregation low flows
ldt=30

RetPerGPD=c()
RetPerGEV=c()
RetLevGEV=c()
RetLevGPD=c()
parlist=c()
peaklist=c()
catlist=c()
IRES=c()

dset="sim"
for (s in 1:length(ValidSYl_out$outl2)){
  id_HERA=s
  Station=ValidSYl_out$V1[s]
  # Station="6139261"
  # s=which(Station_data_IDs==Station)
  print(s)
  id_obs=match(Station,station_ido)
  id_sim=match(Station,Station_data_IDs)
  station_ido[id_obs]
  Station_data_IDs[id_sim]
  uparea=ValidSYl$upa[which(ValidSYl$V1==Station)]
  woo=ValidSYl[which(ValidSYl$V1==Station),]

  #Select the first station
  if (dset=="obs"){
    Q_loc=as.numeric(Q_obs_val[,id_obs+1])
    Q=Q_loc[-1]
  }else if (dset=="sim"){
    Q_loc=as.numeric(Q_sim_val[,id_sim+1])
    Q=Q_loc[-1]
  }
  
  txx=time2
  start_time <- Sys.time()
  catch=as.numeric(ValidSYl_out$outl2[s])
  timeStamps=txx
  thresh=thresh_vec[which(thresh_vec$cid==catch),]
  thresh=thresh$th
  frosttime=NA
  series=data.frame(txx,Q)
  
  names(series)=c("date","Qs")
  rmv=which(as.integer(format(series$date, "%Y"))==1950)
  if (length(rmv)>0){
    series=series[-rmv,]
  }
  #remove first day to avoid bug
  series=series[-1,]
  if (haz=="drought"){
    trans="rev"
    #seasonal split
    catmat=Catf7[which(Catf7$outlets==catch),]
    Tcatmat=mycat[which(mycat$HYBAS_ID==catmat$HYBAS_ID),]
    Tcatchment=which(colnames(frostcat)==Tcatmat$pointid)
    
    intermit=interid(series,trans,WindowSize=ldt)
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
  if(length(na.omit(series$dis))>1 & interflag<1 & nv>15){
    
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
    Nonstat<-TsEvaNs(timeAndSeries, timeWindow, transfType='trendPeaks',ciPercentile = ciPercentile, minPeakDistanceInDays = minPeakDistanceInDays,lowdt=ldt,trans=trans,tail = tail,TrendTh = thresh)
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
      #####
      RLgev=RLevs100$ReturnLevels[2]
      RLgpd=RLevs100$ReturnLevels[3]
      ERgev=RLevs100$ReturnLevels[4]
      ERgpd=RLevs100$ReturnLevels[5]
      
      nRPgev=nRPgpd=10
      params=c()
      for (t in 2:length(Impdates)){
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
    cat(paste0("\n No values in this pixel ",Station," \n or intermittent river (flag = ",interflag,")"))
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
      oops=data.frame(RP=rep(NA,length(year(Impdates))))
      #print(oops)
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
    #parlist=as.data.frame(mapply(c,parlist,params))
    peaklist=as.data.frame((rbind(peaklist,pikos)))
    
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
  
}

# parlist=as.data.frame(parlist)
# peaklist=as.data.frame(peaklist)

Results_sim=list(parameters=parlist,RetLevGEV=RetLevGEV,RetLevGPD=RetLevGPD,RetPerGEV=RetPerGEV,RetPerGPD=RetPerGPD,Peaks=peaklist,catrest=data.frame(catlist,IRES))


save(Results_sim, file=paste0(hydroDir,"/TrendAnalysis/ResValid_",haz,"_",dset,"_1951_2020_30d.Rdata"))

haz="drought"
dset="sim"
load(file=paste0(hydroDir,"/TrendAnalysis/ResValid_",haz,"_",dset,"_1951_2020_v1.Rdata"))
dset="obs"
load(file=paste0(hydroDir,"/TrendAnalysis/ResValid_",haz,"_",dset,"_1951_2020_v1.Rdata"))
#save(parlist,file=paste0(workDir,"TSEVA/out/",sce,"/",haz,"/",season,"/parCat6h_",outlets,haz,"_",Nsq,Quarter,"_1950_2020.Rdata"))

#save(Results, file=paste0(workDir,"TSEVA/out/",runid,"/",haz,"/",season,"/ResCat6h_",outlets,haz,"_Arno_1951_2020_v2.Rdata"))


#Ok now I plot the results

#If I am in drought mode I hav to put a minus sign

if (haz=="drought"){
  IRES=Results_obs$catrest
  IR_locs=which(IRES$IRES==1)
  RLGPD_obs=-Results_obs$RetLevGPD
  for (id in c(1:70)){
    RLGPD_obs[which(RLGPD_obs[,id]<0),id]=0
    RLGPD_obs[which(is.infinite(RLGPD_obs[,id])),id]=NA
  }
  
  RLIRobs<- RLGPD_obs[IR_locs,]
  RLGPD_obs[IR_locs,c(1:70)]<-NA
  
  Peaks_obs=-Results_obs$Peaks$value
  Peaks_sim=-Results_sim$Peaks$value
  
  IRES=Results_sim$catrest
  IR_locs=which(IRES$IRES==1)
  RLGPD_sim=-Results_sim$RetLevGPD
  for (id in c(1:70)){
    RLGPD_sim[which(RLGPD_sim[,id]<0),id]=0
    RLGPD_sim[which(is.infinite(RLGPD_sim[,id])),id]=NA
  }
  
  RLIRsim<- RLGPD_sim[IR_locs,]
  RLGPD_sim[IR_locs,c(1:70)]<-NA
}else{
  RLGPD_obs=Results_obs$RetLevGPD
  RLGPD_sim=Results_sim$RetLevGPD
}
cn=colnames(RLGPD_obs)
RLGPD_trend=RLGPD_obs-RLGPD_obs[,15]
#year 2005 
yi=which(cn=="2005")
RLGPD_trend=data.frame(ValidSYl_out,RLGPD_trend[,yi])


cn=colnames(RLGPD_sim)
RLGPD_trendsim=RLGPD_sim-RLGPD_sim[,15]
#year 2005 
yi=which(cn=="2005")
RLGPD_trendsim=data.frame(ValidSYl_out,RLGPD_trendsim[,yi])


#Now plot this
plot(RLGPD_trend$RLGPD_trend...yi.,RLGPD_trendsim$RLGPD_trendsim...yi.,ylim=c(-100,100),xlim=c(-100,100))

Recaplus <- st_as_sf(RLGPD_trend, coords = c("Var1.x", "Var2.x"), crs = 4326)
Recaplus <- st_transform(Recaplus, crs = 3035)

names(Recaplus)[22]="trend"
Recaplus$trend_qsp=Recaplus$trend/Recaplus$upa*1000

Recaput <- st_as_sf(RLGPD_trendsim, coords = c("Var1.x", "Var2.x"), crs = 4326)
Recaput <- st_transform(Recaput, crs = 3035)

names(Recaput)[22]="trend"
Recaplus$trendsim=Recaput$trend

Recaplus$trendsim_qsp=Recaput$trend/Recaput$upa*1000
#Do a division by mean discharge!!!
msim=aggregate(list(q=Peaks_sim),
               by = list(Station_ID=Results_sim$Peaks$catch),
               FUN = function(x) c(mean= mean(x,na.rm=T),median=median(x,na.rm=T),max=max(x,na.rm=T)))
msim <- do.call(data.frame, msim)

msim2=rowMeans(RLGPD_sim[,c(1:70)],na.rm=T)
msim2=data.frame(catlist,msim2)
mx=match(msim$Station_ID,msim2$catlist)
plot(msim$q.mean,msim2$msim2[mx])
mobs=aggregate(list(q=Peaks_obs),
               by = list(Station_ID=Results_obs$Peaks$catch),
               FUN = function(x) c(mean= mean(x,na.rm=T),median=median(x,na.rm=T),max=max(x,na.rm=T)))
mobs <- do.call(data.frame, mobs)

mobs=rowMeans(RLGPD_obs[,c(1:70)],na.rm=T)

cor.test(mobs/Recaplus$upa,msim/Recaplus$upa)

trend_oqv<-c()
trend_sqv<-c()
for (st in 1:length(Recaplus$outl2)){
  sta=Recaplus$outl2[st]
  mat=match(sta,catlist)
  #qmean=msim$q.mean[which(msim$Station_ID==sta)]
  qmean=msim2$msim2[mat]
  trend_sqmean=Recaplus$trendsim[st]/qmean*100
  trend_sqv=c(trend_sqv,trend_sqmean)
  
  #qmean=mobs$q.mean[which(mobs$Station_ID==sta)]
  qmean=mobs[mat]
  trend_oqmean=Recaplus$trend[st]/qmean*100
  trend_oqv=c(trend_oqv,trend_oqmean)
  
}

plot(trend_oqv,trend_sqv,ylim=c(-100,100),xlim=c(-100,100))

cor.test(Recaplus$trend_qsp,Recaplus$trendsim_qsp)

plot(Recaplus$trend,Recaplus$trendsim)
cor.test(trend_oqv,trend_sqv)
plot(trend_oqv,trend_sqv)

Recaplus$trend_oqv=trend_oqv
Recaplus$trend_sqv=trend_sqv

#################################################################

#plot strange results

craps=Recaplus[which(abs(Recaplus$trend_sqv)>3e2),]

Station=craps$V1[2]
id=craps$outl2[2]
s=which(Station_data_IDs==Station)
print(ds)
id_obs=match(Station,station_ido)
id_sim=match(Station,Station_data_IDs)
station_ido[id_obs]
Station_data_IDs[id_sim]
uparea=ValidSYl$upa[which(ValidSYl$V1==Station)]
woo=ValidSYl[which(ValidSYl$V1==Station),]
dset="sim"
#Select the first station
if (dset=="obs"){
  Q_loc=as.numeric(Q_obs_val[,id_obs+1])
  Q=Q_loc[-1]
}else if (dset=="sim"){
  Q_loc=as.numeric(Q_sim_val[,id_sim+1])
  Q=Q_loc[-1]
}

txx=time2
catch=as.numeric(ValidSYl_out$outl2[s])
timeStamps=txx
thresh=thresh_vec[which(thresh_vec$cid==catch),]
thresh=thresh$th
frosttime=NA
series=data.frame(txx,Q)

plot(series)

#get the result curve for that location
mid=match(Station,ValidSYl_out$V1)
merd=as.numeric(RLGPD_trendsim[mid,])
dass=as.numeric(Results_sim$RetLevGPD[mid,])

pdas=(Results_sim$parameters[which(Results_sim$parameters[,1]==id),])



#blss=Results_sim$parameters[match(craps$outl2,Results_sim$parameters[,1]),]

#comes from the shape parameter (lool)

blss=which(Results_sim$parameters[,10]>(1))
shp_rm=unique(Results_sim$parameters[blss,1])

blss2=which(Results_obs$parameters[,10]>(1))
shp_rm2=unique(Results_obs$parameters[blss2,1])

Recaplus$trend_oqv[match(shp_rm2,Recaplus$outl2)]=NA
Recaplus$trend_sqv[match(shp_rm,Recaplus$outl2)]=NA
plot(Recaplus$trend_oqv,Recaplus$trend_sqv)
cor.test(Recaplus$trend_oqv,Recaplus$trend_sqv,na.rm=T)
map=ggplot(basemap) +
  geom_sf(fill="gray95", color=NA) +
  geom_sf(data=Recaplus,aes(geometry=geometry,fill=trend_oqv,size=upa),color="transparent",alpha=.9,shape=21,stroke=0)+ 
  geom_sf(data=Recaplus,aes(geometry=geometry,size=upa),col="grey20",alpha=1,stroke=0.05,shape=1)+ 
  geom_sf(fill=NA, color="grey20") +
  scale_x_continuous(breaks=seq(-30,40, by=5)) +
  scale_size(range = c(1, 3), trans="sqrt")+
  scale_fill_gradientn(
    colors=palet,n.breaks=5,oob = scales::squish, limits=c(-100,100), name="Trend (% of mean flood peak)") +
  coord_sf(xlim = c(min(nco[,1]),max(nco[,1])), ylim = c(min(nco[,2]),max(nco[,2])))+
  labs(x="Longitude", y = "Latitude") +
  guides(fill = guide_colourbar(barwidth = 0.5, barheight = 12,reverse=F),size="none")+
  theme(axis.title=element_text(size=tsize),
        panel.background = element_rect(fill = "aliceblue", colour = "grey1"),
        panel.border = element_rect(linetype = "solid", fill = NA, colour="black"),
        legend.title = element_text(size=tsize),
        legend.text = element_text(size=osize),
        legend.position = "right",
        panel.grid.major = element_line(colour = "grey85",linetype="dashed"),
        panel.grid.minor = element_line(colour = "grey90"),
        legend.key = element_rect(fill = "transparent", colour = "transparent"),
        legend.key.size = unit(.8, "cm"))


ggsave("D:/tilloal/Documents/LFRuns_utils/TAplots/trend_obs_1965_2005_relative_dr.jpg", map,width=18, height=20, units=c("cm"),dpi=1000)

map

Recapsame=Recaplus[which(Recaplus$Rlen>21900),]

map=ggplot(basemap) +
  geom_sf(fill="gray95", color=NA) +
  geom_sf(data=Recapsame,aes(geometry=geometry,fill=trend_oqv,size=upa),color="transparent",alpha=.9,shape=21,stroke=0)+ 
  geom_sf(data=Recapsame,aes(geometry=geometry,size=upa),col="grey20",alpha=1,stroke=0.05,shape=1)+ 
  geom_sf(fill=NA, color="grey20") +
  scale_x_continuous(breaks=seq(-30,40, by=5)) +
  scale_size(range = c(1, 3), trans="sqrt")+
  scale_fill_gradientn(
    colors=palet,n.breaks=10,oob = scales::squish, limits=c(-50,50), name="Trend (% of mean peak)") +
  coord_sf(xlim = c(min(nco[,1]),max(nco[,1])), ylim = c(min(nco[,2]),max(nco[,2])))+
  labs(x="Longitude", y = "Latitude") +
  guides(fill = guide_colourbar(      title.position = "top",  # Position the title above the bar
                                      title.hjust = 0.5,       # Center-align the title
                                      barheight = 0.5, barwidth = 30,reverse=F),size="none")+
  theme(axis.title=element_text(size=tsize),
        panel.background = element_rect(fill = "aliceblue", colour = "grey1"),
        panel.border = element_rect(linetype = "solid", fill = NA, colour="black"),
        legend.title = element_text(size=tsize),
        legend.text = element_text(size=osize),
        legend.position = "bottom",
        panel.grid.major = element_line(colour = "grey85",linetype="dashed"),
        panel.grid.minor = element_line(colour = "grey90"),
        legend.key = element_rect(fill = "transparent", colour = "transparent"),
        legend.key.size = unit(.8, "cm"))


ggsave("D:/tilloal/Documents/LFRuns_utils/TAplots/trend_obs_1965_2005_pcmpeak_dr.jpg", map,width=18, height=20, units=c("cm"),dpi=1000)


# Part 3: Trend comparisons------------------------------------

max(Recaplus$Rlen)
Recaplus$trendiff=Recaplus$trend_oqv-Recaplus$trend_sqv
min_rec=60*365.25
Recapsame=Recaplus[which(Recaplus$Rlen>min_rec),]
#Recapsame=Recaplus
#Recapsame$trend_modqsp=Recapsame$trendsim_qsp
if(length(which(is.na(Recapsame$trend_oqv)))>0) Recapsame=Recapsame[-which(is.na(Recapsame$trend_oqv)),]
if(length(which(is.na(Recapsame$trend_sqv)))>0) Recapsame=Recapsame[-which(is.na(Recapsame$trend_sqv)),]

Recapsame$tdcat="A"
Recapsame$tdcat[which((Recapsame$trend_oqv*Recapsame$trend_sqv)>0)]="B"
Recapsame$abtd=abs(Recapsame$trendiff)
Recapsame$abtd[which(Recapsame$tdcat=="A")]=Recapsame$abtd[which(Recapsame$tdcat=="A")]

Rsa=Recapsame[which(Recapsame$tdcat=="A"),]
Rsb=Recapsame[which(Recapsame$tdcat=="B"),]

length(Rsb$upa)/length(Recapsame$upa)
paletb=c(hcl.colors(9, palette = "Blues", alpha = NULL, rev = T, fixup = TRUE))
paleta=c(hcl.colors(9, palette = "Reds", alpha = NULL, rev = T, fixup = TRUE))

library(ggnewscale)
map=ggplot(basemap) +
  geom_sf(fill="gray95", color=NA) +
  geom_sf(data=Rsa,aes(geometry=geometry,fill=abtd,size=upa),color="transparent",alpha=.9,shape=21,stroke=0)+ 
  scale_fill_gradientn(
    colors=paleta,n.breaks=7,oob = scales::squish,  labels=NULL, limits=c(0,100), name="                      | different sign      ",guide = "coloursteps") +
  
  #guides(fill = guide_colourbar(barwidth = 15, barheight = 1,reverse=F),size="none")+
  new_scale_fill()+
  geom_sf(data=Rsb,aes(geometry=geometry,fill=abtd,size=upa),color="transparent",alpha=.9,shape=21,stroke=0)+ 
  scale_fill_gradientn(
      colors=paletb,n.breaks=7,oob = scales::squish,limits=c(0,100), name="Absolute difference | same sign  ",guide = "coloursteps") +
  #guides(fill = guide_coloursteps(barwidth = 15, barheight = unit(.5, "lines"), reverse = FALSE)) +
  geom_sf(data=Recapsame,aes(geometry=geometry,size=upa),col="grey20",alpha=1,stroke=0.05,shape=1)+ 
  geom_sf(fill=NA, color="grey20") +
  scale_x_continuous(breaks=seq(-30,40, by=5)) +
  scale_size(range = c(1, 3), trans="sqrt")+
  coord_sf(xlim = c(min(nco[,1]),max(nco[,1])), ylim = c(min(nco[,2]),max(nco[,2])))+
  labs(x="Longitude", y = "Latitude") +
  guides(size="none")+
  theme(axis.title=element_text(size=tsize),
        panel.background = element_rect(fill = "aliceblue", colour = "grey1"),
        panel.border = element_rect(linetype = "solid", fill = NA, colour="black"),
        legend.text = element_text(size=osize),
        legend.key.width  = unit(3, "lines"),
        legend.key.height = unit(.5, "lines"),
        legend.position = "bottom",
        #legend.title = element_text(size = osize),
        legend.spacing.y = unit(0.1, "cm"),
        legend.title = element_text(size = osize, margin = margin(t = 2, r = 0, b = 2, l = 4)),
        legend.box = "vertical",  # Stacks legends vertically
        panel.grid.major = element_line(colour = "grey85",linetype="dashed"),
        panel.grid.minor = element_line(colour = "grey90"),
        legend.key = element_rect(fill = "transparent", colour = "transparent"),
        legend.key.size = unit(.8, "cm"))+
  ggtitle("difference in 10-y RL flood change between observed and simulated flow")

map


ggsave("D:/tilloal/Documents/LFRuns_utils/TAplots/trend_difference_1965_2005_relativ_dr.jpg", map,width=18, height=20, units=c("cm"),dpi=1000)


#agreement in sign of change

lp1=length(which(Recapsame$trend_oqv>0 & Recapsame$trend_sqv >0))
ln1=length(which(Recapsame$trend_oqv<=0 & Recapsame$trend_sqv<=0))

lp2=length(which(Recapsame$trend_oqv>0 & Recapsame$trend_sqv<=0))
ln2=length(which(Recapsame$trend_oqv<=0 &Recapsame$trend_sqv>0))

lt=lp1+lp2+ln1+ln2
truesign=round((lp1+ln1)/lt*100,1)


# Recaplus=Recaplus[-which(is.na(Recaplus$sen_simlf)),]
# Recaplus=Recaplus[-which(is.na(Recaplus$sen_obs)),]
lmod=(lm(trend_oqv ~ trend_sqv, data=Recapsame))
R2lf=round(summary(lm(trend_oqv ~ trend_sqv, data=Recapsame))$r.squared,2)
sqrt(R2lf)
rmself=round(sqrt(mean((Recapsame$trend_oqv - Recapsame$trend_sqv)^2,na.rm=T)),1)

wallou=expression(paste("R",(km^2)))
palet=c(hcl.colors(11, palette = "RdYlBu", alpha = NULL, rev = F, fixup = TRUE))
tsize=16
osize=16
lims=c(-40,40)
if(length(which(is.na(Recapsame$trend_oqv)))>0) Recapsame=Recapsame[-which(is.na(Recapsame$trend_oqv)),]
if(length(which(is.na(Recapsame$trend_sqv)))>0) Recapsame=Recapsame[-which(is.na(Recapsame$trend_sqv)),]

vlim=100

Recapsame$density <- get_density(Recapsame$trend_oqv, Recapsame$trend_sqv, n = 1000)
ps<-ggplot(data=Recapsame, aes(x=trend_oqv, y=trend_sqv, size=upa, color=density)) +
  annotate("rect", xmin=-Inf, xmax=0, ymin=-Inf, ymax=0, fill="lightskyblue", alpha=0.4) +

  annotate("rect", xmin=0, xmax=Inf, ymin=-Inf, ymax=0, fill="lightsalmon", alpha=0.4) +

  annotate("rect", xmin=0, xmax=Inf, ymin=0, ymax=Inf, fill="lightskyblue", alpha=0.4) +

  annotate("rect", xmin=-Inf, xmax=0, ymin=0, ymax=Inf, fill="lightsalmon", alpha=0.4) +
  #geom_abline(slope = coef(lmod)[["trend_sqv"]], 
  #            intercept = coef(lmod)[["(Intercept)"]])+

  geom_point(alpha=0.6) +
  geom_smooth(method='lm', color="black", fill="grey", fullrange=TRUE)+
  #geom_point(data= Rsig, aes(x=ObsChange, y=SimChange),fill="transparent", color="gray25",shape=21,size=2)+
  geom_abline(slope=1,intercept=0,col="gray25",lwd=1.5,alpha=.5)+
  scale_color_viridis(option="F")+
  # 
  # annotate("label", x=-40, y=50, label= paste0("N = ",ln2),size=5)+
  # annotate("label", x=40, y=50, label= paste0("N = ",lp1),size=5)+
  # annotate("label", x=40, y=-50, label= paste0("N = ",lp2),size=5)+
  # annotate("label", x=-40, y=-50, label= paste0("N = ",ln1),size=5)+
  
  annotate("label", x=-60, y=80, label= paste0("r = ",round(sqrt(R2lf),2),"\nAgreement = ",truesign,"%"),size=5)+
  scale_size(range = c(2, 8), trans="sqrt",name= expression(paste("Upstream area ", (km^2),
                                                                  sep = " ")),
             breaks=c(101,1000,10000,100000,500000), labels=c("100","1000", "10 000", "100 000", "500 000"))+
  scale_x_continuous(name="Observed change (% of mean peak discharge)",breaks = seq(-vlim,vlim,by=25),limits=c(-vlim,vlim),expand = c(0, 0))+
  scale_y_continuous(name="Simulated change (% of mean peak discharge)",breaks = seq(-vlim,vlim,by=25),limits=c(-vlim,vlim),expand = c(0, 0))+
  # scale_color_gradientn(
  #   colors=palet, limits=c(0,1),oob = scales::squish,
  #   name="temporal correlation",breaks=seq(0,1, by=0.2))+
  # scale_alpha_continuous(name="trend significance",range = c(0.5, 1))+
  guides( size= "none",col="none")+
  theme(axis.title=element_text(size=tsize),
        panel.background = element_rect(fill = "white", colour = "grey1"),
        panel.border = element_rect(linetype = "solid", fill = NA, colour="black"),
        legend.title = element_text(size=tsize),
        legend.text = element_text(size=osize),
        axis.text =element_text(size=tsize),
        legend.position = "right",
        panel.grid.major = element_line(colour = "grey70"),
        panel.grid.minor = element_line(colour = "grey90"),
        legend.key = element_rect(fill = "transparent", colour = "transparent"),
        legend.key.size = unit(.8, "cm"))+
  ggtitle(paste0("HERA - Observed and simulated trends of\n10-y RL ",haz," (1965-2005) (n = ",lt," stations)"))# adjust transparency range if needed

ps

ggsave("D:/tilloal/Documents/LFRuns_utils/TAplots/scatter2_correlation_HERAvsobs_19652005_Ltseva2.jpg", ps, width=20, height=20, units=c("cm"),dpi=300)

