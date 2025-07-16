##########################################################################################
############   THIS SCRPT IS FOR RUNNING 1 SQUARE OF THE HERA DOMAIN ON A HPC ############
##########################################################################################

#Library importation ------------------
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
suppressWarnings(suppressMessages(library(evd)))
suppressWarnings(suppressMessages(library(POT)))
suppressWarnings(suppressMessages(library(RtsEva)))

#  FUNCTIONS   #################################################


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

# Arguments importation #############################
args <- commandArgs(TRUE)
argus=as.vector(unlist(strsplit(args, split = " ")))
Nsq = argus[1]
haz = argus [2]
outlets = argus [3]
startid = argus [4]
endid = argus [5]
Quarter = argus [6]
sce = argus [7]
foldin= argus [8]
season= argus [9]

#hard coded arguments
# Nsq = 43
# haz = "flood"
# outlets = "RNetwork"
# startid = 1
# endid = 5
# Quarter = 1
# sce = "calibrated"
# foldin= "Calout"
# season= "nonfrost"

var = "dis"

#default: tail is high
tail="high"

workDir = "/BGFS/CLIMEX/tilloal/HydroMeteo/"
setwd(workDir)

rspace= read.csv(paste0(workDir,"subspace_efas.csv"))
rspace=rspace[,-1]
nrspace=rspace[Nsq,]
print(nrspace)

outletname="efas_rnet_100km_01min"
nameout="UCRnet"
outhybas=outletopen(workDir,outletname,nrspace)
hydroDir<-paste0("/BGFS/CLIMEX/tilloal/HydroMeteo/Timeseries/dis6_",foldin)
Idstart=as.numeric(Nsq)*100000
if (length(outhybas$outlets)>0){
  outhybas$outlets=seq((Idstart+1),(Idstart+length(outhybas$outlets)))
}

unikout=outhybas$outlets
outhybas$latlong=paste(round(outhybas$Var1,4),round(outhybas$Var2,4),sep=" ")
file_out=file=paste0(workDir,"TSEVA/out/",sce,"/",haz,"/",season,"/ResCat6h_",outlets,haz,"_",Nsq,Quarter,"_1951_2020.Rdata")
fex=file.exists(file_out)
if (fex==TRUE){
  f.size=as.numeric(fs::file_size(file_out))
}else{
  f.size=0
}

#Scenario differentiation
if (sce=="Histo") code="h"
if (sce=="SCF") code="scf"
if (sce=="calibrated") code=""
if (sce=="calibratedCCF") code="c"
if (sce=="calibratedRCF") code="r"
if (sce=="calibratedRWCF") code="r"
if (sce=="calibratedSCF") code=""
if (sce=="WStat") code="wcf"
if (sce=="RWStat") code="rwcf"
dirout=paste0(workDir,"TSEVA/out/",sce,"/",haz,"/",season)
fex2=file.exists(dirout)
if (fex2==TRUE){
  print("output folder exists")
}else{
  print("output folder does not exist")
}

#Loop on all pixels within squarq
if (length(outhybas$outlets)>0 & f.size<1500){
  #Load the file
  #loading the files as netcdf (needs to be checked offline)
  if (code==""){
	filename=paste0("dis_",Nsq,"_1950_2020_",code,"cf")
  }
  if (code=="c"){
	filename=paste0("dis_",Nsq,"_1951_2020_",code,"cf")
  }
  if (code=="r"){
	filename=paste0("dis_",Nsq,"_1951_2020_",code,"")
  }
  if (code=="h"){
	filename=paste0("dis_",Nsq,"_1951_2020_",code,"")
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
  dists=disNcopenloc(filename,hydroDir,outhybas,1)
  df.dis=dists 
  print(paste0("hazard: ",haz," opening square ", Nsq, " /88"))
  timeStamps=(as.Date(df.dis$time,origin="1979-01-01"))
  timeStamps=as.POSIXct(timeStamps-1/24)
  txx=timeStamps
  df.dis$timeStamps=txx
  
  names(df.dis)[c(1,2)]=c("dis","outlets")
  
  #loading the frost file for drought
  if (haz=="drought"){
    load(file=paste0(workDir,"TSEVA/in/catchment_frost.Rdata"))
	rmv=which(year(frostcat$time)==1950)
	frostcat=frostcat[-rmv,]
	#remove first day
	frostcat=frostcat[-1,]
    Catchmentrivers7=read.csv(paste0(workDir,"TSEVA/in/from_hybas_eu_onlyid.csv"),encoding = "UTF-8", header = T, stringsAsFactors = F)
    outletname="outletsv8_hybas07_01min"
    outhyb07=outletopen(workDir,outletname,nrspace)
    catmatch=match(outhyb07$outlets,Catchmentrivers7$pointid)
    mycat=Catchmentrivers7[catmatch,]
    
    hybas07 <- read_sf(dsn = paste0(workDir,"TSEVA/in/hybas_eu_lev07_v1c.shp"))
    hybasf7=fortify(hybas07) 
    Catamere07=inner_join(hybasf7,Catchmentrivers7,by= "HYBAS_ID")
    Catamere07$llcoord=paste(round(Catamere07$POINT_X,4),round(Catamere07$POINT_Y,4),sep=" ") 
    Catf7=inner_join(Catamere07,outhybas,by= c("llcoord"="latlong"))
    st_geometry(Catf7)=NULL	
	tail="low"
    
  }
 
ThDir<-paste0(workDir,"TSEVA/out/Thresholds")
TH1=read.csv(paste0(ThDir,"/trenTH_Histo_",tail,"_",Nsq,".csv"))
TH2=read.csv(paste0(ThDir,"/trenTH_SCF_",tail,"_",Nsq,".csv"))
TH3=inner_join(TH1,TH2,by="cid")

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
      Nonstat<-TsEvaNs(timeAndSeries, timeWindow, transfType='trendPeaks',ciPercentile = ciPercentile, minPeakDistanceInDays = minPeakDistanceInDays,lowdt=7,trans=trans,tail = tail,TrendTh = thresh)
      nonStationaryEvaParams=Nonstat[[1]]
      stationaryTransformData=Nonstat[[2]]
	
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
  
  
  Results=list(parameters=parlist,RetLevGEV=RetLevGEV,RetLevGPD=RetLevGPD,RetPerGEV=RetPerGEV,RetPerGPD=RetPerGPD,Peaks=peaklist,catrest=data.frame(catlist,IRES))
  save(Results, file=paste0(workDir,"TSEVA/out/",sce,"/",haz,"/",season,"/ResCat6h_",outlets,haz,"_",Nsq,Quarter,"_1951_2020.Rdata"))
  
  
}else if (length(outhybas$outlets)==0){
  print(paste0("There is no river pixel in the box ",Nsq))
}else if (f.size>=1500){
  print(paste0("Output file already exists for box ",Nsq))
}


