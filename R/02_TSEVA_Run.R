##########################################################################################
############   THIS SCRPT IS FOR RUNNING 1 SQUARE OF THE HERA DOMAIN ON A HPC ############
##########################################################################################

setwd(dirname(rstudioapi::getSourceEditorContext()$path))
getwd()
source("functions_trends.R")

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

startid=1
endid=length(unikout)
endid=3

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




