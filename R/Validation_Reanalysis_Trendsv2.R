library(raster)
library(rgdal)
library(sf)
library(tidyverse)
library(ncdf4)
library(sp)
library(reshape)
library(MASS)
library(ggplot2)
library(viridis)
library(rnaturalearth)
library(hydroGOF)
library(readxl)
source("~/LFRuns_utils/TrendAnalysis/R/function_loading.R")

get_density <- function(x, y, ...) {
  dens <- MASS::kde2d(x, y, ...)
  ix <- findInterval(x, dens$x)
  iy <- findInterval(y, dens$y)
  ii <- cbind(ix, iy)
  return(dens$z[ii])
}

# Files paths -------------------------------------------------------------
main_path = 'D:/tilloal/Documents/06_Floodrivers/'
valid_path = paste0(main_path,'DataPaper/')
dis_path<-paste0(main_path,'dis/calibrated/filtered/Histo/')
setwd(valid_path)
outletname = "outletsv8_hybas07_01min"
hydroDir<-("D:/tilloal/Documents/LFRuns_utils/data")

# Part 1: Create the Valid station file -------------------------------------------

# Part 2: TSEVA on observed values------------------------------------
# TSEVA on observed values
ValidSf=read.csv(file="Stations/Stations_ValidationF.csv")[,-1]
ValidSY=ValidSf[which(ValidSf$removal!="YES"),]
ValidSY$Rlenyr=ValidSY$Rlen/365

#Here I need at least 40 years betwee 1960 and 2010
#ValidSY=ValidSY[-which(ValidSY$Rlenyr<40),]
#load frost days
# load(file=paste0(hydroDir,"/Drought/catchment_frost.Rdata"))
# saveRDS(frostcat,paste0(hydroDir,"/Drought/catchment_frost.RDS"))
#frostcat=readRDS(paste0(hydroDir,"/Drought/catchment_frost.RDS"))
#Hybas07
#Catchmentrivers7=read.csv(paste0(hydroDir,"/Catchments/from_hybas_eu_onlyid.csv"),encoding = "UTF-8", header = T, stringsAsFactors = F)




## MHM loading ---------------
mHM_loc=read.table(paste0(valid_path,"Revisions/mHM_EU/mHM_locations_EFAS_grid.txt"),sep=",")

#transalting to R indexing after python
mHM_loc$V2=as.numeric(mHM_loc$V2)+1
mHM_loc$V3=as.numeric(mHM_loc$V3)+1

#load the xls file for catchment area
# Load the Excel file into a tibble
data <- read_excel(paste0(valid_path,"Revisions/mHM_EU/european_catchments_with_filename.xlsx"))


mHM_sta=full_join(data,mHM_loc, by=c("filename"="V1"))


#load the mHM discharge
mHM_dis=read.table(paste0(valid_path,"Revisions/mHM_EU/Q_mHM_filename.txt"),sep=" ")
mHM_dis=as.data.frame(mHM_dis)
colnames(mHM_dis)=mHM_dis[1,]
mHM_dis=mHM_dis[-1,]
mHM_dis[] <- lapply(mHM_dis, as.numeric)
name_vector=colnames(mHM_dis)
#extract mhm_station that are in the dataset

matm=na.omit(match(name_vector,mHM_sta$filename))
mHM_sta=mHM_sta[matm,]
#find the date
date1=seq(as.Date("1960-01-01"),as.Date("2010-12-31"),by="days")



## Observation loading --------------------
Q_data <- read.csv(paste0('Q_19502020.csv'), header = T)  # CSVs with observations
Station_data_IDs <- as.vector(t(Q_data[1, -c(1,2)]))
#remove first days to match with HERA

Q_data=Q_data[-c(1:4),]

## HERA loading ---------------------

HERA_data<-read.csv(file="out/HERA_Val2_19502020.csv")
HERA_cordata<-read.csv(file="out/HERA_CorStat_19502020.csv")

replax=match(as.numeric(HERA_cordata[1,]),as.numeric(HERA_data[1,]))[-1]
HERA_data[,replax]=HERA_cordata[,-1]

date2=seq(as.Date("1950-01-04"),as.Date("2020-12-31"),by="days")


#Sample Q_data and HERA only for the 1960-2010 period
# mtime=match(date1,date2)
# Q_data=Q_data[mtime,-1]
######

# lR=c()
# for (id in 1:length(Station_data_IDs)){
#   Qs=Q_data[,id+1]
#   xl=length(Qs)
#   l=length(which(!is.na(Qs)))
#   lR=c(lR,l)
# }
# 
# #match station data with validSY
# val_lr=match(ValidSY$V1,Station_data_IDs)
# 
# ValidSY$Rlen_trend=lR[val_lr]/365.25

#keep only stations with more than 30yrs
ValidST=ValidSY[-which(ValidSY$Rlenyr<=30),]



#now extract the HERA discharge for catchments in ValidSLong
station_HERA<-HERA_data[1,-1]
SelectHera=which(station_HERA %in% ValidST$V1)
station_ids=data.frame(id=t(station_HERA[SelectHera]))
HERA_dates<-date2
HERA_data2<-HERA_data[-1,-1]
#HERA_data2=HERA_data2[mtime,]
HERA_comp<-HERA_data2[,SelectHera]

station_Q<-Station_data_IDs
SelectQ=which(station_Q %in% ValidST$V1)
station_ido=data.frame(id=station_Q[SelectQ])

Q_dates<-date1
Q_data2<-Q_data[,-c(1,2)]
Q_comp<-Q_data2[,SelectQ]
Q_comp=Q_comp[,order(station_ido$id)]
HERA_comp<-HERA_comp[,order(station_ids$X1)]
station_ido=station_ido$id[order(station_ido$id)]
station_ids=station_ids$X1[order(station_ids$X1)]

#match station_ids with mHM
## Distance between "official gauges" and efas points -----------------------
hera=ValidST[,c(1,2,3)]
heraloc=st_as_sf(hera, coords = c("Var1", "Var2"), crs = 4326)
mhml=mHM_sta[,c(3,2,4)]
mhmloc=st_as_sf(mhml, coords = c("LON", "LAT"), crs = 4326)
dist=c()
mlm=c()
for (r in 1:length(heraloc$upa)){
  cat(paste0(r,"\n"))
  v1=heraloc[r,]
  v2=mhmloc
  oula=st_distance(v1,v2)
  oula=oula/1000
  ziz=which.min(oula)
  guez=oula[ziz]
  dist=c(dist,guez)
  mlm=c(mlm,ziz)
}

hist(dist)

ValidST$dist2mhm=dist
mhm_candidate=mHM_sta[mlm,]

finalcom=data.frame(ValidST,mhm_candidate)
finalcom$UpA=finalcom$upa
plot(finalcom$UpA,finalcom$Area_given_km2)
difupa=(finalcom$UpA)/(finalcom$Area_given_km2)

plot(finalcom$dist2mhm[order(finalcom$dist2mhm)], ylim=c(0,10))
hist(difupa,breaks=20000,xlim=c(0,4))
finalcom_upaclean=finalcom[which(difupa<=1.25 & difupa>=0.85 & dist<5),]

un_st=unique(finalcom_upaclean$filename)

#loop to remove double stations
rmf=c()
for (s in 1:length(un_st)){
  st=un_st[s]
  matsta=which(!is.na(match(finalcom_upaclean$filename,st)))
  if (length(matsta)>1){
    print("multi match")
    print(length(matsta))
    print(st)
    torm=which.min(finalcom_upaclean$dist2mhm[matsta])
    rmt=matsta[-torm]
    rmf=c(rmf,rmt)
  }
}

finalcom_upaclean=finalcom_upaclean[-rmf,]



plot(finalcom_upaclean$UpA,finalcom_upaclean$Area_given_km2)




# Fit TSEVA on observed values ------------------------
dates_Hera <- seq(as.Date("1950-01-04"), as.Date("2020-12-31"), by = "days")
dates_Q <- seq(as.Date("1950-01-01"), as.Date("2020-12-31"), by = "days")
id=56
station_ids[id]
station_ido[id]
plot(Q_comp[,id])
plot(HERA_comp[,id-1])
#load TSEVA functions
source("~/LFRuns_utils/TSEVA_demo/demo_functions.R")
library(xts)


#Extract 1 timeseries out of HERA_comp
#HERA_comp is a matrix with the first row being the dates
#and the other rows the discharge values for each station
haz="flood"
season="nonfrost"
dset="obs"
#start the loop
RetPerGPD=c()
RetPerGEV=c()
RetLevGEV=c()
RetLevGPD=c()
parlist=c()
peaklist=c()
catlist=c()
IRES=c()
Reservoir_i=c()
#loop over all stations
for(id in seq_along(station_ids$X1)){
  # id=23
  start_time <- Sys.time()
  Reservoir_alteration=0
  station=station_idplus$X1[id]
  print(paste0(id,"/",length(station_ids$X1)))

  #the hybas catchment for frost days
  catch=station_idplus$pointid.y[id]
  #Select the first station
  if (dset=="obs"){
    Q=Q_comp[,id]
    txx=dates_Q
  }else if (dset=="sim"){
    Q=HERA_comp[,id]
    txx=dates_Hera
  }
  Q=as.numeric(Q)
  df.dis=Q

  rmv=which(year(txx)==1950)
  df.dis=data.frame(txx,df.dis)
  names(df.dis)[c(1,2)]=c("time","dis")
  #remove 1950 which is not reliable
  df.dis=df.dis[-rmv,]
  timeStamps=txx
  dt = difftime(timeStamps[2],timeStamps[1],units="days")
  dt= as.numeric(dt)

  if (haz=="drought"){
    #seasonality divide: frost vs non frost
    percentile=95
    Tcatchment=which(colnames(frostcat)==catch)
    data=df.dis
    names(data)=c("date","Qs")
    intermit=interid(data, WindowSize=7)
    interflag=intermit$flags[2]
    if (!exists("trans")){trans="rev"}
    print(paste0(trans," transformation used for low flows"))
    series=data.frame(txx[-rmv],intermit$trdis$Q7)
    #remove frost timesteps, this can be modified to do the anlysis only on frost moments
    if (length(Tcatchment)>0){
      frostserie=data.frame(frostcat[,1],frostcat[,Tcatchment])
      names(frostserie)=c("time","Ta")
      frostserie=min_daily_value(frostserie)
      rmv2=which(year(frostserie$date)==1950)
      frostserie=frostserie[-rmv2,]
      frostserie=frostserie[-rmv,]
      frosttime=which(frostserie[,2]<0)
    }else{
      frosttime=NA
    }
  }else if (haz=="flood"){
    data=df.dis
    catmat=Catf7[which(Catf7$pointid==catch),]
    names(data)=c("date","Qs")
    # Extract the maximum daily value
    series <- max_daily_value(data)
    # series=data.frame(txx,df.disX$outlets)
    interflag=0
    frosttime=NA
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
  tindexes=match(Impdates,timeDays)
  names(series)=c("timestamp","dis")
  nv=length(unique(series$dis))
  timeAndSeries=series
  names(timeAndSeries)=c("timestamp","data")

  if(length(na.omit(series$dis))>1 & interflag<1 & nv>15){

    if (haz=="drought" && length(!is.na(frosttime))>1) {
      if (season=="nonfrost") {
        print("nonfrost season")
        timeAndSeries$data[frosttime]=NA
      } else if (season=="frost") {
        print("frost season")
        timeAndSeries$data[-frosttime]=NA
      } else {
        print("season must be frost or nonfrost")
      }
    }

    tsm=1/dt
    series=timeAndSeries[,2]
    timeWindow = 365.25*30; #time windows in days, the correction is done within the functions
    windowSize=366
    minPeakDistanceInDays=30
    timeStamps=timeAndSeries$timestamp
    # Choose transformation
    transftypes=c("ori","rev","inv","lninv")
    trendtypes=c("trend","trendPeaks")
    trendtrans=expand.grid(transftypes,trendtypes)
    #choose transformation
    if (haz=="flood") tt=5
    if (haz=="drought") tt=6

    # Apply TsEvaNs function
    Nonstat <- TsEvaNs(timeAndSeries, timeWindow, transfType=trendtrans[tt,2], 
    ciPercentile= 90, minPeakDistanceInDays = minPeakDistanceInDays, 
    mode=haz,TrendTh=0.85, trans=trendtrans[tt,1])

    nonStationaryEvaParams=Nonstat[[1]]
    stationaryTransformData=Nonstat[[2]]

    # Define range for plotting
    ExRange= c(min(nonStationaryEvaParams$potObj$parameters$peaks),max(nonStationaryEvaParams$potObj$parameters$peaks))
    if (haz=="flood") wr2 <- c(seq(min(ExRange),max(ExRange),length.out=700))
    if (haz=="drought") wr2 <- c(seq(1.1*min(ExRange),0.1*max(ExRange),length.out=700))

    # Apply tsEvaPlotGPDImageScFromAnalysisObj function
    #Plot1= tsEvaPlotGPDImageScFromAnalysisObj(wr2, nonStationaryEvaParams, stationaryTransformData, minYear = '1950',trans="rev")  
    timeIndex=1
    Plot2 = tsEvaPlotReturnLevelsGPDFromAnalysisObj(nonStationaryEvaParams, stationaryTransformData, timeIndex, trans=trendtrans[tt,1],ylabel="Discharge (m3/s)",ope=T)
    #print(Plot2)
    ggsave(paste0(hydroDir,"/TSEVA_hybas/valid/GPDBeam_",haz,"_",dset,"_Cat_",station,".jpg"),Plot2, width=24, height=20, units=c("cm"),dpi=300)
    #Plot3 = tsEvaPlotReturnLevelsGEVFromAnalysisObj(nonStationaryEvaParams, stationaryTransformData, timeIndex, trans="rev",ylabel="Discharge (m3/s)")
    
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
    RLevs100=ComputeReturnLevels(nonStationaryEvaParams, RPgoal, timeIndex,trans=trendtrans[tt,1])

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
        RLevs100i=ComputeReturnLevels(nonStationaryEvaParams, RPgoal, timeIndex,trans=trendtrans[tt,1])
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
    cat(paste0("\n No values in this pixel \n or intermittent river (flag = ",interflag,")"))
    if (is.na(interflag)) interflag=-9999
    if (interflag>0){
      filling=intermit$flags[3]
      datex=yday(intermit$zeroRP$time)
      dtect=c(diff(datex),-1)
      last_days <- intermit$zeroRP$time[which(dtect<0)]
      tindexes=match(last_days,intermit$zeroRP$time)
      oops=intermit$zeroRP[tindexes,]
      # verif=series$timestamp[which(!is.na(match(year(series$timestamp),year(Impdates))))]
      # dfrep=data.frame(year(verif),oops)
      # yagg=aggregate(list(m=dfrep$oops),
      #                by = list(year=dfrep$year.verif.),
      #                FUN = function(x) c(max=max(x)))

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
    #parlist=as.data.frame(mapply(c,parlist,params))
    peaklist=as.data.frame((rbind(peaklist,pikos)))

  }

  #Saving main outputs
  catlist=c(catlist,catch)
  IRES=c(IRES,interflag)
  Reservoir_i=c(Reservoir_i,Reservoir_alteration)
  RetLevGEV=rbind(RetLevGEV,RLgev)
  RetLevGPD=rbind(RetLevGPD,RLgpd)
  RetPerGEV=rbind(RetPerGEV,nRPgev)
  RetPerGPD=rbind(RetPerGPD,nRPgpd)
  end_time <- Sys.time()
  cat(paste0("\nloop duration: ",round(end_time-start_time,2)," seconds\n"))
}


#saving the outputs
parlist=as.data.frame(parlist)
Results=list(parameters=parlist,RetLevGEV=RetLevGEV,RetLevGPD=RetLevGPD,RetPerGEV=RetPerGEV,RetPerGPD=RetPerGPD,Peaks=pikos,catrest=data.frame(catlist,Reservoir_i,IRES))


save(parlist,file=paste0(hydroDir,"/TSEVA_hybas/outputs/parValid_",haz,"_",dset,"_1950_2020.Rdata"))
save(Results, file=paste0(hydroDir,"/TSEVA_hybas/outputs/ResValid_",haz,"_",dset,"_1950_2020.Rdata"))


# Part 2.2: Trend computation using Bloshl (2019) method------------------------------------


#Extract annual maxima for Observed and HERA
dates_Hera <- seq(as.Date("1950-01-04"), as.Date("2020-12-31"), by = "days")
dates_Q <- seq(as.Date("1950-01-01"), as.Date("2020-12-31"), by = "days")
id=56
station_ids[id]
station_ido[id]
plot(Q_comp[,id])
plot(HERA_comp[,id])

Q_sim1=data.frame(time=dates_Hera,Q=HERA_comp[,id])
Amax_HERA=computeAnnualMaxima(Q_sim1)

Q_obs=data.frame(time=dates_Hera,Q=Q_comp[,id])
Amax_obs=computeAnnualMaxima(Q_obs)


station_fx=finalcom_upaclean$V1

station_f=station_ido
#keep only station_f discharge data 4all

vecto=which(!is.na((match(station_ido,finalcom_upaclean$V1))))

station_f=station_ido[vecto]
Q_comp=Q_comp[,vecto]
HERA_comp<-HERA_comp[,vecto]
crap=finalcom_upaclean$filename[order(finalcom_upaclean$V1)]
vecto=which(!is.na((match(name_vector,crap))))

name_vector2=name_vector[vecto]

woow=match(crap,name_vector2)

verif=data.frame(name_vector2[woow],crap)

mHM_comp<-mHM_dis[,vecto][woow]
##################################

timeAndSeries=Q_obs
computeAnnualMean<-function(timeAndSeries) {
  timeStamps <- timeAndSeries[,1]
  srs <- timeAndSeries[,2]
  
  tmvec <- as.Date(timeStamps)
  years <- format(tmvec, "%Y")
  
  findMean <- function(subIndxs) {
    meanX <- mean(srs[subIndxs])

  }
  annualMean <- tapply(1:length(srs), years, findMean)
  random_days <- timeAndSeries %>%
    group_by(year = year(time)) %>%
    sample_n(1)
  annualMeanDate <- random_days$time
  annualMeanIndx <- match(annualMeanDate,timeStamps)

  
  return(list(annualMean = annualMean, annualMeanDate = annualMeanDate, annualMeanIndx = annualMeanIndx))
}
# Q_comp=Q_comp[,order(station_f)]
# HERA_comp=HERA_comp[,order(station_f)]
### 1- Anmax for every station-----------------------------
hit_out=c()
rs_obs=c()
rs_sim=c()
#rs_sim2=c()
for (id in 1:length(station_f)){
  print(id)
  s=station_ids[id]
  
  Q_sim1=data.frame(time=date2,Q=HERA_comp[,id])
  #Q_sim2=data.frame(time=date1,Q=mHM_comp[,id])
  Q_obs=data.frame(time=date2,Q=Q_comp[,id])
  rmv=which(as.integer(format(date2, "%Y"))==1950)
  Q_sim1=Q_sim1[-rmv,]
  Q_obs=Q_obs[-rmv,]
  # Q_sim1$Q=tsEvaNanRunningMean(Q_sim1$Q,30)
  # Q_obs$Q=tsEvaNanRunningMean(Q_obs$Q,30)
  iy=which(is.na(Q_obs$Q))
  hit=NA
  if (length(iy)>0){
    Q_obs=Q_obs[-iy,]
    Q_sim1=Q_sim1[-iy,]
    #Q_sim2=Q_sim2[-iy,]
    
  }
  if (length(Q_obs$Q>0)){
    AMAX_obs <- computeAnnualMaxima(Q_obs)
    AMAX_sim <- computeAnnualMaxima(Q_sim1)
    
    # AMAX_obs <- computeAnnualMinima(Q_obs)
    # AMAX_sim <- computeAnnualMinima(Q_sim1)
    # 
    # AMAX_obs <- computeAnnualMean(Q_obs)
    # AMAX_sim <- computeAnnualMean(Q_sim1)
    # 
    # AMAX_obs <- (Q_obs)
    # AMAX_sim <- computeAnnualMinima(Q_sim1)
    #AMAX_sim2 <- computeAnnualMaxima(Q_sim2)
    if (length(AMAX_sim[[3]])==length(AMAX_obs[[3]])){
      AMAX_diff=AMAX_obs[[3]]-AMAX_sim[[3]]
      hit=length(which(abs(AMAX_diff)<8))/length(AMAX_diff)
    }
    data=data.frame(Q_obs,sim=Q_sim1$Q)
    names(data)=c("date","Q","Qs")
    # data=data.frame(Q_obs,sim=Q_sim1$Q,sim2=Q_sim2$Q)
    # names(data)=c("date","Q","Qs","Qs2")
    
    res_obs=data.frame(Station=rep(s,length(AMAX_obs[[1]])),
                       AnMax=AMAX_obs[[1]],AnMaxDate=AMAX_obs[[2]],
                       Year=as.integer(format(AMAX_obs[[2]], "%Y")))
    res_sim=data.frame(Station=rep(s,length(AMAX_sim[[1]])),
                       AnMax=AMAX_sim[[1]],AnMaxDate=AMAX_sim[[2]], 
                       Year=as.integer(format(AMAX_sim[[2]], "%Y")))
    # res_sim2=data.frame(Station=rep(s,length(AMAX_sim2$annualMax)),
    #                    AnMax=AMAX_sim2$annualMax,AnMaxDate=AMAX_sim2$annualMaxDate, 
    #                    Year=as.integer(format(AMAX_sim2$annualMaxDate, "%Y")))
    
    rs_obs=rbind(rs_obs,res_obs)
    rs_sim=rbind(rs_sim,res_sim)
   # rs_sim2=rbind(rs_sim2,res_sim2)
  }else{print("no observations")}
  hit_out=c(hit_out,hit)
}

rs_obsx=as_tibble(rs_obs[,-3])
rs_obsx$Year=as.character(rs_obsx$Year)
rs_obs2=rs_obsx %>% pivot_wider(names_from=Year,values_from=AnMax)

rs_simx=as_tibble(rs_sim[,-3])
rs_simx$Year=as.character(rs_simx$Year)
rs_simlf=rs_simx %>% pivot_wider(names_from=Year,values_from=AnMax)

# rs_simo=as_tibble(rs_sim2[,-3])
# rs_simo$Year=as.character(rs_simo$Year)
# rs_simhh=rs_simo %>% pivot_wider(names_from=Year,values_from=AnMax)

#reading in the data; path needs to be changed to file-location on computer
# setwd("~/06_Floodrivers/DataPaper/European_floods/europe_floods-master")
# data_europe = read.csv2("europe_data.csv",
#                         header=T,stringsAsFactors = F,dec=".")
# 


##############################################################################
#fitting the trend and evaluating the significance
##############################################################################
# Fig. 1 & EDF 2a |Sen-slope and Mann-Kendall

#two-sided test for trend
mann_kendall_trend <- function(dframe,discharge_name="Qmax",year_name="year",continuity=TRUE) {
  
  #Remove NAs
  if(length(which(is.na(dframe[,discharge_name])))>0) {
    dframe <- dframe[-which(is.na(dframe[,discharge_name])),]
  }
  
  #mann-kendall
  mk <- NA
  for(i in c(1:(nrow(dframe)-1) )) {
    for(j in c((i+1):nrow(dframe))) {
      mk <- c(mk,sign(dframe[j,discharge_name]-dframe[i,discharge_name]) )
    }
  }
  mk <- mk[-1]
  mk <- sum(mk,na.rm=T)
  
  #for variance, need number of ties
  series <- dframe[,discharge_name]
  n <- length(series)
  if(any(duplicated(series))) {
    #indices of elements, that occur multiple times
    duplicate_indices <- which(duplicated(series))
    #number of elements with multiple occurences
    p <- length(unique(series[duplicate_indices]))
    #first occurrence of element with multiple occurrences
    t_index <- sapply(unique(series[duplicate_indices]), function(k) min(which(series==k)) )
    #number of occurrences of a duplicated element
    t <- sapply(unique(series[duplicate_indices]), function(k) length(which(series==k)))
    
    var_mk <- (n*(n-1)*(2*n+5) - sum(sapply(t, function(j) j*(j-1)*(2*j+5))) )/18
  }
  #no ties
  else if(!any(duplicated(series))) {
    var_mk <- (n*(n-1)*(2*n+5))/18
  }
  #actual test-statistic
  if(continuity==FALSE) {
    Z <- mk/sqrt(var_mk)
  }
  else if(continuity==TRUE) {
    Z <- sign(mk)*(abs(mk)-1)/sqrt(var_mk)
  }
  #p-value with normal distribution, assuming two-sided test here
  p_val <- (1-pnorm(abs(Z)))*2
  
  return(c(sign_sum=mk,test_stat=Z,p_val=p_val))
}
#calculating sen-slope on original scale and percentage of mean per decade
sen_slope <- function(dframe,discharge_name="Qmax",year_name="year") {
  
  #check for duplicated values 
  if(length(which(duplicated(dframe)))>0) {
    dframe <- dframe[-which(duplicated(dframe)),]
  }
  #check for NAs
  if(length(which(is.na(dframe[,discharge_name])))>0) {
    dframe <- dframe[-which(is.na(dframe[,discharge_name])),]
  }
  #slopes
  mk <- 0
  for(i in c(1:(nrow(dframe)-1) )) {
    for(j in c((i+1):nrow(dframe))) {
      mk <- c(mk,ifelse(!is.na((dframe[j,discharge_name]-dframe[i,discharge_name])/(dframe[j,year_name]-dframe[i,year_name]) ),
                        (dframe[j,discharge_name]-dframe[i,discharge_name])/(dframe[j,year_name]-dframe[i,year_name]) ,
                        NA))
      
    }
  }
  mk <- mk[-1]
  sen <- median(mk,na.rm=TRUE)
  
  mean_flow <- mean(dframe[,discharge_name],na.rm=TRUE)
  #slope as percentage of mean per decade
  sen_pct_of_mean_per_dec = (sen/mean_flow)*100*10
  return(c(sen=sen,sen_pct_of_mean_per_dec=sen_pct_of_mean_per_dec))
}


#applying to data
rs_obs2$mk_pval = rep(NA,nrow(rs_obs2))
rs_obs2$sen_original = rep(NA,nrow(rs_obs2))
rs_obs2$pct_of_mean_per_decade = rep(NA,nrow(rs_obs2))

rs_simlf$mk_pval = rep(NA,nrow(rs_simlf))
rs_simlf$sen_original = rep(NA,nrow(rs_simlf))
rs_simlf$pct_of_mean_per_decade = rep(NA,nrow(rs_simlf))
rs_simlf$cor = rep(NA,nrow(rs_simlf))

rs_simhh$mk_pval = rep(NA,nrow(rs_simhh))
rs_simhh$sen_original = rep(NA,nrow(rs_simhh))
rs_simhh$pct_of_mean_per_decade = rep(NA,nrow(rs_simhh))
rs_simhh$cor = rep(NA,nrow(rs_simhh))

ylen=71
for(j in (1:length(rs_obs2$Station)) ) {
  print(j)
  #p-value of mann-kendall test
  rs_obs2$mk_pval[j] = mann_kendall_trend(dframe=data.frame(year=c(1951:2020),Qmax=as.numeric(rs_obs2[j,c(2:ylen)]) ) )[3] 
  #sen-slope on scale of data and as percentage of mean per decade
  sen_tmp = sen_slope(dframe=data.frame(year=c(1951:2020),Qmax=as.numeric(rs_obs2[j,c(2:ylen)]) ) )
  rs_obs2$sen_original[j] = sen_tmp[1]
  rs_obs2$pct_of_mean_per_decade[j] = sen_tmp[2]
  
  #p-value of mann-kendall test
  rs_simlf$mk_pval[j] = mann_kendall_trend(dframe=data.frame(year=c(1951:2020),Qmax=as.numeric(rs_simlf[j,c(2:ylen)]) ) )[3] 
  #sen-slope on scale of data and as percentage of mean per decade
  sen_tmp = sen_slope(dframe=data.frame(year=c(1951:2020),Qmax=as.numeric(rs_simlf[j,c(2:ylen)]) ) )
  rs_simlf$sen_original[j] = sen_tmp[1]
  rs_simlf$pct_of_mean_per_decade[j] = sen_tmp[2]
  
  #p-value of mann-kendall test
  # rs_simhh$mk_pval[j] = mann_kendall_trend(dframe=data.frame(year=c(1960:2010),Qmax=as.numeric(rs_simhh[j,c(2:52)]) ) )[3] 
  # #sen-slope on scale of data and as percentage of mean per decade
  # sen_tmp = sen_slope(dframe=data.frame(year=c(1960:2010),Qmax=as.numeric(rs_simhh[j,c(2:52)]) ) )
  # rs_simhh$sen_original[j] = sen_tmp[1]
  # rs_simhh$pct_of_mean_per_decade[j] = sen_tmp[2]
  
  Qmax_simlf=(as.numeric(rs_simlf[j,c(2:ylen)]))
  #Qmax_simhh=(as.numeric(rs_simhh[j,c(2:52)]))
  Qmax_obs=(as.numeric(rs_obs2[j,c(2:ylen)]))
  if (length(which(is.na(Qmax_simlf)))>0){
    Qmax_obs=Qmax_obs[-which(is.na(Qmax_simlf))]
   # Qmax_simhh=Qmax_simhh[-which(is.na(Qmax_simlf))]
    Qmax_simlf=Qmax_simlf[-which(is.na(Qmax_simlf))]
  }
  # Kor=cor(Qmax_simhh,Qmax_obs)
  # print(Kor)
  # rs_simhh$cor[j]=Kor
  
  Kor=cor(Qmax_simlf,Qmax_obs)
  print(Kor)
  rs_simlf$cor[j]=Kor
  
}

plot(rs_simlf$cor,rs_simhh$cor)
abline(a=0,b=1)
mean(rs_simlf$cor)
mean(rs_simhh$cor)
#Now I have to compare the trends

# Recap=data.frame(stat=rs_obs2$Station, pval_obs=rs_obs2$mk_pval,pval_simlf=rs_simlf$mk_pval,pval_simhh=rs_simhh$mk_pval,
#                  sen_obs=rs_obs2$pct_of_mean_per_decade,sen_simlf=rs_simlf$pct_of_mean_per_decade,
#                  sen_simhh=rs_simhh$pct_of_mean_per_decade, corr_lf=rs_simlf$cor, corr_hh=rs_simhh$cor)

Recap=data.frame(stat=rs_obs2$Station, pval_obs=rs_obs2$mk_pval,pval_simlf=rs_simlf$mk_pval,
                 sen_obs=rs_obs2$pct_of_mean_per_decade,sen_simlf=rs_simlf$pct_of_mean_per_decade,
                  corr_lf=rs_simlf$cor)
#Recap=Recap[-which(is.na(Recap$sen_simlf)),]
Rsig=Recap[which(Recap$pval_obs<0.1),]


#I need to match the coordinates with Recap

Recaplus=inner_join(Recap, ValidST,by = c("stat"="V1"))

# Part 3: Trend comparisons------------------------------------


#agreement in sign of change
lp1=length(which(Recap$sen_obs>0 & Recap$sen_simlf>0))
ln1=length(which(Recap$sen_obs<=0 & Recap$sen_simlf<=0))

lp2=length(which(Recap$sen_obs>0 & Recap$sen_simlf<=0))
ln2=length(which(Recap$sen_obs<=0 & Recap$sen_simlf>0))

lt=lp1+lp2+ln1+ln2
truesign=round((lp1+ln1)/lt*100)

lp1=length(which(Recap$sen_obs>0 & Recap$sen_simhh>0))
ln1=length(which(Recap$sen_obs<=0 & Recap$sen_simhh<=0))

lp2=length(which(Recap$sen_obs>0 & Recap$sen_simhh<=0))
ln2=length(which(Recap$sen_obs<=0 & Recap$sen_simhh>0))

lt=lp1+lp2+ln1+ln2
truesign=round((lp1+ln1)/lt*100)

# Recaplus=Recaplus[-which(is.na(Recaplus$sen_simlf)),]
# Recaplus=Recaplus[-which(is.na(Recaplus$sen_obs)),]

R2lf=round(summary(lm(sen_obs ~ sen_simlf, data=Recap))$r.squared,2)
R2hh=round(summary(lm(sen_obs ~ sen_simhh, data=Recap))$r.squared,2)

rmself=round(sqrt(mean((Recaplus$sen_obs - Recaplus$sen_simlf)^2)),1)
rmsehh=round(sqrt(mean((Recaplus$sen_obs - Recaplus$sen_simhh)^2)),1)

wallou=expression(paste("R",(km^2)))
palet=c(hcl.colors(11, palette = "RdYlBu", alpha = NULL, rev = F, fixup = TRUE))
tsize=16
osize=16
lims=c(-40,40)

Recaplus$density <- get_density(Recaplus$sen_obs, Recaplus$sen_sim, n = 100)
ps<-ggplot() +
  # annotate("rect", xmin=-Inf, xmax=0, ymin=-Inf, ymax=0, fill="lightskyblue", alpha=0.4) +
  # 
  # annotate("rect", xmin=0, xmax=Inf, ymin=-Inf, ymax=0, fill="lightsalmon", alpha=0.4) +
  # 
  # annotate("rect", xmin=0, xmax=Inf, ymin=0, ymax=Inf, fill="lightskyblue", alpha=0.4) +
  # 
  # annotate("rect", xmin=-Inf, xmax=0, ymin=0, ymax=Inf, fill="lightsalmon", alpha=0.4) +
  # 
  geom_point(data=Recaplus, aes(x=sen_obs, y=sen_simlf, size=upa, col=density),alpha=0.6) +
  #geom_point(data= Rsig, aes(x=ObsChange, y=SimChange),fill="transparent", color="gray25",shape=21,size=2)+
  geom_abline(slope=1,intercept=0,col="gray25",lwd=1.5,alpha=.5)+
  scale_color_viridis(option="F")+
  # 
  # annotate("label", x=-40, y=50, label= paste0("N = ",ln2),size=5)+
  # annotate("label", x=40, y=50, label= paste0("N = ",lp1),size=5)+
  # annotate("label", x=40, y=-50, label= paste0("N = ",lp2),size=5)+
  # annotate("label", x=-40, y=-50, label= paste0("N = ",ln1),size=5)+
  
  annotate("label", x=-32, y=36, label= paste0("R2 = ",R2lf,"\nRMSE = ",rmself,"%"),size=5)+
  scale_size(range = c(2, 8), trans="sqrt",name= expression(paste("Upstream area ", (km^2),
                                                                  sep = " ")),
             breaks=c(101,1000,10000,100000,500000), labels=c("100","1000", "10 000", "100 000", "500 000"))+
  scale_x_continuous(name="Observed change (% of mean annual maximum discharge per decade)",breaks = seq(-50,50,by=10), limits = lims)+
  scale_y_continuous(name="Simulated change (% of mean annual maximum discharge per decade)",breaks = seq(-100,100,by=10), limits = lims)+
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
  ggtitle(paste0("HERA - Observed and simulated trends of\nannual maximum discharge (1951-2020) (n = ",lt,")"))# adjust transparency range if needed

ps

ggsave("Revisions/scatter2_correlation_HERAvsobs_19512020_AMAX.jpg", ps, width=20, height=20, units=c("cm"),dpi=300)




#Now spatialized map
##############################################################################
#interpolating the sen-slopes via kriging
##############################################################################

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


#in order to perform kriging we need a different coordinate system than long-lat and a grid

laea_proj="+proj=laea +lat_0=52 +lon_0=10 +x_0=4321000 +y_0=3210000 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs +type=crs"
lambert_proj="+proj=laea +lat_0=55 +lon_0=20 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=km +no_defs"

spatial_krig = SpatialPointsDataFrame(coords=Recaplus[,c("Var1", "Var2")], data=Recaplus[,c("sen_obs","pval_obs")],proj4string=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"))

spatial_krig_trafo = spTransform(spatial_krig,CRS(laea_proj))

# spatial_krig = SpatialPointsDataFrame(coords=data_europe[,c("LON", "LAT")], data=data_europe[,c("pct_of_mean_per_decade","mk_pval")],proj4string=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"))
# spatial_krig_trafo = spTransform(spatial_krig,CRS("+proj=laea +lat_0=55 +lon_0=20 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=km +no_defs"))
# data_europe$lambert_x = spatial_krig_trafo@coords[,1]
# data_europe$lambert_y = spatial_krig_trafo@coords[,2]

Recaplus$laea_x = spatial_krig_trafo@coords[,1]
Recaplus$laea_y = spatial_krig_trafo@coords[,2]

europe_krig = SpatialPointsDataFrame(coords=europe_frame[,c("long", "lat")], data=europe_frame[,c(1:4)],proj4string=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"))
europe_krig_trafo = spTransform(europe_krig,CRS(laea_proj))

europe_frame$laea_x = europe_krig_trafo@coords[,1]
europe_frame$laea_y = europe_krig_trafo@coords[,2]


#cover of europe for spatial grid
europe_cover = data.frame(x=c(-2700,-2700,2000,2000),y=c(-2500,2000,2000,-2500),data=rep("data",4))
#another europe cover
lims=c(xmin=2665274, xmax=6528705, ymin=1464363, ymax=5360738)
europe_cover = data.frame(x=c(2000000,2000000,7000000,7000000),y=c(1000000,6000000,6000000,1000000),data=rep("data",4))

coordinates(europe_cover) =~x+y
my_big_grid = makegrid(europe_cover,cellsize=20000)
keep_point = rep(0,nrow(my_big_grid))
relevant_countries = unique(europe_frame$plot_group)
relevant_countries <- relevant_countries[-grep("Rus|Aze|Turk|Ukr|Isra", relevant_countries, ignore.case = TRUE)]

#for the grid, we check which points are contained in any country
for(country in relevant_countries)  {
  tmp_country = europe_frame[which(europe_frame$plot_group==country),]
  pt_in_country_tmp = point.in.polygon(point.x=my_big_grid[,1], point.y=my_big_grid[,2], pol.x=tmp_country[,"laea_x"], pol.y=tmp_country[,"laea_y"], mode.checked=FALSE)
  if(length(which(pt_in_country_tmp==1))>0) {
    keep_point[which(pt_in_country_tmp==1)] = 1
  }
  print(which(country==relevant_countries)/length(relevant_countries))
}
hist(keep_point)
my_grid = my_big_grid[which(keep_point==1),]
names(my_grid) = c("lambert_x","lambert_y")

grid_krig = SpatialPointsDataFrame(coords=my_grid[,c(1,2)], data=my_grid[,c(1,2)],proj4string=CRS(laea_proj))

#grid_krig = SpatialPointsDataFrame(coords=my_grid[,c(1,2)], data=my_grid[,c(1,2)],proj4string=CRS("+proj=laea +lat_0=52 +lon_0=10 +x_0=4321000 +y_0=3210000 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs +type=crs"))

#kriging
library("automap")
kriging_result <- autoKrige(sen_obs~1, input_data=spatial_krig_trafo, new_data = grid_krig)


# Convert kriging results to a data frame
kriging_data <- as.data.frame(kriging_result$krige_output)

# Extract the coordinates from the SpatialPointsDataFrame
coords <- coordinates(kriging_result$krige_output)

# Combine the coordinates with the kriging data
kriging_data$lon <- coords[, 1]
kriging_data$lat <- coords[, 2]


library(dplyr)
# Transform back to longlat projection if necessary
kriging_data_out = SpatialPointsDataFrame(coords=kriging_data[,c("lon", "lat")], data=kriging_data[,c(3:5)],proj4string=CRS(laea_proj))
kriging_data_longlat <- spTransform(kriging_data_out, CRS("+proj=longlat +datum=WGS84"))

kriging_data_longlat=as.data.frame(kriging_data_out)


# Create a ggplot map of the kriging results
m3<-kriging_data_longlat
m4 <- m3 %>%
  # create a new variable from count
  mutate(countv1=cut(var1.pred, breaks=c(-24,-12,-5,-2, 0,2, 5, 12),
                     labels=c( "-24 - -12", "-12 - -5", "-5 - -2","-2 - 0","0 - 2", "0 - 5", "5 - 12"))) %>%
  # change level order
  mutate(countv1=factor(as.character(countv1), levels=rev(levels(countv1))))


outll=outletopen(hydroDir,outletname)
cord.dec=outll[,c(2,3)]
cord.dec = SpatialPoints(cord.dec, proj4string=CRS("+proj=longlat"))
cord.UTM <- spTransform(cord.dec, CRS("+init=epsg:3035"))
nco=cord.UTM@coords
coord.lamb=spTransform(cord.dec, CRS(laea_proj))
nci=coord.lamb@coords

coord.laea=spTransform(cord.dec, CRS(laea_proj))
nci=coord.laea@coords

world <- ne_countries(scale = "medium", returnclass = "sf")
Europe <- world[which(world$continent == "Europe"),]
e2=st_transform(Europe,  crs=3035)
w2=st_transform(world,  crs=3035)
basemap=w2
basemap=st_transform(basemap, CRS(laea_proj))

Recapoint <- st_as_sf(Recaplus, coords = c("Var1", "Var2"), crs = 4326)
Recapoint <- st_transform(Recapoint, crs = 3035)
mp <- Recapoint %>%
  # create a new variable from count
  mutate(countv1=cut(sen_obs, breaks=c(-100,-24,-12,-5,-2, 0,2, 5, 12,100),
                     labels=c("< -24", "-24 - -12", "-12 - -5", "-5 - -2","-2 - 0","0 - 2", "0 - 5", "5 - 12","> 12"))) %>%
  # change level order
  mutate(countv1=factor(as.character(countv1), levels=rev(levels(countv1))))


tsize=12
osize=10

paletc=c(hcl.colors(9, palette = "RdYlBu", alpha = NULL, rev = T, fixup = TRUE))
paletf=c(hcl.colors(7, palette = "RdYlBu", alpha = NULL, rev = T, fixup = TRUE))
pl3<-ggplot(basemap) +
  geom_sf(fill="white")+
  geom_tile(data=m4, aes(x = lon, y = lat, fill = countv1)) +
  #metR::geom_contour_fill(data=kriging_data_longlat, aes(x = lon, y = lat, z = var1.pred),binwidth=2,alpha=1) +
  geom_sf(fill=NA, color="grey") +
  geom_sf(data=mp,aes(geometry=geometry,color=countv1),size=1.5,alpha=.5,shape=16,stroke=0.5)+ 
  geom_sf(data=mp,aes(geometry=geometry),color="black",size=1.5,alpha=.5,shape=21,stroke=0.5)+ 
  #scale_fill_distiller(palette = "RdYlBu",limits=c(-12,12),oob = scales::squish,guide="coloursteps",direction=1)+
  #scale_fill_gradientn(colors=paletx,limits=c(-12,12),oob = scales::squish) +
  scale_fill_manual(values=paletf, na.value = "grey90")+
  scale_color_manual(values=paletc, na.value = "grey90")+
  coord_sf(xlim = c(min(nci[,1]),max(nci[,1])), ylim = c(min(nci[,2]),max(nci[,2])))+
  labs(x="Longitude", y = "Latitude")+
  # guides(fill = guide_coloursteps(barwidth = 1, barheight = 10))+
  labs(fill = "change in mean annual flood per decade (%)") +
  guides( color= "none",
          fill = guide_legend(theme = theme(
            legend.title.position = "right")))+
  theme(axis.title=element_text(size=tsize),
        panel.background = element_rect(fill = "aliceblue", colour = "grey1"),
        panel.border = element_rect(linetype = "solid", fill = NA, colour="black"),
        legend.title = element_text(angle = 90 , size=tsize),
        legend.text = element_text(size=osize),
        legend.key.height= unit(1.6, 'cm'),
        legend.key.width= unit(.2, 'cm'),
        legend.position = "right",
        legend.text.position = "left",
        panel.grid.major = element_line(colour = "grey70"),
        panel.grid.minor = element_line(colour = "grey90"),
        legend.key = element_rect(fill = "transparent", colour = "transparent"),
        legend.key.size = unit(.8, "cm"))+
  ggtitle("Flood trends (observed) (1960-2010)")
pl3

ggsave("Revisions/Obs_floodtrends_1960-2010.jpg", pl3, width=20, height=20, units=c("cm"),dpi=300)






setwd("~/LFRuns_utils/TrendAnalysis")
#Plotting function
Diff.plots.points=function(basemap,sppoints, lims=c(-50,50),trans="identity", brk=NULL, scale="diverging",title=" ",name= " "){
  if (scale=="diverging")  palet=c(hcl.colors(8, palette = "RdYlBu", alpha = NULL, rev = TRUE, fixup = TRUE))
  if (scale=="sequential")  palet=c(hcl.colors(8, palette = "YlGnBu", alpha = NULL, rev = TRUE, fixup = TRUE))
  ggplot(basemap) +
    geom_sf(fill="gray85")+
    geom_sf(fill=NA, color="grey") +
    geom_sf(data=sppoints,aes(color=diffplot,geometry=geometry),alpha=1,size=0.1,stroke=0,shape=15)+ 
    coord_sf(xlim = c(min(nco[,1]),max(nco[,1])), ylim = c(min(nco[,2]),max(nco[,2])))+
    scale_color_gradientn(
      colors=palet, limits=lims,oob = scales::squish,
      name=name, trans=trans,breaks=brk)   +
    labs(x="Longitude", y = "Latitude")+
    guides(fill = guide_colourbar(barwidth = 15, barheight = .8))+
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
}
Diff.plots.catch=function(basemap,spcatch, lims=c(-50,50),trans="identity", brk=NULL, scale="diverging",title=" ",name= " "){
  
  if (scale=="diverging")  palet=c(hcl.colors(8, palette = "RdYlBu", alpha = NULL, rev = FALSE, fixup = TRUE))
  if (scale=="sequential")  palet=c(hcl.colors(8, palette = "YlGnBu", alpha = NULL, rev = TRUE, fixup = TRUE))
  ggplot(basemap) +
    geom_sf(fill="gray85")+
    geom_sf(data=spcatch,aes(fill=diffplot,geometry=geometry),color="transparent",alpha=1,size=0.1,stroke=0,shape=15)+ 
    geom_sf(fill=NA, color="grey") +
    coord_sf(xlim = c(min(nco[,1]),max(nco[,1])), ylim = c(min(nco[,2]),max(nco[,2])))+
    scale_fill_gradientn(
      colors=palet, limits=lims, breaks=brk, trans=trans, oob = scales::squish,name=name)   +
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
    ggtitle (title)
  
}


#observed
load(paste0(hydroDir,"/TrendAnalysis/ResValid_flood_obs_1950_2020.Rdata"))
dataObs=Results

#simulated
load(paste0(hydroDir,"/TrendAnalysis/ResValid_flood_sim_1950_2020.Rdata"))
dataSim=Results


plot(dataObs$catrest$catlist,dataSim$catrest$catlist)
#add latitude and longitude to input
ValidSf=read.csv(file=paste0(valid_path,"Stations/Stations_Validation.csv"))[,-1]
ValidSY=ValidSf[which(ValidSf$removal!="YES"),]
#compute Rls for modelled and observed

dataObs$catrest

#use station_idplus to get the correponding station

Catnames=data.frame(dataObs$catrest,station_idplus$X1)

RLObs=data.frame(dataObs$RetLevGPD,unikout=as.numeric(station_ido$X2))
RLSim=data.frame(dataSim$RetLevGPD,unikout=station_ids$X1)
order(RLObs$unikout)
RLSim=RLSim[order(RLSim$unikout),]
order(RLSim$unikout)
#generate the changes
period=c(1965,2005)
years=c(period[1]:period[2])
valuenames=paste0("X",years)

datacol=names(RLObs)
valcol=match(valuenames,datacol)
datatwin=RLObs
valcol2=valcol
datatwin=as.data.frame(datatwin)

cref=paste0("X",period[1])
crefloc=match(cref,dtc)
finalperiod=paste0("X",period[2])
colsel=match(finalperiod,datacol)

palet=c(hcl.colors(11, palette = "RdYlBu", alpha = NULL, rev = F, fixup = TRUE))
# title=paste0("Relative change in 10-years ",haz," Return Level between ", period[1], " and ", period[2])
legend="Relative Change (%)"
tmpval=(datatwin[,valcol2])/(datatwin[,crefloc])*100-100
datatwin[,valcol2]=tmpval
RLObs_ch=datatwin


datacol=names(RLSim)
valcol=match(valuenames,datacol)
datatwin=RLSim
valcol2=valcol
datatwin=as.data.frame(datatwin)

cref=paste0("X",period[1])
crefloc=match(cref,dtc)
finalperiod=paste0("X",period[2])
colsel=match(finalperiod,datacol)

tmpval=(datatwin[,valcol2])/(datatwin[,crefloc])*100-100
datatwin[,valcol2]=tmpval
RLSim_ch=datatwin



br=c(-100,-80,-60,-40,-20,0,20,40,60,80,100)
limi=c(-100,100)
trans=scales::modulus_trans(.8)
colNA="green"

data_compObsSim=inner_join(RLObs,RLSim,by="unikout")
plot(data_compObsSim$X2005.y,data_compObsSim$X2005.x,log="xy")
abline(a=0,b=1,col=2,lwd=2)

data_compObsSim=inner_join(RLObs_ch,RLSim_ch,by="unikout")


#make a plot with four squares showing the success areas
plot(data_compObsSim$X2005.y,data_compObsSim$X2005.x)
abline(a=0,b=1,col=2,lwd=2)
model <- cor.test(data_compObsSim$X2005.y,data_compObsSim$X2005.x)
model


data_compObsSim$diffplot=data_compObsSim$X2005.y - data_compObsSim$X2005.x

plot(data_compObsSim$diffplot)

#legend="Kendall tau"
mkoa=c()
mksa=c()
corv=c()
for (it in 1:length(RLSim[,1])){
  if (it%%1000==0) print(it)
  mks=NA
  mkt=NA
  miniTSim=as.numeric(RLSim[it,-71])
  miniTObs=as.numeric(RLObs[it,-71])
  Kor=cor(miniTSim,miniTObs)
  if (!is.na(miniTS[2])){
    #mk=MannKendall(miniTS[-1])
    mko1=mmkh(miniTObs[-1],ci=0.95)
    mks1=mmkh(miniTSim[-1],ci=0.95)
    mko2=data.frame(sens=mko1[7],sl=mko1[2])
    mks2=data.frame(sens=mks1[7],sl=mks1[2])
    #compute trend as well with sen.slope
  }else{
    mko2=data.frame(sens=NA,sl=NA)
    mks2=data.frame(sens=NA,sl=NA)

  }
  mkoa=rbind(mkoa,mko2)
  mksa=rbind(mksa,mks2)
  corv=c(corv,Kor)
}

hist(corv,breaks=50)
mkoa$siglvl=0
mkoa$siglvl[which(mkoa$sl<=0.1)]=1
length(mkoa$sens[which(mkoa$siglvl==1)])

mksa$siglvl=0
mksa$siglvl[which(mksa$sl<=0.1)]=1
length(mksa$sens[which(mksa$siglvl==1)])

Recap=cbind(mkoa,mksa,corv)
colnames(Recap)=c("senslope_o","pval_o","siglevel_o","senslope_s","pval_s","siglevel_s","correlation")
Recap$doublesig=0
Recap$doublesig[which(Recap$siglevel_o==1 & Recap$siglevel_s==1)]=1
plot(Recap$senslope_o,Recap$senslope_s,col=Recap$siglevel_o+3)


#Add the 1965-2005 changes to Recap
Recap$ObsChange=data_compObsSim$X2005.x
Recap$SimChange=data_compObsSim$X2005.y
Recap$alpha=1/(round(Recap$pval_o,2)+1e-2)/100
plot(Recap$ObsChange,Recap$SimChange)
abline(a=0,b=1,col=2,lwd=2)


#agreement in sign of change
lp1=length(which(Recap$ObsChange>0 & Recap$SimChange>0))
ln1=length(which(Recap$ObsChange<=0 & Recap$SimChange<=0))

lp2=length(which(Recap$ObsChange>0 & Recap$SimChange<=0))
ln2=length(which(Recap$ObsChange<=0 & Recap$SimChange>0))

lt=lp1+lp2+ln1+ln2
truesign=(lp1+ln1)/lt


Rsig=Recap[which(Recap$siglevel_o==1),]
lp1=length(which(Rsig$ObsChange>0 & Rsig$SimChange>0))
ln1=length(which(Rsig$ObsChange<=0 & Rsig$SimChange<=0))

lp2=length(which(Rsig$ObsChange>0 & Rsig$SimChange<=0))
ln2=length(which(Rsig$ObsChange<=0 & Rsig$SimChange>0))
median(Rsig$correlation)
r=cor(Recap$ObsChange,Recap$SimChange)
r2=round(r^2,3)
# Assuming 'data' is your data frame, 'x' and 'y' are your variables, and 'z' is the variable for transparency
ps<-ggplot() +
  annotate("rect", xmin=-Inf, xmax=0, ymin=-Inf, ymax=0, fill="lightskyblue", alpha=0.4) +

  annotate("rect", xmin=0, xmax=Inf, ymin=-Inf, ymax=0, fill="lightsalmon", alpha=0.4) +

  annotate("rect", xmin=0, xmax=Inf, ymin=0, ymax=Inf, fill="lightskyblue", alpha=0.4) +

  annotate("rect", xmin=-Inf, xmax=0, ymin=0, ymax=Inf, fill="lightsalmon", alpha=0.4) +

  geom_point(data=Recap, aes(x=ObsChange, y=SimChange, color=correlation),alpha=0.7,size=2) +
  geom_point(data= Rsig, aes(x=ObsChange, y=SimChange),fill="transparent", color="gray25",shape=21,size=2)+
  geom_abline(slope=1,intercept=0,col="gray25",lwd=1.5,alpha=.5)+

  annotate("label", x=-90, y=100, label= paste0("N = ",ln2),size=5)+
  annotate("label", x=90, y=100, label= paste0("N = ",lp1),size=5)+
  annotate("label", x=90, y=-100, label= paste0("N = ",lp2),size=5)+
  annotate("label", x=-90, y=-100, label= paste0("N = ",ln1),size=5)+
  
  scale_x_continuous(name="Observed change (%)",breaks = seq(-100,100,by=25), limits = c(-100,100))+
  scale_y_continuous(name="Simulated change (%)",breaks = seq(-100,100,by=25), limits = c(-100,100))+
  scale_color_gradientn(
    colors=palet, limits=c(-1,1),oob = scales::squish,
    name="temporal correlation",breaks=seq(-1,1, by=0.2))+
  # scale_alpha_continuous(name="trend significance",range = c(0.5, 1))+
  theme(axis.title=element_text(size=tsize),
        panel.background = element_rect(fill = "white", colour = "grey1"),
        panel.border = element_rect(linetype = "solid", fill = NA, colour="black"),
        legend.title = element_text(size=tsize),
        legend.text = element_text(size=osize),
        legend.position = "right",
        panel.grid.major = element_line(colour = "grey70"),
        panel.grid.minor = element_line(colour = "grey90"),
        legend.key = element_rect(fill = "transparent", colour = "transparent"),
        legend.key.size = unit(.8, "cm"))+
  ggtitle("Observed and simulated relative changes in 10y-RL flood between 1965 and 2005")# adjust transparency range if needed

ggsave("plots/scatter_correlation_floodtrends.jpg", ps, width=20, height=15, units=c("cm"),dpi=1500)



#Histogram of correlation like in Validation
Recap$bins <-(cut(Recap$correlation, breaks = 20))
Recag=Recap %>%
  group_by(bins) %>%
  summarise(mean = mean(correlation), n = length(correlation))
plot(Recap$correlation)
colors <- c(hcl.colors(19, palette = "RdYlBu", alpha = NULL, rev = F, fixup = TRUE))
length(Recap$correlation[which(Recap$correlation>0)])
length(Recap$correlation[which(Recap$correlation>0.5)])/length(Recap$correlation)
length(Recap$correlation[which(Recap$SimChange>0.5 & Recap$siglevel_o==1)])
length(Recap$correlation[which(Recap$siglevel_o==1)])
name="Temporal correlation"

p<-ggplot(data=Recap, aes(x=correlation)) + 
  geom_histogram(fill=colors,color="grey20",alpha=0.9,lwd=1,bins=20) +
  #geom_bar( data=Recap, aes(x=binN, fill=binN, group=bins),color="grey20",alpha=0.9,lwd=1) +
  geom_point(x=median(Recap$correlation,na.rm=T),y=0,pch=21,size=4,stroke=2,fill="royalblue")+
  geom_vline(xintercept=1,col=2,lwd=2)+
  # scale_fill_gradientn(
  #   colors=palet,oob = scales::squish,
  #   name="temporal correlation",breaks=seq(-1,1, by=0.2))+
  scale_y_continuous(name="Number of stations")+
  scale_x_continuous(name=name,breaks=seq(-1,1,by=0.2), limit=c(-1,1)) +
  #scale_x_discrete(labels= c("< -0.41","-0.41-0.2", "0.2-0.5","0.5-0.7","0.7-0.8", ">0.8"), name="KGE")+
  theme(axis.title=element_text(size=16, face="bold"),
        axis.text = element_text(size=13),
        panel.background = element_rect(fill = "white", colour = "white"),
        panel.grid = element_blank(),
        panel.border = element_rect(linetype = "solid", fill = NA, colour="black"),
        legend.title = element_text(size=14),
        legend.text = element_text(size=12),
        legend.position = "none",
        panel.grid.major = element_line(colour = "grey80"),
        panel.grid.minor.x = element_line(colour = "grey90",linetype="dashed"),
        legend.key = element_rect(fill = "transparent", colour = "transparent"),
        legend.key.size = unit(.8, "cm"))

p
ggsave("plots/histo_correlation_floodtrends.jpg", p, width=20, height=15, units=c("cm"),dpi=1500)
#TP TN etc plots


data_plot=data.frame(loc=data_compObsSim$unikout,var=Recap$correlation, sigo=Recap$siglevel_o)
data_plot=inner_join(data_plot,ValidSY,by=c("loc"="V1"))
ppl <- st_as_sf(data_plot, coords = c("Var1", "Var2"), crs = 4326)
ppl <- st_transform(ppl, crs = 3035)

ppl2=ppl[which(ppl$sigo==1),]
br=c(-1,-0.75,-0.5,-0.25,0,0.25,0.50,0.75,1)
lims=c(-1,1)
tsize=12
osize=12
title="Correlation between simulated and observed 10y-RL flood trends between (1965-2005)"
palet=c(hcl.colors(8, palette = "RdYlBu", alpha = NULL, rev = F, fixup = TRUE))
plot(ppl$var)

pm<-ggplot(basemap) +
  geom_sf(fill="gray85")+
  geom_sf(fill=NA, color="grey") +
  geom_sf(data=ppl,aes(color=var,geometry=geometry,size=UpA),alpha=.9,shape=19,stroke=0)+ 
  geom_sf(data=ppl2,aes(geometry=geometry,size=UpA),color="black",alpha=.9,shape=21,stroke=0.1,fill=NA)+ 
  coord_sf(xlim = c(min(nco[,1]),max(nco[,1])), ylim = c(min(nco[,2]),max(nco[,2])))+
  scale_size(range = c(1, 4), trans="sqrt",name= expression(paste("Upstream area ", (km^2),
                                                                  sep = " ")),
             breaks=c(101,1000,10000,100000,500000), labels=c("100","1000", "10 000", "100 000", "500 000"))+
  scale_color_gradientn(
    colors=palet,oob = scales::squish,
    name="r",breaks=br, labels=br)   +
  labs(x="Longitude", y = "Latitude")+
  guides(fill = guide_colourbar(barwidth = 1, barheight = 12))+
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

ggsave("plots/map_correlation_floodtrends.jpg", pm, width=20, height=15, units=c("cm"),dpi=1500)




t1="Comparison of changes in 10-years RL flood between simulated and observed data"
n1="Cal - Uncal \n relative difference (%)"
dp1=Diff.plots.points(basemap,data_diffCalUncal,c(-50,50), brk=c(-50,-25,0,25,50),scale="diverging",t1,n1)
dp1
#Catchment level
data_changecal=Plot.change[[4]]
data_changeuncal=Plot.change2[[4]]
st_geometry(data_changecal)=NULL
st_geometry(data_changeuncal)=NULL

data_diffCalUncal=inner_join(data_changecal,data_changeuncal,by="HYBAS_ID")
data_diffCalUncal$diffplot=data_diffCalUncal$Rchange.mean.x - data_diffCalUncal$Rchange.mean.y
data_diffCalUncal=inner_join(hybasf7,data_diffCalUncal,by= "HYBAS_ID")

t2="Comparison of changes in 100-years RL flood between uncalibrated and calibrated runs"
n2="Cal - Uncal \n relative difference (%)"
dp2=Diff.plots.catch(basemap,data_diffCalUncal,c(-50,50), brk=c(-50,-25,0,25,50),scale="diverging",t2,n2)
dp2
#load the simulated discharge

#compute the changes, plot the changes, save the changes, compare






cord.dec=ValidSY[,c(1,2)]
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
palet=c(hcl.colors(9, palette = "viridis", alpha = NULL, rev = T, fixup = TRUE))
#palet = c('#f61f0f','#f4cccc','#f6ebeb','#ffffff','#e4edf5','#91bfdb','#1c87d8'))

ppl <- st_as_sf(ValidSY, coords = c("Var1", "Var2"), crs = 4326)
ppl <- st_transform(ppl, crs = 3035)

ggplot(basemap) +
  geom_sf(fill="gray95", color=NA) +
  geom_sf(data=ppl,aes(geometry=geometry,fill=Rlen,size=UpA),color="transparent",alpha=.9,shape=21,stroke=0)+ 
  geom_sf(fill=NA, color="grey20") +
  #geom_sf(data=ppl,aes(geometry=geometry,size=UpA),col="grey",alpha=1,stroke=0.1,shape=1)+ 
  scale_x_continuous(breaks=seq(-30,40, by=5)) +
  scale_size(range = c(1, 4), trans="sqrt",name= expression(paste("Upstream area ", (km^2),
                                                                 sep = " ")),
             breaks=c(101,1000,10000,100000,500000), labels=c("100","1000", "10 000", "100 000", "500 000"))+
  scale_fill_gradientn(
    colors=palet,oob = scales::squish, name="record length (years)", trans="sqrt",
    breaks=c(365,1825,3650,7300,14600,21900), labels=c(1,5,10,20,40,60)) +
  coord_sf(xlim = c(min(nco[,1]),max(nco[,1])), ylim = c(min(nco[,2]),max(nco[,2])))+
  labs(x="Longitude", y = "Latitude") +
  guides(fill = guide_colourbar(barwidth = 20, barheight = 0.5,reverse=F),
         size= guide_legend(override.aes = list(fill = "grey50")))+
  theme(axis.title=element_text(size=tsize),
        panel.background = element_rect(fill = "aliceblue", colour = "grey1"),
        panel.border = element_rect(linetype = "solid", fill = NA, colour="black"),
        legend.title = element_text(size=tsize),
        legend.text = element_text(size=osize),
        legend.position = "bottom",
        legend.box = "vertical",
        panel.grid.major = element_line(colour = "grey85",linetype="dashed"),
        panel.grid.minor = element_line(colour = "grey90"),
        legend.key = element_rect(fill = "transparent", colour = "transparent"),
        legend.key.size = unit(.8, "cm"))
ggsave("Plots/ValidStations.jpg", width=20, height=24, units=c("cm"),dpi=1500)



rat1=abs((dat$sim-dat$obs))/dat$obs
rat2=abs((dat$obs-dat$sim))/dat$sim
rmv=which(rat1>3 | rat2>3)
dat$rat1=rat1
dat$rat2=rat2
dat$flag1=NA
dat$flag1[which(dat$rat1>3 | dat$rat2>3)]=1
dat$flag1[which(dat$obs>10 & dat$flag1==1)]=2
dat$flag1[which(dat$rat1>6 | dat$rat2>6)]=3

#loading Dominiks initial database
Q_stations <- st_read(paste0(valid_path,'Europe_daily_combined_2022_v2_WGS84_NUTS.shp')) 

#Writing result: merging dat and validSta
ValidS=inner_join(ValidSta,dat,by=c("V1"="Station_ID"))
ValidS=inner_join(ValidS,Q_stations,by = c("V1"="StationID"))

#reorganise the data
ValidSf=ValidS[,c(1,2,7,19,4,5,6,14,10,11,12,13,15,18,22)]

## Identify stations with low KGE ---------------
kgefile="out/EFAS_obs_kgeAY.csv"
kge=SpatialSkillPlot(ValidSf,"kge",kgefile)
ValidSf=kge[[1]][,-16]
#Isolate problematic stations stations:
ValidSf$flag2=NA
ValidSf$flag2[which(ValidSf$skill<(-0.0))]=1
ValidSf$removal = ""
ValidSf$comment = ""
#write.csv(rest,file="stations_manual_check.csv")

## Individual check of the remaining stations ------------------
#load remaining stations checked
Mancheck=read.csv(file="Stations/stations_manual_check.csv")
Mancheck=Mancheck[,-1]
manmatch=match(Mancheck$V1,ValidSf$V1)
Vcheck=Mancheck[c(25,26)]
colnames(Vcheck)
colnames(ValidSf)[c(19,20)]
ValidSf[manmatch,c(19,20)]=Vcheck

#Other manual Checks
#Gagnieres 6139070
ValidSf$comment[which(ValidSf$V1==6139070)]="wrong river"
ValidSf$removal[which(ValidSf$V1==6139070)]="YES"
#6233203
ValidSf$comment[which(ValidSf$V1==6233203)]="dubious observations"
ValidSf$removal[which(ValidSf$V1==6233203)]="YES"
#6118175
ValidSf$comment[which(ValidSf$V1==6118175)]="wrong location"
ValidSf$removal[which(ValidSf$V1==6118175)]="YES"
#6124440
ValidSf$comment[which(ValidSf$V1==6124440)]="wrong location"
ValidSf$removal[which(ValidSf$V1==6124440)]="YES"
#6340215
ValidSf$comment[which(ValidSf$V1==6340215)]="wrong location"
ValidSf$removal[which(ValidSf$V1==6340215)]="YES"

ValidSf$comment[which(ValidSf$flag1==2 & ValidSf$removal=="")]="mismatch between obs and sim Qmean"
ValidSf$removal[which(ValidSf$flag1==2 & ValidSf$removal=="")]="YES"

ValidSf$comment[which(ValidSf$flag1==3 & ValidSf$removal=="")]="mismatch between obs and sim Qmean"
ValidSf$removal[which(ValidSf$flag1==3 & ValidSf$removal=="")]="YES"

## Distance between "official gauges" and efas points -----------------------
points=ValidSf[,c(1,2,3)]
Vsfloc=st_as_sf(points, coords = c("Var1", "Var2"), crs = 4326)
Statloc=Q_stations
Sfloc=inner_join(points,Statloc,by=c("V1"="StationID"))
dist=c()
for (r in Vsfloc$V1){
  cat(paste0(r,"\n"))
  v1=Vsfloc[which(Vsfloc$V1==r),]
  v2=Statloc[which(Statloc$StationID==r),]
  oula=st_distance(v1,v2)
  oula=oula/1000
  dist=c(dist,oula)
}

ValidSf$distance=dist
ValidSf$flag2[which(ValidSf$distance>2.5 & ValidSf$flag2==1 & ValidSf$removal=="")]=2
ValidSf$removal[which(ValidSf$flag2==2)]="YES"
ValidSf$comment[which(ValidSf$flag2==2)]="too far from original station"

length(ValidSf$comment[which(ValidSf$comment=="too far from original station")])
length(ValidSf$comment[which(ValidSf$comment=="mismatch between obs and sim Qmean")])
konar=(ValidSf[which(ValidSf$removal=="YES"),])
length(which(ValidSf$removal=="YES"))
length(ValidSf$comment[which(ValidSf$csource=="EFAS")])

#Writing of the final file
#write.csv(ValidSf,file="Stations/Stations_Validation.csv")
StatCheck=ValidSf[which(ValidSf$flag2==1),]
#write.csv(StatCheck,file="Stations/Stations_lowKGE.csv")


# Part 2: Diagnostic plots------------------------------------

ValidSf=read.csv(file="Stations/Stations_Validation.csv")[,-1]
ValidSY=ValidSf[which(ValidSf$removal!="YES"),]
length(ValidSY$comment[which(ValidSY$UpA<=250)])/2901
length(ValidSY$comment[which(ValidSY$skill>-0.41)])
median(ValidSY$skill)

vEFAS=ValidSY[which(ValidSY$csource=="EFAS"),]
median(vEFAS$skill)
#Map of upstream area
#Plot parameters

cord.dec=ValidSY[,c(1,2)]
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
palet=c(hcl.colors(9, palette = "viridis", alpha = NULL, rev = T, fixup = TRUE))
#palet = c('#f61f0f','#f4cccc','#f6ebeb','#ffffff','#e4edf5','#91bfdb','#1c87d8'))

ppl <- st_as_sf(ValidSY, coords = c("Var1", "Var2"), crs = 4326)
ppl <- st_transform(ppl, crs = 3035)

ggplot(basemap) +
  geom_sf(fill="gray95", color=NA) +
  geom_sf(data=ppl,aes(geometry=geometry,fill=Rlen,size=UpA),color="transparent",alpha=.9,shape=21,stroke=0)+ 
  geom_sf(fill=NA, color="grey20") +
  #geom_sf(data=ppl,aes(geometry=geometry,size=UpA),col="grey",alpha=1,stroke=0.1,shape=1)+ 
  scale_x_continuous(breaks=seq(-30,40, by=5)) +
  scale_size(range = c(1, 4), trans="sqrt",name= expression(paste("Upstream area ", (km^2),
                                                                 sep = " ")),
             breaks=c(101,1000,10000,100000,500000), labels=c("100","1000", "10 000", "100 000", "500 000"))+
  scale_fill_gradientn(
    colors=palet,oob = scales::squish, name="record length (years)", trans="sqrt",
    breaks=c(365,1825,3650,7300,14600,21900), labels=c(1,5,10,20,40,60)) +
  coord_sf(xlim = c(min(nco[,1]),max(nco[,1])), ylim = c(min(nco[,2]),max(nco[,2])))+
  labs(x="Longitude", y = "Latitude") +
  guides(fill = guide_colourbar(barwidth = 20, barheight = 0.5,reverse=F),
         size= guide_legend(override.aes = list(fill = "grey50")))+
  theme(axis.title=element_text(size=tsize),
        panel.background = element_rect(fill = "aliceblue", colour = "grey1"),
        panel.border = element_rect(linetype = "solid", fill = NA, colour="black"),
        legend.title = element_text(size=tsize),
        legend.text = element_text(size=osize),
        legend.position = "bottom",
        legend.box = "vertical",
        panel.grid.major = element_line(colour = "grey85",linetype="dashed"),
        panel.grid.minor = element_line(colour = "grey90"),
        legend.key = element_rect(fill = "transparent", colour = "transparent"),
        legend.key.size = unit(.8, "cm"))
ggsave("Plots/ValidStations.jpg", width=20, height=24, units=c("cm"),dpi=1500)

## Scatter plot of quantiles -------------
fileobs="out/obs_percentilesAY.csv"
filesim="out/EFAS_percentilesAY.csv"
obs=read.csv(paste0(valid_path,fileobs))
sim=read.csv(paste0(valid_path,filesim))
names(obs)[1] = names(sim)[1]="Station_ID"
names(obs)[c(2:100)] = names(sim)[c(2:100)]=paste('q',1:99,sep="_")

#quantile selection:
qtil=05
myq=paste0("q_",qtil)
if (qtil<10) qtil=paste0(0,qtil)
nobs=paste0("Q",qtil,"_obs")
nsim=paste0("Q",qtil,"_sim")
scplot=paste0("scatterplot_Q",qtil)
qloc=match(myq,colnames(obs))

  scale_x_log10(name=expression(paste("Upstream area ", (km^2),sep = " ")),
                breaks=c(100,1000,10000,100000), minor_breaks = log10_minor_break(),
                labels=c("100","1 000","10 000","100 000")) +
  theme(axis.title=element_text(size=16, face="bold"),
        axis.text = element_text(size=16),
        panel.background = element_rect(fill = "white", colour = "white"),
        panel.grid = element_blank(),
        panel.border = element_rect(linetype = "solid", fill = NA, colour="black"),
        legend.title = element_text(size=14),
        legend.text = element_text(size=12),
        panel.grid.major = element_line(colour = "grey80"),
        panel.grid.minor.x = element_line(colour = "grey90",linetype="dashed"),
        legend.key = element_rect(fill = "transparent", colour = "transparent"),
        legend.key.size = unit(.8, "cm"))+
  annotate("label", x=450000, y=500, label= paste0("n = ",length(ppl$UpA)),size=6)

p
ggsave("histo_stations_f.jpg", p, width=20, height=15, units=c("cm"),dpi=1500)

#Map of station used in validation (number of years on record) -Decommisionned

palet=c(hcl.colors(9, palette = "viridis", alpha = NULL, rev = T, fixup = TRUE))

ggplot(basemap) +
  geom_sf(fill="gray95", color=NA) +
  geom_sf(data=parpl,aes(geometry=geometry,col=years),alpha=.7,size=1.4,stroke=0,shape=16)+ 
  geom_sf(fill=NA, color="grey20") +
  geom_sf(data=parpef,aes(geometry=geometry),col="black",alpha=1,size=1.4,stroke=0.5,shape=1)+ 
  scale_x_continuous(breaks=seq(-30,40, by=5)) +
  scale_color_gradientn(
    colors=palet, name="length of record (Years)",breaks=c(0,10,20,30,40,50,60,70),oob = scales::squish)  +
  labs(x="Longitude", y = "Latitude")+
  coord_sf(xlim = c(min(nco[,1]),max(nco[,1])), ylim = c(min(nco[,2]),max(nco[,2])))+
  # scale_color_manual(values=c("EFAS" ="blue","FloodDrivers"="red","Both"="green"),name="River Gauges")   +
  labs(x="Longitude", y = "Latitude")+
  guides(color = guide_bins(barwidth = 10, barheight = 2,reverse=F,override.aes = list(size = 5, shape=15)))+
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

ggsave("Validation_gauges.jpg", width=20, height=20, units=c("cm"),dpi=1500)


### Boxplot of performance by decade -----
#groups= "1950-1960","1960-1970","1970-1980","1980-1990","1990-2000","2000-2010","2010-2020"

kgefile="out/EFAS_obs_kge.csv"
skill=read.csv(paste0(valid_path,kgefile))
if (length(skill[1,])>2){
  years=c(1950:2020)
  skill=skill[-73]
  names(skill)[c(2:72)]=years
}
names(skill)[1]="Station_ID"


#Loading all txt files with matching locations: I want to have all matches ==> I can get the area
#load matching coordinates
Sloc=read.table(paste0(valid_path,"out/Stations_locations_EFAS_grid_1950.txt"),sep=",")
yrlist=c(1951:2020)
chd=seq(1960,2020,by=10)
Slolos=list()
ix=0
lslo=rep(NA,length(chd))
for (yi in yrlist){
  print(yi)
  if (length(which(yi==chd))>0){
    ix=ix+1
    Sloc_final=Sloc[which(Sloc$V2!=0),]
    onlyV=which(!is.na(match(Sloc_final$V1,ValidSY$V1)))
    Sloc_fx=Sloc_final[onlyV,]
    Slolos[[ix]]=Sloc_fx
    lslo[ix]=length(Sloc_fx$V1)
    Sloc=read.table(paste0(valid_path,"out/Stations_locations_EFAS_grid_",yi,".txt"),sep=",")
  }else{
    Slocy=read.table(paste0(valid_path,"out/Stations_locations_EFAS_grid_",yi,".txt"),sep=",")
    Slocrep=Slocy[which(Sloc$V2==0),]
    Sloc[which(Sloc$V2==0),]=Slocrep
  }
}

Sloc_final=Sloc[which(Sloc$V2!=0),]
Sloc_final$csource="SpatialQMatch"

#Only keep valid stations
skillm<-melt(skill,id.vars ="Station_ID")
skillm$year=as.numeric(as.character(skillm$variable))
skillm$decade=1950
skillm$decade[which(skillm$year>1959 & skillm$year<1970)]=1960
skillm$decade[which(skillm$year>1969 & skillm$year<1980)]=1970
skillm$decade[which(skillm$year>1979 & skillm$year<1990)]=1980
skillm$decade[which(skillm$year>1989 & skillm$year<2000)]=1990
skillm$decade[which(skillm$year>1999 & skillm$year<2010)]=2000
skillm$decade[which(skillm$year>2009)]=2010
#average
skillm$value[which(is.infinite(skillm$value))]=NA
dat_skill=aggregate(list(skill=skillm$value),
                      by = list(Station_ID=skillm$Station_ID,decade=skillm$decade),
                      FUN = function(x) c(mean= mean(x,na.rm=T),median=median(x,na.rm=T),max=max(x,na.rm=T)))
dat_skill <- do.call(data.frame, dat_skill)
onlyV=which(!is.na(match(dat_skill$Station_ID,ValidSY$V1)))
dat_skill=dat_skill[onlyV,]

plot(dat_skill$skill.median)

palet=c(hcl.colors(9, palette = "Blues", alpha = NULL, rev = T, fixup = TRUE))

meds <- c(by(dat_skill$skill.median, dat_skill$decade, median,na.rm=T))
qshit <- c(by(dat_skill$skill.median, dat_skill$decade, quantile ,na.rm=T))
merdecol=match(dat_skill$decade,names(meds))
dat_skill$col=meds[merdecol]

labelsX=c("1950-1960","1960-1970","1970-1980","1980-1990","1990-2000","2000-2010","2010-2020")

p3<-ggplot(dat_skill, aes(x=factor(decade), y=skill.median)) +
  # labs(title=tit1,x="Warmimg levels (?C)", y = y1)+
  geom_boxplot(notch=F,position=position_dodge(.9),alpha=.8,aes(fill=col),linewidth=0.8,outlier.alpha = 0.1)+
  #geom_violin(aes(fill=col),scale="width",draw_quantiles = c(0.25, 0.5, 0.75),position=position_dodge(.9),alpha=0.7,lwd=1)+
  scale_fill_gradientn(
    colors=palet, n.breaks=6,limits=c(0.4,0.6)) +
  scale_y_continuous(limits = c(-0.5,1),name="KGE",breaks = seq(-1,1,by=0.5),minor_breaks = seq(-1,1,0.1))+
  scale_x_discrete(name="Decade",labels=labelsX)+
  geom_text(data=data.frame(), aes(x=names(meds), y=meds-0.05, label=round(lslo,2)), col='black', size=4,fontface="bold")+
  theme(axis.title=element_text(size=16, face="bold"),
        axis.text = element_text(size=12),
        panel.background = element_rect(fill = "white", colour = "white"),
        panel.grid = element_blank(),
        panel.border = element_rect(linetype = "solid", fill = NA, colour="black"),
        legend.title = element_text(size=14),
        legend.text = element_text(size=12),
        legend.position="none",
        panel.grid.major = element_line(colour = "grey60"),
        panel.grid.minor.y = element_line(colour = "grey80",linetype="dashed"),
        legend.key = element_rect(fill = "transparent", colour = "transparent"),
        legend.key.size = unit(.8, "cm"))
p3
ggsave("Plots/boxp_KGEtime.jpg", p3, width=20, height=15, units=c("cm"),dpi=1500)



### Boxplot of performances according to the presence or not of reservoirs ------------------------

#load reservoir ratio map
res_path<-("D:/tilloal/Documents/LFRuns_utils/data/reservoirs/")
outletname="res_ratio_European_01min.nc"
dir=res_path
ResData=resOpen(res_path,outletname,ValidSY)  
rsx=ResData[[2]]
length(rsx$res.ratio[which(rsx$res.ratio<0.5)])/length(rsx$res.ratio)
ValidSYp=ResData
#do two boxplots, one for res.ratio>0.5 and one for others

ValidSYp$res.group=0
ValidSYp$res.group[which(ValidSYp$res.ratio>=0.5 & ValidSYp$res.ratio<=1)]=1
ValidSYp$res.group[which(ValidSYp$res.ratio>1)]=2

meds <- c(by(ValidSYp$skill, ValidSYp$res.group, median))
len <- c(by(ValidSYp$skill, ValidSYp$res.group, length))
merdecol=match(ValidSYp$res.group,names(meds))
ValidSYp$col=meds[merdecol]
palet=c(hcl.colors(9, palette = "Blues", alpha = NULL, rev = T, fixup = TRUE))
p4<-ggplot(ValidSYp, aes(x=factor(res.group), y=skill)) +
  # scale_fill_manual(values = cols, name = "Socioeconomic \n Pathways")+
  # labs(title=tit1,x="Warmimg levels (?C)", y = y1)+
  geom_boxplot(notch=F,position=position_dodge(.9),alpha=.8,aes(fill=col),linewidth=0.8,outlier.alpha = 0.4)+
  scale_y_continuous(limits = c(-0.5,1),name="KGE",breaks = seq(-1,1,by=0.5),minor_breaks = seq(-1,1,0.1))+
  scale_x_discrete(labels=c("0" = "c < 0.5 \n (low reservoir impact)", "1" = "0.5 <= c <= 1  \n (medium reservoir impact)",
                            "2" = "c > 1 \n (high reservoir impact)"),name="c ratio (reservoir volume to mean annual streamflow) ")+
  scale_fill_gradientn(
    colors=palet, n.breaks=6,limits=c(0.2,0.8)) +
  geom_text(data=data.frame(), aes(x=names(meds), y=meds-0.05, label=len), col='black', size=4,fontface="bold")+
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
p4
ggsave("Plots/boxp_KGE_reservoirs.jpg", p4, width=20, height=15, units=c("cm"),dpi=1500)

## Performance on annual maxima and minima -----------------------
Q_data <- read.csv(paste0('Q_19502020.csv'), header = F)  # CSVs with observations
Q_sim <- read.csv(file="out/EFAS_19502020.csv")
Station_data_IDo <- as.numeric(as.vector(t(Q_data[1, -c(1)])))
Station_data_IDs <- as.vector(t(Q_sim[1, ]))
# kind of hit rate: 2 boxplots
hit_out=c()
rs_obs=c()
rs_sim=c()
for (id in 1:length(Station_data_IDo)){
  print(id)
  s=Station_data_IDo[id]
  # s=6574362
  ms=Q_data[-c(1,2,3,4),c(1,id+1)]
  ms$V2=as.numeric(ms$V2)
  ms$time=as.Date(ms$V2-1,origin="0000-01-01")
  ms[,2]=as.numeric(ms[,2])
  ms=data.frame(time=ms$time,data=ms[,2])
  #extract peak moments
  sl=which(Station_data_IDs==s)
  mss=data.frame(time=ms$time,data=Q_sim[-1,sl])
  iy=which(is.na(ms$data))
  hit=NA
  if (length(iy)>0){
    ms=ms[-iy,]
    mss=mss[-iy,]
  }
  if (length(ms$data>0)){
    AMAX_obs <- computeAnnualMaxima(ms)
    AMAX_sim <- computeAnnualMaxima(mss)
    if (length(AMAX_sim$annualMaxIndx)==length(AMAX_obs$annualMaxIndx)){
      AMAX_diff=AMAX_obs$annualMaxIndx-AMAX_sim$annualMaxIndx
      hit=length(which(abs(AMAX_diff)<8))/length(AMAX_diff)
    }
    data=data.frame(ms,sim=mss$data)
    names(data)=c("date","Q","Qs")
    res_obs=data.frame(Station=rep(s,length(AMAX_obs$annualMax)),AnMax=AMAX_obs$annualMax,AnMaxDate=AMAX_obs$annualMaxDate)
    res_sim=data.frame(Station=rep(s,length(AMAX_sim$annualMax)),AnMax=AMAX_sim$annualMax,AnMaxDate=AMAX_sim$annualMaxDate)
    # for (ev in 1:length(AMAX_obs$annualMaxIndx)){
    #   #print(ev)
    #   test=StatEpisode(data, date.deb = data$date[AMAX_obs$annualMaxIndx[ev]], win=30, qx.min = 100, mode = "FL",obs=TRUE)
    #   res=rbind(res,test)
    # }
    # save=c(Station_ID=s,rmean=mean(res[,2],na.rm=T),errmean=mean(res[,3],na.rm=T))
    rs_obs=rbind(rs_obs,res_obs)
    rs_sim=rbind(rs_sim,res_sim)
  }else{print("no observations")}
  hit_out=c(hit_out,hit)
}

hit_out=data.frame(hit_out,Station_data_IDo)
hitm=match(ValidSY$V1,hit_out$Station_data_IDo)
hit_out=hit_out[hitm,]
boxplot(hit_out$hit_out)
#1- mean on Anmax for every station
catagg_obs=aggregate(list(Q=rs_obs$AnMax),
                   by = list(Station=rs_obs$Station),
                   FUN = function(x) c(med=median(x,na.rm=T),mean=mean(x,na.rm=T),dev=sd(x,na.rm=T),len=length(x)))
catagg_obs <- do.call(data.frame, catagg_obs)

catagg_sim=aggregate(list(Q=rs_sim$AnMax),
                     by = list(Station=rs_sim$Station),
                     FUN = function(x) c(med=median(x,na.rm=T),mean=mean(x,na.rm=T),dev=sd(x,na.rm=T),len=length(x)))
catagg_sim <- do.call(data.frame, catagg_sim)

catagg_obssim=inner_join(catagg_obs,catagg_sim, by="Station",suffix = c("obs","sim"))

plot(catagg_obssim$Q.meanobs,catagg_obssim$Q.meansim)

#reataining only good matches
matV=match(ValidSY$V1,catagg_obssim$Station)

Catagg_f=catagg_obssim[matV,]
Catagg_f$dPeaks=(Catagg_f$Q.meansim-Catagg_f$Q.meanobs)/Catagg_f$Q.meanobs*100
plot(Catagg_f$Q.meanobs,Catagg_f$Q.meansim)


#Spatial plot, where is dpeak shit?
palet=c(hcl.colors(9, palette = "RdBu", alpha = NULL, rev = F, fixup = TRUE))
points=inner_join(ValidSY,Catagg_f,by=c("V1"="Station"))
parpl <- st_as_sf(points, coords = c("Var1", "Var2"), crs = 4326)
parpl <- st_transform(parpl, crs = 3035)
cord.dec=points[,c(1,2)]
cord.dec = SpatialPoints(cord.dec, proj4string=CRS("+proj=longlat"))
cord.UTM <- spTransform(cord.dec, CRS("+init=epsg:3035"))
nco=cord.UTM@coords
world <- ne_countries(scale = "medium", returnclass = "sf")
Europe <- world[which(world$continent == "Europe"),]
e2=st_transform(Europe,  crs=3035)
w2=st_transform(world,  crs=3035)
tsize=12
osize=12
limi=c(-100,100)
basemap=w2
ggplot(basemap) +
  geom_sf(fill="gray95", color=NA) +
  geom_sf(data=parpl,aes(geometry=geometry,fill=dPeaks,size=UpA),color="transparent",alpha=.9,shape=21,stroke=0)+ 
  geom_sf(fill=NA, color="grey20") +
  geom_sf(data=parpl,aes(geometry=geometry,size=UpA),col="black",alpha=1,stroke=0.1,shape=1)+ 
  scale_x_continuous(breaks=seq(-30,40, by=5)) +
  scale_size(range = c(1, 3), trans="sqrt")+
  scale_fill_gradientn(
    colors=palet,n.breaks=5,oob = scales::squish, limits=limi, name="mean Qpeak difference") +
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
ggsave("Plots/Spatial_Qpeakdiff.jpg", width=20, height=15, units=c("cm"),dpi=1500)

median(Catagg_f$dPeaks)
boxplot(Catagg_f$dPeaks,ylim=c(-100,100))

matO=which(!is.na(match(rs_obs$Station,ValidSY$V1)))
matS=which(!is.na(match(rs_sim$Station,ValidSY$V1)))

rs_obsF=rs_obs[matO,]
rs_simF=rs_sim[matS,]

rs_obsF$yday=yday(rs_obsF$AnMaxDate)
seasonMax_obs=aggregate(list(Q=rs_obsF$yday),
                        by = list(Station=rs_obsF$Station),
                        FUN = function(x) c(season=season1(x),len=length(x)))
seasonMax_obs <- do.call(data.frame, seasonMax_obs)

rs_simF$yday=yday(rs_simF$AnMaxDate)
seasonMax_sim=aggregate(list(Q=rs_simF$yday),
                        by = list(Station=rs_simF$Station),
                        FUN = function(x) c(season=season1(x),len=length(x)))
seasonMax_sim <- do.call(data.frame, seasonMax_sim)

seasonMax_obs$dseason=seasonMax_obs$Q.season-seasonMax_sim$Q.season
seasonMax_obs$dseason[which(seasonMax_obs$dseason>182)]=seasonMax_obs$dseason[which(seasonMax_obs$dseason>182)]-365.25
seasonMax_obs$dseason[which(seasonMax_obs$dseason<=-182)]=365.25+seasonMax_obs$dseason[which(seasonMax_obs$dseason<=-182)]

hist(seasonMax_obs$dseason)
palet=c(hcl.colors(9, palette = "RdBu", alpha = NULL, rev = F, fixup = TRUE))
points=inner_join(ValidSY,seasonMax_obs,by=c("V1"="Station"))
parpl <- st_as_sf(points, coords = c("Var1", "Var2"), crs = 4326)
parpl <- st_transform(parpl, crs = 3035)
cord.dec=points[,c(1,2)]
cord.dec = SpatialPoints(cord.dec, proj4string=CRS("+proj=longlat"))
cord.UTM <- spTransform(cord.dec, CRS("+init=epsg:3035"))
nco=cord.UTM@coords
world <- ne_countries(scale = "medium", returnclass = "sf")
Europe <- world[which(world$continent == "Europe"),]
e2=st_transform(Europe,  crs=3035)
w2=st_transform(world,  crs=3035)
tsize=12
osize=12
limi=c(-100,100)
basemap=w2
ggplot(basemap) +
  geom_sf(fill="gray95", color=NA) +
  geom_sf(data=parpl,aes(geometry=geometry,fill=dseason,size=UpA),color="transparent",alpha=.9,shape=21,stroke=0)+ 
  geom_sf(fill=NA, color="grey20") +
  geom_sf(data=parpl,aes(geometry=geometry,size=UpA),col="black",alpha=1,stroke=0.1,shape=1)+ 
  scale_x_continuous(breaks=seq(-30,40, by=5)) +
  scale_size(range = c(1, 3), trans="sqrt")+
  scale_fill_gradientn(
    colors=palet,n.breaks=5,oob = scales::squish, limits=limi, name="Qpeak season difference") +
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
ggsave("Plots/Spatial_Speakdiff.jpg", width=20, height=15, units=c("cm"),dpi=1500)


#Same for low flows
hit_outl=c()
rsl_obs=c()
rsl_sim=c()
for (id in 1:length(Station_data_IDo)){
  print(id)
  s=Station_data_IDo[id]
  # s=6574362
  ms=Q_data[-c(1,2,3,4),c(1,id+1)]
  ms$V2=as.numeric(ms$V2)
  ms$time=as.Date(ms$V2-1,origin="0000-01-01")
  ms[,2]=as.numeric(ms[,2])
  ms=data.frame(time=ms$time,data=ms[,2])
  #extract peak moments
  sl=which(Station_data_IDs==s)
  mss=data.frame(time=ms$time,data=Q_sim[-1,sl])
  iy=which(is.na(ms$data))
  hit=NA
  if (length(iy)>0){
    ms=ms[-iy,]
    mss=mss[-iy,]
  }
  if (length(ms$data>0)){
    AMIN_obs <- computeAnnualMinima(ms)
    AMIN_sim <- computeAnnualMinima(mss)
    if (length(AMIN_sim$annualMinIndx)==length(AMIN_obs$annualMinIndx)){
      AMIN_diff=AMIN_obs$annualMinIndx-AMIN_sim$annualMinIndx
      hit=length(which(abs(AMIN_diff)<15))/length(AMIN_diff)
    }
    data=data.frame(ms,sim=mss$data)
    names(data)=c("date","Q","Qs")
    res_obs=data.frame(Station=rep(s,length(AMIN_obs$annualMin)),AnMin=AMIN_obs$annualMin,AnMinDate=AMIN_obs$annualMinDate)
    res_sim=data.frame(Station=rep(s,length(AMIN_sim$annualMin)),AnMin=AMIN_sim$annualMin,AnMinDate=AMIN_sim$annualMinDate)
    # for (ev in 1:length(AMAX_obs$annualMaxIndx)){
    #   #print(ev)
    #   test=StatEpisode(data, date.deb = data$date[AMAX_obs$annualMaxIndx[ev]], win=30, qx.min = 100, mode = "FL",obs=TRUE)
    #   res=rbind(res,test)
    # }
    # save=c(Station_ID=s,rmean=mean(res[,2],na.rm=T),errmean=mean(res[,3],na.rm=T))
    rsl_obs=rbind(rsl_obs,res_obs)
    rsl_sim=rbind(rsl_sim,res_sim)
  }else{print("no observations")}
  hit_outl=c(hit_outl,hit)
}
hit_outl=data.frame(hit_outl,Station_data_IDo)
hitm=match(ValidSY$V1,hit_outl$Station_data_IDo)
hit_outl=hit_outl[hitm,]
boxplot(hit_outl$hit_out)
#1- mean on Anmax for every station
cataggl_obs=aggregate(list(Q=rsl_obs$AnMin),
                     by = list(Station=rsl_obs$Station),
                     FUN = function(x) c(med=median(x,na.rm=T),mean=mean(x,na.rm=T),dev=sd(x,na.rm=T),len=length(x)))
cataggl_obs <- do.call(data.frame, cataggl_obs)

cataggl_sim=aggregate(list(Q=rsl_sim$AnMin),
                     by = list(Station=rsl_sim$Station),
                     FUN = function(x) c(med=median(x,na.rm=T),mean=mean(x,na.rm=T),dev=sd(x,na.rm=T),len=length(x)))
cataggl_sim <- do.call(data.frame, cataggl_sim)

cataggl_obssim=inner_join(cataggl_obs,cataggl_sim, by="Station",suffix = c("obs","sim"))

plot(cataggl_obssim$Q.meanobs,cataggl_obssim$Q.meansim)

#reataining only good matches
matV=match(ValidSY$V1,cataggl_obssim$Station)

Cataggl_f=cataggl_obssim[matV,]
Cataggl_f$dPeaks=(Cataggl_f$Q.meansim-Cataggl_f$Q.meanobs)/Cataggl_f$Q.meanobs*100
plot(Cataggl_f$Q.meanobs,Cataggl_f$Q.meansim)


#Spatial plot, where is dpeak shit?
palet=c(hcl.colors(9, palette = "RdBu", alpha = NULL, rev = F, fixup = TRUE))
points=inner_join(ValidSY,Cataggl_f,by=c("V1"="Station"))
parpl <- st_as_sf(points, coords = c("Var1", "Var2"), crs = 4326)
parpl <- st_transform(parpl, crs = 3035)
cord.dec=points[,c(1,2)]
cord.dec = SpatialPoints(cord.dec, proj4string=CRS("+proj=longlat"))
cord.UTM <- spTransform(cord.dec, CRS("+init=epsg:3035"))
nco=cord.UTM@coords
world <- ne_countries(scale = "medium", returnclass = "sf")
Europe <- world[which(world$continent == "Europe"),]
e2=st_transform(Europe,  crs=3035)
w2=st_transform(world,  crs=3035)
tsize=12
osize=12
limi=c(-100,100)
basemap=w2
ggplot(basemap) +
  geom_sf(fill="gray95", color=NA) +
  geom_sf(data=parpl,aes(geometry=geometry,fill=dMin,size=UpA),color="transparent",alpha=.9,shape=21,stroke=0)+ 
  geom_sf(fill=NA, color="grey20") +
  geom_sf(data=parpl,aes(geometry=geometry,size=UpA),col="black",alpha=1,stroke=0.1,shape=1)+ 
  scale_x_continuous(breaks=seq(-30,40, by=5)) +
  scale_size(range = c(1, 3), trans="sqrt")+
  scale_fill_gradientn(
    colors=palet,n.breaks=5,oob = scales::squish, limits=limi, name="mean Qpeak difference") +
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
ggsave("Plots/Spatial_QMindiff.jpg", width=20, height=15, units=c("cm"),dpi=1500)

matO=which(!is.na(match(rsl_obs$Station,ValidSY$V1)))
matS=which(!is.na(match(rsl_sim$Station,ValidSY$V1)))

rsl_obsF=rsl_obs[matO,]
rsl_simF=rsl_sim[matS,]

rsl_obsF$yday=yday(rsl_obsF$AnMinDate)
seasonMin_obs=aggregate(list(Q=rsl_obsF$yday),
                        by = list(Station=rsl_obsF$Station),
                        FUN = function(x) c(season=season1(x),len=length(x)))
seasonMin_obs <- do.call(data.frame, seasonMin_obs)

rsl_simF$yday=yday(rsl_simF$AnMinDate)
seasonMin_sim=aggregate(list(Q=rsl_simF$yday),
                        by = list(Station=rsl_simF$Station),
                        FUN = function(x) c(season=season1(x),len=length(x)))
seasonMin_sim <- do.call(data.frame, seasonMin_sim)

seasonMin_obs$dseason=seasonMin_obs$Q.season-seasonMin_sim$Q.season
seasonMin_obs$dseason[which(seasonMin_obs$dseason>182)]=seasonMin_obs$dseason[which(seasonMin_obs$dseason>182)]-365.25
seasonMin_obs$dseason[which(seasonMin_obs$dseason<=-182)]=365.25+seasonMin_obs$dseason[which(seasonMin_obs$dseason<=-182)]

hist(seasonMin_obs$dseason)
palet=c(hcl.colors(9, palette = "RdBu", alpha = NULL, rev = F, fixup = TRUE))
points=inner_join(ValidSY,seasonMin_obs,by=c("V1"="Station"))
parpl <- st_as_sf(points, coords = c("Var1", "Var2"), crs = 4326)
parpl <- st_transform(parpl, crs = 3035)
cord.dec=points[,c(1,2)]
cord.dec = SpatialPoints(cord.dec, proj4string=CRS("+proj=longlat"))
cord.UTM <- spTransform(cord.dec, CRS("+init=epsg:3035"))
nco=cord.UTM@coords
world <- ne_countries(scale = "medium", returnclass = "sf")
Europe <- world[which(world$continent == "Europe"),]
e2=st_transform(Europe,  crs=3035)
w2=st_transform(world,  crs=3035)
tsize=12
osize=12
limi=c(-100,100)
basemap=w2
ggplot(basemap) +
  geom_sf(fill="gray95", color=NA) +
  geom_sf(data=parpl,aes(geometry=geometry,fill=dseason,size=UpA),color="transparent",alpha=.9,shape=21,stroke=0)+ 
  geom_sf(fill=NA, color="grey20") +
  geom_sf(data=parpl,aes(geometry=geometry,size=UpA),col="black",alpha=1,stroke=0.1,shape=1)+ 
  scale_x_continuous(breaks=seq(-30,40, by=5)) +
  scale_size(range = c(1, 3), trans="sqrt")+
  scale_fill_gradientn(
    colors=palet,n.breaks=5,oob = scales::squish, limits=limi, name="Qpeak season difference") +
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
ggsave("Plots/Spatial_SMindiff.jpg", width=20, height=15, units=c("cm"),dpi=1500)

names(hit_out)=names(hit_outl)
Perfres=rbind(hit_out,hit_outl)
Perfres$hit_outl=Perfres$hit_outl*100
Perfres$group=1
Perfres$group[c(2902:5802)]=2
ggplot(Perfres, aes(x=factor(group), y=hit_outl)) +
  # labs(title=tit1,x="Warmimg levels (?C)", y = y1)+
  geom_boxplot(notch=F,position=position_dodge(.9),alpha=.8,aes(fill=factor(group)),linewidth=0.8,outlier.alpha = 0.4)+
  scale_y_continuous(limits = c(0,100),name="Hit rate (%)",breaks = seq(0,100,by=10),minor_breaks = seq(-100,100,5))+
  scale_x_discrete(labels=c("1" = "Annual maxs", "2" = "Annual mins"),name="")+
  scale_fill_manual(values = c("1"="royalblue","2" ="tomato"), name = "")+
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
ggsave("Plots/Bxplot_Hitrate.jpg", width=20, height=15, units=c("cm"),dpi=1500)



Perfres=rbind(seasonMax_obs[,c(1,4)],seasonMin_obs[,c(1,4)])
Perfres$group=1
Perfres$group[c(2902:5802)]=2

ggplot(Perfres, aes(x=factor(group), y=dseason)) +
  # labs(title=tit1,x="Warmimg levels (?C)", y = y1)+
  geom_violin(aes(fill=factor(group)),position=position_dodge(.9),alpha=0.8,size=0.9)+
  geom_boxplot(width=0.05,notch=F,size=0.8)+
  coord_flip()+
  #geom_boxplot(notch=F,position=position_dodge(.9),alpha=.8,aes(fill=factor(group)),linewidth=0.8,outlier.alpha = 0.4)+
  scale_y_continuous(limits = c(-60,60),name="Error (days)",breaks = seq(-60,60,by=10),minor_breaks = seq(-100,100,5))+
  scale_x_discrete(labels=c("1" = "Annual maxs", "2" = "Annual mins"),name="")+
  scale_fill_manual(values = c("1"="#9B98F6","2" ="#FF8888"), name = "")+
  theme(axis.title=element_text(size=16),
        axis.text = element_text(size=14),
        panel.background = element_rect(fill = "white", colour = "white"),
        panel.grid = element_blank(),
        panel.border = element_rect(linetype = "solid", fill = NA, colour="black"),
        legend.title = element_text(size=14),
        legend.text = element_text(size=12),
        legend.position = "none",
        panel.grid.major = element_line(colour = "grey60"),
        panel.grid.minor = element_line(colour = "grey80",linetype="dashed"),
        legend.key = element_rect(fill = "transparent", colour = "transparent"),
        legend.key.size = unit(.8, "cm")) 

ggsave("Plots/Bxplot_seasonF.jpg", width=20, height=7, units=c("cm"),dpi=1500)



Perfres=rbind(hit_out,Cataggl_f[,c(1,10)])
Perfres$group=1
Perfres$group[c(2902:5802)]=2
ggplot(Perfres, aes(x=factor(group), y=dPeaks)) +
  # labs(title=tit1,x="Warmimg levels (?C)", y = y1)+
  #geom_violin(scale="width",draw_quantiles = c(0.25, 0.5, 0.75),position=position_dodge(.9),alpha=0.7)
  geom_boxplot(notch=F,position=position_dodge(.9),alpha=.8,aes(fill=factor(group)),linewidth=0.8,outlier.alpha = 0.4)+
  scale_y_continuous(limits = c(-100,100),name="Error (%)",breaks = seq(-100,100,by=20),minor_breaks = seq(-100,100,10))+
  scale_x_discrete(labels=c("1" = "Flood peaks", "2" = "Drought lows"),name="")+
  scale_fill_manual(values = c("1"="royalblue","2" ="tomato"), name = "")+
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
ggsave("Plots/Bxplot_peakdiff.jpg", width=20, height=15, units=c("cm"),dpi=1500)

#Now I just need to redu the scatterplot for different quantiles and goodbye

obsinr=rbind(Cataggl_f[,c(1,2,3)])
siminr=rbind(Cataggl_f[,c(1,6,7)])
names(obsinr)=names(siminr)=names(obsin)
obsinr=obsinr[which(obsinr$value>0),]
siminr=siminr[which(siminr$value>0),]
nobs="AnMin_obs"
nsim="AnMin_sim"
scplot="scatterplot_AnMin"
Myplot2=ScatterValid(obsinr,nobs,siminr,nsim,scplot,valid_path,ValidSY)

Myplot2
ggsave(paste0("Plots/",scplot,".jpg"), width=20, height=15, units=c("cm"),dpi=1500)



outf=c()
for( Nsq in 1:88){
  print(Nsq)
  rspace= read.csv(paste0(hydroDir,"/subspace_efas.csv"))
  rspace=rspace[,-1]
  nrspace=rspace[Nsq,]
  outletname="efas_rnet_100km_01min"
  
  outhybas=outletopen(hydroDir,outletname,nrspace)
  Idstart=as.numeric(Nsq)*100000
  if (length(outhybas$outlets)>0){
    outhybas$outlets=seq((Idstart+1),(Idstart+length(outhybas$outlets)))
    outhybas$latlong=paste(round(outhybas$Var1,4),round(outhybas$Var2,4),sep=" ")
    # outcut=which(!is.na(match(outhybas$outlets,Paramsfl$catchment)))
    # zebi=unique(Paramsfl$catchment)
    # outhloc=outhybas[outcut,]
    outf=rbind(outf,outhybas)
  }
  
  
}


### Hybas07 ----
Catchmentrivers7=read.csv(paste0(hydroDir,"/Catchments/from_hybas_eu_onlyid.csv"),encoding = "UTF-8", header = T, stringsAsFactors = F)
hybas07 <- read_sf(dsn = paste0(hydroDir,"/Catchments/hydrosheds/hybas_eu_lev07_v1c.shp"))
hybasf7=fortify(hybas07) 
Catamere07=inner_join(hybasf7,Catchmentrivers7,by= "HYBAS_ID")
Catamere07$llcoord=paste(round(Catamere07$POINT_X,4),round(Catamere07$POINT_Y,4),sep=" ") 
#maybe useless
Catf7=inner_join(Catamere07,outf,by= c("llcoord"="latlong"))
cst7=st_transform(Catf7,  crs=3035)


data2=inner_join(ValidStar,cst7,by=c("Var1","Var2"))
# points <- st_as_sf(data2, coords = c("Var1", "Var2"), crs = 4326)
# points <- st_transform(points, crs = 3035)

#Aggregate by catchment
pointagg=aggregate(list(KGE=data2$kge.mean),
                   by = list(HYBAS_ID=data2$HYBAS_ID),
                   FUN = function(x) c(med=median(x,na.rm=T),mean=mean(x,na.rm=T),dev=sd(x,na.rm=T),len=length(x)))
pointagg <- do.call(data.frame, pointagg)


colorz=c("red","orange","hotpink","darkorchid1","royalblue","darkblue","black")


pointagg$kge=pointagg$KGE.med
pointagg$kge[which(is.infinite(pointagg$kge))]=NA
pointagg$kgecode=0
pointagg$kgecode[which(pointagg$kge>-0.41 & pointagg$kge<=0.2)]=1
pointagg$kgecode[which(pointagg$kge>0.2 & pointagg$kge<=0.5)]=2
pointagg$kgecode[which(pointagg$kge>0.5 & pointagg$kge<=0.7)]=3
pointagg$kgecode[which(pointagg$kge>0.7 & pointagg$kge<=0.8)]=4
pointagg$kgecode[which(pointagg$kge>0.8 & pointagg$kge<=0.9)]=5
pointagg$kgecode[which(pointagg$kge>0.9)]=6

pointplot=inner_join(hybasf7,pointagg,by= "HYBAS_ID")

tsize=16
osize=16
pl2=ggplot(basemap) +
  geom_sf(fill="gray95")+
  geom_sf(data=point,aes(fill=factor(kgecode),geometry=geometry),color="transparent")+ 
  geom_sf(fill=NA, color="grey") +
  geom_sf(data=parpl,aes(geometry=geometry),color="grey60",alpha=.7,size=1.4,stroke=0,shape=16)+
  geom_sf(data=parpef,aes(geometry=geometry),col="black",alpha=.7,size=1.4,stroke=0.5,shape=1)+
  coord_sf(xlim = c(min(nco[,1]),max(nco[,1])), ylim = c(min(nco[,2]),max(nco[,2])))+
  # coord_sf(xlim = c(min(nco[,1]),max(nco[,1])), ylim = c(min(nco[,2]),max(nco[,2])))+
  # geom_point(data=event,aes(x=lon,y=lat,col=dr))+
  # scale_color_distiller(palette = "Spectral",
  #                       direction = -1, limit=c(0,max(eventp$dr)),oob = scales::squish)+
  scale_fill_manual(values = colorz, labels= c("< -0.41","-0.41 - 0.2", "0,2 - 0.5","0.5 - 0.7","0.7 - 0.8","0.8 - 0.9", ">0.9"), name="KGE")   +
  # scale_fill_gradientn(
  #   colors=palet,
#     breaks=br,limits=limi,trans=trans,
#     oob = scales::squish,na.value=colNA, name=legend)   +
  labs(x="Longitude", y = "Latitude")+
  guides(colour = guide_legend(override.aes = list(size = 4)))+
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
pl2
ggsave("Validation_KGEmedcatch.jpg", width=20, height=20, units=c("cm"),dpi=1500)
#Histogram of all considered station upstream area




matp100=matpix[which(matpix$DrainingArea.km2.LDD>100),]
df=matp100
p<-ggplot(df, aes(x=DrainingArea.km2.LDD)) + 
  geom_histogram(color="steelblue", fill="slategray1",bins=15,alpha=0.9,lwd=1)+
  scale_y_continuous(name="Number of stations")+
  scale_x_log10(name="Upstream area (km2)",
                breaks=c(100,1000,10000,100000), minor_breaks = log10_minor_break(),
                labels=c("100","1 000","10 000","100 000")) +
  theme(axis.title=element_text(size=16, face="bold"),
        axis.text = element_text(size=16),
        panel.background = element_rect(fill = "white", colour = "white"),
        panel.grid = element_blank(),
        panel.border = element_rect(linetype = "solid", fill = NA, colour="black"),
        legend.title = element_text(size=14),
        legend.text = element_text(size=12),
        panel.grid.major = element_line(colour = "grey80"),
        panel.grid.minor.x = element_line(colour = "grey90",linetype="dashed"),
        legend.key = element_rect(fill = "transparent", colour = "transparent"),
        legend.key.size = unit(.8, "cm"))+
  annotate("label", x=250000, y=140, label= paste0("n = ",length(df$upAc)),size=6)

p
ggsave("histo_stations.jpg", p, width=20, height=15, units=c("cm"),dpi=1500)





#################OLD shitty station detection########################
shittystations=dat_stat[which(dat_stat$rat1>5 | dat_stat$rat2>5),]

#Remove all strange 
shitty_ori=inner_join(shittystations,dat,by=c("station"="Station_ID"))



shitty_rmv=shitty_ori[which(shitty_ori$qm>10),]
shitty_rmv$uniki=paste0(shitty_rmv$station,shitty_rmv$year)
dat$uniki=paste0(dat$Station_ID,dat$year)
rmv=(match(shitty_rmv$uniki,dat$uniki))
do=dat[rmv,]
dat=dat[-rmv,]



#load stations from EFAS to check if there are faulty matches

efas_stations=read.csv(paste0(valid_path,"hres_calib_stations.csv"))
matefas=match(flefas$StationID,shittystations$station)
#First I join with my matching stations to have all fieds
ptn=inner_join(flefas,efas_stations,by=c("National_Station_Identifier"))

shittystations_dat=inner_join(ptn,shittystations, by=c("StationID"="station"))
shittystations_keep=shittystations_dat[which(!is.na(shittystations_dat$Cal_area)),]

ggplot(shittystations_dat) + geom_point(aes(x=obs, y=sim,col=factor(Station_ID)),alpha=0.3) +
  scale_x_log10(name="observed",
                breaks=c(1,10,100,1000), minor_breaks = log10_minor_break(),
                labels=c("1","10","100","1000"), limits=c(0.001,tmax)) +
  scale_y_log10(name="simulated",
                breaks=c(1,10,100,1000), minor_breaks = log10_minor_break(),
                labels=c("1","10","100","1000"), limits=c(0.001,tmax)) +
  geom_abline(slope=1, intercept=0, lwd =1, alpha=0.5)+
  # scale_x_continuous(name=nobs,breaks=seq(0,tmax, by=bi))+
  # scale_y_continuous(name=nsim,breaks=seq(0,tmax, by=bi))+
  annotate("label", x=800, y=max(dat$sim), label= paste0("R2 = ",r2),size=5)+
  scale_color_discrete()+
  theme(axis.title=element_text(size=16, face="bold"),
        axis.text = element_text(size=16),
        panel.background = element_rect(fill = "white", colour = "white"),
        panel.grid = element_blank(),
        panel.border = element_rect(linetype = "solid", fill = NA, colour="black"),
        legend.title = element_text(size=14),
        legend.text = element_text(size=12),
        legend.position = "none",
        panel.grid.major = element_line(colour = "grey80"),
        panel.grid.minor = element_line(colour = "grey90",linetype="dashed"),
        legend.key = element_rect(fill = "transparent", colour = "transparent"),
        legend.key.size = unit(.8, "cm"))


#Remove all strange 
shitty_ori=inner_join(shittystations,dat,by=c("station"="Station_ID"))
shitty_rmv=shitty_ori[which(shitty_ori$qm>10),]
shitty_rmv$uniki=paste0(shitty_rmv$station,shitty_rmv$year)
dat$uniki=paste0(dat$Station_ID,dat$year)
rmv=(match(shitty_rmv$uniki,dat$uniki))
do=dat[rmv,]
dat=dat[-rmv,]




