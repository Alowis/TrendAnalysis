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


# Load inputs from HPC computation ----------------------------------------

RPGPDfl=c()
RLGPDfl=c()
naallfl=c()
Paramsfl=c()
parlist=c()
Peaks=c()
Peaksave=c()
outf=c()
hydroDir<-("D:/tilloal/Documents/LFRuns_utils/data")
outlets="Rnet"
hazard="Drought"
season="nonfrost"

lf=list.files(path = paste0(hydroDir,"/TSEVA/Histo/",hazard,"/",season), full.names = TRUE, recursive = TRUE)
lf
cal=T
for (file in lf){
  load(file)
  fils=sub(".*/", "", file)
  tt=unlist(strsplit(fils, "[_]"))
  Nsq=as.numeric(tt[3])
  if (cal==F){
    if (Nsq>88) Nsq=floor(Nsq/10)
  }else if (cal==T){
    Nsq=floor(Nsq/10)
  }
  
  out.type=sub("\\_.*", "", fils)
  print(Nsq)
  if (out.type=="ResCat6h"){
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
    
    if (cal==F){
      RetLevGPD=Results$RetLevGPD[which(!is.na(Results$RetLevGPD[,1])),]
      RetPerGPD=Results$RetPerGPD[which(!is.na(Results$RetLevGPD[,1])),]
      RetLevGEV=Results$RetLevGEV[which(!is.na(Results$RetLevGPD[,1])),]
      RetPerGEV=Results$RetPerGEV[which(!is.na(Results$RetLevGPD[,1])),]
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
  }
}


# Saving outputs of the loop
saveout=T
if (saveout==T){
  save(RLGPDfl,file=paste0(hydroDir,"/Flood/RL100.flood.RWCF4.Rdata"))
  save(Paramsfl,file=paste0(hydroDir,"/Flood/params.flood.RWCF4.Rdata"))
  save(Peaksave,file=paste0(hydroDir,"/Flood/peaks.flood.RWCF4.Rdata"))
}

##do the same for socCF and start the comparison

## CALIBRATED DROUGHT ----
RPGPDdr=c()
RLGPDdr=c()
naalldr=c()
Paramsdr=c()
Peaksave=c()
outf=c()
hydroDir<-("D:/tilloal/Documents/LFRuns_utils/data")
outlets="Rnet"
hazard="Drought"
season="year"
lf=list.files(path = paste0(hydroDir,"/TSEVA/SocCF/",hazard,"/",season), full.names = TRUE, recursive = TRUE)
lf
cal=T
for (file in lf){
  file=lf[46]
  load(file)
  fils=sub(".*/", "", file)
  tt=unlist(strsplit(fils, "[_]"))
  Nsq=as.numeric(tt[3])
  if (cal==F){
    if (Nsq>88) Nsq=floor(Nsq/10)
  }else if (cal==T){
    Nsq=floor(Nsq/10)
  }
  
  out.type=sub("\\_.*", "", fils)
  print(Nsq)
  if (out.type=="ResCat6h"){
    parlist=Results$parameters
    parlist=as.data.frame(parlist)
    Results$catrest
    
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
    if (length(outhybas$outlets)>0){
      outhybas$outlets=seq((Idstart+1),(Idstart+length(outhybas$outlets)))
      outhybas$outl2=seq((Idstart2+1),(Idstart2+length(outhybas$outlets)))
      outhybas$latlong=paste(round(outhybas$Var1,4),round(outhybas$Var2,4),sep=" ")
      #outcut=which(!is.na(match(outhybas$outlets,parlist$catchment)))
      zebi=seq(parlist$catchment[1],parlist$catchment[length(parlist$catchment)])
      outcut=which(!is.na(match(outhybas$outlets,zebi)))
      outhloc=outhybas[outcut,]
      outf=rbind(outf,outhloc)
    }
    unikout=unique(outhloc$outl2)
    
    pc=c()
    
    RetLevGPD=Results$RetLevGPD
    RetPerGPD=Results$RetPerGPD
    RetLevGEV=Results$RetLevGEV
    RetPerGEV=Results$RetPerGEV
    
    if (cal==F){
      RetLevGPD=Results$RetLevGPD[which(!is.na(Results$RetLevGPD[,1])),]
      RetPerGPD=Results$RetPerGPD[which(!is.na(Results$RetLevGPD[,1])),]
      RetLevGEV=Results$RetLevGEV[which(!is.na(Results$RetLevGPD[,1])),]
      RetPerGEV=Results$RetPerGEV[which(!is.na(Results$RetLevGPD[,1])),]
    }
    
    Impdates=names(RetLevGPD)
    
    names(RetLevGPD)=paste0("Y",Impdates)
    names(RetPerGPD)=paste0("Y",Impdates)
    RetLevGPD=data.frame(RetLevGPD,unikout)
    RetPerGPD=data.frame(RetPerGPD,unikout)
    
    names(RetLevGEV)=paste0("Y",Impdates)
    names(RetPerGEV)=paste0("Y",Impdates)
    RetLevGEV=data.frame(RetLevGEV,unikout)
    RetPerGEV=data.frame(RetPerGEV,unikout)
    
    RPGPDdr=rbind(RPGPDdr,RetPerGPD)
    RLGPDdr=rbind(RLGPDdr,RetLevGPD)
    Paramsdr=rbind(Paramsdr,parlist)
    Peakfoir=Results$Peaks
    myIDs=paste(Peakfoir$timeID,Peakfoir$catch,sep=" ")
    unid=unique(myIDs)
    locu=match(unid,myIDs)
    Peaks=Peakfoir[locu,]
    Peaks$time=as.POSIXct(Peaks$time,origin="1970-01-01 00:00:00")
    Peaks$d=yday(Peaks$time)
    Peaks$year=year(Peaks$time)
    Peaksave=rbind(Peaksave,Peaks)
  }
}

saveout=F
if (saveout==T){
  save(RLGPDdr,file=paste0(hydroDir,"/Drought/RL100.drought.Soccf1.Rdata"))
  save(Paramsdr,file=paste0(hydroDir,"/Drought/params.drought.Soccf1.Rdata"))
  save(Peaksave,file=paste0(hydroDir,"/Drought/peaks.drought.Socf1.Rdata"))
}

