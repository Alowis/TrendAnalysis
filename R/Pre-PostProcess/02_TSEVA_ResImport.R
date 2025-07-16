


setwd(dirname(rstudioapi::getSourceEditorContext()$path))
getwd()
source("functions_trends2.R")
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
scenario="Histo"
hazard="Flood"
season="year"
mmx=""
lf=list.files(path = paste0(hydroDir,"/TSEVA/",scenario,"/",hazard,"/",season,mmx), full.names = TRUE, recursive = TRUE)
lf
nfiles=length(lf)
cal=T
if (nfiles<175){
  print("files missing")
}
for (file in lf){
  #file=lf[1]
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
  parlist=Results$parameters
  parlist=data.table(parlist)
  rspace= read.csv(paste0(hydroDir,"/subspace_efas.csv"))
  rspace=rspace[,-1]
  nrspace=rspace[Nsq,]
  
  #specify outlets
  outletname="efas_rnet_100km_01min"
  
  outhybas=outletopen(hydroDir,outletname,nrspace)
  Idstart=as.numeric(Nsq)*10000
  Idstart2=as.numeric(Nsq)*100000
  if (length(outhybas$outlets)>0){
    outhybas$outlets=seq((Idstart+1),(Idstart+length(outhybas$outlets)))
    outhybas$outl2=seq((Idstart2+1),(Idstart2+length(outhybas$outlets)))
    outhybas$latlong=paste(round(outhybas$Var1,4),round(outhybas$Var2,4),sep=" ")
    out2cat=seq(parlist$catchment[1],parlist$catchment[length(parlist$catchment)])
    outcut=which(!is.na(match(outhybas$outl2,out2cat)))
    outhloc=outhybas[outcut,]
    outf=rbind(outf,outhloc)
  }
  unikout=unique(outhloc$outl2)
  
  pc=c()
  
  RetLevGPD=Results$RetLevGPD
  RetPerGPD=Results$RetPerGPD
  RetLevGEV=Results$RetLevGEV
  RetPerGEV=Results$RetPerGEV
  
  #remove values for year 1950 is necessary
  rm50=which(colnames(RetLevGPD)=="1950")
  
  if (length(rm50)>0){
    print("remove")
    RetLevGPD=RetLevGPD[,-rm50]
    RetPerGPD=RetPerGPD[,-rm50]
    RetLevGEV=RetLevGEV[,-rm50]
    RetPerGEV=RetPerGEV[,-rm50]
  }
  
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

print(paste0(hydroDir,"/",hazard,"/RL100.",hazard,".",season,".",scenario,mmx,".Rdata"))
# Saving outputs of the loop
saveout=T
if (saveout==T){
  
  save(RLGPDfl,file=paste0(hydroDir,"/",hazard,"/RL100.",hazard,".",season,".",scenario,mmx,"2.Rdata"))
  save(Paramsfl,file=paste0(hydroDir,"/",hazard,"/params.",hazard,".",season,".",scenario,mmx,"2.Rdata"))
  save(Peaksave,file=paste0(hydroDir,"/",hazard,"/peaks.",hazard,".",season,".",scenario,mmx,"2.Rdata"))
}

gc()


## Intermittent Rivers and Ephemeral Streams ----

IRES_save=c()
hydroDir<-("D:/tilloal/Documents/LFRuns_utils/data")
outlets="Rnet"
scenario="SocCF"
hazard="Drought"
season="nonfrost"
mmx=""
lf=list.files(path = paste0(hydroDir,"/TSEVA/",scenario,"/",hazard,"/",season,mmx), full.names = TRUE, recursive = TRUE)
lf
nfiles=length(lf)
cal=T
if (nfiles<175){
  print("files missing")
}
for (file in lf){
  #file=lf[1]
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
  parlist=Results$parameters
  parlist=data.table(parlist)
  
  # unikpar=unique(parlist$catchment)
  # dparmerde=diff(unikpar)
  # plot(dparmerde)
  
  rspace= read.csv(paste0(hydroDir,"/subspace_efas.csv"))
  rspace=rspace[,-1]
  nrspace=rspace[Nsq,]
  
  
  
  IRES=Results$catrest
 
  IRES_save=rbind(IRES_save,IRES)
  
}

length(which(IRES_save$IRES>1))

print(paste0(hydroDir,"/",hazard,"/IRES.",season,".",scenario,mmx,"_new2.Rdata"))

save(IRES_save,file=paste0(hydroDir,"/",hazard,"/IRES.",season,".",scenario,mmx,"_new2.Rdata"))


