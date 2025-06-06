
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
library(matrixStats)


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



# Example usage of TsEvaNs function
timeAndSeries <- ArdecheStMartin
#go from six-hourly values to daily max
timeAndSeries <- max_daily_value(timeAndSeries)
#keep only the 30 last years
yrs <- as.integer(format(timeAndSeries$date, "%Y"))
tokeep <- which(yrs>=1990)
timeAndSeries <- timeAndSeries[tokeep,]
timeWindow <- 10*365 # 10 years
result <- TsEvaNs(timeAndSeries, timeWindow,
                  transfType = 'trendPeaks',tail = 'high')
mierda=result$nonStationaryEvaParams$potObj

#2 Pre-loaded results -----------
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
basemap=w2



##2.2 Loading saved results in .Rdata ---------------------------

###load UpArea -----
#load upstream area
main_path = 'D:/tilloal/Documents/06_Floodrivers/'
valid_path = paste0(main_path,'DataPaper/')
outletname="/GIS/upArea_European_01min.nc"
dir=valid_path
outf$idlalo=paste(outf$idlo, outf$idla, sep=" ")
UpArea=UpAopen(valid_path,outletname,outf)
head(UpArea)

###load historical run -----

haz="Drought"

if (haz == "Drought") namefile="Drought.nonfrost.Histo22"
if (haz == "Flood") namefile="flood.year.Histo22"

load(file=paste0(hydroDir,"/",haz,"/params.",namefile,".Rdata"))
load(file=paste0(hydroDir,"/",haz,"/RL100.",namefile,".Rdata"))
load(file=paste0(hydroDir,"/",haz,"/peaks.",namefile,".Rdata"))
gc()

Paramsfl=(Paramsfl[,-c(4:9,17)])
ParamsflH=Paramsfl
PeakH=Peaksave
RLGPDflH=RLGPDfl
rm(Paramsfl,RLGPDfl)
gc()

###load Socio-CF run -----
if (haz == "Drought") namefile="Drought.nonfrost.SocCF22"
if (haz == "Flood") namefile="Flood.year.socCF22"

load(file=paste0(hydroDir,"/",haz,"/params.",namefile,".Rdata"))
load(file=paste0(hydroDir,"/",haz,"/RL100.",namefile,".Rdata"))
load(file=paste0(hydroDir,"/",haz,"/peaks.",namefile,".Rdata"))

RLGPDflSCF=RLGPDfl
Paramsfl=Paramsfl[,-c(4:9,17)]
ParamsflSCF=data.table(Paramsfl)
PeakSCF=Peaksave
rm(Paramsfl,RLGPDfl)
gc()

###load results from Res+WU CF run -----
if (haz == "Drought") namefile="Drought.nonfrost.RWCF22"
if (haz == "Flood") namefile="flood.year.RWCF22"
load(file=paste0(hydroDir,"/",haz,"/params.",namefile,".Rdata"))
load(file=paste0(hydroDir,"/",haz,"/RL100.",namefile,".Rdata"))
load(file=paste0(hydroDir,"/",haz,"/peaks.",namefile,".Rdata"))
RLGPDflRWCF=RLGPDfl
Paramsfl=Paramsfl[,-c(4:9,17)]
ParamsflRWCF=data.table(Paramsfl)
PeakRWCF=Peaksave
rm(Paramsfl,RLGPDfl)
gc()

###load results from Water CF run -----
if (haz == "Drought") namefile="Drought.nonfrost.WCF2"
if (haz == "Flood") namefile="flood.year.WCF2"
load(file=paste0(hydroDir,"/",haz,"/params.",namefile,".Rdata"))
load(file=paste0(hydroDir,"/",haz,"/RL100.",namefile,".Rdata"))
load(file=paste0(hydroDir,"/",haz,"/peaks.",namefile,".Rdata"))

RLGPDflWCF=RLGPDfl
Paramsfl=Paramsfl[,-c(4:9,17)]
ParamsflWCF=data.table(Paramsfl)
PeakWCF=Peaksave
rm(Paramsfl,RLGPDfl)
gc()

rm(catmap)
gc()

ParamsflRWCF[which(ParamsflRWCF$catchment==4103819),]
RLGPDflWCF[which(RLGPDflWCF$unikout==4103819),]

hist(ParamsflSCF$epsilonGPD,xlim=c(-2,1),breaks=10000)
quantile(ParamsflSCF$epsilonGPD,0.025,na.rm=T)






###Confidence interval for change, needs to me moved on the HPC -------------------
library(RtsEva)
#test on CI for change extimate
tabsave=c()
for (ix in 1:100){
  pix=RLGPDflH$unikout[ix]
  paramD=ParamsflSCF[which(ParamsflSCF$catchment==pix),]
  idb=which(paramD$Year==2015)
  epsilon <- paramD$epsilonGPD[idb]
  sigma <- paramD$sigmaGPD[c(4,idb)]
  threshold <- paramD$thresholdGPD[c(4,idb)]
  timeHorizonInYears <- 70
  
  nPeaks <- paramD$nPeaks[1]
  epsilonStdErr <- paramD$epsilonStdErrGPD[c(4,idb)]
  sigmaStdErr <- paramD$sigmaStdErrGPD[c(4,idb)]
  thresholdStdErr <- paramD$thresholdStdErrGPD[c(4,idb)]
  
  returnPeriods <- c(10)
  
  returnLevels1 <- tsEvaComputeReturnLevelsGPD(epsilon, sigma, 
                                               threshold, epsilonStdErr, sigmaStdErr, thresholdStdErr, 
                                               nPeaks = nPeaks, sampleTimeHorizon = timeHorizonInYears, 
                                               returnPeriods = returnPeriods)
  
  RL1=-returnLevels1$returnLevels[1]
  supRLCI1 <- -(returnLevels1$returnLevels[1] - returnLevels1$returnLevelsErr[1])
  infRLCI1 <- -(returnLevels1$returnLevels[1] + returnLevels1$returnLevelsErr[1])
  Rlba=data.frame(ts="base",rl=RL1,sci=supRLCI1,ici=infRLCI1)
  
  RL2=-returnLevels1$returnLevels[2]
  supRLCI2 <- -(returnLevels1$returnLevels[2] - returnLevels1$returnLevelsErr[2])
  infRLCI2 <- -(returnLevels1$returnLevels[2] + returnLevels1$returnLevelsErr[2])
  Rlpre=data.frame(ts="pre",rl=RL2,sci=supRLCI2,ici=infRLCI2)
  
  if(!is.nan(supRLCI1)){
    if(RL2<supRLCI1) {
      col="green"
      sig=1
    }
    if(RL2>supRLCI1) {
      col="red"
      sig=0.5*0.16
    }
    if(RL2<infRLCI1) {
      col="red"
      sig=0.5*0.16
    }
    if(supRLCI2<infRLCI1) {
      col="purple"
      sig=0.16*0.16
    }
    if(infRLCI2>supRLCI1) {
      col="purple"
      sig=0.16*0.16
    }
  }else(
    sig=NA
  )
  

  
  tavbl=rbind(Rlba,Rlpre)
  tavbl$id=pix
  tavbl$sig=sig
  tavbl$col=col
  colx=c("black",col)
  tabsave=rbind(tabsave,tavbl)

}

#plot any pixel
tabvl=tabsave[which(tabsave$id==pix),]
# Create a plot with rl as a point and a line from ici to sci
colx=c("black",tabvl$col[1])
ggplot(tavbl, aes(x = ts)) + 
  geom_segment(aes(x = ts, xend = ts, y = ici, yend = sci), size = 1,alpha=0.6) + 
  geom_point(aes(y = rl),color=colx,size=4) + 
  labs(title = "Comparison of Base and Pre", x = "Category", y = "Value")





# Define the function that performs the calculations for a single row
pix=RLGPDflH$unikout[2]
calculate_mcci <- function(pix, ParamsflSCF, timeHorizonInYears = 70, returnPeriods = c(10)) {
  paramD <- ParamsflSCF[which(ParamsflSCF$catchment == pix), ]
  idb <- which(paramD$Year==2015)
  
  epsilon <- paramD$epsilonGPD[idb]
  sigma <- paramD$sigmaGPD[c(4, idb)]
  threshold <- paramD$thresholdGPD[c(4, idb)]
  
  nPeaks <- paramD$nPeaks[1]
  epsilonStdErr <- paramD$epsilonStdErrGPD[c(4, idb)]
  sigmaStdErr <- paramD$sigmaStdErrGPD[c(4, idb)]
  thresholdStdErr <- paramD$thresholdStdErrGPD[c(4, idb)]
  
  returnLevels1 <- tsEvaComputeReturnLevelsGPD(
    epsilon, sigma, threshold, epsilonStdErr, sigmaStdErr, thresholdStdErr,
    nPeaks = nPeaks, sampleTimeHorizon = timeHorizonInYears, returnPeriods = returnPeriods
  )
  
  supRLCI1 <- returnLevels1$returnLevels[1] + returnLevels1$returnLevelsErr[1]
  infRLCI1 <- returnLevels1$returnLevels[1] - returnLevels1$returnLevelsErr[1]
  
  supRLCI2 <- returnLevels1$returnLevels[2] + returnLevels1$returnLevelsErr[2]
  infRLCI2 <- returnLevels1$returnLevels[2] - returnLevels1$returnLevelsErr[2]
  
  #((returnLevels1$returnLevels[2] - returnLevels1$returnLevels[1]) / returnLevels1$returnLevels[1] * 100)
  maxCI <- ((supRLCI2 - supRLCI1) / returnLevels1$returnLevels[1] * 100)
  minCI <- ((infRLCI2 - infRLCI1) / returnLevels1$returnLevels[1] * 100)
  
  myCI <- abs(maxCI - minCI)
  return(myCI)
}

# Use sapply to apply the function to each element in RLGPDflH$unikout
library(RtsEva)


time_taken <- system.time({
  mcci <- sapply(RLGPDflH$unikout[1:10], calculate_mcci, ParamsflSCF = ParamsflSCF)
})
# Print the time taken
print(time_taken)
hist(mcci, breaks=10)
median(mcci,na.rm=T)


# Load necessary libraries
library(parallel)
library(RtsEva)

#recuperate the wight I used to assess changes for low flows for consistence with other analysis

# Define the function that performs the calculations for a single row
calculate_mcci <- function(paramD, timeHorizonInYears = 70, returnPeriods = c(10)) {
  idb <- which(paramD$Year==2015)
  epsilon <- paramD$epsilonGPD[1]
  sigma <- paramD$sigmaGPD[c(4, idb)]
  threshold <- paramD$thresholdGPD[c(4, idb)]
  
  nPeaks <- paramD$nPeaks[1]
  epsilonStdErr <- paramD$epsilonStdErrGPD[c(4, idb)]
  sigmaStdErr <- paramD$sigmaStdErrGPD[c(4, idb)]
  thresholdStdErr <- paramD$thresholdStdErrGPD[c(4, idb)]
  
  returnLevels1 <- tsEvaComputeReturnLevelsGPD(
    epsilon, sigma, threshold, epsilonStdErr, sigmaStdErr, thresholdStdErr,
    nPeaks = nPeaks, sampleTimeHorizon = timeHorizonInYears, returnPeriods = returnPeriods
  )
  
  
  supRLCI1 <- returnLevels1$returnLevels[1] + returnLevels1$returnLevelsErr[1]
  infRLCI1 <- returnLevels1$returnLevels[1] - returnLevels1$returnLevelsErr[1]
  
  #CI1=(returnLevels1$returnLevelsErr[2]-returnLevels1$returnLevelsErr[1])/returnLevels1$returnLevels[1] * 100
  supRLCI2 <- returnLevels1$returnLevels[2] + returnLevels1$returnLevelsErr[2]
  infRLCI2 <- returnLevels1$returnLevels[2] - returnLevels1$returnLevelsErr[2]
  meanCI<- (returnLevels1$returnLevels[2] - returnLevels1$returnLevels[1]) 
  maxCI <- (supRLCI2 - supRLCI1) 
  minCI <- (infRLCI2 - infRLCI1) 

  initRL=returnLevels1$returnLevels[1]
  myCI <- c(maxCI,meanCI,minCI,initRL)
  return(myCI)

}

pix=RLGPDflH$unikout[100]
merde=calculate_mcci(paramD=ParamsflSCF[which(ParamsflSCF$catchment == pix),])

calculate_tbs <- function(paramD, timeHorizonInYears = 70, returnPeriods = c(10)) {
  idb <- which(paramD$Year==2015)
  epsilon <- paramD$epsilonGPD[1]
  sigma <- paramD$sigmaGPD[c(4, idb)]
  threshold <- paramD$thresholdGPD[c(4, idb)]
  
  nPeaks <- paramD$nPeaks[1]
  epsilonStdErr <- paramD$epsilonStdErrGPD[c(4, idb)]
  sigmaStdErr <- paramD$sigmaStdErrGPD[c(4, idb)]
  thresholdStdErr <- paramD$thresholdStdErrGPD[c(4, idb)]
  
  returnLevels1 <- tsEvaComputeReturnLevelsGPD(
    epsilon, sigma, threshold, epsilonStdErr, sigmaStdErr, thresholdStdErr,
    nPeaks = nPeaks, sampleTimeHorizon = timeHorizonInYears, returnPeriods = returnPeriods
  )
  
  RL1=-returnLevels1$returnLevels[1]
  supRLCI1 <- -(returnLevels1$returnLevels[1] - returnLevels1$returnLevelsErr[1])
  infRLCI1 <- -(returnLevels1$returnLevels[1] + returnLevels1$returnLevelsErr[1])
  
  RL2=-returnLevels1$returnLevels[2]
  supRLCI2 <- -(returnLevels1$returnLevels[2] - returnLevels1$returnLevelsErr[2])
  infRLCI2 <- -(returnLevels1$returnLevels[2] + returnLevels1$returnLevelsErr[2])
  
  if(!is.nan(supRLCI1)){
    if(RL2<supRLCI1) {
      col=0
      sig=1
    }
    if(RL2>supRLCI1) {
      col=1
      sig=0.5*0.16
    }
    if(RL2<infRLCI1) {
      col=1
      sig=0.5*0.16
    }
    if(supRLCI2<infRLCI1) {
      col=2
      sig=0.16*0.16
    }
    if(infRLCI2>supRLCI1) {
      col=2
      sig=0.16*0.16
    }
  }else{
    sig=NA
    supRLCI1=NA
    supRLCI2=NA
    infRLCI1=NA
    infRLCI2=NA
    col=NA
  }
  Rlpre=c(RL2,supRLCI2,infRLCI2)
  #Rlpre=data.frame(ts="pre",rl=RL2,sci=supRLCI2,ici=infRLCI2)
  
  Rlba=c(RL1,supRLCI1,infRLCI1)
  #Rlba=data.frame(ts="base",rl=RL1,sci=supRLCI1,ici=infRLCI1)
 
  tavbl=c(Rlba,Rlpre,paramD$catchment[1],sig,col) 
  # tavbl=rbind(Rlba,Rlpre)
  # tavbl$id=paramD$catchment[1]
  # tavbl$sig=sig
  # tavbl$col=col
  return(tavbl)
  
}
merde=calculate_tbs(paramD=ParamsflSCF[which(ParamsflSCF$catchment == pix),])
# Use parallel processing to apply the function
pentry=ParamsflSCF[1:1000]
time_taken <- system.time({
  mcci <- sapply(RLGPDflH$unikout[1:100], calculate_tbs, paramD = ParamsflSCF)
})

cl <- makeCluster(detectCores() - 3)  # Use all but one core
clusterExport(cl, "tsEvaComputeReturnLevelsGPD")  # Export the function to the cluster
clusterExport(cl, "calculate_tbs")  # Export the function to the cluster
clusterExport(cl, "ParamsflSCF") 
time_taken <- system.time({
  tbsave <- parSapply(cl, RLGPDflH$unikout[1:1000], function(pix) {
    paramD <- ParamsflSCF[which(ParamsflSCF$catchment == pix), ]
    if (!is.null(paramD)) {
      calculate_tbs(paramD)
    } else {
      NA  # Handle cases where pix is not found in param_list_named
    }
  })
})

# Stop the cluster
stopCluster(cl)

# Print the time taken
print(time_taken)

tbsave2=data.frame(t(tbsave))

# Combine the results
mcci_df <- do.call(rbind, mcci)
rownames(mcci_df) <- RLGPDflH$unikout[1:100]

plot(mcci_df[,2],type="o")
lines(mcci_df[,1],col=2)
lines(mcci_df[,3],col=3)
plot(abs(mcci_df[,1]-mcci_df[,2])/abs(mcci_df[,2])*100,col="red", lwd=2,type="o",log="y")

relative_uncertainty=abs(mcci_df[,1]-mcci_df[,2])/abs(mcci_df[,2])*100
median(relative_uncertainty,na.rm=T)
hist(relative_uncertainty, breaks=3000, xlim=c(0,1000))
lines((mcci_df[,3]-mcci_df[,2]),col="blue", lwd=2)
mccx=unlist(mccx)
names(mccx)=NA

hist(mcci,breaks=100)

#To be brought on the HPC for super fast runs


function (epsilon, sigma, threshold, epsilonStdErr, sigmaStdErr, 
          thresholdStdErr, nPeaks, sampleTimeHorizon, returnPeriods) 
{
  X0 <- nPeaks/sampleTimeHorizon
  XX <- X0 * returnPeriods
  npars <- length(sigma)
  nt <- length(returnPeriods)
  XX_ <- rep(XX, each = npars)
  sigma_ <- rep(sigma, times = nt)
  sigmaStdErr_ <- rep(sigmaStdErr, times = nt)
  threshold_ <- rep(threshold, times = nt)
  thresholdStdErr_ <- rep(thresholdStdErr, times = nt)
  if (epsilon != 0) {
    returnLevels <- threshold_ + sigma_/epsilon * ((XX_)^epsilon - 
                                                     1)
    dxm_u <- 1
    dxm_sigma <- 1/epsilon * (XX_^epsilon - 1)
    dxm_epsilon <- -sigma_/epsilon^2 * ((XX_)^epsilon - 1) + 
      sigma_/epsilon * log(XX_) * XX_^epsilon
    returnLevelsErr <- sqrt((dxm_u * thresholdStdErr_)^2 + 
                              (dxm_sigma * sigmaStdErr_)^2 + (dxm_epsilon * epsilonStdErr)^2)
  }
  else {
    returnLevels <- threshold_ + sigma_ * log(XX_)
    dxm_u <- 1
    dxm_sigma <- log(XX_)
    returnLevelsErr <- ((dxm_u * thresholdStdErr_)^2 + (dxm_sigma * 
                                                          sigmaStdErr_)^2)^0.5
  }
  return(list(returnLevels = returnLevels, returnLevelsErr = returnLevelsErr))
}
