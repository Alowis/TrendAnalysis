#This script is dedicated to the analysis of the timeseries outputs of LISFLOOD
#Each output corresponds to a watershed In Europe
#The aim os to extract seasonnalities of snwmelt, snow cover, precipitation and soil moisture
setwd("D:/tilloal/Documents/LFRuns_utils/TS-EVA")


source("functions.R")
setwd("D:/tilloal/Documents/LFRuns_utils")
suppressWarnings(suppressMessages(library(ncdf4)))
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)
library(rgeos)
library(reshape)
library(raster)
library(dplyr)
library(ggplot2)
library(viridis)
library(hrbrthemes)
library(tidyverse)
library(raster)

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

domainopen=function(dir,outletname,nrspace=rep(NA,5)){
  ncbassin=paste0(dir,"/",outletname,".nc")
  ncb=nc_open(ncbassin)
  name.vb=names(ncb[['var']])
  namev=name.vb[1]
  if ("Band1"%in% name.vb)namev="Band1"
  #time <- ncvar_get(ncb,"time")
  
  #timestamp corretion
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
  
  # outll=outll[which(!is.na(outlets)),]
  # outlets=outlets[which(!is.na(outlets))]
  outll=data.frame(outlets,outll)
  return (outll)
}

hydroDir<-("D:/tilloal/Documents/LFRuns_utils/data")
#dataDir<-("D:/tilloal/Documents/LFRuns_utils/data/tss")
dataDir<-("D:/tilloal/Documents/LFRuns_utils/LFPostProcess/Diagnostics")
world <- ne_countries(scale = "medium", returnclass = "sf")
Europe <- world[which(world$continent == "Europe"),]
#lf=list.files(path = paste0(dataDir,"/timeseries/"), full.names = TRUE, recursive = TRUE)
lf=list.files(path = paste0(dataDir), full.names = TRUE, recursive = TRUE)
lf
#lf2=list.files(path = paste0(dataDir,"/timeseries/"), full.names = TRUE, recursive = TRUE)
seasony=function(x){
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
  return(c(Di,R))
}


outletname="outletsv8_hybas07_01min"

outhybas=outletopen(hydroDir,outletname)
outhybas$latlong=paste(round(outhybas$Var1,4),round(outhybas$Var2,4),sep=" ")
#Hybas07
Catchmentrivers7=read.csv(paste0(hydroDir,"/Catchments/from_hybas_eu_onlyid.csv"),encoding = "UTF-8", header = T, stringsAsFactors = F)
hybas07 <- read_sf(dsn = paste0(hydroDir,"/Catchments/hydrosheds/hybas_eu_lev07_v1c.shp"))
hybasf7=fortify(hybas07) 
Catamere07=inner_join(hybasf7,Catchmentrivers7,by= "HYBAS_ID")
Catamere07$llcoord=paste(round(Catamere07$POINT_X,4),round(Catamere07$POINT_Y,4),sep=" ") 
Catf7=inner_join(Catamere07,outhybas,by= c("llcoord"="latlong"))

#Aggregate at catchment level and plot
#Plot parameters
palet=c(hcl.colors(9, palette = "BuPu", alpha = NULL, rev = TRUE, fixup = TRUE))
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

cst7=st_transform(Catf7,  crs=3035)
basemap=w2

#tavg<-read.csv(paste0(hydroDir,"/tss/Tavg/tAvgUps_1954_2020.csv"))
tavg<-read.csv(paste0(hydroDir,"/tss/tAvgUpsv2_1950_2020.csv"))



#reoder tavg according to time
time=as.POSIXct(tavg$X)
tavg=tavg[order(time),]

#if 6-hourly timestep
#timeStampx=as.POSIXct(time*21600-3600,origin = "1950-01-01 00:00:00")
timeStampx=time[order(time)]

#if daily timestep
#timeStampx=as.POSIXct(time*86400-3600,origin = "1950-01-01 00:00:00")

cadegage=which(timeStampx<timeStampx[1])
#timeStamps=as.Date(timeStampx[-cadegage])
tstamp=data.frame(timeStampx)

cname=colnames(tavg)
cname[1]="time"
cname=gsub("X","",cname)
colnames(tavg)=cname
# tavg$tstamp=timeStampx
save=tavg[1,]

tavg$tday=as.Date(tavg$time)
#loop on each pixel to create a file similar my input

tid=seq(1,length(timeStampx))
dsel=hour(timeStampx)
tsday=as.Date(timeStampx[which(dsel==12 | dsel==13)])

zdays=matrix(-9999,ncol=length(cname[-1]),nrow=length(timeStampx))
#zdays2=matrix(-9999,ncol=3,nrow=length(timeStampx)*length(cname[-1]))
mavg=matrix(-9999,ncol=3,nrow=length(tsday)*length(cname[-1]))
#aggregate to daily time scale, this is too much
id=0
for (cn in cname[-1])
{
  id=id+1
  col=which(colnames(tavg)==cn)
  print(col)
  tavg_out=tavg[,col]
  # plot(tavg_out[1:100],type="o")
  #30 days running average
  # tavg_rnm=round(tsEvaNanRunningMean(tavg_out,4),2)
  # points(tavg_rnm[1:100],col=2,type="o")
  
  # tavg_mo1=round(tsEvaNanRunningMean(tavg_out,120),2)
  # points(tavg_mo1[1:100],col=3,type="o")
  
  tavg_mo2=round(tsEvaNanRunningMean(tavg_out,120),2)
  # points(tavg_mo2[1:100],col=4,type="o")
  
  #extract daily values
  
  dataD=tavg_mo2[which(dsel==12 | dsel==13)]
  zdt=tid[which(dsel==12 | dsel==13)]
  lo=length(dataD)
  zdn=as.numeric(rep(cn,length(zdt)))

  td=tid
  zddf=cbind(zdt,zdn,dataD)
  mavg[c(((id-1)*lo+1):(id*lo)),]=zddf
  
  zdays[,id]=tavg_mo2
}

zdays=data.frame(timeStampx, zdays)

names(zdays) = cname

mavg=as.data.frame(mavg)

mavg$time=as.Date(timeStampx[mavg$V1])
mavg$year=year(mavg$time)
AvgFrostAgg=aggregate(list(temp=mavg$V3),
                      by = list(outlet=mavg$V2,yr=mavg$year),
                      FUN = function(x) c(frost=length(which(x<0)),len=length(x),mean=mean(x)))
AvgFrostAgg <- do.call(data.frame, AvgFrostAgg)

AvgFrostAgg$ratpix=(AvgFrostAgg$temp.frost/AvgFrostAgg$temp.len)*100

AvgFAgg=aggregate(list(tf=AvgFrostAgg$temp.frost),
                      by = list(outlet=AvgFrostAgg$outlet),
                      FUN = function(x) c(frost=length(which(x>0)),len=length(x),mean=mean(x)))
AvgFAgg <- do.call(data.frame, AvgFAgg)


names(AvgFAgg)[1]="outlets"
AvgFAgg$frost=AvgFAgg$tf.mean
AvgFAgg$frost[which(AvgFAgg$tf.frost<50)]=0
countries <- read_sf(dsn = paste0(hydroDir,"/Countries/Countries_EFAS.shp"))
catchEFAS <- st_intersection(cst7, countries)
Frostplot=inner_join(AvgFAgg, catchEFAS,by="outlets")

ggplot(basemap) +
  geom_sf(fill="gray85")+
  geom_sf(data=Frostplot,aes(fill=frost,geometry=geometry),color="transparent")+ 
  geom_sf(fill=NA, color="grey") +
  coord_sf(xlim = c(min(nco[,1]),max(nco[,1])), ylim = c(min(nco[,2]),max(nco[,2])))+
  # coord_sf(xlim = c(min(nco[,1]),max(nco[,1])), ylim = c(min(nco[,2]),max(nco[,2])))+
  # geom_point(data=event,aes(x=lon,y=lat,col=dr))+
  # scale_color_distiller(palette = "Spectral",
  #                       direction = -1, limit=c(0,max(eventp$dr)),oob = scales::squish)+
  scale_fill_gradientn(
    colors=palet, limits=c(0,182),
    oob = scales::squish,na.value="transparent", name="days")   +
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
  ggtitle("Average number of frost days per year")


ggsave("D:/tilloal/Documents/LFRuns_utils/Figures/FrostdaysPerYr.jpg", width=30, height=21, units=c("cm"),dpi=1000)
#Now I have to link catchment which have not enough frost seasons 
noFrost=AvgFAgg$outlets[which(AvgFAgg$frost==0)]

colrep=match(noFrost,cname)

#remove the columns without frost
frostcat=zdays[,-colrep]

save(frostcat, file=paste0(hydroDir,"/Drought/catchment_frost.Rdata"))



outlets=unique(AvgFrostAgg$outlet)
plot(AvgFrostAgg$yr[which(AvgFrostAgg$outlet==1029)],AvgFrostAgg$temp.frost[which(AvgFrostAgg$outlet==1029)],type="l")
#



#30 days running averaging for temperature 
id=0
for (cn in cname[-1])
{
  id=id+1
  col=which(colnames(tavg)==cn)
  print(col)
  tavg_out=tavg[,col]
  #30 days running average
  tavg_day=round(tsEvaNanRunningMean(tavg_out,120),2)
  td=tid
  zdn=as.numeric(rep(cn,length(tavg_rnm)))
  lo=length(tavg_rnm)
  
  zdays[,id]=zday_out
}


mavg=as.data.frame(mavg)

mavg$time=timeStampx[mavg$V1]
mavg$year=year(mavg$time)
AvgFrostAgg=aggregate(list(temp=mavg$V3),
                   by = list(outlet=mavg$V2,yr=mavg$year),
                   FUN = function(x) c(frost=length(which(x<0)),len=length(x),mean=mean(x)))
AvgFrostAgg <- do.call(data.frame, AvgFrostAgg)

AvgFrostAgg$ratpix=(AvgFrostAgg$temp.frost/AvgFrostAgg$temp.len)*100
AvgFrostAgg$ratpix[which(AvgFrostAgg$temp.frost<120)]=NA


names(AvgFrostAgg)[1]="outlets"
Frostplot=inner_join(AvgFrostAgg, cst7,by="outlets")


ggplot(basemap) +
  geom_sf(fill="gray85")+
  geom_sf(data=Frostplot,aes(fill=ratpix,geometry=geometry),color="transparent")+ 
  geom_sf(fill=NA, color="grey") +
  coord_sf(xlim = c(min(nco[,1]),max(nco[,1])), ylim = c(min(nco[,2]),max(nco[,2])))+
  # coord_sf(xlim = c(min(nco[,1]),max(nco[,1])), ylim = c(min(nco[,2]),max(nco[,2])))+
  # geom_point(data=event,aes(x=lon,y=lat,col=dr))+
  # scale_color_distiller(palette = "Spectral",
  #                       direction = -1, limit=c(0,max(eventp$dr)),oob = scales::squish)+
  scale_fill_gradientn(
    colors=palet, limits=c(0,50),
    oob = scales::squish,na.value="grey", name="timesteps with frost (%)")   +
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



palet=c(hcl.colors(9, palette = "YlOrRd", alpha = NULL, rev = TRUE, fixup = TRUE))
ggplot(basemap) +
  geom_sf(fill="gray85")+
  geom_sf(data=Frostplot,aes(fill=temp.mean,geometry=geometry),color="transparent")+ 
  geom_sf(fill=NA, color="grey") +
  coord_sf(xlim = c(min(nco[,1]),max(nco[,1])), ylim = c(min(nco[,2]),max(nco[,2])))+
  # coord_sf(xlim = c(min(nco[,1]),max(nco[,1])), ylim = c(min(nco[,2]),max(nco[,2])))+
  # geom_point(data=event,aes(x=lon,y=lat,col=dr))+
  # scale_color_distiller(palette = "Spectral",
  #                       direction = -1, limit=c(0,max(eventp$dr)),oob = scales::squish)+
  scale_fill_gradientn(
    colors=palet, limits=c(-5,20),
    oob = scales::squish,na.value="grey", name="mean 30 days averaged temperature")   +
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




for (cn in cname[-1])
{
  id=id+1
  col=which(colnames(tavg)==cn)
  print(col)
  tavg_out=tavg[,col]
  zday_out=rep(0,length(tavg_out))
  zday_out[which(tavg_out<0)]=1
  
  
  zdt=tid
  zdn=as.numeric(rep(cn,length(zday_out)))
  lo=length(zday_out)
  
  zddf=cbind(zdt,zdn,zday_out)
  
  zdays[,id]=zday_out
  zdays2[c(((id-1)*lo+1):(id*lo)),]=zddf
}

zdays=as.data.frame(zdays)
zdays2=as.data.frame(zdays2)

FrostAgg=aggregate(list(temp=zdays2$V3),
                   by = list(outlet=zdays2$V2),
                   FUN = function(x) c(frost=length(which(x==0)),len=length(x)))
FrostAgg <- do.call(data.frame, FrostAgg)

FrostAgg$ratpix=(1-FrostAgg$temp.frost/FrostAgg$temp.len)*100

outletname="outletsv8_hybas07_01min"

outhybas=outletopen(hydroDir,outletname)
outhybas$latlong=paste(round(outhybas$Var1,4),round(outhybas$Var2,4),sep=" ")
#Hybas07
Catchmentrivers7=read.csv(paste0(hydroDir,"/Catchments/from_hybas_eu_onlyid.csv"),encoding = "UTF-8", header = T, stringsAsFactors = F)
hybas07 <- read_sf(dsn = paste0(hydroDir,"/Catchments/hydrosheds/hybas_eu_lev07_v1c.shp"))
hybasf7=fortify(hybas07) 
Catamere07=inner_join(hybasf7,Catchmentrivers7,by= "HYBAS_ID")
Catamere07$llcoord=paste(round(Catamere07$POINT_X,4),round(Catamere07$POINT_Y,4),sep=" ") 
Catf7=inner_join(Catamere07,outhybas,by= c("llcoord"="latlong"))

#Aggregate at catchment level and plot
#Plot parameters
palet=c(hcl.colors(9, palette = "BuPu", alpha = NULL, rev = TRUE, fixup = TRUE))
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

cst7=st_transform(Catf7,  crs=3035)
basemap=w2
names(FrostAgg)[1]="outlets"
Frostplot=inner_join(FrostAgg, cst7,by="outlets")


ggplot(basemap) +
  geom_sf(fill="gray85")+
  geom_sf(data=Frostplot,aes(fill=ratpix,geometry=geometry),color="transparent")+ 
  geom_sf(fill=NA, color="grey") +
  coord_sf(xlim = c(min(nco[,1]),max(nco[,1])), ylim = c(min(nco[,2]),max(nco[,2])))+
  # coord_sf(xlim = c(min(nco[,1]),max(nco[,1])), ylim = c(min(nco[,2]),max(nco[,2])))+
  # geom_point(data=event,aes(x=lon,y=lat,col=dr))+
  # scale_color_distiller(palette = "Spectral",
  #                       direction = -1, limit=c(0,max(eventp$dr)),oob = scales::squish)+
  scale_fill_gradientn(
    colors=palet, limits=c(0,50),
    oob = scales::squish,na.value="grey", name="timesteps with frost(%)")   +
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





unikY=unique(year(timeStampx))
TavgAggp=c()
for (yrsel in unikY){
  print(yrsel)
  tavgyears=c(year(timeStampx))
  tavg_yr=tavg[which(tavgyears==yrsel),]
  outlets=as.numeric(colnames(tavg_yr)[-1])
  time=tavg_yr$time
  tid=seq(1,length(time))
  #monthly aggregation by catchments
  tavd=data.frame(t(tavg_yr[,-1]))
  tavd=data.frame(outlets,tavd)
  colnames(tavd)=c("outlets",tid)
  
  lo=length(tavd$outlets)
  newTavgY<-melt(tavd,id.vars ="outlets")
  newTavgY$variable=as.numeric(as.character(newTavgY$variable))
  newTavgY$time=time[newTavgY$variable]
  
  newTavgY$belowZ=0
  newTavgY$belowZ[which(newTavgY$value<0)]=1
  
  
  #newTavgY$time=as.POSIXct(newTavgY$variable*21600-3600,origin = "1950-01-01 00:00:00")
  
  newTavgY$month=month(newTavgY$time)
  TavgAgg=aggregate(list(temp=newTavgY$value),
                     by = list(outlet=newTavgY$outlets,month=newTavgY$month),
                     FUN = function(x) c(mean=mean(x,na.rm=T),dev=sd(x,na.rm=T),len=length(x)))
  TavgAgg <- do.call(data.frame, TavgAgg)

  #Only retain catchments with at least 1 month of T<0
  TavgAgg$Frost=0
  TavgAgg$Frost[which(TavgAgg$temp.mean<0)]=1
  TavgAgg$Year=yrsel
  
  TavgAggp=rbind(TavgAggp,TavgAgg)
}
#save(TavgAggp,file=paste0(hydroDir,"/tss/Tavg/Tavg_agg.Rdata"))

#link to hybas07 catchments

load(file=paste0(hydroDir,"/tss/Tavg/Tavg_agg.Rdata"))
outletname="outletsv8_hybas07_01min"

outhybas=outletopen(hydroDir,outletname)
outhybas$latlong=paste(round(outhybas$Var1,4),round(outhybas$Var2,4),sep=" ")
#Hybas07
Catchmentrivers7=read.csv(paste0(hydroDir,"/Catchments/from_hybas_eu_onlyid.csv"),encoding = "UTF-8", header = T, stringsAsFactors = F)
hybas07 <- read_sf(dsn = paste0(hydroDir,"/Catchments/hydrosheds/hybas_eu_lev07_v1c.shp"))
hybasf7=fortify(hybas07) 
Catamere07=inner_join(hybasf7,Catchmentrivers7,by= "HYBAS_ID")
Catamere07$llcoord=paste(round(Catamere07$POINT_X,4),round(Catamere07$POINT_Y,4),sep=" ") 

Catf7=inner_join(Catamere07,outhybas,by= c("llcoord"="latlong"))


#for streams I have the hybas07_raster file now
Gridhybas07=raster( paste0(hydroDir,"/hybas07_raster2.tif"))
Ghybas=as.data.frame(Gridhybas07,xy=T)
Ghybas=Ghybas[which(!is.na(Ghybas[,3])),]
Ghybas$llcoord=paste(round(Ghybas$x,4),round(Ghybas$y,4),sep=" ") 
Ghybas_riv=inner_join(Ghybas,outhybas,by= c("llcoord"="latlong"))


TavgAggF=aggregate(list(frost=TavgAggp$Frost),
                   by = list(outlets=TavgAggp$outlet),
                   FUN = function(x) c(sum=sum(x)))


TavgAggF <- do.call(data.frame, TavgAggF)

TavgAggF$frostpY=TavgAggF$frost/length(unikY)
TavgAggF$Isfrost=TavgAggF$frostpY
TavgAggF$Isfrost[which(TavgAggF$frostpY<1)]=0

TavgAgg_FP=aggregate(list(frost=TavgAggp$Frost),
                   by = list(outlets=TavgAggp$outlet, month=TavgAggp$month),
                   FUN = function(x) c(sum=sum(x), mean=mean(x)))


TavgAgg_FP <- do.call(data.frame, TavgAgg_FP)
TavgAgg_FP$frostp=0
TavgAgg_FP$frostp[which(TavgAgg_FP$frost.mean>0.75)]=1


#save output for the 
save(TavgAgg_FP,file="SeasonFrostDrought.Rdata")


Tagcat=TavgAgg_FP[which(TavgAgg_FP$outlets==1029),]
plot(Tagcat$month,Tagcat$frost.mean,type="l")
hist(TavgAgg_FP$frost.mean)



TavgAgg_jan=TavgAgg_FP[which(TavgAgg_FP$month==1),]

TavgAgg_jan$frostp=0
TavgAgg_jan$frostp[which(TavgAgg_jan$frost.mean>0.75)]=1



TavgAgg_Years=aggregate(list(frost=TavgAggp$Frost),
                     by = list(outlets=TavgAggp$outlet, year=TavgAggp$Year),
                     FUN = function(x) c(sum=sum(x), mean=mean(x)))

TavgAgg_Years <- do.call(data.frame, TavgAgg_Years)

TavgAgg_Years$out=1
TavgAgg_Years$out[which(TavgAgg_Years$frost.sum<1)]=0

TavgAgg_fr=TavgAgg_Years[which(TavgAgg_Years$out==0),]


TavgAgg_frf=aggregate(list(frost=TavgAgg_Years$out),
                        by = list(outlets=TavgAgg_Years$outlet),
                        FUN = function(x) c(sum=sum(x)))

TavgAgg_frf$frostplus=TavgAgg_frf$frost
TavgAgg_frf$frostplus[which(TavgAgg_frf$frost<50)]=0


#trend on number of catchments experiencin frost per year
TavgAgg_trend=aggregate(list(frost=TavgAgg_Years$out),
                      by = list(outlets=TavgAgg_Years$year),
                      FUN = function(x) c(sum=sum(x)))

TavgAgg_trend=TavgAgg_trend[-67,]
plot(TavgAgg_trend,type="l")
abline(lm(frost ~ outlets, data=TavgAgg_trend),col="red")

lareg=lm(frost ~ outlets, data=TavgAgg_trend)
summary(lareg)
#Plot parameters
palet=c(hcl.colors(9, palette = "BuPu", alpha = NULL, rev = TRUE, fixup = TRUE))
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

cst7=st_transform(Catf7,  crs=3035)
basemap=w2

Frostplot=inner_join(TavgAgg_frf, cst7,by="outlets")


ggplot(basemap) +
  geom_sf(fill="gray85")+
  geom_sf(data=Frostplot,aes(fill=frostplus,geometry=geometry),color="transparent")+ 
  geom_sf(fill=NA, color="grey") +
  coord_sf(xlim = c(min(nco[,1]),max(nco[,1])), ylim = c(min(nco[,2]),max(nco[,2])))+
  # coord_sf(xlim = c(min(nco[,1]),max(nco[,1])), ylim = c(min(nco[,2]),max(nco[,2])))+
  # geom_point(data=event,aes(x=lon,y=lat,col=dr))+
  # scale_color_distiller(palette = "Spectral",
  #                       direction = -1, limit=c(0,max(eventp$dr)),oob = scales::squish)+
  scale_fill_gradientn(
    colors=palet, limits=c(50,67),
    oob = scales::squish,na.value="grey", name="# of years with frost")   +
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



Frostplot=inner_join(TavgAggF, cst7,by="outlets")
ggplot(basemap) +
  geom_sf(fill="gray85")+
  geom_sf(data=Frostplot,aes(fill=Isfrost,geometry=geometry),color="transparent")+ 
  geom_sf(fill=NA, color="grey") +
  coord_sf(xlim = c(min(nco[,1]),max(nco[,1])), ylim = c(min(nco[,2]),max(nco[,2])))+
  # coord_sf(xlim = c(min(nco[,1]),max(nco[,1])), ylim = c(min(nco[,2]),max(nco[,2])))+
  # geom_point(data=event,aes(x=lon,y=lat,col=dr))+
  # scale_color_distiller(palette = "Spectral",
  #                       direction = -1, limit=c(0,max(eventp$dr)),oob = scales::squish)+
  scale_fill_gradientn(
    colors=palet, limits=c(0,6),
    oob = scales::squish,na.value="grey", name="average # of frost months ")   +
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

ggsave("FrostMontsperYear.jpg", width=30, height=21, units=c("cm"),dpi=1000)
#create netcdfs containing frost extend *boolean for each month of the year)






Catpixel07=Catamere07[,c(1,13,16,17,18,19,20)]
st_geometry(Catpixel07) <- NULL

wtfhybas=data.frame(as.numeric(unique(Ghybas$hybas07_raster)))
gaaa=data.frame(cst7$HYBAS_ID)




Frostplot=inner_join(TavgAggF, catchEFAS,by="outlets")
ggplot(basemap) +
  geom_sf(fill="gray85")+
  geom_sf(data=Frostplot,aes(fill=Isfrost,geometry=geometry),color="transparent")+ 
  geom_sf(fill=NA, color="grey") +
  coord_sf(xlim = c(min(nco[,1]),max(nco[,1])), ylim = c(min(nco[,2]),max(nco[,2])))+
  # coord_sf(xlim = c(min(nco[,1]),max(nco[,1])), ylim = c(min(nco[,2]),max(nco[,2])))+
  # geom_point(data=event,aes(x=lon,y=lat,col=dr))+
  # scale_color_distiller(palette = "Spectral",
  #                       direction = -1, limit=c(0,max(eventp$dr)),oob = scales::squish)+
  scale_fill_gradientn(
    colors=palet, limits=c(0,6),
    oob = scales::squish,na.value="grey", name="average # of frost months ")   +
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

ggsave("FrostMontsperYear2.jpg", width=30, height=21, units=c("cm"),dpi=1000)



#Monthly maps

Frostplot=inner_join(TavgAgg_jan, catchEFAS,by="outlets")
ggplot(basemap) +
  geom_sf(fill="gray85")+
  geom_sf(data=Frostplot,aes(fill=frostp,geometry=geometry),color="transparent")+ 
  geom_sf(fill=NA, color="grey") +
  coord_sf(xlim = c(min(nco[,1]),max(nco[,1])), ylim = c(min(nco[,2]),max(nco[,2])))+
  # coord_sf(xlim = c(min(nco[,1]),max(nco[,1])), ylim = c(min(nco[,2]),max(nco[,2])))+
  # geom_point(data=event,aes(x=lon,y=lat,col=dr))+
  # scale_color_distiller(palette = "Spectral",
  #                       direction = -1, limit=c(0,max(eventp$dr)),oob = scales::squish)+
  scale_fill_gradientn(
    colors=palet, limits=c(0,1),
    oob = scales::squish,na.value="grey", name="average # of frost months ")   +
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


Catpixel07F=inner_join(Frostplot,Ghybas,by=c("SORT"= "hybas07_raster2"))

Catpixel07F=Catpixel07F[,-17]
cp07 <- st_as_sf(Catpixel07F, coords = c("x", "y"), crs = 4326)
cp07 <- st_transform(cp07, crs = 3035)

#catchEFAS <- st_intersection(cp07, countries)

#this is working, now I have to save this shit as netcdf

tsize=14
osize=12
ggplot(basemap) +
  geom_sf(fill="gray90")+
  geom_sf(data=cp07,aes(col=frostp,geometry=geometry),alpha=1,size=.1,stroke=0,shape=15)+
  geom_sf(fill=NA, color="grey") +
  coord_sf(xlim = c(min(nco[,1]),max(nco[,1])), ylim = c(min(nco[,2]),max(nco[,2])))+
  scale_color_gradientn(
    colors=palet, limits=c(0,1),
    oob = scales::squish,na.value="grey", name="# of years with frost")+
  theme(axis.title=element_text(size=tsize),
        panel.background = element_rect(fill = "aliceblue", colour = "grey1"),
        panel.border = element_rect(linetype = "solid", fill = NA, colour="black"),
        legend.title = element_text(size=tsize),
        legend.text = element_text(size=osize),
        legend.position = "right",
        panel.grid.major = element_line(colour = "grey70"),
        panel.grid.minor = element_line(colour = "grey90"),
        legend.key = element_rect(fill = "transparent", colour = "transparent"),
        legend.key.size = unit(.8, "cm"))


TavgAgg_FP$frostp=0
TavgAgg_FP$frostp[which(TavgAgg_FP$frost.mean>0.75)]=1

TavgAgg_Final=aggregate(list(frost=TavgAgg_FP$frostp),
                     by = list(outlets=TavgAgg_FP$outlet),
                     FUN = function(x) c(sum=sum(x), mean=mean(x)))


TavgAgg_Final <- do.call(data.frame, TavgAgg_Final)



Frostplot=inner_join(TavgAgg_Final, catchEFAS,by="outlets")
ggplot(basemap) +
  geom_sf(fill="gray85")+
  geom_sf(data=Frostplot,aes(fill=frost.sum,geometry=geometry),color="transparent")+ 
  geom_sf(fill=NA, color="grey") +
  coord_sf(xlim = c(min(nco[,1]),max(nco[,1])), ylim = c(min(nco[,2]),max(nco[,2])))+
  # coord_sf(xlim = c(min(nco[,1]),max(nco[,1])), ylim = c(min(nco[,2]),max(nco[,2])))+
  # geom_point(data=event,aes(x=lon,y=lat,col=dr))+
  # scale_color_distiller(palette = "Spectral",
  #                       direction = -1, limit=c(0,max(eventp$dr)),oob = scales::squish)+
  scale_fill_gradientn(
    colors=palet, limits=c(0,6),breaks=c(6:0), labels=rev(c("0","1","2","3","4","5",">=6")),
    oob = scales::squish,na.value="grey", name="Number of frost months ", guide="legend")   +
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

ggsave("Figures/FrostMonths_ana_th0.75.jpg", width=30, height=21, units=c("cm"),dpi=1000)




dom=domainopen(hydroDir,outletname)
dom$latlong=paste(round(dom$Var1,4),round(dom$Var2,4),sep=" ")
ncfile_dom = paste0(hydroDir,"/upArea_European_01min.nc")
nc=nc_open(ncfile_dom)

#1. I load one subregion
name.lon="lon"
name.lat="lat"
londat = ncvar_get(nc,name.lon) 
latdat = ncvar_get(nc,name.lat) 



nlon=length(unique(dom$idlo))
nlat=length(unique(dom$idla))
timedim <- ncdim_def("time","Months",as.double(c(1:12)))
londim <- ncdim_def("lon","degrees_east",as.double(londat)) 
latdim <- ncdim_def("lat","degrees_north",as.double(latdat)) 
tmp_def <- ncvar_def("frost","boolean",list(londim,latdim,timedim),-9999,"frost months",prec="short",compression=4)
proj <-ncvar_def("wgs_1984","1",NULL,NULL,longname="wgs_1984",prec="integer")
cat(paste0("\ Initializing netcdf \n"))
ncname <- paste0("FrostSeasons_Europe.nc")
ncpath<-hydroDir
ncfname <- paste(ncpath, ncname, sep="/")
ncout <- nc_create(ncfname,list(tmp_def,proj),force_v4=TRUE)  

for (mo in 1:12){
  print(mo)
  start_time <- Sys.time()
  TavgAgg_mon=TavgAgg_FP[which(TavgAgg_FP$month==mo),]
  TavgAgg_mon$frostp=0
  TavgAgg_mon$frostp[which(TavgAgg_mon$frost.mean>0.75)]=1
  Frostplot=inner_join(TavgAgg_mon, catchEFAS,by="outlets")
  Frostplot=Frostplot[,-41]
  Catpixel07F=inner_join(Frostplot,Ghybas,by=c("SORT"= "hybas07_raster2"))
  #cp07 <- st_as_sf(Catpixel07F, coords = c("x", "y"), crs = 4326)
  fixshit=match(dom$latlong,Catpixel07F$llcoord.y)
  coordstuff=dom$latlong[fixshit]
  dom$frost=Catpixel07F$frostp[fixshit]
  
  # frosty=full_join(Catpixel07F,dom,by=c("llcoord.x"="latlong"))
  # jtemerde=inner_join(frosty,dom,by=c("llcoord.x"="latlong"))
  # 
  # frosty=frosty[with(frosty, order(idlo.y,idla.y)),]
  frost_array <- array(dom$frost, dim=c(nlon,nlat))
  sts=c(1,1,mo)
  cts=c(nlon,nlat,1)
  ncvar_put(ncout,tmp_def, frost_array, start = sts, count=cts)
  end_time <- Sys.time()
  cat(paste0("\ duration: ",round(end_time-start_time,4)," seconds\n"))
}


# library(lattice)
# levelplot(frost_array)

ncatt_put(ncout,"lon","axis","X") #,verbose=FALSE) #,definemode=FALSE)
ncatt_put(ncout,"lat","axis","Y")

ncatt_put(ncout,tmp_def,"standard_name","frost days")
ncatt_put(ncout,tmp_def,"grid_mapping","wgs_1984")
ncatt_put(ncout,tmp_def,"esri_pe_string",'GEOGCS["WGS 84",DATUM["WGS_1984",SPHEROID["WGS 84",6378137,298.257223563,AUTHORITY["EPSG","7030"]],AUTHORITY["EPSG","6326"]],PRIMEM["Greenwich",0],UNIT["degree",0.0174532925199433],AUTHORITY["EPSG","4326"]]')
ncatt_put(ncout,tmp_def,"missing_value",-9999)



# put the CRS attributes
projname <- "wgs_1984"
ncatt_put(ncout,proj,"name",projname)
ncatt_put(ncout,proj,"long_name",projname)
ncatt_put(ncout,proj,"grid_mapping_name","latitude_longitude")
ncatt_put(ncout,proj,"semi_major_axis", 6378137)
ncatt_put(ncout,proj,"inverse_flattening", 298.257223563)
ncatt_put(ncout,proj,"proj4_params", "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")
ncatt_put(ncout,proj,"EPSG_code","EPSG:4326")


history = 'Created Sept 2023' #####
Conventions = 'CF-1.6'
Source_Software = 'R netCDF4'
reference = 'JRC Climate risk team'  #####
title = 'Frost monthly maps (1954-2019)'
keywords = 'Lisflood, Global'
source = 'ERA5-land'
institution = 'European Commission - Economics of climate change Unit (JRC.C.6) : https://ec.europa.eu/jrc/en/research-topic/climate-change'
comment = 'no.'

# add global attributes
ncatt_put(ncout,0,"title",title)
ncatt_put(ncout,0,"institution",institution)
ncatt_put(ncout,0,"source",source)
ncatt_put(ncout,0,"references",reference)
history <- paste("A.M. Tilloy", date(), sep=", ")
ncatt_put(ncout,0,"history",history)
ncatt_put(ncout,0,"Conventions",Conventions)

# Get a summary of the created file:
ncout

nc_close(ncout)




outletCNT=pixfrance[,c(1:4)]
st_geometry(outletCNT)=NULL
#loop over all the squares to get the ID corresponding to the pixels
rspace= read.csv(paste0(hydroDir,"/subspace_efas.csv"))
rspace=rspace[,-1]



#This netcdf will be loaded during TSEVA







timeseq=data.frame(time1)

#Loop to generate yearly netcdf files for all the domain
for (idt in 1:length(years)){
  
  cat(paste0("\ year ",years[idt],"\n"))
  timeY=seq(as.Date(paste0(years[idt],"-01-01")),as.Date(paste0(years[idt],"-12-31")),by="day")
  timeY=data.frame(time=timeY)
  nt=length(timeY$time)
  td=as.numeric(timeY$time-as.Date("1950-01-01"))
  timedim <- ncdim_def("time",tunitsf,as.double(td))
  londim <- ncdim_def("lon","degrees_east",as.double(londat)) 
  latdim <- ncdim_def("lat","degrees_north",as.double(latdat)) 
  tmp_def <- ncvar_def("flow","boolean",list(londim,latdim,timedim),-9999,"flow days",prec="double",compression=4)
  proj <-ncvar_def("wgs_1984","1",NULL,NULL,longname="wgs_1984",prec="integer")
  cat(paste0("\ Initializing netcdf \n"))
  ncname <- paste0("Intermittent_Rivers_Bool_",years[idt],".nc")
  ncpath<-hydroDir
  ncfname <- paste(ncpath, ncname, sep="/")
  ncout <- nc_create(ncfname,list(tmp_def,proj),force_v4=TRUE)  
  
  zeryear=Zertot[which(Zertot$yr==years[idt]),]
  #fill the netcdf at a daily timescale
  for (idtx in 1:nt){
    start_time <- Sys.time()
    
    cat(paste0("\ day ",timeY$time[idtx],"\n"))
    if (length(zeryear>0)){
      zerday=zeryear[which(zeryear$time==timeY$time[idtx]),]
      zerday=full_join(zerday,outf,by=c("catch"="outlets"))
      zerday$flow=1
      zerday$flow[which(!is.na(zerday$zerod))]=0
    }else{
      zerday=outhybas
      zerday$flow=1
    }
    
    
    zerfinal=full_join(zerday,dom,by=c("latlong"="latlong"))
    zerfinal=zerfinal[with(zerfinal, order(idla.y,idlo.y)),]
    
    nlon=length(unique(zerfinal$idlo.y))
    nlat=length(unique(zerfinal$idla.y))
    zer_array <- array(zerfinal$flow, dim=c(nlon,nlat))
    sts=c(1,1,idtx)
    cts=c(nlon,nlat,1)
    ncvar_put(ncout,tmp_def, zer_array, start = sts, count=cts)
    end_time <- Sys.time()
    cat(paste0("\ duration: ",round(end_time-start_time,4)," seconds\n"))
    
  }
  
  
  ncatt_put(ncout,"lon","axis","X") #,verbose=FALSE) #,definemode=FALSE)
  ncatt_put(ncout,"lat","axis","Y")
  
  ncatt_put(ncout,tmp_def,"standard_name","flow days")
  ncatt_put(ncout,tmp_def,"grid_mapping","wgs_1984")
  ncatt_put(ncout,tmp_def,"esri_pe_string",'GEOGCS["WGS 84",DATUM["WGS_1984",SPHEROID["WGS 84",6378137,298.257223563,AUTHORITY["EPSG","7030"]],AUTHORITY["EPSG","6326"]],PRIMEM["Greenwich",0],UNIT["degree",0.0174532925199433],AUTHORITY["EPSG","4326"]]')
  ncatt_put(ncout,tmp_def,"missing_value",-9999)
  
  
  
  # put the CRS attributes
  projname <- "wgs_1984"
  ncatt_put(ncout,proj,"name",projname)
  ncatt_put(ncout,proj,"long_name",projname)
  ncatt_put(ncout,proj,"grid_mapping_name","latitude_longitude")
  ncatt_put(ncout,proj,"semi_major_axis", 6378137)
  ncatt_put(ncout,proj,"inverse_flattening", 298.257223563)
  ncatt_put(ncout,proj,"proj4_params", "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")
  ncatt_put(ncout,proj,"EPSG_code","EPSG:4326")
  
  
  history = 'Created Feb 2023' #####
  Conventions = 'CF-1.6'
  Source_Software = 'R netCDF4'
  reference = 'JRC Climate risk team'  #####
  title = 'Riverflow maps'
  keywords = 'Lisflood, Global'
  source = 'ERA5-land'
  institution = 'European Commission - Economics of climate change Unit (JRC.C.6) : https://ec.europa.eu/jrc/en/research-topic/climate-change'
  comment = 'no.'
  
  # add global attributes
  ncatt_put(ncout,0,"title",title)
  ncatt_put(ncout,0,"institution",institution)
  ncatt_put(ncout,0,"source",source)
  ncatt_put(ncout,0,"references",reference)
  history <- paste("A.M. Tilloy", date(), sep=", ")
  ncatt_put(ncout,0,"history",history)
  ncatt_put(ncout,0,"Conventions",Conventions)
  
  # Get a summary of the created file:
  ncout
  
  nc_close(ncout)
  
  # for (id in 1:length(unikzero)){
  #   print(id)
  #   Rpix=Zerall[which(Zerall$catch == unikzero[id] ),]
  #   Rpix$time=as.Date(Rpix$td,origin="1951-01-01")
  #   Rpix=full_join(Rpix,timeseq,by="time")
  #   Rpix$catch=unikzero[id]
  #   Rpix$flow=1
  #   Rpix$flow[which(!is.na(Rpix$zerod))]=0
  #   Rpix=Rpix[,c(5,6,7)]
  #   #pixdat=inner_join(Rpix,outhybas,by=c("catch"="outlets"))
  #   pixlite=outhybas[match(Rpix$catch,outhybas$outlets),]
  #   lonp=unique(pixlite$Var1)
  #   latp=unique(pixlite$Var2)
  #   
  #   locx=which(!is.na(match(londat,lonp)))
  #   latx=which(!is.na(match(latdat,latp)))
  #   sts=c(locx,latx,1)
  #   cts=c(1,1,length(time))
  #   datfill=Rpix$flow
  #   datfillb=as.integer(datfill)
  #   Rpg=rbind(Rpg,Rpix)
  #   start_time <- Sys.time()
  #   ncvar_put(ncout,tmp_def, datfill, start = sts, count=cts)
  #   end_time <- Sys.time()
  #   cat(paste0("\ duration: ",round(end_time-start_time,4)," seconds\n"))
  # }
}






