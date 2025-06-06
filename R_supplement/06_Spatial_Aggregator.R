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



#Functions
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



#Set data directory
hydroDir<-("D:/tilloal/Documents/LFRuns_utils/data")

#load outf

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

#Load my shapefile on which to aggregate

### Hybas07 ----
Catchmentrivers7=read.csv(paste0(hydroDir,"/Catchments/from_hybas_eu_onlyid.csv"),encoding = "UTF-8", header = T, stringsAsFactors = F)
hybas07 <- read_sf(dsn = paste0(hydroDir,"/Catchments/hydrosheds/hybas_eu_lev07_v1c.shp"))
hybasf7=fortify(hybas07)
Catamere07=inner_join(hybasf7,Catchmentrivers7,by= "HYBAS_ID")
Catamere07$llcoord=paste(round(Catamere07$POINT_X,4),round(Catamere07$POINT_Y,4),sep=" ")
cst7=right_join(Catamere07,outf,by= c("llcoord"="latlong"))
GNF=cst7
st_geometry(GNF)=NULL

UnHY=unique(GNF$HYBAS_ID)
# 
#cst7=st_transform(cst7,  crs=3035)


### HydroRegions ----

GridHR=raster( paste0(hydroDir,"/HydroRegions_raster_WGS84.tif"))
GHR=as.data.frame(GridHR,xy=T)
GHR=GHR[which(!is.na(GHR[,3])),]
GHR$llcoord=paste(round(GHR$x,4),round(GHR$y,4),sep=" ")
GHR_riv=inner_join(GHR,outf,by= c("llcoord"="latlong"))
GHshpp <- read_sf(dsn ="Z:/ClimateRun4/nahaUsers/tilloal/HydroRegions/her_all_adjusted.shp")
HydroRsf=fortify(GHshpp)

UnHY=unique(GHR_riv$HydroRegions_raster_WGS84)



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

#load data to be aggregated



### Plot parameters ----
palet2=c(hcl.colors(9, palette = "Blues", alpha = NULL, rev = TRUE, fixup = TRUE))
# Plot of ordered change by region, can be important
palet=c(hcl.colors(11, palette = "RdYlBu", alpha = NULL, rev = F, fixup = TRUE))
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

#Land use
library(exactextractr)
rastforest=raster( "D:/tilloal/Documents/06_Floodrivers/landuse/fracforest_ch20201951.tif")

rastsealed=raster( "D:/tilloal/Documents/06_Floodrivers/landuse/fracsealed_ch20201951.tif")

rastirrigated=raster( "D:/tilloal/Documents/06_Floodrivers/landuse/fracirrigated_ch20201951.tif")

rastother=raster( "D:/tilloal/Documents/06_Floodrivers/landuse/fracother_ch20201951.tif")

rastrice=raster( "D:/tilloal/Documents/06_Floodrivers/landuse/fracrice_ch20201951.tif")

rastwater=raster( "D:/tilloal/Documents/06_Floodrivers/landuse/fracwater_ch20201951.tif")

#I could do a relative change as well

forestchange<- exact_extract(rastforest, GHshpp, 'mean')
GHshpp$forestchange=forestchange
GHshpp$sealedchange <- exact_extract(rastsealed, GHshpp, 'mean')
GHshpp$irrigatedchange <- exact_extract(rastirrigated, GHshpp, 'mean')
GHshpp$otherchange <- exact_extract(rastother, GHshpp, 'mean')
GHshpp$ricechange <- exact_extract(rastrice, GHshpp, 'mean')
GHshpp$waterchange <- exact_extract(rastwater, GHshpp, 'mean')

mhh=match(UnHY,GHshpp$Id)
GHshppH=GHshpp[mhh,]
df_GHshppH=data.frame(GHshppH)
st_geometry(df_GHshppH)<-NULL


lucmap=list()
luclass=c("forest","sealed","irrigated","other","rice","water")

for (li in 1:length(luclass))
{
  lu=luclass[li]
  print(lu)
  GHshppH$fill=as.numeric(df_GHshppH[,7+li])
  
  flims=(quantile(GHshppH$fill,c(0.01,0.99),na.rm=T))
  lims=c(-round(max(abs(flims)),1),round(max(abs(flims)),1))
  if (diff(lims)==0){
    lims=c(-round(max(abs(flims)),2),round(max(abs(flims)),2))
  }
  if (diff(lims)==0){
    lims=c(-round(max(abs(flims)),4),round(max(abs(flims)),4))
  }
  #lims=c(-25,25)
  lucmap<-ggplot(basemap) +
    geom_sf(fill="gray95",color="transparent",size=0.5)+
    geom_sf(data=GHshppH,aes(fill=fill*100,geometry=geometry),alpha=1,color="transparent")+
    geom_sf(fill="transparent",color="gray30",size=0.5)+
    coord_sf(xlim = c(min(nco[,1]),max(nco[,1])), ylim = c(min(nco[,2]),max(nco[,2])))+
    scale_fill_gradientn(
      colors=palet,
      limits=lims*100,oob = scales::squish,
      name=paste0("Change in ",lu,"(%)"))   +
    labs(x="Longitude", y = "Latitude")+
    theme(axis.title=element_text(size=tsize),
          title = element_text(size=16),
          axis.text=element_text(size=osize),
          panel.background = element_rect(fill = "aliceblue", colour = "grey1"),
          panel.border = element_rect(linetype = "solid", fill = NA, colour="black"),
          legend.title = element_text(size=tsize),
          legend.text = element_text(size=osize),
          legend.position = "right",
          panel.grid.major = element_line(colour = "grey70"),
          panel.grid.minor = element_line(colour = "grey90"),
          legend.key = element_rect(fill = "transparent", colour = "transparent"),
          legend.key.size = unit(1, "cm"))+
    ggtitle(lu)
  ggsave(paste0("D:/tilloal/Documents/LFRuns_utils/TAplots/luchange_",lu,"_20201951_HR.jpg"), lucmap,width=30, height=20, units=c("cm"),dpi=300)

}



#Same for water demand change



#divide by area
workDir<-("D:/tilloal/Documents/06_Floodrivers/")
ncpix= paste0(workDir,"/mapscal/pixarea_European_01min.nc")
ncp=nc_open(ncpix)
t=ncp$var[[2]]
tsize<-t$varsize
pixarea=ncvar_get(ncp,names(ncp[['var']])[2]) 
plon=ncvar_get(ncp,"lon") 
plat=ncvar_get(ncp,"lat") 
pixarea=as.matrix(t(pixarea))

rast_totwd<-raster( "D:/tilloal/Documents/06_Floodrivers/wateruse/wateruse_sums/all_demands_2020.tif")
rast_totwd<-raster( "D:/tilloal/Documents/06_Floodrivers/wateruse/wateruse_sums/all_demands_1951.tif")
#Convert RasterLayer to matrix
rast_tmat <- (as.matrix(rast_totwd))
#multiply by 30.4 to have the real sum values
rast_tmat=rast_tmat*30.4
#m3/m2 to m3
rast_tmat=rast_tmat*pixarea
#convert mm to m3
rast_tmat=rast_tmat/1e3
#m3 to km3
rast_tmat=rast_tmat/1e9
rast_tmout <- raster(nrows=nrow(rast_totwd), ncols=ncol(rast_totwd), ext=extent(rast_totwd))
crs(rast_tmout) <- crs(rast_tmat)
values(rast_tmout) <- rast_tmat
rast_totwd=rast_tmout

HRwgs84 <- st_transform(GHshpp, crs = 4326)
HRwgs84$totalwd2020 <- exact_extract(rast_totwd, HRwgs84, 'sum')
sum(HRwgs84$totalwd2020)
HRwgs84$wdpkm2=HRwgs84$totalwd2020/HRwgs84$SURF_KM2*1000*1000




mhh=match(UnHY,GHshpp$HYBAS_ID)
hybas07$totalwd2020 <- exact_extract(rast_totwd, hybas07, 'sum')
hybas07H=hybas07[mhh,]
sum(hybas07H$totalwd2020,na.rm=T)
hybas07H$wdpkm2=hybas07H$totalwd2020/hybas07H$UP_AREA*1000*1000

#unit is now in mm
wd="total water demand in 1951"
tsize=12
flims=(quantile(HRwgs84$wdpkm2,c(0.1,0.95),na.rm=T))
lims=c(0,round(max(abs(flims)),1))
lims=c(0,120)
wdmap<-ggplot(basemap) +
  geom_sf(fill="gray95",color="transparent",size=0.5)+
  geom_sf(data=HRwgs84,aes(fill=wdpkm2,geometry=geometry),alpha=1,color="transparent")+
  geom_sf(fill="transparent",color="gray30",size=0.5)+
  coord_sf(xlim = c(min(nco[,1]),max(nco[,1])), ylim = c(min(nco[,2]),max(nco[,2])))+
  scale_fill_gradientn(
    colors=palet2,
    limits=lims,oob = scales::squish, trans="sqrt", breaks=c(1,5,20,50,100,200),
    name=paste0("(mm/year)"))   +
  labs(x="Longitude", y = "Latitude")+
  theme(axis.title=element_text(size=tsize),
        title = element_text(size=16),
        axis.text=element_text(size=osize),
        panel.background = element_rect(fill = "aliceblue", colour = "grey1"),
        panel.border = element_rect(linetype = "solid", fill = NA, colour="black"),
        legend.title = element_text(size=tsize),
        legend.text = element_text(size=osize),
        legend.position = "right",
        panel.grid.major = element_line(colour = "grey70"),
        panel.grid.minor = element_line(colour = "grey90"),
        legend.key = element_rect(fill = "transparent", colour = "transparent"),
        legend.key.size = unit(1, "cm"))+
  ggtitle(wd)

wdmap
ggsave(paste0("D:/tilloal/Documents/LFRuns_utils/TAplots/TotalWaterDemand_1951.jpg"), wdmap,width=30, height=20, units=c("cm"),dpi=300)



rast_ene=raster( "D:/tilloal/Documents/06_Floodrivers/wateruse/ene_ysum_ch20201951.tif")
#Convert RasterLayer to matrix
rast_enemat <- (as.matrix(rast_ene))
#multiply by 30.4 to have the real sum values
rast_enemat=rast_enemat*30.4
#m3/m2 to m3
rast_enemat=rast_enemat*pixarea
#convert mm to m3
rast_enemat=rast_enemat/1e3
#m3 to km3
rast_enemat=rast_enemat/1e9
rast_eneout <- raster(nrows=nrow(rast_ene), ncols=ncol(rast_ene), ext=extent(rast_ene))
crs(rast_eneout) <- crs(rast_ene)
values(rast_eneout) <- rast_enemat
rast_ene=rast_eneout

rast_dom=raster( "D:/tilloal/Documents/06_Floodrivers/wateruse/dom_ysum_ch20201951.tif")
#Convert RasterLayer to matrix
rast_dommat <- (as.matrix(rast_dom))
#multiply by 30.4 to have the real sum values
rast_dommat=rast_dommat*30.4
#m3/m2 to m3
rast_dommat=rast_dommat*pixarea
#convert mm to m3
rast_dommat=rast_dommat/1e3
#m3 to km3
rast_dommat=rast_dommat/1e9
rast_domout <- raster(nrows=nrow(rast_dom), ncols=ncol(rast_dom), ext=extent(rast_dom))
crs(rast_domout) <- crs(rast_dom)
values(rast_domout) <- rast_dommat
rast_dom=rast_domout
rast_liv=raster( "D:/tilloal/Documents/06_Floodrivers/wateruse/liv_ysum_ch20201951.tif")
#Convert RasterLayer to matrix
rast_livmat <- (as.matrix(rast_liv))
#multiply by 30.4 to have the real sum values
rast_livmat=rast_livmat*30.4
#m3/m2 to m3
rast_livmat=rast_livmat*pixarea
#convert mm to m3
rast_livmat=rast_livmat/1e3
#m3 to km3
rast_livmat=rast_livmat/1e9
rast_livout <- raster(nrows=nrow(rast_liv), ncols=ncol(rast_liv), ext=extent(rast_liv))
crs(rast_livout) <- crs(rast_liv)
values(rast_livout) <- rast_livmat
rast_liv=rast_livout

rast_ind=raster( "D:/tilloal/Documents/06_Floodrivers/wateruse/ind_ysum_ch20201951.tif")
#Convert RasterLayer to matrix
rast_indmat <- (as.matrix(rast_ind))
#multiply by 30.4 to have the real sum values
rast_indmat=rast_indmat*30.4
#m3/m2 to m3
rast_indmat=rast_indmat*pixarea
#convert mm to m3
rast_indmat=rast_indmat/1e3
#m3 to km3
rast_indmat=rast_indmat/1e9
rast_indout <- raster(nrows=nrow(rast_ind), ncols=ncol(rast_ind), ext=extent(rast_ind))
crs(rast_indout) <- crs(rast_ind)
values(rast_indout) <- rast_indmat
rast_ind=rast_indout

rast_total=raster( "D:/tilloal/Documents/06_Floodrivers/wateruse/all_ysum_ch20201951.tif")
#Convert RasterLayer to matrix
rast_totalmat <- (as.matrix(rast_total))
#multiply by 30.4 to have the real sum values
rast_totalmat=rast_totalmat*30.4
#m3/m2 to m3
rast_totalmat=rast_totalmat*pixarea
#convert mm to m3
rast_totalmat=rast_totalmat/1e3
#m3 to km3
rast_totalmat=rast_totalmat/1e9
rast_totalout <- raster(nrows=nrow(rast_total), ncols=ncol(rast_total), ext=extent(rast_total))
crs(rast_totalout) <- crs(rast_total)
values(rast_totalout) <- rast_totalmat
rast_total=rast_totalout


HRwgs84$enechange <- exact_extract(rast_ene, HRwgs84, 'sum')
HRwgs84$domchange <- exact_extract(rast_dom, HRwgs84, 'sum')
HRwgs84$livchange <- exact_extract(rast_liv, HRwgs84, 'sum')
HRwgs84$indchange <- exact_extract(rast_ind, HRwgs84, 'sum')
HRwgs84$totalchange <- exact_extract(rast_total, HRwgs84, 'sum')


mhh=match(UnHY,HRwgs84$Id)
HRwgs84h=HRwgs84[mhh,]
df_HRwgs84h=data.frame(HRwgs84h)
st_geometry(df_HRwgs84h)<-NULL

# mhh=match(UnHY,hybas07$HYBAS_ID)
# hybas07H=hybas07[mhh,]
# df_hybas07H=data.frame(hybas07H)
# st_geometry(df_hybas07H)<-NULL

wdclass=c("ene","dom","liv","ind","total")
wdmap=list()
for (li in 1:length(wdclass))
{
  wd=wdclass[li]
  print(wd)
  HRwgs84h$fill=as.numeric(df_HRwgs84h[,15+li])
  
  flims=(quantile(HRwgs84h$fill,c(0.01,0.99),na.rm=T))
  lims=c(-round(max(abs(flims)),1),round(max(abs(flims)),1))
  if (diff(lims)==0){
    lims=c(-round(max(abs(flims)),2),round(max(abs(flims)),2))
  }
  if (diff(lims)==0){
    lims=c(-round(max(abs(flims)),4),round(max(abs(flims)),4))
  }
  tsize=12
  wdmap[[li]]<-ggplot(basemap) +
    geom_sf(fill="gray95",color="transparent",size=0.5)+
    geom_sf(data=HRwgs84h,aes(fill=fill,geometry=geometry),alpha=1,color="transparent")+
    geom_sf(fill="transparent",color="gray30",size=0.5)+
    coord_sf(xlim = c(min(nco[,1]),max(nco[,1])), ylim = c(min(nco[,2]),max(nco[,2])))+
    scale_fill_gradientn(
      colors=palet,
      limits=lims,oob = scales::squish,
      name=paste0("Change (km3/year)"))   +
    labs(x="Longitude", y = "Latitude")+
    theme(axis.title=element_text(size=tsize),
          title = element_text(size=16),
          axis.text=element_text(size=osize),
          panel.background = element_rect(fill = "aliceblue", colour = "grey1"),
          panel.border = element_rect(linetype = "solid", fill = NA, colour="black"),
          legend.title = element_text(size=tsize),
          legend.text = element_text(size=osize),
          legend.position = "right",
          panel.grid.major = element_line(colour = "grey70"),
          panel.grid.minor = element_line(colour = "grey90"),
          legend.key = element_rect(fill = "transparent", colour = "transparent"),
          legend.key.size = unit(1, "cm"))+
    ggtitle(wd)
  ggsave(paste0("D:/tilloal/Documents/LFRuns_utils/TAplots/wdchangeHR_",wd,"_20201951.jpg"),wdmap[[li]],width=30, height=20, units=c("cm"),dpi=300)
  
}

wdmap[[1]]

###################################################################
#detailed analysis of water demand (seasonal)

#load netcdf of water demand

print(year)
df_wd=c()
for (w in wdclass[-5]){
  
  nc <- nc_open(paste0("D:/tilloal/Documents/06_Floodrivers/wateruse/wateruse_histo/",w,"_1950_2020.nc"))
  print(w)
  # Extract the dimensions of the file
  lon_dim <- nc$dim[["lon"]]
  lat_dim <- nc$dim[["lat"]]
  time_dim <- nc$dim[["time"]]
  name.vb=names(nc[['var']])
  namev=name.vb[1]
  time <- ncvar_get(nc,"time")
  timestamp=as.Date(time,origin="1950-01-01")
  lt=length(time)
  
  
  # Extract the longitude and latitude values
  lon <- ncvar_get(nc, "lon")
  lat <- ncvar_get(nc, "lat")
  
  outll=expand.grid(lon,lat)
  mvs=c()
  for (d in 1:lt){
    print(d/lt*100)
    mon_val <- ncvar_get(nc, namev, 
                           start = c(1, 1, d), 
                           count = c(lon_dim$len, lat_dim$len, 1))
    md=sum(mon_val,na.rm=T)
    
    mvs=c(mvs,md)
  }

  df_wd=rbind(df_wd,mvs)
}

df_wd=as.data.frame(t(df_wd))
names(df_wd)=wdclass[-5]
df_wd$sum=df_wd$ind+df_wd$liv+df_wd$dom+df_wd$ene
df_wd$time=timestamp

save(df_wd,file=paste0("D:/tilloal/Documents/06_Floodrivers/wateruse/wd_domain.Rdata"))


# Create a data frame from the first day
dx <- data.frame(lon = as.vector(outll$Var1), 
                 lat = as.vector(outll$Var2), 
                 value = as.vector(first_day))

####################################################################
#load results from drought catchment analysis 

#floods
load(file=paste0(hydroDir,"/TSEVA/output_plots/outputs_flood_year_qsp.Rdata"))
#droughts
load(file=paste0(hydroDir,"/TSEVA/output_plots/outputs_drought_nonfrost_qsp.Rdata"))

#I extract the trend at MUTS3 level fist

FloodTrends=Output_fl_year$TrendOutlets
# mhh=match(UnHY,FloodTrends$HydroR)
FloodTrends=FloodTrends[which(FloodTrends$driver=="Landuse"),]
length(which(!is.na(FloodTrends$Y2020)))
Output_dr_nfrost=Output_dr_year
DroughtTrends=Output_dr_nfrost$TrendOutlets
DroughtTrends=DroughtTrends[which(DroughtTrends$driver=="Landuse"),]


#join flood and land use changes
#Flood_xplain=inner_join(FloodTrends,df_hybas07H,by=c("HydroR"="HYBAS_ID"))
Flood_xplain=inner_join(FloodTrends,df_hybas07H,by=c("HYBAS_ID"))

plot(Flood_xplain$Rchange.Y2020,Flood_xplain$sealedchange)
corsealed=cor(Flood_xplain$Rchange.Y2020,Flood_xplain$sealedchange)

plot(Flood_xplain$Rchange.Y2020,Flood_xplain$forestchange)
corforest=cor(Flood_xplain$Rchange.Y2020,Flood_xplain$forestchange)

plot(Flood_xplain$Rchange.Y2020,Flood_xplain$sealedchange)
corsealed=cor(Flood_xplain$Rchange.Y2020,Flood_xplain$sealedchange)

plot(Flood_xplain$Rchange.Y2020,Flood_xplain$sealedchange)
corsealed=cor(Flood_xplain$Rchange.Y2020,Flood_xplain$sealedchange)

Flood_xplainLU=Flood_xplain[,c(71,86,87,88,89,90,91)]
Flood_xplainLU=Flood_xplain[,c(95,112,113,114,115,116,117)]


library(tidyverse)
library(ggpubr)



# Reshape the data from wide to long format
long_data <- Flood_xplainLU %>%
  gather(key = "variable", value = "value", -Y2020)
long_data=long_data[-which(abs(long_data$value)<=1e-3),]
# Create the base plot with facets
p <- long_data %>%
  ggplot(aes(x = value*100, y = Y2020)) +
  facet_wrap(~ variable, scales = "free") +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +  # Add a linear regression line
  # Calculate the Pearson correlation for each facet and add it as text
  geom_text(data = long_data %>% 
              group_by(variable) %>%
              summarize(correlation = cor(value*100, Y2020, use = "complete.obs")) %>%
              mutate(label = paste0("r = ", round(correlation, 2))),
            aes(label = label, x = Inf, y = Inf),
            hjust = 1.1, vjust = 1.1, inherit.aes = FALSE, check_overlap = TRUE)+
labs(title = "change in 10y flood (l/s/km2) attributed to LUC vs changes in land use fractions (%)")

# Print the final plot
print(p)
ggsave(paste0("D:/tilloal/Documents/LFRuns_utils/TAplots/LUdriversFLOOD_20201951_outlets.jpg"), p,width=30, height=20, units=c("cm"),dpi=300)


#Same for drought


#join flood and land use changes
Drought_xplain=inner_join(DroughtTrends,df_hybas07H,by=c("HYBAS_ID"))

#Drought_xplainLU=Drought_xplain[,c(71,86,87,88,89,90,91)]
Drought_xplainLU=Drought_xplain[,c(95,112,113,114,115,116,117)]

library(tidyverse)
library(ggpubr)



# Reshape the data from wide to long format
long_data <- Drought_xplainLU %>%
  gather(key = "variable", value = "value", -Y2020)
long_data=long_data[-which(abs(long_data$value)<=1e-3),]
# Create the base plot with facets
p <- long_data %>%
  ggplot(aes(x = value*100, y = Y2020)) +
  facet_wrap(~ variable, scales = "free") +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +  # Add a linear regression line
  # Calculate the Pearson correlation for each facet and add it as text
  geom_text(data = long_data %>% 
              group_by(variable) %>%
              summarize(correlation = cor(value*100, Y2020, use = "complete.obs")) %>%
              mutate(label = paste0("r = ", round(correlation, 2))),
            aes(label = label, x = Inf, y = Inf),
            hjust = 1.1, vjust = 1.1, inherit.aes = FALSE, check_overlap = TRUE)+
  labs(title = "change in 10y drought (l/s/km2) attributed to LUC vs changes in land use fractions (%)")

# Print the final plot
print(p)
ggsave(paste0("D:/tilloal/Documents/LFRuns_utils/TAplots/LUdriversDROUGHT_20201951_outlets.jpg"), p,width=30, height=20, units=c("cm"),dpi=300)





#Water demand
DroughtTrends=Output_dr_nfrost$TrendOutlets
DroughtTrends=DroughtTrends[which(DroughtTrends$driver=="Wateruse"),]

Drought_xplain=inner_join(DroughtTrends,df_hybas07H,by=c("HYBAS_ID"))

Drought_xplainWD=Drought_xplain[,c(71,92,93,94,95,96)]
Drought_xplainWD=Drought_xplain[,c(95,118,119,120,121,122)]



# Reshape the data from wide to long format
long_data <- Drought_xplainWD %>%
  gather(key = "variable", value = "value", -Y2020)
long_data=long_data[-which(abs(long_data$value)<=1e-3),]
# Create the base plot with facets
p <- long_data %>%
  ggplot(aes(x = value, y = Y2020)) +
  facet_wrap(~ variable, scales = "free") +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +  # Add a linear regression line
  # Calculate the Pearson correlation for each facet and add it as text
  geom_text(data = long_data %>% 
              group_by(variable) %>%
              summarize(correlation = cor(value, Y2020, use = "complete.obs")) %>%
              mutate(label = paste0("R = ", round(correlation, 2))),
            aes(label = label, x = Inf, y = Inf),
            hjust = 1.1, vjust = 1.1, inherit.aes = FALSE, check_overlap = TRUE)+
  labs(title = "change in 10y drought (% of 1951 RL) attributed to water demand vs changes in water demand (km3)")

# Print the final plot
print(p)
ggsave(paste0("D:/tilloal/Documents/LFRuns_utils/TAplots/WDdriversDROUGHT_20201951_outlets.jpg"), p,width=30, height=20, units=c("cm"),dpi=300)


#Reservoirs


