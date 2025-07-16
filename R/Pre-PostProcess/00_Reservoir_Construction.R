#Identification of pre and post 1951 reservoirs in EFAS-----
# Library calling --------------------------------------------------
library(rgdal)
library(raster)
library(rgdal)
library(raster)
library(ncdf4)
library(lubridate)
library(ggplot2)
library(rasterVis)
library(tidyverse)
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)
library(rgeos)

dataDir<-("D:/tilloal/Documents/LFRuns_utils/data")
# Load the CSV file into a data frame
data <- read.csv(paste0(dataDir,"/Catchments/from_hybas_eu_allATTR.csv"), header = TRUE)

# Convert the HYBAS_ID column to a factor
data$HYBAS_ID <- as.character(data$HYBAS_ID)

data$llcoord=paste(round(data$POINT_X,4),round(data$POINT_Y,4),sep =" ")
# Create an empty list to store the filtered data frames
dir=hydroDir

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
#Import reservoir locations from netcdf
resOpen=function(dir,outletname){
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
  outlets=as.vector(outlets)
  outll=expand.grid(londat,latdat)
  lonlatloop=expand.grid(c(1:llo),c(1:lla))
  outll$res=outlets
  outll$idlo=lonlatloop$Var1
  outll$idla=lonlatloop$Var2
  
  outll$idlalo=paste(outll$idlo,outll$idla,sep=" ")
  outfinal=outll[which(!is.na(outll$res)),]
  return (outfinal)
}

hydroDir<-("D:/tilloal/Documents/LFRuns_utils/data")
#comparison file between 2020 and 1951

res2020=resOpen(hydroDir,"/reservoirs/reservoirs_volumes_2020_Domain2.nc")
res2020$idla=2970-res2020$idla+1
res2020$idlalo=paste(res2020$idlo,res2020$idla,sep=" ")
res1951=resOpen(hydroDir,"/reservoirs/reservoirs_volumes_1951.nc")

max(res2020$res)
matres=na.omit(match(res1951$idlalo,res2020$idlalo))
res_old=res2020[matres,]
res_new=res2020[-matres,]

res_comp=left_join(res_old,res1951,by="idlalo")


### [Plot] - Figure S4 - Old and new reservoirs included in HERA -----

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

colorz = c("new"="orange","old" ='tomato4')
lab1=c("post-1951","pre-1951")
res_new$status="new"
res_old$status="old"

res_f=rbind(res_new,res_old)
pointout <- st_as_sf(res_f, coords = c("Var1", "Var2"), crs = 4326)
pointout <- st_transform(pointout, crs = 3035)

rnetwork=st_as_sf(outll, coords = c("Var1", "Var2"), crs = 4326)
rnetwork <- st_transform(rnetwork, crs = 3035)



damm<-ggplot(basemap) +
  geom_sf(fill="gray95",color="gray10",size=0.5)+
  geom_sf(data=pointout,aes(geometry=geometry,size=res,col=status),alpha=.9,stroke=0,shape=16)+
  geom_sf(data=rnetwork,aes(geometry=geometry),col="royalblue",size=0.1,alpha=.9,stroke=0,shape=15)+
  coord_sf(xlim = c(min(nco[,1]),max(nco[,1])), ylim = c(min(nco[,2]),max(nco[,2])))+
  scale_colour_manual(values = colorz, name="Dam construction", labels=lab1) +
  # scale_fill_manual(values = colorz, name="Largest change driver",labels=lab1) +
  scale_size(range = c(0.1, 5), trans="sqrt",name= expression(paste("Reservoir volume ", (m^3),sep = " ")),
             breaks=c(1e+5,1e+6,1e+7,1e+8,1e+9,1e+10))+
  #ggtitle(title)+
  labs(x="Longitude", y = "Latitude")+
  guides(colour = guide_legend(override.aes = list(size = 10)))+
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

ggsave(paste0("D:/tilloal/Documents/LFRuns_utils/TrendAnalysis/plots/DamKontraktion.jpg"), damm, width=20, height=20, units=c("cm"),dpi=1000) 
