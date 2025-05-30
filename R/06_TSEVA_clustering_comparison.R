
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


source("D:/tilloal/Documents/LFRuns_utils/TrendAnalysis/R/functions_trends2.R")

# The objective of this scipt is to compare the quality of different spatial aggregations:

# HydroRegions
# Hybas levels 7 and lower
# NUTS3
# NUTS2


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

##2.1 Spatial data for catchments ----

### Hybas07 ----
Catchmentrivers7=read.csv(paste0(hydroDir,"/Catchments/from_hybas_eu_onlyid.csv"),encoding = "UTF-8", header = T, stringsAsFactors = F)
hybas07 <- read_sf(dsn = paste0(hydroDir,"/Catchments/hydrosheds/hybas_eu_lev07_v1c.shp"))
hybasf7=fortify(hybas07)
Catamere07=inner_join(hybasf7,Catchmentrivers7,by= "HYBAS_ID")
Catamere07$llcoord=paste(round(Catamere07$POINT_X,4),round(Catamere07$POINT_Y,4),sep=" ")
cst7=right_join(Catamere07,outf,by= c("llcoord"="latlong"))
GNF=cst7
st_geometry(GNF)=NULL
rm(Catamere07)
#cst7=st_transform(cst7,  crs=3035)
outlethybas07="outletsv8_hybas07_01min"
outhybas07=outletopen(hydroDir,outlethybas07)

#matching outlets with pixel Ids
outhybas07$latlong=paste(round(outhybas07$Var1,4),round(outhybas07$Var2,4),sep=" ")
mhy=match(outhybas07$latlong,outf$latlong)
outhybas07$outID=outf$outl2[mhy]

### European Biogeo regions ----

biogeo <- read_sf(dsn = paste0(hydroDir,"/eea_3035_biogeo-regions_2016/BiogeoRegions2016_wag84.shp"))
biogeof=fortify(biogeo)
st_geometry(biogeof)<-NULL
biogeoregions=raster( paste0(hydroDir,"/eea_3035_biogeo-regions_2016/Biogeo_rasterized_wsg84.tif"))
Gbiogeoregions=as.data.frame(biogeoregions,xy=T)
biogeomatch=inner_join(biogeof,Gbiogeoregions,by= c("PK_UID"="Biogeo_rasterized_wsg84"))
biogeomatch$latlong=paste(round(biogeomatch$x,4),round(biogeomatch$y,4),sep=" ")
biogeo_rivers=right_join(biogeomatch,outf, by="latlong")

### HydroRegions ----

GridHR=raster( paste0(hydroDir,"/HydroRegions_raster_WGS84.tif"))
GHR=as.data.frame(GridHR,xy=T)
GHR=GHR[which(!is.na(GHR[,3])),]
GHR$llcoord=paste(round(GHR$x,4),round(GHR$y,4),sep=" ")
GHR_riv=inner_join(GHR,outf,by= c("llcoord"="latlong"))
GHshpp <- read_sf(dsn ="Z:/ClimateRun4/nahaUsers/tilloal/HydroRegions/her_all_adjusted.shp")
HydroRsf=fortify(GHshpp)

### NUTS3 ----

NUTS3 <- read_sf(dsn = paste0(hydroDir,"/Countries/NUTS3/NUTS3_Extended_domain.shp"))
NUTS3$N3ID=c(1:length(NUTS3$NUTS_ID))
N2ID=unique(NUTS3$NUTS2_ID)
N2IDn=c(1:length(N2ID))
mati=match(NUTS3$NUTS2_ID,N2ID)
NUTS3$N2ID=N2IDn[mati]
#st_write(NUTS3, paste0(hydroDir,"/Countries/NUTS3/NUTS3_modified.shp"), driver = "ESRI Shapefile")

NUTS3 <- read_sf(dsn = paste0(hydroDir,"/Countries/NUTS3/NUTS3_modified.shp"))
GridNUTS3=raster( paste0(hydroDir,"/Countries/NUTS3/NUTS3_Raster3ID.tif"))
GN3=as.data.frame(GridNUTS3,xy=T)
GN3=GN3[which(!is.na(GN3[,3])),]
GN3$llcoord=paste(round(GN3$x,4),round(GN3$y,4),sep=" ")
GN3_riv=right_join(GN3,outf,by= c("llcoord"="latlong"))
# 
NUTS3 <- read_sf(dsn = paste0(hydroDir,"/Countries/NUTS3/NUTS3_modified.shp"))
GridNUTS2=raster( paste0(hydroDir,"/Countries/NUTS3/NUTS3_Raster2ID.tif"))
GN2=as.data.frame(GridNUTS2,xy=T)
GN2=GN2[which(!is.na(GN2[,3])),]
GN2$llcoord=paste(round(GN2$x,4),round(GN2$y,4),sep=" ")
GN2_riv=right_join(GN2,outf,by= c("llcoord"="latlong"))


NUTS1 <- read_sf(dsn = paste0(hydroDir,"/Countries/NUTS_RG_20M_2021_4326.shp"))
NUTS12 <- read_sf(dsn = paste0(hydroDir,"/Countries/NUTS_RG_20M_2024_4326.shp"))
NUTS1=NUTS1[which(NUTS1$LEVL_CODE==1),]
NUTS12=NUTS12[which(NUTS12$LEVL_CODE==1),]

NUTS1=rbind(NUTS1,NUTS12)
unik=unique(NUTS1$NUTS_ID)
muk=match(unik,NUTS1$NUTS_ID)
NUTS1=NUTS1[muk,]
NUTS1$NID=c(1:length(NUTS1$NUTS_ID))
# 
# st_write(NUTS1, paste0(hydroDir,"/Countries/NUTS1a_2021_4326.shp"), driver = "ESRI Shapefile")
GridNUTS1=raster( paste0(hydroDir,"/Countries/NUTS3/NUTS1_Raster1ID.tif"))
GN1=as.data.frame(GridNUTS1,xy=T)
GN1=GN1[which(!is.na(GN1[,3])),]
GN1$llcoord=paste(round(GN1$x,4),round(GN1$y,4),sep=" ")
GN1_riv=right_join(GN1,outf,by= c("llcoord"="latlong"))
# 
# GNF=right_join(GN3,GN2_riv,by="llcoord")
# 
GNUTS1sf=fortify(NUTS1) 
# 
# GNFx=GNF[which(is.na(GNF$NUTS3_Raster3ID)),]

Regio=GNUTS3sf
yrname=colnames(DataSfC)[12:81]
RegioP=GN3_riv
RegioName=RegioName
RegionAggregate<- function(Drivertrend, RegioP, Regio, RegioName) {
  
  # Data preparation
  
  # Aggregation choice
  data=Drivertrend
  DataC <- right_join(RegioP, data, by = c("outl2"))
  #DataC <- right_join(GNF, data, by = c("outl2" = "unikout"))
  
  if (RegioName=="H7"){
    DataC$IDR=DataC$HYBAS_ID
    Regio$IDR=Regio$HYBAS_ID
  }
  if (RegioName=="H6"){
    DataC$IDR=DataC$HYBAS_ID
    Regio$IDR=Regio$HYBAS_ID
  }
  if (RegioName=="HR"){
    DataC$IDR=DataC$HER
    Regio$IDR=Regio$CODEB
  }
  if (RegioName=="N3"){
    DataC$IDR=DataC$NUTS3_Raster3ID
    Regio$IDR=Regio$N3ID
  }
  if (RegioName=="N2"){
    DataC$IDR=DataC$NUTS3_Raster2ID
    Regio$IDR=Regio$N2ID
  }
  if (RegioName=="N1"){
    DataC$IDR=DataC$NUTS1_Raster1ID
    Regio$IDR=Regio$NID
  }
  
  DataC=data.frame(DataC)
  length(which(is.na(DataC$Var2)))
  # Matching Regio Ids
  HRM <- match(DataC$IDR, Regio$IDR)
  
  # Adding NUTS3 and NUTS2 IDs to the data
  DataC$Regio_id <- Regio$IDR[HRM]
  
  
  # Point aggregation by NUTS3
  # Here I aggregate relative changes, while at pixel level I show absolute changes
  pointagg <- aggregate(list(Rchange_rel = DataC$Y2020),
                        by = list(Region = DataC$Regio_id),
                        FUN = function(x) c(mean = mean(x, na.rm = TRUE), 
                                            dev = sd(x, na.rm = TRUE), 
                                            len = length(x), 
                                            med = median(x, na.rm = TRUE)))
  pointD <- do.call(data.frame, pointagg)
  
  
  # Return the final trendClim as the main output
  return(pointD)
}

###load UpArea -----
#load upstream area
main_path = 'D:/tilloal/Documents/06_Floodrivers/'
valid_path = paste0(main_path,'DataPaper/')
outletname="/GIS/upArea_European_01min.nc"
dir=valid_path
outf$idlalo=paste(outf$idlo, outf$idla, sep=" ")
UpArea=UpAopen(valid_path,outletname,outf)
head(UpArea)

#floods
load(file=paste0(hydroDir,"/TSEVA/output_plots/Flood_pixChange.Rdata"))
DataSflood=DataSave
#droughts
load(file=paste0(hydroDir,"/TSEVA/output_plots/Drought_pixChange.Rdata"))
DataSfdrought=DataSave

#1 Analysis for flood: looking at climate driven change and total changes

DataSfC=DataSflood$Total
#Aggregation at different levels


#first is biogeoregions
# Point aggregation by NUTS3
# Here I aggregate relative changes, while at pixel level I show absolute changes
pointagg <- aggregate(list(Rchange_rel = DataSfC$Y2020),
                      by = list(Region = DataSfC$Biogeo_id),
                      FUN = function(x) c(mean = mean(x, na.rm = TRUE), 
                                          dev = sd(x, na.rm = TRUE), 
                                          len = length(x), 
                                          med = median(x, na.rm = TRUE)))
pointD <- do.call(data.frame, pointagg)

k=length(pointD$Region)
bgvar=var(pointD$Rchange_rel.mean)
disBG=(pointD$Rchange_rel.len/sum(pointD$Rchange_rel.len)*pointD$Rchange_rel.dev^2)
mean(disHR)
k=length(pointD$Region)
n=length(DataSfC$x)


CH_indexBG=(bgvar/(k-1))/(sum(disBG)/(n-k))



Regio=HydroRsf
RegioP=GHR_riv
RegioName="HR"
HER_agg=RegionAggregate(DataSfC, RegioP, Regio, RegioName)

Regio=GNUTS3sf
RegioP=GN3_riv
RegioName="N3"
NUTS3_agg=RegionAggregate(DataSfC, RegioP, Regio, RegioName) 

Regio=GNUTS3sf
RegioP=GN2_riv
RegioName="N2"
NUTS2_agg=RegionAggregate(DataSfC, RegioP, Regio, RegioName) 

Regio=GNUTS1sf
RegioP=GN1_riv
RegioName="N1"
NUTS1_agg=RegionAggregate(DataSfC, RegioP, Regio, RegioName) 


Regio=hybasf7
RegioP=GNF
RegioName="H7"
Hybas7_agg=RegionAggregate(DataSfC, RegioP, Regio, RegioName) 


totalVar=var(DataC$Y2020,na.rm=T)
#dispersion measure
k=length(HER_agg$Region)
bcvar=var(HER_agg$Rchange_rel.mean*HER_agg$Rchange_rel.len/sum(HER_agg$Rchange_rel.len))
disHR=(HER_agg$Rchange_rel.len/sum(HER_agg$Rchange_rel.len)*HER_agg$Rchange_rel.dev^2)
mean(disHR)
k=length(HER_agg$Region)
n=length(DataC$x.x)


CH_index=(bcvar/(k-1))/(sum(disHR)/(n-k))


f_statistic <- ( sum(disHR)/ (n-k))

n2var=sd(NUTS2_agg$Rchange_rel.mean)^2
disN2=(NUTS2_agg$Rchange_rel.len/sum(NUTS2_agg$Rchange_rel.len)*NUTS2_agg$Rchange_rel.dev^2)
mean(disN2)
disN2=(NUTS2_agg$Rchange_rel.dev^2)
k=length(NUTS2_agg$Region)
mean(disN2)*k
n=length(DataC$x.x)
f_statistic2 <-  ( sum(disN2)/ (n-k))


n3var=sd(NUTS3_agg$Rchange_rel.mean,na.rm=T)^2
NUTS3_agg$Rchange_rel.dev[which(is.na(NUTS3_agg$Rchange_rel.dev))]=0
disN3=sum(NUTS3_agg$Rchange_rel.len/sum(NUTS3_agg$Rchange_rel.len)*NUTS3_agg$Rchange_rel.dev^2)
k=length(NUTS3_agg$Region)
n=length(DataC$x.x)
f_statistic3 <-  k*( disN3/ (n-k))

k=length(NUTS1_agg$Region)
n1var=var(NUTS1_agg$Rchange_rel.mean*NUTS1_agg$Rchange_rel.len/sum(NUTS1_agg$Rchange_rel.len))
#NUTS1_agg$Rchange_rel.dev[which(is.na(NUTS1_agg$Rchange_rel.dev))]=0
disN1=(NUTS1_agg$Rchange_rel.len/sum(NUTS1_agg$Rchange_rel.len)*NUTS1_agg$Rchange_rel.dev^2)
mean(disN1)
NUTS1_agg$variability=disN1

n=length(DataC$x.x)
f_statisticn1 <-  ( sum(disN1)/ (n-k))
CH_index1=(n1var/(k-1))/(sum(disN1)/(n-k))


h7var=sd(Hybas7_agg$Rchange_rel.mean,na.rm=T)^2
Hybas7_agg$Rchange_rel.dev[which(is.na(Hybas7_agg$Rchange_rel.dev))]=0
disH7=sum(Hybas7_agg$Rchange_rel.len/sum(Hybas7_agg$Rchange_rel.len)*Hybas7_agg$Rchange_rel.dev^2)
k=length(Hybas7_agg$Region)
n=length(DataC$x.x)
f_statistic4 <- k*( disH7/ (n-k))

#2 Analysis for drought: looking at climate driven change and total changes





#Standard error plot



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
catmap=cst7
rm(cst7)
basemap=w2

st_geometry(databiclim) <- NULL
databic2=inner_join(HER_agg,GHshpp,by=c("Region"="CODEB"))

NUTS1_tr=st_transform(NUTS1,crs=3035)
databic2=inner_join(NUTS1_agg,NUTS1_tr,by=c("Region"="NID"))

palet=c(hcl.colors(11, palette = "OrRd", alpha = NULL, rev = T, fixup = TRUE))

map <- ggplot(basemap) +
  geom_sf(fill="white")+
  geom_sf(data = databic2, mapping = aes(fill = variability, geometry=geometry), alpha=1) +
  geom_sf(fill=NA, color="gray42") +
  # bi_scale_fill(pal = "BlueOr", dim = 3, na.value=colNA ) +
  # scale_fill_manual(values = colorp) +
  scale_fill_gradientn(
    colors=palet, name="Intra-cluster variability", limits=c(0,10), oob = scales::squish)   +
  coord_sf(xlim = c(min(nco[,1]),max(nco[,1])), ylim = c(min(nco[,2]),max(nco[,2])))+
  #bi_scale_color(pal = "BlueOr", dim = 3, na.value=colNA ) +
  # scale_color_manual(values = colorp) +
  ggtitle("Changes in drought flows")+
  labs()+
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


map

map <- ggplot(basemap) +
  geom_sf(fill="white")+
  geom_sf(data = databic2, mapping = aes(fill = xfl,geometry=geometry ), alpha=0.7, color = "transparent", size = 0.01) +
  geom_sf(fill=NA, color="gray42") +
  # bi_scale_fill(pal = "BlueOr", dim = 3, na.value=colNA ) +
  # scale_fill_manual(values = colorp) +
  scale_fill_gradientn(
    colors=palet, name="95% CI of mean", limits=c(0,5), oob = scales::squish)   +
  coord_sf(xlim = c(min(nco[,1]),max(nco[,1])), ylim = c(min(nco[,2]),max(nco[,2])))+
  #bi_scale_color(pal = "BlueOr", dim = 3, na.value=colNA ) +
  # scale_color_manual(values = colorp) +
  labs()+
  ggtitle("Changes in flood flows")+
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


map


#I could compare how namy region have significant changes or not with different aggregations




