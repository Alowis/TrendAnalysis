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

library(patchwork)

lucmap  <- list()
luclass <- c("forest", "sealed", "irrigated", "other", "rice", "water")
panel_letters <- letters[1:length(luclass)]

# ── 1. Define break points and labels ─────────────────────────
breaks <- c(-Inf, -10, -5, -1, -0.1, 0, 0.1, 1, 5, 10, Inf)
labels <- c("< -10", "-10 - -5", "-5 - -1", "-1 - -0.1","-0.1 - 0",
            "0 - 0.1","0.1 - 1",  "1 - 5",   "5 - 10",  "> 10")

# ── 2. Define a diverging palette (8 categories) ───────────────
# Reds for negative, blues for positive, adjust to your palet style
cat_colors <- c(
  "#8B0000",   # < -10%      dark red
  "#D73027",   # -10 to -5%  red
  "#FC8D59",   # -5 to -1%   light red
  "#FEE010",   # -1 to 0%    pale yellow-red
  "#FFFCCC",   # -1 to 0%    pale yellow-red
  "#E0F3F1",   # 0 to 1%     pale blue
  "#CCFFFF",   # 0 to 1%     pale blue
  "#91BFDB",   # 1 to 5%     light blue
  "#4575B4",   # 5 to 10%    blue
  "#1A237E"    # > 10%       dark blue
)
# ── Build plots WITHOUT any legend ────────────────────────────
lucmap <- list()

# ── 2. Pre-cut ALL columns with the SAME factor levels ─────────
for (li in 1:length(luclass)) {
  GHshppH[[paste0("fill_cat_", li)]] <- factor(
    cut(as.numeric(df_GHshppH[, 7 + li]) * 100,
        breaks = breaks, labels = labels,
        include.lowest = TRUE, right = TRUE),
    levels = labels    # <-- force identical levels in every panel
  )
}
for (li in 1:length(luclass)) {
  lu    <- luclass[li]
  label <- panel_letters[li]
  
  GHshppH$fill_cat <- GHshppH[[paste0("fill_cat_", li)]]
  
  lucmap[[li]] <- ggplot(basemap) +
    geom_sf(fill = "gray95", color = "transparent", size = 0.5) +
    geom_sf(data = GHshppH, aes(fill = fill_cat, geometry = geometry),
            alpha = 1, color = "transparent") +
    geom_sf(fill = "transparent", color = "gray30", size = 0.5) +
    coord_sf(xlim = c(min(nco[, 1]), max(nco[, 1])),
             ylim = c(min(nco[, 2]), max(nco[, 2]))) +
    scale_fill_manual(
      values   = setNames(cat_colors, labels),
      name     = "Change (%)",
      drop     = FALSE,
      na.value = "grey80",
      limits   = labels
    ) +
    labs(x = NULL, y = NULL) +
    ggtitle(paste0("(",label, ") ", lu)) +
    theme(
      plot.title       = element_text(size = 16, face = "bold", hjust = 0),
      axis.text        = element_text(size = osize),
      panel.background = element_rect(fill = "aliceblue", colour = "grey1"),
      panel.border     = element_rect(linetype = "solid", fill = NA, colour = "black"),
      legend.position  = "none",    # <-- NO legend on any panel
      panel.grid.major = element_line(colour = "grey70"),
      panel.grid.minor = element_line(colour = "grey90"),
      legend.key       = element_rect(fill = "transparent", colour = "transparent"),
      legend.key.size  = unit(1, "cm")
    )
}

# ── Build a standalone legend plot ────────────────────────────
legend_plot <- ggplot(
  data.frame(cat = factor(labels, levels = labels)),
  aes(x = 1, y = cat, fill = cat)) +
  geom_tile() +
  scale_fill_manual(
    values   = setNames(cat_colors, labels),
    name     = "Change (%)",
    drop     = FALSE,
    limits   = labels
  ) +
  theme_void() +
  theme(
    legend.position  = "right",
    legend.title     = element_text(size = 20, face = "bold"),
    legend.text      = element_text(size = 16),
    legend.key.size  = unit(2, "cm")
  ) +
  guides(fill = guide_legend(ncol = 1))

# ── Extract legend only using cowplot ─────────────────────────
library(cowplot)
standalone_legend <- get_legend(legend_plot)

# ── Combine maps + legend ─────────────────────────────────────
maps_combined <- wrap_plots(lucmap, ncol = 3)

combined <- plot_grid(
  maps_combined,
  standalone_legend,
  ncol   = 2,
  rel_widths = c(1, 0.16)   # tweak 0.12 to make legend wider/narrower
)


ggsave("D:/tilloal/Documents/LFRuns_utils/TAplots/luchange_all_20201951_HR.jpg",
       combined, width = 60, height = 35, units = "cm", dpi = 300)


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



library(patchwork)

wdclass       <- c("ene", "dom", "liv", "ind", "total")
panel_letters <- letters[1:length(wdclass)]

# ── 1. Define breaks/labels/colors ONCE ───────────────────────
# In km3/year — adjust thresholds to your data range
breaks_wd <- c(-Inf, -0.5, -0.1, -0.01, 0, 0.01, 0.1, 0.5, Inf)
labels_wd <- c("< -0.5", "-0.5 to -0.1", "-0.1 to -0.01", "-0.01 to 0",
               "0 to 0.01", "0.01 to 0.1", "0.1 to 0.5", "> 0.5")

cat_colors_wd <- c(
  "#8B0000",   # < -0.5      dark red
  "#D73027",   # -0.5 to -0.1
  "#FC8D59",   # -0.1 to -0.01
  "#FEE090",   # -0.01 to 0
  "#E0F3F8",   # 0 to 0.01
  "#91BFDB",   # 0.01 to 0.1
  "#4575B4",   # 0.1 to 0.5
  "#1A237E"    # > 0.5       dark blue
)

# ── 2. Pre-cut ALL columns with identical factor levels ────────
for (li in 1:length(wdclass)) {
  HRwgs84h[[paste0("fill_cat_", li)]] <- factor(
    cut(as.numeric(df_HRwgs84h[, 15 + li]),
        breaks = breaks_wd, labels = labels_wd,
        include.lowest = TRUE, right = TRUE),
    levels = labels_wd    # force identical levels in every panel
  )
}
# ── 3. Build maps WITHOUT legend ──────────────────────────────
wdmap <- list()
tsize <- 12

for (li in 1:length(wdclass)) {
  wd    <- wdclass[li]
  label <- panel_letters[li]
  
  HRwgs84h$fill_cat <- HRwgs84h[[paste0("fill_cat_", li)]]
  
  wdmap[[li]] <- ggplot(basemap) +
    geom_sf(fill = "gray95", color = "transparent", size = 0.5) +
    geom_sf(data = HRwgs84h, aes(fill = fill_cat, geometry = geometry),
            alpha = 1, color = "transparent") +
    geom_sf(fill = "transparent", color = "gray30", size = 0.5) +
    coord_sf(xlim = c(min(nco[, 1]), max(nco[, 1])),
             ylim = c(min(nco[, 2]), max(nco[, 2]))) +
    scale_fill_manual(
      values   = setNames(cat_colors_wd, labels_wd),
      name     = "Change (km³/year)",
      drop     = FALSE,
      na.value = "grey80",
      limits   = labels_wd
    ) +
    labs(x = NULL, y = NULL) +
    ggtitle(paste0("(",label, ") ", wd)) +
    theme(
      plot.title       = element_text(size = 16, face = "bold", hjust = 0),
      axis.text        = element_text(size = osize),
      panel.background = element_rect(fill = "aliceblue", colour = "grey1"),
      panel.border     = element_rect(linetype = "solid", fill = NA, colour = "black"),
      legend.position  = "none",    # no legend on panels
      panel.grid.major = element_line(colour = "grey70"),
      panel.grid.minor = element_line(colour = "grey90"),
      legend.key.size  = unit(1, "cm")
    )
}

# ── 4. Standalone legend ──────────────────────────────────────
legend_plot_wd <- ggplot(
  data.frame(cat = factor(labels_wd, levels = labels_wd)),
  aes(x = 1, y = cat, fill = cat)) +
  geom_tile() +
  scale_fill_manual(
    values = setNames(cat_colors_wd, labels_wd),
    name   = "Change (km³/year)",
    drop   = FALSE,
    limits = labels_wd
  ) +
  theme_void() +
  theme(
    legend.position = "right",
    legend.title    = element_text(size = 20, face = "bold"),
    legend.text     = element_text(size = 16),
    legend.key.size = unit(2, "cm")
  ) +
  guides(fill = guide_legend(ncol = 1))

standalone_legend_wd <- cowplot::get_legend(legend_plot_wd)

# ── 5. Combine: 5 maps in 2 rows + legend ─────────────────────
maps_combined_wd <- wrap_plots(wdmap, ncol = 3)

combined_wd <- plot_grid(
  maps_combined_wd,
  standalone_legend_wd,
  ncol       = 2,
  rel_widths = c(1, 0.16)
)

ggsave("D:/tilloal/Documents/LFRuns_utils/TAplots/wdchange_all_20201951.jpg",
       combined_wd, width = 55, height = 30, units = "cm", dpi = 300)
