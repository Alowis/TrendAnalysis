# Library calling --------------------------------------------------
setwd(dirname(rstudioapi::getSourceEditorContext()$path))
getwd()
source("functions_trends.R")

# Pre-loaded results -----------

#Set data directory
hydroDir<-("D:/tilloal/Documents/LFRuns_utils/ChangingHydroExtremes/data")

###create the outlets file outf ----
if (!exists("outf")){
  outf=c()
  for( Nsq in 1:88){
    print(Nsq)
    rspace= read.csv(paste0(hydroDir,"/subspace_efas.csv"))
    rspace=rspace[,-1]
    nrspace=rspace[Nsq,]
    #outletname="outletsv8_hybas07_01min"
    #outletname="outlets_hybas09_01min"
    outletname="GeoData/efas_rnet_100km_01min"
    
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

###load results from previous script ----

#### floods ----
load(file=paste0(hydroDir,"/Flood/outputs_flood_year_relxHR.Rdata"))
#### droughts ----
load(file=paste0(hydroDir,"/Drought/outputs_drought_nonfrost_relxHR.Rdata"))

#Trend at regional level
FloodTrends=Output_fl_year$TrendRegio
DroughtTrends=Output_dr_nonfrost$TrendRegio


#1. Bivariate results plot ----

Totrend=FloodTrends[which(FloodTrends$driver=="Total"),]
colnames(Totrend)[c(2:71)]=c(1951:2020)
trtest <- suppressWarnings(melt(Totrend, id.vars = "HydroR", variable.name = "variable", value.name = "value"))

# Convert variable to character and extract the year
trtest$variable <- as.character(trtest$variable)
trtest$yr <- as.numeric(trtest$variable)
trtest$value=as.numeric(trtest$value)
#Temporal aggregation over the domain
trtime_global <- aggregate(list(value = trtest$value),
                           by = list(yr = trtest$yr),
                           FUN = function(x) c(mean = mean(x, na.rm = TRUE),
                                               l = length(x),
                                               med= median(x, na.rm=T),
                                               ql = quantile(x, 0.25, na.rm = TRUE),
                                               qh = quantile(x, 0.75, na.rm = TRUE),
                                               w1 = quantile(x, 0.025, na.rm = TRUE),
                                               w2 = quantile(x, 0.975, na.rm = TRUE)))

# Convert to a data frame
FloodTh <- do.call(data.frame, trtime_global)

#extract results by decade
decades=c(1955,1965,1975,1985,1995,2005,2015)
valcol=which(!is.na(match(trtest$yr,decades)))
trtest=trtest[valcol,]
trtest$value=as.numeric(trtest$value)
trtest$decad=trtest$yr

# Aggregate only by decade (without location)
trtime_global <- aggregate(list(value = trtest$value),
                           by = list(yr = trtest$decad),
                           FUN = function(x) c(mean = mean(x, na.rm = TRUE),
                                               l = length(x),
                                               med= median(x, na.rm=T),
                                               ql = quantile(x, 0.25, na.rm = TRUE),
                                               qh = quantile(x, 0.75, na.rm = TRUE),
                                               w1 = quantile(x, 0.025, na.rm = TRUE),
                                               w2 = quantile(x, 0.975, na.rm = TRUE)))

# Convert to a data frame
FloodGlobal <- do.call(data.frame, trtime_global)

#Similar procedure for drought changes
Totrend=DroughtTrends[which(DroughtTrends$driver=="Total"),]
colnames(Totrend)[c(2:71)]=c(1951:2020)
trtest <- suppressWarnings(melt(Totrend, id.vars = "HydroR", variable.name = "variable", value.name = "value"))

# Convert variable to character and extract the year
trtest$variable <- as.character(trtest$variable)
trtest$yr <- as.numeric(trtest$variable)
trtest$value=as.numeric(trtest$value)

trtime_global <- aggregate(list(value = trtest$value),
                           by = list(yr = trtest$yr),
                           FUN = function(x) c(mean = mean(x, na.rm = TRUE),
                                               l = length(x),
                                               med= median(x, na.rm=T),
                                               ql = quantile(x, 0.25, na.rm = TRUE),
                                               qh = quantile(x, 0.75, na.rm = TRUE),
                                               w1 = quantile(x, 0.025, na.rm = TRUE),
                                               w2 = quantile(x, 0.975, na.rm = TRUE)))

# Convert to a data frame
DroughtTh <- do.call(data.frame, trtime_global)

#extract results by decade
decades=c(1955,1965,1975,1985,1995,2005,2015)
valcol=which(!is.na(match(trtest$yr,decades)))
trtest=trtest[valcol,]
trtest$value=as.numeric(trtest$value)
trtest$decad=trtest$yr

# Aggregate only by decade (without location)
trtime_global <- aggregate(list(value = trtest$value),
                           by = list(yr = trtest$decad),
                           FUN = function(x) c(mean = mean(x, na.rm = TRUE),
                                               l = length(x),
                                               med= median(x, na.rm=T),
                                               ql = quantile(x, 0.25, na.rm = TRUE),
                                               qh = quantile(x, 0.75, na.rm = TRUE),
                                               w1 = quantile(x, 0.025, na.rm = TRUE),
                                               w2 = quantile(x, 0.975, na.rm = TRUE)))

# Convert to a data frame
DroughtGlobal <- do.call(data.frame, trtime_global)


## [Plot] - Figure 1 - total changes in time for flood and drought hazards -----

FloodTh$haz="Flood"
DroughtTh$haz="Drought"
FDChanges=rbind(FloodTh,DroughtTh)
colorz = c("Flood" ='darkblue',"Drought"="darkorange")
colorn = c("Flood" ='darkblue',"Drought"="orange")
fac=1
xlabs=seq(1950,2020,5)
clabels=c("Drought","Flood")
nplot="Change (% 10yRL)"
br=c(seq(-100,-20,20),c(-10,-5,-2,0,2,5,10),seq(20,100,20))

ggplot() +
  geom_hline(yintercept = 0, color = "black", size = 1) +
  # IQR represented as rectangles
  geom_line(data=FDChanges,aes(x=yr, y=value.mean,color = (haz),group=(haz)),
                lwd=2) +
geom_ribbon(data = FDChanges, aes(x = yr, ymin = value.w1.2.5., ymax = value.w2.97.5., fill = haz, group = haz), 
            alpha = 0.2) +
  scale_color_manual(values = colorn, name = "Hazard", labels = clabels) +
  # Y-axis settings
  scale_y_continuous(name = nplot, breaks = br, trans=scales::modulus_trans(.4)) +
  
  # X-axis settings
  scale_x_continuous(breaks = xlabs, labels = xlabs, name = "Years",limits=c(1955,2020),
                     minor_breaks = seq(1955,2005,10), expand = c(.001,0.001)) +
  
  # Manual fill and color scales
  scale_fill_manual(values = colorn, name = "Hazard", labels = clabels) +

  # Theme customization
  theme(
    axis.title = element_text(size = 18, face = "bold"),
    title = element_text(size = 22, face = "bold"),
    axis.text = element_text(size = 16),
    axis.text.x = element_text(size = 16, face = "bold"),
    panel.background = element_rect(fill = "white", colour = "white"),
    panel.grid = element_blank(),
    panel.border = element_rect(linetype = "solid", fill = NA, colour = "black"),
    legend.title = element_text(size = 20, face = "bold"),
    legend.text = element_text(size = 16),
    axis.ticks.y = element_blank(),
    panel.grid.major.y = element_line(color = "lightgray",linetype = "dashed"),
    legend.position = "right",
    legend.key = element_rect(fill = "transparent", colour = "transparent"),
    legend.key.size = unit(0.8, "cm")
  ) +
  # Add title
  ggtitle("Europe")

ggsave(paste0("D:/tilloal/Documents/LFRuns_utils/TrendAnalysis/plots/temporalCh_clim_fldr23.jpg"), width=30, height=10, units=c("cm"),dpi=1000) 

##[Plot] - Supplement figure - Boxplot of changes in flood and drought hazards by decade -----

trtF=rbind(FloodGlobal,DroughtGlobal)
trtF$year=rep(seq(1950,2010,10),2)
trtF$hazard=c(rep("Flood",7),rep("Drought",7))

names(trtF)[c(2,4,5,6,7,8)]=c("changeC","med","cq1","cq2","w1","w2")

colorz = c("Flood" ='darkblue',"Drought"="darkorange")
colorn = c("Flood" ='darkblue',"Drought"="orange")

xlabs=seq(1950,2010,10)
clabels=c("Drought","Flood")

nplot="Mean change (% 10yRL)"
br=seq(-50,100,10)
br=seq(-50,300,10)
br=c(seq(-100,-10,10),seq(-5,5,5),seq(10,100,10))
ggplot() +
  # IQR represented as rectangles
  geom_linerange(data=trtF,aes(x=year, ymin=fac*(w1),ymax=fac*w2,color = (hazard),group=(hazard)),
                 position = position_dodge2(width = 9),lwd=1,alpha=0.8) +
  scale_color_manual(values = colorn, name = "Hazard", labels = clabels) +
  new_scale_color()+
  
  geom_rect(data=trtF, aes(xmin = year - 4.5, xmax = year + 4.5, 
                           ymin = cq1, ymax = cq2, fill = (hazard), group = (hazard)), 
            alpha = 0.5, position = position_dodge(width = 9)) +
  
  # Mean as points over the IQR
  geom_point(data=trtF, aes(x = year, y = changeC, color = (hazard), group = (hazard)),
             position = position_dodge(width = 9), size = 3) +
  
  # Median as a horizontal line within the IQR rectangle
  geom_rect(data=trtF, aes(xmin = year - 4.5, xmax = year + 4.5, 
                           ymin = med-1e-1, ymax = med+1e-1, fill = (hazard), group = (hazard)), 
            alpha = 1, position = position_dodge(width = 9)) +
  
  # Y-axis settings
  scale_y_continuous(name = nplot, breaks = br, trans=scales::modulus_trans(.6)) +
  
  # X-axis settings
  scale_x_continuous(breaks = xlabs, labels = xlabs, name = "Decades",
                     minor_breaks = seq(1955,2005,10), expand = c(.01,0.01)) +
  
  # Manual fill and color scales
  scale_fill_manual(values = colorn, name = "Hazard", labels = clabels) +
  scale_color_manual(values = colorz, name = "Hazard", labels = clabels) +

  # Theme customization
  theme(
    axis.title = element_text(size = 18, face = "bold"),
    title = element_text(size = 22, face = "bold"),
    axis.text = element_text(size = 16),
    axis.text.x = element_text(size = 16, face = "bold"),
    panel.background = element_rect(fill = "white", colour = "white"),
    panel.grid = element_blank(),
    panel.border = element_rect(linetype = "solid", fill = NA, colour = "black"),
    legend.title = element_text(size = 20, face = "bold"),
    legend.text = element_text(size = 16),
    axis.ticks.y = element_blank(),
    panel.grid.major.y = element_line(color = "lightgray",linetype = "dashed"),
    legend.position = "right",
    panel.grid.minor.x = element_line(colour = "grey23",linetype = "dashed"),
    legend.key = element_rect(fill = "transparent", colour = "transparent"),
    legend.key.size = unit(0.8, "cm")
  ) +
  
  # Add title
  ggtitle("Europe")

ggsave(paste0("D:/tilloal/Documents/LFRuns_utils/TrendAnalysis/plots/bxplot_fldr5.jpg"), width=30, height=20, units=c("cm"),dpi=1000) 



driver=unique(DroughtTrends$driver)

## load Biogeographic regions ----
biogeo <- read_sf(dsn = paste0(hydroDir,"/eea_3035_biogeo-regions_2016/BiogeoRegions2016_wag84.shp"))
biogeof=fortify(biogeo)
st_geometry(biogeof)<-NULL
biogeoregions=raster( paste0(hydroDir,"/eea_3035_biogeo-regions_2016/Biogeo_rasterized_wsg84.tif"))
Gbiogeoregions=as.data.frame(biogeoregions,xy=T)
biogeomatch=inner_join(biogeof,Gbiogeoregions,by= c("PK_UID"="Biogeo_rasterized_wsg84"))
biogeomatch$latlong=paste(round(biogeomatch$x,4),round(biogeomatch$y,4),sep=" ")
biogeo_rivers=right_join(biogeomatch,outf, by="latlong")

### load Hybas07 ----
Catchmentrivers7=read.csv(paste0(hydroDir,"/Catchments/from_hybas_eu_onlyid.csv"),encoding = "UTF-8", header = T, stringsAsFactors = F)
hybas07 <- read_sf(dsn = paste0(hydroDir,"/Catchments/hydrosheds/hybas_eu_lev07_v1c.shp"))
hybasf7=fortify(hybas07)
Catamere07=inner_join(hybasf7,Catchmentrivers7,by= "HYBAS_ID")
Catamere07$llcoord=paste(round(Catamere07$POINT_X,4),round(Catamere07$POINT_Y,4),sep=" ")
cst7=right_join(Catamere07,outf,by= c("llcoord"="latlong"))
GNF=cst7
st_geometry(GNF)=NULL
UnHY=unique(GNF$HYBAS_ID)

### load HydroRegions ----

GridHR=raster( paste0(hydroDir,"/HydroRegions_raster_WGS84.tif"))
GHR=as.data.frame(GridHR,xy=T)
GHR=GHR[which(!is.na(GHR[,3])),]
a1=(unique(GHR$HydroRegions_raster_WGS84))
GHR$llcoord=paste(round(GHR$x,4),round(GHR$y,4),sep=" ")
GHR_riv=inner_join(GHR,outf,by= c("llcoord"="latlong"))
a2=(unique(GHR_riv$HydroRegions_raster_WGS84))

ups=which(is.na(match(a1,a2)))
mierda=a1[ups]
GHshpp <- read_sf(dsn ="Z:/ClimateRun4/nahaUsers/tilloal/HydroRegions/her_all_adjusted_cor.shp")
HydroRsf=fortify(GHshpp)
length(unique(GHshpp$IRST_NAMEB))

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
biobase=st_transform(biogeo, crs=3035)
tsize=12
osize=12
Impdates=seq(1950,2020,by=10)
valuenames=paste0("Y",Impdates)
catmap=cst7
basemap=w2


#recomputing trend significance
trendplot=Output_fl_year$TrendRegio[which(Output_fl_year$TrendRegio$driver=="Clim"),]
median(trendplot$Rchange.Y2015)
trendpfas=t(trendplot[which(trendplot$HydroR==23),-c(1,72)])
plot(trendpfas)
trendpagg=trendplot[,c(1,71,72)]
trendFlood=calculateTrendSig(trendPlot=trendplot,pointagg=trendpagg)

length(which(trendFlood$Rchange.Y2020>0))
length(which(trendFlood$change<(-1)))

trendplot=Output_dr_nonfrost$TrendRegio[which(Output_dr_nonfrost$TrendRegio$driver=="Clim"),]
median(trendplot$Rchange.Y2015)
trendpagg=trendplot[,c(1,71,72)]
trendDrought=calculateTrendSig(trendPlot=trendplot,pointagg=trendpagg)
length(which(trendDrought$Rchange.Y2020<0))
length(which(trendDrought$change<(0)))

FloodTrends=Output_fl_year$TrendRegio
DroughtTrends=Output_dr_nonfrost$TrendRegio

#other analysis
DroughtTP=Output_dr_nonfrost$TrendPix
DroughtTP=DroughtTP[which(DroughtTP$driver=="Clim"),]
length(DroughtTP$Y2015[which(DroughtTP$Y2015>0)])/length(DroughtTP$Y2015[which(!is.na(DroughtTP$Y2015))])

FloodTP=Output_fl_year$TrendPix
FloodTP=FloodTP[which(FloodTP$driver=="Clim"),]
length(FloodTP$Y2015[which(FloodTP$Y2015>=15)])/length(FloodTP$Y2015[which(!is.na(FloodTP$Y2015))])

#1. TOTAL Trend ------------

##1.1 Bivariate categorization at catchment level ---------------
TotalFloodTrend=FloodTrends[which(FloodTrends$driver==driver[5]),c(2:71)]
TotalDroughtTrend=DroughtTrends[which(DroughtTrends$driver==driver[5]),c(2:71)]
TotalFloodTrend$HydroR=FloodTrends[which(FloodTrends$driver==driver[1]),1]
TotalDroughtTrend$HydroR=DroughtTrends[which(DroughtTrends$driver==driver[1]),1]

FlPoint=inner_join(GHshpp,TotalFloodTrend,by=c("CODEB"="HydroR"))
st_geometry(FlPoint)<-NULL


FlpointT=TotalFloodTrend
colnames(FlpointT)[c(1:70)]=c(1951:2020)
trtest <- suppressWarnings(melt(FlpointT, id.vars = "HydroR", variable.name = "variable", value.name = "value"))
# Convert variable to character and extract the year
trtest$variable <- as.character(trtest$variable)
trtest$variable<-as.numeric(trtest$variable)

#Plot trajectories of different regions
for (idr in FlpointT$HydroR){
  trp=trtest[which(trtest$HydroR==idr),]
  plot(trp$variable,trp$value,main=idr)
}
ggplot(trtest,aes(x=variable,y=value,color=HydroR,group=HydroR))+
  geom_line()



Flplot=inner_join(HydroRsf,FlPoint,by= c("Id"))
period=c(1951,2020)
haz="flood"

Flplot <- st_transform(Flplot, crs = 3035)
DrPoint=inner_join(GHshpp,TotalDroughtTrend,by=c("CODEB"="HydroR"))

st_geometry(DrPoint)<-NULL
Drplot=inner_join(HydroRsf,DrPoint,by= c("Id"))
period=c(1951,2020)
haz="flood"
Drplot <- st_transform(Drplot, crs = 3035)

Flplot$d2015=Drplot$Rchange.Y2015
Flplot$f2015=Flplot$Rchange.Y2015
FlplotTot=Flplot


#normalize flood change possibility
#here I do not normalize
a=b=1
databitot=FlplotTot
databitot$x=databitot$f2015/a
databitot$y=databitot$d2015/b

# breakers for category change: 0% and 5% absolute change
breaker1=0
breaker2=5

### flood -----------------
alterclass=data.frame(databitot$x)
alterclass$class=NA
alterclass$class[which(alterclass[,1]<=(-breaker2))]=1
alterclass$class[which(alterclass[,1]>(-breaker2) & alterclass[,1]<(breaker1))]=2
alterclass$class[which(alterclass[,1]>=(breaker1) & alterclass[,1]<(breaker2))]=3
alterclass$class[which(alterclass[,1]>=breaker2)]=4
c1=alterclass$class

### drought ----------------
alterclass=data.frame(databitot$y)
alterclass$class=NA
alterclass$class[which(alterclass[,1]<=(-breaker2))]=1
alterclass$class[which(alterclass[,1]>(-breaker2) & alterclass[,1]<(breaker1))]=2
alterclass$class[which(alterclass[,1]>=(breaker1) & alterclass[,1]<(breaker2))]=3
alterclass$class[which(alterclass[,1]>=breaker2)]=4
c2=alterclass$class

cx=paste(c2,c1,sep="-")
databitot$bi_class=cx
databitot$combined_category <-databitot$bi_class
##2.2 Bivariate categorization at pixel level ---------------

FloodTrendsP=Output_fl_year$TrendPix
DroughtTrendsP=Output_dr_nonfrost$TrendPix

TotalFloodTrendPix=FloodTrendsP[which(FloodTrendsP$driver==driver[5]),]
TotalDroughtTrendPix=DroughtTrendsP[which(DroughtTrendsP$driver==driver[5]),]

ClimFloodTrendPix=FloodTrendsP[which(FloodTrendsP$driver==driver[1]),]
ClimDroughtTrendPix=DroughtTrendsP[which(DroughtTrendsP$driver==driver[1]),]

LuFloodTrendPix=FloodTrendsP[which(FloodTrendsP$driver==driver[2]),]
LuDroughtTrendPix=DroughtTrendsP[which(DroughtTrendsP$driver==driver[2]),]

ResFloodTrendPix=FloodTrendsP[which(FloodTrendsP$driver==driver[3]),]
ResDroughtTrendPix=DroughtTrendsP[which(DroughtTrendsP$driver==driver[3]),]

WuFloodTrendPix=FloodTrendsP[which(FloodTrendsP$driver==driver[4]),]
WuDroughtTrendPix=DroughtTrendsP[which(DroughtTrendsP$driver==driver[4]),]

# Remove NA rows 
TotalFloodTrendPix = TotalFloodTrendPix[which(!is.na(TotalFloodTrendPix$Var1)),]
TotalDroughtTrendPix = TotalDroughtTrendPix[which(!is.na(TotalDroughtTrendPix$Var1)),]

ClimFloodTrendPix = ClimFloodTrendPix[which(!is.na(ClimFloodTrendPix$Var1)),]
ClimDroughtTrendPix = ClimDroughtTrendPix[which(!is.na(ClimDroughtTrendPix$Var1)),]

LuFloodTrendPix = LuFloodTrendPix[which(!is.na(LuFloodTrendPix$Var1)),]
LuDroughtTrendPix = LuDroughtTrendPix[which(!is.na(LuDroughtTrendPix$Var1)),]

ResFloodTrendPix = ResFloodTrendPix[which(!is.na(ResFloodTrendPix$Var1)),]
ResDroughtTrendPix = ResDroughtTrendPix[which(!is.na(ResDroughtTrendPix$Var1)),]

WuFloodTrendPix = WuFloodTrendPix[which(!is.na(WuFloodTrendPix$Var1)),]
WuDroughtTrendPix = WuDroughtTrendPix[which(!is.na(WuDroughtTrendPix$Var1)),]


points <- st_as_sf(TotalFloodTrendPix, coords = c("Var1", "Var2"), crs = 4326)
points <- st_transform(points, crs = 3035)
Flpixplot=points
Flpixplot$d2015=TotalDroughtTrendPix$Y2015

a=b=1
databipi=Flpixplot
databipi$x=databipi$Y2015/a
databipi$y=databipi$d2015/b

breaker1=0
breaker2=5

### floods --------------------
alterclass=data.frame(databipi$x)
alterclass$class=NA
alterclass$class[which(alterclass[,1]<=(-breaker2))]=1
alterclass$class[which(alterclass[,1]>(-breaker2) & alterclass[,1]<(breaker1))]=2
alterclass$class[which(alterclass[,1]>=(breaker1) & alterclass[,1]<(breaker2))]=3
alterclass$class[which(alterclass[,1]>=breaker2)]=4
c1=alterclass$class

### droughts --------------------
alterclass=data.frame(databipi$y)
alterclass$class=NA
alterclass$class[which(alterclass[,1]<=(-breaker2))]=1
alterclass$class[which(alterclass[,1]>(-breaker2) & alterclass[,1]<(breaker1))]=2
alterclass$class[which(alterclass[,1]>=(breaker1) & alterclass[,1]<(breaker2))]=3
alterclass$class[which(alterclass[,1]>=breaker2)]=4
c2=alterclass$class

cx=paste(c2,c1,sep="-")
databipi$bi_class=cx

##1.3 [Plot] - Figure 4 - Bivariate hydrological extremes change ---------

colors <- c(
  "1-1" = "#dd6a29",  # high x, low y
  "2-1" = "#d9926a",  # medium-high x, low y
  "3-1" = "#d6b3a0",  # medium-low x, low y
  "4-1" = "#d3d3d3",  # low x, low y

  "1-2" = "#a36229",  # high x, medium-low y
  "2-2" = "#a08769",  # medium-high x, medium-low y
  "3-2" = "#9ea69f",  # medium-low x, medium-low y
  "4-2" = "#9cc4d2",  # low x, medium-low y
  
  "1-3" = "#635929",  # high x, medium-high y
  "2-3" = "#617b69",  # medium-high x, medium-high y
  "3-3" = "#60979f",  # medium-low x, medium-high y
  "4-3" = "#5fb2d1",  # low x, medium-high y
  
  "1-4" = "#174f28",  # high x, high y
  "2-4" = "#166d68",  # medium-high x, high y
  "3-4" = "#16869e",  # medium-low x, high y
  "4-4" = "#169dd0"  # low x, high y

)

colorp <- c(
  "1-1" = "#dd6a40",  # high x, low y
  "2-1" = "#d9926a",  # medium-high x, low y
  "3-1" = "#DEB887",  # medium-low x, low y
  "4-1" = "#FFD39B",  # low x, low y
  
  "1-2" = "#a36229",  # high x, medium-low y
  "2-2" = "#999999",  # medium-high x, medium-low y
  "3-2" = "#999999",  # medium-low x, medium-low y
  "4-2" = "#9cc4d2",  # low x, medium-low y
  
  "1-3" = "#635929",  # high x, medium-high y
  "2-3" = "#999999",  # medium-high x, medium-high y
  "3-3" = "#999999",  # medium-low x, medium-high y
  "4-3" = "#5fb2d1",  # low x, medium-high y
  
  "1-4" = "#174f28",  # high x, high y
  "2-4" = "#166d68",  # medium-high x, high y
  "3-4" = "#16869e",  # medium-low x, high y
  "4-4" = "#169dd0"  # low x, high y
  
)

loscolors=c("Accelerating" = "#174f28","Drying" = "#dd6a29","Stable"="gray60","Wetting" = "#169dd0","Decelerating" = "burlywood")

bi_pal(pal = colors, dim = 4)

#removal of pixel with NA direction
databipi$bi_class[which(databipi$bi_class=="NA-NA")]=NA
databipi$bi_class[which(is.na(databipi$Y2015))]=NA
databipi$bi_class[which(is.na(databipi$d2015))]=NA

map <- ggplot(basemap) +
  geom_sf(fill="white")+
  geom_sf(data = databitot, mapping = aes(fill = combined_category), alpha=0.7, color = "transparent", size = 0.01,show.legend = F) +
  geom_sf(data = databipi, mapping = aes(col = bi_class,geometry=geometry,size=upa), alpha=1,stroke=0,shape=15, show.legend = FALSE) +
  geom_sf(fill=NA, color="gray42") +
  scale_fill_manual(values = colorp) +
  coord_sf(xlim = c(min(nco[,1]),max(nco[,1])), ylim = c(min(nco[,2]),max(nco[,2])))+
  scale_size(range = c(0.08, 0.4), trans="sqrt",name= expression(paste("Upstream area ", (km^2),
                                                                       sep = " ")),
             breaks=c(101,1000,10000,100000,500000), labels=c("100","1000", "10 000", "100 000", "500 000"),
             guide = "none")+
  scale_color_manual(values = colorp,na.value=NA) +
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



legend <- bi_legend(pal = colorp,
                    dim = 4,
                    xlab = "  +  Drought intensity  -  ",
                    ylab = "  -  Flood intensity  +  ",
                    size = 16,
                    arrows = FALSE)

pl=ggarrange(map, legend,
             labels = c("Map", "Key"),
             ncol = 2, nrow = 1,widths = c(2,1), heights=c(1,1), vjust=-1)


ggsave(paste0("D:/tilloal/Documents/LFRuns_utils/TrendAnalysis/plots/totchange_bvPIX_HR_n10.jpg"), pl, width=20, height=20, units=c("cm"),dpi=500)


#1.4 Aggregation by Biogeoregion ----

#Color scale
colord <- c(
  "1-1" = 12,  # high x, low y
  "2-1" = 11,  # medium-high x, low y
  "3-1" = 10,  # medium-low x, low y
  "4-1" = 9,  # low x, low y
  
  "1-2" = 13,  # high x, medium-low y
  "2-2" = 4,  # medium-high x, medium-low y
  "3-2" = 3,  # medium-low x, medium-low y
  "4-2" = 8,  # low x, medium-low y
  
  "1-3" = 14,  # high x, medium-high y
  "2-3" = 1,  # medium-high x, medium-high y
  "3-3" = 2,  # medium-low x, medium-high y
  "4-3" = 7,  # low x, medium-high y
  
  "1-4" = 15,  # high x, high y
  "2-4" = 16,  # medium-high x, high y
  "3-4" = 5,  # medium-low x, high y
  "4-4" = 6  # low x, high y
  
)

zob=unique(databipi$bi_class)
zob=zob[-which(is.na(zob))]

databipi2=databipi

databipi2$trcat="Wetting"
databipi2$trcat[which(databipi2$bi_class=="2-2" | databipi2$bi_class=="3-2" |
                    databipi2$bi_class=="2-3" | databipi2$bi_class=="3-3")]="Stable"
databipi2$trcat[which(databipi2$bi_class=="1-1" | databipi2$bi_class=="1-2" |
                    databipi2$bi_class=="2-1")]="Drying"
databipi2$trcat[which(databipi2$bi_class=="1-4" | databipi2$bi_class=="2-4" |
                    databipi2$bi_class=="1-3")]="Accelerating"
databipi2$trcat[which(databipi2$bi_class=="4-1" | databipi2$bi_class=="4-2" |
                    databipi2$bi_class=="3-1")]="Decelerating"

databipi2$Biogeo_id[which(databipi2$Biogeo_id=="Pannonian")]="Continental"

if (length(which(is.na(databipi2$bi_class)))>0){
  databipi3=databipi2[-which(is.na(databipi2$bi_class)),]
}else{
  databipi3=databipi2
}


#Aggregation over each change trajectory subclasses
EuAgg = aggregate(list(val=databipi2$upa),
                 by = list(traj=databipi2$bi_class),
                 FUN = function(x) c(len=length(x)))
EuAgg$trcat="Wetting"
EuAgg$trcat[which(EuAgg$traj=="2-2" | EuAgg$traj=="3-2" |
                   EuAgg$traj=="2-3" | EuAgg$traj=="3-3")]="Stable"
EuAgg$trcat[which(EuAgg$traj=="1-1" | EuAgg$traj=="1-2" |
                   EuAgg$traj=="2-1")]="Drying"
EuAgg$trcat[which(EuAgg$traj=="1-4" | EuAgg$traj=="2-4" |
                   EuAgg$traj=="1-3")]="Accelerating"
EuAgg$trcat[which(EuAgg$traj=="4-1" | EuAgg$traj=="4-2" |
                   EuAgg$traj=="3-1")]="Decelerating"

#Aggregation over each change trajectory
EuAgg2 = aggregate(list(val=EuAgg$val),
                  by = list(traj=EuAgg$trcat),
                  FUN = function(x) c(sum=sum(x)))
SumEuPix=sum(EuAgg2$val)
EuAgg2$rel.val=EuAgg2$val/SumEuPix*100


RegioSize = aggregate(list(val=databipi3$upa),
                      by = list(region=databipi3$Biogeo_id),
                      FUN = function(x) c(len=length(x)))

sum(RegioSize$val)

#Aggregate by region and trajectory
magg = aggregate(list(val=databipi2$upa),
                 by = list(reg=databipi2$Biogeo_id ,traj=databipi2$bi_class),
                 FUN = function(x) c(len=length(x)))
magg <- do.call(data.frame, magg)

#Remove peripheral regions
magg=magg[-which(magg$reg=="Steppic" | magg$reg=="BlackSea" | magg$reg=="Arctic"),]
magg$trcat="Wetting"
magg$trcat[which(magg$traj=="2-2" | magg$traj=="3-2" |
                 magg$traj=="2-3" | magg$traj=="3-3")]="Stable"
magg$trcat[which(magg$traj=="1-1" | magg$traj=="1-2" |
                 magg$traj=="2-1")]="Drying"
magg$trcat[which(magg$traj=="1-4" | magg$traj=="2-4" |
                   magg$traj=="1-3")]="Accelerating"
magg$trcat[which(magg$traj=="4-1" | magg$traj=="4-2" |
                   magg$traj=="3-1")]="Decelerating"

#Plot preparation
magg$order=1
magg$order[which(magg$trcat=="Wetting")]=2
magg$order[which(magg$trcat=="Decelerating")]=3
magg$order[which(magg$trcat=="Drying")]=4
magg$order[which(magg$trcat=="Accelerating")]=5

matm=match(magg$traj,names(colord))
magg$ordf=colord[matm]
magg <- magg %>% 
  mutate(traj = reorder(traj, ordf, FUN = mean))


magg2 = aggregate(list(val=magg$val),
                   by = list(traj=magg$trcat,regio=magg$reg),
                   FUN = function(x) c(sum=sum(x)))

rsize=(match(magg2$reg,RegioSize$region))
magg2$rsize=NA
magg2$rsize=RegioSize$val[rsize]
magg2$ratio=magg2$val/magg2$rsize*100
sum(magg$val[which(magg$reg=="Mediterranean")])


### [Plot] - Figure 4 - Stacked barplot of change trajectory by Biogeoregion ----
loscolors=c("Accelerating" = "#174f28","Drying" = "#dd6a29","Stable"="gray60","Wetting" = "#169dd0","Decelerating" = "burlywood")

tsize=15
osize=18
ggplot(magg, aes(x = reg, y = val, fill = traj)) +
  geom_bar(stat = "identity", position = "fill") +
  scale_fill_manual(values = colorp) +
  scale_x_discrete(name="Biogeoregions")+
  scale_y_continuous(label=c(0,25,50,75,100),name="River pixels (%)")+
  theme(axis.title=element_text(size=tsize),
        axis.text=element_text(size=tsize),
        panel.background = element_rect(fill = "white", colour = "grey1"),
        panel.border = element_rect(linetype = "solid", fill = NA, colour="black"),
        legend.title = element_text(size=tsize),
        legend.text = element_text(size=osize),
        legend.position = "none",
        legend.key = element_rect(fill = "transparent", colour = "transparent"),
        legend.key.size = unit(.8, "cm"))


ggsave(paste0("D:/tilloal/Documents/LFRuns_utils/TrendAnalysis/plots/Regional_trajectories5.jpg"), width=20, height=40, units=c("cm"),dpi=800)


#2. CLIMATE TREND ------------------

##2.1 Bivariate categorization at catchment level ---------------
ClimateFloodTrend=FloodTrends[which(FloodTrends$driver==driver[1]),]
ClimateDroughtTrend=DroughtTrends[which(DroughtTrends$driver==driver[1]),]
ClimateFloodTrend$HydroR=FloodTrends[which(FloodTrends$driver==driver[1]),1]
ClimateDroughtTrend$HydroR=DroughtTrends[which(DroughtTrends$driver==driver[1]),1]


FlPoint=full_join(GHshpp,ClimateFloodTrend,by=c("CODEB"="HydroR"))
st_geometry(FlPoint)<-NULL
Flplot=inner_join(HydroRsf,FlPoint,by= c("Id"))
period=c(1951,2020)
haz="flood"

Flplot <- st_transform(Flplot, crs = 3035)

DrPoint=full_join(GHshpp,ClimateDroughtTrend,by=c("CODEB"="HydroR"))
st_geometry(DrPoint)<-NULL
period=c(1951,2020)
haz="flood"
Drplot=inner_join(HydroRsf,DrPoint,by= c("Id"))

mv=match(Drplot$Id,Flplot$Id)
Flplot$d2015=Drplot$Rchange.Y2015
Flplot$f2015=Flplot$Rchange.Y2015

colNA="transparent"
tsize=16
osize=12
legend="Change in Q (%)"
palet=c(hcl.colors(11, palette = "RdYlBu", alpha = NULL, rev = F, fixup = TRUE))
### [SPlot] - SUpplementary plot - Univariate plot for verification ----
uplot=F
if (uplot==T){
  br=seq(-100,100,2)
  labels=br
  limi=c(-20,20)
  ggplot(basemap) +
    geom_sf(fill="white")+
    geom_sf(data = Flplot, mapping = aes(fill = d2015), alpha=0.9, color = "transparent", size = 0.01, show.legend = TRUE) +
    geom_sf(fill=NA, color="gray42") +
    coord_sf(xlim = c(min(nco[,1]),max(nco[,1])), ylim = c(min(nco[,2]),max(nco[,2])))+
    scale_fill_gradientn(
      colors=palet,
      breaks=br,limits=limi,
      oob = scales::squish,na.value=colNA, name=legend)  +
    labs()+
    theme(axis.title=element_text(size=tsize),
          panel.background = element_rect(fill = "aliceblue", colour = "grey1"),
          panel.border = element_rect(linetype = "solid", fill = NA, colour="black"),
          legend.title = element_text(size=tsize),
          legend.text = element_text(size=osize),
          legend.position = "bottom",
          panel.grid.major = element_line(colour = "grey70"),
          panel.grid.minor = element_line(colour = "grey90"),
          legend.key = element_rect(colour = "transparent"),
          legend.key.size = unit(.8, "cm"))
}

#normalize change
FlplotClim=Flplot
a=b=1
databiclim=FlplotClim
databiclim$x=databiclim$f2015/a
databiclim$y=databiclim$d2015/b

breaker1=0
breaker2=5

#### floods --------------------
alterclass=data.frame(databiclim$x)
alterclass$class=NA
alterclass$class[which(alterclass[,1]<=(-breaker2))]=1
alterclass$class[which(alterclass[,1]>(-breaker2) & alterclass[,1]<(breaker1))]=2
alterclass$class[which(alterclass[,1]>=(breaker1) & alterclass[,1]<(breaker2))]=3
alterclass$class[which(alterclass[,1]>=breaker2)]=4
c1=alterclass$class

#### droughts --------------------
alterclass=data.frame(databiclim$y)
alterclass$class=NA
alterclass$class[which(alterclass[,1]<=(-breaker2))]=1
alterclass$class[which(alterclass[,1]>(-breaker2) & alterclass[,1]<(breaker1))]=2
alterclass$class[which(alterclass[,1]>=(breaker1) & alterclass[,1]<(breaker2))]=3
alterclass$class[which(alterclass[,1]>=breaker2)]=4
c2=alterclass$class

cx=paste(c2,c1,sep="-")
databiclim$bi_class=cx
# Combine the main category and subcategory to make a label
databiclim$combined_category <-databiclim$bi_class

##1.2 Bivariate categorization at pixel level ---------------
a=sd(ClimFloodTrendPix$Y2015,na.rm=T)
b=sd(ClimDroughtTrendPix$Y2015,na.rm=T)
a=b=1
databipic=Flpixplot
databipic$x=ClimFloodTrendPix$Y2015/a
databipic$y=ClimDroughtTrendPix$Y2015/b

sd(databipic$x,na.rm=T)
breaker1=0
breaker2=5

alterclass=data.frame(databipic$x)
alterclass$class=NA
alterclass$class[which(alterclass[,1]<=(-breaker2))]=1
alterclass$class[which(alterclass[,1]>(-breaker2) & alterclass[,1]<(breaker1))]=2
alterclass$class[which(alterclass[,1]>=(breaker1) & alterclass[,1]<(breaker2))]=3
alterclass$class[which(alterclass[,1]>=breaker2)]=4

c1=alterclass$class
alterclass=data.frame(databipic$y)
alterclass$class=NA
alterclass$class[which(alterclass[,1]<=(-breaker2))]=1
alterclass$class[which(alterclass[,1]>(-breaker2) & alterclass[,1]<(breaker1))]=2
alterclass$class[which(alterclass[,1]>=(breaker1) & alterclass[,1]<(breaker2))]=3
alterclass$class[which(alterclass[,1]>=breaker2)]=4

c2=alterclass$class

cx=paste(c2,c1,sep="-")

databipic$bi_class=cx

mbicli=data.frame(x=mean(databiclim$d2015,na.rm=T),y=mean(databiclim$f2015,na.rm=T),xsd=sd(databiclim$d2015,na.rm=T),ysd=sd(databiclim$f2015,na.rm=T))

mbicli=data.frame(x=mean(databiclim$d2015,na.rm=T),y=mean(databiclim$f2015,na.rm=T),
                 xq1=quantile(databiclim$d2015,0.05,na.rm=T),yq1=quantile(databiclim$f2015,0.05,na.rm=T),
                 xq2=quantile(databiclim$d2015,0.95,na.rm=T),yq2=quantile(databiclim$f2015,0.95,na.rm=T))

map <- ggplot(basemap) +
  geom_sf(fill="white")+
  geom_sf(data = databiclim, mapping = aes(fill = combined_category), alpha=0.7, color = "transparent", size = 0.01,show.legend = F) +
  geom_sf(data = databipic, mapping = aes(col = bi_class,geometry=geometry,size=upa), alpha=1,stroke=0,shape=15, show.legend = FALSE) +
  geom_sf(fill=NA, color="gray42") +
  scale_fill_manual(values = colorp) +
  coord_sf(xlim = c(min(nco[,1]),max(nco[,1])), ylim = c(min(nco[,2]),max(nco[,2])))+
  scale_size(range = c(0.08, 0.4), trans="sqrt",name= expression(paste("Upstream area ", (km^2),
                                                                       sep = " ")),
             breaks=c(101,1000,10000,100000,500000), labels=c("100","1000", "10 000", "100 000", "500 000"),
             guide = "none")+
  scale_color_manual(values = colorp, na.value=NA) +
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



legend <- bi_legend(pal = colorp,
                    dim = 4,
                    xlab = "  -  Drought flows  +  ",
                    ylab = "  -  Flood flows  +  ",
                    size = 16,
                    arrows = FALSE)

pl=ggarrange(map, legend,
             labels = c("Map", "Key"),
             ncol = 2, nrow = 1,widths = c(2,1), heights=c(1,1), vjust=-1)



ggsave(paste0("D:/tilloal/Documents/LFRuns_utils/TrendAnalysis/plots/climchange_HR_9.jpg"), pl, width=20, height=20, units=c("cm"),dpi=800)


#3. LAND USE TREND --------------

##3.1 Bivariate categorization at catchment level ---------------

LanduseFloodTrend=FloodTrends[which(FloodTrends$driver==driver[2]),]
LanduseDroughtTrend=DroughtTrends[which(DroughtTrends$driver==driver[2]),]
LanduseFloodTrend$HydroR=FloodTrends[which(FloodTrends$driver==driver[1]),1]
LanduseDroughtTrend$HydroR=DroughtTrends[which(DroughtTrends$driver==driver[1]),1]

LanduseFloodTrend$f2015=LanduseFloodTrend$Rchange.Y2015
LanduseDroughtTrend$f2015=LanduseDroughtTrend$Rchange.Y2015

FlPoint=full_join(GHshpp,LanduseFloodTrend,by=c("CODEB"="HydroR"))
st_geometry(FlPoint)<-NULL
Flplot=inner_join(HydroRsf,FlPoint,by= c("Id"))
period=c(1951,2020)
haz="flood"

Flplot <- st_transform(Flplot, crs = 3035)
DrPoint=full_join(GHshpp,LanduseDroughtTrend,by=c("CODEB"="HydroR"))
st_geometry(DrPoint)<-NULL
Drplot=inner_join(HydroRsf,DrPoint,by= c("Id"))
period=c(1951,2020)
haz="flood"
Drplot <- st_transform(Drplot, crs = 3035)
Flplot$d2015=Drplot$Rchange.Y2015
Flplot$f2015=Flplot$Rchange.Y2015

FlplotLu=Flplot
a=b=1
databilu=FlplotLu
databilu$x=databilu$f2015/a
databilu$y=databilu$d2015/b

breaker1=0
breaker2=5

#### floods --------------------
alterclass=data.frame(databilu$x)
alterclass$class=NA
alterclass$class[which(alterclass[,1]<=(-breaker2))]=1
alterclass$class[which(alterclass[,1]>(-breaker2) & alterclass[,1]<(breaker1))]=2
alterclass$class[which(alterclass[,1]>=(breaker1) & alterclass[,1]<(breaker2))]=3
alterclass$class[which(alterclass[,1]>=breaker2)]=4
c1=alterclass$class

#### droughts --------------------
alterclass=data.frame(databilu$y)
alterclass$class=NA
alterclass$class[which(alterclass[,1]<=(-breaker2))]=1
alterclass$class[which(alterclass[,1]>(-breaker2) & alterclass[,1]<(breaker1))]=2
alterclass$class[which(alterclass[,1]>=(breaker1) & alterclass[,1]<(breaker2))]=3
alterclass$class[which(alterclass[,1]>=breaker2)]=4
c2=alterclass$class

cx=paste(c2,c1,sep="-")
databilu$bi_class=cx
# Combine the main category and subcategory to make a label
databilu$combined_category <-databilu$bi_class

##3.2 Bivariate categorization at pixel level ---------------

luflood=LuFloodTrendPix$Y2015
ludrought=LuDroughtTrendPix$Y2015
a=b=1
databipilu=Flpixplot
databipilu$x=luflood/a
databipilu$y=ludrought/b


sd(databipilu$x,na.rm=T)

alterclass=data.frame(databipilu$x)
alterclass$class=NA
alterclass$class[which(alterclass[,1]<=(-breaker2))]=1
alterclass$class[which(alterclass[,1]>(-breaker2) & alterclass[,1]<(breaker1))]=2
alterclass$class[which(alterclass[,1]>=(breaker1) & alterclass[,1]<(breaker2))]=3
alterclass$class[which(alterclass[,1]>=breaker2)]=4

c1=alterclass$class
alterclass=data.frame(databipilu$y)
alterclass$class=NA
alterclass$class[which(alterclass[,1]<=(-breaker2))]=1
alterclass$class[which(alterclass[,1]>(-breaker2) & alterclass[,1]<(breaker1))]=2
alterclass$class[which(alterclass[,1]>=(breaker1) & alterclass[,1]<(breaker2))]=3
alterclass$class[which(alterclass[,1]>=breaker2)]=4

c2=alterclass$class

cx=paste(c2,c1,sep="-")

databipilu$bi_class=cx

mbilu=data.frame(x=mean(databilu$d2015,na.rm=T),y=mean(databilu$f2015,na.rm=T),xsd=sd(databilu$d2015,na.rm=T),ysd=sd(databilu$f2015,na.rm=T))


mbilu=data.frame(x=mean(databilu$d2015,na.rm=T),y=mean(databilu$f2015,na.rm=T),
                  xq1=quantile(databilu$d2015,0.05,na.rm=T),yq1=quantile(databilu$f2015,0.05,na.rm=T),
                  xq2=quantile(databilu$d2015,0.95,na.rm=T),yq2=quantile(databilu$f2015,0.95,na.rm=T))


#4. RESERVOIR trend --------------

##4.1 Bivariate categorization at HER level level ---------------
ReservoirFloodTrend=FloodTrends[which(FloodTrends$driver==driver[3]),c(2:71)]
ReservoirDroughtTrend=DroughtTrends[which(DroughtTrends$driver==driver[3]),c(2:71)]
ReservoirFloodTrend$HydroR=FloodTrends[which(FloodTrends$driver==driver[1]),1]
ReservoirDroughtTrend$HydroR=DroughtTrends[which(DroughtTrends$driver==driver[1]),1]

FlPoint=full_join(GHshpp,ReservoirFloodTrend,by=c("CODEB"="HydroR"))
st_geometry(FlPoint)<-NULL
Flplot=inner_join(HydroRsf,FlPoint,by= c("Id"))
period=c(1951,2020)
haz="flood"

Flplot <- st_transform(Flplot, crs = 3035)

DrPoint=full_join(GHshpp,ReservoirDroughtTrend,by=c("CODEB"="HydroR"))
st_geometry(DrPoint)<-NULL
Drplot=inner_join(HydroRsf,DrPoint,by= c("Id"))
period=c(1951,2020)
haz="flood"
Drplot <- st_transform(Drplot, crs = 3035)

Flplot$d2015=Drplot$Rchange.Y2015
Flplot$f2015=Flplot$Rchange.Y2015

FlplotRes=Flplot
a=b=1
databire=FlplotRes
databire$x=databire$f2015/a
databire$y=databire$d2015/b

alterclass=data.frame(databire$x)
alterclass$class=NA
alterclass$class[which(alterclass[,1]<=(-breaker2))]=1
alterclass$class[which(alterclass[,1]>(-breaker2) & alterclass[,1]<(breaker1))]=2
alterclass$class[which(alterclass[,1]>=(breaker1) & alterclass[,1]<(breaker2))]=3
alterclass$class[which(alterclass[,1]>=breaker2)]=4

c1=alterclass$class
alterclass=data.frame(databire$y)
alterclass$class=NA
alterclass$class[which(alterclass[,1]<=(-breaker2))]=1
alterclass$class[which(alterclass[,1]>(-breaker2) & alterclass[,1]<(breaker1))]=2
alterclass$class[which(alterclass[,1]>=(breaker1) & alterclass[,1]<(breaker2))]=3
alterclass$class[which(alterclass[,1]>=breaker2)]=4

c2=alterclass$class

cx=paste(c2,c1,sep="-")

databire$bi_class=cx

# Combine the main category and subcategory to make a label
databire$combined_category <-databire$bi_class

##4.2 Bivariate categorization at pixel level ---------------
crsiflood=ResFloodTrendPix$Y2015
crsidrought=ResDroughtTrendPix$Y2015

a=b=1
databipire=Flpixplot
databipire$x=crsiflood/a
databipire$y=crsidrought/b




alterclass=data.frame(databipire$x)
alterclass$class=NA
alterclass$class[which(alterclass[,1]<=(-breaker2))]=1
alterclass$class[which(alterclass[,1]>(-breaker2) & alterclass[,1]<(breaker1))]=2
alterclass$class[which(alterclass[,1]>=(breaker1) & alterclass[,1]<(breaker2))]=3
alterclass$class[which(alterclass[,1]>=breaker2)]=4

c1=alterclass$class
alterclass=data.frame(databipire$y)
alterclass$class=NA
alterclass$class[which(alterclass[,1]<=(-breaker2))]=1
alterclass$class[which(alterclass[,1]>(-breaker2) & alterclass[,1]<(breaker1))]=2
alterclass$class[which(alterclass[,1]>=(breaker1) & alterclass[,1]<(breaker2))]=3
alterclass$class[which(alterclass[,1]>=breaker2)]=4

c2=alterclass$class

cx=paste(c2,c1,sep="-")

databipire$bi_class=cx

mbires=data.frame(x=mean(databire$d2015,na.rm=T),y=mean(databire$f2015,na.rm=T),
                  xq1=quantile(databire$d2015,0.05,na.rm=T),yq1=quantile(databire$f2015,0.05,na.rm=T),
                  xq2=quantile(databire$d2015,0.95,na.rm=T),yq2=quantile(databire$f2015,0.95,na.rm=T))


#5. Water Demand signal ------------------

##5.1 Bivariate categorization at HER level ---------------
WaterDemandFloodTrend=FloodTrends[which(FloodTrends$driver==driver[4]),c(2:71)]
WaterDemandDroughtTrend=DroughtTrends[which(DroughtTrends$driver==driver[4]),c(2:71)]
WaterDemandFloodTrend$HydroR=FloodTrends[which(FloodTrends$driver==driver[4]),1]
WaterDemandDroughtTrend$HydroR=DroughtTrends[which(DroughtTrends$driver==driver[4]),1]

FlPoint=full_join(GHshpp,WaterDemandFloodTrend,by=c("CODEB"="HydroR"))
st_geometry(FlPoint)<-NULL
Flplot=inner_join(HydroRsf,FlPoint,by= c("Id"))
period=c(1951,2020)
haz="flood"

Flplot <- st_transform(Flplot, crs = 3035)

DrPoint=full_join(GHshpp,WaterDemandDroughtTrend,by=c("CODEB"="HydroR"))
st_geometry(DrPoint)<-NULL
Drplot=inner_join(HydroRsf,DrPoint,by= c("Id"))
period=c(1951,2020)
haz="flood"
Drplot <- st_transform(Drplot, crs = 3035)

Flplot$d2015=Drplot$Rchange.Y2015
Flplot$f2015=Flplot$Rchange.Y2015

FlplotWD=Flplot
a=b=1
databiwd=FlplotWD
databiwd$x=databiwd$f2015/a
databiwd$y=databiwd$d2015/b

alterclass=data.frame(databiwd$x)
alterclass$class=NA
alterclass$class[which(alterclass[,1]<=(-breaker2))]=1
alterclass$class[which(alterclass[,1]>(-breaker2) & alterclass[,1]<(breaker1))]=2
alterclass$class[which(alterclass[,1]>=(breaker1) & alterclass[,1]<(breaker2))]=3
alterclass$class[which(alterclass[,1]>=breaker2)]=4

c1=alterclass$class
alterclass=data.frame(databiwd$y)
alterclass$class=NA
alterclass$class[which(alterclass[,1]<=(-breaker2))]=1
alterclass$class[which(alterclass[,1]>(-breaker2) & alterclass[,1]<(breaker1))]=2
alterclass$class[which(alterclass[,1]>=(breaker1) & alterclass[,1]<(breaker2))]=3
alterclass$class[which(alterclass[,1]>=breaker2)]=4

c2=alterclass$class

cx=paste(c2,c1,sep="-")

databiwd$bi_class=cx
# Combine the main category and subcategory to make a label
databiwd$combined_category <-databiwd$bi_class


##5.2 Bivariate categorization at pixel level ---------------
crsiflood=WuFloodTrendPix$Y2015
crsidrought=WuDroughtTrendPix$Y2015
a=sd(crsiflood,na.rm=T)
b=sd(crsidrought,na.rm=T)
a=b=1
databipiwd=Flpixplot
databipiwd$x=crsiflood/a
databipiwd$y=crsidrought/b


sd(databipiwd$x,na.rm=T)

alterclass=data.frame(databipiwd$x)
alterclass$class=NA
alterclass$class[which(alterclass[,1]<=(-breaker2))]=1
alterclass$class[which(alterclass[,1]>(-breaker2) & alterclass[,1]<(breaker1))]=2
alterclass$class[which(alterclass[,1]>=(breaker1) & alterclass[,1]<(breaker2))]=3
alterclass$class[which(alterclass[,1]>=breaker2)]=4

c1=alterclass$class
alterclass=data.frame(databipiwd$y)
alterclass$class=NA
alterclass$class[which(alterclass[,1]<=(-breaker2))]=1
alterclass$class[which(alterclass[,1]>(-breaker2) & alterclass[,1]<(breaker1))]=2
alterclass$class[which(alterclass[,1]>=(breaker1) & alterclass[,1]<(breaker2))]=3
alterclass$class[which(alterclass[,1]>=breaker2)]=4

c2=alterclass$class

cx=paste(c2,c1,sep="-")

databipiwd$bi_class=cx

mbiwd=data.frame(x=mean(databiwd$d2015,na.rm=T),y=mean(databiwd$f2015,na.rm=T),
                  xq1=quantile(databiwd$d2015,0.05,na.rm=T),yq1=quantile(databiwd$f2015,0.05,na.rm=T),
                  xq2=quantile(databiwd$d2015,0.95,na.rm=T),yq2=quantile(databiwd$f2015,0.95,na.rm=T))


# Bivariate plots of mean contribution ----------------------

databipicl=data.frame(databipic[,c(1,2,80:86)])
databipicl$driver="Climate"

databipilul=data.frame(databipilu[,c(1,2,80:86)])
databipilul$driver="LandUse"

databipirel=data.frame(databipire[,c(1,2,80:86)])
databipirel$driver="Reservoir"

databipiwdl=data.frame(databipiwd[,c(1,2,80:86)])
databipiwdl$driver="WaterDemand"

databipD=rbind(databipicl,databipilul,
               databipirel,databipiwdl)

databipD$bi_class[which(databipD$bi_class=="NA-NA")]=NA
databipD$bi_class[which(is.na(databipD$x))]=NA
databipD$bi_class[which(is.na(databipD$y))]=NA
unique(databipD$bi_class)
magg2 = aggregate(list(val=databipD$upa),
                 by = list(reg=databipD$driver ,traj=databipD$bi_class),
                 FUN = function(x) c(len=length(x)))
magg2 <- do.call(data.frame, magg2)
# magg2=magg2[-which(magg2$traj=="NA-1" | magg2$traj=="NA-2" |
#                   magg2$traj=="NA-3"),]
magg2$trcat="Wetting"
magg2$trcat[which(magg2$traj=="2-2" | magg2$traj=="3-2" |
                   magg2$traj=="2-3" | magg2$traj=="3-3")]="Stable"
magg2$trcat[which(magg2$traj=="1-1" | magg2$traj=="1-2" |
                   magg2$traj=="2-1")]="Drying"
magg2$trcat[which(magg2$traj=="1-4" | magg2$traj=="2-4" |
                   magg2$traj=="1-3")]="Accelerating"
magg2$trcat[which(magg2$traj=="4-1" | magg2$traj=="4-2" |
                   magg2$traj=="3-1")]="Decelerating"
magg2$order=1
magg2$order[which(magg2$trcat=="Wetting")]=2
magg2$order[which(magg2$trcat=="Decelerating")]=3
magg2$order[which(magg2$trcat=="Drying")]=4
magg2$order[which(magg2$trcat=="Accelerating")]=5

matshit=match(magg2$traj,names(colord))
magg2$ordf=colord[matshit]
magg2 <- magg2 %>% 
  mutate(traj = reorder(traj, ordf, FUN = mean))

loscolors=c("Accelerating" = "#174f28","Drying" = "#dd6a29","Stable"="gray60","Wetting" = "#169dd0","Decelerating" = "burlywood")
magg3=magg2[which(magg2$reg=="LandUse"),]
ggplot(magg3, aes(x = trcat, y = val, fill = trcat)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = loscolors) +
  scale_x_discrete(name="Trajectory")+
  #scale_y_continuous(label=c(0,25,50,75,100),name="River pixels (%)")+
  theme(axis.title=element_text(size=tsize),
        panel.background = element_rect(fill = "white", colour = "grey1"),
        panel.border = element_rect(linetype = "solid", fill = NA, colour="black"),
        legend.title = element_text(size=tsize),
        legend.text = element_text(size=osize),
        legend.position = "none",
        legend.key = element_rect(fill = "transparent", colour = "transparent"),
        legend.key.size = unit(.8, "cm"))



### match biogeoregions with HRs ----
databitotxHR=databipi
RegioRLi=aggregate(list(val=databitotxHR$y),
                   by = list(HydroR=databitotxHR$HydroRegions_raster_WGS84),
                   FUN = function(x) c(mean=mean(x,na.rm=T),dev=sd(x,na.rm=T),len=length(x),med=median(x,na.rm=T),q1=quantile(x, 0.05, na.rm=T),q3=quantile(x, 0.95, na.rm=T)))
RegioRLi <- do.call(data.frame, RegioRLi)

HRM=na.omit((match(RegioRLi$HydroR,HydroRsf$Id)))
HydroRsf_dom=HydroRsf[HRM,]


layerMatch=match(GHR_riv$outl2,biogeo_rivers$outl2)

GHR_riv$biogeoR=biogeo_rivers$code[layerMatch]
unique(GHR_riv$HydroRegions_raster_WGS84)
verif=(match(GHR_riv$HydroRegions_raster_WGS84,HydroRsf_dom$Id))
unique(verif)
unique(HydroRsf_dom$Id[verif])
GHR_riv$HydrRName=HydroRsf_dom$IRST_NAMEB[verif]


bhp_m=na.omit(unique(match(HydroRsf_dom$Id,GHR_riv$HydroRegions_raster_WGS84)))
HydroRsf_dom$biogeoreg=GHR_riv$biogeoR[bhp_m]

sample=GHR_riv[bhp_m,]


RegioAg = aggregate(list(oc=GHR_riv$HydrRName),
                    by = list(BR=GHR_riv$biogeoR,HR=GHR_riv$HydroRegions_raster_WGS84),
                    FUN = function(x) c(len=length(x)))
RegioAg <- do.call(data.frame, RegioAg)

ur=unique(RegioAg$HR)

### matching biogeoregions and pixel level estimate -----
databiwdxBG=inner_join(biogeo_rivers,databipiwd,by=c("outl2"))
databiluxBG=inner_join(biogeo_rivers,databipilu,by=c("outl2"))
databiclixBG=inner_join(biogeo_rivers,databipic,by=c("outl2"))
databiresxBG=inner_join(biogeo_rivers,databipire,by=c("outl2"))
databitotxBG=inner_join(biogeo_rivers,databipi,by=c("outl2"))


### matching hydroregions and pixel level estimates ----
databiwdxHR=inner_join(GHR_riv,databipiwd,by=c("outl2"))
databiluxHR=inner_join(GHR_riv,databipilu,by=c("outl2"))
databiclixHR=inner_join(GHR_riv,databipic,by=c("outl2"))
databiresxHR=inner_join(GHR_riv,databipire,by=c("outl2"))
databitotxHR=inner_join(GHR_riv,databipi,by=c("outl2"))


mbf=rbind(mbicli,mbires,mbilu,mbiwd)
mbf$names=c('Climate',
            'Reservoirs',
            'Landuse',
            'WaterDemand')


#### [SPlot] - Supplementary plot - Bivariate change aggregate -----
loscolors=c("Accelerating" = "#174f28","Drying" = "#dd6a29","Stable"="gray60","Wetting" = "#169dd0","Decelerating" = "burlywood")
clabels=c("Climate","Land use","Reservoirs", "Water demand")
colorn = c("WaterDemand" ='limegreen',"Reservoirs" ='tomato4',"Landuse" ='orange',"Climate" ='royalblue')
tsize=osize=20
bpc=ggplot(mbf, aes(x = x, y = y, fill = names)) +
  geom_linerange(data=mbf,aes(x=x, ymin=yq1,ymax=yq2,color = names,group=factor(names)),lwd=1,alpha=0.8) +
  geom_linerange(data=mbf,aes(y=y, xmin=xq1,xmax=xq2,color = names,group=factor(names)),lwd=1,alpha=0.8) +
  geom_point(shape = 21, color = "black", size = 5, show.legend = TRUE) +
  scale_fill_manual(values = colorn, name = "Drivers", labels = clabels) +
  scale_color_manual(values = colorn, name = "Drivers", labels = clabels) +
  coord_cartesian(xlim=c(-100,100),ylim=c(-50,50)) +
  annotate("text", x = 50, y = 20, label = "Wetting", color = "#169dd0", size = 7, fontface = "bold") +
  annotate("text", x = -50, y = 20, label = "Accelerating", color = "#174f28", size = 7, fontface = "bold") +
  annotate("text", x = -50, y = -20, label = "Drying", color = "#dd6a29", size = 7, fontface = "bold") +
  annotate("text", x = 50, y = -20, label = "Decelerating", color = "burlywood", size = 7, fontface = "bold") +
  labs(x = "Change in drought (%)", 
       y = "Change in flood (%)") +
  scale_x_continuous(trans=scales::modulus_trans(.5),
                     breaks=seq(-100,100,50),
                     minor_breaks = seq(-100,100,10)) +
  scale_y_continuous(trans=scales::modulus_trans(.5),
                     breaks=seq(-100,100,25),
                     minor_breaks = seq(-100,100,10)) +
  guides( fill = guide_legend(override.aes = list(size = 6)))+
  theme(
    axis.title = element_text(size = 18, face = "bold"),
    title = element_text(size = 22, face = "bold"),
    axis.text = element_text(size = 16, face = "bold"),
    panel.background = element_rect(fill = "white", colour = "white"),
    panel.grid = element_blank(),
    panel.border = element_rect(linetype = "solid", fill = NA, colour = "black"),
    legend.title = element_text(size = 20, face = "bold"),
    legend.text = element_text(size = 16),
    #legend.key = element_rect(fill = "lightgray", colour = "lightgray"),
    legend.position = "right",
    panel.grid.major = element_line(colour = "grey60"),
    panel.grid.minor= element_line(colour = "grey70", linetype = "dashed"),
    #legend.key = element_rect(fill = "white", colour = "white"),
    #legend.key.size = unit(1, "cm")
  )

bpc
ggsave(paste0("D:/tilloal/Documents/LFRuns_utils/TrendAnalysis/plots/Contribution_bv.jpg"),bpc, width=25, height=20, units=c("cm"),dpi=300) 



#6. Bivariate change by driver Aggregated to HER level ----

##6.1 Aggregation of pixel level results to HER level ----

###6.1.1 Aggregation of drought changes ----
mbiwd1=aggregate(list(val=databiwdxHR$y.y),
                 by = list(HR=databiwdxHR$HydrRName),
                 FUN = function(x) c(dmean=mean(x,na.rm=T),dmed=median(x,na.rm=T),dq1=quantile(x,0.05,na.rm=T),dq2=quantile(x,0.95,na.rm=T),l=length(x),sd=sd(x,na.rm=T)))
mbiwd1 <- do.call(data.frame, mbiwd1)

mbicli1=aggregate(list(val=databiclixHR$y.y),
                  by = list(HR=databiclixHR$HydrRName),
                  FUN = function(x) c(dmean=mean(x,na.rm=T),dmed=median(x,na.rm=T),dq1=quantile(x,0.05,na.rm=T),dq2=quantile(x,0.95,na.rm=T),l=length(x),sd=sd(x,na.rm=T)))
mbicli1 <- do.call(data.frame, mbicli1)

mbires1=aggregate(list(val=databiresxHR$y.y),
                  by = list(HR=databiresxHR$HydrRName),
                  FUN = function(x) c(dmean=mean(x,na.rm=T),dmed=median(x,na.rm=T),dq1=quantile(x,0.05,na.rm=T),dq2=quantile(x,0.95,na.rm=T),l=length(x),sd=sd(x,na.rm=T)))
mbires1 <- do.call(data.frame, mbires1)

mbilu1=aggregate(list(val=databiluxHR$y.y),
                 by = list(HR=databiluxHR$HydrRName),
                 FUN = function(x) c(dmean=mean(x,na.rm=T),dmed=median(x,na.rm=T),dq1=quantile(x,0.05,na.rm=T),dq2=quantile(x,0.95,na.rm=T),l=length(x),sd=sd(x,na.rm=T)))
mbilu1 <- do.call(data.frame, mbilu1)

mball1=aggregate(list(val=databitotxHR$y.y),
                 by = list(HR=databitotxHR$HydrRName),
                 FUN = function(x) c(dmean=mean(x,na.rm=T),dmed=median(x,na.rm=T),dq1=quantile(x,0.05,na.rm=T),dq2=quantile(x,0.95,na.rm=T),l=length(x),sd=sd(x,na.rm=T)))
mball1 <- do.call(data.frame, mball1)

###6.1.2 Aggregation of flood changes ----
mbiwd2=aggregate(list(val=databiwdxHR$x.y),
                 by = list(HR=databiwdxHR$HydrRName),
                 FUN = function(x) c(fmean=mean(x,na.rm=T),fq1=quantile(x,0.05,na.rm=T),fq2=quantile(x,0.95,na.rm=T),l=length(x),sd=sd(x,na.rm=T)))
mbiwd2 <- do.call(data.frame, mbiwd2)

mbicli2=aggregate(list(val=databiclixHR$x.y),
                  by = list(HR=databiclixHR$HydrRName),
                  FUN = function(x) c(fmean=mean(x,na.rm=T),fq1=quantile(x,0.05,na.rm=T),fq2=quantile(x,0.95,na.rm=T),l=length(x),sd=sd(x,na.rm=T)))
mbicli2 <- do.call(data.frame, mbicli2)

mbires2=aggregate(list(val=databiresxHR$x.y),
                  by = list(HR=databiresxHR$HydrRName),
                  FUN = function(x) c(fmean=mean(x,na.rm=T),fq1=quantile(x,0.05,na.rm=T),fq2=quantile(x,0.95,na.rm=T),l=length(x),sd=sd(x,na.rm=T)))
mbires2 <- do.call(data.frame, mbires2)

mbilu2=aggregate(list(val=databiluxHR$x.y),
                 by = list(HR=databiluxHR$HydrRName),
                 FUN = function(x) c(fmean=mean(x,na.rm=T),fq1=quantile(x,0.05,na.rm=T),fq2=quantile(x,0.95,na.rm=T),l=length(x),sd=sd(x,na.rm=T)))
mbilu2 <- do.call(data.frame, mbilu2)

mball2=aggregate(list(val=databitotxHR$x.y),
                 by = list(HR=databitotxHR$HydrRName),
                 FUN = function(x) c(fmean=mean(x,na.rm=T),fq1=quantile(x,0.05,na.rm=T),fq2=quantile(x,0.95,na.rm=T),l=length(x),sd=sd(x,na.rm=T)))
mball2 <- do.call(data.frame, mball2)


##6.2 Joint changes by region ----

####merge flood and drought changes ----
mbiwd3=inner_join(mbiwd1,mbiwd2,by="HR")
mbires3=inner_join(mbires1,mbires2,by="HR")
mbilu3=inner_join(mbilu1,mbilu2,by="HR")
mbicli3=inner_join(mbicli1,mbicli2,by="HR")
mball3=inner_join(mball1,mball2,by="HR")

RegionName=inner_join(mbicli2,HydroRsf,by=c("HR"="IRST_NAMEB"))
l1=length(mbicli2$HR)

####create one data.frame containing all changes from different drivers ----
mbfX=rbind(mbicli3,mbires3,mbilu3,mbiwd3,mball3)
mbfX$names=c(rep('Climate',l1),
            rep('Reservoirs',l1),
            rep('Landuse',l1),
            rep('WaterDemand',l1),
            rep('Historical',l1))

zeb=match(mbfX$HR,RegionName$HR)
mbfX$regioname=RegionName$IRST_NAMEB[zeb]

####correlation between change in drought and flood of several drivers ----
cor.test(mbicli3$val.dmean, mbicli3$val.fmean)
cor.test(mbires3$val.dmean, mbires3$val.fmean)
cor.test(mbilu3$val.dmean, mbilu3$val.fmean)


colnames(mbfX)[c(2,4,5,8,9,10)]=c("x","xq1","xq2","y","yq1","yq2")
mbfX$xl=((mbfX$xq2-mbfX$xq1)+(mbfX$yq2-mbfX$yq1))/2

####only data for Historical changes ----
mbfH=mbfX[which(mbfX$names=="Historical"),]
####only data for  changes ----
mbfX=mbfX[-which(mbfX$names=="Historical"),]

loscolors=c("Accelerating" = "#174f28","Drying" = "#dd6a29","Stable"="gray60","Wetting" = "#169dd0","Decelerating" = "burlywood")
clabels=c("Climate","Land use","Reservoirs", "Water demand")
colorn = c("WaterDemand" ='limegreen',"Reservoirs" ='tomato4',"Landuse" ='orange',"Climate" ='royalblue')
shapex=c("Alpine"=0,"Atlantic"=1,"Boreal"=2, "Continental"=3,"Mediterranean"=4)

tsize=osize=20
lalim=70

## [Plot] - Figure 5 - Bivariate plot of driver contribution to hydro extreme changes ---- 
bpc=ggplot() +
  geom_vline(xintercept = 0, 
  color = "black", size = 1.5) +
  geom_hline(yintercept = 0, 
             color = "black", size = 1.5) +
  geom_point(data=mbfX, aes(x = x, y = y, fill = names, size=val.l.x), shape = 21, alpha=.5, stroke=0, show.legend = TRUE) +
  geom_point(data=mbfX, aes(x = x, y = y, color = names, size=val.l.x), shape = 1, stroke=.8, show.legend = FALSE, alpha=.6) +
  scale_fill_manual(values = colorn, name = "Drivers", labels = clabels) +
  scale_color_manual(values = colorn, name = "Drivers", labels = clabels) +
  scale_size(range = c(2, 10))+
  coord_cartesian(xlim=c(-lalim,lalim),ylim=c(-lalim,lalim)) +
  annotate("text", x = 35, y = 45, label = "Wetting", color = "#169dd0", size = 6, fontface = "bold") +
  annotate("text", x = -35, y = 45, label = "Accelerating", color = "#174f28", size = 6, fontface = "bold") +
  annotate("text", x = -35, y = -45, label = "Drying", color = "#dd6a29", size = 6, fontface = "bold") +
  annotate("text", x = 35, y = -45, label = "Decelerating", color = "burlywood", size = 6, fontface = "bold") +
  labs(x = "Change in drought intensity (s%)", 
       y = "Change in flood intensity (s%)") +
  scale_x_continuous(trans=scales::modulus_trans(.5),
                     breaks=c(-200,-100,-50,-20,-10,-5,0,5,10,20,50,100,200),
                     labels=rev(c(-200,-100,-50,-20,-10,-5,0,5,10,20,50,100,200)),
                     minor_breaks = seq(-100,100,10)) +
  scale_y_continuous(trans=scales::modulus_trans(.5),
                     breaks=c(-200,-100,-50,-20,-10,-5,0,5,10,20,50,100,200),
                     minor_breaks = seq(-100,100,10), expand = c(0,0)) +
  guides( color = guide_legend(override.aes = list(size = 6)),size="none")+
  theme(
    axis.title = element_text(size = 18, face = "bold"),
    title = element_text(size = 22, face = "bold"),
    axis.text = element_text(size = 16, face = "bold"),
    panel.background = element_rect(fill = "white", colour = "white"),
    panel.grid = element_blank(),
    panel.border = element_rect(linetype = "solid", fill = NA, colour = "black"),
    legend.title = element_text(size = 20, face = "bold"),
    legend.text = element_text(size = 16),
    legend.position = "bottom",
    panel.grid.major = element_line(colour = "grey60"),
    panel.grid.minor= element_line(colour = "grey70", linetype = "dashed"),
  )

bpc
ggsave(paste0("D:/tilloal/Documents/LFRuns_utils/TrendAnalysis/plots/Contribution_HR_11.jpg"),bpc, width=20, height=20, units=c("cm"),dpi=500) 



#Temporal evolution per region
mbicli=aggregate(list(val=databiclixHR$y.y),
                  by = list(HR=databiclixHR$HydrRName),
                  FUN = function(x) c(dmean=mean(x,na.rm=T),dmed=median(x,na.rm=T),dq1=quantile(x,0.05,na.rm=T),dq2=quantile(x,0.95,na.rm=T),l=length(x),sd=sd(x,na.rm=T)))
mbicli1 <- do.call(data.frame, mbicli1)


mb1=mbfX[which(mbfX$bg=="Mediterranean"),]
mb1=mbfX[which(mbfX$bg=="Continental"),]

length(unique(GHshpp$IRST_NAMEB))

#plot by driver for different trajectories

mbfX$traj="none"
mbfX$traj[which(mbfX$y>=0 & mbfX$x>0)]="Wetting"

mbfX$traj[which(mbfX$y<=0 & mbfX$x>=0)]="Decelerating"

mbfX$traj[which(mbfX$y<=0 & mbfX$x<=0)]="Drying"

mbfX$traj[which(mbfX$y>=0 & mbfX$x<=0)]="Accelerating"


#Aggregation of trajectories by regions
magg2 = aggregate(list(val=mbfX$HR),
                  by = list(reg=mbfX$names,traj=mbfX$traj),
                  FUN = function(x) c(len=length(x)))
magg2 <- do.call(data.frame, magg2)
magg2$perc=magg2$val/s1*100

clabels=c("Climate","Landuse","Reservoirs", "WaterDemand")
atraj=c("Wetting","Decelerating","Drying","Accelerating")

## [Plot] - Figure 5.b - Proportion of trajectories by drivers -----
for (driv in atraj){
  magg3=magg2[which(magg2$traj==driv),]

  mierda=match(clabels,magg3$reg)
  addup=clabels[which(is.na(mierda))]
  if (length(addup)>0){
    mad=c(addup,driv,0)
    magg3=rbind(magg3,mad)
    magg3$val=as.numeric(magg3$val)
  }
  magg3$valp=magg3$val/length(mbfH$HR)*100
  ggplot(magg3, aes(x = reg, y = valp, fill = reg)) +
    geom_bar(stat = "identity") +
    scale_fill_manual(values = colorn) +
    scale_x_discrete(name=" ")+
    scale_y_continuous("HydroRegion (%)", limits=c(0,100))+
    theme(axis.title=element_text(size=16),
          axis.text = element_text(size=12),
          panel.background = element_rect(fill = "white", colour = "grey1"),
          panel.border = element_rect(linetype = "solid", fill = NA, colour="black"),
          legend.title = element_text(size=tsize),
          legend.text = element_text(size=osize),
          legend.position = "none",
          plot.title=element_text(size=20),
          legend.key = element_rect(fill = "transparent", colour = "transparent"),
          legend.key.size = unit(.8, "cm"))+
    ggtitle(driv)
  
  ggsave(paste0("D:/tilloal/Documents/LFRuns_utils/TrendAnalysis/plots/driverxtraj3_",driv,".jpg"), width=12, height=8, units=c("cm"),dpi=300) 
  

}


##6.3 plot of historical change with colors by biogeoregions ----
Hydroplot=c()
for (id in 1:length(ur)){
  myr=ur[id]
  maty=which(!is.na((match(RegioAg$HR,myr))))
  lar=RegioAg[maty,]
  if(length(lar$oc)>1){
    sumi=sum(lar$oc)
    lar=lar[which.max(lar$oc),]
    lar$oc=lar$oc/sumi*100
  }else{  lar$oc=100}
  Hydroplot=rbind(Hydroplot,lar)
} 

bioplot <- inner_join(HydroRsf_dom,Hydroplot, by=c("Id"="HR"))
bioplot$BR[which(bioplot$BR=="Pannonian")]="Continental"
bioplot$BR[which(bioplot$BR=="Steppic")]="Continental"
br=c(-50,-20,-10,0,10,20,50)
labels=br
limi=c(-60,60)
colx=c("Alpine"="deeppink","Atlantic"="forestgreen","Boreal"="darkviolet", "Continental"="orange3","Mediterranean"="gold")

###[Plot] - Extended Figure 1 - Bioregions and HydroRegions map ----
ggplot(basemap) +
  geom_sf(fill="gray95")+
  geom_sf(fill=NA, color="grey") +
  geom_sf(data=bioplot,aes(fill=factor(BR),geometry=geometry),color="black",alpha=.7,size=0.0001)+ 
  coord_sf(xlim = c(min(nco[,1]),max(nco[,1])), ylim = c(min(nco[,2]),max(nco[,2])))+
  scale_fill_manual(values = colx, name="BioGeoRegions") +
  labs(x="Longitude", y = "Latitude")+
  guides(color = guide_colourbar(barwidth = 12, barheight = 1))+
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

ggsave(paste0("D:/tilloal/Documents/LFRuns_utils/TrendAnalysis/plots/Map_HRxBG2.jpg"), width=25, height=20, units=c("cm"),dpi=300) 

h2p=match(mbfH$HR,bioplot$IRST_NAMEB)
mbfH$biogeo=bioplot$BR[h2p]
mbfH$CODEB=bioplot$CODEB[h2p]


#Significance of trends - not used -
n2t=match(trendFlood$HydroR ,bioplot$CODEB)
trendFlood$Hname=bioplot$IRST_NAMEB[n2t]
s2p=match(mbfH$HR,trendFlood$Hname)
mbfH$fsig=trendFlood$change[s2p]

length(trendDrought$HydroR[which(trendDrought$change==2)])/length(trendDrought$HydroR)
n2t=match(trendDrought$HydroR ,bioplot$CODEB)
trendDrought$Hname=bioplot$IRST_NAMEB[n2t]
s2p=match(mbfH$HR,trendDrought$Hname)
mbfH$dsig=trendDrought$change[s2p]

mbfH$fdsig=0
mbfH$fdsig=abs(mbfH$fsig)+abs(mbfH$dsig)

mbfHsig=mbfH[which(mbfH$fdsig==4),]


mbfHuns=mbfH[which(mbfH$fdsig==3),]

###[Splot] - Total Joint changes by HER and BGR ----
colx=c("Alpine"="deeppink","Atlantic"="forestgreen","Boreal"="darkviolet", "Continental"="orange3","Mediterranean"="gold")
bpc=ggplot() +
  geom_vline(xintercept = 0, 
             color = "black", size = 1.5) +
  geom_hline(yintercept = 0, 
             color = "black", size = 1.5) +

  geom_point(data=mbfH, aes(x = x, y = y, color=biogeo,size=val.l.x), show.legend = TRUE, stroke=0, alpha=.7) +
  geom_text(data=mbfH, aes(x=x, y=y,label = CODEB), size = 3, color = "black",fontface = "bold") +
  geom_point(data=mbfHsig, aes(x = x, y = y,size=val.l.x), fill=NA,show.legend = F, shape=1, stroke=1, col="black", alpha=.8) +
  scale_color_manual(values = colx, name = "Regions") +
  scale_size(range = c(6, 14))+
  coord_cartesian(xlim=c(-lalim,lalim+10),ylim=c(-lalim,lalim)) +
  annotate("text", x = 40, y = 45, label = "Wetting", color = "#169dd0", size = 6, fontface = "bold") +
  annotate("text", x = -35, y = 45, label = "Accelerating", color = "#174f28", size = 6, fontface = "bold") +
  annotate("text", x = -35, y = -45, label = "Drying", color = "#dd6a29", size = 6, fontface = "bold") +
  annotate("text", x = 40, y = -45, label = "Decelerating", color = "burlywood", size = 6, fontface = "bold") +
  labs(x = "Change in drought intensity (%)", 
       y = "Change in flood intensity (%)") +
  scale_x_continuous(trans=scales::modulus_trans(.5),
                     breaks=c(-200,-100,-50,-20,-10,-5,0,5,10,20,50,100,200),
                     labels=rev(c(-200,-100,-50,-20,-10,-5,0,5,10,20,50,100,200)),
                     minor_breaks = seq(-100,100,10)) +
  scale_y_continuous(trans=scales::modulus_trans(.5),
                     breaks=c(-200,-100,-50,-10,-5,0,5,10,50,100,200),
                     minor_breaks = seq(-100,100,10), expand = c(0,0)) +
  guides( color = guide_legend(override.aes = list(size = 6)),size="none")+
  theme(
    axis.title = element_text(size = 18, face = "bold"),
    title = element_text(size = 22, face = "bold"),
    axis.text = element_text(size = 16, face = "bold"),
    panel.background = element_rect(fill = "white", colour = "white"),
    panel.grid = element_blank(),
    panel.border = element_rect(linetype = "solid", fill = NA, colour = "black"),
    legend.title = element_text(size = 20, face = "bold"),
    legend.text = element_text(size = 16),
    legend.position = "right",
    panel.grid.major = element_line(colour = "grey60"),
    panel.grid.minor= element_line(colour = "grey70", linetype = "dashed") )

bpc
ggsave(paste0("D:/tilloal/Documents/LFRuns_utils/TrendAnalysis/plots/Contribution_HRxBG_7.jpg"),bpc, width=25, height=20, units=c("cm"),dpi=300) 

###[SPlot] - histogram per BGR of count of HER per categories ----
mbfH$fdg=mbfH$fsig+mbfH$dsig
mbfH$traj=0
#wetting
mbfH$traj[which(mbfH$fsig>0 & mbfH$dsig>0)]=1
#signifcant wetting
mbfH$traj[which(mbfH$traj==1 & mbfH$fdsig==4)]=2

#wetting
mbfH$traj[which(mbfH$fsig>0 & mbfH$dsig>0)]=1
#signifcant wetting
mbfH$traj[which(mbfH$traj==1 & mbfH$fdsig==4)]=2

#decelerating
mbfH$traj[which(mbfH$fsig<0 & mbfH$dsig>0)]=3
#signifcant decelerating
mbfH$traj[which(mbfH$traj==3 & mbfH$fdsig==4)]=4

#drying
mbfH$traj[which(mbfH$fsig<0 & mbfH$dsig<0)]=5
#signifcant drying
mbfH$traj[which(mbfH$traj==5 & mbfH$fdsig==4)]=6

#accelerating
mbfH$traj[which(mbfH$fsig>0 & mbfH$dsig<0)]=7
#signifcant accelerating
mbfH$traj[which(mbfH$traj==7 & mbfH$fdsig==4)]=8
#aggregation
magg = aggregate(list(val=mbfH$x),
                 by = list(reg=mbfH$biogeo,traj=mbfH$traj),
                 FUN = function(x) c(len=length(x)))
magg <- do.call(data.frame, magg)

mas = aggregate(list(val=mbfH$x),
                 by = list(reg=mbfH$biogeo),
                 FUN = function(x) c(len=length(x)))
mas <- do.call(data.frame, mas)
mmt=match(magg$reg,mas$reg)
magg$rlen=mas$val[mmt]
magg$val2=magg$val/magg$rlen*100

paet=c("lightblue","royalblue", "beige","lightgoldenrod",
       "orange","darkorange","lightgreen","forestgreen")
labeli=c("wetting","sig. wetting","decelerating","sig. decelerating",
         "drying","sig. drying","accelerating","sig. accelerating")
ggplot(magg, aes(x = reg, y = val2, fill = factor(traj))) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(
    values=paet, labels=labeli)



#Percedntage of regions in each category
v=0
#wetting
length(which(mbfH$x>v & mbfH$y>v))/length(mbfH$HR)*100
#accelerating
length(which(mbfH$x<(-v) & mbfH$y>v))/length(mbfH$HR)*100
#decelerating
length(which(mbfH$x>v & mbfH$y<(-v)))/length(mbfH$HR)*100
#drying
length(which(mbfH$x<(-v) & mbfH$y<(-v)))/length(mbfH$HR)*100

mbfC=mbfX[which(mbfX$names=="Climate"),]

#wetting
length(which(mbfC$x>v & mbfC$y>v))/length(mbfC$HR)*100
#accelerating
length(which(mbfC$x<(-v) & mbfC$y>v))/length(mbfC$HR)*100
#decelerating
length(which(mbfC$x>v & mbfC$y<(-v)))/length(mbfC$HR)*100
#drying
length(which(mbfC$x<(-v) & mbfC$y<(-v)))/length(mbfC$HR)*100

mbfL=mbfX[which(mbfX$names=="Landuse"),]
length(which(mbfL$x<v & mbfL$y>v))/length(mbfL$HR)*100

#wetting
length(mbfHsig$HR)/length(mbfH$HR)*100
length(which(mbfHsig$x>v & mbfHsig$y>v))/length(mbfH$HR)*100
#accelerating
length(which(mbfHsig$x<(-v) & mbfHsig$y>v))/length(mbfH$HR)*100
#decelerating
length(which(mbfHsig$x>v & mbfHsig$y<(-v)))/length(mbfH$HR)*100
#drying
length(which(mbfHsig$x<(-v) & mbfHsig$y<(-v)))/length(mbfH$HR)*100

mbfR=mbfX[which(mbfX$names=="Reservoirs"),]
length(which(mbfL$x>v & mbfL$y>v))/length(mbfL$HR)*100


##6.4 Saving outputs at HER level ----
write.csv(mbfH,file=paste0("D:/tilloal/Documents/LFRuns_utils/TrendAnalysis/Histo_res_bvHR_6.csv"))
write.csv(mbfX,file=paste0("D:/tilloal/Documents/LFRuns_utils/TrendAnalysis/Drivers_res_bvHR_6.csv"))




#Standard error plot
st_geometry(databiclim) <- NULL
databic2=inner_join(mbfH,GHshpp,by=c("HR"="IRST_NAMEB"))
#databic2=full_join(databic2,mbfH,by=c("IRST_NAMEB.y"="HR"))
# databic2$xdr=databic2$val.l.x/length(databiclixBG$PK_UID)*(databic2$val.sd.x)^2
hist(abs(databic2$x))
databic2$xdr=databic2$val.sd.x
databic2$xfl=2*(databic2$val.sd.x/sqrt(databic2$val.l.x))*1.96
plot(databic2$xdr)
palet=c(hcl.colors(11, palette = "YlOrRd", alpha = NULL, rev = T, fixup = TRUE))

map <- ggplot(basemap) +
  geom_sf(fill="white")+
  geom_sf(data = databic2, mapping = aes(fill = xdr, geometry=geometry), alpha=1) +
  geom_sf(fill=NA, color="gray42") +
  # bi_scale_fill(pal = "BlueOr", dim = 3, na.value=colNA ) +
  # scale_fill_manual(values = colorp) +
  scale_fill_gradientn(
    colors=palet, name="Standard Deviation", limits=c(0,50), oob = scales::squish)   +
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
        legend.position = "right",
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

###[Plot]- Extended Figure 1 - hydroregions plot
color = grDevices::colors()[grep('gr(a|e)y|white|black', grDevices::colors(), invert = T)]
coHR=sample(color, length(unique(databic2$HR)))
map <- ggplot(basemap) +
  geom_sf(fill="white")+
  geom_sf(data = databic2, mapping = aes(fill = factor(HR),geometry=geometry), alpha=0.7, color = "black", size = 0.001,show.legend = FALSE) +
  geom_sf(fill=NA, color="gray42") +
  # bi_scale_fill(pal = "BlueOr", dim = 3, na.value=colNA ) +
  # scale_fill_gradientn(
  #   colors=palet, name="95% CI of mean", limits=c(0,10), oob = scales::squish)   +
  coord_sf(xlim = c(min(nco[,1]),max(nco[,1])), ylim = c(min(nco[,2]),max(nco[,2])))+
  #bi_scale_color(pal = "BlueOr", dim = 3, na.value=colNA ) +
  scale_fill_manual(values = coHR) +
  labs()+
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

ggsave(paste0("D:/tilloal/Documents/LFRuns_utils/TrendAnalysis/plots/mapreg_HR.jpg"), map, width=30, height=20, units=c("cm"),dpi=300) 


#7 Combination of drivers -----

##7.1 HER Level ----

###7.1.1 clim+reservoirs ----
databicr=databiclim
databicr$x= databire$x + databiclim$x
databicr$y= databire$y + databiclim$y

breaker1=0
breaker2=5

alterclass=data.frame(databicr$x)
alterclass$class=NA
alterclass$class[which(alterclass[,1]<=(-breaker2))]=1
alterclass$class[which(alterclass[,1]>(-breaker2) & alterclass[,1]<(breaker1))]=2
alterclass$class[which(alterclass[,1]>=(breaker1) & alterclass[,1]<(breaker2))]=3
alterclass$class[which(alterclass[,1]>=breaker2)]=4
c1=alterclass$class

alterclass=data.frame(databicr$y)
alterclass$class=NA
alterclass$class[which(alterclass[,1]<=(-breaker2))]=1
alterclass$class[which(alterclass[,1]>(-breaker2) & alterclass[,1]<(breaker1))]=2
alterclass$class[which(alterclass[,1]>=(breaker1) & alterclass[,1]<(breaker2))]=3
alterclass$class[which(alterclass[,1]>=breaker2)]=4
c2=alterclass$class

cx=paste(c2,c1,sep="-")
databicr$bi_class=cx

# Combine the main category and subcategory to make a label
databicr$combined_category <-databicr$bi_class
databicrl=databiclim

###7.1.2 clim+reservoirs+lu ----
databicrl$x=databilu$x + databire$x + databiclim$x
databicrl$y=databilu$y + databire$y + databiclim$y

alterclass=data.frame(databicrl$x)
alterclass$class=NA
alterclass$class[which(alterclass[,1]<=(-breaker2))]=1
alterclass$class[which(alterclass[,1]>(-breaker2) & alterclass[,1]<(breaker1))]=2
alterclass$class[which(alterclass[,1]>=(breaker1) & alterclass[,1]<(breaker2))]=3
alterclass$class[which(alterclass[,1]>=breaker2)]=4
c1=alterclass$class

alterclass=data.frame(databicrl$y)
alterclass$class=NA
alterclass$class[which(alterclass[,1]<=(-breaker2))]=1
alterclass$class[which(alterclass[,1]>(-breaker2) & alterclass[,1]<(breaker1))]=2
alterclass$class[which(alterclass[,1]>=(breaker1) & alterclass[,1]<(breaker2))]=3
alterclass$class[which(alterclass[,1]>=breaker2)]=4
c2=alterclass$class
cx=paste(c2,c1,sep="-")
databicrl$bi_class=cx
# Combine the main category and subcategory to make a label
databicrl$combined_category <-databicrl$bi_class


###7.1.3 clim+reservoirs+lu+wd ----
databicrlw=databiclim

databicrlw$x=databilu$x + databire$x + databiclim$x + databiwd$x
databicrlw$y=databilu$y + databire$y + databiclim$y + databiwd$y


alterclass=data.frame(databicrlw$x)
alterclass$class=NA
alterclass$class[which(alterclass[,1]<=(-breaker2))]=1
alterclass$class[which(alterclass[,1]>(-breaker2) & alterclass[,1]<(breaker1))]=2
alterclass$class[which(alterclass[,1]>=(breaker1) & alterclass[,1]<(breaker2))]=3
alterclass$class[which(alterclass[,1]>=breaker2)]=4
c1=alterclass$class

alterclass=data.frame(databicrlw$y)
alterclass$class=NA
alterclass$class[which(alterclass[,1]<=(-breaker2))]=1
alterclass$class[which(alterclass[,1]>(-breaker2) & alterclass[,1]<(breaker1))]=2
alterclass$class[which(alterclass[,1]>=(breaker1) & alterclass[,1]<(breaker2))]=3
alterclass$class[which(alterclass[,1]>=breaker2)]=4
c2=alterclass$class
cx=paste(c2,c1,sep="-")
databicrlw$bi_class=cx
# Combine the main category and subcategory to make a label
databicrlw$combined_category <-databicrlw$bi_class

#databicrlw=databicrlw[-which(is.na(databicrlw$d2015)),]
#plot(databicrlw$y,databitot$y)

##7.2 Pixel level -----------

###7.2.1 clim+reservoirs ----
databipicr=databipic
databipicr$x= databipire$x + databipic$x
databipicr$y= databipire$y + databipic$y

alterclass=data.frame(databipicr$x)
alterclass$class=NA
alterclass$class[which(alterclass[,1]<=(-breaker2))]=1
alterclass$class[which(alterclass[,1]>(-breaker2) & alterclass[,1]<(breaker1))]=2
alterclass$class[which(alterclass[,1]>=(breaker1) & alterclass[,1]<(breaker2))]=3
alterclass$class[which(alterclass[,1]>=breaker2)]=4
c1=alterclass$class

alterclass=data.frame(databipicr$y)
alterclass$class=NA
alterclass$class[which(alterclass[,1]<=(-breaker2))]=1
alterclass$class[which(alterclass[,1]>(-breaker2) & alterclass[,1]<(breaker1))]=2
alterclass$class[which(alterclass[,1]>=(breaker1) & alterclass[,1]<(breaker2))]=3
alterclass$class[which(alterclass[,1]>=breaker2)]=4
c2=alterclass$class

cx=paste(c2,c1,sep="-")
databipicr$bi_class=cx
# Combine the main category and subcategory to make a label
databipicr$combined_category <-databipicr$bi_class

###7.2.2 clim+reservoirs+lu ----
databipicrl=databipic

databipicrl$x=databipilu$x + databipire$x + databipic$x
databipicrl$y=databipilu$y + databipire$y + databipic$y

alterclass=data.frame(databipicrl$x)
alterclass$class=NA
alterclass$class[which(alterclass[,1]<=(-breaker2))]=1
alterclass$class[which(alterclass[,1]>(-breaker2) & alterclass[,1]<(breaker1))]=2
alterclass$class[which(alterclass[,1]>=(breaker1) & alterclass[,1]<(breaker2))]=3
alterclass$class[which(alterclass[,1]>=breaker2)]=4
c1=alterclass$class

alterclass=data.frame(databipicrl$y)
alterclass$class=NA
alterclass$class[which(alterclass[,1]<=(-breaker2))]=1
alterclass$class[which(alterclass[,1]>(-breaker2) & alterclass[,1]<(breaker1))]=2
alterclass$class[which(alterclass[,1]>=(breaker1) & alterclass[,1]<(breaker2))]=3
alterclass$class[which(alterclass[,1]>=breaker2)]=4
c2=alterclass$class
cx=paste(c2,c1,sep="-")
databipicrl$bi_class=cx
# Combine the main category and subcategory to make a label
databipicrl$combined_category <-databipicrl$bi_class

###7.2.3 clim+reservoirs+lu+wd ----
databipicrlw=databipilu

databipicrlw$x=databipilu$x + databipire$x + databipic$x + databipiwd$x
databipicrlw$y=databipilu$y + databipire$y + databipic$y + databipiwd$y

alterclass=data.frame(databipicrlw$x)
alterclass$class=NA
alterclass$class[which(alterclass[,1]<=(-breaker2))]=1
alterclass$class[which(alterclass[,1]>(-breaker2) & alterclass[,1]<(breaker1))]=2
alterclass$class[which(alterclass[,1]>=(breaker1) & alterclass[,1]<(breaker2))]=3
alterclass$class[which(alterclass[,1]>=breaker2)]=4
c1=alterclass$class

alterclass=data.frame(databipicrlw$y)
alterclass$class=NA
alterclass$class[which(alterclass[,1]<=(-breaker2))]=1
alterclass$class[which(alterclass[,1]>(-breaker2) & alterclass[,1]<(breaker1))]=2
alterclass$class[which(alterclass[,1]>=(breaker1) & alterclass[,1]<(breaker2))]=3
alterclass$class[which(alterclass[,1]>=breaker2)]=4
c2=alterclass$class
cx=paste(c2,c1,sep="-")
databipicrlw$bi_class=cx

# Combine the main category and subcategory to make a label
databipicrlw$combined_category <-databipicrlw$bi_class



# 7.3 SANKEY AT HER LEVEL ------------

#Sankey diagram of transfers between classes

databiclim$maxicat="Wetting"
databiclim$maxicat[which(databiclim$combined_category == "1-1" |
                         databiclim$combined_category == "1-2" |
                         databiclim$combined_category == "2-1" )]="Drying"

databiclim$maxicat[which(databiclim$combined_category == "3-1" |
                           databiclim$combined_category == "4-2" |
                           databiclim$combined_category == "4-1" )]="Decelerating"

databiclim$maxicat[which(databiclim$combined_category == "1-3" |
                           databiclim$combined_category == "1-4" |
                           databiclim$combined_category == "2-4" )]="Accelerating"

databiclim$maxicat[which(databiclim$combined_category == "3-3" |
                          databiclim$combined_category == "3-2" |
                          databiclim$combined_category == "2-3" |
                          databiclim$combined_category == "2-2" )]="Stable"

databicr$maxicat="Wetting"
databicr$maxicat[which(databicr$combined_category == "1-1" |
                           databicr$combined_category == "1-2" |
                           databicr$combined_category == "2-1")]="Drying"

databicr$maxicat[which(databicr$combined_category == "3-1" |
                           databicr$combined_category == "4-2" |
                           databicr$combined_category == "4-1" )]="Decelerating"

databicr$maxicat[which(databicr$combined_category == "1-3" |
                           databicr$combined_category == "1-4" |
                           databicr$combined_category == "2-4" )]="Accelerating"

databicr$maxicat[which(databicr$combined_category == "3-3" |
                           databicr$combined_category == "3-2" |
                           databicr$combined_category == "2-3" |
                           databicr$combined_category == "2-2" )]="Stable"

databicrl$maxicat="Wetting"
databicrl$maxicat[which(databicrl$combined_category == "1-1" |
                           databicrl$combined_category == "1-2" |
                           databicrl$combined_category == "2-1" )]="Drying"

databicrl$maxicat[which(databicrl$combined_category == "3-1" |
                           databicrl$combined_category == "4-2" |
                           databicrl$combined_category == "4-1" )]="Decelerating"

databicrl$maxicat[which(databicrl$combined_category == "1-3" |
                           databicrl$combined_category == "1-4" |
                           databicrl$combined_category == "2-4" )]="Accelerating"

databicrl$maxicat[which(databicrl$combined_category == "3-3" |
                           databicrl$combined_category == "3-2" |
                           databicrl$combined_category == "2-3" |
                           databicrl$combined_category == "2-2" )]="Stable"


databicrlw$maxicat="Wetting"
databicrlw$maxicat[which(databicrlw$combined_category == "1-1" |
                           databicrlw$combined_category == "1-2" |
                           databicrlw$combined_category == "2-1" )]="Drying"

databicrlw$maxicat[which(databicrlw$combined_category == "3-1" |
                           databicrlw$combined_category == "4-2" |
                           databicrlw$combined_category == "4-1")]="Decelerating"

databicrlw$maxicat[which(databicrlw$combined_category == "1-3" |
                           databicrlw$combined_category == "1-4" |
                           databicrlw$combined_category == "2-4" )]="Accelerating"

databicrlw$maxicat[which(databicrlw$combined_category == "3-3" |
                           databicrlw$combined_category == "3-2" |
                           databicrlw$combined_category == "2-3" |
                           databicrlw$combined_category == "2-2" )]="Stable"

class1=databiclim$maxicat
class2=databicr$maxicat
class3=databicrl$maxicat
class4=databicrlw$maxicat

ltot=length(databiclim$bi_class)
loscolors=c("Accelerating" = "#174f28","Drying" = "#dd6a29","Stable"="gray60","Wetting" = "#169dd0","Decelerating" = "burlywood")
#8 categories 
# Define the nodes
brk=unique(class1)
library(ggsankey)


d <- data.frame(cbind(class1,class2,class3,class4))
names(d) <- c('OnlyClimate',
              'ClimateReservoirs',
              'ClimateReservoirLanduse',
              'AllDrivers')

df <- d%>%
  make_long(OnlyClimate, 
            ClimateReservoirs,
            ClimateReservoirLanduse,
            AllDrivers)


df$node <- factor(df$node,levels = brk)
df$next_node <- factor(df$next_node,levels = brk)

df$idnext=c(1:length(df$next_node))
reagg <- df%>%
  dplyr::group_by(node,x)%>%  # Here we are grouping the data by node and then we are taking the frequency of it 
  tally()

df2 <- merge(df, 
             reagg, 
             by.x = 'node', 
             by.y = 'node', 
             all.x = FALSE)

df2=df2[which(df2$x.x==df2$x.y),]
df2$np=round(df2$n/ltot*100)
pl <- ggplot(df2, aes(x = x.x,                        
                      next_x = next_x,                                     
                      node = node,
                      next_node = next_node,        
                      fill = factor(node),
                      label= paste0(node,"\n", np,"%")))


###[SPlot] - Sankey diagram at HER level ----

pl <- pl +geom_sankey(flow.alpha = 0.5,          #This Creates the transparency of your node 
                      node.color = "black",     # This is your node color        
                      show.legend = TRUE)        # This determines if you want your legend to show

pl <- pl + geom_sankey_label(size = 3, 
                             color = "black", 
                             fill = "white")


pl <- pl + theme_sankey(base_size = 18) 
pl <- pl + theme(legend.position = 'none')
pl <- pl + theme(axis.title = element_blank(),
                 axis.text.y = element_blank(),
                 axis.ticks = element_blank(),
                 panel.grid = element_blank())

pl=pl + scale_fill_manual(values = loscolors) +
  ggtitle("Change in the Europe's hydrological cycle between the 1950s and the 2010s")

pl


ggsave(paste0("D:/tilloal/Documents/LFRuns_utils/TrendAnalysis/plots/Sankey_changes_catNew.jpg"), pl, width=30, height=20, units=c("cm"),dpi=300) 

mean_lf=c(mean(databitot$d2015,na.rm=T),mean(databiclim$d2015,na.rm=T),mean(databicr$d2015,na.rm=T),mean(databicrl$d2015,na.rm=T))
mean_scen=data.frame()

ggplot() +
  geom_point(data=databitot, aes(x = d2015, y = f2015, fill = combined_category),
             shape = 21, color = "transparent", fill="blue",size = 3,alpha=0.1,show.legend = FALSE) +
  geom_point(data=databiclim, aes(x = d2015, y = f2015, fill = combined_category),
             shape = 21, color = "transparent", fill= "red",alpha=0.1, size = 3,show.legend = FALSE) +
  geom_point(data=databicr, aes(x = d2015, y = f2015, fill = combined_category),
             shape = 21, color = "transparent", fill= "darkgreen",alpha=0.1,size = 3,show.legend = FALSE) +
  geom_point(data=databicrl, aes(x = d2015, y = f2015, fill = combined_category),
             shape = 21, color = "transparent", fill= "purple",alpha=0.1,size = 3,show.legend = FALSE) +

  coord_cartesian(xlim=c(-100,100),ylim=c(-50,50)) +
  
  # Add main quadrant labels
  annotate("text", x = 1.8, y = 40, label = "Wetting", color = "gray12", size = 5, fontface = "bold") +
  annotate("text", x = -1.8, y = 40, label = "Accelerating", color = "gray12", size = 5, fontface = "bold") +
  annotate("text", x = -1.8, y = -40, label = "Drying", color = "gray12", size = 5, fontface = "bold") +
  annotate("text", x = 1.8, y = -40, label = "Decelerating", color = "gray12", size = 5, fontface = "bold") +
  
  # Set axis labels
  labs(x = "Change in drought flows (l/s/km2)", 
       y = "Change in flood flows (l/s/km2)") +
  
  # Customize the theme
  theme(axis.title=element_text(size=tsize),
        panel.background = element_rect(fill = "transparent", colour = "grey1"),
        panel.border = element_rect(linetype = "solid", fill = NA, colour="black"),
        legend.title = element_text(size=tsize),
        legend.text = element_text(size=osize),
        legend.position = "bottom",
        panel.grid.major = element_line(colour = "grey70"),
        panel.grid.minor = element_line(colour = "grey90"),
        legend.key = element_rect(fill = "transparent", colour = "transparent"),
        legend.key.size = unit(.8, "cm"))

##7.3 SANKEY for pixels #################


mean(databipic$y,na.rm=T)
mean(databipicr$y,na.rm=T)
mean(databipicrl$y,na.rm=T)

plot(databipi$x[1:1000]-databipicrlw$x[1:1000],na.rm=T)
median(databipicrlw$y,na.rm=T)

cor.test(databipi$y,databipicrlw$y,na.rm=T)
#Sankey diagram of transers between classes
databipi=databipi[-which(is.na(databipicrl$x)),]
databipic=databipic[-which(is.na(databipicrl$x)),]
databipicr=databipicr[-which(is.na(databipicrl$x)),]
databipicrlw=databipicrlw[-which(is.na(databipicrl$x)),]
databipicrl=databipicrl[-which(is.na(databipicrl$x)),]

databipic$maxicat="Wetting"
databipic$maxicat[which(databipic$bi_class == "1-1" |
                           databipic$bi_class == "1-2" |
                           databipic$bi_class == "2-1" )]="Drying"

databipic$maxicat[which(databipic$bi_class == "3-1" |
                           databipic$bi_class == "4-2" |
                           databipic$bi_class == "4-1" )]="Decelerating"

databipic$maxicat[which(databipic$bi_class == "1-3" |
                           databipic$bi_class == "1-4" |
                           databipic$bi_class == "2-4" )]="Accelerating"

databipic$maxicat[which(databipic$bi_class == "3-3" |
                           databipic$bi_class == "3-2" |
                           databipic$bi_class == "2-3" |
                           databipic$bi_class == "2-2" )]="Stable"

databipicr$maxicat="Wetting"
databipicr$maxicat[which(databipicr$bi_class == "1-1" |
                          databipicr$bi_class == "1-2" |
                          databipicr$bi_class == "2-1")]="Drying"

databipicr$maxicat[which(databipicr$bi_class == "3-1" |
                          databipicr$bi_class == "4-2" |
                          databipicr$bi_class == "4-1" )]="Decelerating"

databipicr$maxicat[which(databipicr$bi_class == "1-3" |
                          databipicr$bi_class == "1-4" |
                          databipicr$bi_class == "2-4" )]="Accelerating"

databipicr$maxicat[which(databipicr$bi_class == "3-3" |
                          databipicr$bi_class == "3-2" |
                          databipicr$bi_class == "2-3" |
                          databipicr$bi_class == "2-2" )]="Stable"

databipicrl$maxicat="Wetting"
databipicrl$maxicat[which(databipicrl$bi_class == "1-1" |
                           databipicrl$bi_class == "1-2" |
                           databipicrl$bi_class == "2-1" )]="Drying"

databipicrl$maxicat[which(databipicrl$bi_class == "3-1" |
                           databipicrl$bi_class == "4-2" |
                           databipicrl$bi_class == "4-1" )]="Decelerating"

databipicrl$maxicat[which(databipicrl$bi_class == "1-3" |
                           databipicrl$bi_class == "1-4" |
                           databipicrl$bi_class == "2-4" )]="Accelerating"

databipicrl$maxicat[which(databipicrl$bi_class == "3-3" |
                           databipicrl$bi_class == "3-2" |
                           databipicrl$bi_class == "2-3" |
                           databipicrl$bi_class == "2-2" )]="Stable"

databipicrlw$maxicat="Wetting"
databipicrlw$maxicat[which(databipicrlw$bi_class == "1-1" |
                            databipicrlw$bi_class == "1-2" |
                            databipicrlw$bi_class == "2-1" )]="Drying"

databipicrlw$maxicat[which(databipicrlw$bi_class == "3-1" |
                            databipicrlw$bi_class == "4-2" |
                            databipicrlw$bi_class == "4-1" )]="Decelerating"

databipicrlw$maxicat[which(databipicrlw$bi_class == "1-3" |
                            databipicrlw$bi_class == "1-4" |
                            databipicrlw$bi_class == "2-4" )]="Accelerating"

databipicrlw$maxicat[which(databipicrlw$bi_class == "3-3" |
                            databipicrlw$bi_class == "3-2" |
                            databipicrlw$bi_class == "2-3" |
                            databipicrlw$bi_class == "2-2" )]="Stable"


databipi$maxicat="Wetting"
databipi$maxicat[which(databipi$bi_class == "1-1" |
                          databipi$bi_class == "1-2" |
                          databipi$bi_class == "2-1" )]="Drying"

databipi$maxicat[which(databipi$bi_class == "3-1" |
                          databipi$bi_class == "4-2" |
                          databipi$bi_class == "4-1")]="Decelerating"

databipi$maxicat[which(databipi$bi_class == "1-3" |
                          databipi$bi_class == "1-4" |
                          databipi$bi_class == "2-4" )]="Accelerating"

databipi$maxicat[which(databipi$bi_class == "3-3" |
                          databipi$bi_class == "3-2" |
                          databipi$bi_class == "2-3" |
                          databipi$bi_class == "2-2" )]="Stable"



###7.3.2 Side analysis -----
#how many pixel are accelerating, deceleration, etch in each biogeoregion?

bgname=c("Alpine","Atlantic","Boreal", "Continental","Mediterranean")
categorix=c("Accelerating","Drying","Stable","Wetting" ,"Decelerating")

#only climate driver
results_df <- data.frame(matrix(0, nrow = length(bgname), ncol = length(categorix)))
rownames(results_df) <- bgname
colnames(results_df) <- categorix

# Loop through each 'neu' in 'bgname'
for (neu in bgname) {
  databipicBG <- databipic[which(databipic$Biogeo_id == neu), ]
  
  # Loop through each 'cat' in 'categorix'
  for (cat in categorix) {
    lol <- length(which(databipicBG$maxicat == cat)) / length(databipicBG$maxicat) * 100
    
    # Fill the data frame with the percentage for the current 'neu' and 'cat'
    results_df[neu, cat] <- round(lol)
  }
}

# Print the results as a table
print(results_df)

#All drivers
# Initialize an empty data frame with 'bgname' as row names and 'categorix' as column names
results_df_all <- data.frame(matrix(0, nrow = length(bgname), ncol = length(categorix)))
rownames(results_df_all) <- bgname
colnames(results_df_all) <- categorix

# Loop through each 'neu' in 'bgname'
for (neu in bgname) {
  databipicBG <- databipi[which(databipi$Biogeo_id == neu), ]
  
  # Loop through each 'cat' in 'categorix'
  for (cat in categorix) {
    lol <- length(which(databipicBG$maxicat == cat)) / length(databipicBG$maxicat) * 100
    
    # Fill the data frame with the percentage for the current 'neu' and 'cat'
    results_df_all[neu, cat] <- round(lol)
  }
}

# Print the results as a table
print(results_df_all)


# Loop through each 'cat' in 'categorix'
results_eu=c()
for (cat in categorix) {
  lol <- length(which(databipic$maxicat == cat)) / length(databipic$maxicat) * 100
  
  # Fill the data frame with the percentage for the current 'neu' and 'cat'
  results_eu <- c(results_eu,round(lol))
}

# Print the results as a table
print(results_eu)


####[SPlot] - Sankey diagram at pixel level ----
class1=databipic$maxicat
class2=databipicr$maxicat
class3=databipicrl$maxicat
class4=databipicrlw$maxicat

ltot=length(databipi$bi_class)

loscolors=c("Accelerating" = "#174f28","Drying" = "#dd6a29","Stable"="gray60","Wetting" = "#169dd0","Decelerating" = "burlywood")
#8 categories 
# Define the nodes
brk=unique(class1)


d <- data.frame(cbind(class1,class2,class3,class4))
names(d) <- c('OnlyClimate',
              'ClimateReservoirs',
              'ClimateReservoirLanduse',
              'AllDrivers')

df <- d%>%
  make_long(OnlyClimate, 
            ClimateReservoirs,
            ClimateReservoirLanduse,
            AllDrivers)


df$node <- factor(df$node,levels = brk)
df$next_node <- factor(df$next_node,levels = brk)

df$idnext=c(1:length(df$next_node))
reagg <- df%>%
  dplyr::group_by(node,x)%>%  # Here we are grouping the data by node and then we are taking the frequency of it 
  tally()

df2 <- merge(df, 
             reagg, 
             by.x = 'node', 
             by.y = 'node', 
             all.x = FALSE)
df2=df2[which(df2$x.x==df2$x.y),]

df2$np=round(df2$n/ltot*100)
pl <- ggplot(df2, aes(x = x.x,                        
                      next_x = next_x,                                     
                      node = node,
                      next_node = next_node,        
                      fill = factor(node),
                      label= paste0(node,"\n", np,"%")))



pl <- pl +geom_sankey(flow.alpha = 0.5,          #This Creates the transparency of your node 
                      node.color = "black",     # This is your node color        
                      show.legend = TRUE)        # This determines if you want your legend to show

pl <- pl + geom_sankey_label(size = 3, 
                             color = "black", 
                             fill = "white")


pl <- pl + theme_sankey(base_size = 18) 
pl <- pl + theme(legend.position = 'none')
pl <- pl + theme(axis.title = element_blank(),
                 axis.text.y = element_blank(),
                 axis.ticks = element_blank(),
                 panel.grid = element_blank())

pl=pl + scale_fill_manual(values = loscolors) +
  ggtitle("Change in the Europe's hydrological cycle between the 1950s and the 2010s")

pl


ggsave(paste0("D:/tilloal/Documents/LFRuns_utils/TrendAnalysis/plots/Sankey_changes_pixels5.jpg"), pl, width=30, height=10, units=c("cm"),dpi=300) 
