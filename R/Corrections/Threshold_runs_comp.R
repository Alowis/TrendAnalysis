#Analysing the impact of thresholds


setwd(dirname(rstudioapi::getSourceEditorContext()$path))
getwd()
source("D:/tilloal/Documents/LFRuns_utils/TrendAnalysis/R/functions_trends.R")

#2 Pre-loaded results -----------
#Set data directory
hydroDir<-("D:/tilloal/Documents/LFRuns_utils/ChangingHydroExtremes/data")
hydroDir<-("D:/tilloal/Documents/LFRuns_utils/ChangingHydroExtremes/data")

#outlets file outf
if (!exists("outf")){
  outf=c()
  for( Nsq in 1:88){
    print(Nsq)
    rspace= read.csv(paste0(hydroDir,"/subspace_efas.csv"))
    rspace=rspace[,-1]
    nrspace=rspace[Nsq,]
    outletname="GeoData/efas_rnet_100km_01min"
    
    outhybas=outletopen(hydroDir,outletname,nrspace)
    Idstart=as.numeric(Nsq)*10000
    Idstart2=as.numeric(Nsq)*100000
    if (length(outhybas$outlets)>0){
      outhybas$outlets=seq((Idstart+1),(Idstart+length(outhybas$outlets)))
      outhybas$outl2=seq((Idstart2+1),(Idstart2+length(outhybas$outlets)))
      outhybas$latlong=paste(round(outhybas$Var1,4),round(outhybas$Var2,4),sep=" ")
      outhloc=outhybas
      outf=rbind(outf,outhloc)
    }
  }
}

##2.1 Spatial data for catchments ----

### Hybas07 ----
Catchmentrivers7=read.csv(paste0(hydroDir,"/GeoData/HYBAS07/from_hybas_eu_onlyid.csv"),encoding = "UTF-8", header = T, stringsAsFactors = F)
hybas07 <- read_sf(dsn = paste0(hydroDir,"/GeoData/HYBAS07/hybas_eu_lev07_v1c.shp"))
hybasf7=fortify(hybas07)
Catamere07=inner_join(hybasf7,Catchmentrivers7,by= "HYBAS_ID")
Catamere07$llcoord=paste(round(Catamere07$POINT_X,4),round(Catamere07$POINT_Y,4),sep=" ")
cst7=right_join(Catamere07,outf,by= c("llcoord"="latlong"))
GNF=cst7
length(unique(GNF$HYBAS_ID))
st_geometry(GNF)=NULL
rm(Catamere07)
outlethybas07="/GeoData/HYBAS07/outletsv8_hybas07_01min"
outhybas07=outletopen(hydroDir,outlethybas07)

#matching outlets with pixel Ids
outhybas07$latlong=paste(round(outhybas07$Var1,4),round(outhybas07$Var2,4),sep=" ")
mhy=match(outhybas07$latlong,outf$latlong)
outhybas07$outID=outf$outl2[mhy]

### European Biogeo regions ----
biogeo <- read_sf(dsn = paste0(hydroDir,"/GeoData/eea_3035_biogeo-regions_2016/BiogeoRegions2016_wag84.shp"))
biogeof=fortify(biogeo)
st_geometry(biogeof)<-NULL
biogeoregions=raster( paste0(hydroDir,"/GeoData/eea_3035_biogeo-regions_2016/Biogeo_rasterized_wsg84.tif"))
Gbiogeoregions=as.data.frame(biogeoregions,xy=T)
biogeomatch=inner_join(biogeof,Gbiogeoregions,by= c("PK_UID"="Biogeo_rasterized_wsg84"))
biogeomatch$latlong=paste(round(biogeomatch$x,4),round(biogeomatch$y,4),sep=" ")
biogeo_rivers=right_join(biogeomatch,outf, by="latlong")

### HydroRegions ----
GridHR=raster( paste0(hydroDir,"/GeoData/HER/HydroRegions_raster_WGS84.tif"))
GHR=as.data.frame(GridHR,xy=T)
GHR2=GHR
GHR2$llcoord=paste(round(GHR2$x,4),round(GHR2$y,4),sep=" ")
GHR=GHR[which(!is.na(GHR[,3])),]
GHR$llcoord=paste(round(GHR$x,4),round(GHR$y,4),sep=" ")
GHR_riv=inner_join(GHR,outf,by= c("llcoord"="latlong"))
GHshpp <- read_sf(dsn =paste0(hydroDir,"/GeoData/HER/her_all_adjusted.shp"))
HydroRsf=fortify(GHshpp)

### Plot parameters ----
palet2=c(hcl.colors(9, palette = "Blues", alpha = NULL, rev = TRUE, fixup = TRUE))
outletname="/GeoData/efas_rnet_100km_01min"
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

##2.2 Loading saved results in .Rdata ---------------------------

###load UpArea -----
#load upstream area
# main_path = 'D:/tilloal/Documents/06_Floodrivers/'
# valid_path = paste0(main_path,'DataPaper/')
outletname="/GeoData/upArea_European_01min.nc"
#dir=valid_path
outf$idlalo=paste(outf$idlo, outf$idla, sep=" ")
UpArea=UpAopen(hydroDir,outletname,outf)
head(UpArea)

hist(log(UpArea$upa))
n1=round(length(UpArea$upa[which(UpArea$upa<45000)])/length(UpArea$upa)*100,1)
n2=round(length(UpArea$upa[which(UpArea$upa>=45000 & UpArea$upa<100000)])/length(UpArea$upa)*100,1)
n3=round(length(UpArea$upa[which(UpArea$upa>=100000)])/length(UpArea$upa)*100,1)

#' Generate Minor Logarithmic Breaks for Plotting
#'
#' @description
#' The `log10_minor_break` function generates a sequence of minor logarithmic breaks
#' for use in plotting functions, particularly useful when customizing the minor
#' breaks on a logarithmic scale.
#'
#' @param ... Additional arguments (currently not used).
#'
#' @return A function that takes a numeric vector `x` and returns a numeric vector
#'   of minor breaks in the original data scale (not log-transformed). These breaks
#'   are calculated between the major logarithmic breaks of `x`.
#'
#' @examples
#' # Generate a set of logarithmic minor breaks for plotting
#' x_data <- c(0.1, 1, 10, 100, 1000)
#' minor_breaks_func <- log10_minor_break()
#' minor_breaks <- minor_breaks_func(x_data)
#'
#' # Use the minor breaks in a plot (for example, with ggplot2)
#' # ggplot(data, aes(x, y)) +
#' #   geom_point() +
#' #   scale_x_log10(minor_breaks = minor_breaks)
#'
#' @export
#'
#' @seealso
#' Plotting functions in packages like `ggplot2`, which allow for customization
#' of minor breaks on a logarithmic scale.
#'
log10_minor_break = function (...){
  function(x) {
    minx         = floor(min(log10(x), na.rm=T))-1;
    maxx         = ceiling(max(log10(x), na.rm=T))+1;
    n_major      = maxx-minx+1;
    major_breaks = seq(minx, maxx, by=1)
    minor_breaks = 
      rep(log10(seq(1, 9, by=1)), times = n_major)+
      rep(major_breaks, each = 9)
    return(10^(minor_breaks))
  }
}

###load UpArea -----
#load upstream area
# main_path = 'D:/tilloal/Documents/06_Floodrivers/'
# valid_path = paste0(main_path,'DataPaper/')
outletname="/GeoData/upArea_European_01min.nc"
#dir=valid_path
outf$idlalo=paste(outf$idlo, outf$idla, sep=" ")
UpArea=UpAopen(hydroDir,outletname,outf)
head(UpArea)


p<-ggplot(UpArea, aes(x=upa)) + 
  geom_histogram(color="steelblue", fill="slategray1",bins=15,alpha=0.9,lwd=1)+
  scale_y_continuous(breaks=seq(0,100000, by=10000),name="Number of pixels")+
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
  annotate("label", x=300000, y=50000, label= paste0("n_small = ",n1,"%\n n_medium = ",n2,"% \n n_large = ",n3,"%"),size=6)

p
ggsave("D:/tilloal/Documents/01_Projects/FlooDrough/Figures/Revisions/histo_pixels_revisions.jpg", p, width=20, height=15, units=c("cm"),dpi=1500)


#Loading fitting results from the 4 runs and for all 282 000 river pixels in
# in the domain. Requires at least 20 GB of free RAM.

###load historical run -----
haz="Flood"
if (haz == "Drought") namefile="Drought.nonfrost.Histo"
if (haz == "Flood") namefile="flood.year.Histo"

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

epsilonH=ParamsflH$epsilonGPD[which(ParamsflH$Year==1955)]

epsilon=ParamsflH$epsilonGPD[which(ParamsflH$Year==1955)]
sigma=ParamsflH$sigmaGPD[which(ParamsflH$Year==1955)]
threshold=ParamsflH$thresholdGPD[which(ParamsflH$Year==1955)]
epsilonStdErrGPD=ParamsflH$epsilonStdErrGPD[which(ParamsflH$Year==1955)]
sigmaStdErrGPD=ParamsflH$sigmaStdErrGPD[which(ParamsflH$Year==1955)]
thresholdStdErrGPD=ParamsflH$thresholdStdErrGPD[which(ParamsflH$Year==1955)]
nPeaks=ParamsflH$nPeaks[which(ParamsflH$Year==1955)]

RPgoal=100

calcGPDReturnLevel_Single <- function(epsilon, sigma, threshold, nPeaks, sampleTimeHorizon, returnPeriod) {
  
  # Calculate expected number of events in T years
  XX <- (nPeaks / sampleTimeHorizon) * returnPeriod
  
  # Initialize result vector with NAs
  returnLevel <- rep(NA_real_, length(sigma))
  
  # 1. Handle GPD Case: epsilon is NOT NA and NOT zero
  idx_gpd <- !is.na(epsilon) & abs(epsilon) > 1e-9
  if (any(idx_gpd)) {
    returnLevel[idx_gpd] <- threshold[idx_gpd] + 
      (sigma[idx_gpd] / epsilon[idx_gpd]) * (XX[idx_gpd]^epsilon[idx_gpd] - 1)
  }
  
  # 2. Handle Gumbel Case: epsilon is NOT NA and IS zero
  idx_gum <- !is.na(epsilon) & abs(epsilon) <= 1e-9
  if (any(idx_gum)) {
    returnLevel[idx_gum] <- threshold[idx_gum] + 
      sigma[idx_gum] * log(XX[idx_gum])
  }
  
  # Note: Rows where epsilon/sigma/threshold were originally NA will remain NA
  return(returnLevel)
}


RL100s=calcGPDReturnLevel_Single(epsilon, sigma, threshold,
                                   nPeaks, sampleTimeHorizon=70, returnPeriod = RPgoal)

###load Socio-CF run -----
if (haz == "Drought") namefile="Drought.nonfrost.SocCF"
if (haz == "Flood") namefile="Flood.year.socCF"

haz="Flood"
load(file=paste0(hydroDir,"/",haz,"/params.",namefile,".Rdata"))
load(file=paste0(hydroDir,"/",haz,"/RL100.",namefile,".Rdata"))
load(file=paste0(hydroDir,"/",haz,"/peaks.",namefile,".Rdata"))

RLGPDflSCF=RLGPDfl
Paramsfl=Paramsfl[,-c(4:9,17)]
ParamsflSCF=data.table(Paramsfl)
PeakSCF=Peaksave
rm(Paramsfl,RLGPDfl)
gc()

epsilonS=ParamsflSCF$epsilonGPD[which(ParamsflSCF$Year==1955)]

###load results from Res+WU CF run -----
if (haz == "Drought") namefile="Drought.nonfrost.RWCF"
if (haz == "Flood") namefile="flood.year.RWCF"
load(file=paste0(hydroDir,"/",haz,"/params.",namefile,".Rdata"))
load(file=paste0(hydroDir,"/",haz,"/RL100.",namefile,".Rdata"))
load(file=paste0(hydroDir,"/",haz,"/peaks.",namefile,".Rdata"))
RLGPDflRWCF=RLGPDfl
Paramsfl=Paramsfl[,-c(4:9,17)]
ParamsflRWCF=data.table(Paramsfl)
PeakRWCF=Peaksave
rm(Paramsfl,RLGPDfl)
gc()

epsilonRW=ParamsflRWCF$epsilonGPD[which(ParamsflRWCF$Year==1955)]

###load results from Water CF run -----
if (haz == "Drought") namefile="Drought.nonfrost.WCF"
if (haz == "Flood") namefile="flood.year.WCF"
load(file=paste0(hydroDir,"/",haz,"/params.",namefile,".Rdata"))
load(file=paste0(hydroDir,"/",haz,"/RL100.",namefile,".Rdata"))
load(file=paste0(hydroDir,"/",haz,"/peaks.",namefile,".Rdata"))

RLGPDflWCF=RLGPDfl
Paramsfl=Paramsfl[,-c(4:9,17)]
ParamsflWCF=data.table(Paramsfl)
PeakWCF=Peaksave
rm(Paramsfl,RLGPDfl)
gc()

epsilonW=ParamsflWCF$epsilonGPD[which(ParamsflWCF$Year==1955)]
rm(catmap)
gc()

#compare shape parameters

EpsilonTab=data.frame(catchment=ParamsflWCF$catchment[which(ParamsflWCF$Year==1955)],epsilonH,epsilonS,epsilonRW,epsilonW)
EpsilonTab=inner_join(EpsilonTab,UpArea,by=c("catchment"="outl2"))

EpsilonTab$Wuchange=EpsilonTab$epsilonW-EpsilonTab$epsilonH
EpsilonTab$SEchange=EpsilonTab$epsilonS-EpsilonTab$epsilonH
EpsilonTab$Luchange=EpsilonTab$epsilonS-EpsilonTab$epsilonRW
quantile(abs(EpsilonTab$Luchange),0.5)
points <- st_as_sf(EpsilonTab, coords = c("Var1.x", "Var2.x"), crs = 4326)
points <- st_transform(points, crs = 3035)

# palet=c(hcl.colors(11, palette = "YlOrRd", alpha = NULL, rev = T, fixup = TRUE))

palet=c(hcl.colors(11, palette = "RdYlBu", alpha = NULL, rev = T, fixup = TRUE))
cazzo=ggplot(basemap) +
  geom_sf(fill="gray95",color="gray10",size=0.5)+
  geom_sf(data=points,aes(col=Luchange,geometry=geometry,size=upa),alpha=.9,stroke=0,shape=15)+
  coord_sf(xlim = c(min(nco[,1]),max(nco[,1])), ylim = c(min(nco[,2]),max(nco[,2])))+
  # scale_colour_manual(values = colIR, name="IR", labels=c("0"="perennial","1"="casi-perennial","2"="IRES")) +
  scale_size(range = c(0.08, 0.4), trans="sqrt",name= expression(paste("Upstream area ", (km^2),
                                                                       sep = " ")),
             breaks=c(101,1000,10000,100000,500000), labels=c("100","1000", "10 000", "100 000", "500 000"),
             guide = "none")+
  scale_color_gradientn(
    colors=palet,
    breaks=seq(-0.1,0.1,by=0.02), limits=c(-0.1,0.1),
    oob = scales::squish)   +
  labs(x="Longitude", y = "Latitude")+
  # guides(colour = guide_legend(override.aes = list(size = 10)))+
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

ggsave(cazzo,filename=paste0("D:/tilloal/Documents/LFRuns_utils/TrendAnalysis/plots/maps_shapes_Lu",haz,".jpg"), width=23, height=20, units=c("cm"),dpi=1000) 


#compare shape parameters

SigmaTab=data.frame(catchment=ParamsflWCF$catchment[which(ParamsflWCF$Year==1955)],
                    SigmaH=ParamsflH$sigmaGPD[which(ParamsflH$Year==1955)],
                    SigmaW=ParamsflWCF$sigmaGPD[which(ParamsflH$Year==1955)],
                    SigmaRW=ParamsflRWCF$sigmaGPD[which(ParamsflH$Year==1955)],
                    SigmaS=ParamsflSCF$sigmaGPD[which(ParamsflH$Year==1955)])
SigmaTab=inner_join(SigmaTab,UpArea,by=c("catchment"="outl2"))

SigmaTab$Wuchange=SigmaTab$SigmaW-SigmaTab$SigmaH
points <- st_as_sf(SigmaTab, coords = c("Var1.x", "Var2.x"), crs = 4326)
points <- st_transform(points, crs = 3035)

# palet=c(hcl.colors(11, palette = "YlOrRd", alpha = NULL, rev = T, fixup = TRUE))

palet=c(hcl.colors(11, palette = "RdYlBu", alpha = NULL, rev = T, fixup = TRUE))
cazzo=ggplot(basemap) +
  geom_sf(fill="gray95",color="gray10",size=0.5)+
  geom_sf(data=points,aes(col=Wuchange,geometry=geometry,size=upa),alpha=.9,stroke=0,shape=15)+
  coord_sf(xlim = c(min(nco[,1]),max(nco[,1])), ylim = c(min(nco[,2]),max(nco[,2])))+
  # scale_colour_manual(values = colIR, name="IR", labels=c("0"="perennial","1"="casi-perennial","2"="IRES")) +
  scale_size(range = c(0.08, 0.4), trans="sqrt",name= expression(paste("Upstream area ", (km^2),
                                                                       sep = " ")),
             breaks=c(101,1000,10000,100000,500000), labels=c("100","1000", "10 000", "100 000", "500 000"),
             guide = "none")+
  scale_color_gradientn(
    colors=palet,
    breaks=seq(-0.1,0.1,by=0.02), limits=c(-0.1,0.1),
    oob = scales::squish)   +
  labs(x="Longitude", y = "Latitude")+
  # guides(colour = guide_legend(override.aes = list(size = 10)))+
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

ggsave(cazzo,filename=paste0("D:/tilloal/Documents/LFRuns_utils/TrendAnalysis/plots/maps_scale_wu",haz,".jpg"), width=23, height=20, units=c("cm"),dpi=1000) 







#compare threshold values for different runs

head(PeakH)

library(dplyr)

# ── 1. Add inter-event time (days between consecutive peaks per catchment) ──
PeakH2 <- PeakH %>%
  arrange(catch, time) %>%
  group_by(catch) %>%
  mutate(
    time        = as.POSIXct(time),
    dt_days     = as.numeric(difftime(time, lag(time), units = "days"))
  ) %>%
  ungroup()

# ── 3. Ljung-Box test per catchment ───────────────────────────────────
lb_results <- PeakH2 %>%
  arrange(catch, time) %>%
  group_by(catch) %>%
  summarise(
    n_peaks  = n(),
    lb_stat  = Box.test(value, lag = min(10, floor(n()/5)),
                        type = "Ljung-Box")$statistic,
    lb_pval  = Box.test(value, lag = min(10, floor(n()/5)),
                        type = "Ljung-Box")$p.value,
    .groups  = "drop"
  ) %>%
  mutate(autocorrelated = lb_pval < 0.05)

# how many catchments show significant autocorrelation?
table(lb_results$autocorrelated)

# ── 4. Lag-1 correlation between consecutive peak values ──────────────
lag1_results <- PeakH2 %>%
  arrange(catch, time) %>%
  group_by(catch) %>%
  summarise(
    lag1_cor = cor(value[-n()], value[-1], use = "complete.obs"),
    .groups  = "drop"
  )

hist(lag1_results$lag1_cor, breaks = 30,
     main = "Lag-1 correlation across catchments",
     xlab = "Pearson r")
abline(v = 0, col = "red", lty = 2)

# ── 5. Check independence of inter-event times (exponential = Poisson process) ──
# For a Poisson process, inter-arrival times should be exponential and independent
iet_results <- PeakH %>%
  arrange(catch, time) %>%
  group_by(catch) %>%
  summarise(
    iet_lb_pval = tryCatch(
      Box.test(diff(as.numeric(time)), lag = 5,
               type = "Ljung-Box")$p.value,
      error = function(e) NA_real_
    ),
    .groups = "drop"
  )


library(randtests)

PeakH3=PeakH2[c(1:10000),]


# tIDstart and tIDend define the event window — check no overlap between consecutive events
overlap_check <- PeakH2 %>%
  arrange(catch, time) %>%
  group_by(catch) %>%
  mutate(
    prev_end    = lag(tIDend),
    overlapping = tIDstart < prev_end   # TRUE = overlap with previous event
  ) %>%
  ungroup()

# any overlaps?
overlap_check %>% 
  filter(overlapping, !is.na(overlapping)) %>%
  nrow()

# per catchment summary
final_check=overlap_check %>%
  group_by(catch) %>%
  summarise(n_overlaps = sum(overlapping, na.rm = TRUE)) 


cat("Catchments with non-overlapping windows:   ",
    sum(final_check$n_overlaps, na.rm = TRUE), "/", nrow(final_check), "\n")


#Now plot thresholds

###[Plot] Figure S. - Standard deviation of thresholds for all pixels ----
head(ParamsflH)

threshold_stats <- ParamsflH %>%
  group_by(catchment) %>%
  summarise(
    mean_threshold = mean(thresholdGPD, na.rm = TRUE),
    sd_threshold   = sd(thresholdGPD,   na.rm = TRUE),
    .groups = "drop"
  )

threshold_statsSCF <- ParamsflSCF %>%
  group_by(catchment) %>%
  summarise(
    mean_threshold = mean(thresholdGPD, na.rm = TRUE),
    sd_threshold   = sd(thresholdGPD,   na.rm = TRUE),
    .groups = "drop"
  )
threshold_statsRWCF <- ParamsflRWCF%>%
  group_by(catchment) %>%
  summarise(
    mean_threshold = mean(thresholdGPD, na.rm = TRUE),
    sd_threshold   = sd(thresholdGPD,   na.rm = TRUE),
    .groups = "drop"
  )
threshold_statsWCF <- ParamsflWCF %>%
  group_by(catchment) %>%
  summarise(
    mean_threshold = mean(thresholdGPD, na.rm = TRUE),
    sd_threshold   = sd(thresholdGPD,   na.rm = TRUE),
    .groups = "drop"
  )


library(dplyr)


library(data.table)

# 1. Convert to data.table (in-place to save memory)
setDT(PeakH); setDT(PeakRWCF); setDT(PeakSCF); setDT(PeakWCF)

# 2. Create a "Master List" of unique peaks from File 1
# This will be our reference point for the windows
master_peaks <- PeakH[, .(catch, time)] %>% unique()
master_peaks[, `:=`(win_min = time - 3, win_max = time + 3)]

# 3. Create unique lookup tables for the other files
p2 <- PeakRWCF[, .(catch, time_p2 = time)] %>% unique()
p3 <- PeakSCF[, .(catch, time_p3 = time)] %>% unique()
p4 <- PeakWCF[, .(catch, time_p4 = time)] %>% unique()

# 4. Perform Sequential Non-Equi Joins
# This looks for matches within the +/- 3 day window for each file
shared_peaks <- master_peaks[p2, on = .(catch, win_min <= time_p2, win_max >= time_p2), nomatch = 0]
# We use unique() here to avoid double-counting if one p1 peak matches two p2 peaks
shared_peaks <- unique(shared_peaks[, .(catch, time, win_min, win_max)])

shared_peaks <- shared_peaks[p3, on = .(catch, win_min <= time_p3, win_max >= time_p3), nomatch = 0]
shared_peaks <- unique(shared_peaks[, .(catch, time, win_min, win_max)])

shared_peaks <- shared_peaks[p4, on = .(catch, win_min <= time_p4, win_max >= time_p4), nomatch = 0]
shared_peaks <- unique(shared_peaks[, .(catch, time)])

# 5. Calculate the Final Comparison
# Total peaks in File 1 vs Peaks shared by all 4 within the window
final_comparison <- PeakH[, .(total_peaks_f1 = uniqueN(time)), by = catch]

# Count how many of those f1 peaks survived the joins
shared_counts <- shared_peaks[, .(shared_count = .N), by = catch]

result <- merge(final_comparison, shared_counts, by = "catch", all.x = TRUE)
result[is.na(shared_count), shared_count := 0]
result[, percent_match := (shared_count / total_peaks_f1) * 100]

print(result)

# 3. Create unique lookup tables for the other files
p1 <- PeakH[, .(catch, time_p1 = time)] %>% unique()
p2 <- PeakRWCF[, .(catch, time_p2 = time)] %>% unique()
p3 <- PeakWCF[, .(catch, time_p3 = time)] %>% unique()
p4 <- PeakSCF[, .(catch, time_p4 = time)] %>% unique()
# 1. Prepare each DT with its own window columns

library(data.table)

# Function to clean and unique each run
prep_exact <- function(dt, new_name) {
  setDT(dt)
  # Force columns to be simple vectors to avoid 'list' errors
  clean_dt <- unique(data.table(
    catch = unlist(dt$catch), 
    time = unlist(dt$time)
  ))
  return(clean_dt)
}

# Apply to your 4 dataframes with your new names
dt_list <- list(
  Histo = prep_exact(p1), 
  CF3 = prep_exact(p3), 
  CF2 = prep_exact(p2), 
  CF1 = prep_exact(p4)
)


run_names <- names(dt_list)
sim_matrix_exact <- matrix(1, nrow = 4, ncol = 4, dimnames = list(run_names, run_names))

for (i in 1:3) {
  for (j in (i+1):4) {
    A <- dt_list[[i]]
    B <- dt_list[[j]]
    
    # Exact match join
    common_count <- nrow(merge(A, B, by = c("catch", "time")))
    
    # Union (Total unique peaks across both)
    total_unique <- nrow(unique(rbind(A, B)))
    
    sim_score <- common_count / total_unique
    sim_matrix_exact[i, j] <- sim_score
    sim_matrix_exact[j, i] <- sim_score
  }
}

print(sim_matrix_exact)

library(corrplot)
my_palette <- colorRampPalette(c("#D73027", "white", "#4575B4"))(200)
corrplot(sim_matrix_exact, type = "lower", method = "color", 
         tl.col = "black",is.corr=F,
         # --- Adding the Values ---
         addCoef.col = "black",     # Color of the text inside boxes
         number.digits=2,
         number.cex = 0.7,          # Font size of the text
         col.lim=c(0.8,1),
                 # Text label color (axis labels)
         tl.srt = 45 
         )# Rotate axis labels for readability)




Peakspoints=inner_join(result,UpArea,by=c("catch"="outl2"))




#colIR=c("0"="royalblueResults#colIR=c("0"="royalblue","1"="lightblue","2"="orangered","3"="tomato","4"="purple")
points <- st_as_sf(Peakspoints, coords = c("Var1.x", "Var2.x"), crs = 4326)
points <- st_transform(points, crs = 3035)

#palet=c(hcl.colors(11, palette = "YlOrRd", alpha = NULL, rev = T, fixup = TRUE))

pp<-ggplot(Peakspoints, aes(x=percent_match)) + 
  geom_histogram(color="steelblue", fill="slategray1",bins=100,alpha=0.9,lwd=1)+
  scale_y_continuous(breaks=seq(0,200000, by=10000),name="Number of pixels")+
  scale_x_continuous(breaks=seq(0,100, by=10),name="number of shared peaks between the four runs")+
  # scale_x_sqrt(name=expression(paste("Upstream area ", (km^2),sep = " ")),
  #               breaks=c(100,1000,10000,100000), minor_breaks = log10_minor_break(),
  #               labels=c("100","1 000","10 000","100 000")) +
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
        legend.key.size = unit(.8, "cm"))
#annotate("label", x=300000, y=50000, label= paste0("n_small = ",n1,"%\n n_medium = ",n2,"% \n n_large = ",n3,"%"),size=6)

pp
ggsave(paste0("D:/tilloal/Documents/01_Projects/FlooDrough/Figures/Revisions/histo_shared_peaks_",haz,".jpg"), pp, width=20, height=15, units=c("cm"),dpi=1500)
90/282
length(which(result$percent_match>90))/length(result$percent_match)
palet=c(hcl.colors(11, palette = "RdYlBu", alpha = NULL, rev = F, fixup = TRUE))
peaksh=ggplot(basemap) +
  geom_sf(fill="gray95",color="gray10",size=0.5)+
  geom_sf(data=points,aes(col=percent_match,geometry=geometry,size=upa),alpha=.9,stroke=0,shape=15)+
  coord_sf(xlim = c(min(nco[,1]),max(nco[,1])), ylim = c(min(nco[,2]),max(nco[,2])))+
  # scale_colour_manual(values = colIR, name="IR", labels=c("0"="perennial","1"="casi-perennial","2"="IRES")) +
  scale_size(range = c(0.08, 0.4), trans="sqrt",name= expression(paste("Upstream area ", (km^2),
                                                                       sep = " ")),
             breaks=c(101,1000,10000,100000,500000), labels=c("100","1000", "10 000", "100 000", "500 000"),
             guide = "none")+
  scale_color_gradientn(
    colors=palet,
    breaks=seq(0,100,by=10), limits=c(0,100),
    oob = scales::squish)   +
  labs(x="Longitude", y = "Latitude")+
  # guides(colour = guide_legend(override.aes = list(size = 10)))+
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

ggsave(peaksh,filename=paste0("D:/tilloal/Documents/LFRuns_utils/TrendAnalysis/plots/map_shared_peaks",haz,".jpg"), width=23, height=20, units=c("cm"),dpi=1000) 




ParamsflH[which(ParamsflH$catchment==4400031),]
PeakH[which(PeakH$catch==4400031),]
threshold_stats=ParamsflH[which(ParamsflH$Year==1955),]
threshold_statsSCF=ParamsflSCF[which(ParamsflSCF$Year==1955),]
threshold_statsRWCF=ParamsflRWCF[which(ParamsflRWCF$Year==1955),]
threshold_statsWCF=ParamsflWCF[which(ParamsflWCF$Year==1955),]

final_threshold=data.frame(catchment=threshold_stats$catchment,
                      means=rowMeans(cbind(-threshold_stats$thresholdGPD,
                      -threshold_statsSCF$thresholdGPD,
                      -threshold_statsWCF$thresholdGPD,
                      -threshold_statsRWCF$thresholdGPD)),
                      sds=rowSds(cbind(threshold_stats$thresholdGPD,
                           threshold_statsSCF$thresholdGPD,
                           threshold_statsWCF$thresholdGPD,
                           threshold_statsRWCF$thresholdGPD)))

final_threshold$means[which(final_threshold$means<0)]=NA

ranges=rowRanges(cbind(threshold_stats$thresholdGPD,
                 threshold_statsSCF$thresholdGPD,
                 threshold_statsWCF$thresholdGPD,
                 threshold_statsRWCF$thresholdGPD))
row_diffs=ranges[, 2] - ranges[, 1]
final_threshold$rel_range <- row_diffs / final_threshold$means*100

final_threshold$cv=final_threshold$sds / final_threshold$means * 100

final_npeaks=data.frame(catchment=threshold_stats$catchment,
                           means=rowMeans(cbind(threshold_stats$nPeaks,
                                                threshold_statsSCF$nPeaks,
                                                threshold_statsWCF$nPeaks,
                                                threshold_statsRWCF$nPeaks)),
                           sds=rowSds(cbind(threshold_stats$nPeaks,
                                            threshold_statsSCF$nPeaks,
                                            threshold_statsWCF$nPeaks,
                                            threshold_statsRWCF$nPeaks)))

ranges=rowRanges(cbind(threshold_stats$nPeaks,
                       threshold_statsSCF$nPeaks,
                       threshold_statsWCF$nPeaks,
                       threshold_statsRWCF$nPeaks))
row_diffs=ranges[, 2] - ranges[, 1]
final_npeaks$rel_range <- row_diffs / final_npeaks$means*100

final_npeaks$cv=final_npeaks$sds / final_npeaks$means * 100
length(which(final_npeaks$cv<=10))/length(final_npeaks$cv)
length(which(final_npeaks$rel_range<=10))/length(final_npeaks$cv)
hist(final_npeaks$cv)
#final_npeaks$cv[which(final_npeaks$means<0)]=NA
hist(threshold_stats$nPeaks,breaks=1000,xlim=c(0,100))
threshold_stats[which.max(final_npeaks$cv),]
threshold_statsSCF[which.max(final_npeaks$cv),]
threshold_stats[which.max(final_npeaks$cv),]
threshold_stats[which.max(final_npeaks$cv),]
hist(final_npeaks$cv,breaks=10000,xlim=c(0,100),ylim=c(0,10))


Thpoints=inner_join(final_npeaks,UpArea,by=c("catchment"="outl2"))
head(Thpoints)
#colIR=c("0"="royalblue","1"="lightblue","2"="orangered","3"="tomato","4"="purple")
points <- st_as_sf(Thpoints, coords = c("Var1.x", "Var2.x"), crs = 4326)
points <- st_transform(points, crs = 3035)

# palet=c(hcl.colors(11, palette = "YlOrRd", alpha = NULL, rev = T, fixup = TRUE))

psds<-ggplot(Thpoints, aes(x=rel_range)) + 
  geom_histogram(color="steelblue", fill="slategray1",bins=51,alpha=0.9,lwd=1)+
  scale_y_continuous(breaks=seq(0,200000, by=10000),name="Number of pixels")+
  scale_x_sqrt(breaks=seq(0,100, by=10),limits=c(-0.1,50),name="relative range of npeaks")+
  # scale_x_sqrt(name=expression(paste("Upstream area ", (km^2),sep = " ")),
  #               breaks=c(100,1000,10000,100000), minor_breaks = log10_minor_break(),
  #               labels=c("100","1 000","10 000","100 000")) +
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
        legend.key.size = unit(.8, "cm"))
#annotate("label", x=300000, y=50000, label= paste0("n_small = ",n1,"%\n n_medium = ",n2,"% \n n_large = ",n3,"%"),size=6)

psds
ggsave(paste0("D:/tilloal/Documents/01_Projects/FlooDrough/Figures/Revisions/histo_rrnpy",haz,".jpg"), psds, width=20, height=15, units=c("cm"),dpi=1500)

palet=c(hcl.colors(11, palette = "RdYlBu", alpha = NULL, rev = T, fixup = TRUE))
cazzo=ggplot(basemap) +
  geom_sf(fill="gray95",color="gray10",size=0.5)+
  geom_sf(data=points,aes(col=cv,geometry=geometry,size=upa),alpha=.9,stroke=0,shape=15)+
  coord_sf(xlim = c(min(nco[,1]),max(nco[,1])), ylim = c(min(nco[,2]),max(nco[,2])))+
  # scale_colour_manual(values = colIR, name="IR", labels=c("0"="perennial","1"="casi-perennial","2"="IRES")) +
  scale_size(range = c(0.08, 0.4), trans="sqrt",name= expression(paste("Upstream area ", (km^2),
                                                                       sep = " ")),
             breaks=c(101,1000,10000,100000,500000), labels=c("100","1000", "10 000", "100 000", "500 000"),
             guide = "none")+
  scale_color_gradientn(
    colors=palet,
    breaks=seq(0,20,by=2), limits=c(0,20),
    oob = scales::squish)   +
  labs(x="Longitude", y = "Latitude")+
  # guides(colour = guide_legend(override.aes = list(size = 10)))+
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

ggsave(cazzo,filename=paste0("D:/tilloal/Documents/LFRuns_utils/TrendAnalysis/plots/maps_cv_npeaks",haz,".jpg"), width=23, height=20, units=c("cm"),dpi=1000) 

#palet=c(hcl.colors(6, palette = "Spectral", alpha = NULL, rev = T, fixup = TRUE))
Thpoints$means[which(Thpoints$means<30)]=NA
max(Thpoints$means,na.rm=T)/70


pnpy<-ggplot(Thpoints, aes(x=means/70)) + 
  geom_histogram(color="steelblue", fill="slategray1",bins=50,alpha=0.9,lwd=1)+
  scale_y_continuous(breaks=seq(0,100000, by=10000),name="Number of pixels")+
  scale_x_continuous(breaks=seq(0,5, by=0.5),name="mean peaks per year")+
  # scale_x_log10(name=expression(paste("Upstream area ", (km^2),sep = " ")),
  #               breaks=c(100,1000,10000,100000), minor_breaks = log10_minor_break(),
  #               labels=c("100","1 000","10 000","100 000")) +
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
        legend.key.size = unit(.8, "cm"))
  #annotate("label", x=300000, y=50000, label= paste0("n_small = ",n1,"%\n n_medium = ",n2,"% \n n_large = ",n3,"%"),size=6)

pnpy
ggsave(paste0("D:/tilloal/Documents/01_Projects/FlooDrough/Figures/Revisions/histo_ppy_",haz,".jpg"), pnpy, width=20, height=15, units=c("cm"),dpi=1500)

points <- st_as_sf(Thpoints, coords = c("Var1.x", "Var2.x"), crs = 4326)
points <- st_transform(points, crs = 3035)

ppy=ggplot(basemap) +
  geom_sf(fill="gray95",color="gray10",size=0.5)+
  geom_sf(data=points,aes(col=means/70,geometry=geometry,size=upa),alpha=.9,stroke=0,shape=15)+
  coord_sf(xlim = c(min(nco[,1]),max(nco[,1])), ylim = c(min(nco[,2]),max(nco[,2])))+
  # scale_colour_manual(values = colIR, name="IR", labels=c("0"="perennial","1"="casi-perennial","2"="IRES")) +
  scale_size(range = c(0.08, 0.4), trans="sqrt",name= expression(paste("Upstream area ", (km^2),
                                                                       sep = " ")),
             breaks=c(101,1000,10000,100000,500000), labels=c("100","1000", "10 000", "100 000", "500 000"),
             guide = "none")+
  scale_color_gradientn(
    colors=palet,
    breaks=seq(0,3,by=.5), limits=c(0,3),
    oob = scales::squish,name="mean peaks per year")   +
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

ggsave(ppy,filename=paste0("D:/tilloal/Documents/LFRuns_utils/TrendAnalysis/plots/map_mean_ppy",haz,".jpg"), width=23, height=20, units=c("cm"),dpi=1000) 



length(which(final_threshold$rel_range<=50))/length(final_threshold$rel_range)
length(which(final_threshold$cv<=50))/length(final_threshold$rel_range)
range(final_threshold$means,na.rm=T)
Thpoints=inner_join(final_threshold,UpArea,by=c("catchment"="outl2"))
head(Thpoints)
#colIR=c("0"="royalblue","1"="lightblue","2"="orangered","3"="tomato","4"="purple")
points <- st_as_sf(Thpoints, coords = c("Var1.x", "Var2.x"), crs = 4326)
points <- st_transform(points, crs = 3035)


psds<-ggplot(Thpoints, aes(x=rel_range)) + 
  geom_histogram(color="steelblue", fill="slategray1",bins=50,alpha=0.9,lwd=1)+
  scale_y_continuous(breaks=seq(0,200000, by=10000),name="Number of pixels")+
  scale_x_continuous(breaks=seq(0,50, by=10),limits=c(-1,50),name="coefficient of variation of threshold")+
  # scale_x_sqrt(name=expression(paste("Upstream area ", (km^2),sep = " ")),
  #               breaks=c(100,1000,10000,100000), minor_breaks = log10_minor_break(),
  #               labels=c("100","1 000","10 000","100 000")) +
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
        legend.key.size = unit(.8, "cm"))
#annotate("label", x=300000, y=50000, label= paste0("n_small = ",n1,"%\n n_medium = ",n2,"% \n n_large = ",n3,"%"),size=6)

psds
ggsave(paste0("D:/tilloal/Documents/01_Projects/FlooDrough/Figures/Revisions/histo_cvth_",haz,".jpg"), psds, width=20, height=15, units=c("cm"),dpi=1500)


cazzo=ggplot(basemap) +
  geom_sf(fill="gray95",color="gray10",size=0.5)+
  geom_sf(data=points,aes(col=cv,geometry=geometry,size=upa),alpha=.9,stroke=0,shape=15)+
  coord_sf(xlim = c(min(nco[,1]),max(nco[,1])), ylim = c(min(nco[,2]),max(nco[,2])))+
  # scale_colour_manual(values = colIR, name="IR", labels=c("0"="perennial","1"="casi-perennial","2"="IRES")) +
  scale_size(range = c(0.08, 0.4), trans="sqrt",name= expression(paste("Upstream area ", (km^2),
                                                                       sep = " ")),
             breaks=c(101,1000,10000,100000,500000), labels=c("100","1000", "10 000", "100 000", "500 000"),
             guide = "none")+
  scale_color_gradientn(
    colors=palet,
    breaks=seq(0,20,by=2), limits=c(0,20),
    oob = scales::squish)   +
  labs(x="Longitude", y = "Latitude")+
  # guides(colour = guide_legend(override.aes = list(size = 10)))+
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

ggsave(cazzo,filename=paste0("D:/tilloal/Documents/LFRuns_utils/TrendAnalysis/plots/map_cv_threshold",haz,".jpg"), width=23, height=20, units=c("cm"),dpi=1000) 


#NOW COMPUTE DIFFERENT rlS

calcGPDReturnLevel_Single <- function(epsilon, sigma, threshold, nPeaks, sampleTimeHorizon, returnPeriod) {
  
  # Calculate expected number of events in T years
  XX <- (nPeaks / sampleTimeHorizon) * returnPeriod
  
  # Initialize result vector with NAs
  returnLevel <- rep(NA_real_, length(sigma))
  
  # 1. Handle GPD Case: epsilon is NOT NA and NOT zero
  idx_gpd <- !is.na(epsilon) & abs(epsilon) > 1e-9
  if (any(idx_gpd)) {
    returnLevel[idx_gpd] <- threshold[idx_gpd] + 
      (sigma[idx_gpd] / epsilon[idx_gpd]) * (XX[idx_gpd]^epsilon[idx_gpd] - 1)
  }
  
  # 2. Handle Gumbel Case: epsilon is NOT NA and IS zero
  idx_gum <- !is.na(epsilon) & abs(epsilon) <= 1e-9
  if (any(idx_gum)) {
    returnLevel[idx_gum] <- threshold[idx_gum] + 
      sigma[idx_gum] * log(XX[idx_gum])
  }
  
  # Note: Rows where epsilon/sigma/threshold were originally NA will remain NA
  return(returnLevel)
}

  ###load historical run -----
  haz="Drought"
for (haz in c("Flood","Drought")){
  
  if (haz == "Drought") namefile="Drought.nonfrost.Histo"
  if (haz == "Flood") namefile="flood.year.Histo"
  print(namefile)
  load(file=paste0(hydroDir,"/",haz,"/params.",namefile,".Rdata"))
  
  Paramsfl=(Paramsfl[,-c(4:9,17)])
  ParamsflH=Paramsfl
  rm(Paramsfl)
  gc()
  
  # 1. Pre-calculate the groups once (massive speed boost over 'which')
  # This creates a list where each element is the data for one year
  params_split <- split(ParamsflH, ParamsflH$Year)
  Years <- names(params_split)
  
  # 2. Pre-allocate the result matrix (prevents memory fragmentation)
  # Rows = number of spatial points/observations, Cols = number of years
  num_rows_per_year <- nrow(params_split[[1]]) 
  RL100mat <- matrix(NA_real_, nrow = num_rows_per_year, ncol = length(Years))
  colnames(RL100mat) <- Years
  
  # 3. The Optimized Loop
  for (i in seq_along(Years)) {
    year_label <- Years[i]
    cat("Processing Year:", year_label, "\n") # faster than print()
    
    # Extract the data for this year from our pre-split list
    df_year <- params_split[[year_label]]
    
    # Run the vectorized function (using the NA-safe version we discussed)
    # This calculates all 300k (or whatever your spatial count is) rows at once
    RL100s <- calcGPDReturnLevel_Single(
      epsilon = df_year$epsilonGPD, 
      sigma = df_year$sigmaGPD, 
      threshold = df_year$thresholdGPD,
      nPeaks = df_year$nPeaks, 
      sampleTimeHorizon = 70, 
      returnPeriod = 100
    )
    
    # Assign directly into the pre-allocated column (very fast)
    RL100mat[, i] <- RL100s
  }
  
  save(RL100mat,file=paste0(hydroDir,"/",haz,"/RL100x.",namefile,".Rdata"))
  ###load Socio-CF run -----
  if (haz == "Drought") namefile="Drought.nonfrost.SocCF"
  if (haz == "Flood") namefile="Flood.year.socCF"
  print(namefile)
  load(file=paste0(hydroDir,"/",haz,"/params.",namefile,".Rdata"))
  
  Paramsfl=Paramsfl[,-c(4:9,17)]
  ParamsflSCF=data.table(Paramsfl)
  rm(Paramsfl)
  gc()
  
  
  # 1. Pre-calculate the groups once (massive speed boost over 'which')
  # This creates a list where each element is the data for one year
  params_split <- split(ParamsflSCF, ParamsflH$Year)
  Years <- names(params_split)
  
  # 2. Pre-allocate the result matrix (prevents memory fragmentation)
  # Rows = number of spatial points/observations, Cols = number of years
  num_rows_per_year <- nrow(params_split[[1]]) 
  RL100mat <- matrix(NA_real_, nrow = num_rows_per_year, ncol = length(Years))
  colnames(RL100mat) <- Years
  
  # 3. The Optimized Loop
  for (i in seq_along(Years)) {
    year_label <- Years[i]
    cat("Processing Year:", year_label, "\n") # faster than print()
    
    # Extract the data for this year from our pre-split list
    df_year <- params_split[[year_label]]
    
    # Run the vectorized function (using the NA-safe version we discussed)
    # This calculates all 300k (or whatever your spatial count is) rows at once
    RL100s <- calcGPDReturnLevel_Single(
      epsilon = df_year$epsilonGPD, 
      sigma = df_year$sigmaGPD, 
      threshold = df_year$thresholdGPD,
      nPeaks = df_year$nPeaks, 
      sampleTimeHorizon = 70, 
      returnPeriod = 100
    )
    
    # Assign directly into the pre-allocated column (very fast)
    RL100mat[, i] <- RL100s
  }
  
  save(RL100mat,file=paste0(hydroDir,"/",haz,"/RL100x.",namefile,".Rdata"))
  
  ###load results from Res+WU CF run -----
  if (haz == "Drought") namefile="Drought.nonfrost.RWCF"
  if (haz == "Flood") namefile="flood.year.RWCF"
  print(namefile)
  load(file=paste0(hydroDir,"/",haz,"/params.",namefile,".Rdata"))
  Paramsfl=Paramsfl[,-c(4:9,17)]
  ParamsflRWCF=data.table(Paramsfl)
  rm(Paramsfl)
  gc()
  
  # 1. Pre-calculate the groups once (massive speed boost over 'which')
  # This creates a list where each element is the data for one year
  params_split <- split(ParamsflRWCF, ParamsflH$Year)
  Years <- names(params_split)
  
  # 2. Pre-allocate the result matrix (prevents memory fragmentation)
  # Rows = number of spatial points/observations, Cols = number of years
  num_rows_per_year <- nrow(params_split[[1]]) 
  RL100mat <- matrix(NA_real_, nrow = num_rows_per_year, ncol = length(Years))
  colnames(RL100mat) <- Years
  
  # 3. The Optimized Loop
  for (i in seq_along(Years)) {
    year_label <- Years[i]
    cat("Processing Year:", year_label, "\n") # faster than print()
    
    # Extract the data for this year from our pre-split list
    df_year <- params_split[[year_label]]
    
    # Run the vectorized function (using the NA-safe version we discussed)
    # This calculates all 300k (or whatever your spatial count is) rows at once
    RL100s <- calcGPDReturnLevel_Single(
      epsilon = df_year$epsilonGPD, 
      sigma = df_year$sigmaGPD, 
      threshold = df_year$thresholdGPD,
      nPeaks = df_year$nPeaks, 
      sampleTimeHorizon = 70, 
      returnPeriod = 100
    )
    
    # Assign directly into the pre-allocated column (very fast)
    RL100mat[, i] <- RL100s
  }
  
  save(RL100mat,file=paste0(hydroDir,"/",haz,"/RL100x.",namefile,".Rdata"))
  
  ###load results from Water CF run -----
  if (haz == "Drought") namefile="Drought.nonfrost.WCF"
  if (haz == "Flood") namefile="flood.year.WCF"
  print(namefile)
  load(file=paste0(hydroDir,"/",haz,"/params.",namefile,".Rdata"))
  
  
  Paramsfl=Paramsfl[,-c(4:9,17)]
  ParamsflWCF=data.table(Paramsfl)
  rm(Paramsfl)
  gc()
  
  # 1. Pre-calculate the groups once (massive speed boost over 'which')
  # This creates a list where each element is the data for one year
  params_split <- split(ParamsflWCF, ParamsflH$Year)
  Years <- names(params_split)
  
  # 2. Pre-allocate the result matrix (prevents memory fragmentation)
  # Rows = number of spatial points/observations, Cols = number of years
  num_rows_per_year <- nrow(params_split[[1]]) 
  RL100mat <- matrix(NA_real_, nrow = num_rows_per_year, ncol = length(Years))
  colnames(RL100mat) <- Years
  
  # 3. The Optimized Loop
  for (i in seq_along(Years)) {
    year_label <- Years[i]
    cat("Processing Year:", year_label, "\n") # faster than print()
    
    # Extract the data for this year from our pre-split list
    df_year <- params_split[[year_label]]
    
    # Run the vectorized function (using the NA-safe version we discussed)
    # This calculates all 300k (or whatever your spatial count is) rows at once
    RL100s <- calcGPDReturnLevel_Single(
      epsilon = df_year$epsilonGPD, 
      sigma = df_year$sigmaGPD, 
      threshold = df_year$thresholdGPD,
      nPeaks = df_year$nPeaks, 
      sampleTimeHorizon = 70, 
      returnPeriod = 100
    )
    
    # Assign directly into the pre-allocated column (very fast)
    RL100mat[, i] <- RL100s
  }
  
  save(RL100mat,file=paste0(hydroDir,"/",haz,"/RL100x.",namefile,".Rdata"))
}
