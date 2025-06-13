#Code to identify trend threshold for each catchment, 
#if threshold is different, I need a flag and I check manually what is the best

getwd()
setwd("D:/tilloal/Documents/LFRuns_utils/TrendAnalysis/")
source("R/functions_trends.R")
library(RtsEva)
library(lubridate)
library(zoo)
suppressWarnings(suppressMessages(library(ncdf4)))
suppressWarnings(suppressMessages(library(sf)))
suppressWarnings(suppressMessages(library(rnaturalearth)))
suppressWarnings(suppressMessages(library(rnaturalearthdata)))
suppressWarnings(suppressMessages(library(rgeos)))
suppressWarnings(suppressMessages(library(dplyr)))
main_path = "D:/tilloal/Documents/06_Floodrivers/"
valid_path = paste0(main_path)

dis_path<-paste0(main_path,'HERA/')
load(file=paste0(dis_path,"HERA6h_H07_19502020.Rdata"))
Station_data_IDs <- as.vector(t(Q_sim[1, ]))[-1]
tqx=as.POSIXct((Q_sim[,1]*3600*24), origin="1979-01-01 00:00:00")
txx=tqx[-1]
Q_simH=Q_sim

rm(Q_sim)


#extract values for 1 location

#remove locations without discharge
wolala=which(is.na(Q_simH[2,]))-1
Q_simH=Q_simH[,-which(is.na(Q_simH[2,]))]
Station_data_IDs=Station_data_IDs[-wolala]

#match HYBAS_id with station_IDs
# hybmatch=match(Station_data_IDs,cst7$outlets)
# hybas_IDs=cst7$HYBAS_ID[hybmatch]
txx=tqx[-1]
rmv=which(year(txx)==1950)
txx=txx[-rmv]
# Q_H=Q_H[-rmv]
# Q_SCF=Q_SCF[-rmv]

check_timeserie2=function(timeseries,yro){
  ts_years <- as.integer((lubridate::year(timeseries)))
  year_check <- yro %in% ts_years
  runs <- rle(year_check)
  rf=which(runs$values==FALSE)
  if (any(runs$lengths[rf] >= 2)) {
    return(FALSE)
  }
  else {
    return(TRUE)
  }
}

tsEvaFindTrendThreshold2<-function(series, timeStamps, timeWindow){
  ptn = timeStamps[which(!is.na(series))]
  bounds = unique(lubridate::year(ptn))
  nr <- rep(1, length(series))
  nr = nr + rnorm(length(series), 0, 1e-05)
  sts <- c()
  lnegs = c()
  pctd = c()
  pcts <- seq(0.4, 0.95, by = 0.05)
  for (iter in 1:length(pcts)) {
    thrsdt <- quantile(series, pcts[iter], na.rm = TRUE)
    series_no_na <- series
    series_no_na[which(is.na(series_no_na))] <- -9999
    serieb <- series_no_na
    timeb = timeStamps
    timeb = timeb[-which(serieb < thrsdt)]
    serieb[which(serieb < thrsdt)] <- NA
    checkY=check_timeserie2(timeb,bounds)
    if (checkY == FALSE) {
      print(paste0("not all years - q= ",pcts[iter]))
      break
    }
    rs <- tsEvaDetrendTimeSeries(timeStamps, serieb, timeWindow, 
                                 fast = T)

    norm_trend <- rs@trendSeries/mean(rs@trendSeries, na.rm = TRUE)
    dtr1 = serieb - rs@trendSeries
    lneg = length(which(dtr1 < 0))
    stab <- cor(nr, norm_trend, use = "pairwise.complete.obs")
    if (iter == 1) 
       nr <- norm_trend
    if (lneg >= 1) 
      lnegs = c(lnegs, lneg)
      sts <- c(sts, stab)
      pctd = c(pctd, pcts[iter])
    }
  dow=abs(diff(sts))[-1]
  
  rval = pctd[length(pctd)]
  
  if(max(dow)>0.2){
    print("breaking point")
    rval = pctd[which.max(dow)+1]
  }
  if (sum(lnegs) > 1) {
    rval = pctd[which.min(lnegs)]
  }
  
  return(rval)
}


#modify the function by using findpeaks instead
tsEvaFindTrendThreshold3<-function(series, timeStamps){
  #compute monthly maxs
  tserie=data.frame(series,timeStamps)
  
  monthly_max=tserie %>% 
    mutate(month = floor_date(timeStamps, "month")) %>% 
    group_by(month) %>% 
    summarise(max_value = max(series))
  
  #use findpeaks
  med=median(series,na.rm=T)
  peakox=findpeaks(tserie$series,threshold = med)
  peakos=data.frame(time=timeStamps[peakox[,2]],values=peakox[,1])

  
  return(list(peaks=peakos,monmax=monthly_max))
}

trendPeaks=tsEvaFindTrendThreshold3(series,timeStamps)
timeWindow = 365.25*30; #time windows in days, the correction is done within the functions

library(foreach)
library(doParallel)
library(RtsEva)
# Register parallel backend
cl <- makeCluster(4)
registerDoParallel(cl)
Trth_H_list=c()
tail="low"
Trth_H_list <- foreach(id = 1:length(Station_data_IDs), .combine = rbind, .packages=c("RtsEva","zoo")) %dopar% {
  #print(id)
  stid <- Station_data_IDs[id]
  catch <- stid
  Q_H <- Q_simH[,id+1][-1]
  Q_H <- Q_H[-rmv]
  
  data <- data.frame(txx,Q_H)
  names(data) <- c("date","Qs")
  if (tail=="high"){
    # Extract the maximum daily value
    series <- max_daily_value(data)
    timeAndSeriesH=series
    names(timeAndSeriesH)=c("timestamp","data")
  } else if (tail=="low"){
    
    minPeakDistanceInDays=30
    #7 day average flow for drought analysis
    WindowSize=7
    names(data)=c("date","Qs")
    dt1=min(diff(data$date),na.rm=T)
    dt=as.numeric(dt1)
    tdim=attributes(dt1)$units
    if (tdim=="hours") dt=dt/24
    nRunMn = ceiling(WindowSize/dt);
    colnames(data)
    data$Q7=rollmean(data$Qs,nRunMn, align = "right", fill=NA)
    timeStamps=data$date
    series=data$Q7
    start_index=1
    indices_to_extract <- seq(from = start_index, to = length(series), by = WindowSize/dt)
    series=series[indices_to_extract]
    series=-1*series
    timeStamps=timeStamps[indices_to_extract]
    timeAndSeriesH=data.frame(timeStamps,series)
    names(timeAndSeriesH)=c("timestamp","data")
  }
  
  TrendTh_H1 <- tsEvaFindTrendThreshold(series=timeAndSeriesH$data, timeStamps=timeAndSeriesH$timestamp, timeWindow)
  
  if (is.null(TrendTh_H1))TrendTh_H1=NA
  
  TrendTh_H2 <- tsEvaFindTrendThreshold2(series=timeAndSeriesH$data, timeStamps=timeAndSeriesH$timestamp, timeWindow)
  m1=median(timeAndSeriesH$data)
  o=findpeaks(timeAndSeriesH$data, threshold = m1)
  plot(o[,1])
  mean(o[,1])
  m1=mean(timeAndSeriesH$data)
  if (is.null(TrendTh_H2))TrendTh_H2=NA
  Trth_H <- c(catch,TrendTh_H1,TrendTh_H2)
  
  Trth_H
}

stopCluster(cl)
write.csv(Trth_H_list,file="trenTH_Hlow.csv")

rm(Q_simH)
dis_path<-paste0(main_path,'HERA_RWstat/')
load(file=paste0(dis_path,"HERA_CFR6h_H07_19502020.Rdata"))
Q_simCFR=Q_sim
Q_simCFR=Q_simCFR[,-which(is.na(Q_simCFR[2,]))]
rm(Q_sim)
cl <- makeCluster(4)
registerDoParallel(cl)
Trth_R_list <- foreach(id = 1:length(Station_data_IDs), .combine = rbind,.packages=c("RtsEva","zoo")) %dopar% {
  
  #print(id)
  stid <- Station_data_IDs[id]
  catch <- stid

  Q_RCF <- Q_simCFR[,id+1][-1]

  #Q_RCF <- Q_RCF[-rmv]

  data <- data.frame(txx,Q_RCF)
  names(data) <- c("date","Qs")
  
  if (tail=="high"){
    # Extract the maximum daily value
    series <- max_daily_value(data)
    timeAndSeriesR=series
    names(timeAndSeriesR)=c("timestamp","data")
  } else if (tail=="low"){
    
    minPeakDistanceInDays=30
    #7 day average flow for drought analysis
    WindowSize=7
    names(data)=c("date","Qs")
    dt1=min(diff(data$date),na.rm=T)
    dt=as.numeric(dt1)
    tdim=attributes(dt1)$units
    if (tdim=="hours") dt=dt/24
    nRunMn = ceiling(WindowSize/dt);
    colnames(data)
    data$Q7=rollmean(data$Qs,nRunMn, align = "right", fill=NA)
    timeStamps=data$date
    series=data$Q7
    start_index=1
    indices_to_extract <- seq(from = start_index, to = length(series), by = WindowSize/dt)
    series=series[indices_to_extract]
    series=-1*series
    timeStamps=timeStamps[indices_to_extract]
    timeAndSeriesR=data.frame(timeStamps,series)
    names(timeAndSeriesR)=c("timestamp","data")
  }
  # 
  # data <- data.frame(txx,Q_SCF)
  # names(data) <- c("date","Qs")
  # series <- max_daily_value(data)
  # timeAndSeriesS <- series
  # names(timeAndSeriesS) <- c("timestamp","data")
  

  TrendTh_R1 <- tsEvaFindTrendThreshold(series=timeAndSeriesR$data, timeStamps=timeAndSeriesR$timestamp, timeWindow)
  if (is.null(TrendTh_R1))TrendTh_R1=NA
  TrendTh_R2 <- tsEvaFindTrendThreshold2(series=timeAndSeriesR$data, timeStamps=timeAndSeriesR$timestamp, timeWindow)
  if (is.null(TrendTh_R2))TrendTh_R2=NA
  Trth_R <- c(catch,TrendTh_R1,TrendTh_R2)

  
  Trth_R
  
}

stopCluster(cl)

write.csv(Trth_R_list,file="trenTH_Rlow.csv")
rm(Q_simCFR)
dis_path<-paste0(main_path,'HERA_CFS/')
load(file=paste0(dis_path,"HERA_CFS6h_H07_19502020.Rdata"))
Q_simCFS=Q_sim
Q_simCFS=Q_simCFS[,-which(is.na(Q_simCFS[2,]))]

cl <- makeCluster(4)
registerDoParallel(cl)
Trth_S_list <- foreach(id = 1:length(Station_data_IDs), .combine = rbind, .packages=c("RtsEva","zoo")) %dopar% {
  
  
  #print(id)
  stid <- Station_data_IDs[id]
  catch <- stid
  
  Q_SCF <- Q_simCFS[,id+1][-1]
  Q_SCF <- Q_SCF[-rmv]

  data <- data.frame(txx,Q_SCF)
  names(data) <- c("date","Qs")
  if (tail=="high"){
    # Extract the maximum daily value
    series <- max_daily_value(data)
    timeAndSeriesS=series
    names(timeAndSeriesS)=c("timestamp","data")
  } else if (tail=="low"){
    
    minPeakDistanceInDays=30
    #7 day average flow for drought analysis
    WindowSize=7
    names(data)=c("date","Qs")
    dt1=min(diff(data$date),na.rm=T)
    dt=as.numeric(dt1)
    tdim=attributes(dt1)$units
    if (tdim=="hours") dt=dt/24
    nRunMn = ceiling(WindowSize/dt);
    colnames(data)
    data$Q7=rollmean(data$Qs,nRunMn, align = "right", fill=NA)
    timeStamps=data$date
    series=data$Q7
    start_index=1
    indices_to_extract <- seq(from = start_index, to = length(series), by = WindowSize/dt)
    series=series[indices_to_extract]
    series=-1*series
    timeStamps=timeStamps[indices_to_extract]
    timeAndSeriesS=data.frame(timeStamps,series)
    names(timeAndSeriesS)=c("timestamp","data")
  }
  
  TrendTh_S1 <- tsEvaFindTrendThreshold(series=timeAndSeriesS$data, timeStamps=timeAndSeriesS$timestamp, timeWindow)
  if (is.null(TrendTh_S1))TrendTh_S1=NA
  TrendTh_S2 <- tsEvaFindTrendThreshold2(series=timeAndSeriesS$data, timeStamps=timeAndSeriesS$timestamp, timeWindow)
  if (is.null(TrendTh_S2))TrendTh_S2=NA
  Trth_S <- c(catch,TrendTh_S1,TrendTh_S2)
  
  
  Trth_S
  
}
stopCluster(cl)
write.csv(Trth_S_list,file="trenTH_Slow.csv")





plot(Trth_H_list[,3], Trth_S_list[,3],ylim=c(0.5,1),xlim=c(0.5,1))


dx=data.frame(hist=Trth_H_list[,3],scf=Trth_S_list[,3], catch=Trth_H_list[,1])
dx2=dx[which(dx$hist!=dx$scf),]
dx$common=dx$hist
dx$common[which(dx$hist>1)]=0.5

save(dx,file="thresholds_low_catchments.rdata")
#the rule could be to use the soccf th if no value for reservoir and the reservoir run instead

#large reduction of unmatching threshold with new method




#match dx2 with catchments and plot



## Spatial data for catchments ----
#Set data directory
hydroDir<-("D:/tilloal/Documents/LFRuns_utils/data")


#load my shitty rdata
load(file=paste0(hydroDir,"/dataHybas07_ids.Rdata"))
outletname="outletsv8_hybas07_01min"
nameout="UCH07"

### Hybas07 ----
Catchmentrivers7=read.csv(paste0(hydroDir,"/Catchments/from_hybas_eu_onlyid.csv"),encoding = "UTF-8", header = T, stringsAsFactors = F)
hybas07 <- read_sf(dsn = paste0(hydroDir,"/Catchments/hydrosheds/hybas_eu_lev07_v1c.shp"))
hybasf7=fortify(hybas07) 
Catamere07=inner_join(hybasf7,Catchmentrivers7,by= "HYBAS_ID")
outf=outletopen(hydroDir,outletname)
outf$latlong=paste(round(outf$Var1,4),round(outf$Var2,4),sep=" ")
Catamere07$llcoord=paste(round(Catamere07$POINT_X,4),round(Catamere07$POINT_Y,4),sep=" ") 
Catf7=inner_join(Catamere07,outf,by= c("llcoord"="latlong"))
cst7=st_transform(Catf7,  crs=3035)


reservoirdog=outletopen(paste0(hydroDir,"/reservoirs"),"res_ratio_diff_2020-1951")
reservoirdog$llcoord=paste(round(reservoirdog$Var1,4),round(reservoirdog$Var2,4),sep=" ") 


dx=dx[-which(dx$scf>1),]
dx=dx[-which(dx$hist>1),]

#match dx2 with catchments


#other method to extract right catchments
#efas rnt
out1=outletopen(hydroDir,"efas_rnet_100km_01min")
out1$latlong=paste(round(out1$Var1,4),round(out1$Var2,4),sep=" ")
out2=inner_join(out1,outf, by="latlong")


### Plot parameters ----
palet2=c(hcl.colors(9, palette = "Blues", alpha = NULL, rev = TRUE, fixup = TRUE))
outll=outletopen(hydroDir,outletname)
cord.dec=outll[,c(2,3)]
cord.dec = SpatialPoints(cord.dec, proj4string=CRS("+proj=longlat"))
cord.UTM <- spTransform(cord.dec, CRS("+init=epsg:3035"))
nco=cord.UTM@coords
coord.lamb=spTransform(cord.dec, CRS("+proj=laea +lat_0=55 +lon_0=20 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=km +no_defs"))
nci=coord.lamb@coords
world <- ne_countries(scale = "medium", returnclass = "sf")
Europe <- world[which(world$continent == "Europe"),]
e2=st_transform(Europe,  crs=3035)
w2=st_transform(world,  crs=3035)
tsize=12
osize=12
catchmentproj=cst7
Impdates=seq(1951,2020,by=1)
valuenames=paste0("X",Impdates)
catmap=catchmentproj
basemap=w2




#load reservoir influence
reservoirdog=outletopen(paste0(hydroDir,"/reservoirs"),"res_ratio_diff_2020-1951")
reservoirdog$llcoord=paste(round(reservoirdog$Var1,4),round(reservoirdog$Var2,4),sep=" ") 
reservoir_bg=inner_join(reservoirdog,outf,by=c("llcoord"="latlong"))


hybmatch=match(dx2$catch,cst7$outlets)
dx2$hybasID=cst7$HYBAS_ID[hybmatch]


data=right_join(catmap,dx2,by = c("outlets"="catch"))

palet=c(hcl.colors(11, palette = "RdYlBu", alpha = NULL, rev = F, fixup = TRUE))
br=c(-100,-80,-50,-30,-20,-10,-8,-6,-4,-2,0,2,4,6,8,10,20,30,50,80,100)
br=c(-100,-80,-50,-20,-10,0,10,20,50,80,100)
limi=c(-50,50)
trans=scales::modulus_trans(.8)
colNA="darkgrey"

ggplot(basemap) +
  geom_sf(fill="white")+
  geom_sf(data=data,aes(geometry=geometry),color="transparent", fill="red")+ 
  geom_sf(fill=NA, color="grey") +
  coord_sf(xlim = c(min(nco[,1]),max(nco[,1])), ylim = c(min(nco[,2]),max(nco[,2])))+
  labs(x="Longitude", y = "Latitude")+
  guides(fill = guide_colourbar(barwidth = 20, barheight = .8))+
  theme(axis.title=element_text(size=tsize),
        panel.background = element_rect(fill = "aliceblue", colour = "grey1"),
        panel.border = element_rect(linetype = "solid", fill = NA, colour="black"),
        legend.title = element_text(size=tsize),
        legend.text = element_text(size=osize),
        legend.position = "bottom",
        legend.title.position = "top",
        panel.grid.major = element_line(colour = "grey70"),
        panel.grid.minor = element_line(colour = "grey90"),
        legend.key = element_rect(fill = "transparent", colour = "transparent"),
        legend.key.size = unit(.8, "cm"))

#Check if it is correlated with reservoirs

#match reservoirdog with catchment and link the differences with reservoir status

reservoir_bg=inner_join(reservoirdog,outf,by=c("llcoord"="latlong"))

#match this with dx

mrx=match(dx$catch,reservoir_bg$outlets.y)

dx$reservoir=reservoir_bg$outlets.x[mrx]


dx2=dx[which(dx$hist!=dx$scf),]
cres=dx2$catch[which(dx2$hist>1)]

dx3=dx2[which(dx2$hist>1),]
dx4=dx3[2,]
catch =cres[1]

id=which(Station_data_IDs==catch)

Q_H=Q_simH[,id+1][-1]
Q_SCF=Q_simCFS[,id+1][-1]
Q_RCF=Q_simCFR[,id+1][-1]

Q_H=Q_H[-rmv]
Q_SCF=Q_SCF[-rmv]
Q_RCF=Q_RCF[-rmv]

data=data.frame(txx,Q_H)
names(data)=c("date","Qs")
# Extract the maximum daily value
series <- max_daily_value(data)
timeAndSeriesH=series
names(timeAndSeriesH)=c("timestamp","data")

data=data.frame(txx,Q_RCF)
names(data)=c("date","Qs")
# Extract the maximum daily value
series <- max_daily_value(data)
timeAndSeriesR=series
names(timeAndSeriesR)=c("timestamp","data")

data=data.frame(txx,Q_SCF)
names(data)=c("date","Qs")
# Extract the maximum daily value
series <- max_daily_value(data)
timeAndSeriesS=series
names(timeAndSeriesS)=c("timestamp","data")
print("coucou")

plot(timeAndSeriesH)
points(timeAndSeriesS,col="blue")

TrendTh_H1=tsEvaFindTrendThreshold2(series=timeAndSeriesH$data, timeStamps=timeAndSeriesH$timestamp, timeWindow)
TrendTh_R1=tsEvaFindTrendThreshold2(series=timeAndSeriesR$data, timeStamps=timeAndSeriesR$timestamp, timeWindow)
TrendTh_S1=tsEvaFindTrendThreshold2(series=timeAndSeriesS$data, timeStamps=timeAndSeriesS$timestamp, timeWindow)
#check all locations where no threshold was found for reservoir runs

#try to compute the trend

timeStamps=timeAndSeriesH$timestamp
series=timeAndSeriesH$data
timeWindow
TrendTh=0.5
trasfData = tsEvaTransformSeriesToStationaryPeakTrend( timeStamps, series, timeWindow, TrendTh);

message(paste0("trend threshold= ", TrendTh))
qd <- quantile(series, TrendTh, na.rm = T)
serieb <- series
serieb[which(serieb < qd)] <- NA
plot(serieb)
rs <- tsEvaDetrendTimeSeries(timeStamps, serieb, timeWindow)


plot(trasfData$trendSeries)
dxx=dx[which(dx$hist>1),]







data=right_join(catmap,dx4,by = c("outlets"="catch"))

palet=c(hcl.colors(11, palette = "RdYlBu", alpha = NULL, rev = F, fixup = TRUE))
br=c(-100,-80,-50,-30,-20,-10,-8,-6,-4,-2,0,2,4,6,8,10,20,30,50,80,100)
br=c(-100,-80,-50,-20,-10,0,10,20,50,80,100)
limi=c(-50,50)
trans=scales::modulus_trans(.8)
colNA="darkgrey"

ggplot(basemap) +
  geom_sf(fill="white")+
  geom_sf(data=data,aes(geometry=geometry),color="transparent", fill="red")+ 
  geom_sf(fill=NA, color="grey") +
  coord_sf(xlim = c(min(nco[,1]),max(nco[,1])), ylim = c(min(nco[,2]),max(nco[,2])))+
  labs(x="Longitude", y = "Latitude")+
  guides(fill = guide_colourbar(barwidth = 20, barheight = .8))+
  theme(axis.title=element_text(size=tsize),
        panel.background = element_rect(fill = "aliceblue", colour = "grey1"),
        panel.border = element_rect(linetype = "solid", fill = NA, colour="black"),
        legend.title = element_text(size=tsize),
        legend.text = element_text(size=osize),
        legend.position = "bottom",
        legend.title.position = "top",
        panel.grid.major = element_line(colour = "grey70"),
        panel.grid.minor = element_line(colour = "grey90"),
        legend.key = element_rect(fill = "transparent", colour = "transparent"),
        legend.key.size = unit(.8, "cm"))

