setwd(dirname(rstudioapi::getSourceEditorContext()$path))
getwd()
source("functions_trends.R")

#2 Pre-loaded results -----------
#Set data directory
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

###load Socio-CF run -----
if (haz == "Drought") namefile="Drought.nonfrost.SocCF"
if (haz == "Flood") namefile="Flood.year.socCF"

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
if (haz == "Drought") namefile="Drought.nonfrost.RWCF"
if (haz == "Flood") namefile="flood.year.RWCF"
load(file=paste0(hydroDir,"/data/",haz,"/params.",namefile,".Rdata"))
load(file=paste0(hydroDir,"/data/",haz,"/RL100.",namefile,".Rdata"))
load(file=paste0(hydroDir,"/data/",haz,"/peaks.",namefile,".Rdata"))
RLGPDflRWCF=RLGPDfl
Paramsfl=Paramsfl[,-c(4:9,17)]
ParamsflRWCF=data.table(Paramsfl)
PeakRWCF=Peaksave
rm(Paramsfl,RLGPDfl)
gc()

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

rm(catmap)
gc()


#histogram of shape files
hist(ParamsflH$epsilonGPD,xlim=c(-1,2),breaks=10000)
max(ParamsflH$epsilonGPD,na.rm=T)
quantile(ParamsflH$epsilonGPD,0.5,na.rm=T)

#3 Data cleaning -----
RLGPDflSCF=as.data.frame(RLGPDflSCF)
RLGPDflH=as.data.frame(RLGPDflH)
RLGPDflWCF=as.data.frame(RLGPDflWCF)
RLGPDflRWCF=as.data.frame(RLGPDflRWCF)

##3.1 Identify locations where last years RL were not computed -----
for (run in c(1,2,3,4)){
  print(run)
  if (run==1) RLcheck=RLGPDflSCF
  if (run==2) RLcheck=RLGPDflRWCF
  if (run==3) RLcheck=RLGPDflWCF
  if (run==4) RLcheck=RLGPDflH
  
  pbloc=RLcheck$unikout[which(is.nan((RLcheck$Y2020)))]
  if (length(pbloc)>1){
    for (f in 1:length(pbloc)){
      fi=pbloc[f]
      print(fi)
      c=RLcheck[which(RLcheck$unikout==fi),]
      last=which(is.nan(as.numeric(c)))[1]-1
      c[which(is.nan(as.numeric(c)))]=c[last]
      RLcheck[which(RLcheck$unikout==fi),]=c
    }
    if (run==1) RLGPDflSCF=RLcheck
    if (run==2) RLGPDflRWCF=RLcheck
    if (run==3) RLGPDWCF=RLcheck
    if (run==4) RLGPDflH=RLcheck
  }
}


##3.2 load IRES status -----
if (haz=="Drought"){
  hazard="Drought"
  season="nonfrost"
  mmx=""
  load(file=paste0(hydroDir,"/",hazard,"/IRES.",season,".Histo",mmx,".Rdata"))
  # In loade data:
  ##IRES=1 correspond to casi-perrenial rivers
  ##IRES=2 correspond to river with a constant floor low flow
  ##IRES=3 correspond to an IRES
  
  #Bring class 1 and 2 together here
  if (length(which(IRES_save$IRES==2))){
  IRES_save$IRES[which(IRES_save$IRES==2)]=1
  }
  IRES_save$IRES[which(IRES_save$IRES==3)]=2
  IRES_Histo=IRES_save
  
  load(file=paste0(hydroDir,"/",hazard,"/IRES.",season,".WCF",mmx,".Rdata"))
  if (length(which(IRES_save$IRES==2))){
    IRES_save$IRES[which(IRES_save$IRES==2)]=1
  }
  IRES_save$IRES[which(IRES_save$IRES==3)]=2
  IRES_WCF=IRES_save
  
  load(file=paste0(hydroDir,"/",hazard,"/IRES.",season,".RWCF",mmx,".Rdata"))
  if (length(which(IRES_save$IRES==2))){
    IRES_save$IRES[which(IRES_save$IRES==2)]=1
  }
  IRES_save$IRES[which(IRES_save$IRES==3)]=2
  IRES_RWCF=IRES_save
  
  load(file=paste0(hydroDir,"/",hazard,"/IRES.",season,".SocCF",mmx,".Rdata"))
  if (length(which(IRES_save$IRES==2))){
    IRES_save$IRES[which(IRES_save$IRES==2)]=1
  }
  IRES_save$IRES[which(IRES_save$IRES==3)]=2
  IRES_SocCF=IRES_save
  
  #identification of missing pixels (no fit) in SocCF
  msoc=match(IRES_Histo$catlist,IRES_SocCF$catlist)
  IRES_Histo$catlist[which(is.na(msoc))]
  
  #combination of all runs to detect any river pixel which is IRES
  # in at least one run
  IRES_comb=full_join(IRES_Histo,IRES_SocCF,by="catlist")
  IRES_comb=full_join(IRES_comb,IRES_WCF,by="catlist")
  IRES_comb=full_join(IRES_comb,IRES_RWCF,by="catlist")
  IRES_comb$IRES.y[which(is.na(IRES_comb$IRES.y))]=0
  IRES_comb$IRES.x.x[which(is.na(IRES_comb$IRES.x.x))]=0
  IRES_comb$IRES.y.y[which(is.na(IRES_comb$IRES.y.y))]=0
  
  #conservative aproach: if the pixel is an IRES in at least one of the 4 runs, 
  #it is considered as IRES
  IRES_comb$gen_IR=ceiling((IRES_comb$IRES.x+IRES_comb$IRES.y+IRES_comb$IRES.x.x+IRES_comb$IRES.y.y)/4)
  IRES_comb$dtectblss=(IRES_comb$IRES.x+IRES_comb$IRES.y+IRES_comb$IRES.x.x+IRES_comb$IRES.y.y)
  
  names(IRES_comb)[c(2,3,4,5)]=c("Histo","SCF","WCF","RWCF")
  IR_locs=which(IRES_comb$gen_IR>=2)
  length(which(IRES_comb$gen_IR==1))/length(IRES_comb$catlist)
  length(which(IRES_comb$gen_IR>=2))/length(IRES_comb$catlist)
  
  ###[Plot] Figure S.6 - Intermittent rivers plot ----
  
  IRpoints=inner_join(IRES_comb,UpArea,by=c("catlist"="outl2"))
  
  colIR=c("0"="royalblue","1"="lightblue","2"="orangered","3"="tomato","4"="purple")
  points <- st_as_sf(IRpoints, coords = c("Var1.x", "Var2.x"), crs = 4326)
  points <- st_transform(points, crs = 3035)
  
  IRmap=ggplot(basemap) +
    geom_sf(fill="gray95",color="gray10",size=0.5)+
    geom_sf(data=points,aes(col=factor(gen_IR),geometry=geometry,size=upa),alpha=.9,stroke=0,shape=15)+
    coord_sf(xlim = c(min(nco[,1]),max(nco[,1])), ylim = c(min(nco[,2]),max(nco[,2])))+
    scale_colour_manual(values = colIR, name="IR", labels=c("0"="perennial","1"="casi-perennial","2"="IRES")) +
    scale_size(range = c(0.08, 0.4), trans="sqrt",name= expression(paste("Upstream area ", (km^2),
                                                                         sep = " ")),
               breaks=c(101,1000,10000,100000,500000), labels=c("100","1000", "10 000", "100 000", "500 000"),
               guide = "none")+
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
  
  #ggsave(paste0("D:/tilloal/Documents/LFRuns_utils/TrendAnalysis/plots/IRES_allrunsNew19.jpg"), IRmap, width=20, height=20, units=c("cm"),dpi=1000) 
}

##3.3 drought RL corrections ----------------
if (haz=="Drought"){
  RLIRH<- RLGPDflH[IR_locs,]
  #putting NAs on IRES pixels
  RLGPDflH[IR_locs,c(1:70)]<-NA
  #reversing RL values as the rev tranformation was used for low flow
  RLGPDflH[,c(1:70)]=-RLGPDflH[,c(1:70)]
  #drought RL that are negative are set to 0
  for (id in c(1:70)){
    RLGPDflH[which(RLGPDflH[,id]<0),id]=0
    RLGPDflH[which(is.infinite(RLGPDflH[,id])),id]=NA
  }
  
  RLIRWCF<- RLGPDflWCF[IR_locs,]
  #putting NAs on IRES pixels
  RLGPDflWCF[IR_locs,c(1:70)]<-NA
  #reversing RL values as the rev tranformation was used for low flow
  RLGPDflWCF[,c(1:70)]=-RLGPDflWCF[,c(1:70)]
  #drought RL that are negative are set to 0
  for (id in c(1:70)){
    RLGPDflWCF[which(RLGPDflWCF[,id]<0),id]=0
    RLGPDflWCF[which(is.infinite(RLGPDflWCF[,id])),id]=NA
  }
  
  RLIRSCF<- RLGPDflSCF[IR_locs,]
  #putting NAs on IRES pixels
  RLGPDflSCF[IR_locs,c(1:70)]<-NA
  #reversing RL values as the rev tranformation was used for low flow
  RLGPDflSCF[,c(1:70)]=-RLGPDflSCF[,c(1:70)]
  #drought RL that are negative are set to 0
  for (id in c(1:70)){
    RLGPDflSCF[which(RLGPDflSCF[,id]<0),id]=0
    RLGPDflSCF[which(is.infinite(RLGPDflSCF[,id])),id]=NA
  }
  
  
  RLIRRWCF<- RLGPDflRWCF[IR_locs,]
  #putting NAs on IRES pixels
  RLGPDflRWCF[IR_locs,c(1:70)]<-NA
  #reversing RL values as the rev tranformation was used for low flow
  RLGPDflRWCF[,c(1:70)]=-RLGPDflRWCF[,c(1:70)]
  #drought RL that are negative are set to 0
  for (id in c(1:70)){
    RLGPDflRWCF[which(RLGPDflRWCF[,id]<0),id]=0
    RLGPDflRWCF[which(is.infinite(RLGPDflRWCF[,id])),id]=NA
  }
}

##3.4 Irrealistic shape parameter removal

if (haz=="Flood") {
  na="10-y high flow \n(l/s/km2)"
  shp_bnd=c(1,-0.5)
}
if (haz=="Drought") {
  na="10-y low flows \n(l/s/km2)"
  shp_bnd=c(0,-1.5)
}

#shape parameter vector to remove pixels with crazy shape parameter
Shapepar1=ParamsflH[,c(1,2,4)]
Shapepar2=ParamsflSCF[,c(1,2,4)]
Shapepar3=ParamsflRWCF[,c(1,2,4)]
Shapepar4=ParamsflWCF[,c(1,2,4)]
Shapeparf=cbind(Shapepar2[which(Shapepar1$Year==2015),],Shapepar3$epsilonGPD[which(Shapepar2$Year==2015)],
                Shapepar4$epsilonGPD[which(Shapepar3$Year==2015)],Shapepar1$epsilonGPD[which(Shapepar4$Year==2015)])

rmpixs1=which(Shapepar1$epsilonGPD>shp_bnd[1] | Shapepar1$epsilonGPD<shp_bnd[2])
length(rmpixs1)/length(Shapepar1$catchment)
rmpixs2=which(Shapepar2$epsilonGPD>shp_bnd[1] | Shapepar2$epsilonGPD<shp_bnd[2])
rmpixs3=which(Shapepar3$epsilonGPD>shp_bnd[1] | Shapepar3$epsilonGPD<shp_bnd[2])
rmpixs4=which(Shapepar4$epsilonGPD>shp_bnd[1] | Shapepar4$epsilonGPD<shp_bnd[2])
rmpixs=unique(c(rmpixs1,rmpixs2,rmpixs3,rmpixs4))
length(rmpixs)/length(Shapepar1$catchment)*100
rmp2=na.omit(unique(ParamsflSCF$catchment[rmpixs]))
ShapeparRM=Shapeparf[(match(rmp2,Shapeparf$catchment)),]
Shapeparf=Shapeparf[-(match(rmp2,Shapeparf$catchment)),]

Shapeparf$mean=(Shapeparf$epsilonGPD+Shapeparf$V2+Shapeparf$V3+Shapeparf$V4)/4
Shapeparf$sd=sqrt(((Shapeparf$epsilonGPD-Shapeparf$mean)^2+(Shapeparf$V2-Shapeparf$mean)^2+
  (Shapeparf$V3-Shapeparf$mean)^2+(Shapeparf$V4-Shapeparf$mean)^2)/4)
##3.4 Large scale error computation --------------

RlevErrtH=c()
for (yr in (1951:2020)){
  print(yr)
  ParamSpecial=ParamsflSCF[which(ParamsflSCF$Year==yr),]
  ParamSpecial<-ParamSpecial[-match(rmp2,ParamSpecial$catchment),]
  RlevErri <- calculate_return_levels(ParamSpecial,ci=2)
  RlevErri$year=yr
  mapu=na.omit((match(ParamSpecial$catchment,UpArea$outl2)))
  ParamSpecial$upa=UpArea$upa[mapu]
  #RlevErri=RlevErri/ParamSpecial$upa*1000
  mr1=mean(RlevErri$returnLevels,na.rm=T)
  mr1er=mean(RlevErri$returnLevelErr,na.rm=T)
  RLEr=c(mr1,mr1er,yr)
  RlevErrtH=rbind(RLEr,RlevErrtH)
  
}
RlevErrtH=data.frame(RlevErrtH)
if(haz=="Drought") RlevErrtH$X1=-RlevErrtH$X1
#base 100 in 1955
idb=which(RlevErrtH$X3==1955)
RlevErrtH$X4=RlevErrtH$X1/RlevErrtH$X1[idb]*100
RlevErrtH$X5=RlevErrtH$X2/RlevErrtH$X1[idb]*100
xlabs=seq(1950,2020,5)
clabels=c("Drought","Flood")


###[Plot] - Supplement - Mean change in Hazard intensity with CI ----
nplot="Change in 10yRL"
br=seq(80,120,2)
ggplot() +
  # IQR represented as rectangles
  geom_line(data=RlevErrtH,aes(x=X3, y=X4),
            lwd=2) +
  geom_ribbon(data = RlevErrtH, aes(x = X3, ymin = X4-X5, ymax = X4+X5),fill="blue", 
              alpha = 0.2) +

  # Y-axis settings
  scale_y_continuous(name = nplot, breaks = br, trans=scales::modulus_trans(.6)) +
  # 
  # X-axis settings
  scale_x_continuous(breaks = xlabs, labels = xlabs, name = "Years",limits=c(1955,2020),
                     minor_breaks = seq(1955,2005,10), expand = c(.001,0.001)) +
  
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
  ggtitle("Europe")



### Error on the estimation of RL in 1955 -----
# SocCF scenario
unish=unique(ParamsflSCF$catchment)
paraU=ParamsflSCF[which(ParamsflSCF$Year==1955),]

ParamSpecial<-paraU[-match(rmp2,paraU$catchment),]

Paramf=ParamsflH[which(ParamsflH$Year==2015),]
Paramf<-Paramf[-match(rmp2,Paramf$catchment),]
#Paramf<-Paramf[-which(is.na(ParamSpecial$epsilonGPD)),]
RlevErrf <- calculate_return_levels(Paramf,ci=1)
RlevErrf$year=2015

#ParamSpecial<-ParamSpecial[-which(is.na(ParamSpecial$epsilonGPD)),]
RlevErri <- calculate_return_levels(ParamSpecial,ci=1)
RlevErri$year=1955


if(haz=="Drought"){
  RlevErri$returnLevels= -RlevErri$returnLevels
  RlevErri$returnLevels[which(RlevErri$returnLevels<0)]=RlevErri$returnLevelErr[which(RlevErri$returnLevels<0)]
  
  RlevErrf$returnLevels= -RlevErrf$returnLevels
  RlevErrf$returnLevels[which(RlevErrf$returnLevels<0)]=RlevErrf$returnLevelErr[which(RlevErrf$returnLevels<0)]
}
mapu=na.omit((match(ParamSpecial$catchment,UpArea$outl2)))
ParamSpecial$upa=UpArea$upa[mapu]
RlevErri$catchment=ParamSpecial$catchment
RlevErri$epsilon=ParamSpecial$epsilonGPD
RlevErri$relErr=(RlevErri$returnLevelErr/RlevErri$returnLevels+1e-4)*100
RlevErri$rlf=RlevErrf$returnLevels
RlevErri$change=abs(RlevErri$rlf-RlevErri$returnLevels)
RlevErri$relcf=(RlevErri$rlf/RlevErri$returnLevels+1e-4)*100
RlevErri$ErrVsCh=RlevErri$rlf-RlevErri$returnLevelErr
RlevErri$S2n="Err > Change"
RlevErri$S2n[which(RlevErri$ErrVsCh>=0)]="Change > Err"
RlevErri$S2n[which(is.nan(RlevErri$ErrVsCh))]="Unstable"
length(which(RlevErri$S2n=="Err > Change"))
length(which(RlevErri$S2n=="Unstable"))/length(IRpoints$SCF)
ParamPlot=inner_join(RlevErri,UpArea,by=c("catchment"="outl2"))
max(RlevErri$returnLevelErr,na.rm=T)
Paraplot <- st_as_sf(ParamPlot, coords = c("Var1.x", "Var2.x"), crs = 4326)
Paraplot <- st_transform(Paraplot, crs = 3035)
#Paraplot=Paraplot[-which(is.na(Paraplot$returnLevels)),]

Unstable_locs=ParamPlot[which(ParamPlot$S2n=="Unstable"),]
Unstable_locs <- right_join(GHR_riv, Unstable_locs, by = c("outl2" = "catchment"))
#### [SPlot] - Supplement map - Relative of 10y-RL error in 1955 -----
br=seq(0,200,by=20)
labels=br
tsize=22
osize=16
colNA="transparent"
titleX=paste0("10-year RL error for ",haz)
paletS=c(hcl.colors(11, palette = "YlGnBu", alpha = NULL, rev = T, fixup = TRUE))
m1=ggplot(basemap) +
  geom_sf(fill="gray90",color="darkgrey",size=0.5)+
  geom_sf(data=Paraplot,aes(col=relErr,geometry=geometry,size=upa),alpha=1,stroke=0,shape=15)+ 
  
  scale_size(range = c(0.1, 0.5), trans="sqrt",name= expression(paste("Upstream area ", (km^2),
                                                                      sep = " ")),
             breaks=c(101,1000,10000,100000,500000), labels=c("100","1000", "10 000", "100 000", "500 000"),
             guide = "none")+
  coord_sf(xlim = c(min(nco[,1]),max(nco[,1])), ylim = c(min(nco[,2]),max(nco[,2])))+
  scale_color_gradientn(
    colors=paletS,
    breaks=br,labels=labels, limits=c(0,200),
    oob = scales::squish,na.value="white", name="RL(1955) RelErr (%)")   +
  labs(x="Longitude", y = "Latitude")+
  guides(colour = guide_colourbar(barheight = 16, barwidth = 1.5),
         fill = guide_colourbar(barheight = 16, barwidth = .6))+
  theme(axis.title=element_text(size=tsize),
        title = element_text(size=osize),
        axis.text=element_text(size=osize),
        panel.background = element_rect(fill = "aliceblue", colour = "grey1"),
        panel.border = element_rect(linetype = "solid", fill = NA, colour="black"),
        legend.text = element_text(size=osize),
        legend.title = element_text(size = osize, margin = margin(t = 2, r = 2, b = 6, l = 0)),
        legend.spacing.x = unit(0.2, "cm"),
        legend.position = "right",
        legend.box = "horizontal",  # Stack legends vertically
        panel.grid.major = element_line(colour = "grey70"),
        panel.grid.minor = element_line(colour = "grey90"),
        legend.key = element_rect(fill = "transparent", colour = "transparent"),
        legend.key.size = unit(1, "cm"))

ggsave(paste0("D:/tilloal/Documents/LFRuns_utils/TrendAnalysis/plots/Error_relative_",haz,".jpg"), m1, width=23, height=20, units=c("cm"),dpi=1000) 

#### [SPlot] - Supplement map - Comparison of Error and 2015-1955 changes -----
m2=ggplot(basemap) +
  geom_sf(fill="gray90",color="darkgrey",size=0.5)+
  geom_sf(data=Paraplot,aes(col=S2n,geometry=geometry,size=upa),alpha=1,stroke=0,shape=15)+ 
  
  scale_size(range = c(0.1, 0.5), trans="sqrt",name= expression(paste("Upstream area ", (km^2),
                                                                      sep = " ")),
             breaks=c(101,1000,10000,100000,500000), labels=c("100","1000", "10 000", "100 000", "500 000"),
             guide = "none")+
  coord_sf(xlim = c(min(nco[,1]),max(nco[,1])), ylim = c(min(nco[,2]),max(nco[,2])))+
  scale_color_manual(
    values=c("Unstable"="purple","Err > Change"="darkred","Change > Err"="royalblue"),
    name="Error classes")   +
  labs(x="Longitude", y = "Latitude")+
  guides(colour = guide_legend(override.aes = list(size = 10)))+
  # guides(colour = guide_colourbar(barheight = 16, barwidth = 1.5),
  #        fill = guide_colourbar(barheight = 16, barwidth = .6))+
  theme(axis.title=element_text(size=tsize),
        title = element_text(size=osize),
        axis.text=element_text(size=osize),
        panel.background = element_rect(fill = "aliceblue", colour = "grey1"),
        panel.border = element_rect(linetype = "solid", fill = NA, colour="black"),
        legend.text = element_text(size=osize),
        legend.title = element_text(size = osize, margin = margin(t = 2, r = 2, b = 6, l = 0)),
        legend.spacing.x = unit(0.2, "cm"),
        legend.position = "right",
        legend.box = "horizontal",  # Stack legends vertically
        panel.grid.major = element_line(colour = "grey70"),
        panel.grid.minor = element_line(colour = "grey90"),
        legend.key = element_rect(fill = "transparent", colour = "transparent"),
        legend.key.size = unit(1, "cm"))

ggsave(paste0("D:/tilloal/Documents/LFRuns_utils/TrendAnalysis/plots/Error_class_",haz,".jpg"), m2, width=23, height=20, units=c("cm"),dpi=1000) 



RlevErrip=inner_join(RlevErri,Shapeparf,by="catchment")



###5.1.1 Spatial smoothing for drought to remove noise from unstable GPD fits ----
if (haz=="Drought"){
  
  #Spatial smooting of LUAgg
  #tweak DataL to identify potentially problematic pixels
  IRpoints <- right_join(GHR_riv, IRpoints, by = c("llcoord" = "latlong"))
  matRLi=match(IRpoints$outl2,RLGPDflSCF$unikout)
  IRpoints$RLi=RLGPDflSCF$Y2015[matRLi]
  length(which(IRpoints$RLi>10))

  
  # Pixels on which the spatial smoothing will be applied on
  unikR=unique(HydroRsf$Id)
  DataRSmooth=c()
  for (rid in 1:length(unikR)){
    print(rid)
    regioid=unikR[rid]
    DataReI=IRpoints[which(IRpoints$HydroRegions_raster_WGS84==regioid),]
    DataRefn=DataReI[which(is.na(DataReI$Y2015)),]
    #spatial smoothing applied to pixels with unstable fits and pixels surrounding them
    
    DataReUn=Unstable_locs[which(Unstable_locs$HydroRegions_raster_WGS84==regioid),]
    matchuns=which(!is.na(match(DataReI$outl2,DataReUn$outl2)))
    DataRefix=DataReI[matchuns,]
    if (length(DataRefix$x>0)){
      # 
      points <- data.frame(x=DataReI$x, y=DataReI$y, z=DataReI$outl2)
      pointX<- data.frame(x=DataRefix$x, y=DataRefix$y, z=DataRefix$outl2)
      coordinates(points) <- ~x + y
      coordinates(pointX) <- ~x + y
      
      max_distance <- 5  # Define the maximum distance for weighting
      neighbours <- sapply(1:length(pointX), function(i) {
        point <- pointX[i, ]
        neighbour_finder(point, points,max_distance)
      })
      neighbours=unique(as.vector(unlist(neighbours)))
      matchd=which(!is.na(match(DataReI$outl2,neighbours)))
      DataRefix2=DataReI[matchd,]
      #remove large rivers
      if (length(which(DataRefix2$RLi>10))>0){
        DataRefix2=DataRefix2[-which(DataRefix2$RLi>10),]
      }
      #DataRefix2=rbind(DataRefix2,DataRefn)
      DataRSmooth=rbind(DataRSmooth,DataRefix2)
    }else{ print("no pixel to smooth")}
  }
  
  
  #plot pixel status
  mld=match(DataRSmooth$outl2,IRpoints$outl2)
  length(mld)/length(IRpoints$x)
  mli=which(IRpoints$gen_IR==2)
  mlo=match(rmp2,IRpoints$outl2)
  IRpoints$status="original"
  IRpoints$status[mld]="smoothed"
  IRpoints$status[mli]="IRES"
  IRpoints$status[mlo]="off-limit parameter"
  
  length(which(IRpoints$status=="original"))/length(IRpoints$status)+
    length(which(IRpoints$status=="IRES"))/length(IRpoints$status)+
    length(which(IRpoints$status=="smoothed"))/length(IRpoints$status)
  
  
  length(which(IRpoints$status=="off-limit parameter"))/length(IRpoints$status)*100
  Paraplot <- st_as_sf(IRpoints, coords = c("Var1.x", "Var2.x"), crs = 4326)
  Paraplot <- st_transform(Paraplot, crs = 3035)
  
  #### [SPlot] - Supplement map - Relative of 10y-RL error in 1955 -----
  br=seq(0,200,by=20)
  labels=br
  tsize=22
  osize=16
  colNA="transparent"
  titleX=paste0("10-year RL error for ",haz)
  paletS=c(hcl.colors(11, palette = "YlGnBu", alpha = NULL, rev = T, fixup = TRUE))
  m1=ggplot(basemap) +
    geom_sf(fill="gray90",color="darkgrey",size=0.5)+
    geom_sf(data=Paraplot,aes(geometry=geometry,size=upa,col=status),alpha=1,stroke=0,shape=15)+ 
    scale_color_manual(values=c("original"="royalblue","smoothed"="darkorange","IRES"="darkred","off-limit parameter"="black"),name=" ")+
    
    scale_size(range = c(0.1, 0.5), trans="sqrt",name= expression(paste("Upstream area ", (km^2),
                                                                        sep = " ")),
               breaks=c(101,1000,10000,100000,500000), labels=c("100","1000", "10 000", "100 000", "500 000"),
               guide = "none")+
    coord_sf(xlim = c(min(nco[,1]),max(nco[,1])), ylim = c(min(nco[,2]),max(nco[,2])))+
    guides(colour = guide_legend(override.aes = list(size = 10)))+
    labs(x="Longitude", y = "Latitude")+
    theme(axis.title=element_text(size=tsize),
          title = element_text(size=osize),
          axis.text=element_text(size=osize),
          panel.background = element_rect(fill = "aliceblue", colour = "grey1"),
          panel.border = element_rect(linetype = "solid", fill = NA, colour="black"),
          legend.text = element_text(size=osize),
          legend.title = element_text(size = osize, margin = margin(t = 2, r = 2, b = 6, l = 0)),
          legend.spacing.x = unit(0.2, "cm"),
          legend.position = "right",
          legend.box = "horizontal",  # Stack legends vertically
          panel.grid.major = element_line(colour = "grey70"),
          panel.grid.minor = element_line(colour = "grey90"),
          legend.key = element_rect(fill = "transparent", colour = "transparent"),
          legend.key.size = unit(1, "cm"))
  
  ggsave(paste0("D:/tilloal/Documents/LFRuns_utils/TrendAnalysis/plots/smoothfit_",haz,".jpg"), m1, width=23, height=20, units=c("cm"),dpi=1000) 
  length(DataRSmooth$x)
}


###[Plot] - Supplement - 10y RL in 1950 ----
plotmean=F
if (plotmean==T){
 
  #mean 10y RL over the period for each pixel
  #compute mean Return level
  #Socio CF
  data = data.frame(RLGPDflSCF)
  #Historical
  data2 = data.frame(RLGPDflH)
  
  mdat=rowMeans(data[,c(1:70)],na.rm=T)
  mdat2=rowMeans(data2[,c(1:70)],na.rm=T)
  
  mdat[(which(mdat==0))]=1
  hist(mdat,xlim=c(0,1000),breaks=10)
  
  RLev1=inner_join(RLGPDflH,UpArea,by=c("unikout"="outl2"))
  RLev1 <- st_as_sf(RLev1, coords = c("Var1.x", "Var2.x"), crs = 4326)
  RLev1 <- st_transform(RLev1, crs = 3035)
  RLev1$Y1951=RLev1$Y1951+1e-4
  RLev1$mq=mdat+1e-4
  
  rmp2=na.omit(unique(ParamsflSCF$catchment[rmpixs]))
  RLev1=RLev1[-match(rmp2,RLev1$unikout),]
 
   #compute specific discharge
  RLev1$Mqsp=RLev1$mq/RLev1$upa*1000
  min(RLev1$Mqsp,na.rm=T)
  max(RLev1$Mqsp,na.rm=T)
  plotmean=F
  
  ###[Plot] - Supplement - Mean 10-Y RL in specific discharge  ----
  br=c(1e-4,0.001,0.01,.1,1,10,100,1000)
  labels=br
  tsize=22
  osize=16
  colNA="transparent"
  titleX=paste0("10-year Return level for ",haz," in 1951")
  
  legend="Q (l/s/km2)"
  #legend="Q (m3/s)"
  palet=c(hcl.colors(11, palette = "YlGnBu", alpha = NULL, rev = T, fixup = TRUE))
  
  m0=ggplot(basemap) +
    geom_sf(fill="gray90",color="darkgrey",size=0.5)+
    geom_sf(data=RLev1,aes(col=Mqsp,geometry=geometry,size=upa),alpha=1,stroke=0,shape=15)+ 
    
    scale_size(range = c(0.1, 0.5), trans="sqrt",name= expression(paste("Upstream area ", (km^2),
                                                                        sep = " ")),
               breaks=c(101,1000,10000,100000,500000), labels=c("100","1000", "10 000", "100 000", "500 000"),
               guide = "none")+
    coord_sf(xlim = c(min(nco[,1]),max(nco[,1])), ylim = c(min(nco[,2]),max(nco[,2])))+
    scale_color_gradientn(
      colors=palet,
      breaks=br,trans="log", limits=c(1e-3,10),
      oob = scales::squish,na.value=colNA, name=na)   +
    labs(x="Longitude", y = "Latitude")+
    guides(colour = guide_colourbar(barheight = 16, barwidth = 1.5),
           fill = guide_colourbar(barheight = 16, barwidth = .6))+
    theme(axis.title=element_text(size=tsize),
          title = element_text(size=osize),
          axis.text=element_text(size=osize),
          panel.background = element_rect(fill = "aliceblue", colour = "grey1"),
          panel.border = element_rect(linetype = "solid", fill = NA, colour="black"),
          legend.text = element_text(size=osize),
          legend.title = element_text(size = osize, margin = margin(t = 2, r = 2, b = 6, l = 0)),
          legend.spacing.x = unit(0.2, "cm"),
          legend.position = "right",
          legend.box = "horizontal",  # Stack legends vertically
          panel.grid.major = element_line(colour = "grey70"),
          panel.grid.minor = element_line(colour = "grey90"),
          legend.key = element_rect(fill = "transparent", colour = "transparent"),
          legend.key.size = unit(1, "cm"))
  
  ggsave(paste0("D:/tilloal/Documents/LFRuns_utils/TrendAnalysis/plots/RL10_",haz,"_mean_scfxNew2.jpg"), m0, width=23, height=20, units=c("cm"),dpi=1000) 
 
  ####[SPlot] - Supplement - Bound on the GPD in 1955 -----
  # SocCF scenario
  paraU$bound=paraU$thresholdGPD-(paraU$sigmaGPD/paraU$epsilonGPD)
  paraU$bound[which(abs(paraU$epsilonGPD)<1e-3)]=NA
  hist(paraU$bound,xlim=c(-10,0),breaks=100000)
  paraU$bound[which(paraU$bound>0)]=0
  paraU$bound=-paraU$bound
  paraU$bound[which(paraU$bound<1e-4)]=1e-4
  paraU=as.data.frame(paraU)
  ParamPlot=inner_join(paraU,UpArea,by=c("catchment"="outl2"))
  Paraplot <- st_as_sf(ParamPlot, coords = c("Var1.x", "Var2.x"), crs = 4326)
  Paraplot <- st_transform(Paraplot, crs = 3035)
  
  br=c(1e-2,1e-1,1,10,100,1000)
  labels=c("0.01","0.1","1","10","100","1000")
  tsize=22
  osize=16
  colNA="transparent"
  titleX=paste0("Lower bound of fitted GPD for ",haz)
  paletS=c(hcl.colors(11, palette = "YlGnBu", alpha = NULL, rev = T, fixup = TRUE))
  
  m3=ggplot(basemap) +
    geom_sf(fill="gray90",color="darkgrey",size=0.5)+
    geom_sf(data=Paraplot,aes(col=bound,geometry=geometry,size=upa),alpha=1,stroke=0,shape=15)+ 
    
    scale_size(range = c(0.1, 0.5), trans="sqrt",name= expression(paste("Upstream area ", (km^2),
                                                                        sep = " ")),
               breaks=c(101,1000,10000,100000,500000), labels=c("100","1000", "10 000", "100 000", "500 000"),
               guide = "none")+
    coord_sf(xlim = c(min(nco[,1]),max(nco[,1])), ylim = c(min(nco[,2]),max(nco[,2])))+
    scale_color_gradientn(
      colors=paletS,
      breaks=br,labels=labels, limits=c(1e-2,1000),trans="log",
      oob = scales::squish,na.value="white", name="Q7min (m3/s)")   +
    labs(x="Longitude", y = "Latitude")+
    guides(colour = guide_colourbar(barheight = 16, barwidth = 1.5),
           fill = guide_colourbar(barheight = 16, barwidth = .6))+
    theme(axis.title=element_text(size=tsize),
          title = element_text(size=osize),
          axis.text=element_text(size=osize),
          panel.background = element_rect(fill = "aliceblue", colour = "grey1"),
          panel.border = element_rect(linetype = "solid", fill = NA, colour="black"),
          legend.text = element_text(size=osize),
          legend.title = element_text(size = osize, margin = margin(t = 2, r = 2, b = 6, l = 0)),
          legend.spacing.x = unit(0.2, "cm"),
          legend.position = "right",
          legend.box = "horizontal",  # Stack legends vertically
          panel.grid.major = element_line(colour = "grey70"),
          panel.grid.minor = element_line(colour = "grey90"),
          legend.key = element_rect(fill = "transparent", colour = "transparent"),
          legend.key.size = unit(1, "cm"))
  
  ggsave(paste0("D:/tilloal/Documents/LFRuns_utils/TrendAnalysis/plots/lowerbound_",haz,".jpg"), m3, width=23, height=20, units=c("cm"),dpi=1000) 
  
  ####[SPlot] - Supplement - Shape parameter instability check at regional level----
  #regional plot comparing shape parameters between runs
  regioid=149
  paraS=ParamsflSCF[which(ParamsflSCF$Year==2015),]
  matPar=match(GHR_riv$outl2,paraS$catchment)
  plot(matPar)
  paraS$hydroR=0
  paraS$hydroR[matPar]=GHR_riv$HydroRegions_raster_WGS84
  paraS1=paraS[which(paraS$hydroR==regioid),]
  
  paraRW=ParamsflRWCF[which(ParamsflSCF$Year==2015),]
  matPar=match(GHR_riv$outl2,paraRW$catchment)
  plot(matPar)
  paraRW$hydroR=0
  paraRW$hydroR[matPar]=GHR_riv$HydroRegions_raster_WGS84
  paraRW1=paraRW[which(paraRW$hydroR==regioid),]
  
  paraS1$epsilon_diff=paraRW1$epsilonGPD - paraS1$epsilonGPD
  paraS1$sigma_diff=paraRW1$sigmaGPD - paraS1$sigmaGPD
  paraS1$th_diff=paraRW1$thresholdGPD - paraS1$thresholdGPD
  ParamPlot=inner_join(paraS1,UpArea,by=c("catchment"="outl2"))
  
  
  cord.dec=ParamPlot[,c(22,23)]
  cord.dec = SpatialPoints(cord.dec, proj4string=CRS("+proj=longlat"))
  cord.UTM <- spTransform(cord.dec, CRS("+init=epsg:3035"))
  nco2=cord.UTM@coords
  Paraplot <- st_as_sf(ParamPlot, coords = c("Var1.x", "Var2.x"), crs = 4326)
  Paraplot <- st_transform(Paraplot, crs = 3035)
  
  br=c(seq(-0.5,0.5,0.1))
  labels=br
  tsize=22
  osize=16
  colNA="transparent"

  paletS=c(hcl.colors(11, palette = "RdYlBu", alpha = NULL, rev = T, fixup = TRUE))
  m4=ggplot(basemap) +
    geom_sf(fill="white",color="darkgrey",size=0.5)+
    geom_sf(data=Paraplot,aes(col=epsilon_diff,geometry=geometry,size=upa),alpha=1,stroke=0,shape=15)+ 
    
    scale_size(range = c(1, 5), trans="sqrt",name= expression(paste("Upstream area ", (km^2),
                                                                        sep = " ")),
               breaks=c(101,1000,10000,100000,500000), labels=c("100","1000", "10 000", "100 000", "500 000"),
               guide = "none")+
    coord_sf(xlim = c(min(nco2[,1]),max(nco2[,1])), ylim = c(min(nco2[,2]),max(nco2[,2])))+
    scale_color_gradientn(
      colors=paletS,
      breaks=br,labels=labels, limits=c(-.5,.5),
      oob = scales::squish,na.value="white", name="Parameter")   +
    labs(x="Longitude", y = "Latitude")+
    guides(colour = guide_colourbar(barheight = 16, barwidth = 1.5),
           fill = guide_colourbar(barheight = 16, barwidth = .6))+
    theme(axis.title=element_text(size=tsize),
          title = element_text(size=osize),
          axis.text=element_text(size=osize),
          panel.background = element_rect(fill = "aliceblue", colour = "grey1"),
          panel.border = element_rect(linetype = "solid", fill = NA, colour="black"),
          legend.text = element_text(size=osize),
          legend.title = element_text(size = osize, margin = margin(t = 2, r = 2, b = 6, l = 0)),
          legend.spacing.x = unit(0.2, "cm"),
          legend.position = "right",
          legend.box = "horizontal",  # Stack legends vertically
          panel.grid.major = element_line(colour = "grey70"),
          panel.grid.minor = element_line(colour = "grey90"),
          legend.key = element_rect(fill = "transparent", colour = "transparent"),
          legend.key.size = unit(1, "cm"))
  
  m4
  
  paraH=ParamsflH[which(ParamsflSCF$Year==2015),]
  paraRW=ParamsflRWCF[which(ParamsflSCF$Year==2015),]
  paraW=ParamsflWCF[which(ParamsflSCF$Year==2015),]
  
  #### [SPlot] - Supplement map - Shape parameter plot -----
  br=seq(-1,0, by=0.2)
  labels=br
  tsize=22
  osize=16
  colNA="transparent"
  titleX=paste0("Shape parameter for ",haz)
  paletS=c(hcl.colors(11, palette = "YlGnBu", alpha = NULL, rev = T, fixup = TRUE))
  
  m5=ggplot(basemap) +
    geom_sf(fill="gray90",color="darkgrey",size=0.5)+
    geom_sf(data=Paraplot,aes(col=epsilonGPD,geometry=geometry,size=upa),alpha=1,stroke=0,shape=15)+ 
    
    scale_size(range = c(0.1, 0.5), trans="sqrt",name= expression(paste("Upstream area ", (km^2),
                                                                        sep = " ")),
               breaks=c(101,1000,10000,100000,500000), labels=c("100","1000", "10 000", "100 000", "500 000"),
               guide = "none")+
    coord_sf(xlim = c(min(nco[,1]),max(nco[,1])), ylim = c(min(nco[,2]),max(nco[,2])))+
    scale_color_gradientn(
      colors=paletS,
      breaks=br, limits=c(-1,0),
      oob = scales::squish,na.value=colNA, name="shape")   +
    labs(x="Longitude", y = "Latitude")+
    guides(colour = guide_colourbar(barheight = 16, barwidth = 1.5),
           fill = guide_colourbar(barheight = 16, barwidth = .6))+
    theme(axis.title=element_text(size=tsize),
          title = element_text(size=osize),
          axis.text=element_text(size=osize),
          panel.background = element_rect(fill = "aliceblue", colour = "grey1"),
          panel.border = element_rect(linetype = "solid", fill = NA, colour="black"),
          legend.text = element_text(size=osize),
          legend.title = element_text(size = osize, margin = margin(t = 2, r = 2, b = 6, l = 0)),
          legend.spacing.x = unit(0.2, "cm"),
          legend.position = "right",
          legend.box = "horizontal",  # Stack legends vertically
          panel.grid.major = element_line(colour = "grey70"),
          panel.grid.minor = element_line(colour = "grey90"),
          legend.key = element_rect(fill = "transparent", colour = "transparent"),
          legend.key.size = unit(1, "cm"))
  
  ggsave(paste0("D:/tilloal/Documents/LFRuns_utils/TrendAnalysis/plots/ShapePar_",haz,".jpg"), m5, width=23, height=20, units=c("cm"),dpi=1000) 
  
}  

#4 Change attribution-----
RLGPDflSCF=as.data.frame(RLGPDflSCF)
RLGPDflH=as.data.frame(RLGPDflH)
RLGPDflWCF=as.data.frame(RLGPDflWCF)
RLGPDflRWCF=as.data.frame(RLGPDflRWCF)

#Before change attribution, I compute total change to compare with final results
#mato=match(RLGPDflSCF[,71],RLGPDflH[,71])
colnames(RLGPDflSCF)
cd=as.numeric(which(colnames(RLGPDflSCF)=="Y1955"))

##4.1 change from climate--------------------
Climtrend=RLGPDflSCF[,c(1:70)]-(RLGPDflSCF[,cd])

##4.2 change from socoeconomy -----------

###4.2.1 land use -----------
Soctrend=(RLGPDflRWCF[,-71]-RLGPDflSCF[,-71])


###4.2.2 reservoir trend ----------
Restrend=(RLGPDflWCF[,-71]-RLGPDflRWCF[,-71])

###4.2.3 water demand trend ----------
Wutrend=RLGPDflH[,-71]-RLGPDflWCF[,-71]

##4.4 Total trend -----
#Totaltrend=Climtrend+Soctrend+Wutrend+Restrend
Totaltrend=RLGPDflH[,c(1:70)]-(RLGPDflSCF[,cd])

#Upstream area vector
UpAvec=UpArea[,c(12,3)]

##4.5 remove pixels with crazy shape parameter --------
#Aggregation by Regions of initial RL
data = data.frame(RLGPDflSCF)
data=data[-match(rmp2,data$unikout),]


##4.6 Spatial aggregation to desired regions: HydroRegion ----------


# DataI=right_join(GHR_riv,data,by = c("outl2"="unikout"))
# HRM=match(DataI$HydroRegions_raster_WGS84,HydroRsf$Id)
# DataI$Hid=HydroRsf$Id[HRM]
# DataI$mdatS=DataI$mdat * 1000 / (DataI$upa)
# DataI$mdatmaxS=DataI$mdatmax * 1000 / (DataI$upa)
# mean(DataI$mdatS,na.rm=T)
# RegioRLi=aggregate(list(RL10=DataI$mdat),
#                    by = list(HydroR=DataI$Hid),
#                    FUN = function(x) c(mean=mean(x,na.rm=T),dev=sd(x,na.rm=T),len=length(x),med=median(x,na.rm=T),q1=quantile(x, 0.05, na.rm=T),q3=quantile(x, 0.95, na.rm=T)))
# RegioRLi <- do.call(data.frame, RegioRLi)

#matching biogeoregions and hybas07
bio_hybas=inner_join(biogeo_rivers,GNF,by=c("outl2"))
HRM=na.omit(unique((match(bio_hybas$HYBAS_ID,hybasf7$HYBAS_ID))))
hybasf7_dom=hybasf7[HRM,]

bhp_m=na.omit(unique(((match(hybasf7_dom$HYBAS_ID,bio_hybas$HYBAS_ID)))))
hybasf7_dom$biogeoreg=bio_hybas$code[bhp_m]

#5 Plots & Aggregation by Regions -----
herziz=na.omit(match(GHR_riv$HydroRegions_raster_WGS84,GHshpp$Id))
GHR_riv$HER=GHshpp$CODEB[herziz]

#Postprocessing of Data for water demand and land use change
### parameter definition
yrname=colnames(Climtrend)
unikout=RLGPDflSCF[, 71]
parameters=data.frame(ParamsflSCF)
rmpixels=rmpixs
Regio=hybasf7
Regio=HydroRsf
Drivertrend=Climtrend

#All changes relative to the same total change between 1955 socCF
# and Historical

DataI=data.frame(RLGPDflH)
c1=as.numeric(which(colnames(DataI)=="Y1955"))
irange=match(yrname,colnames(DataI))
c2=irange
#symetric difference
DataI$Init=RLGPDflSCF[,c1]
DataI[,irange]=(DataI$Init+DataI[,irange])/2




##5.1 Land use change data processing -----------
change="socio"
DataL=ComputeChange(Drivertrend=Soctrend, unikout, DataI,
                      outhybas07,parameters, rmpixels, UpAvec, GHR_riv, HydroRsf, 
                    yrname,change="socio",eps=0.1)
skw=which.min(DataL$Y1985)
yrange=match(yrname,colnames(DataL))

palet=c(hcl.colors(11, palette = "RdYlBu", alpha = NULL, rev = T, fixup = TRUE))
###5.1.1 Spatial smoothing for drought to remove noise from unstable GPD fits ----
if (haz=="Drought"){
  Change_smoothed=c()
  loc_smoothed=c()
  batar=c()
  fdp=c()
  unikR=unique(HydroRsf$Id)
  for (rid in 1:length(unikR)){
    print(rid)
    regioid=unikR[rid]
    unique(DataL$HydroRegions_raster_WGS84)
    DataReI=DataL[which(DataL$Regio_id==regioid),]
    DataRefix=DataRSmooth[which(DataRSmooth$HydroRegions_raster_WGS84==regioid),]
    dm=match(DataRefix$outl2,DataReI$outl2)
    DataRefix$ch2015=DataReI$Y2015[dm]
    #I have to avoid the fitting on  pixels with low changes
    DataRefix=DataRefix[which(abs(DataRefix$ch2015)>0.5 ),]
    #DataRetiny=DataReI[which(DataReI$RLi<.1),]
    #DataRefn=DataReI[which(is.na(DataReI$Y2015)),]
    #different rules to identify unstable pixels on which a spatial smoothing will 
    #be applied
    # DataRefix=DataReI[which(DataReI$epsD>(0.1)),]
    # DataRefix=DataReI[which( DataReI$RLi<10),]
    # DataRefix=DataRefix[which(DataRefix$epsD>(0.1) | DataRefix$RLi<1),]
    # DataRefix=DataRefix[which(abs(DataRefix$Y2015)>0.5 ),]
    # DataRefix=rbind(DataRefix,DataRefn)
    # 
    # DataReUn=Unstable_locs[which(Unstable_locs$HydroRegions_raster_WGS84==regioid),]
    # matchuns=which(!is.na(match(DataReI$outlets,DataReUn$outlets.x)))
    # DataRefix=DataReI[matchuns,]
    # DataRefix=rbind(DataRefix,DataRefn)
    # 
    # points <- data.frame(x=DataReI$x, y=DataReI$y, z=DataReI$outlets)
    # pointX<- data.frame(x=DataRefix$x, y=DataRefix$y, z=DataRefix$outlets)
    # coordinates(points) <- ~x + y
    # coordinates(pointX) <- ~x + y
    # 
    # neighbours=unique(as.vector(unlist(neighbours)))
    # matchd=which(!is.na(match(DataReI$outlets,neighbours)))
    # DataRefix=DataReI[matchd,]
    # DataRefix=rbind(DataRefix,DataRefn)
      
    if (length(DataRefix$x)>1){
      points <- data.frame(x=DataReI$x, y=DataReI$y, z=DataReI$Y2015)
      pointX<- data.frame(x=DataRefix$x, y=DataRefix$y, z=DataRefix$ch2015)
      # rivmask=data.frame(x=DataReI$x, y=DataReI$y, llcoord=DataReI$llcoord)
      coordinates(points) <- ~x + y
      coordinates(pointX) <- ~x + y
      # Apply the weighted average function to each point
      
      max_distance <- 15  # Define the maximum distance for weighting
      smoothed_values <- sapply(1:length(pointX), function(i) {
        point <- pointX[i, ]
        weighted_average(r, point, points, rivermask,max_distance)
        
      })
      
      rlocs=match(DataRefix$llcoord,DataReI$llcoord)
      DataReI$Y2015b=DataReI$Y2015
      DataReI$Y2015b[rlocs]=smoothed_values
      Change_smoothed=c(Change_smoothed,smoothed_values)
      loc_smoothed=c(loc_smoothed, DataRefix$llcoord)
      batar=c(batar,DataRefix$ch2015)
      fdp=c(fdp,DataRefix$HydroRegions_raster_WGS84)
    }
    
    # ggplot()+
    #   geom_point(data=as.data.frame(DataReI),aes(x=x,y=y, col=epsD))+
    #   scale_color_gradientn(colors=palet,limits=c(-.1,.1))+
    #   geom_point(data=DataRefix,aes(x=x,y=y),col="black")
    # 
    # 
    # ggplot(as.data.frame(DataReI),aes(x=x,y=y, col=Y2015b))+
    #   geom_point()+
    #   scale_color_gradientn(colors=palet,limits=c(-10,10),oob=scales::squish)
    # 
    # ggplot(as.data.frame(DataReI),aes(x=x,y=y, col=Y2015))+
    #   geom_point()+
    #   scale_color_gradientn(colors=palet,limits=c(-10,10),oob=scales::squish)
    # #

  }
  
  #DataL=DataLs
  plot(batar,Change_smoothed)
  mean(batar,na.rm=T)
  cmL=Change_smoothed
  btL=batar

  lmL=loc_smoothed
  hist(cmL)
  ratio=cmL/btL
  ratio[which(is.nan(ratio))]=1
  ratio[which(is.infinite(ratio))]=1
  ratio[which(is.na(btL))]=NA
  plot(ratio)
  rlocs=match(lmL,DataL$llcoord)
  DataL$Y2015b=DataL$Y2015
  DataL$Y2015b[rlocs]=cmL
  DataL2=DataL
  DataL2[rlocs,yrange]=DataL[rlocs,yrange]*ratio
  
  #verification of consistency of the averages
  DataLs=DataL
  rlocp=which(is.na(DataL2$Y2015[rlocs]))
  rlocx=rlocs
  
  #initial mean change in affected pixels
  print(mean(DataLs$Y2015[rlocx],na.rm=T))
  #interpolated mean change
  print(mean(DataLs$Y2015b[rlocx],na.rm=T))
  #corrected mean change
  print(mean(DataL2$Y2015[rlocx],na.rm=T))
  DataL=DataL2
  DataL=DataL[,-c(85:89)]
}

###5.1.2 Aggregation of relative changes ----
pointagg <- aggregate(list(Rchange_rel = DataL$Y2015),
                      by = list(HydroR = DataL$HER),
                      FUN = function(x) c(mean = mean(x, na.rm = TRUE), 
                                          dev = sd(x, na.rm = TRUE), 
                                          len = length(x), 
                                          med = median(x, na.rm = TRUE), 
                                          q1 = quantile(x, 0.025, na.rm = TRUE), 
                                          q3 = quantile(x, 0.975, na.rm = TRUE)))
pointSoc <- do.call(data.frame, pointagg)
head(pointSoc)

trendagg <- aggregate(list(Rchange = DataL[, yrange]),
                      by = list(HydroR = DataL$HER),
                      FUN = function(x) c(mean = mean(x, na.rm = TRUE)))
trendSoc <- do.call(data.frame, trendagg)


##5.2 Water demand change data processing ------
DataW=ComputeChange(Drivertrend=Wutrend, unikout,DataI,
                    outhybas07,parameters, rmpixels, UpAvec, GHR_riv, HydroRsf, 
                    yrname,change="socio",eps=0.1)

###5.2.1 Spatial smoothing for drought to remove noise from unstable GPD fits ----
if (haz=="Drought"){
  Change_smoothed=c()
  loc_smoothed=c()
  batar=c()
  fdp=c()
  for (rid in 1:length(unikR)){
    print(rid)
    regioid=unikR[rid]
    unique(DataW$HydroRegions_raster_WGS84)
    DataReI=DataW[which(DataW$Regio_id==regioid),]
    DataRefix=DataRSmooth[which(DataRSmooth$HydroRegions_raster_WGS84==regioid),]
    dm=match(DataRefix$outl2,DataReI$outl2)
    DataRefix$ch2015=DataReI$Y2015[dm]
    #I have to avoid the fitting on  pixels with low changes
    DataRefix=DataRefix[which(abs(DataRefix$ch2015)>0.5 ),]
    if (length(DataRefix$x)>1){
      points <- data.frame(x=DataReI$x, y=DataReI$y, z=DataReI$Y2015)
      pointX<- data.frame(x=DataRefix$x, y=DataRefix$y, z=DataRefix$ch2015)
      # rivmask=data.frame(x=DataReI$x, y=DataReI$y, llcoord=DataReI$llcoord)
      coordinates(points) <- ~x + y
      coordinates(pointX) <- ~x + y
      # Apply the weighted average function to each point
      
      max_distance <- 15  # Define the maximum distance for weighting
      smoothed_values <- sapply(1:length(pointX), function(i) {
        point <- pointX[i, ]
        weighted_average(r, point, points, rivermask,max_distance)
        
      })
      
      rlocs=match(DataRefix$llcoord,DataReI$llcoord)
      DataReI$Y2015b=DataReI$Y2015
      DataReI$Y2015b[rlocs]=smoothed_values
      Change_smoothed=c(Change_smoothed,smoothed_values)
      loc_smoothed=c(loc_smoothed, DataRefix$llcoord)
      batar=c(batar,DataRefix$ch2015)
      fdp=c(fdp,DataRefix$HydroRegions_raster_WGS84)
    }
    
    # ggplot()+
    #   geom_point(data=as.data.frame(DataReI),aes(x=x,y=y, col=epsD))+
    #   scale_color_gradientn(colors=palet,limits=c(-.1,.1))+
    #   geom_point(data=DataRefix,aes(x=x,y=y),col="black")
    # 
    # 
    # ggplot(as.data.frame(DataReI),aes(x=x,y=y, col=Y2015b))+
    #   geom_point()+
    #   scale_color_gradientn(colors=palet,limits=c(-10,10),oob=scales::squish)
    # 
    # ggplot(as.data.frame(DataReI),aes(x=x,y=y, col=Y2015))+
    #   geom_point()+
    #   scale_color_gradientn(colors=palet,limits=c(-10,10),oob=scales::squish)
    # #
    
  }
  
  btW=batar
  cmW=Change_smoothed
  lmW=loc_smoothed
  ratio=cmW/btW
  ratio[which(is.nan(ratio))]=1
  ratio[which(is.infinite(ratio))]=1
  ratio[which(is.na(btL))]=NA
  rlocs=match(lmW,DataW$llcoord)
  DataW$Y2015b=DataW$Y2015
  DataW$Y2015b[rlocs]=cmW
  DataW2=DataW
  DataW2[rlocs,yrange]=DataW2[rlocs,yrange]*ratio
  DataWs=DataW
  #rlocp=which(is.na(DataW2$Y2015[rlocs]))
  rlocx=rlocs
  #initial mean change in affected pixels

  # aa= which(is.nan(DataW$Y2015))
  # DataW[aa,yrange]=NA
  print(mean(DataW2$Y2015[rlocx],na.rm=T))
  print(mean(DataWs$Y2015[rlocx],na.rm=T))
  #interpolated mean change
  print(mean(DataWs$Y2015b[rlocx],na.rm=T))
  #corrected mean change
  print(mean(DataW2$Y2015[rlocx],na.rm=T))

  DataW2=DataW2[,-c(85:89)]
  DataW=DataW2
}

###5.2.2 Aggregation of relative changes ----
pointagg <- aggregate(list(Rchange_rel = DataW$Y2015),
                      by = list(HydroR = DataW$HER),
                      FUN = function(x) c(mean = mean(x, na.rm = TRUE), 
                                          dev = sd(x, na.rm = TRUE), 
                                          len = length(x), 
                                          med = median(x, na.rm = TRUE), 
                                          q1 = quantile(x, 0.025, na.rm = TRUE), 
                                          q3 = quantile(x, 0.975, na.rm = TRUE)))
pointWu <- do.call(data.frame, pointagg)

yrange=match(yrname,colnames(DataL))


trendagg <- aggregate(list(Rchange = DataW[, yrange]),
                      by = list(HydroR = DataW$HER),
                      FUN = function(x) c(mean = mean(x, na.rm = TRUE)))
trendWu <- do.call(data.frame, trendagg)


##5.3 Total change----

change="total"
DataT=ComputeChange(Drivertrend=Totaltrend, unikout, DataI,
                    outhybas07,parameters, rmpixels, UpAvec, GHR_riv, HydroRsf,
                    yrname,change,eps=0.1)

###5.3.1 Spatial smoothing for drought to remove noise from unstable GPD fits ----
if (haz=="Drought"){
  matRLi=match(DataT$outl2,RLGPDflWCF$unikout)
  DataT$RLi=RLGPDflH$Y2015[matRLi]
  Param1=ParamsflSCF[which(ParamsflSCF$Year==2015),]
  Param2=ParamsflH[which(ParamsflH$Year==2015),]
  matP=match(DataT$outl2,Param1$catchment)
  DataT$epsilon1=Param1$epsilonGPD[matP]
  DataT$epsilon2=Param2$epsilonGPD[matP]
  DataT$epsD=abs(DataT$epsilon1-DataT$epsilon2)
  median(DataT$RLi,na.rm=T)
  length(which(DataT$RLi<10 & (DataT$epsD)>(0.01) & abs(DataT$Y2015)>0.5))
  length(which(DataT$RLi<10 & (DataT$epsilon1)<(-.51) & abs(DataT$Y2015)>0.5))
  length(which(DataT$upa<10000 & abs(DataT$Y2015)>1))
  
  Change_smoothed=c()
  loc_smoothed=c()
  batar=c()
  fdp=c()
  for (rid in 1:length(unikR)){
    print(rid)
    regioid=unikR[rid]
    unique(DataT$HydroRegions_raster_WGS84)
    DataReI=DataT[which(DataT$Regio_id==regioid),]
    DataRefix=DataRSmooth[which(DataRSmooth$HydroRegions_raster_WGS84==regioid),]
    dm=match(DataRefix$outl2,DataReI$outl2)
    DataRefix$ch2015=DataReI$Y2015[dm]
    #I have to avoid the fitting on  pixels with low changes
    DataRefix=DataRefix[which(abs(DataRefix$ch2015)>0.5 ),]
    if (length(DataRefix$x)>1){
      points <- data.frame(x=DataReI$x, y=DataReI$y, z=DataReI$Y2015)
      pointX<- data.frame(x=DataRefix$x, y=DataRefix$y, z=DataRefix$ch2015)
      # rivmask=data.frame(x=DataReI$x, y=DataReI$y, llcoord=DataReI$llcoord)
      coordinates(points) <- ~x + y
      coordinates(pointX) <- ~x + y
      # Apply the weighted average function to each point
      
      max_distance <- 15  # Define the maximum distance for weighting
      smoothed_values <- sapply(1:length(pointX), function(i) {
        point <- pointX[i, ]
        weighted_average(r, point, points, rivermask,max_distance)
        
      })
      
      rlocs=match(DataRefix$llcoord,DataReI$llcoord)
      DataReI$Y2015b=DataReI$Y2015
      DataReI$Y2015b[rlocs]=smoothed_values
      Change_smoothed=c(Change_smoothed,smoothed_values)
      loc_smoothed=c(loc_smoothed, DataRefix$llcoord)
      batar=c(batar,DataRefix$ch2015)
      fdp=c(fdp,DataRefix$HydroRegions_raster_WGS84)
    }
    
    # ggplot()+
    #   geom_point(data=as.data.frame(DataReI),aes(x=x,y=y, col=epsD))+
    #   scale_color_gradientn(colors=palet,limits=c(-.1,.1))+
    #   geom_point(data=DataRefix,aes(x=x,y=y),col="black")
    # 
    # 
    # ggplot(as.data.frame(DataReI),aes(x=x,y=y, col=Y2015b))+
    #   geom_point()+
    #   scale_color_gradientn(colors=palet,limits=c(-10,10),oob=scales::squish)
    # 
    # ggplot(as.data.frame(DataReI),aes(x=x,y=y, col=Y2015))+
    #   geom_point()+
    #   scale_color_gradientn(colors=palet,limits=c(-10,10),oob=scales::squish)
    # #
    
  }
  btW=batar
  cmW=Change_smoothed
  lmW=loc_smoothed
  plot(btW,cmW)
  ratio=cmW/btW
  ratio[which(is.nan(ratio))]=1
  ratio[which(is.infinite(ratio))]=1
  ratio[which(is.na(btL))]=NA
  rlocs=match(lmW,DataT$llcoord)
  DataTs=DataT
  rlocp=which(is.na(DataW2$Y2015[rlocs]))
  rlocx=rlocs[-rlocp]
  DataT$Y2015b=DataT$Y2015
  DataT$Y2015b[rlocs]=cmW
  DataT2=DataT
  DataT2[rlocs,yrange]=DataT2[rlocs,yrange]*ratio
  
  #initial mean change in affected pixels
  print(mean(DataT$Y2015[rlocx],na.rm=T))
  #interpolated mean change
  print(mean(DataT$Y2015b[rlocx],na.rm=T))
  #corrected mean change
  print(mean(DataT2$Y2015[rlocx],na.rm=T))
  
  DataT2=DataT2[,-c(85:89)]
  DataT=DataT2
}

###5.3.2 Aggregation of relative changes ----
pointagg <- aggregate(list(Rchange_rel = DataT$Y2015),
                      by = list(HydroR = DataT$HER),
                      FUN = function(x) c(mean = mean(x, na.rm = TRUE), 
                                          dev = sd(x, na.rm = TRUE), 
                                          len = length(x), 
                                          med = median(x, na.rm = TRUE), 
                                          q1 = quantile(x, 0.025, na.rm = TRUE), 
                                          q3 = quantile(x, 0.975, na.rm = TRUE)))
pointTot <- do.call(data.frame, pointagg)

yrange=match(yrname,colnames(DataT))

trendagg <- aggregate(list(Rchange = DataT[, yrange]),
                      by = list(HydroR = DataT$HER),
                      FUN = function(x) c(mean = mean(x, na.rm = TRUE)))
trendTot <- do.call(data.frame, trendagg)


##5.4 Clim change----
change="clim"
DataC=ComputeChange(Drivertrend=Climtrend, unikout, DataI,
                    outhybas07,parameters, rmpixels, UpAvec, GHR_riv, HydroRsf, 
                    yrname,change, eps=0.1)


###5.4.1 Aggregation of relative changes ----
pointagg <- aggregate(list(Rchange_rel = DataC$Y2015),
                      by = list(HydroR = DataC$HER),
                      FUN = function(x) c(mean = mean(x, na.rm = TRUE), 
                                          dev = sd(x, na.rm = TRUE), 
                                          len = length(x), 
                                          med = median(x, na.rm = TRUE), 
                                          q1 = quantile(x, 0.025, na.rm = TRUE), 
                                          q3 = quantile(x, 0.975, na.rm = TRUE)))
pointClim <- do.call(data.frame, pointagg)

yrange=match(yrname,colnames(DataC))

trendagg <- aggregate(list(Rchange = DataC[, yrange]),
                      by = list(HydroR = DataC$HER),
                      FUN = function(x) c(mean = mean(x, na.rm = TRUE)))
trendClim <- do.call(data.frame, trendagg)


## 5.5 Reservoir change------

DataR=ComputeChange(Drivertrend=Restrend, unikout, DataI, outhybas07,
                    parameters, rmpixels, UpAvec, GHR_riv, HydroRsf,
                    yrname,change="socio",eps=0.1)

###5.5.1 Aggregation of relative changes ----
pointagg <- aggregate(list(Rchange_rel = DataR$Y2015),
                      by = list(HydroR = DataR$HER),
                      FUN = function(x) c(mean = mean(x, na.rm = TRUE), 
                                          dev = sd(x, na.rm = TRUE), 
                                          len = length(x), 
                                          med = median(x, na.rm = TRUE), 
                                          q1 = quantile(x, 0.025, na.rm = TRUE), 
                                          q3 = quantile(x, 0.975, na.rm = TRUE)))
pointRes <- do.call(data.frame, pointagg)

yrange=match(yrname,colnames(DataC))

trendagg <- aggregate(list(Rchange = DataR[, yrange]),
                      by = list(HydroR = DataR$HER),
                      FUN = function(x) c(mean = mean(x, na.rm = TRUE)))
trendRes <- do.call(data.frame, trendagg)


###5.5.2 Correction of ResTrend  ----
# removal of large change in pixels non-affected by reservoirs 
# connected to parameter instability
RestrendCor <- data.frame(Restrend, unikout = RLGPDflH$unikout)
r_path = paste0(hydroDir,'/reservoirs/')
outletname="res_ratio_diff_2020-1951.nc"
new_reservoirs=ReservoirOpen(r_path,outletname,outf)
inf_reservoirs=new_reservoirs[which(new_reservoirs$upa>0),]
length(inf_reservoirs$Var1.x)/length(DataC$x)
rmat=match(inf_reservoirs$outl2,DataR$outl2)
Rcrap=DataR[-rmat,]

#percentage of pixels with dicarded changes
length(which(abs(Rcrap$Y2015)>2))/length(Rcrap$Y2015)*100
Rcrap=Rcrap[which(abs(Rcrap$Y2015)>2),]
Rcrap=Rcrap[-which(is.na(Rcrap$x)),]


##5.6 Extra analysis ----

#save the data 
savedat=FALSE
if (savedat==TRUE){
 DataSave=list("Total"=DataT,"Climate"=DataC,"LandUse"=DataL,"Reservoirs"=DataR,"WaterDemand"=DataW)
 save(DataSave,file=paste0(hydroDir,"/TSEVA/output_plots/Drought_pixChange_v4.Rdata"))
}

#Suorva Dam example
suorva=5701511
SuorvaC=DataC[which(DataC$outl2==suorva),]
SuorvaL=DataL[which(DataL$outl2==suorva),]
SuorvaR=DataR[which(DataR$outl2==suorva),]
SuorvaC=DataC[which(DataC$outl2==suorva),]



#match dataT with Biogeoregions
period=c(1955,2015)
br=seq(-100,100,by=25)
labels=br
colNA="transparent"
palet=c(hcl.colors(11, palette = "RdYlBu", alpha = NULL, rev = F, fixup = TRUE))
title=paste0("Change in 10-years ",haz," Return Level between ", period[1], " and ", period[2])
legend="Change (%)"
pointP=pointClim
pointP=pointP[order(pointP$Rchange_rel.mean),]
pointP$id=c(1:length(pointP$Rchange_rel.mean))

#mean of all HER
mean(pointP$Rchange_rel.mean)
mean(pointP$Rchange_rel.dev)

####Question: How many catchment see an increase with Climate?-----

#Increase: If at least 90% of rivers have an increase in the catchment
climI=which(pointP$Rchange_rel.q1.2.5.>=0)
length(climI)/length(pointP$HydroR)


climT=which(pointP$Rchange_rel.mean<0)
length(climT)/length(pointP$HydroR)

#pixel level
climIp=length(which(DataC$Y2015>=15))/length(DataC$HydroRegions_raster_WGS84)
climIp

climS=which(pointP$Rchange_rel.q1.2.5.<=0 & pointP$Rchange_rel.q3.97.5.>=0)
length(climS)/length(pointP$HydroR)

climD=which(pointP$Rchange_rel.q3.97.5.<=0)
length(climD)/length(pointP$HydroR)



##5.7 Main drivers at different levels ------

###5.7.1 Main driver at HER level----
PointAI=data.frame(clim=abs(pointClim$Rchange_rel.mean),lu=abs(pointSoc$Rchange_rel.mean),
                   res=abs(pointRes$Rchange_rel.mean),wu=abs(pointWu$Rchange_rel.mean))


max_col_numbers <- apply(PointAI, 1, function(x) which.max(x))
max_col_numbers<-as.numeric(max_col_numbers)

colorz = c("4"="limegreen","3" ='tomato4',"2" ='orange',"1" ='royalblue')


pointap=full_join(GHshpp,pointClim,by=c("CODEB"="HydroR"))
st_geometry(pointap)<-NULL

cmat=match(Regio$Id,pointClim$HydroR)
pointap=pointap[which(!is.na(cmat)),]
cmat=cmat[which(!is.na(cmat))]
pointap$maxcol=max_col_numbers[cmat]
pointap=pointap[which(!is.na(pointap$maxcol)),]

pointplot=left_join(HydroRsf,pointap,by= c("CODEB"="Id"))
pointplot=pointplot[which(!is.na(pointplot$maxcol)),]
period=c(1951,2020)

nutplot <- st_transform(pointplot, crs = 3035)
br=c(-50,-20,-10,0,10,20,50)
labels=br
limi=c(-50,50)

#### [Plot] - Supplement map - Largest driver of change at the HER level----
pl1=ggplot(basemap) +
  geom_sf(fill="gray95")+
  geom_sf(fill=NA, color="grey") +
  geom_sf(data=nutplot,aes(fill=factor(maxcol),geometry=geometry),color="transparent",alpha=1,size=0.25,stroke=0,shape=15)+ 
  coord_sf(xlim = c(min(nco[,1]),max(nco[,1])), ylim = c(min(nco[,2]),max(nco[,2])))+
  scale_fill_manual(values = colorz, name=" ") +
  labs(x="Longitude", y = "Latitude")+
  guides(color = guide_colourbar(barwidth = 12, barheight = 1))+
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
pl1

ggsave(paste0("D:/tilloal/Documents/LFRuns_utils/TrendAnalysis/plots/ShapePar_",haz,".jpg"), pl1, width=23, height=20, units=c("cm"),dpi=1000) 

###Main driver at pixel level----

#pixel level drivers 
idc=match(c("Y2020","upa"),colnames(DataC))
idl=c(3,idc)
DataC1=DataC[,idl]
DataL1=DataL[,idl]
DataW1=DataW[,idl]
DataR1=DataR[,idl]


DataAI=data.frame(clim=abs(DataC1$Y2020),lu=abs(DataL1$Y2020),
                   resw=abs(DataR1$Y2020),wu=abs(DataW1$Y2020))

max_col_numbers <- apply(DataAI, 1, function(x) which.max(x))
max_col_numbers=as.numeric(max_col_numbers)

#### [Plot] - Supplement map - Largest driver of change at the pixel level-----
datap=DataC
datap$maxcol=max_col_numbers
datap=datap[-which(is.na(datap$Var1)),]
points <- st_as_sf(datap, coords = c("Var1", "Var2"), crs = 4326)
points <- st_transform(points, crs = 3035)
points=points[-which(is.na(points$maxcol)),]
points[which(points$maxcol==4),]
if (length(which(is.na(nutplot$maxcol)))>0) nutplot=nutplot[-which(is.na(nutplot$maxcol)),]
lab1=c("Climate","Land use","Reservoirs","Water demand")
tsize=14
osize=12
pl2<-ggplot(basemap) +
  geom_sf(fill="gray95",color="gray10",size=0.5)+
  geom_sf(data=nutplot,aes(fill=factor(maxcol),geometry=geometry),color="transparent",alpha=.6,size=0.25,stroke=0,shape=15)+ 
  geom_sf(data=points,aes(col=factor(maxcol),geometry=geometry,size=upa),alpha=.9,stroke=0,shape=15)+
  coord_sf(xlim = c(min(nco[,1]),max(nco[,1])), ylim = c(min(nco[,2]),max(nco[,2])))+
  scale_colour_manual(values = colorz, name="Largest change driver", labels=lab1) +
  scale_fill_manual(values = colorz, name="Largest change driver",labels=lab1) +
  scale_size(range = c(0.08, 0.4), trans="sqrt",name= expression(paste("Upstream area ", (km^2),
                                                                  sep = " ")),
             breaks=c(101,1000,10000,100000,500000), labels=c("100","1000", "10 000", "100 000", "500 000"),
             guide = "none")+
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

ggsave(paste0("D:/tilloal/Documents/LFRuns_utils/TrendAnalysis/plots/ShapePar_",haz,".jpg"), pl2, width=23, height=20, units=c("cm"),dpi=1000) 

#count pixel for each dominant driver 

#climate
length(which(points$maxcol==1))/length(points$x)
length(which(nutplot$maxcol==1))/length(nutplot$Id)

#landuse
length(which(points$maxcol==2))/length(points$x)
length(which(nutplot$maxcol==2))/length(nutplot$Id)

#reservoirs
length(which(points$maxcol==3))/length(points$x)
length(which(nutplot$maxcol==3))/length(nutplot$Id)

#waterdemand
length(which(points$maxcol==4))/length(points$x)
length(which(nutplot$maxcol==4))/length(nutplot$Id)


# length(which(Rcrap$maxcol==3))/length(Rcrap$HYBAS_ID)

#What is happening in the 1980s?
pointagg <- aggregate(list(Rchange_rel = DataL$Y1985),
                      by = list(HydroR = DataL$HER),
                      FUN = function(x) c(mean = mean(x, na.rm = TRUE), 
                                          dev = sd(x, na.rm = TRUE), 
                                          len = length(x), 
                                          med = median(x, na.rm = TRUE), 
                                          q1 = quantile(x, 0.025, na.rm = TRUE), 
                                          q3 = quantile(x, 0.975, na.rm = TRUE)))
pointSoc80s <- do.call(data.frame, pointagg)

trendS130=t(trendClim[which(trendClim$HydroR==130),])[-1]
plot(trendS130)



ggsave(paste0("D:/tilloal/Documents/LFRuns_utils/TrendAnalysis/plots/largestdriver_HR",haz,"23.jpg"), ok, width=20, height=20, units=c("cm"),dpi=1000) 


##6. Map generation -------

# values for plot saving
#method for trend computation (th is thresold or trendPeaks)
mmx="th"
#iteration for plot ID
it=27
###6.1 large loop for change from different drivers ----

###6.1 [Plot] - Figure 2 and Figure 3 ----
yrlist=c(1951:2020)
driverlist=c("climate","landuse","reservoirs","wateruse","all")
driver=driverlist[2]
Dchangelist=list()
for (driver in driverlist){
  print(driver)
  if (driver=="climate"){
    trendPlot=trendClim
    datap=DataC
    pointagg=pointClim
  }
  if (driver=="landuse"){
    trendPlot=trendSoc
    datap=DataL
    pointagg=pointSoc
  }
  if (driver=="reservoirs"){
    trendPlot=trendRes
    datap=DataR
    datap$Y2015[match(Rcrap$outl2,datap$outl2)]=0
    pointagg=pointRes
  }
  if (driver=="wateruse"){
    trendPlot=trendWu
    datap=DataW
    pointagg=pointWu
  }
  if (driver=="all"){
    trendPlot=trendTot
    datap=DataT
    pointagg=pointTot
  }
  Pplot=calculatePoints(trendPlot, yrlist, pointagg, Regio, GHshpp, datap)
  
  save=T
  
  if (save==T){ 
    colNA="transparent"
    
    if (driver=="climate"){
      ####[Plot] - supplement- ordered change aggregated at the HER level ----
      pointP=Pplot$PagD
      uhi=unique(pointP$CODEB)
      pointP=pointP[match(uhi,pointP$CODEB),]
      pointP=pointP[order(pointP$Rchange_rel.mean),]
      pointP$id=c(1:length(pointP$Rchange_rel.mean))
      
      hist(pointP$Rchange_rel.len,breaks=100)
      hist(GHshpp$SURF_KM2,breaks=100)
      mean(GHshpp$SURF_KM2)
      median(pointP$Rchange_rel.len)
      print(length(which(pointP$change==-2)))/length(pointP$Id)
      print(length(which(pointP$change==2)))/length(pointP$Id)
      manualcol=c("-2"="#A51122","-1"="#F1C363", "1"="#ACD2BB","2"= "#324DA0")
      manualab=c("sig. decrease","decrease", "increase","sig. increase")
      brl=c(-200,200)
      by=50
 
      br=c(brl[1],-100,-50,-10,0,10,50,100,brl[2])
      limi=c(-5000,5000)
      
      ggplot() +
        coord_cartesian(ylim=c(brl[1],brl[2]))+
        geom_hline(yintercept = 0,lwd=1, col="black")+
        geom_segment(data=pointP,aes(x=id, xend=id, y=Rchange_rel.q1.2.5.,
                                     yend=Rchange_rel.q3.97.5.,color=factor(change),size = Rchange_rel.len), alpha=.99)+
        scale_color_manual(values = manualcol,breaks=manualab,labels=manualab,name="")+
        geom_point(data=pointP, aes(x=id, y=Rchange_rel.mean), pch=21, fill="white",colour="gray3",size=3,stroke=1,alpha=1) + 
        geom_text(data=pointP, aes(x=id, y=Rchange_rel.mean,label = CODEB), size=1.5, color = "black",fontface = "bold") +
        scale_size(range = c(0.8, 4),trans="sqrt",
                   guide = "none")+
        scale_y_continuous(limits=limi,breaks=br,name="Change (%)",trans=scales::modulus_trans(.3))+
        scale_x_continuous(name="HER",expand=c(.01,.01),breaks=c(-100,200))+
        guides(colour = guide_legend(override.aes = list(size = 10)))+
        theme(axis.title=element_text(size=20, face="bold",color="black"),
              axis.text = element_text(size=18,color="black"),
              panel.background = element_rect(fill = "white", colour = "white"),
              panel.grid = element_blank(),
              panel.border = element_rect(linetype = "solid", fill = NA, colour="black",linewidth=2),
              legend.title = element_text(size=20),
              legend.text = element_text(size=18,color="black"),
              legend.position = "none",
              panel.grid.major = element_line(colour = "transparent"),
              panel.grid.minor.y = element_line(colour = "transparent",linetype="dashed"),
              legend.key = element_rect(fill = "transparent", colour = "transparent"),
              legend.key.size = unit(.8, "cm"))
      
      ggsave(paste0("D:/tilloal/Documents/LFRuns_utils/TrendAnalysis/plots/ordered_HR",it,"_",haz,".jpg"),width=40, height=8, units=c("cm"),dpi=400) 
      
    }
    
    #differentiated plot schemes for flood and drought
    if (haz=="Flood"){
      br=c(-50,-20,-10,-5,0,5,10,20,50)
      labels=br
      limi=c(-50,50)
      tsize=16
      osize=12
      legend2="Change (%)    "
      palet=c(hcl.colors(11, palette = "RdYlBu", alpha = NULL, rev = F, fixup = TRUE))
      paletf=c(hcl.colors(11, palette = "RdBu", alpha = NULL, rev = F, fixup = TRUE))
      points=Pplot$points
      pag=Pplot$PagD
      pointsInside=Pplot$psp
      ####[Plot] - Figure 2b - Map of changes in flood 10-Y RL driven by climatic changes ----
      if (driver=="climate" | driver=="all"){
        
        titleX=paste0("Change in 10-year ",haz," attributed \nto ",driver," changes (% of  10y flood) -  1955-2015")
        fmap<-ggplot(basemap) +
          geom_sf(fill="white",color="darkgrey",size=0.5)+
          geom_sf(data=pag,aes(fill=Rchange_rel.mean,geometry=geometry),alpha=0.2,color="transparent")+
          geom_sf(data=points,aes(col=Y2015,geometry=geometry,size=upa),alpha=.9,stroke=0,shape=15)+ 
          
          scale_size(range = c(0.08, 0.4), trans="sqrt",name= expression(paste("Upstream area ", (km^2),
                                                                               sep = " ")),
                     breaks=c(101,1000,10000,100000,500000), labels=c("100","1000", "10 000", "100 000", "500 000"),
                     guide = "none")+
          scale_fill_gradientn(
            colors=paletf,
            breaks=br,limits=limi,trans=scales::modulus_trans(.3),
            oob = scales::squish,na.value=colNA, name=legend2)   +
          guides(fill = "none")+
          new_scale_fill()+
          geom_sf(data=pointsInside,aes(geometry=geometry, fill=sign),alpha=.6,size=.5,stroke=0,shape=21,color="black")+
          scale_fill_manual(values=c("tomato4","steelblue4"), name="Significant trends")+
          coord_sf(xlim = c(min(nco[,1]),max(nco[,1])), ylim = c(min(nco[,2]),max(nco[,2])))+
          scale_color_gradientn(
            colors=palet,
            breaks=br,limits=limi,trans=scales::modulus_trans(.3),
            oob = scales::squish,na.value=colNA, name=legend2)   +
          labs(x="Longitude", y = "Latitude")+
          guides(colour = guide_colourbar(barwidth = 22, barheight = 1),
                  fill = guide_legend(override.aes = list(size = 10)))+
          theme(axis.title=element_text(size=tsize),
                title = element_text(size=osize),
                axis.text=element_text(size=osize),
                panel.background = element_rect(fill = "aliceblue", colour = "grey1"),
                panel.border = element_rect(linetype = "solid", fill = NA, colour="black"),
                legend.title = element_text(size=tsize),
                legend.text = element_text(size=osize),
                legend.position = "bottom",
                legend.box = "vertical",  # Stack legends vertically
                panel.grid.major = element_line(colour = "grey70"),
                panel.grid.minor = element_line(colour = "grey90"),
                legend.key = element_rect(fill = "transparent", colour = "transparent"),
                legend.key.size = unit(1, "cm"))+
          ggtitle(titleX)
        
        
        #frequency of significant change regions
        ls=length(unique(pointsInside$Id))
        lsp=length(unique(pointsInside$Id[which(pointsInside$sign=="increase")]))
        lsn=length(unique(pointsInside$Id[which(pointsInside$sign=="decrease")]))
        lx=length(pag$Id)
        lsp/lx
        lsn/lx
        ggsave(paste0("D:/tilloal/Documents/LFRuns_utils/TrendAnalysis/plots/mapF_",driver,"_",haz,mmx,"rel_New2",it,".jpg"), fmap, width=22, height=20, units=c("cm"),dpi=1000) 
        
      }else{
          ####[Plot] - Figure 3 - Map of changes in flood 10-Y RL driven by socioeconomic changes ----
          titleX=paste0("Change in 10-year ",haz," attributed \nto ",driver," changes (% of 10y flood) -  1955-2015")
          legend2="Change (%)"
          ocrap<-ggplot(basemap) +
            geom_sf(fill="white",color="darkgrey",size=0.5)+
            geom_sf(data=pag,aes(fill=Rchange_rel.mean,geometry=geometry),alpha=0.6,color="transparent")+
            geom_sf(data=points,aes(col=Y2015,geometry=geometry,size=upa),alpha=.9,stroke=0,shape=15)+ 
            
            scale_size(range = c(0.1, 0.5), trans="sqrt",name= expression(paste("Upstream area ", (km^2),
                                                                                 sep = " ")),
                       breaks=c(101,1000,10000,100000,500000), labels=c("100","1000", "10 000", "100 000", "500 000"),
                       guide = "none")+
            scale_fill_gradientn(
              colors=palet,
              breaks=br,limits=limi,trans=scales::modulus_trans(.3),
              oob = scales::squish,na.value=colNA, name="Change (%)")   +
            coord_sf(xlim = c(min(nco[,1]),max(nco[,1])), ylim = c(min(nco[,2]),max(nco[,2])))+
            scale_color_gradientn(
              colors=palet,
              breaks=br,limits=limi,trans=scales::modulus_trans(.3),
              oob = scales::squish,na.value="transparent", name="Change (%)")   +
            labs(x="Longitude", y = "Latitude")+
            guides(colour = guide_colourbar(barwidth = 1.5, barheight = 14),
                   fill = guide_colourbar(barwidth = 1.5, barheight = 14))+
            theme(axis.title=element_text(size=tsize),
                  title = element_text(size=osize),
                  axis.text=element_text(size=osize),
                  panel.background = element_rect(fill = "aliceblue", colour = "grey1"),
                  panel.border = element_rect(linetype = "solid", fill = NA, colour="black"),
                  legend.text = element_text(size=8),
                  legend.title = element_text(size = osize, margin = margin(t = 2, r = 2, b = 6, l = 0)),
                  legend.spacing.x = unit(0.2, "cm"),
                  legend.position = "bottom",
                  legend.box = "vertical",  # Stack legends vertically
                  panel.grid.major = element_line(colour = "grey70"),
                  panel.grid.minor = element_line(colour = "grey90"),
                  legend.key = element_rect(fill = "transparent", colour = "transparent"),
                  legend.key.size = unit(1, "cm"))+
            ggtitle(titleX)
          
          ggsave(paste0("D:/tilloal/Documents/LFRuns_utils/TrendAnalysis/plots/maF_",driver,"_",haz,"_",it,".jpg"), ocrap, width=23, height=20, units=c("cm"),dpi=1000) 
          
      } 
    }else if(haz=="Drought"){
        br=c(-50,-20,-10,-5,0,5,10,20,50)
        labels=br
        limi=c(-50,50)
        tsize=16
        osize=12
        legend2="Change (%)"
        palet=c(hcl.colors(11, palette = "RdYlBu", alpha = NULL, rev = F, fixup = TRUE))
        paletf=c(hcl.colors(11, palette = "RdYlBu", alpha = NULL, rev = F, fixup = TRUE))
        points=Pplot$points
        pag=Pplot$PagD
        pointsInside=Pplot$psp
          ####[Plot] - Figure 2a - Map of changes in drought 10-Y RL driven by climatic changes ----
        if (driver=="climate" | driver=="all"){
          #specifically designed plot for climate
          titleX=paste0("Change in 10-year ",haz," attributed \n to ",driver, "changes (% of 10y drought) - 1955-2015")
          points=points[-which(is.na(points$Y2015)),]
          ocrap<-ggplot(basemap) +
            geom_sf(fill="white",color="darkgrey",size=0.5)+
            geom_sf(data=pag,aes(fill=Rchange_rel.mean,geometry=geometry),alpha=0.2,color="transparent")+
            geom_sf(data=points,aes(col=Y2015,geometry=geometry,size=upa),alpha=.9,stroke=0,shape=15)+ 
            
            scale_size(range = c(0.08, 0.4), trans="sqrt",name= expression(paste("Upstream area ", (km^2),
                                                                                 sep = " ")),
                       breaks=c(101,1000,10000,100000,500000), labels=c("100","1000", "10 000", "100 000", "500 000"),
                       guide = "none")+
            scale_fill_gradientn(
              colors=paletf,
              breaks=br,limits=limi,labels = labels, trans=scales::modulus_trans(.3),
              oob = scales::squish,na.value=colNA, name=legend2)   +
            guides(fill = "none")+
            new_scale_fill()+
            geom_sf(data=pointsInside,aes(geometry=geometry, fill=sign),alpha=.6,size=.5,stroke=0,shape=21,color="black")+
            scale_fill_manual(values=c("tomato4","steelblue4"), name="Significant trends")+
            coord_sf(xlim = c(min(nco[,1]),max(nco[,1])), ylim = c(min(nco[,2]),max(nco[,2])))+
            scale_color_gradientn(
              colors=palet,
              breaks=br,limits=limi,labels = labels,trans=scales::modulus_trans(.3),
              oob = scales::squish,na.value="transparent", name=legend2)   +
            labs(x="Longitude", y = "Latitude")+
            guides(colour = guide_colourbar(barwidth = 1.5, barheight = 14),
                   fill = guide_legend(override.aes = list(size = 10)))+
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
                  legend.key.size = unit(1, "cm"))+
            ggtitle(titleX)
          
          ggsave(paste0("D:/tilloal/Documents/LFRuns_utils/TrendAnalysis/plots/mapF_",driver,"_",haz,"_",it,".jpg"), ocrap, width=22, height=20, units=c("cm"),dpi=1000) 
          
        }else{
          ####[Plot] - Figure 3 - Map of changes in drought 10-Y RL driven by socioeconomic changes ----
          titleX=paste0("Change in 10-year ",haz," attributed \nto ", driver," changes (% of  10y drought) - 1955-2015")
          ocrap<-ggplot(basemap) +
            geom_sf(fill="white",color="darkgrey",size=0.5)+
            geom_sf(data=pag,aes(fill=Rchange_rel.mean,geometry=geometry),alpha=0.6,color="transparent")+
            geom_sf(data=points,aes(col=Y2015,geometry=geometry,size=upa),alpha=.9,stroke=0,shape=15)+ 
            
            scale_size(range = c(0.1, 0.5), trans="sqrt",name= expression(paste("Upstream area ", (km^2),
                                                                                sep = " ")),
                       breaks=c(101,1000,10000,100000,500000), labels=c("100","1000", "10 000", "100 000", "500 000"),
                       guide = "none")+
            scale_fill_gradientn(
              colors=paletf,
              breaks=br,limits=limi,labels = labels,trans=scales::modulus_trans(.3),
              oob = scales::squish,na.value=colNA, name="Change (%)")   +
            coord_sf(xlim = c(min(nco[,1]),max(nco[,1])), ylim = c(min(nco[,2]),max(nco[,2])))+
            scale_color_gradientn(
              colors=paletf,
              breaks=br,limits=limi,labels = labels, trans=scales::modulus_trans(.3),
              oob = scales::squish,na.value=colNA, name="Change (%)")   +
            labs(x="Longitude", y = "Latitude")+
            guides(colour = guide_colourbar(barwidth = 1.5, barheight = 14),
                   fill = guide_colourbar(barwidth = 1.5, barheight = 14))+
            theme(axis.title=element_text(size=tsize),
                  title = element_text(size=osize),
                  axis.text=element_text(size=osize),
                  panel.background = element_rect(fill = "aliceblue", colour = "grey1"),
                  panel.border = element_rect(linetype = "solid", fill = NA, colour="black"),
                  legend.text = element_text(size=8),
                  legend.title = element_text(size = osize, margin = margin(t = 2, r = 2, b = 6, l = 0)),
                  legend.spacing.x = unit(0.2, "cm"),
                  legend.position = "right",
                  legend.box = "horizontal",  # Stack legends vertically
                  panel.grid.major = element_line(colour = "grey70"),
                  panel.grid.minor = element_line(colour = "grey90"),
                  legend.key = element_rect(fill = "transparent", colour = "transparent"),
                  legend.key.size = unit(1, "cm"))+
            ggtitle(titleX)
          
          ggsave(paste0("D:/tilloal/Documents/LFRuns_utils/TrendAnalysis/plots/mapF_",driver,"_",haz,"_",it,".jpg"), ocrap, width=23, height=20, units=c("cm"),dpi=1000) 
          
        }
      }
  }
  Dchangelist=c(Dchangelist,list(Pplot))
}


#create a barplot with 1 category per driver and number of locations with increasing|decreasing trend
lm=match(pointsag,RegioRLi$HydroR)
RegioRLi$HydroR[lm]

paggC=Dchangelist[[1]]$PagD
paggL=Dchangelist[[2]]$PagD
paggR=Dchangelist[[3]]$PagD
paggW=Dchangelist[[4]]$PagD
paggC$driver="Clim"
paggL$driver="Lu"
paggR$driver="Res"
paggW$driver="Wu"
colnames(paggC)
colnames(paggR)
pointsAD=rbind(paggC,paggL,paggR,paggW)

length(which(DataC$Y2020>15))/length(DataC$Y2020)
length(which(pointSoc$Rchange_rel.mean<0))/length(pointSoc$Rchange_rel.mean)



TdataClim=processTrendData(trendData =  trendClim, DataTr = DataC,
                           id_var = "HydroR")
TdataLuse=processTrendData(trendSoc,DataL,id_var = "HydroR")
TdataRes=processTrendData(trendRes,DataR,id_var = "HydroR")
TdataWuse=processTrendData(trendWu,DataW,id_var = "HydroR")

bio_names=unique(biogeo$code)

#I keep only some bioregions
bio_names=bio_names[c(1,3,4,6,7,9,11)]

# combine trends from all drivers
trtF=rbind(TdataClim$tGlobal,TdataLuse$tGlobal,TdataRes$tGlobal,TdataWuse$tGlobal)
trtF$year=rep(seq(1950,2010,10),4)
trtF$driver=c(rep("Clim",7),rep("LUC",7),rep("Res",7),rep("WU",7))

names(trtF)[c(2,4,5,6,7,8)]=c("changeC","med","cq1","cq2","w1","w2")

colorz = c("Clim" ='dodgerblue4',"LUC" ='gold4',"Res" ='firebrick4',"WU"="olivedrab")
colorn = c("WU" ='limegreen',"Res" ='tomato4',"LUC" ='orange',"Clim" ='royalblue')


#### [Plot] - Figure 1 - Boxplot in time ----

fac=1
xlabs=seq(1950,2010,10)
clabels=c("Climate","Land use","Reservoirs", "Water demand")
nplot="Mean change (% 10yRL)"
br=c(seq(-100,-10,10),seq(-5,5,5),seq(10,100,10))
ggplot() +
  # IQR represented as rectangles
  geom_linerange(data=trtF,aes(x=year, ymin=fac*(w1),ymax=fac*w2,color = factor(driver),group=factor(driver)),
                 position = position_dodge2(width = 9),lwd=1,alpha=0.8) +
  scale_color_manual(values = colorn, name = "Drivers", labels = clabels) +
  new_scale_color()+
  
  geom_rect(data=trtF, aes(xmin = year - 4.5, xmax = year + 4.5, 
                           ymin = cq1, ymax = cq2, fill = factor(driver), group = factor(driver)), 
            alpha = 0.5, position = position_dodge(width = 9)) +
  
  # Mean as points over the IQR
  geom_point(data=trtF, aes(x = year, y = changeC, color = factor(driver), group = factor(driver)),
             position = position_dodge(width = 9), size = 3) +
  
  # Median as a horizontal line within the IQR rectangle
  geom_rect(data=trtF, aes(xmin = year - 4.5, xmax = year + 4.5, 
                           ymin = med-1e-1, ymax = med+1e-1, fill = factor(driver), group = factor(driver)), 
            alpha = 1, position = position_dodge(width = 9)) +
  
  # Y-axis settings
  scale_y_continuous(name = nplot, breaks = br, trans=scales::modulus_trans(.6)) +
  
  # X-axis settings
  scale_x_continuous(breaks = xlabs, labels = xlabs, name = "Decades",
                     minor_breaks = seq(1955,2005,10), expand = c(.01,0.01)) +
  
  # Manual fill and color scales
  scale_fill_manual(values = colorn, name = "Drivers", labels = clabels) +
  scale_color_manual(values = colorz, name = "Drivers", labels = clabels) +
  
  # Customize legend
  guides(color = guide_legend(override.aes = list(color = colorz))) +
  
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
    #panel.grid.minor.y = element_line(color = "lightgray"),
    legend.position = "right",
    #panel.grid.major = element_line(colour = "grey80"),
    panel.grid.minor.x = element_line(colour = "grey23",linetype = "dashed"),
    legend.key = element_rect(fill = "transparent", colour = "transparent"),
    legend.key.size = unit(0.8, "cm")
  ) +
  
  # Add title
  ggtitle("Europe")

ggsave(paste0("D:/tilloal/Documents/LFRuns_utils/TrendAnalysis/plots/bxp_in_time_EU7_",haz,"_23.jpg"), width=30, height=20, units=c("cm"),dpi=1000) 

#### [Plot] - Supplements - Boxplot in time by BGR ----
for (bn in bio_names){
  print(bn)
  trtF=rbind(TdataClim$BgData[which(TdataClim$BgData$loc==bn),],TdataLuse$BgData[which(TdataLuse$BgData$loc==bn),],
             TdataRes$BgData[which(TdataRes$BgData$loc==bn),],TdataWuse$BgData[which(TdataWuse$BgData$loc==bn),])
  
  
  trtF$year=rep(seq(1950,2010,10),4)
  trtF$driver=c(rep("Clim",7),rep("LUC",7),rep("Res",7),rep("WU",7))
  
  names(trtF)[c(3,5,6,7,8,9)]=c("changeC","med","cq1","cq2","w1","w2")
  
  colorz = c("Clim" ='dodgerblue4',"LUC" ='gold4',"Res" ='firebrick4',"WU"="olivedrab")
  colorn = c("WU" ='limegreen',"Res" ='tomato4',"LUC" ='orange',"Clim" ='royalblue')
  
  fac=1
  xlabs=seq(1950,2010,10)
  clabels=c("Climate","Land use","Reservoirs", "Water demand")
  
  nplot="Mean change (% of  10Y RL)"
  br=c(seq(-200,-10,10),seq(-5,5,5),seq(10,200,10))
  
  ggplot() +
    # IQR represented as rectangles
    geom_linerange(data=trtF,aes(x=year, ymin=fac*(w1),ymax=fac*w2,color = factor(driver),group=factor(driver)),
                   position = position_dodge2(width = 9),lwd=1,alpha=0.8) +
    scale_color_manual(values = colorn, name = "Drivers", labels = clabels) +
    new_scale_color()+
    
    geom_rect(data=trtF, aes(xmin = year - 4.5, xmax = year + 4.5, 
                             ymin = cq1, ymax = cq2, fill = factor(driver), group = factor(driver)), 
              alpha = 0.5, position = position_dodge(width = 9)) +
    
    # Mean as points over the IQR
    geom_point(data=trtF, aes(x = year, y = changeC, color = factor(driver), group = factor(driver)),
               position = position_dodge(width = 9), size = 3) +
    
    # Median as a horizontal line within the IQR rectangle
    geom_rect(data=trtF, aes(xmin = year - 4.5, xmax = year + 4.5, 
                             ymin = med-1e-1, ymax = med+1e-1, fill = factor(driver), group = factor(driver)), 
              alpha = 1, position = position_dodge(width = 9)) +
    
    # Y-axis settings
    scale_y_continuous(name = nplot, breaks = br,trans=scales::modulus_trans(.6)) +
    
    # X-axis settings
    scale_x_continuous(breaks = xlabs, labels = xlabs, name = "Decades",
                       minor_breaks = seq(1955,2005,10), expand = c(.01,0.01)) +
    
    # Manual fill and color scales
    scale_fill_manual(values = colorn, name = "Drivers", labels = clabels) +
    scale_color_manual(values = colorz, name = "Drivers", labels = clabels) +
    
    # Customize legend
    guides(color = guide_legend(override.aes = list(color = colorz))) +
    
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
     # panel.grid.minor.y = element_line(color = "lightgray"),
      legend.position = "right",
      #panel.grid.major = element_line(colour = "grey80"),
      panel.grid.minor.x = element_line(colour = "grey23", linetype = "dashed"),
      legend.key = element_rect(fill = "transparent", colour = "transparent"),
      legend.key.size = unit(0.8, "cm")
    ) +
    
    # Add title
    ggtitle(bn)
  
   
  ggsave(paste0("D:/tilloal/Documents/LFRuns_utils/TrendAnalysis/plots/bxp_in_time",bn,"_",haz,"_23.jpg"), width=30, height=20, units=c("cm"),dpi=1000) 
  
}

#[Plot] - Supplement - box plot by catchment size
Cdata=DataC
Cdata$upagroup=1
Cdata$upagroup[which(Cdata$upa>200 & Cdata$upa<=500)]=2
Cdata$upagroup[which(Cdata$upa>500 & Cdata$upa<=1000)]=3
Cdata$upagroup[which(Cdata$upa>1000 & Cdata$upa<=10000)]=4
Cdata$upagroup[which(Cdata$upa>10000)]=5

ClimUpag<-UpATrendData (DataC)
LuseUpag<-UpATrendData (DataL)
ResUpag<-UpATrendData (DataR)
WdemUpag<-UpATrendData (DataW)

UpaF=rbind(ClimUpag,LuseUpag,ResUpag,WdemUpag)
UpaF$driver=c(rep("Clim",5),rep("LUC",5),rep("Res",5),rep("WU",5))
UpaF=UpaF[,-1]
names(UpaF)[c(2,4,5,6,7,8)]=c("changeC","med","cq1","cq2","w1","w2")

colorz = c("Clim" ='dodgerblue4',"LUC" ='gold4',"Res" ='firebrick4',"WU"="olivedrab")
colorn = c("WU" ='limegreen',"Res" ='tomato4',"LUC" ='orange',"Clim" ='royalblue')

#use pointap for bars
fac=1
xlabs=seq(1950,2010,10)
clabels=c("Climate","Land use","Reservoirs", "Water demand")

unique(UpaF$loc)
nplot="Mean change (% of mean 10yRL)"
br=c(seq(-300,-20,20),seq(-20,20,10),seq(20,300,20))
ggplot() +
  # IQR represented as rectangles
  geom_linerange(data=UpaF,aes(x=loc, ymin=fac*(w1),ymax=fac*w2,color = factor(driver),group=factor(driver)),
                 position = position_dodge2(width = .9),lwd=1,alpha=0.8) +
  scale_color_manual(values = colorn, name = "Drivers", labels = clabels) +
  new_scale_color()+
  
  geom_rect(data=UpaF, aes(xmin = loc - 0.45, xmax = loc + 0.45, 
                           ymin = cq1, ymax = cq2, fill = factor(driver), group = factor(driver)), 
            alpha = 0.5, position = position_dodge(width = .9)) +
  
  # Median as points over the IQR
  geom_point(data=UpaF, aes(x = loc, y = changeC, color = factor(driver), group = factor(driver)),
             position = position_dodge(width = .9), size = 3) +
  
  # Median as a horizontal line within the IQR rectangle
  geom_rect(data=UpaF, aes(xmin = loc - .45, xmax = loc + .45, 
                           ymin = med-1e-1, ymax = med+1e-1, fill = factor(driver), group = factor(driver)), 
            alpha = 1, position = position_dodge(width = .9)) +
  
  # Y-axis settings
  scale_y_continuous(name = nplot, breaks = br,trans=scales::modulus_trans(.6)) +
  
  # X-axis settings
  scale_x_continuous(breaks=c(1,2,3,4,5),labels=c( "100-200",  "200-500",
                             "500-1000", "1000-10 000",
                              ">10 000"),name="Catchment Area (km2)")+
  # Manual fill and color scales
  scale_fill_manual(values = colorn, name = "Drivers", labels = clabels) +
  scale_color_manual(values = colorz, name = "Drivers", labels = clabels) +
  
  # Customize legend
  guides(color = guide_legend(override.aes = list(color = colorz))) +
  
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
    panel.grid.major.y = element_line(color = "lightgray"),
    panel.grid.minor.y = element_line(color = "lightgray"),
    legend.position = "right",
    panel.grid.major = element_line(colour = "grey80"),
    panel.grid.minor.x = element_line(colour = "grey90", linetype = "dashed"),
    legend.key = element_rect(fill = "transparent", colour = "transparent"),
    legend.key.size = unit(0.8, "cm")
  ) +
  
  # Add title
  ggtitle("Europe")

ggsave(paste0("D:/tilloal/Documents/LFRuns_utils/TrendAnalysis/plots/bxp_in_UpA_EU_",haz,"_23.jpg"), width=30, height=20, units=c("cm"),dpi=1000) 



##7. SAVING OUTPUTS ###################

#save important outputs for next script
DataT$driver="Total"
DataC$driver="Clim"
DataR$driver="Reservoirs"

DataL$driver="Landuse"
DataW$driver="Wateruse"


Alltrend=rbind(DataC,DataL,DataR,DataW,DataT)

trendTot$driver="Total"
trendClim$driver="Clim"
trendRes$driver="Reservoirs"
trendSoc$driver="Landuse"
trendWu$driver="Wateruse"

trendRegio=rbind(trendClim,trendSoc,trendRes,trendWu,trendTot)

if(haz=="Flood"){
  Output_fl_year=list(TrendPix=Alltrend,TrendRegio=trendRegio,Out2020=pointsAD,DataI=DataI)
  save(Output_fl_year,file=paste0(hydroDir,"/TSEVA/output_plots/outputs_flood_year_relxHR_8.Rdata"))
}
if (haz=="Drought")
{
  Output_dr_nonfrost=list(TrendPix=Alltrend,TrendRegio=trendRegio,Out2020=pointsAD,DataI=DataI)
  save(Output_dr_nonfrost,file=paste0(hydroDir,"/TSEVA/output_plots/outputs_drought_nonfrost_relxHR_9.Rdata"))
}
