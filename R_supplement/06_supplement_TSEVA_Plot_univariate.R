
#plot shape parameter


max(Shapepar$epsilonGPD)
hist(Shapepar$epsilonGPD,breaks=100)
mean(Shapepar$epsilonGPD)

#now plot the shape parameter
Shapeplot=inner_join(Shapepar,UpArea,by=c("catchment"= "outl2"))
points <- st_as_sf(Shapeplot, coords = c("Var1.x", "Var2.x"), crs = 4326)
points <- st_transform(points, crs = 3035)

pls=ggplot(basemap) +
  geom_sf(fill="gray85")+
  geom_sf(fill=NA, color="grey") +
  geom_sf(data=points,aes(col=epsilonGPD,geometry=geometry),alpha=1,size=0.25,stroke=0,shape=15)+ 
  coord_sf(xlim = c(min(nco[,1]),max(nco[,1])), ylim = c(min(nco[,2]),max(nco[,2])))+
  scale_color_gradientn(
    colors=palet,
    breaks=c(-0.5,-0.25,0,0.25,0.5,0.75,1),limits=c(-0.5,0.5),
    oob = scales::squish,na.value=colNA, name=legend)   +
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

pls




#identify dominant driver for each region


ziz=full_join(GHshpp,pointagg,by=c("CODEB"="HydroR"))
st_geometry(ziz)<-NULL
#natch pointagg with hybasf
pointplot=inner_join(HydroRsf,ziz,by= c("Id"="Id"))

# legend="10y flood specific discharge (l/s/km2)"
# title=paste0("10-years ",haz," Return Level")
# br=c(1,5,10,25,100,250,1000)
# labels=br
# limi=c(1,1200)

haz="flood"
palet=c(hcl.colors(11, palette = "RdYlBu", alpha = NULL, rev = F, fixup = TRUE))
title=paste0("Change in 10-years ",haz," Return Level between ", period[1], " and ", period[2])
legend="Change in specific discharge (l/s/km2)"
if (haz=="drought"){
  data[,valcol]=tmpval*100
  br=c(-150,-100,-50,0,50,100,150)
  labels=br/100
  limi=c(-150,150)
}
if (haz=="flood"){
  # data[,valcol]=tmpval
  br=c(-20,-15,-10,5,0,5,10,15,20)
  labels=br
  limi=c(-20,20)
}
pl2=ggplot(basemap) +
  geom_sf(fill="white")+
  geom_sf(data=pointplot,aes(fill=Rchange.med,geometry=geometry),color="transparent",alpha=.5)+ 
  geom_sf(data=points,aes(col=Y2020,geometry=geometry),alpha=1,size=0.25,stroke=0,shape=15)+ 
  #geom_sf(fill=NA, color="grey") +
  geom_sf(data=pointplot,aes(geometry=geometry),color="darkgrey",fill=NA)+ 
  scale_fill_gradientn(
    colors=palet,
    breaks=br,labels=labels, limits=limi,
    oob = scales::squish,na.value=colNA, name=legend)   +
  scale_color_gradientn(
    colors=palet,
    breaks=br,limits=limi,
    oob = scales::squish,na.value=colNA, name=legend)   +
  #geom_sf(data=dt3,aes(geometry=geometry,fill=factor(plot)),color="gray12")+ 
  #scale_fill_manual(values=c("1"="purple","-1"="tan2"), name="IRES trend")+
  coord_sf(xlim = c(min(nco[,1]),max(nco[,1])), ylim = c(min(nco[,2]),max(nco[,2])))+
  labs(x="Longitude", y = "Latitude")+
  guides(color = guide_colourbar(barwidth = 20, barheight = .8), fill="none")+
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
        legend.key.size = unit(.8, "cm"))+
  ggtitle(title)

min(points$fillplot)
ggsave("D:/tilloal/Documents/LFRuns_utils/TrendAnalysis/plots/mapHR_abs2.jpg", pl2, width=20, height=20, units=c("cm"),dpi=1000) 



pointagg=aggregate(list(Rchange=dataX[,valcol+8]),
                   by = list(HydroR=dataX$HR_id),
                   FUN = function(x) c(mean=mean(x,na.rm=T)))
pointagg <- do.call(data.frame, pointagg)


pointagp=aggregate(list(Rchange=dataX[,valcol+8]),
                   by = list(HydroR=dataX$HR_id),
                   FUN = function(x) c(mean=mean(x,na.rm=T),sd=sd(x,na.rm=T),q10=quantile(x,0.1, na.rm=T),q90=quantile(x,0.9,na.rm=T)))
pointagp <- do.call(data.frame, pointagp)


pap=data.frame(t(pointagp))

tpag=data.frame(t(pointagg))

colnames(tpag)=tpag[1,]
tpag=tpag[-1,]


# Compute the correlation matrix
cor_mat <- cor(tpag[ -c(1,2),])

# Perform hierarchical clustering
hclust_res <- hclust(as.dist(1 - cor_mat))
# Plot the dendrogram
plot(hclust_res, main = "Hierarchical Clustering Dendrogram")


# Extract the cluster assignments


print(dunn_index)
wss <- vector()
dunn_index <- vector()
for (i in 1:10) {
  cluster_assignments <- cutree(hclust_res, k = i)
  wss[i] <- sum(cor_mat[cluster_assignments == cluster_assignments[1]]^2)
  dunn_index[i] <- dunn(as.dist(1 - cor_mat),cluster_assignments)
}

# Plot the WSS as a function of the number of clusters
plot(1:10, wss, type = "b", xlab = "Number of Clusters", ylab = "WSS")

plot(2:10, dunn_index[2:10], type = "b", xlab = "Number of Clusters", ylab = "WSS")

cluster_assignments <- cutree(hclust_res, k = 4)

#legend="Kendall tau"
legend="Relative change (%)"
tmpval=(datatwin[,valcol2])
mkta=c()
mksa=c()
for (it in 1:length(pointagg[,1])){
  print(it)
  mks=NA
  mkt=NA
  miniTS=as.numeric(pointagg[it,])
  if (!is.na(miniTS[2])){
    #mk=MannKendall(miniTS[-1])
    mk2=mmkh(miniTS[-1],ci=0.95)
    mk=data.frame(tau=mk2[6],sl=mk2[2])
    
    #compute trend as well with sen.slope
    mkt=mk$tau
    mks=mk$sl
    print(mkt)
  }else{
    mkt=NA
    mks=NA
  }
  mkta=c(mkta,mkt)
  mksa=c(mksa,mks)
}
# br=c(-1,-.80,-.60,-.40,-.20,0,.20,.40,.60,.80,1)
# limi=c(-1,1)
trans=scales::modulus_trans(.6)
colNA="darkgrey"
#title="Trends in 100 years Re flood magnitude"


#plot all the changes through all catchments
pa_l <- reshape2::melt(pointagg, id.vars = "HydroR", variable.name = "Year", value.name = "value")
pa_l$value=as.numeric(pa_l$value)

tpmean=aggregate(list(rl=pa_l$value),
                 by = list(HydroR=pa_l$HydroR),
                 FUN = function(x) c(mean=mean(x,na.rm=T)))
tpmean <- do.call(data.frame, tpmean)

pa_l$value2=pa_l$value/tp
# Add the cluster assignments to the data frame
pa_l$cluster <- factor(cluster_assignments)

pointagg$mkta=mkta
pointagg$sl=mksa
meds <- c(by(pa_l$value, pa_l$HydroR, mean))

plot(meds)
pa2=pa_l[which(pa_l$cluster==4),]
ggplot(pa2, aes(x = Year, y = value, group = HydroR,col=cluster)) +
  geom_line() +
  scale_y_continuous(limits=c(-10,10))+
  theme_minimal()

plot(as.numeric(pointagg[1,-1]))

pointagg$cluster=factor(cluster_assignments)

#natch pointagg with hybasf
pointsag=inner_join(GHshpp,pointagg,by=c("CODEB"="HydroR"))



br=c(-30,-20,-10,0,10,20,30)
labels=br
limi=c(-30,30)
tsize=16
osize=12
legend="Change in Qsp \n(l/s/km2)"

ggplot(basemap) +
  geom_sf(fill="gray85",color="darkgrey",size=0.5)+
  geom_sf(data=pointsag,aes(fill=cluster,geometry=geometry),alpha=0.8,color="transparent")+
  #geom_sf(data=pointsag,aes(color=mkta,geometry=geometry),alpha=0.4,fill="transparent")+
  scale_color_gradientn(
    colors=palet,
    breaks=seq(-1,1,by=0.2), limits=c(-1,1),
    oob = scales::squish,na.value=colNA, name=legend) +
  coord_sf(xlim = c(min(nco[,1]),max(nco[,1])), ylim = c(min(nco[,2]),max(nco[,2])))+
  # geom_sf(data=datasl1,aes(col=fillplot,geometry=geometry),alpha=0.8,size=0.05,stroke=0,shape=1)+ 
  # geom_sf(data=datasl2,aes(col=fillplot,geometry=geometry),alpha=0.9,size=0.1,stroke=0,shape=1)+ 
  #geom_sf(data=data,aes(fill=fillplot,geometry=geometry),color="transparent")+ 
  #geom_sf(data=datasig_f,aes(size=siglvl),fill="grey",color="transparent",shape=21,alpha=0.1)+ 
  #scale_color_manual(values=c("blue" ="darkblue","red"="darkred"))+
  #geom_col_pattern(data=pointsag,aes(fill=siglvl,geometry=geometry),colour='black', pattern = 'circle') +
  #scale_fill_manual(values=c("tomato4","steelblue4"), name="Significant trends")+
  ggtitle(title)+
  # coord_sf(xlim = c(min(nco[,1]),max(nco[,1])), ylim = c(min(nco[,2]),max(nco[,2])))+
  # geom_point(data=event,aes(x=lon,y=lat,col=dr))+
  # scale_color_distiller(palette = "Spectral",
  #                       direction = -1, limit=c(0,max(eventp$dr)),oob = scales::squish)+
  # scale_color_gradientn(
  #   colors=palet,
  #   breaks=br,limits=limi,trans=trans,
  #   oob = scales::squish,na.value=colNA, name=legend)   +
  labs(x="Longitude", y = "Latitude")+
  guides(colour = guide_colourbar(barwidth = 1.5, barheight = 10),
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
        legend.key.size = unit(1, "cm"))




pointsag$siglvl=0
pointsag$siglvl[which(pointsag$sl<=0.1)]=1

catsig=pointsag[which(pointsag$siglvl>0),]
catsig <- st_transform(catsig, crs = 3035)
## create a grid of points
grdpts <- sf::st_make_grid(catsig, what = "centers",cellsize = 40000)
my.points <- sf::st_sf(grdpts)

sf_use_s2(FALSE)
pointsInside <- sf::st_join(x = my.points, y = catsig, left = FALSE)
pointsInside$sign="positive"
pointsInside$sign[which(pointsInside$mkta<0)]="negative"
pointsag$siglvl=factor(pointsag$siglvl)


#create a time vector with years
dates <- seq.Date(ymd("1951-06-01"), ymd("2020-06-01"), by = "year")
dates= as.POSIXct(dates, format = "%Y-%m-%d")
dfall=c()
plotOut=FALSE
tail="high"


# Plotting the outputs ----------------------------------------------------

## 1. selection of data to be plotted ----

### drought ----

datar1=RLGPDdr
#problem with 2020, to be solved
datar1$Y2020=datar1$Y2019
length(which(!is.na(datar1$Y2020)))
hist(datar1$Y2020,xlim=c(0,1500), breaks=1000)
#Inverting return period values
datar1[,c(1:71)]=-datar1[,c(1:71)]
datar=data.table(datar1)



### flood ----
RLGPDflSCF=RLGPDfl


datarSCF=RLGPDflSCF
datarSCF$Y2020=datarSCF$Y2019
#add latitude and longitude to input



length((datar$unikout))




datariSCF=inner_join(outf,datarSCF,by = c("outl2"="unikout"))
#datar$Y2020=datar$Y2019
unikout=datar$unikout


# datar=RLGPDfl
# datariH=inner_join(outf,datar,by = c("outlets"="unikout"))
# datar$Y2020=datar$Y2019
# unikout=datar$unikout
# 
# 
# #add latitude and longitude to input
# datari=inner_join(outf,datar,by = c("outlets"="unikout"))
# 
# datariSCF=inner_join(outf,datar,by = c("outlets"="unikout"))
# #join with catchment data
dataricat=inner_join(Catamere07,datariSCF,by = c("llcoord"="latlong"))

#checkpars=inner_join(datari,Paramsfl[which(Paramsfl$Year==2020),],by=c("outlets"="catchment"))

Impdates=seq(1950,2020,by=10)
valuenames=paste0("Y",Impdates)

## 2. different plots of the results ----

### Historical distribution of changes ----
totalch=plotHistoDates(outf,datariSCF,law="GPD",type,period=c(1951,2020),parlist=Paramsfl,valuenames)
totalch[[2]]

#trick to also correct
# v1=Paramsfl[which(Paramsfl$Year==2019),]
# v1=v1[match(unique(v1$catchment),v1$catchment),]
# v1$Year=2020
# Paramsfl[which(Paramsfl$Year==2020),]=v1

### Changes in RLs or RPs at pixel and catchment levels ----
Plot.change=plotchangemapix(basemap,catmap=cst7,datariSCF, law="GPD",type="RLchange",period=c(1950,2020),Paramsfl,hybasf = HydroRsf,haz="flood")

Plot.change[[1]]
ggsave("D:/tilloal/Documents/LFRuns_utils/Figures/RP_change_floodCal3.jpg", width=20, height=15, units=c("cm"),dpi=1500)

Plot.change[[2]]
ggsave("D:/tilloal/Documents/LFRuns_utils/Figures/RL_change_hybas07.jpg", width=20, height=15, units=c("cm"),dpi=1000)

### Changes in RLs with mk test at catchment level ----
Plot.change.sig=plotTrendSipix(basemap,dataricat,period=c(1951,2020),hybasf = hybasf7,valuenames,nco)
Plot.change.sig
ggsave("D:/tilloal/Documents/LFRuns_utils/Figures/RLchange_flood_sig19502020.jpg",Plot.change.sig, width=20, height=15, units=c("cm"),dpi=1500)



#load upstream area
main_path = 'D:/tilloal/Documents/06_Floodrivers/'
valid_path = paste0(main_path,'DataPaper/')
outletname="/GIS/upArea_European_01min.nc"
dir=valid_path
outf$idlalo=paste(outf$idlo, outf$idla, sep=" ")
UpArea=UpAopen(valid_path,outletname,outf)
head(UpArea)

Impdates=seq(1951,2020,by=1)
valuenames=paste0("Y",Impdates)
period=c(1951,2020)
haz="drought"




# TO be improved
plotchangemapix_qspU=function(basemap,catmap,datar,upArea, GHR_riv, HydroRsf, law="GPD",type,period=c(1951,2020),parlist,valuenames,haz){
  
  datar=datariSCF
  data=datar
  # data=right_join(GHR_riv,datar,by = c("outlets"="unikout"))
  # data=right_join(cst7,data,by = c("outlets"="outlets"))
  #names(data)[28]="HydroR"
  upag=match(data$outlets,UpArea$outlets)
  data$uparea=UpArea$upa[upag]
  data$outlets
  data$pointid
  datacol=names(data)
  valcol=match(valuenames,datacol)
  # data=data[which(data$IRES==1),]
  datatwin=data
  st_geometry(datatwin) <- NULL
  valcol2=valcol
  
  datatwin=as.data.frame(datatwin)
  class(datatwin)
  dtc=names(datatwin)
  # mcor=unique(match(datar$unikout,Paramsdr$catchment))
  # 
  # if (length(which(is.na(mcor)))>0) datar=datar[-which(is.na(mcor)),]
  # 
  cref=paste0("Y",period[1])
  crefloc=match(cref,dtc)
  finalperiod=paste0("Y",period[2])
  colsel=match(finalperiod,datacol)
  
  
  palet=c(hcl.colors(11, palette = "RdYlBu", alpha = NULL, rev = F, fixup = TRUE))
  title=paste0("Change in 10-years ",haz," Return Level between ", period[1], " and ", period[2])
  legend="Change in specific discharge (l/s/km2)"
  # datathin=datatwin[,c(valcol2)]
  #value_L_s_m2 <- (depth_mm_day / 1000) / 86400 * 1000
  tmpval=(datatwin[,valcol2]-datatwin[,crefloc])
  
  #not change but raw values
  #tmpval=(datatwin[,valcol2])
  tmpval=tmpval*1000/(data$uparea)
  if (haz=="drought"){
    data[,valcol]=tmpval*100
    br=c(-150,-100,-50,0,50,100,150)
    labels=br/100
    limi=c(-150,150)
  }
  if (haz=="flood"){
    data[,valcol]=tmpval
    br=c(-30,-20,-10,0,10,20,30)
    labels=br
    limi=c(-30,30)
  }
  trans=scales::modulus_trans(.8)
  colNA="gray10"
  
  names(data)[colsel]="fillplot"
  data[which(data$IRES==1),colsel]=NA
  
  
  
  data2=inner_join(data,catmap[,c(1:20)],by=c("latlong"="llcoord"))
  points <- st_as_sf(data2[,c(1:80)], coords = c("Var1", "Var2"), crs = 4326)
  points <- st_transform(points, crs = 3035)
  
  pl1=ggplot(basemap) +
    geom_sf(fill="gray85")+
    geom_sf(fill=NA, color="grey") +
    geom_sf(data=points,aes(col=fillplot,geometry=geometry),alpha=1,size=0.25,stroke=0,shape=15)+ 
    coord_sf(xlim = c(min(nco[,1]),max(nco[,1])), ylim = c(min(nco[,2]),max(nco[,2])))+
    scale_color_gradientn(
      colors=palet,
      breaks=br,limits=limi,trans=trans,
      oob = scales::squish,na.value=colNA, name=legend)   +
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
          legend.key.size = unit(.8, "cm"))+
    ggtitle(title)
  
  pl1
  
  ggsave("D:/tilloal/Documents/LFRuns_utils/TrendAnalysis/plots/map_verifx.jpg", pl1, width=20, height=20, units=c("cm"),dpi=1000) 
  # savecrap=inner_join(data,Catamere07,by= c("HYBAS_ID"))
  # st_write(data, paste0(hydroDir,"/disasterdrought3.shp"))
  # 
  
  #  points[68270,]
  #  colnames(Paramsfl)
  #  wtff=Paramsfl[which(Paramsfl$catchment==4300027 ),]
  #  
  # wtf= datatwin[which(datatwin$outlets==4300027),]
  #plot(as.numeric(wtf[1,c(7:77)]))
  #Now aggregate by Hydroregions
  
  rmfuckers=unique(ParamsflSCF$catchment[which(ParamsflSCF$epsilonGPD>1.5)])
  
  data=data[-match(rmfuckers,data$outl2),]
  dataX=right_join(GHR_riv,data,by = c("outl2"="outl2"))
  dataX$uparea
  HRM=match(dataX$HydroRegions_raster_WGS84,GHshpp$Id)
  
  dataX$HR_id=GHshpp$CODEB[HRM]
  dataX$upaHR=GHshpp$SURF_KM2[HRM]
  
  dataXL=dataX[which(dataX$uparea>dataX$upaHR),]
  dataX=dataX[-which(dataX$uparea>dataX$upaHR)]
  pointagg=aggregate(list(Rchange=dataX$fillplot),
                     by = list(HydroR=dataX$HR_id),
                     FUN = function(x) c(mean=median(x,na.rm=T),dev=sd(x,na.rm=T),len=length(x),q1=quantile(x, 0.25, na.rm=T),q3=quantile(x, 0.75, na.rm=T)))
  pointagg <- do.call(data.frame, pointagg)
  
  ziz=full_join(GHshpp,pointagg,by=c("CODEB"="HydroR"))
  st_geometry(ziz)<-NULL
  #natch pointagg with hybasf
  pointplot=inner_join(HydroRsf,ziz,by= c("Id"="Id"))
  
  legend="10y flood specific discharge (l/s/km2)"
  title=paste0("10-years ",haz," Return Level")
  br=c(1,5,10,25,100,250,1000)
  labels=br
  limi=c(1,1200)
  
  pl2=ggplot(basemap) +
    geom_sf(fill="white")+
    geom_sf(data=pointplot,aes(fill=Rchange.mean,geometry=geometry),color="transparent",alpha=.5)+ 
    geom_sf(data=points,aes(col=fillplot,geometry=geometry),alpha=1,size=0.25,stroke=0,shape=15)+ 
    #geom_sf(fill=NA, color="grey") +
    geom_sf(data=pointplot,aes(geometry=geometry),color="darkgrey",fill=NA)+ 
    scale_fill_gradientn(
      colors=palet,
      breaks=br,labels=labels, limits=limi,trans="log",
      oob = scales::squish,na.value=colNA, name=legend)   +
    scale_color_gradientn(
      colors=palet,
      breaks=br,limits=limi,trans="log",
      oob = scales::squish,na.value=colNA, name=legend)   +
    #geom_sf(data=dt3,aes(geometry=geometry,fill=factor(plot)),color="gray12")+ 
    #scale_fill_manual(values=c("1"="purple","-1"="tan2"), name="IRES trend")+
    coord_sf(xlim = c(min(nco[,1]),max(nco[,1])), ylim = c(min(nco[,2]),max(nco[,2])))+
    labs(x="Longitude", y = "Latitude")+
    guides(color = guide_colourbar(barwidth = 20, barheight = .8), fill="none")+
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
          legend.key.size = unit(.8, "cm"))+
    ggtitle(title)
  
  min(points$fillplot)
  ggsave("D:/tilloal/Documents/LFRuns_utils/TrendAnalysis/plots/mapHR_abs_flscf.jpg", pl2, width=20, height=20, units=c("cm"),dpi=1000) 
  
  
  
  pointagg=aggregate(list(Rchange=dataX[,valcol+10]),
                     by = list(HydroR=dataX$HR_id),
                     FUN = function(x) c(mean=mean(x,na.rm=T)))
  pointagg <- do.call(data.frame, pointagg)
  
  
  pointagp=aggregate(list(Rchange=dataX[,valcol+10]),
                     by = list(HydroR=dataX$HR_id),
                     FUN = function(x) c(mean=mean(x,na.rm=T),sd=sd(x,na.rm=T),q10=quantile(x,0.1, na.rm=T),q90=quantile(x,0.9,na.rm=T)))
  pointagp <- do.call(data.frame, pointagp)
  
  
  pap=data.frame(t(pointagp))
  
  tpag=data.frame(t(pointagg))
  
  colnames(tpag)=tpag[1,]
  tpag=tpag[-1,]
  
  
  # Compute the correlation matrix
  cor_mat <- cor(tpag[ -c(1,2),])
  
  # Perform hierarchical clustering
  hclust_res <- hclust(as.dist(1 - cor_mat))
  # Plot the dendrogram
  plot(hclust_res, main = "Hierarchical Clustering Dendrogram")
  
  
  # Extract the cluster assignments
  
  library(clValid)
  print(dunn_index)
  wss <- vector()
  dunn_index <- vector()
  for (i in 1:10) {
    cluster_assignments <- cutree(hclust_res, k = i)
    wss[i] <- sum(cor_mat[cluster_assignments == cluster_assignments[1]]^2)
    dunn_index[i] <- dunn(as.dist(1 - cor_mat),cluster_assignments)
  }
  
  # Plot the WSS as a function of the number of clusters
  plot(1:10, wss, type = "b", xlab = "Number of Clusters", ylab = "WSS")
  
  plot(2:10, dunn_index[2:10], type = "b", xlab = "Number of Clusters", ylab = "WSS")
  
  cluster_assignments <- cutree(hclust_res, k = 3)
  
  #legend="Kendall tau"
  legend="Relative change (%)"
  tmpval=(datatwin[,valcol2])
  mkta=c()
  mksa=c()
  for (it in 1:length(pointagg[,1])){
    if (it%%1000==0) print(it)
    mks=NA
    mkt=NA
    miniTS=as.numeric(pointagg[it,])
    print(it)
    if (!is.na(miniTS[2])){
      #mk=MannKendall(miniTS[-1])
      mk2=mmkh(miniTS[-1],ci=0.95)
      mk=data.frame(tau=mk2[6],sl=mk2[2])
      
      #compute trend as well with sen.slope
      mkt=mk$tau
      mks=mk$sl
    }else{
      mkt=NA
      mks=NA
    }
    mkta=c(mkta,mkt)
    mksa=c(mksa,mks)
  }
  # br=c(-1,-.80,-.60,-.40,-.20,0,.20,.40,.60,.80,1)
  # limi=c(-1,1)
  trans=scales::modulus_trans(.6)
  colNA="darkgrey"
  #title="Trends in 100 years Re flood magnitude"
  
  
  #plot all the changes through all catchments
  pa_l <- reshape2::melt(pointagg, id.vars = "HydroR", variable.name = "Year", value.name = "value")
  
  
  tpmean=aggregate(list(rl=pa_l$value),
                   by = list(HydroR=pa_l$HydroR),
                   FUN = function(x) c(mean=mean(x,na.rm=T)))
  tpmean <- do.call(data.frame, tpmean)
  
  pa_l$value2=pa_l$value/tpmean$rl
  # Add the cluster assignments to the data frame
  pa_l$cluster <- factor(cluster_assignments)
  
  pointagg$mkta=mkta
  pointagg$sl=mksa
  meds <- c(by(pa_l$value, pa_l$HydroR, mean))
  
  plot(meds)
  pa2=pa_l[which(pa_l$cluster=="4"),]
  ggplot(pa2, aes(x = Year, y = value2, group = HydroR,col=cluster)) +
    geom_line() +
    scale_y_continuous(limits=c(0.5,1.5))+
    theme_minimal()
  
  plot(as.numeric(pointagg[1,-1]))
  
  pointagg$cluster=factor(cluster_assignments)
  
  #natch pointagg with hybasf
  pointsag=inner_join(GHshpp,pointagg,by=c("CODEB"="HydroR"))
  
  
  
  br=c(-30,-20,-10,0,10,20,30)
  labels=br
  limi=c(-30,30)
  tsize=16
  osize=12
  legend="Change in Qsp \n(l/s/km2)"
  
  ggplot(basemap) +
    geom_sf(fill="gray85",color="darkgrey",size=0.5)+
    geom_sf(data=pointsag,aes(fill=cluster,geometry=geometry),alpha=0.8,color="transparent")+
    #geom_sf(data=pointsag,aes(color=mkta,geometry=geometry),alpha=0.4,fill="transparent")+
    scale_color_gradientn(
      colors=palet,
      breaks=seq(-1,1,by=0.2), limits=c(-1,1),
      oob = scales::squish,na.value=colNA, name=legend) +
    coord_sf(xlim = c(min(nco[,1]),max(nco[,1])), ylim = c(min(nco[,2]),max(nco[,2])))+
    # geom_sf(data=datasl1,aes(col=fillplot,geometry=geometry),alpha=0.8,size=0.05,stroke=0,shape=1)+ 
    # geom_sf(data=datasl2,aes(col=fillplot,geometry=geometry),alpha=0.9,size=0.1,stroke=0,shape=1)+ 
    #geom_sf(data=data,aes(fill=fillplot,geometry=geometry),color="transparent")+ 
    #geom_sf(data=datasig_f,aes(size=siglvl),fill="grey",color="transparent",shape=21,alpha=0.1)+ 
    #scale_color_manual(values=c("blue" ="darkblue","red"="darkred"))+
    #geom_col_pattern(data=pointsag,aes(fill=siglvl,geometry=geometry),colour='black', pattern = 'circle') +
    #scale_fill_manual(values=c("tomato4","steelblue4"), name="Significant trends")+
    ggtitle(title)+
    # coord_sf(xlim = c(min(nco[,1]),max(nco[,1])), ylim = c(min(nco[,2]),max(nco[,2])))+
    # geom_point(data=event,aes(x=lon,y=lat,col=dr))+
    # scale_color_distiller(palette = "Spectral",
    #                       direction = -1, limit=c(0,max(eventp$dr)),oob = scales::squish)+
    # scale_color_gradientn(
    #   colors=palet,
    #   breaks=br,limits=limi,trans=trans,
    #   oob = scales::squish,na.value=colNA, name=legend)   +
    labs(x="Longitude", y = "Latitude")+
    guides(colour = guide_colourbar(barwidth = 1.5, barheight = 10),
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
          legend.key.size = unit(1, "cm"))
  
  
  
  
  pointsag$siglvl=0
  pointsag$siglvl[which(pointsag$sl<=0.1)]=1
  
  catsig=pointsag[which(pointsag$siglvl>0),]
  catsig <- st_transform(catsig, crs = 3035)
  ## create a grid of points
  grdpts <- sf::st_make_grid(catsig, what = "centers",cellsize = 40000)
  my.points <- sf::st_sf(grdpts)
  
  sf_use_s2(FALSE)
  pointsInside <- sf::st_join(x = my.points, y = catsig, left = FALSE)
  pointsInside$sign="positive"
  pointsInside$sign[which(pointsInside$mkta<0)]="negative"
  pointsag$siglvl=factor(pointsag$siglvl)
  
  # #TBC with the river pixels and better colors
  #   ggplot(basemap) +
  #   geom_sf(data=pointsInside,aes(geometry=geometry, fill=sign),alpha=.7,size=2,stroke=0,shape=21,color="black")+
  #   coord_sf(xlim = c(min(nco[,1]),max(nco[,1])), ylim = c(min(nco[,2]),max(nco[,2])))+
  #   scale_fill_manual(values=c("tomato4","steelblue4"), name="Significant trends")
  #   
  # ggplot(regiod)+
  # geom_sf(mapping=aes(geometry=geometry,group=name,fill=name))
  
  library(ggnewscale)
  br=c(-30,-20,-10,0,10,20,30)
  labels=br
  limi=c(-30,30)
  tsize=16
  osize=12
  legend="Change in Qsp \n(l/s/km2)"
  ocrap<-ggplot(basemap) +
    geom_sf(fill="gray85",color="darkgrey",size=0.5)+
    geom_sf(data=pointsag,aes(fill=mkta,geometry=geometry),alpha=0.4,color="transparent")+
    geom_sf(data=points,aes(col=fillplot,geometry=geometry),alpha=.9,size=0.15,stroke=0,shape=15)+ 
    #geom_sf(data=pointplot,aes(geometry=geometry),color="darkgrey",fill=NA)+ 
    scale_fill_gradientn(
      colors=palet,
      breaks=seq(-1,1,by=0.2), limits=c(-1,1),
      oob = scales::squish,na.value=colNA, name=legend) +
    guides(fill = "none")+
    new_scale_fill()+
    geom_sf(data=pointsInside,aes(geometry=geometry, fill=sign),alpha=.8,size=.6,stroke=0,shape=21,color="black")+
    #geom_sf(fill=NA, color="grey") +
    #geom_sf(regiod,mapping=aes(geometry=geometry,group=name),fill=NA,color="orange")+
    coord_sf(xlim = c(min(nco[,1]),max(nco[,1])), ylim = c(min(nco[,2]),max(nco[,2])))+
    # geom_sf(data=datasl1,aes(col=fillplot,geometry=geometry),alpha=0.8,size=0.05,stroke=0,shape=1)+ 
    # geom_sf(data=datasl2,aes(col=fillplot,geometry=geometry),alpha=0.9,size=0.1,stroke=0,shape=1)+ 
    #geom_sf(data=data,aes(fill=fillplot,geometry=geometry),color="transparent")+ 
    #geom_sf(data=datasig_f,aes(size=siglvl),fill="grey",color="transparent",shape=21,alpha=0.1)+ 
    #scale_color_manual(values=c("blue" ="darkblue","red"="darkred"))+
    #geom_col_pattern(data=pointsag,aes(fill=siglvl,geometry=geometry),colour='black', pattern = 'circle') +
    scale_fill_manual(values=c("tomato4","steelblue4"), name="Significant trends")+
    ggtitle(title)+
    # coord_sf(xlim = c(min(nco[,1]),max(nco[,1])), ylim = c(min(nco[,2]),max(nco[,2])))+
    # geom_point(data=event,aes(x=lon,y=lat,col=dr))+
    # scale_color_distiller(palette = "Spectral",
    #                       direction = -1, limit=c(0,max(eventp$dr)),oob = scales::squish)+
    scale_color_gradientn(
      colors=palet,
      breaks=br,limits=limi,trans=trans,
      oob = scales::squish,na.value=colNA, name=legend)   +
    labs(x="Longitude", y = "Latitude")+
    guides(colour = guide_colourbar(barwidth = 1.5, barheight = 10),
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
          legend.key.size = unit(1, "cm"))
  
  ggsave("D:/tilloal/Documents/LFRuns_utils/TrendAnalysis/plots/mapHR_final_flscf8020_sig.jpg", ocrap, width=20, height=20, units=c("cm"),dpi=1000) 
  
  
  #ocrap
  
  return(list(pl1, pl2))
}














### Comparison of two scenarios with mk test at catchment level ----

### Beam plot of change for Europe ----
plotbeam=FascPlotChange(datarSCF,outf,period=c(1950,2020),lims=c(-60, 60),q1=0.25,q2=0.75)

ggsave("D:/tilloal/Documents/LFRuns_utils/Figures/Change_alleurope_faisceau_drought.jpg",plotbeam, width=40, height=15, units=c("cm"),dpi=1500)


### Plot by biogeographic regions ----
biogeo <- read_sf(dsn = paste0(hydroDir,"/eea_3035_biogeo-regions_2016/BiogeoRegions2016_wag84.shp"))
biogeof=fortify(biogeo)
st_geometry(biogeof)<-NULL
biogeoregions=raster( paste0(hydroDir,"/eea_3035_biogeo-regions_2016/Biogeo_rasterized_wsg84.tif"))
Gbiogeoregions=as.data.frame(biogeoregions,xy=T)
biogeomatch=inner_join(biogeof,Gbiogeoregions,by= c("PK_UID"="Biogeo_rasterized_wsg84"))
biogeomatch$latlong=paste(round(biogeomatch$x,4),round(biogeomatch$y,4),sep=" ")
biogeo_rivers=inner_join(biogeomatch,outf, by="latlong")

unikbiog=unique(biogeo_rivers$PK_UID)
i=0
ridgechange=list()
faschange=list()
statregions=c()
for (b in unikbiog){
  i=i+1
  outf9=biogeo_rivers[which(biogeo_rivers$PK_UID==b),]
  rplots=plotHistoDates(outf9,datari,law="GPD",type,period=c(1950),parlist=Paramsfl,valuenames)
  plotv2=FascPlotChange(datar,outf9,period=c(1951,2019))
  faschange[[i]]=plotv2
  ggsave(paste0("D:/tilloal/Documents/LFRuns_utils/Figures/FaisceauChange_",i,".jpg"), plotv2, width=20, height=15, units=c("cm"),dpi=1000)
  ridgechange[[i]]=rplots[[2]]
  statregions=rbind(statregions,rplots[[1]])
  
}

#All results in a single file
plotkeep=faschange[c(1,2,3,5,6,7,8)]

layout_mat<-rbind(c(1,2),
                  c(3,4),
                  c(5,6),
                  c(8,9))
layout_mat

library(gridExtra)
plots<-marrangeGrob(plotkeep, layout_matrix=layout_mat)

ggsave("D:/tilloal/Documents/LFRuns_utils/Figures/FaisceauChange.jpg", plots, width=50, height=60, units=c("cm"),dpi=1000)
write.csv(statregions,file=paste0(hydroDir,"/Flood/stats_by_bgreg.csv"))


### Plot of the shape parameter in every catchment ----
pc=unique(ParamsflH$catchment)
parcat=ParamsflH[match(pc,ParamsflH$catchment),]
parcat$catchment
tsize=16
osize=16
basemap=w2
parplot=inner_join(parcat,UpArea,by=c("catchment"="outl2"))
parplot=full_join(outf,parcat,by = c("outlets"="catchment"))
parpl <- st_as_sf(parplot, coords = c("Var1", "Var2"), crs = 4326)
parpl <- st_transform(parpl, crs = 3035)


palet=c(hcl.colors(9, palette = "BrBG", alpha = NULL, rev = F, fixup = TRUE))
ggplot(basemap) +
  geom_sf(data=parpl,aes(col=epsilonGPD,geometry=geometry),alpha=1,size=0.1,stroke=0,shape=15)+ 
  geom_sf(fill=NA, color="grey") +
  coord_sf(xlim = c(min(nco[,1]),max(nco[,1])), ylim = c(min(nco[,2]),max(nco[,2])))+
  scale_color_gradientn(
    colors=palet, limits=c(-0.5,1),
    oob = scales::squish,na.value="grey", name="shape parameter")   +
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

ggsave("D:/tilloal/Documents/LFRuns_utils/Figures/ShapeParams.jpg", width=30, height=20, units=c("cm"),dpi=1000)

### Plot the 100Y RL at a given time ----
hazard="drought"
chyr=2000
Ychyr=paste0("Y",chyr)

parplot=inner_join(datar,outf,by = c("unikout"="outlets"))
parpl <- st_as_sf(parplot, coords = c("Var1", "Var2"), crs = 4326)
parpl <- st_transform(parpl, crs = 3035)
colplot=match(Ychyr,names(parpl))
names(parpl)[colplot]="Ysel"
palet=c(hcl.colors(9, palette = "YlGnBu", alpha = NULL, rev = T, fixup = TRUE))
max(parpl$Ysel, na.rm=T)

ggplot(basemap) +
  geom_sf(fill="white", color=NA) +
  geom_sf(data=parpl,aes(col=Ysel,geometry=geometry),alpha=1,size=0.25,stroke=0,shape=15)+ 
  geom_sf(fill=NA, color="grey") +
  coord_sf(xlim = c(min(nco[,1]),max(nco[,1])), ylim = c(min(nco[,2]),max(nco[,2])))+
  scale_color_gradientn(
    colors=palet, limits=c(0.1,200),
    breaks=c(0.1,1,10,100),
    trans="log",
    oob = scales::squish,na.value="darkorange", name="Q (m3/s)")   +
  guides(color = guide_colourbar(barwidth = 10, barheight = 1))+
  labs(x="Longitude", y = "Latitude")+
  theme(axis.title=element_text(size=tsize),
        panel.background = element_rect(fill = "aliceblue", colour = "grey1"),
        panel.border = element_rect(linetype = "solid", fill = NA, colour="black"),
        legend.title = element_text(size=tsize),
        legend.text = element_text(size=osize),
        legend.position = "bottom",
        panel.grid.major = element_line(colour = "grey70"),
        panel.grid.minor = element_line(colour = "grey90"),
        legend.key = element_rect(fill = "transparent", colour = "transparent"),
        legend.key.size = unit(.8, "cm"))+
  ggtitle(paste0("100 years ", hazard ," Return level - ", chyr))

ggsave("D:/tilloal/Documents/LFRuns_utils/Figures/drought_100yRL_2000.jpg", width=30, height=20, units=c("cm"),dpi=1000)



### Look ar change through different RPs ----
RLsave=c()
for (RP in c(5,10,20,50,100))
{
  print(RP)
  RLs=GPDLargeRLs(Paramsfl,RP)
  
  RLsave=cbind(RLsave,RLs$GPD$returnLevels)
}

RLsave=as.data.frame(RLsave)
Paramplus=cbind(Paramsfl,RLsave)

### Maximum year and work on multi-decadal cycle ----

datar=RLGPDfl
datart=datar[,-c(71,72)]
ymax=apply(datart, 1, which.max)
ymin=apply(datart, 1, which.min)

subdata=datart
datrd=subdata/subdata$Y1950*100-100

yp=c(1950:2019)
maxyears=yp[ymax]
minyears=yp[ymin]
hist(maxyears, breaks=7)
hist(minyears, breaks=7)

datart$maxyear=maxyears
datart$minyear=minyears
datari=inner_join(outf,datar,by = c("outlets"="unikout"))

points <- st_as_sf(datari, coords = c("Var1", "Var2"), crs = 4326)
points <- st_transform(points, crs = 3035)

palet=c(hcl.colors(9, palette = "PRGn", alpha = NULL, rev = T, fixup = TRUE))
ggplot(w2) +
  geom_sf(fill="gray85")+
  geom_sf(fill=NA, color="grey") +
  geom_sf(data=points,aes(col=factor(gt2),geometry=geometry),alpha=1,size=0.25,stroke=0,shape=15)+ 
  coord_sf(xlim = c(min(nco[,1]),max(nco[,1])), ylim = c(min(nco[,2]),max(nco[,2])))+
  
  scale_color_manual(values=c("orange","green","blue"))  +
  labs(x="Longitude", y = "Latitude")+
  theme(
    panel.background = element_rect(fill = "aliceblue", colour = "grey1"),
    panel.border = element_rect(linetype = "solid", fill = NA, colour="black"),
    legend.position = "right",
    panel.grid.major = element_line(colour = "grey70"),
    panel.grid.minor = element_line(colour = "grey90"),
    legend.key = element_rect(fill = "transparent", colour = "transparent"),
    legend.key.size = unit(.8, "cm"))

ggsave("D:/tilloal/Documents/LFRuns_utils/Figures/maxgroup.jpg", width=30, height=20, units=c("cm"),dpi=1000)


# Link to GMST ------------------------------------------------------------


#extraction of RL100 values per year
subdata=datari[,-c(1:6,77:81)]
datrd=subdata/subdata$Y1950*100-100

#get quantiles of change per year
dfplot= datrd %>% gather("Year","pixel")
dfplot$year=as.numeric(substr(dfplot$Year,2,5))

dfquantiles=aggregate(list(change=dfplot$pixel),
                      by = list(Year=dfplot$year),
                      FUN = function(x) c(mean=mean(x,na.rm=T),q1=quantile(x, 0.025, na.rm=T),q3=quantile(x, 0.975, na.rm=T)))

dfquantiles <- do.call(data.frame, dfquantiles)
names(dfquantiles)=c("Year", "mean","qlow","qhigh")

dfquantiles$Year[which.max(dfquantiles$mean)]

#Read GMST from HadCRUT
GMST=read.csv(file=paste0(hydroDir,"/timeseries/HadCRUT.5.0.1.0.csv"))

#Preindustrial temperature as defined in IPCC
PreIndustrial=mean(GMST$Anomaly..deg.C.[32:61])

#adjusting change to preIndustrial warming
GMST$Anomaly_PreInd=GMST$Anomaly..deg.C.-PreIndustrial

#30 years moving average of temperature
GMST$runningAnomaly=tsEvaNanRunningMean(GMST$Anomaly_PreInd, 30)
plot(GMST$Time,GMST$Anomaly_PreInd,type="o")
abline(h=GMST_sub$runningAnomaly[1])
lines(GMST$Time,GMST$runningAnomaly,col=2, lwd=2)

GMST$Time
tbound=c(1950,2019)
timestamp=seq(tbound[1],tbound[2])
yRL=dfquantiles$Year

GMST_sub=GMST[match(timestamp,GMST$Time),]
GMST_sub$meanchange=dfquantiles$mean
plot(GMST_sub$Time,GMST_sub$Anomaly_PreInd,type="o")
lines(GMST_sub$Time,GMST_sub$runningAnomaly,col=2, lwd=2)


plot(GMST_sub$runningAnomaly,GMST_sub$meanchange,type="o",ylim=c(-10,10))
lines(GMST_sub$runningAnomaly,dfquantiles$qlow)
lines(GMST_sub$runningAnomaly,dfquantiles$qhigh)


#Development part: approximating values to GW times
peakWL=approx(GMST_sub$Time,GMST_sub$runningAnomaly,xout=peakyears)$y
points(jitter(peakWL),pikos$value,pch=16, col="blue")
GWlevels=c(0.5,0.7,1)
WLyear=round(approx(GMST_sub$runningAnomaly,GMST_sub$Time,xout=GWlevels)$y)


# Seasonal analysis at large scale ----------------------------------------

load(file=paste0(hydroDir,"/Flood/params.flood.cal2.Rdata"))
load(file=paste0(hydroDir,"/Flood/peaks.flood.cal2.Rdata"))

#For visualisation purposes, I run the script on one pixel (loop is run on HPC)

#time settings
#max id of the time vector
tidm=25931
timeStamps6h=as.POSIXct(c(0:(tidm*4))*3600*6 -3600, origin="1950-01-03 12:00:00")
tsm=1/0.25
rmv=c(1:(60*tsm))
timeStamps=timeStamps6h[-rmv]
time=unique(as.Date(timeStamps))
time=time[-length(time)]
timeStamps=timeStamps6h[-rmv]


Nsq=42
print(Nsq)
Idstart=as.numeric(Nsq)*100000
nrspace=rspace[Nsq,]
outhybas=outletopen(hydroDir,outletname,nrspace)
outhybas$outlets=seq((Idstart+1),(Idstart+length(outhybas$outlets)))
outhybas$latlong=paste(round(outhybas$Var1,4),round(outhybas$Var2,4),sep=" ")



Nmat=which(!is.na(match(Peaksave$catch,outhybas$outlets)))
Pikday=Peaksave[Nmat,]
unikzero=unique(Pikday$catch)
id=1
cat(paste0("\n",id))
Rpix=Pikday[which(Pikday$catch == unikzero[id]),]

Rpix$d=yday(Rpix$time)
Rpix$year=year(Rpix$time)

#I don't want to look only at the peak but at all drought times

#vector of all drought days
tsEvents=c()
for (ev in 1:length(Rpix$value)){
  event=Rpix[ev,]
  #check that my new bricoled time stamp is valid
  event$time2=timeStamps[event$timeID]
  if (event$time2==event$time){
    seqdates=seq(event$tIDstart,event$tIDend)
    tev=timeStamps[seqdates]
    evn=rep(ev,length(tev))
    tsEv=data.frame(tev,evn)
    tsEvents=rbind(tsEvents,tsEv)
  }else{
    print("times not matching")
  }
}
tsEvents$day=yday(tsEvents$tev)
tsEvents$date=as.Date(tsEvents$tev)
tsEventsD=aggregate(list(ev=tsEvents$evn),
                    by = list(event=tsEvents$evn, date=tsEvents$date),
                    FUN = function(x) c(len=length(x)))
tsEventsD <- do.call(data.frame, tsEventsD)
tsEventsD$day=yday(tsEventsD$date)
tsEventsD$ISev=1
dayvec=data.frame(day=seq(1,366))
nyears=length(unique(year(time)))
tw=30
windowSize=tw*365.25
timestampsD=data.frame(date=time)

tsEventsF=right_join(tsEventsD,timestampsD,by="date")
tsEventsF=tsEventsF[order(tsEventsF$date),]
tsEventsF$ISev[which(is.na(tsEventsF$ISev))]=0

if (length(which(diff(tsEventsF$date)==0))) tsEventsF=tsEventsF[-which(diff(tsEventsF$date)==0),]

series=tsEventsF$day
rseason=RunningSeason(series, timestampsD$date,windowSize, nyears, 50)

rseason$diff=rseason$rnseas-rseason$rnseas[1]
if (length(which(is.na(rseason$diff)))<1){
  rseason$diff[which(rseason$diff>182)]=rseason$diff[which(rseason$diff>182)]-365.25
  rseason$diff[which(rseason$diff<=-182)]=365.25+rseason$diff[which(rseason$diff<=-182)]
}
rseason$catch=rep(unikzero[id],length(rseason$rnseas))

#save example plot
rseason[[2]]
ggsave("flood_seasonality_pix3.jpg", width=30, height=21, units=c("cm"),dpi=1000)


rseasonT=c()
outf=c()
hydroDir<-("D:/tilloal/Documents/LFRuns_utils/data")
outlets="Rnet"
lf=list.files(path = paste0(hydroDir,"/Flood/Peaks/v2"), full.names = TRUE, recursive = TRUE)
lf
for (file in lf){
  load(file)
  fils=sub(".*/", "", file)
  tt=unlist(strsplit(fils, "[_]"))
  Nsq=as.numeric(tt[2])
  out.type=sub("\\_.*", "", fils)
  print(Nsq)
  
  
  
  rspace= read.csv(paste0(hydroDir,"/subspace_efas.csv"))
  rspace=rspace[,-1]
  nrspace=rspace[Nsq,]
  #outletname="outletsv8_hybas07_01min"
  #outletname="outlets_hybas09_01min"
  outletname="efas_rnet_100km_01min"
  
  outhybas=outletopen(hydroDir,outletname,nrspace)
  Idstart=as.numeric(Nsq)*100000
  if (length(outhybas$outlets)>0){
    outhybas$outlets=seq((Idstart+1),(Idstart+length(outhybas$outlets)))
    outhybas$latlong=paste(round(outhybas$Var1,4),round(outhybas$Var2,4),sep=" ")
    # outcut=which(!is.na(match(outhybas$outlets,parlist$catchment)))
    # zebi=unique(parlist$catchment)
    # outhloc=outhybas[outcut,]
    
    outf=rbind(outf,outhybas)
  }
  
  
  rseasonT=rbind(rseasonT,rseasonF)
  
  
}







#comput mean, median, q025 and q975
summarystat=aggregate(list(diff=rseasonT$rnseas),
                      by = list(year=rseasonT$time),
                      FUN = function(x) c(mean=mean(x,na.rm=T),median=median(x,na.rm=T),q025=quantile(x,0.025,na.rm=T),q975=quantile(x,0.975,na.rm=T)))
summarystat=do.call(data.frame, summarystat)

rseasonT$year=year(rseasonT$time)
rseasonP=rseasonT[which(rseasonT$year==2005),]

seasonplot=full_join(rseasonP,outf, by=c("catch"="outlets"))
seasonplot$sq=round(seasonplot$catch/100000)

seasonP <- st_as_sf(seasonplot, coords = c("Var1", "Var2"), crs = 4326)
seasonP <- st_transform(seasonP, crs = 3035)

cord.dec=outf[,c(2,3)]
cord.dec = SpatialPoints(cord.dec, proj4string=CRS("+proj=longlat"))
cord.UTM <- spTransform(cord.dec, CRS("+init=epsg:3035"))
nco=cord.UTM@coords


####################### Seasonal Plots #######################


#Plots for seasonality

## Mean flood date plot
direction_labeller <- function(x){
  ifelse(x %% 30 == 0, rev(c('Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct',"Nov",'Dec'))[(as.integer(x/30) %% 13)], '')
}
my_colors=(c("#2E22EA","#9E3DFB","#F86BE2","#FCCE7B","#C4E416","#4BBA0F","#447D87","#2C24E9"))
x=seq(360,30,-30)
ptn=direction_labeller(seq(360,30,-30))
ptnn=rev(ptn)
paletx=c(hcl.colors(8, palette = "PurpOr", alpha = NULL, rev = TRUE, fixup = TRUE))
tsize=16
osize=16
ggplot(w2) +
  geom_sf(fill="white")+
  geom_sf(data=seasonP,aes(col=rnp1,geometry=geometry),alpha=1,size=0.25,stroke=0,shape=15)+ 
  geom_sf(fill=NA, color="grey") +
  coord_sf(xlim = c(min(nco[,1]),max(nco[,1])), ylim = c(min(nco[,2]),max(nco[,2])))+
  scale_color_gradientn(
    colors=my_colors,
    breaks=seq(15,345,30),
    label=ptn,
    na.value="grey95",
    limits=c(0,366), name="Date")+
  labs(x="Longitude", y = "Latitude")+
  guides(color = guide_colourbar(barwidth = 2, barheight = 15,reverse=T))+
  theme(axis.title=element_text(size=tsize),
        panel.background = element_rect(fill = "aliceblue", colour = "grey1"),
        panel.border = element_rect(linetype = "solid", fill = NA, colour="black"),
        legend.title = element_text(size=tsize),
        legend.text = element_text(size=osize),
        panel.grid.major = element_line(colour = "grey70"),
        panel.grid.minor = element_line(colour = "grey90"),
        legend.key = element_rect(fill = "transparent", colour = "transparent"),
        legend.key.size = unit(.8, "cm"))+
  ggtitle("mean flood date")
ggsave( filename="D:/tilloal/Documents/LFRuns_utils/Figures/mean_flood_date_1965.jpg", width=16.3, height=15, units="cm",dpi=3000)



## Concentration of flood timing plot
paletx=c(hcl.colors(8, palette = "YlOrRd", alpha = NULL, rev = TRUE, fixup = TRUE))

ggplot(w2) +
  geom_sf(fill="white")+
  geom_sf(data=seasonP,aes(col=rncon,geometry=geometry),alpha=1,size=0.25,stroke=0,shape=15)+ 
  geom_sf(fill=NA, color="grey") +
  coord_sf(xlim = c(min(nco[,1]),max(nco[,1])), ylim = c(min(nco[,2]),max(nco[,2])))+
  scale_color_gradientn(
    colors=paletx,
    limits=c(0,1),
    na.value="grey95",oob = scales::squish, name="R")+
  labs(x="Longitude", y = "Latitude")+
  guides(color = guide_colourbar(barwidth = 1, barheight = 10,reverse=F))+
  theme(axis.title=element_text(size=tsize),
        panel.background = element_rect(fill = "aliceblue", colour = "grey1"),
        panel.border = element_rect(linetype = "solid", fill = NA, colour="black"),
        legend.title = element_text(size=tsize),
        legend.text = element_text(size=osize),
        panel.grid.major = element_line(colour = "grey70"),
        panel.grid.minor = element_line(colour = "grey90"),
        legend.key = element_rect(fill = "transparent", colour = "transparent"),
        legend.key.size = unit(.8, "cm"))+
  ggtitle("Concentration around mean flood date")

ggsave( filename="D:/tilloal/Documents/LFRuns_utils/Figures/concentration_flood_season_2005x.jpg", width=16.3, height=15, units="cm",dpi=3000)



## Change in flood timing: 1950-2020
rseason1950=rseasonT[which(rseasonT$year==1965),]
rseason2020=rseasonT[which(rseasonT$year==2005),]

rseason2020$rnseas1950=rseason1950$rnseas
rseason2020$rnnp1950=rseason1950$rnnp
rseason2020$rnp11950=rseason1950$rnp1
rseason2020$rnp21950=rseason1950$rnp2
rseason2020$diff=rseason2020$rnseas-rseason1950$rnseas

rseason2020$diff[which(rseason2020$diff>182)]=rseason2020$diff[which(rseason2020$diff>182)]-365.25
rseason2020$diff[which(rseason2020$diff<=-182)]=365.25+rseason2020$diff[which(rseason2020$diff<=-182)]


rseason2020=full_join(rseason2020,outf, by=c("catch"="outlets"))
rseason2020$sq=round(rseason2020$catch/100000)
rseason2020P <- st_as_sf(rseason2020, coords = c("Var1", "Var2"), crs = 4326)
rseason2020P <- st_transform(rseason2020P, crs = 3035)

#remove differences that are too big for visualisation
rseason2020P$diff2=rseason2020P$diff
rseason2020P$diff2[which(abs(rseason2020P$diff)>60)]=NA
paletx=c(hcl.colors(8, palette = "Spectral", alpha = NULL, rev = TRUE, fixup = TRUE))

ggplot(w2) +
  geom_sf(fill="white")+
  geom_sf(data=rseason2020P,aes(col=diff,geometry=geometry),alpha=1,size=0.25,stroke=0,shape=15)+ 
  geom_sf(fill=NA, color="grey") +
  coord_sf(xlim = c(min(nco[,1]),max(nco[,1])), ylim = c(min(nco[,2]),max(nco[,2])))+
  scale_color_gradientn(
    colors=paletx,
    limits=c(-60,60),
    na.value="grey95",oob = scales::squish, name=" days")+
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
  ggtitle("change in mean flood timing (2005-1965)")


ggsave( filename="D:/tilloal/Documents/LFRuns_utils/Figures/change_floodTiming_19652005x7.jpg", width=16.3, height=15, units="cm",dpi=3000)


## Aggregate by catchment to decipher a spatial pattern

hist(rseason2020$diff,breaks=365)
mean(rseason2020$diff,na.rm=T)
rseason2020$diff2=rseason2020$diff
rseason2020$diff2[which(abs(rseason2020$diff)>60)]=NA
rseason2020x=inner_join(rseason2020,Catamere07,by=c("latlong"="llcoord"))
rseasonAgg=aggregate(list(Tchange=rseason2020x$diff),
                     by = list(HYBAS_ID=rseason2020x$HYBAS_ID),
                     FUN = function(x) c(med=median(x,na.rm=T),dev=sd(x,na.rm=T),len=length(x)))
rseasonAgg <- do.call(data.frame, rseasonAgg)

rseasonAgg=inner_join(hybasf7,rseasonAgg,by= "HYBAS_ID")

paletx=c(hcl.colors(8, palette = "Spectral", alpha = NULL, rev = TRUE, fixup = TRUE))
ggplot(w2) +
  geom_sf(fill="white")+
  geom_sf(data=rseasonAgg,aes(fill=Tchange.med,geometry=geometry),color="transparent")+ 
  geom_sf(fill=NA, color="gray") +
  coord_sf(xlim = c(min(nco[,1]),max(nco[,1])), ylim = c(min(nco[,2]),max(nco[,2])))+
  scale_fill_gradientn(
    colors=paletx,
    limits=c(-60,60),
    name=" days",
    oob = scales::squish,na.value="gray")   +
  labs(x="Longitude", y = "Latitude")+
  guides(fill = guide_colourbar(barwidth = 2, barheight = 10))+
  theme(axis.title=element_text(size=tsize),
        panel.background = element_rect(fill = "aliceblue", colour = "grey1"),
        panel.border = element_rect(linetype = "solid", fill = NA, colour="black"),
        legend.title = element_text(size=tsize),
        legend.text = element_text(size=osize),
        panel.grid.major = element_line(colour = "grey70"),
        panel.grid.minor = element_line(colour = "grey90"),
        legend.key = element_rect(fill = "transparent", colour = "transparent"),
        legend.key.size = unit(.8, "cm"))

ggsave( filename="D:/tilloal/Documents/LFRuns_utils/Figures/change_floodTiming_19652005_cat.jpg", width=16.3, height=15, units="cm",dpi=3000)


#Change in number of flood season: 1950-2020

rseason2020$chns=rseason2020$rnnp-rseason2020$rnnp1950
rseason2020P <- st_as_sf(rseason2020, coords = c("Var1", "Var2"), crs = 4326)
rseason2020P <- st_transform(rseason2020P, crs = 3035)

paletx=c(hcl.colors(8, palette = "RdBu", alpha = NULL, rev = TRUE, fixup = TRUE))
ggplot(w2) +
  geom_sf(fill="white")+
  geom_sf(data=rseason2020P,aes(col=chns,geometry=geometry),alpha=1,size=0.1,stroke=0,shape=15)+ 
  geom_sf(fill=NA, color="grey") +
  coord_sf(xlim = c(min(nco[,1]),max(nco[,1])), ylim = c(min(nco[,2]),max(nco[,2])))+
  scale_color_gradientn(
    colors=paletx,
    limits=c(-1,1),
    na.value="grey95",oob = scales::squish, name=" change number of flood seasons")+
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
  ggtitle("change in mean flood seasons (2005-1965)")


my_colors2=rev(c("1"="#8fce00","2"="#4d2aa7","3"="#990000","4"="#4d2aa7"))
ggplot(w2) +
  geom_sf(fill="white")+
  geom_sf(data=seasonP,aes(color=factor(rnnp),geometry=geometry),alpha=1,size=0.3,stroke=0,shape=15)+ 
  geom_sf(fill=NA, color="grey") +
  coord_sf(xlim = c(min(nco[,1]),max(nco[,1])), ylim = c(min(nco[,2]),max(nco[,2])))+
  scale_color_manual(
    values=my_colors2,
    na.value="grey",
    breaks=c("4","3","2","1"),
    name=" # flood \n seasons")+
  guides(colour = guide_legend(override.aes = list(size = 10)))+
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

ggsave( filename="D:/tilloal/Documents/LFRuns_utils/Figures/NfloodSeasons_2005v2.jpg", width=16.3, height=15, units="cm",dpi=1600)

#Peaks in seasons (spring, summer autumn, winter) -------------------------


## different approaches to season clustering ----

seasonP$spring=NA
seasonP$summer=NA
seasonP$autumn=NA
seasonP$winter=NA

seasonP$spring[which(seasonP$rnp1>80 & seasonP$rnp1<172)]=1
seasonP$spring[which(seasonP$rnp2>80 & seasonP$rnp2<172)]=1
seasonP$spring[which(seasonP$rnp3>80 & seasonP$rnp3<172)]=1
seasonP$spring[which(seasonP$rnp4>80 & seasonP$rnp4<172)]=1

seasonP$summer[which(seasonP$rnp1>171 & seasonP$rnp1<266)]=1
seasonP$summer[which(seasonP$rnp2>171 & seasonP$rnp2<266)]=1
seasonP$summer[which(seasonP$rnp3>171 & seasonP$rnp3<266)]=1
seasonP$summer[which(seasonP$rnp4>171 & seasonP$rnp4<266)]=1

seasonP$autumn[which(seasonP$rnp1>265 & seasonP$rnp1<356)]=1
seasonP$autumn[which(seasonP$rnp2>265 & seasonP$rnp2<356)]=1
seasonP$autumn[which(seasonP$rnp3>265 & seasonP$rnp3<356)]=1
seasonP$autumn[which(seasonP$rnp4>265 & seasonP$rnp4<356)]=1

seasonP$winter[which(seasonP$rnp1>355 | seasonP$rnp1<81)]=1
seasonP$winter[which(seasonP$rnp2>355 | seasonP$rnp2<81)]=1
seasonP$winter[which(seasonP$rnp3>355 | seasonP$rnp3<81)]=1
seasonP$winter[which(seasonP$rnp4>355 | seasonP$rnp4<81)]=1


rseason1950=rseasonT[which(rseasonT$year==1965),]

rseason2020=rseasonT[which(rseasonT$year==2005),]
rseason1950$spring=NA
rseason1950$summer=NA
rseason1950$autumn=NA
rseason1950$winter=NA
rseason1950$spring[which(rseason1950$rnp1>80 & rseason1950$rnp1<172)]=rseason1950$rnp1[which(rseason1950$rnp1>80 & rseason1950$rnp1<172)]
rseason1950$spring[which(rseason1950$rnp2>80 & rseason1950$rnp2<172)]=rseason1950$rnp2[which(rseason1950$rnp2>80 & rseason1950$rnp2<172)]

rseason1950$summer[which(rseason1950$rnp1>171 & rseason1950$rnp1<266)]=rseason1950$rnp1[which(rseason1950$rnp1>171 & rseason1950$rnp1<266)]
rseason1950$summer[which(rseason1950$rnp2>171 & rseason1950$rnp2<266)]=rseason1950$rnp2[which(rseason1950$rnp2>171 & rseason1950$rnp2<266)]

rseason1950$autumn[which(rseason1950$rnp1>265 & rseason1950$rnp1<356)]=rseason1950$rnp1[which(rseason1950$rnp1>265 & rseason1950$rnp1<356)]
rseason1950$autumn[which(rseason1950$rnp2>265 & rseason1950$rnp2<356)]=rseason1950$rnp2[which(rseason1950$rnp2>265 & rseason1950$rnp2<356)]

rseason1950$winter[which(rseason1950$rnp1>355 | rseason1950$rnp1<81)]=rseason1950$rnp1[which(rseason1950$rnp1>355 | rseason1950$rnp1<81)]
rseason1950$winter[which(rseason1950$rnp2>355 | rseason1950$rnp2<81)]=rseason1950$rnp2[which(rseason1950$rnp2>355 | rseason1950$rnp2<81)]


rseason2020$spring=NA
rseason2020$summer=NA
rseason2020$autumn=NA
rseason2020$winter=NA
rseason2020$spring[which(rseason2020$rnp1>80 & rseason2020$rnp1<172)]=rseason2020$rnp1[which(rseason2020$rnp1>80 & rseason2020$rnp1<172)]
rseason2020$spring[which(rseason2020$rnp2>80 & rseason2020$rnp2<172)]=rseason2020$rnp2[which(rseason2020$rnp2>80 & rseason2020$rnp2<172)]

rseason2020$summer[which(rseason2020$rnp1>171 & rseason2020$rnp1<266)]=rseason2020$rnp1[which(rseason2020$rnp1>171 & rseason2020$rnp1<266)]
rseason2020$summer[which(rseason2020$rnp2>171 & rseason2020$rnp2<266)]=rseason2020$rnp2[which(rseason2020$rnp2>171 & rseason2020$rnp2<266)]

rseason2020$autumn[which(rseason2020$rnp1>265 & rseason2020$rnp1<356)]=rseason2020$rnp1[which(rseason2020$rnp1>265 & rseason2020$rnp1<356)]
rseason2020$autumn[which(rseason2020$rnp2>265 & rseason2020$rnp2<356)]=rseason2020$rnp2[which(rseason2020$rnp2>265 & rseason2020$rnp2<356)]

rseason2020$winter[which(rseason2020$rnp1>355 | rseason2020$rnp1<81)]=rseason2020$rnp1[which(rseason2020$rnp1>355 | rseason2020$rnp1<81)]
rseason2020$winter[which(rseason2020$rnp2>355 | rseason2020$rnp2<81)]=rseason2020$rnp2[which(rseason2020$rnp2>355 | rseason2020$rnp2<81)]


rseason2020$dwint=rseason2020$winter-rseason1950$winter

rseason2020$dspring[which(rseason2020$dspring>182)]=rseason2020$dspring[which(rseason2020$dspring>182)]-365.25
rseason2020$dspring[which(rseason2020$dspring<=-182)]=365.25+rseason2020$dspring[which(rseason2020$dspring<=-182)]

rseason2020$dspring=rseason2020$spring-rseason1950$spring

rseason2020$dwint[which(rseason2020$dwint>182)]=rseason2020$dwint[which(rseason2020$dwint>182)]-365.25
rseason2020$dwint[which(rseason2020$dwint<=-182)]=365.25+rseason2020$dwint[which(rseason2020$dwint<=-182)]

rseason2020$dsum=rseason2020$summer-rseason1950$summer

rseason2020$dsum[which(rseason2020$dsum>182)]=rseason2020$dsum[which(rseason2020$dsum>182)]-365.25
rseason2020$dsum[which(rseason2020$dsum<=-182)]=365.25+rseason2020$dsum[which(rseason2020$dsum<=-182)]

rseason2020$daut=rseason2020$autumn-rseason1950$autumn

rseason2020$daut[which(rseason2020$daut>182)]=rseason2020$daut[which(rseason2020$daut>182)]-365.25
rseason2020$daut[which(rseason2020$daut<=-182)]=365.25+rseason2020$daut[which(rseason2020$daut<=-182)]


seasonplot=full_join(rseason2020,outf, by=c("catch"="outlets"))
seasonplot$sq=round(seasonplot$catch/100000)
seasonP <- st_as_sf(seasonplot, coords = c("Var1", "Var2"), crs = 4326)
seasonP <- st_transform(seasonP, crs = 3035)

### Change in mean flood timing ----
paletx=c(hcl.colors(8, palette = "Spectral", alpha = NULL, rev = TRUE, fixup = TRUE))
ggplot(w2) +
  geom_sf(fill="white")+
  geom_sf(data=seasonP,aes(col=daut,geometry=geometry),alpha=1,size=0.25,stroke=0,shape=15)+ 
  geom_sf(fill=NA, color="grey") +
  coord_sf(xlim = c(min(nco[,1]),max(nco[,1])), ylim = c(min(nco[,2]),max(nco[,2])))+
  scale_color_gradientn(
    colors=paletx,
    limits=c(-60,60),
    na.value="grey95",oob = scales::squish, name=" days")+
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
  ggtitle("change in mean flood timing (2005-1965)")


ggsave( filename="D:/tilloal/Documents/LFRuns_utils/Figures/change_floodTiming_19652005x7.jpg", width=16.3, height=15, units="cm",dpi=3000)

### Aggregation to Hybas catchments ----

seasonPx=inner_join(seasonplot,Catamere07,by=c("latlong"="llcoord"))
rseasonAgg=aggregate(list(Tchange=seasonPx$dsum),
                     by = list(HYBAS_ID=seasonPx$HYBAS_ID),
                     FUN = function(x) c(med=mean(x,na.rm=T),dev=sd(x,na.rm=T),len=length(x)))
rseasonAgg <- do.call(data.frame, rseasonAgg)

rseasonAgg=inner_join(hybasf7,rseasonAgg,by= "HYBAS_ID")

paletx=c(hcl.colors(8, palette = "Spectral", alpha = NULL, rev = TRUE, fixup = TRUE))
ggplot(w2) +
  geom_sf(fill="white")+
  geom_sf(data=rseasonAgg,aes(fill=Tchange.med,geometry=geometry),color="transparent")+ 
  geom_sf(fill=NA, color="gray") +
  coord_sf(xlim = c(min(nco[,1]),max(nco[,1])), ylim = c(min(nco[,2]),max(nco[,2])))+
  scale_fill_gradientn(
    colors=paletx,
    limits=c(-60,60),
    name=" days",
    oob = scales::squish,na.value="gray")   +
  labs(x="Longitude", y = "Latitude")+
  guides(fill = guide_colourbar(barwidth = 2, barheight = 10))+
  theme(axis.title=element_text(size=tsize),
        panel.background = element_rect(fill = "aliceblue", colour = "grey1"),
        panel.border = element_rect(linetype = "solid", fill = NA, colour="black"),
        legend.title = element_text(size=tsize),
        legend.text = element_text(size=osize),
        panel.grid.major = element_line(colour = "grey70"),
        panel.grid.minor = element_line(colour = "grey90"),
        legend.key = element_rect(fill = "transparent", colour = "transparent"),
        legend.key.size = unit(.8, "cm"))

ggsave( filename="D:/tilloal/Documents/LFRuns_utils/Figures/change_floodTiming_19652005_cat.jpg", width=16.3, height=15, units="cm",dpi=3000)


mySplot=plotSeasonFlood(seasonP,nco,tsize,osize)
seasons=c("summer","autumn","winter","spring")
library(ggpubr)
saveplot1=ggarrange(mySplot[[1]],NULL, mySplot[[2]],mySplot[[3]],NULL, mySplot[[4]],
                    ncol = 3, nrow = 2,widths = c(1,-0.3,1,1,-0.3,1))

ggsave( filename="D:/tilloal/Documents/LFRuns_utils/Figures/flood_seasonality6.jpg", saveplot1, width=29.3, height=21, units="cm",dpi=2000)

### Histogram of number of pixels experiencing peaks in each season in period 1 and period 2 ----

seqseason=c("spring","summer","autumn","winter")
lseaf=c()
for (i in (1:4)){
  lseas=data.frame(year=c(1965,2005),season=rep(seqseason[i],2),length=c(length(which(!is.na(rseason1950[,11+i])))/length(rseason1950[,11+i]),length(which(!is.na(rseason2020[,11+i])))/length(rseason2020[,11+i])))
  
  lseaf=rbind(lseaf,lseas)
}
ggplot(lseaf, aes(fill=factor(year), x=factor(season), y=length)) +
  # scale_fill_manual(values = cols, name = "Socioeconomic \n Pathways")+
  # labs(title=tit1,x="Warmimg levels (?C)", y = y1)+
  geom_bar(stat = "identity", position="dodge2")+
  geom_vline(xintercept = 1.5,size=1, col="black",lty="dotted")+
  geom_vline(xintercept = 2.5,size=1, col="black",lty="dotted")+
  geom_vline(xintercept = 3.5,size=1, col="black",lty="dotted")+
  # scale_color_manual(values = cols, name = "Socioeconomic \n Pathways",guide="none")+ 
  theme(axis.title=element_text(size=14),
        panel.background = element_rect(fill = "transparent", colour = "grey10"),
        panel.border = element_rect(linetype = "solid", fill = NA, colour="black"),
        legend.title = element_text(size=14),
        legend.text = element_text(size=12),
        panel.grid.major.y = element_line(colour = "grey70"),
        panel.grid.minor.y = element_line(colour = "grey90"),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        legend.key = element_rect(fill = "transparent", colour = "transparent"),
        legend.key.size = unit(.8, "cm"))

