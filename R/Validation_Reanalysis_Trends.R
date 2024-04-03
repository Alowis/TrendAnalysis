library(raster)
library(rgdal)
library(sf)
library(tidyverse)
library(ncdf4)
library(sp)
library(reshape)
library(MASS)
library(ggplot2)
library(viridis)
library(rnaturalearth)
library(hydroGOF)
source("~/LFRuns_utils/TS-EVA/functions_dvlp.R")

#  Function declaration ---------------------------------------------------
get_density <- function(x, y, ...) {
  dens <- MASS::kde2d(x, y, ...)
  ix <- findInterval(x, dens$x)
  iy <- findInterval(y, dens$y)
  ii <- cbind(ix, iy)
  return(dens$z[ii])
}
log10_minor_break <- function(...){
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
disNcopenloc=function(ncd,outloc,idc){
  name.vb=names(ncd[['var']])
  namev="dis"
  time <- ncvar_get(ncd,"time")
  lt=length(time)
  
  #timestamp corretion
  name.lon = "lon"
  name.lat="lat"
  londat = ncvar_get(ncd,name.lon) 
  llo=length(londat)
  latdat = ncvar_get(ncd,name.lat) 
  lla=length(latdat)
  
  idm=1+(idc-1)*lt
  outloc$llcoord[idc]
  outloc$LisfloodX[idc]
  start=c(outloc$idlo[idc],outloc$idla[idc],1)
  count=c(1,1,lt)
  outlets = ncvar_get(ncd,namev,start = start, count= count)
  outlets=as.vector(outlets)
  plot(outlets)
  print(outloc$StationID[idc])
  outid=rep(outloc$StationID[idc],length(time))
  #lonlatloop=expand.grid(c(1:lla),c(1:lt))
  lon=rep(londat[start[1]],length(time))
  lat=rep(latdat[start[2]],length(time))
  outll=data.frame(outlets,outid,lon,lat,time)
  
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
  
  outll$idlalo=paste(outll$idlo,outll$idla,sep=" ")
  outfinal=inner_join(outll, Sloc_final, by="idlalo")
  return (outfinal)
}
SpatialSkillPlot<-function(points,metric,inputfile){
  #Plot parameters
  cord.dec=points[,c(1,2)]
  cord.dec = SpatialPoints(cord.dec, proj4string=CRS("+proj=longlat"))
  cord.UTM <- spTransform(cord.dec, CRS("+init=epsg:3035"))
  nco=cord.UTM@coords
  world <- ne_countries(scale = "medium", returnclass = "sf")
  Europe <- world[which(world$continent == "Europe"),]
  e2=st_transform(Europe,  crs=3035)
  w2=st_transform(world,  crs=3035)
  tsize=12
  osize=12
  basemap=w2
  
  myfile=inputfile
  skill=read.csv(paste0(valid_path,myfile))
  if (length(skill[1,])>2){
    years=c(1950:2020)
    skill=skill[-73]
    names(skill)[c(2:72)]=years
  }
  names(skill)[1]="Station_ID"
  
  #Only keep valid stations
  valid=match(points$V1, skill$Station_ID)
  skillv=skill[valid,]
  skillm<-melt(skillv,id.vars ="Station_ID")
  skillm2=skillm[which(skillm$value!=0),]
  skillm2$value[which(is.infinite(skillm2$value))]=NA
  
  #average
  if (length(skill[1,])>2){
    dat_skill=aggregate(list(skill=skillm2$value),
                        by = list(Station_ID=skillm2$Station_ID),
                        FUN = function(x) c(mean= mean(x,na.rm=T),median=median(x,na.rm=T),max=max(x,na.rm=T)))
    
    dat_skill <- do.call(data.frame, dat_skill)
  }else{
    dat_skill=skillm2[,c(1,3)]
  }
  
  
  points=inner_join(points,dat_skill,by=c("V1"="Station_ID"))
  if (length(skill[1,])>2){
    points$skill=points$skill.mean
    
  }else{
    points$skill=points$value
  }
  
  #For correlation
  if (metric=="r"){
    palet=c(hcl.colors(9, palette = "YlGnBu", alpha = NULL, rev = T, fixup = TRUE))
    limi=c(0,1)
    var="r"
  }
  #For Bias
  if (metric=="b"){
    palet=c(hcl.colors(9, palette = "RdBu", alpha = NULL, rev = F, fixup = TRUE))
    limi=c(0,2)
    var="b"
  }
  #For Variability
  if (metric=="v"){
    palet=c(hcl.colors(9, palette = "RdBu", alpha = NULL, rev = F, fixup = TRUE))
    #palet = c('#f61f0f','#f4cccc','#f6ebeb','#ffffff','#e4edf5','#91bfdb','#1c87d8')
    limi=c(0,2)
    var="v"
  }
  if (metric=="kge"){
    colorz = c('#d73027','orange','#fee090','lightblue','royalblue',"darkblue")
    kgelabs=c("< -0.41","-0.41 - 0","0 - 0.2", "0.2 - 0.5","0.5 - 0.75", ">0.75")
    points$skill=as.numeric(points$skill)
    points$kgecode=0
    points$kgecode[which(is.na(points$skill))]=NA
    points$kgecode[which(points$skill>(-0.41) & points$skill<=0)]=1
    points$kgecode[which(points$skill>0 & points$skill<=0.2)]=2
    points$kgecode[which(points$skill>0.2 & points$skill<=0.5)]=3
    points$kgecode[which(points$skill>0.5 & points$skill<=0.75)]=4
    points$kgecode[which(points$skill>0.75)]=5
    var="kge"
  }
  parpl <- st_as_sf(points, coords = c("Var1", "Var2"), crs = 4326)
  parpl <- st_transform(parpl, crs = 3035)
  parpef=parpl[which(parpl$calib==TRUE),]
  if (metric =="kge"){
    tsize=size=12
    map=ggplot(basemap) +
      geom_sf(fill="gray95", color=NA) +
      geom_sf(data=parpl,aes(geometry=geometry),col="black",alpha=1,size=1.4,stroke=0.05,shape=1)+ 
      geom_sf(data=parpl,aes(geometry=geometry,colour=factor(kgecode)),alpha=.9,size=1.4,stroke=0,shape=16)+ 
      geom_sf(data=parpef,aes(geometry=geometry),col="black",alpha=1,size=1.4,stroke=0.1,shape=8)+ 
      geom_sf(fill=NA, color="grey20") +
      scale_x_continuous(breaks=seq(-30,40, by=5)) +
      labs(x="Longitude", y = "Latitude")+
      coord_sf(xlim = c(min(nco[,1]),max(nco[,1])), ylim = c(min(nco[,2]),max(nco[,2])))+
      scale_color_manual(values = colorz, labels= kgelabs, name="KGE")   +
      labs(x="Longitude", y = "Latitude") +
      guides(colour = guide_legend(override.aes = list(size = 4)))+
      theme(axis.title=element_text(size=tsize),
            panel.background = element_rect(fill = "aliceblue", colour = "grey1"),
            panel.border = element_rect(linetype = "solid", fill = NA, colour="black"),
            legend.title = element_text(size=tsize),
            legend.text = element_text(size=osize),
            legend.position = "bottom",
            panel.grid.major = element_line(colour = "grey85",linetype="dashed"),
            panel.grid.minor = element_line(colour = "grey90"),
            legend.key = element_rect(fill = "transparent", colour = "transparent"),
            legend.key.size = unit(.8, "cm"))
  }else{
    print(metric)
    #plot for correlation even if called skill
    map=ggplot(basemap) +
      geom_sf(fill="gray95", color=NA) +
      geom_sf(data=parpl,aes(geometry=geometry,fill=skill,size=UpA),color="transparent",alpha=.9,shape=21,stroke=0)+ 
      geom_sf(fill=NA, color="grey20") +
      geom_sf(data=parpl,aes(geometry=geometry,size=UpA),col="black",alpha=1,stroke=0.1,shape=1)+ 
      scale_x_continuous(breaks=seq(-30,40, by=5)) +
      scale_size(range = c(1, 3), trans="sqrt")+
      scale_fill_gradientn(
        colors=palet,n.breaks=5,oob = scales::squish, limits=limi, name=metric) +
      coord_sf(xlim = c(min(nco[,1]),max(nco[,1])), ylim = c(min(nco[,2]),max(nco[,2])))+
      labs(x="Longitude", y = "Latitude") +
      guides(fill = guide_colourbar(barwidth = 0.5, barheight = 12,reverse=F))+
      theme(axis.title=element_text(size=tsize),
            panel.background = element_rect(fill = "aliceblue", colour = "grey1"),
            panel.border = element_rect(linetype = "solid", fill = NA, colour="black"),
            legend.title = element_text(size=tsize),
            legend.text = element_text(size=osize),
            legend.position = "right",
            panel.grid.major = element_line(colour = "grey85",linetype="dashed"),
            panel.grid.minor = element_line(colour = "grey90"),
            legend.key = element_rect(fill = "transparent", colour = "transparent"),
            legend.key.size = unit(.8, "cm"))
  }
  return(list(points,map))
}
ScatterValid=function(obsin,nobs,simin,nsim,scplot,valid_path,ValidSta){
  
  # obs=read.csv(paste0(valid_path,fileobs))
  # sim=read.csv(paste0(valid_path,filesim))
  # names(obs)[1] = names(sim)[1]="Station_ID"
  # 
  # 
  # if (length(obs[1,])>2){
  #   years=c(1950:2020)
  #   obs=obs[-73]
  #   names(obs)[c(2:72)]=years
  #   sim=sim[-73]
  #   names(sim)[c(2:72)]=years
  # }
  # 
  # #Only keep valid stations
  # valid=match(ValidSta$V1, obs$Station_ID)
  # obsv=obs[valid,]
  # obstest<-melt(obsv,id.vars ="Station_ID")
  # obstest2=obstest[which(obstest$value>0),]
  # 
  # length(unique(obstest2$Station_ID))
  # 
  # valid=match(ValidSta$V1, sim$Station_ID)
  # simv=sim[valid,]
  # simtest<-melt(simv,id.vars ="Station_ID")
  # simtest2=simtest[which(simtest$value>0),]
  # 
  # length(unique(simtest2$Station_ID))
  # 
  
  #Output plots
  
  if (length(obsin[1,])>3){
    obs_sim=inner_join(obsin,simin, by=c("Station_ID","variable"))
    obs_sim$value.x[which(is.infinite(obs_sim$value.x))]=NA
    obs_sim$value.y[which(is.infinite(obs_sim$value.y))]=NA
    length(unique(obs_sim$Station_ID))
    names(obs_sim)[c(3,4)]=c("obs","sim")
    dat=as.data.frame(obs_sim[c(1,2,3,4)])
    
    #average
    dat_av=aggregate(list(obs=dat$obs,sim=dat$sim),
                     by = list(station=dat$Station_ID),
                     FUN = function(x) c(mean=mean(x,na.rm=T)))
    
    dat <- do.call(data.frame, dat_av)
    dat=dat[which(!is.na(dat$obs)),]
  }else{
    obsin=obsin[,c(1,3)]
    simin=simin[,c(1,3)]
    obs_sim=inner_join(obsin,simin, by=c("Station_ID"))
    obs_sim$value.x[which(is.infinite(obs_sim$value.x))]=NA
    obs_sim$value.y[which(is.infinite(obs_sim$value.y))]=NA
    length(unique(obs_sim$Station_ID))
    names(obs_sim)[c(2,3)]=c("obs","sim")
    dat=obs_sim
  }
  r=cor(dat$sim,dat$obs)
  r2=round(r^2,3)
  
  mmax=max(max(dat$obs,dat$sim))
  mmin=min(min(dat$obs,dat$sim))
  tmax=round(mmax/1000)*1000
  tmin=max(mmin,1e-3)
  
  bi=2000
  if (tmax<10000) bi=1000
  if (tmax<4000) bi=500
  #dat$year=as.numeric(as.character(dat$variable))
  # dat$obsj=jitter(dat$obs,0.001)+0.0001
  # dat$simj=jitter(dat$sim,0.001)+0.0001
  dat$density <- get_density(log(dat$obs), log(dat$sim), n = 100)
  #dat$density <- get_density(log(dat$obsj), log(dat$simj), n = 100)
  ggplot(dat,environment=environment()) + 
    # geom_point(aes(x=obs, y=sim,col=density),stroke=0,size=3,alpha=0.2,shape=16) +
    # geom_point(aes(x=obs, y=sim,col=density),stroke=0,size=3,alpha=0.2,shape=16) +
    geom_jitter(aes(x=obs, y=sim,col=density),stroke=0,size=3,alpha=0.5,shape=16, height=.1, width=.0)+
    geom_abline(slope=1, intercept=0, lwd =1, alpha=1,col="darkgreen",linetype="dashed")+
    scale_x_log10(name=nobs,
                  breaks=c(0.1,1,10,100,1000,10000), minor_breaks = log10_minor_break(),
                  labels=c("0.1","1","10","100","1000","10000"), limits=c(tmin,tmax)) +
    scale_y_log10(name=nsim,
                  breaks=c(0.1,1,10,100,1000,10000), minor_breaks = log10_minor_break(),
                  labels=c("0.1","1","10","100","1000","10000"), limits=c(tmin,tmax)) +
    # scale_x_continuous(name=nobs,breaks=seq(0,tmax, by=bi))+
    # scale_y_continuous(name=nsim,breaks=seq(0,tmax, by=bi))+
    annotate("label", x=3*tmin, y=max(tmax), label= paste0("R2 = ",r2),size=5)+
    scale_color_viridis(option="A")+
    theme(axis.title=element_text(size=16, face="bold"),
          axis.text = element_text(size=16),
          panel.background = element_rect(fill = "white", colour = "white"),
          panel.grid = element_blank(),
          panel.border = element_rect(linetype = "solid", fill = NA, colour="black"),
          legend.title = element_text(size=14),
          legend.text = element_text(size=12),
          legend.position = "none",
          plot.margin = margin(1,1,1,1, "cm"),
          panel.grid.major = element_line(colour = "grey80"),
          panel.grid.minor = element_line(colour = "grey90",linetype="dashed"),
          legend.key = element_rect(fill = "transparent", colour = "transparent"),
          legend.key.size = unit(.8, "cm"))
  
  #dat$rat=dat$sim/dat$obs
  #return(dat)
  #ggsave(paste0(scplot,".jpg"), width=20, height=15, units=c("cm"),dpi=1500)
}
StatEpisode <- function (data, date.deb = NA, win = 10, qx.min = NULL, mode = "FL",obs=TRUE,plot=FALSE){
  x <- data
  # if (class(x$date)[1] == "POSIXt") {
  # if (inherits(x$date,"POSIXt")){
  #   date.deb <- as.POSIXct(date.deb, tz = "utc")
  # }
  if (mode=="FL"){
    ind <- which(x$date >= date.deb-win/2)
    ip <- which(ind >= min(ind) & ind <= (min(ind) + win))
    ind <- ind[ip]
  }
  if (mode=="DR"){
    ind <- which(x$date >= date.deb-1)
    ip <- which(ind >= min(ind) & ind <= (min(ind) + win+1))
    ind <- ind[ip]
  }
  
  if (length(ind) < 3) {
    cat("Episodes du ", format(date.deb, "%d-%b-%y"), "ignore (n<3)\n")
    return(NULL)
  }
  if (obs==T) xep <- x[ind, c("date", "Q" , "Qs")]
  if (obs==F) {xep <- x[ind, c("date", "Qs")]
  xep$Q=xep$Qs}
  # if (!is.null(fp)) {
  #   tmp <- abat.niv(xep, fp)
  #   xep <- data.frame(xep, Peff = tmp$P)
  # }
  if (is.null(qx.min))
    qx.min <- max(x$Q, na.rm = TRUE)
  #     if (is.null(px.min))
  #         px.min <- max(x$P, na.rm = TRUE)
  if (mode=="FL"){
    q.max <- max(c(xep$Q, xep$Qs, qx.min), na.rm = TRUE)
    q.maxsim<-max(na.omit(xep$Qs))
    if (obs==T) {q.maxobs<-max(na.omit(xep$Q))
    ec.pic = round((q.maxobs-q.maxsim)/q.maxobs*100,1)}
    # p.max <- max(c(xep$P, px.min), na.rm = TRUE)
    q.range <- c(0, 1.5 * q.max)
    # p.range <- c(0, 3 * p.max)
    ec=ec.pic
  }
  if (mode=="DR"){
    q.max <- max(c(xep$Q, xep$Qs, qx.min), na.rm = TRUE)
    q.minsim<-min(na.omit(xep$Qs))
    if (obs==T){ q.minobs<-min(na.omit(xep$Q))
    ec.vol = round((sum(na.omit(xep$Qs)-na.omit(xep$Q)))/sum(na.omit(xep$Q))*100,1)}
    # p.max <- max(c(xep$P, px.min), na.rm = TRUE)
    q.range <- c(0, 1.2 * q.max)
    # p.range <- c(0, 3 * p.max)
    ec=ec.vol
  }
  if (obs==T) kge<-c(kge=KGE(xep$Q,xep$Qs),r=cor(xep$Q,xep$Qs))
  if (plot==T){
    par(mar = c(2, 3, 2, 3), mgp = c(1.7, 0.5, 0))
    # plot.window(range(xep$date))#, p.range)
    plot.new()
    
    #     plot.hyeto(xep$date, xep$P, curve = TRUE, ylim = rev(p.range),
    #         grid = FALSE, col = "lightblue", border = "azure", lwd = 2)
    #     abline(v = axis.POSIXct(1, xep$date), col = "lightgrey",
    #         lty = 3)
    #     if (!is.null(fp))
    #         plot.hyeto(xep$date, xep$Peff, curve = TRUE, ylim = rev(p.range),
    #             grid = FALSE, border = NA)
    #     axis(4, pretty(c(0, p.max), n = 5))
    #     mtext("PrÃ©cipitations (mm)", 4, 1.7, cex = par()$cex)
    
    par(new = TRUE)
    plot(xep$date, xep$Qs, ylim = q.range, type = "n", axes = FALSE,
         xlab = "", ylab = expression(paste("Discharge ", (m^3/s),
                                            sep = " ")))
    h = axis(2, pretty(c(0, q.max*1.5), n = 5))
    axis(1,xep$date, label=format(xep$date,"%d %b"),cex.axis=.8)
    box()
    points(xep$date, xep$Qs, type = "l", col = "red", pch = 2,
           lwd = 2, cex = 0.6)
    
    if (obs==T) {points(xep$date, xep$Q, type = "l", col = "royalblue", lwd = 2)
      legend("left", legend = c("Histo.", "Recal."), 
             col = c("royalblue","red"), pch = c(NA, NA), lwd = 2, bty = "n", cex = 0.7)
    }else{
      
    }
    if (mode =="FL"){
      #abline(h=quantile(data$Q,0.90,na.rm=T),col="darkgrey",lwd=2) 
      if (obs==T){
        txt.lab <- paste("KGE =", round(kge[1],2), "\nEc-pic =",ec.pic," %","\nr =",round(kge[2],3))
        text(mean(xep$date), 1.2*q.max-0.1*q.max, pos=4, lab=txt.lab,cex=0.7)
      }
      titre <- paste("Episode du ", format(xep$date[which.max(xep$Q)],
                                           "%a %d %B %Y"))
      ec=ec.pic
    }
    if (mode =="DR"){
      abline(h=quantile(data$Qs,0.05,na.rm=T),col="darkgrey",lwd=2) 
      if (obs==T){
        txt.lab <- paste("KGE =", round(kge[1],2), "\nEc-vol =",ec.vol," %","\nr =",round(kge[2],3))
        text(mean(xep$date), 1.2*q.max-0.1*q.max, pos=4, lab=txt.lab,cex=0.7)
      }
      titre <- paste("Episode du ", format(xep$date[1],"%a %d %B %Y") ,"au ",
                     format(xep$date[win], "%a %d %B %Y"))
      ec=ec.vol
    }
    mtext(titre,3,font = 2,adj = 0,line = 0.5,cex = .75)
    # legend("right", legend = c("Pluie", "Neige"), pch = 22, lwd = NA,
    # lty = NA, col = "lightblue", pt.bg = c("lightblue", "azure"),
    # pt.lwd = 1, pt.cex = 2, bty = "n", cex = 0.7)
    # legend("right", legend = "Pluie", pch = 22, lwd = NA,
    #         lty = NA, col = "azure", pt.bg = "lightblue",
    #         pt.lwd = 1, pt.cex = 2, bty = "n", cex = 0.7)
  }
  if (obs==T) return(c(kge,error=ec,QmObs=q.maxobs,QmObs=q.maxobs,year=year(xep$date)[1]))
}
computeAnnualMaxima <- function(timeAndSeries) {
  timeStamps <- timeAndSeries[,1]
  srs <- timeAndSeries[,2]
  
  tmvec <- as.Date(timeStamps)
  years <- format(tmvec, "%Y")
  
  findMax <- function(subIndxs) {
    subIndxMaxIndx <- which.max(srs[subIndxs])[1]
    subIndxs[subIndxMaxIndx]
  }
  annualMaxIndx <- tapply(1:length(srs), years, findMax)
  annualMaxIndx <- annualMaxIndx[!is.na(annualMaxIndx)]
  annualMax <- srs[annualMaxIndx]
  annualMaxDate <- timeStamps[annualMaxIndx]
  
  return(list(annualMax = annualMax, annualMaxDate = annualMaxDate, annualMaxIndx = annualMaxIndx))
}
computeAnnualMinima <- function(timeAndSeries) {
  timeStamps <- timeAndSeries[,1]
  srs <- timeAndSeries[,2]
  
  tmvec <- as.Date(timeStamps)
  years <- format(tmvec, "%Y")
  
  findMin <- function(subIndxs) {
    subIndxMinIndx <- which.min(srs[subIndxs])[1]
    subIndxs[subIndxMinIndx]
  }
  annualMinIndx <- tapply(1:length(srs), years, findMin)
  annualMinIndx <- annualMinIndx[!is.na(annualMinIndx)]
  annualMin <- srs[annualMinIndx]
  annualMinDate <- timeStamps[annualMinIndx]
  
  return(list(annualMin = annualMin, annualMinDate = annualMinDate, annualMinIndx = annualMinIndx))
}
#Circular statistics for mean date
seasony=function(x){
  theta=x*(2*pi/365.25)
  #plot(cos(theta),sin(theta),xlim=c(-1,1),ylim=c(-1,1))
  xi=1/(length(theta))*sum(cos(na.omit(theta)))
  yi=1/(length(theta))*sum(sin(na.omit(theta)))
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
season1=function(x){
  l1=length(which(!is.na(x)))
  if(l1>0){
    x=x[which(!is.na(x))]
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
  }else{Di=NA}
  return(Di)
}
season2=function(x){
  l1=length(which(!is.na(x)))
  if(l1>0){
    x=x[which(!is.na(x))]
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
  }else{Di=NA
  R=NA}
  return(R)
}
resOpen=function(dir,outletname,ValidSY){
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
  lonlatloop=expand.grid(c(1:llo),c(lla:1))
  outll$res.ratio=outlets
  outll$idlo=lonlatloop$Var1
  outll$idla=lonlatloop$Var2
  
  outll$idlalo=paste(outll$idlo,outll$idla,sep=" ")
  outlp=outll[which(!is.na(outll$res.ratio)),]
  mf=match(ValidSY$idlalo,outlp$idlalo)
  ValidSY$res.ratio=outlp$res.ratio[mf]
  return (ValidSY)
}

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
  
  mierda=which(!is.na(outlets))
  outll=outll[mierda,]
  outlets=outlets[which(!is.na(outlets))]
  outll=data.frame(outlets,outll)
  return (outll)
}
# Files paths -------------------------------------------------------------
main_path = 'D:/tilloal/Documents/06_Floodrivers/'
valid_path = paste0(main_path,'DataPaper/')
dis_path<-paste0(main_path,'dis/calibrated/filtered/Histo/')
setwd(valid_path)

# Part 1: Create the Valid station file -------------------------------------------
min_daily_value <- function(timeseries) {
  # Convert the 6-hourly time series to daily time series
  daily_timeseries <- apply.daily(timeseries, min)
  tday=unique(as.Date(timeseries[,1]))
  return(data.frame(date=tday,Qmd=daily_timeseries))
}
timeseries=frostserie
#Loading all txt files with matching locations
#load matching coordinates
Sloc=read.table(paste0(valid_path,"out/Stations_locations_EFAS_grid_1950.txt"),sep=",")
yrlist=c(1951:2020)
for (yi in yrlist){
  print(yi)
  Slocy=read.table(paste0(valid_path,"out/Stations_locations_EFAS_grid_",yi,".txt"),sep=",")
  Slocrep=Slocy[which(Sloc$V2==0),]
  Sloc[which(Sloc$V2==0),]=Slocrep
}
Sloc_final=Sloc[which(Sloc$V2!=0),]
Sloc_final$csource="SpatialQMatch"

#transalting to R indexing after python
Sloc_final$V2=Sloc_final$V2+1
Sloc_final$V3=Sloc_final$V3+1

#load file of stations from EFAS
flefas=read.csv(paste0(valid_path,"Stations/efas_flooddriver_match.csv"))
efasmatch=match(flefas$StationID,Sloc_final$V1)
faichier=inner_join(flefas,Sloc_final,by=c("StationID"="V1"))
Sloc_final$csource[efasmatch]="EFAS"
Sloc_final$idlalo=paste(Sloc_final$V3,Sloc_final$V2, sep=" ")

#Now the upstream area
outletname="GIS/upArea_European_01min.nc"
dir=valid_path
UpArea=UpAopen(valid_path,outletname,Sloc_final)
head(UpArea)
ValidSta=UpArea[which(UpArea$upa>100),]


## Flag stations that have very different mean discharges -------------
fileobs="out/obs_meanAY.csv"
filesim="out/EFAS_meanAY.csv"
obs=read.csv(paste0(valid_path,fileobs))
names(obs)=c("Station_ID","mean")
sim=read.csv(paste0(valid_path,filesim))
names(sim)=c("Station_ID","mean")

years=c(1950:2020)
obstest=obs[which(!is.infinite(obs$mean)),]
simtest=sim[which(!is.infinite(sim$mean)),]
obs_sim=inner_join(obstest,simtest, by=c("Station_ID"))
names(obs_sim)[c(2,3)]=c("obs","sim")

#Remove stations that were removed in the first step
rmv0=which(!is.na(match(obs_sim$Station_ID,ValidSta$V1)))
obs_sim=obs_sim[rmv0,]
rmv1=match(obs_sim$Station_ID,ValidSta$V1)
obs_sim$UpA=ValidSta$upa[rmv1]
## Flag EFAS stations that were used for calibration --------
efas_stations_orig=read.csv("Stations/stations_efas_meta.csv",sep=";")
efas_stations_orig$latlong=paste(efas_stations_orig$LisfloodX, efas_stations_orig$LisfloodY, sep=" ")
efas_stations=flefas
efas_stations$latlong=paste(efas_stations$LisfloodX, efas_stations$LisfloodY, sep=" ")
efas_stations2=inner_join(efas_stations, efas_stations_orig,by="latlong")
efas_stations_cal=efas_stations2[which(efas_stations2$EC_calib==1),]
ValidSta$calib=FALSE
matchcal=match(efas_stations_cal$StationID,ValidSta$V1)
ValidSta$calib[matchcal]=TRUE

## Number of years with data ------------
Q_data <- read.csv(paste0('Q_19502020.csv'), header = F)  # CSVs with observations
HERA_data<-read.csv(file="out/EFAS_19502020.csv")
Q_data=Q_data[-1,-1]
Station_data_IDs <- as.vector(t(Q_data[1, -1]))
lR=c()
for (id in 1:length(Station_data_IDs)){
  Qs=Q_data[-1,id+1]
  xl=length(Qs)
  l=length(which(!is.na(Qs)))
  lR=c(lR,l)
}
plot(lR[order(lR)])
RecordLen=data.frame(Station_data_IDs,lR)
rmv2=match(obs_sim$Station_ID,RecordLen[,1])
obs_sim$Rlen=RecordLen$lR[rmv2]
dat=obs_sim

# Part 2: TSEVA on observed values------------------------------------
# TSEVA on observed values
ValidSf=read.csv(file="Stations/Stations_Validation.csv")[,-1]
ValidSY=ValidSf[which(ValidSf$removal!="YES"),]
length(ValidSY$comment[which(ValidSY$UpA<=250)])/2901
length(ValidSY$comment[which(ValidSY$skill>-0.41)])
median(ValidSY$skill)

vEFAS=ValidSY[which(ValidSY$csource=="EFAS"),]
median(vEFAS$skill)

######################################################################################
outletname = "outletsv8_hybas07_01min"
hydroDir<-("D:/tilloal/Documents/LFRuns_utils/data")

#load frost days
# load(file=paste0(hydroDir,"/Drought/catchment_frost.Rdata"))
# saveRDS(frostcat,paste0(hydroDir,"/Drought/catchment_frost.RDS"))

frostcat=readRDS(paste0(hydroDir,"/Drought/catchment_frost.RDS"))
#Hybas07
Catchmentrivers7=read.csv(paste0(hydroDir,"/Catchments/from_hybas_eu_onlyid.csv"),encoding = "UTF-8", header = T, stringsAsFactors = F)

#location of the outlets

outletname = "outletsv8_hybas07_01min"
outhyb07=outletopen(hydroDir,outletname)
outhyb07$latlong=paste(round(outhyb07$Var2,4),round(outhyb07$Var1,4),sep=" ")
head(outhyb07)
UpArea$latlong=paste(round(UpArea$Var2,4),round(UpArea$Var1,4),sep=" ")
Catchmentrivers7$latlong=paste(round(Catchmentrivers7$POINT_Y,4),round(Catchmentrivers7$POINT_X,4),sep=" ")
head(Catchmentrivers7)
CRplus=inner_join(outhyb07,Catchmentrivers7,by="latlong")
outboost=inner_join(UpArea,Catchmentrivers7,by="latlong")
head(CRplus)
head(outboost)

mycat=inner_join(outboost,CRplus, by="HYBAS_ID")
head(mycat)


hybas07 <- read_sf(dsn = paste0(hydroDir,"/Catchments/hydrosheds/hybas_eu_lev07_v1c.shp"))
hybasf7=fortify(hybas07) 
Catamere07=inner_join(hybasf7,outboost,by= "HYBAS_ID")
Catamere07$llcoord=paste(round(Catamere07$POINT_X,4),round(Catamere07$POINT_Y,4),sep=" ") 
Catf7=Catamere07
head(Catf7)
st_geometry(Catf7)=NULL
min(Catf7$pointid)
#Extract only catchments where there is at least 40 years of observation
#70 years in days
yind<-50*365
ValidSLong=ValidSY[which(ValidSY$Rlen>=yind),]
length(ValidSLong$calib)

#now extract the HERA discharge for catchments in ValidSLong
station_HERA<-HERA_data[1,-1]
SelectHera=which(station_HERA %in% ValidSLong$V1)
station_ids=data.frame(id=t(station_HERA[SelectHera]))
HERA_dates<-HERA_data[,1]
HERA_data2<-HERA_data[,-1]
HERA_comp<-HERA_data2[-1,SelectHera]

station_Q<-Q_data[1,-1]
SelectQ=which(station_Q %in% ValidSLong$V1)
Q_dates<-Q_data[,1]
Q_data2<-Q_data[,-1]
Q_comp<-Q_data2[-1,SelectQ]


station_idplus=left_join(station_ids,mycat,by=c("X1"="V1"))

#Fit TSEVA on observed values
dates_Hera <- seq(as.Date("1950-01-04"), as.Date("2020-12-31"), by = "days")
dates_Q <- seq(as.Date("1950-01-01"), as.Date("2020-12-31"), by = "days")
#load TSEVA functions
source("~/LFRuns_utils/TSEVA_demo/demo_functions.R")
library(xts)


#Extract 1 timeseries out of HERA_comp
#HERA_comp is a matrix with the first row being the dates
#and the other rows the discharge values for each station
haz="flood"
season="nonfrost"
dset="obs"
#start the loop
RetPerGPD=c()
RetPerGEV=c()
RetLevGEV=c()
RetLevGPD=c()
parlist=c()
peaklist=c()
catlist=c()
IRES=c()
Reservoir_i=c()

#loop over all stations
for(id in seq_along(station_ids$X1)){
  # id=23
  start_time <- Sys.time()
  Reservoir_alteration=0
  station=station_idplus$X1[id]
  print(paste0(id,"/",length(station_ids$X1)))

  #the hybas catchment for frost days
  catch=station_idplus$pointid.y[id]
  #Select the first station
  if (dset=="obs"){
    Q=Q_comp[,id]
    txx=dates_Q
  }else if (dset=="sim"){
    Q=HERA_comp[,id]
    txx=dates_Hera
  }
  Q=as.numeric(Q)
  df.dis=Q

  rmv=which(year(txx)==1950)
  df.dis=data.frame(txx,df.dis)
  names(df.dis)[c(1,2)]=c("time","dis")
  #remove 1950 which is not reliable
  df.dis=df.dis[-rmv,]
  timeStamps=txx
  dt = difftime(timeStamps[2],timeStamps[1],units="days")
  dt= as.numeric(dt)

  if (haz=="drought"){
    #seasonality divide: frost vs non frost
    percentile=95
    Tcatchment=which(colnames(frostcat)==catch)
    data=df.dis
    names(data)=c("date","Qs")
    intermit=interid(data, WindowSize=7)
    interflag=intermit$flags[2]
    if (!exists("trans")){trans="rev"}
    print(paste0(trans," transformation used for low flows"))
    series=data.frame(txx[-rmv],intermit$trdis$Q7)
    #remove frost timesteps, this can be modified to do the anlysis only on frost moments
    if (length(Tcatchment)>0){
      frostserie=data.frame(frostcat[,1],frostcat[,Tcatchment])
      names(frostserie)=c("time","Ta")
      frostserie=min_daily_value(frostserie)
      rmv2=which(year(frostserie$date)==1950)
      frostserie=frostserie[-rmv2,]
      frostserie=frostserie[-rmv,]
      frosttime=which(frostserie[,2]<0)
    }else{
      frosttime=NA
    }
  }else if (haz=="flood"){
    data=df.dis
    catmat=Catf7[which(Catf7$pointid==catch),]
    names(data)=c("date","Qs")
    # Extract the maximum daily value
    series <- max_daily_value(data)
    # series=data.frame(txx,df.disX$outlets)
    interflag=0
    frosttime=NA
    trans="ori"
  }

  names(series)=c("timestamp","dis")
  dt1=min(diff(series$timestamp),na.rm=T)
  dt=as.numeric(dt1)
  tdim=attributes(dt1)$units
  if (tdim=="hours") dt=dt/24
  if (dt==1){
    timeDays=series$timestamp
  }else{
    timeDays=unique(as.Date(series$timestamp))
  }

  bounds=c(year(timeDays[1]),year(timeDays[length(timeDays)]))
  realbound=bounds
  tbound=c(as.Date(paste0(realbound[1],"-12-31")),as.Date(paste0(realbound[2],"-12-31")))
  Impdates=seq(tbound[1],tbound[2],by="1 year")
  tindexes=match(Impdates,timeDays)
  names(series)=c("timestamp","dis")
  nv=length(unique(series$dis))
  timeAndSeries=series
  names(timeAndSeries)=c("timestamp","data")

  if(length(na.omit(series$dis))>1 & interflag<1 & nv>15){

    if (haz=="drought" && length(!is.na(frosttime))>1) {
      if (season=="nonfrost") {
        print("nonfrost season")
        timeAndSeries$data[frosttime]=NA
      } else if (season=="frost") {
        print("frost season")
        timeAndSeries$data[-frosttime]=NA
      } else {
        print("season must be frost or nonfrost")
      }
    }

    tsm=1/dt
    series=timeAndSeries[,2]
    timeWindow = 365.25*30; #time windows in days, the correction is done within the functions
    windowSize=366
    minPeakDistanceInDays=30
    timeStamps=timeAndSeries$timestamp
    # Choose transformation
    transftypes=c("ori","rev","inv","lninv")
    trendtypes=c("trend","trendPeaks")
    trendtrans=expand.grid(transftypes,trendtypes)
    #choose transformation
    if (haz=="flood") tt=5
    if (haz=="drought") tt=6

    # Apply TsEvaNs function
    Nonstat <- TsEvaNs(timeAndSeries, timeWindow, transfType=trendtrans[tt,2], 
    ciPercentile= 90, minPeakDistanceInDays = minPeakDistanceInDays, 
    mode=haz,TrendTh=0.85, trans=trendtrans[tt,1])

    nonStationaryEvaParams=Nonstat[[1]]
    stationaryTransformData=Nonstat[[2]]

    # Define range for plotting
    ExRange= c(min(nonStationaryEvaParams$potObj$parameters$peaks),max(nonStationaryEvaParams$potObj$parameters$peaks))
    if (haz=="flood") wr2 <- c(seq(min(ExRange),max(ExRange),length.out=700))
    if (haz=="drought") wr2 <- c(seq(1.1*min(ExRange),0.1*max(ExRange),length.out=700))

    # Apply tsEvaPlotGPDImageScFromAnalysisObj function
    #Plot1= tsEvaPlotGPDImageScFromAnalysisObj(wr2, nonStationaryEvaParams, stationaryTransformData, minYear = '1950',trans="rev")  
    timeIndex=1
    Plot2 = tsEvaPlotReturnLevelsGPDFromAnalysisObj(nonStationaryEvaParams, stationaryTransformData, timeIndex, trans=trendtrans[tt,1],ylabel="Discharge (m3/s)",ope=T)
    #print(Plot2)
    ggsave(paste0(hydroDir,"/TSEVA_hybas/valid/GPDBeam_",haz,"_",dset,"_Cat_",station,".jpg"),Plot2, width=24, height=20, units=c("cm"),dpi=300)
    #Plot3 = tsEvaPlotReturnLevelsGEVFromAnalysisObj(nonStationaryEvaParams, stationaryTransformData, timeIndex, trans="rev",ylabel="Discharge (m3/s)")
    
    stationaryTransformData$timeStampsDay=unique(as.Date(stationaryTransformData$timeStamps))
    pikos=data.frame(nonStationaryEvaParams$potObj$parameters$peaks,nonStationaryEvaParams$potObj$parameters$peakID,nonStationaryEvaParams$potObj$parameters$peakST, nonStationaryEvaParams$potObj$parameters$peakEN)
    names(pikos)=c("value","timeID","tIDstart","tIDend")
    pikos$time=timeStamps[pikos$timeID]
    pikos$catch=rep(catch,length(pikos[,1]))
    #Here I need to convert the timeStamp to a daily one if dt is not 1
    dt1=min(diff(timeStamps),na.rm=T)
    dt=as.numeric(dt1)
    tdim=attributes(dt1)$units
    if (tdim=="hours") dt=dt/24
    if (dt==1){
      timeDays=stationaryTransformData$timeStamps
    }else{
      timeDays=stationaryTransformData$timeStampsDay
    }

    bounds=c(year(timeDays[1]),year(timeDays[length(timeDays)]))
    realbound=bounds
    tbound=c(as.Date(paste0(realbound[1],"-12-31")),as.Date(paste0(realbound[2],"-12-31")))
    Impdates=seq(tbound[1],tbound[2],by="1 years")
    datex=yday(timeDays)
    dtect=c(diff(datex),-1)
    last_days <- timeDays[which(dtect<0)]
    tindexes=match(last_days,timeDays)

    # Compute return periods and levels
    RPgoal=10
    timeIndex=tindexes[1]
    RLevs100=ComputeReturnLevels(nonStationaryEvaParams, RPgoal, timeIndex,trans=trendtrans[tt,1])

    if (RLevs100$Fit=="No fit"){

      RLgpd=nonStationaryEvaParams$gevObj$parameters$annualMax
      names(RLgpd)=year(Impdates)

      RLgev=nonStationaryEvaParams$gevObj$parameters$annualMax
      names(RLgev)=year(Impdates)

      nRPgev=rep(NA, length(Impdates))
      names(nRPgev)=year(Impdates)

      nRPgpd=rep(NA, length(Impdates))
      names(nRPgpd)=year(Impdates)

      params=data.frame(matrix(ncol=17,nrow=(length(Impdates)-1)))
      params[,1]=rep(catch,7)
      params[,2]=year(Impdates)[-1]
      params[,4]=rep(interflag,7)
      if (is.null(colnames(parlist))){
        colnames(params)=rep("nom",17)
      }else{
        colnames(params)=colnames(parlist)
      }
    }else{
      #####
      RLgev=RLevs100$ReturnLevels[2]
      RLgpd=RLevs100$ReturnLevels[3]
      ERgev=RLevs100$ReturnLevels[4]
      ERgpd=RLevs100$ReturnLevels[5]
      nRPgev=nRPgpd=10
      params=c()
      for (t in 2:length(Impdates)){
        timeIndex=tindexes[t]
        RLevs100i=ComputeReturnLevels(nonStationaryEvaParams, RPgoal, timeIndex,trans=trendtrans[tt,1])
        params=c(catch,year(Impdates[t]),timeIndex,RLevs100i$Params)
        names(params)[1:3]=c("catchment","Year","timeIndex")

        Rper=RPcalc(params,RPiGEV=RLevs100$ReturnLevels[2],RPiGPD=RLevs100$ReturnLevels[3])
        nRPgpd=c(nRPgpd,Rper[2])
        nRPgev=c(nRPgev,Rper[1])
        RLgev=cbind(RLgev,RLevs100i$ReturnLevels[2])
        RLgpd=cbind(RLgpd,RLevs100i$ReturnLevels[3])
        ERgev=cbind(ERgev,RLevs100i$ReturnLevels[4])
        ERgpd=cbind(ERgpd,RLevs100i$ReturnLevels[5])
        if (length(parlist)>1) colnames(parlist)=names(params)
        parlist=rbind(parlist,params)
      }

      RLgev=as.data.frame(RLgev)
      names(RLgev)=year(Impdates)
      rownames(RLgev)=RPgoal

      RLgpd=as.data.frame(RLgpd)
      names(RLgpd)=year(Impdates)
      rownames(RLgpd)=RPgoal

      nRPgev=as.data.frame(t(nRPgev))
      names(nRPgev)=year(Impdates)

      nRPgpd=as.data.frame(t(nRPgpd))
      names(nRPgpd)=year(Impdates)

      peaklist=rbind(peaklist,pikos)
    }
  }else{
    cat(paste0("\n No values in this pixel \n or intermittent river (flag = ",interflag,")"))
    if (is.na(interflag)) interflag=-9999
    if (interflag>0){
      filling=intermit$flags[3]
      datex=yday(intermit$zeroRP$time)
      dtect=c(diff(datex),-1)
      last_days <- intermit$zeroRP$time[which(dtect<0)]
      tindexes=match(last_days,intermit$zeroRP$time)
      oops=intermit$zeroRP[tindexes,]
      # verif=series$timestamp[which(!is.na(match(year(series$timestamp),year(Impdates))))]
      # dfrep=data.frame(year(verif),oops)
      # yagg=aggregate(list(m=dfrep$oops),
      #                by = list(year=dfrep$year.verif.),
      #                FUN = function(x) c(max=max(x)))

    }else{
      filling=NA
    }
    RLgpd=oops$RP
    names(RLgpd)=year(Impdates)

    RLgev=oops$RP
    names(RLgev)=year(Impdates)

    nRPgev=rep(NA, length(Impdates))
    names(nRPgev)=year(Impdates)

    nRPgpd=rep(NA, length(Impdates))
    names(nRPgpd)=year(Impdates)

    params=data.frame(matrix(ncol=17,nrow=(length(Impdates)-1)))
    params[,1]=rep(catch,length(Impdates)-1)
    params[,2]=year(Impdates)[-1]
    params[,4]=rep(interflag,length(Impdates)-1)
    if (is.null(colnames(parlist))){
      colnames(params)=rep("nom",17)
    }else{
      colnames(params)=colnames(parlist)
    }
    parlist=as.data.frame((rbind(parlist,params)))

    pikos=data.frame(matrix(ncol=6,nrow=1))
    pikos[,1]=NA
    pikos[,6]=catch
    if (is.null(colnames(peaklist))){
      colnames(pikos)=c("value","timeID","tIDstart","tIDend","time","catch")
    }else{
      colnames(pikos)=colnames(peaklist)
    }
    #parlist=as.data.frame(mapply(c,parlist,params))
    peaklist=as.data.frame((rbind(peaklist,pikos)))

  }

  #Saving main outputs
  catlist=c(catlist,catch)
  IRES=c(IRES,interflag)
  Reservoir_i=c(Reservoir_i,Reservoir_alteration)
  RetLevGEV=rbind(RetLevGEV,RLgev)
  RetLevGPD=rbind(RetLevGPD,RLgpd)
  RetPerGEV=rbind(RetPerGEV,nRPgev)
  RetPerGPD=rbind(RetPerGPD,nRPgpd)
  end_time <- Sys.time()
  cat(paste0("\nloop duration: ",round(end_time-start_time,2)," seconds\n"))
}


#saving the outputs
parlist=as.data.frame(parlist)
Results=list(parameters=parlist,RetLevGEV=RetLevGEV,RetLevGPD=RetLevGPD,RetPerGEV=RetPerGEV,RetPerGPD=RetPerGPD,Peaks=pikos,catrest=data.frame(catlist,Reservoir_i,IRES))


save(parlist,file=paste0(hydroDir,"/TSEVA_hybas/outputs/parValid_",haz,"_",dset,"_1950_2020.Rdata"))
save(Results, file=paste0(hydroDir,"/TSEVA_hybas/outputs/ResValid_",haz,"_",dset,"_1950_2020.Rdata"))


#Now I have to compare the trends

#load the simulated discharge

#compute the changes, plot the changes, save the changes, compare






cord.dec=ValidSY[,c(1,2)]
cord.dec = SpatialPoints(cord.dec, proj4string=CRS("+proj=longlat"))
cord.UTM <- spTransform(cord.dec, CRS("+init=epsg:3035"))
nco=cord.UTM@coords
world <- ne_countries(scale = "medium", returnclass = "sf")
Europe <- world[which(world$continent == "Europe"),]
e2=st_transform(Europe,  crs=3035)
w2=st_transform(world,  crs=3035)
tsize=12
osize=12
basemap=w2
palet=c(hcl.colors(9, palette = "viridis", alpha = NULL, rev = T, fixup = TRUE))
#palet = c('#f61f0f','#f4cccc','#f6ebeb','#ffffff','#e4edf5','#91bfdb','#1c87d8'))

ppl <- st_as_sf(ValidSY, coords = c("Var1", "Var2"), crs = 4326)
ppl <- st_transform(ppl, crs = 3035)

ggplot(basemap) +
  geom_sf(fill="gray95", color=NA) +
  geom_sf(data=ppl,aes(geometry=geometry,fill=Rlen,size=UpA),color="transparent",alpha=.9,shape=21,stroke=0)+ 
  geom_sf(fill=NA, color="grey20") +
  #geom_sf(data=ppl,aes(geometry=geometry,size=UpA),col="grey",alpha=1,stroke=0.1,shape=1)+ 
  scale_x_continuous(breaks=seq(-30,40, by=5)) +
  scale_size(range = c(1, 4), trans="sqrt",name= expression(paste("Upstream area ", (km^2),
                                                                 sep = " ")),
             breaks=c(101,1000,10000,100000,500000), labels=c("100","1000", "10 000", "100 000", "500 000"))+
  scale_fill_gradientn(
    colors=palet,oob = scales::squish, name="record length (years)", trans="sqrt",
    breaks=c(365,1825,3650,7300,14600,21900), labels=c(1,5,10,20,40,60)) +
  coord_sf(xlim = c(min(nco[,1]),max(nco[,1])), ylim = c(min(nco[,2]),max(nco[,2])))+
  labs(x="Longitude", y = "Latitude") +
  guides(fill = guide_colourbar(barwidth = 20, barheight = 0.5,reverse=F),
         size= guide_legend(override.aes = list(fill = "grey50")))+
  theme(axis.title=element_text(size=tsize),
        panel.background = element_rect(fill = "aliceblue", colour = "grey1"),
        panel.border = element_rect(linetype = "solid", fill = NA, colour="black"),
        legend.title = element_text(size=tsize),
        legend.text = element_text(size=osize),
        legend.position = "bottom",
        legend.box = "vertical",
        panel.grid.major = element_line(colour = "grey85",linetype="dashed"),
        panel.grid.minor = element_line(colour = "grey90"),
        legend.key = element_rect(fill = "transparent", colour = "transparent"),
        legend.key.size = unit(.8, "cm"))
ggsave("Plots/ValidStations.jpg", width=20, height=24, units=c("cm"),dpi=1500)



rat1=abs((dat$sim-dat$obs))/dat$obs
rat2=abs((dat$obs-dat$sim))/dat$sim
rmv=which(rat1>3 | rat2>3)
dat$rat1=rat1
dat$rat2=rat2
dat$flag1=NA
dat$flag1[which(dat$rat1>3 | dat$rat2>3)]=1
dat$flag1[which(dat$obs>10 & dat$flag1==1)]=2
dat$flag1[which(dat$rat1>6 | dat$rat2>6)]=3

#loading Dominiks initial database
Q_stations <- st_read(paste0(valid_path,'Europe_daily_combined_2022_v2_WGS84_NUTS.shp')) 

#Writing result: merging dat and validSta
ValidS=inner_join(ValidSta,dat,by=c("V1"="Station_ID"))
ValidS=inner_join(ValidS,Q_stations,by = c("V1"="StationID"))

#reorganise the data
ValidSf=ValidS[,c(1,2,7,19,4,5,6,14,10,11,12,13,15,18,22)]

## Identify stations with low KGE ---------------
kgefile="out/EFAS_obs_kgeAY.csv"
kge=SpatialSkillPlot(ValidSf,"kge",kgefile)
ValidSf=kge[[1]][,-16]
#Isolate problematic stations stations:
ValidSf$flag2=NA
ValidSf$flag2[which(ValidSf$skill<(-0.0))]=1
ValidSf$removal = ""
ValidSf$comment = ""
#write.csv(rest,file="stations_manual_check.csv")

## Individual check of the remaining stations ------------------
#load remaining stations checked
Mancheck=read.csv(file="Stations/stations_manual_check.csv")
Mancheck=Mancheck[,-1]
manmatch=match(Mancheck$V1,ValidSf$V1)
Vcheck=Mancheck[c(25,26)]
colnames(Vcheck)
colnames(ValidSf)[c(19,20)]
ValidSf[manmatch,c(19,20)]=Vcheck

#Other manual Checks
#Gagnieres 6139070
ValidSf$comment[which(ValidSf$V1==6139070)]="wrong river"
ValidSf$removal[which(ValidSf$V1==6139070)]="YES"
#6233203
ValidSf$comment[which(ValidSf$V1==6233203)]="dubious observations"
ValidSf$removal[which(ValidSf$V1==6233203)]="YES"
#6118175
ValidSf$comment[which(ValidSf$V1==6118175)]="wrong location"
ValidSf$removal[which(ValidSf$V1==6118175)]="YES"
#6124440
ValidSf$comment[which(ValidSf$V1==6124440)]="wrong location"
ValidSf$removal[which(ValidSf$V1==6124440)]="YES"
#6340215
ValidSf$comment[which(ValidSf$V1==6340215)]="wrong location"
ValidSf$removal[which(ValidSf$V1==6340215)]="YES"

ValidSf$comment[which(ValidSf$flag1==2 & ValidSf$removal=="")]="mismatch between obs and sim Qmean"
ValidSf$removal[which(ValidSf$flag1==2 & ValidSf$removal=="")]="YES"

ValidSf$comment[which(ValidSf$flag1==3 & ValidSf$removal=="")]="mismatch between obs and sim Qmean"
ValidSf$removal[which(ValidSf$flag1==3 & ValidSf$removal=="")]="YES"

## Distance between "official gauges" and efas points -----------------------
points=ValidSf[,c(1,2,3)]
Vsfloc=st_as_sf(points, coords = c("Var1", "Var2"), crs = 4326)
Statloc=Q_stations
Sfloc=inner_join(points,Statloc,by=c("V1"="StationID"))
dist=c()
for (r in Vsfloc$V1){
  cat(paste0(r,"\n"))
  v1=Vsfloc[which(Vsfloc$V1==r),]
  v2=Statloc[which(Statloc$StationID==r),]
  oula=st_distance(v1,v2)
  oula=oula/1000
  dist=c(dist,oula)
}

ValidSf$distance=dist
ValidSf$flag2[which(ValidSf$distance>2.5 & ValidSf$flag2==1 & ValidSf$removal=="")]=2
ValidSf$removal[which(ValidSf$flag2==2)]="YES"
ValidSf$comment[which(ValidSf$flag2==2)]="too far from original station"

length(ValidSf$comment[which(ValidSf$comment=="too far from original station")])
length(ValidSf$comment[which(ValidSf$comment=="mismatch between obs and sim Qmean")])
konar=(ValidSf[which(ValidSf$removal=="YES"),])
length(which(ValidSf$removal=="YES"))
length(ValidSf$comment[which(ValidSf$csource=="EFAS")])

#Writing of the final file
#write.csv(ValidSf,file="Stations/Stations_Validation.csv")
StatCheck=ValidSf[which(ValidSf$flag2==1),]
#write.csv(StatCheck,file="Stations/Stations_lowKGE.csv")


# Part 2: Diagnostic plots------------------------------------

ValidSf=read.csv(file="Stations/Stations_Validation.csv")[,-1]
ValidSY=ValidSf[which(ValidSf$removal!="YES"),]
length(ValidSY$comment[which(ValidSY$UpA<=250)])/2901
length(ValidSY$comment[which(ValidSY$skill>-0.41)])
median(ValidSY$skill)

vEFAS=ValidSY[which(ValidSY$csource=="EFAS"),]
median(vEFAS$skill)
#Map of upstream area
#Plot parameters

cord.dec=ValidSY[,c(1,2)]
cord.dec = SpatialPoints(cord.dec, proj4string=CRS("+proj=longlat"))
cord.UTM <- spTransform(cord.dec, CRS("+init=epsg:3035"))
nco=cord.UTM@coords
world <- ne_countries(scale = "medium", returnclass = "sf")
Europe <- world[which(world$continent == "Europe"),]
e2=st_transform(Europe,  crs=3035)
w2=st_transform(world,  crs=3035)
tsize=12
osize=12
basemap=w2
palet=c(hcl.colors(9, palette = "viridis", alpha = NULL, rev = T, fixup = TRUE))
#palet = c('#f61f0f','#f4cccc','#f6ebeb','#ffffff','#e4edf5','#91bfdb','#1c87d8'))

ppl <- st_as_sf(ValidSY, coords = c("Var1", "Var2"), crs = 4326)
ppl <- st_transform(ppl, crs = 3035)

ggplot(basemap) +
  geom_sf(fill="gray95", color=NA) +
  geom_sf(data=ppl,aes(geometry=geometry,fill=Rlen,size=UpA),color="transparent",alpha=.9,shape=21,stroke=0)+ 
  geom_sf(fill=NA, color="grey20") +
  #geom_sf(data=ppl,aes(geometry=geometry,size=UpA),col="grey",alpha=1,stroke=0.1,shape=1)+ 
  scale_x_continuous(breaks=seq(-30,40, by=5)) +
  scale_size(range = c(1, 4), trans="sqrt",name= expression(paste("Upstream area ", (km^2),
                                                                 sep = " ")),
             breaks=c(101,1000,10000,100000,500000), labels=c("100","1000", "10 000", "100 000", "500 000"))+
  scale_fill_gradientn(
    colors=palet,oob = scales::squish, name="record length (years)", trans="sqrt",
    breaks=c(365,1825,3650,7300,14600,21900), labels=c(1,5,10,20,40,60)) +
  coord_sf(xlim = c(min(nco[,1]),max(nco[,1])), ylim = c(min(nco[,2]),max(nco[,2])))+
  labs(x="Longitude", y = "Latitude") +
  guides(fill = guide_colourbar(barwidth = 20, barheight = 0.5,reverse=F),
         size= guide_legend(override.aes = list(fill = "grey50")))+
  theme(axis.title=element_text(size=tsize),
        panel.background = element_rect(fill = "aliceblue", colour = "grey1"),
        panel.border = element_rect(linetype = "solid", fill = NA, colour="black"),
        legend.title = element_text(size=tsize),
        legend.text = element_text(size=osize),
        legend.position = "bottom",
        legend.box = "vertical",
        panel.grid.major = element_line(colour = "grey85",linetype="dashed"),
        panel.grid.minor = element_line(colour = "grey90"),
        legend.key = element_rect(fill = "transparent", colour = "transparent"),
        legend.key.size = unit(.8, "cm"))
ggsave("Plots/ValidStations.jpg", width=20, height=24, units=c("cm"),dpi=1500)

## Scatter plot of quantiles -------------
fileobs="out/obs_percentilesAY.csv"
filesim="out/EFAS_percentilesAY.csv"
obs=read.csv(paste0(valid_path,fileobs))
sim=read.csv(paste0(valid_path,filesim))
names(obs)[1] = names(sim)[1]="Station_ID"
names(obs)[c(2:100)] = names(sim)[c(2:100)]=paste('q',1:99,sep="_")

#quantile selection:
qtil=05
myq=paste0("q_",qtil)
if (qtil<10) qtil=paste0(0,qtil)
nobs=paste0("Q",qtil,"_obs")
nsim=paste0("Q",qtil,"_sim")
scplot=paste0("scatterplot_Q",qtil)
qloc=match(myq,colnames(obs))
def my_function():
  # This is a comment explaining the purpose of this function
  print("Hello, world!")
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
  annotate("label", x=450000, y=500, label= paste0("n = ",length(ppl$UpA)),size=6)

p
ggsave("histo_stations_f.jpg", p, width=20, height=15, units=c("cm"),dpi=1500)

#Map of station used in validation (number of years on record) -Decommisionned

palet=c(hcl.colors(9, palette = "viridis", alpha = NULL, rev = T, fixup = TRUE))

ggplot(basemap) +
  geom_sf(fill="gray95", color=NA) +
  geom_sf(data=parpl,aes(geometry=geometry,col=years),alpha=.7,size=1.4,stroke=0,shape=16)+ 
  geom_sf(fill=NA, color="grey20") +
  geom_sf(data=parpef,aes(geometry=geometry),col="black",alpha=1,size=1.4,stroke=0.5,shape=1)+ 
  scale_x_continuous(breaks=seq(-30,40, by=5)) +
  scale_color_gradientn(
    colors=palet, name="length of record (Years)",breaks=c(0,10,20,30,40,50,60,70),oob = scales::squish)  +
  labs(x="Longitude", y = "Latitude")+
  coord_sf(xlim = c(min(nco[,1]),max(nco[,1])), ylim = c(min(nco[,2]),max(nco[,2])))+
  # scale_color_manual(values=c("EFAS" ="blue","FloodDrivers"="red","Both"="green"),name="River Gauges")   +
  labs(x="Longitude", y = "Latitude")+
  guides(color = guide_bins(barwidth = 10, barheight = 2,reverse=F,override.aes = list(size = 5, shape=15)))+
  theme(axis.title=element_text(size=tsize),
        panel.background = element_rect(fill = "aliceblue", colour = "grey1"),
        panel.border = element_rect(linetype = "solid", fill = NA, colour="black"),
        legend.title = element_text(size=tsize),
        legend.text = element_text(size=osize),
        legend.position = "bottom",
        panel.grid.major = element_line(colour = "grey85",linetype="dashed"),
        panel.grid.minor = element_line(colour = "grey90"),
        legend.key = element_rect(fill = "transparent", colour = "transparent"),
        legend.key.size = unit(.8, "cm"))

ggsave("Validation_gauges.jpg", width=20, height=20, units=c("cm"),dpi=1500)


### Boxplot of performance by decade -----
#groups= "1950-1960","1960-1970","1970-1980","1980-1990","1990-2000","2000-2010","2010-2020"

kgefile="out/EFAS_obs_kge.csv"
skill=read.csv(paste0(valid_path,kgefile))
if (length(skill[1,])>2){
  years=c(1950:2020)
  skill=skill[-73]
  names(skill)[c(2:72)]=years
}
names(skill)[1]="Station_ID"


#Loading all txt files with matching locations: I want to have all matches ==> I can get the area
#load matching coordinates
Sloc=read.table(paste0(valid_path,"out/Stations_locations_EFAS_grid_1950.txt"),sep=",")
yrlist=c(1951:2020)
chd=seq(1960,2020,by=10)
Slolos=list()
ix=0
lslo=rep(NA,length(chd))
for (yi in yrlist){
  print(yi)
  if (length(which(yi==chd))>0){
    ix=ix+1
    Sloc_final=Sloc[which(Sloc$V2!=0),]
    onlyV=which(!is.na(match(Sloc_final$V1,ValidSY$V1)))
    Sloc_fx=Sloc_final[onlyV,]
    Slolos[[ix]]=Sloc_fx
    lslo[ix]=length(Sloc_fx$V1)
    Sloc=read.table(paste0(valid_path,"out/Stations_locations_EFAS_grid_",yi,".txt"),sep=",")
  }else{
    Slocy=read.table(paste0(valid_path,"out/Stations_locations_EFAS_grid_",yi,".txt"),sep=",")
    Slocrep=Slocy[which(Sloc$V2==0),]
    Sloc[which(Sloc$V2==0),]=Slocrep
  }
}

Sloc_final=Sloc[which(Sloc$V2!=0),]
Sloc_final$csource="SpatialQMatch"

#Only keep valid stations
skillm<-melt(skill,id.vars ="Station_ID")
skillm$year=as.numeric(as.character(skillm$variable))
skillm$decade=1950
skillm$decade[which(skillm$year>1959 & skillm$year<1970)]=1960
skillm$decade[which(skillm$year>1969 & skillm$year<1980)]=1970
skillm$decade[which(skillm$year>1979 & skillm$year<1990)]=1980
skillm$decade[which(skillm$year>1989 & skillm$year<2000)]=1990
skillm$decade[which(skillm$year>1999 & skillm$year<2010)]=2000
skillm$decade[which(skillm$year>2009)]=2010
#average
skillm$value[which(is.infinite(skillm$value))]=NA
dat_skill=aggregate(list(skill=skillm$value),
                      by = list(Station_ID=skillm$Station_ID,decade=skillm$decade),
                      FUN = function(x) c(mean= mean(x,na.rm=T),median=median(x,na.rm=T),max=max(x,na.rm=T)))
dat_skill <- do.call(data.frame, dat_skill)
onlyV=which(!is.na(match(dat_skill$Station_ID,ValidSY$V1)))
dat_skill=dat_skill[onlyV,]

plot(dat_skill$skill.median)

palet=c(hcl.colors(9, palette = "Blues", alpha = NULL, rev = T, fixup = TRUE))

meds <- c(by(dat_skill$skill.median, dat_skill$decade, median,na.rm=T))
qshit <- c(by(dat_skill$skill.median, dat_skill$decade, quantile ,na.rm=T))
merdecol=match(dat_skill$decade,names(meds))
dat_skill$col=meds[merdecol]

labelsX=c("1950-1960","1960-1970","1970-1980","1980-1990","1990-2000","2000-2010","2010-2020")

p3<-ggplot(dat_skill, aes(x=factor(decade), y=skill.median)) +
  # labs(title=tit1,x="Warmimg levels (?C)", y = y1)+
  geom_boxplot(notch=F,position=position_dodge(.9),alpha=.8,aes(fill=col),linewidth=0.8,outlier.alpha = 0.1)+
  #geom_violin(aes(fill=col),scale="width",draw_quantiles = c(0.25, 0.5, 0.75),position=position_dodge(.9),alpha=0.7,lwd=1)+
  scale_fill_gradientn(
    colors=palet, n.breaks=6,limits=c(0.4,0.6)) +
  scale_y_continuous(limits = c(-0.5,1),name="KGE",breaks = seq(-1,1,by=0.5),minor_breaks = seq(-1,1,0.1))+
  scale_x_discrete(name="Decade",labels=labelsX)+
  geom_text(data=data.frame(), aes(x=names(meds), y=meds-0.05, label=round(lslo,2)), col='black', size=4,fontface="bold")+
  theme(axis.title=element_text(size=16, face="bold"),
        axis.text = element_text(size=12),
        panel.background = element_rect(fill = "white", colour = "white"),
        panel.grid = element_blank(),
        panel.border = element_rect(linetype = "solid", fill = NA, colour="black"),
        legend.title = element_text(size=14),
        legend.text = element_text(size=12),
        legend.position="none",
        panel.grid.major = element_line(colour = "grey60"),
        panel.grid.minor.y = element_line(colour = "grey80",linetype="dashed"),
        legend.key = element_rect(fill = "transparent", colour = "transparent"),
        legend.key.size = unit(.8, "cm"))
p3
ggsave("Plots/boxp_KGEtime.jpg", p3, width=20, height=15, units=c("cm"),dpi=1500)



### Boxplot of performances according to the presence or not of reservoirs ------------------------

#load reservoir ratio map
res_path<-("D:/tilloal/Documents/LFRuns_utils/data/reservoirs/")
outletname="res_ratio_European_01min.nc"
dir=res_path
ResData=resOpen(res_path,outletname,ValidSY)  
rsx=ResData[[2]]
length(rsx$res.ratio[which(rsx$res.ratio<0.5)])/length(rsx$res.ratio)
ValidSYp=ResData
#do two boxplots, one for res.ratio>0.5 and one for others

ValidSYp$res.group=0
ValidSYp$res.group[which(ValidSYp$res.ratio>=0.5 & ValidSYp$res.ratio<=1)]=1
ValidSYp$res.group[which(ValidSYp$res.ratio>1)]=2

meds <- c(by(ValidSYp$skill, ValidSYp$res.group, median))
len <- c(by(ValidSYp$skill, ValidSYp$res.group, length))
merdecol=match(ValidSYp$res.group,names(meds))
ValidSYp$col=meds[merdecol]
palet=c(hcl.colors(9, palette = "Blues", alpha = NULL, rev = T, fixup = TRUE))
p4<-ggplot(ValidSYp, aes(x=factor(res.group), y=skill)) +
  # scale_fill_manual(values = cols, name = "Socioeconomic \n Pathways")+
  # labs(title=tit1,x="Warmimg levels (?C)", y = y1)+
  geom_boxplot(notch=F,position=position_dodge(.9),alpha=.8,aes(fill=col),linewidth=0.8,outlier.alpha = 0.4)+
  scale_y_continuous(limits = c(-0.5,1),name="KGE",breaks = seq(-1,1,by=0.5),minor_breaks = seq(-1,1,0.1))+
  scale_x_discrete(labels=c("0" = "c < 0.5 \n (low reservoir impact)", "1" = "0.5 <= c <= 1  \n (medium reservoir impact)",
                            "2" = "c > 1 \n (high reservoir impact)"),name="c ratio (reservoir volume to mean annual streamflow) ")+
  scale_fill_gradientn(
    colors=palet, n.breaks=6,limits=c(0.2,0.8)) +
  geom_text(data=data.frame(), aes(x=names(meds), y=meds-0.05, label=len), col='black', size=4,fontface="bold")+
  theme(axis.title=element_text(size=16, face="bold"),
        axis.text = element_text(size=12),
        panel.background = element_rect(fill = "white", colour = "white"),
        panel.grid = element_blank(),
        panel.border = element_rect(linetype = "solid", fill = NA, colour="black"),
        legend.title = element_text(size=14),
        legend.text = element_text(size=12),
        legend.position = "none",
        panel.grid.major = element_line(colour = "grey60"),
        panel.grid.minor.y = element_line(colour = "grey80",linetype="dashed"),
        legend.key = element_rect(fill = "transparent", colour = "transparent"),
        legend.key.size = unit(.8, "cm"))
p4
ggsave("Plots/boxp_KGE_reservoirs.jpg", p4, width=20, height=15, units=c("cm"),dpi=1500)

## Performance on annual maxima and minima -----------------------
Q_data <- read.csv(paste0('Q_19502020.csv'), header = F)  # CSVs with observations
Q_sim <- read.csv(file="out/EFAS_19502020.csv")
Station_data_IDo <- as.numeric(as.vector(t(Q_data[1, -c(1)])))
Station_data_IDs <- as.vector(t(Q_sim[1, ]))
# kind of hit rate: 2 boxplots
hit_out=c()
rs_obs=c()
rs_sim=c()
for (id in 1:length(Station_data_IDo)){
  print(id)
  s=Station_data_IDo[id]
  # s=6574362
  ms=Q_data[-c(1,2,3,4),c(1,id+1)]
  ms$V2=as.numeric(ms$V2)
  ms$time=as.Date(ms$V2-1,origin="0000-01-01")
  ms[,2]=as.numeric(ms[,2])
  ms=data.frame(time=ms$time,data=ms[,2])
  #extract peak moments
  sl=which(Station_data_IDs==s)
  mss=data.frame(time=ms$time,data=Q_sim[-1,sl])
  iy=which(is.na(ms$data))
  hit=NA
  if (length(iy)>0){
    ms=ms[-iy,]
    mss=mss[-iy,]
  }
  if (length(ms$data>0)){
    AMAX_obs <- computeAnnualMaxima(ms)
    AMAX_sim <- computeAnnualMaxima(mss)
    if (length(AMAX_sim$annualMaxIndx)==length(AMAX_obs$annualMaxIndx)){
      AMAX_diff=AMAX_obs$annualMaxIndx-AMAX_sim$annualMaxIndx
      hit=length(which(abs(AMAX_diff)<8))/length(AMAX_diff)
    }
    data=data.frame(ms,sim=mss$data)
    names(data)=c("date","Q","Qs")
    res_obs=data.frame(Station=rep(s,length(AMAX_obs$annualMax)),AnMax=AMAX_obs$annualMax,AnMaxDate=AMAX_obs$annualMaxDate)
    res_sim=data.frame(Station=rep(s,length(AMAX_sim$annualMax)),AnMax=AMAX_sim$annualMax,AnMaxDate=AMAX_sim$annualMaxDate)
    # for (ev in 1:length(AMAX_obs$annualMaxIndx)){
    #   #print(ev)
    #   test=StatEpisode(data, date.deb = data$date[AMAX_obs$annualMaxIndx[ev]], win=30, qx.min = 100, mode = "FL",obs=TRUE)
    #   res=rbind(res,test)
    # }
    # save=c(Station_ID=s,rmean=mean(res[,2],na.rm=T),errmean=mean(res[,3],na.rm=T))
    rs_obs=rbind(rs_obs,res_obs)
    rs_sim=rbind(rs_sim,res_sim)
  }else{print("no observations")}
  hit_out=c(hit_out,hit)
}

hit_out=data.frame(hit_out,Station_data_IDo)
hitm=match(ValidSY$V1,hit_out$Station_data_IDo)
hit_out=hit_out[hitm,]
boxplot(hit_out$hit_out)
#1- mean on Anmax for every station
catagg_obs=aggregate(list(Q=rs_obs$AnMax),
                   by = list(Station=rs_obs$Station),
                   FUN = function(x) c(med=median(x,na.rm=T),mean=mean(x,na.rm=T),dev=sd(x,na.rm=T),len=length(x)))
catagg_obs <- do.call(data.frame, catagg_obs)

catagg_sim=aggregate(list(Q=rs_sim$AnMax),
                     by = list(Station=rs_sim$Station),
                     FUN = function(x) c(med=median(x,na.rm=T),mean=mean(x,na.rm=T),dev=sd(x,na.rm=T),len=length(x)))
catagg_sim <- do.call(data.frame, catagg_sim)

catagg_obssim=inner_join(catagg_obs,catagg_sim, by="Station",suffix = c("obs","sim"))

plot(catagg_obssim$Q.meanobs,catagg_obssim$Q.meansim)

#reataining only good matches
matV=match(ValidSY$V1,catagg_obssim$Station)

Catagg_f=catagg_obssim[matV,]
Catagg_f$dPeaks=(Catagg_f$Q.meansim-Catagg_f$Q.meanobs)/Catagg_f$Q.meanobs*100
plot(Catagg_f$Q.meanobs,Catagg_f$Q.meansim)


#Spatial plot, where is dpeak shit?
palet=c(hcl.colors(9, palette = "RdBu", alpha = NULL, rev = F, fixup = TRUE))
points=inner_join(ValidSY,Catagg_f,by=c("V1"="Station"))
parpl <- st_as_sf(points, coords = c("Var1", "Var2"), crs = 4326)
parpl <- st_transform(parpl, crs = 3035)
cord.dec=points[,c(1,2)]
cord.dec = SpatialPoints(cord.dec, proj4string=CRS("+proj=longlat"))
cord.UTM <- spTransform(cord.dec, CRS("+init=epsg:3035"))
nco=cord.UTM@coords
world <- ne_countries(scale = "medium", returnclass = "sf")
Europe <- world[which(world$continent == "Europe"),]
e2=st_transform(Europe,  crs=3035)
w2=st_transform(world,  crs=3035)
tsize=12
osize=12
limi=c(-100,100)
basemap=w2
ggplot(basemap) +
  geom_sf(fill="gray95", color=NA) +
  geom_sf(data=parpl,aes(geometry=geometry,fill=dPeaks,size=UpA),color="transparent",alpha=.9,shape=21,stroke=0)+ 
  geom_sf(fill=NA, color="grey20") +
  geom_sf(data=parpl,aes(geometry=geometry,size=UpA),col="black",alpha=1,stroke=0.1,shape=1)+ 
  scale_x_continuous(breaks=seq(-30,40, by=5)) +
  scale_size(range = c(1, 3), trans="sqrt")+
  scale_fill_gradientn(
    colors=palet,n.breaks=5,oob = scales::squish, limits=limi, name="mean Qpeak difference") +
  coord_sf(xlim = c(min(nco[,1]),max(nco[,1])), ylim = c(min(nco[,2]),max(nco[,2])))+
  labs(x="Longitude", y = "Latitude") +
  guides(fill = guide_colourbar(barwidth = 0.5, barheight = 12,reverse=F),size="none")+
  theme(axis.title=element_text(size=tsize),
        panel.background = element_rect(fill = "aliceblue", colour = "grey1"),
        panel.border = element_rect(linetype = "solid", fill = NA, colour="black"),
        legend.title = element_text(size=tsize),
        legend.text = element_text(size=osize),
        legend.position = "right",
        panel.grid.major = element_line(colour = "grey85",linetype="dashed"),
        panel.grid.minor = element_line(colour = "grey90"),
        legend.key = element_rect(fill = "transparent", colour = "transparent"),
        legend.key.size = unit(.8, "cm"))
ggsave("Plots/Spatial_Qpeakdiff.jpg", width=20, height=15, units=c("cm"),dpi=1500)

median(Catagg_f$dPeaks)
boxplot(Catagg_f$dPeaks,ylim=c(-100,100))

matO=which(!is.na(match(rs_obs$Station,ValidSY$V1)))
matS=which(!is.na(match(rs_sim$Station,ValidSY$V1)))

rs_obsF=rs_obs[matO,]
rs_simF=rs_sim[matS,]

rs_obsF$yday=yday(rs_obsF$AnMaxDate)
seasonMax_obs=aggregate(list(Q=rs_obsF$yday),
                        by = list(Station=rs_obsF$Station),
                        FUN = function(x) c(season=season1(x),len=length(x)))
seasonMax_obs <- do.call(data.frame, seasonMax_obs)

rs_simF$yday=yday(rs_simF$AnMaxDate)
seasonMax_sim=aggregate(list(Q=rs_simF$yday),
                        by = list(Station=rs_simF$Station),
                        FUN = function(x) c(season=season1(x),len=length(x)))
seasonMax_sim <- do.call(data.frame, seasonMax_sim)

seasonMax_obs$dseason=seasonMax_obs$Q.season-seasonMax_sim$Q.season
seasonMax_obs$dseason[which(seasonMax_obs$dseason>182)]=seasonMax_obs$dseason[which(seasonMax_obs$dseason>182)]-365.25
seasonMax_obs$dseason[which(seasonMax_obs$dseason<=-182)]=365.25+seasonMax_obs$dseason[which(seasonMax_obs$dseason<=-182)]

hist(seasonMax_obs$dseason)
palet=c(hcl.colors(9, palette = "RdBu", alpha = NULL, rev = F, fixup = TRUE))
points=inner_join(ValidSY,seasonMax_obs,by=c("V1"="Station"))
parpl <- st_as_sf(points, coords = c("Var1", "Var2"), crs = 4326)
parpl <- st_transform(parpl, crs = 3035)
cord.dec=points[,c(1,2)]
cord.dec = SpatialPoints(cord.dec, proj4string=CRS("+proj=longlat"))
cord.UTM <- spTransform(cord.dec, CRS("+init=epsg:3035"))
nco=cord.UTM@coords
world <- ne_countries(scale = "medium", returnclass = "sf")
Europe <- world[which(world$continent == "Europe"),]
e2=st_transform(Europe,  crs=3035)
w2=st_transform(world,  crs=3035)
tsize=12
osize=12
limi=c(-100,100)
basemap=w2
ggplot(basemap) +
  geom_sf(fill="gray95", color=NA) +
  geom_sf(data=parpl,aes(geometry=geometry,fill=dseason,size=UpA),color="transparent",alpha=.9,shape=21,stroke=0)+ 
  geom_sf(fill=NA, color="grey20") +
  geom_sf(data=parpl,aes(geometry=geometry,size=UpA),col="black",alpha=1,stroke=0.1,shape=1)+ 
  scale_x_continuous(breaks=seq(-30,40, by=5)) +
  scale_size(range = c(1, 3), trans="sqrt")+
  scale_fill_gradientn(
    colors=palet,n.breaks=5,oob = scales::squish, limits=limi, name="Qpeak season difference") +
  coord_sf(xlim = c(min(nco[,1]),max(nco[,1])), ylim = c(min(nco[,2]),max(nco[,2])))+
  labs(x="Longitude", y = "Latitude") +
  guides(fill = guide_colourbar(barwidth = 0.5, barheight = 12,reverse=F),size="none")+
  theme(axis.title=element_text(size=tsize),
        panel.background = element_rect(fill = "aliceblue", colour = "grey1"),
        panel.border = element_rect(linetype = "solid", fill = NA, colour="black"),
        legend.title = element_text(size=tsize),
        legend.text = element_text(size=osize),
        legend.position = "right",
        panel.grid.major = element_line(colour = "grey85",linetype="dashed"),
        panel.grid.minor = element_line(colour = "grey90"),
        legend.key = element_rect(fill = "transparent", colour = "transparent"),
        legend.key.size = unit(.8, "cm"))
ggsave("Plots/Spatial_Speakdiff.jpg", width=20, height=15, units=c("cm"),dpi=1500)


#Same for low flows
hit_outl=c()
rsl_obs=c()
rsl_sim=c()
for (id in 1:length(Station_data_IDo)){
  print(id)
  s=Station_data_IDo[id]
  # s=6574362
  ms=Q_data[-c(1,2,3,4),c(1,id+1)]
  ms$V2=as.numeric(ms$V2)
  ms$time=as.Date(ms$V2-1,origin="0000-01-01")
  ms[,2]=as.numeric(ms[,2])
  ms=data.frame(time=ms$time,data=ms[,2])
  #extract peak moments
  sl=which(Station_data_IDs==s)
  mss=data.frame(time=ms$time,data=Q_sim[-1,sl])
  iy=which(is.na(ms$data))
  hit=NA
  if (length(iy)>0){
    ms=ms[-iy,]
    mss=mss[-iy,]
  }
  if (length(ms$data>0)){
    AMIN_obs <- computeAnnualMinima(ms)
    AMIN_sim <- computeAnnualMinima(mss)
    if (length(AMIN_sim$annualMinIndx)==length(AMIN_obs$annualMinIndx)){
      AMIN_diff=AMIN_obs$annualMinIndx-AMIN_sim$annualMinIndx
      hit=length(which(abs(AMIN_diff)<15))/length(AMIN_diff)
    }
    data=data.frame(ms,sim=mss$data)
    names(data)=c("date","Q","Qs")
    res_obs=data.frame(Station=rep(s,length(AMIN_obs$annualMin)),AnMin=AMIN_obs$annualMin,AnMinDate=AMIN_obs$annualMinDate)
    res_sim=data.frame(Station=rep(s,length(AMIN_sim$annualMin)),AnMin=AMIN_sim$annualMin,AnMinDate=AMIN_sim$annualMinDate)
    # for (ev in 1:length(AMAX_obs$annualMaxIndx)){
    #   #print(ev)
    #   test=StatEpisode(data, date.deb = data$date[AMAX_obs$annualMaxIndx[ev]], win=30, qx.min = 100, mode = "FL",obs=TRUE)
    #   res=rbind(res,test)
    # }
    # save=c(Station_ID=s,rmean=mean(res[,2],na.rm=T),errmean=mean(res[,3],na.rm=T))
    rsl_obs=rbind(rsl_obs,res_obs)
    rsl_sim=rbind(rsl_sim,res_sim)
  }else{print("no observations")}
  hit_outl=c(hit_outl,hit)
}
hit_outl=data.frame(hit_outl,Station_data_IDo)
hitm=match(ValidSY$V1,hit_outl$Station_data_IDo)
hit_outl=hit_outl[hitm,]
boxplot(hit_outl$hit_out)
#1- mean on Anmax for every station
cataggl_obs=aggregate(list(Q=rsl_obs$AnMin),
                     by = list(Station=rsl_obs$Station),
                     FUN = function(x) c(med=median(x,na.rm=T),mean=mean(x,na.rm=T),dev=sd(x,na.rm=T),len=length(x)))
cataggl_obs <- do.call(data.frame, cataggl_obs)

cataggl_sim=aggregate(list(Q=rsl_sim$AnMin),
                     by = list(Station=rsl_sim$Station),
                     FUN = function(x) c(med=median(x,na.rm=T),mean=mean(x,na.rm=T),dev=sd(x,na.rm=T),len=length(x)))
cataggl_sim <- do.call(data.frame, cataggl_sim)

cataggl_obssim=inner_join(cataggl_obs,cataggl_sim, by="Station",suffix = c("obs","sim"))

plot(cataggl_obssim$Q.meanobs,cataggl_obssim$Q.meansim)

#reataining only good matches
matV=match(ValidSY$V1,cataggl_obssim$Station)

Cataggl_f=cataggl_obssim[matV,]
Cataggl_f$dPeaks=(Cataggl_f$Q.meansim-Cataggl_f$Q.meanobs)/Cataggl_f$Q.meanobs*100
plot(Cataggl_f$Q.meanobs,Cataggl_f$Q.meansim)


#Spatial plot, where is dpeak shit?
palet=c(hcl.colors(9, palette = "RdBu", alpha = NULL, rev = F, fixup = TRUE))
points=inner_join(ValidSY,Cataggl_f,by=c("V1"="Station"))
parpl <- st_as_sf(points, coords = c("Var1", "Var2"), crs = 4326)
parpl <- st_transform(parpl, crs = 3035)
cord.dec=points[,c(1,2)]
cord.dec = SpatialPoints(cord.dec, proj4string=CRS("+proj=longlat"))
cord.UTM <- spTransform(cord.dec, CRS("+init=epsg:3035"))
nco=cord.UTM@coords
world <- ne_countries(scale = "medium", returnclass = "sf")
Europe <- world[which(world$continent == "Europe"),]
e2=st_transform(Europe,  crs=3035)
w2=st_transform(world,  crs=3035)
tsize=12
osize=12
limi=c(-100,100)
basemap=w2
ggplot(basemap) +
  geom_sf(fill="gray95", color=NA) +
  geom_sf(data=parpl,aes(geometry=geometry,fill=dMin,size=UpA),color="transparent",alpha=.9,shape=21,stroke=0)+ 
  geom_sf(fill=NA, color="grey20") +
  geom_sf(data=parpl,aes(geometry=geometry,size=UpA),col="black",alpha=1,stroke=0.1,shape=1)+ 
  scale_x_continuous(breaks=seq(-30,40, by=5)) +
  scale_size(range = c(1, 3), trans="sqrt")+
  scale_fill_gradientn(
    colors=palet,n.breaks=5,oob = scales::squish, limits=limi, name="mean Qpeak difference") +
  coord_sf(xlim = c(min(nco[,1]),max(nco[,1])), ylim = c(min(nco[,2]),max(nco[,2])))+
  labs(x="Longitude", y = "Latitude") +
  guides(fill = guide_colourbar(barwidth = 0.5, barheight = 12,reverse=F),size="none")+
  theme(axis.title=element_text(size=tsize),
        panel.background = element_rect(fill = "aliceblue", colour = "grey1"),
        panel.border = element_rect(linetype = "solid", fill = NA, colour="black"),
        legend.title = element_text(size=tsize),
        legend.text = element_text(size=osize),
        legend.position = "right",
        panel.grid.major = element_line(colour = "grey85",linetype="dashed"),
        panel.grid.minor = element_line(colour = "grey90"),
        legend.key = element_rect(fill = "transparent", colour = "transparent"),
        legend.key.size = unit(.8, "cm"))
ggsave("Plots/Spatial_QMindiff.jpg", width=20, height=15, units=c("cm"),dpi=1500)

matO=which(!is.na(match(rsl_obs$Station,ValidSY$V1)))
matS=which(!is.na(match(rsl_sim$Station,ValidSY$V1)))

rsl_obsF=rsl_obs[matO,]
rsl_simF=rsl_sim[matS,]

rsl_obsF$yday=yday(rsl_obsF$AnMinDate)
seasonMin_obs=aggregate(list(Q=rsl_obsF$yday),
                        by = list(Station=rsl_obsF$Station),
                        FUN = function(x) c(season=season1(x),len=length(x)))
seasonMin_obs <- do.call(data.frame, seasonMin_obs)

rsl_simF$yday=yday(rsl_simF$AnMinDate)
seasonMin_sim=aggregate(list(Q=rsl_simF$yday),
                        by = list(Station=rsl_simF$Station),
                        FUN = function(x) c(season=season1(x),len=length(x)))
seasonMin_sim <- do.call(data.frame, seasonMin_sim)

seasonMin_obs$dseason=seasonMin_obs$Q.season-seasonMin_sim$Q.season
seasonMin_obs$dseason[which(seasonMin_obs$dseason>182)]=seasonMin_obs$dseason[which(seasonMin_obs$dseason>182)]-365.25
seasonMin_obs$dseason[which(seasonMin_obs$dseason<=-182)]=365.25+seasonMin_obs$dseason[which(seasonMin_obs$dseason<=-182)]

hist(seasonMin_obs$dseason)
palet=c(hcl.colors(9, palette = "RdBu", alpha = NULL, rev = F, fixup = TRUE))
points=inner_join(ValidSY,seasonMin_obs,by=c("V1"="Station"))
parpl <- st_as_sf(points, coords = c("Var1", "Var2"), crs = 4326)
parpl <- st_transform(parpl, crs = 3035)
cord.dec=points[,c(1,2)]
cord.dec = SpatialPoints(cord.dec, proj4string=CRS("+proj=longlat"))
cord.UTM <- spTransform(cord.dec, CRS("+init=epsg:3035"))
nco=cord.UTM@coords
world <- ne_countries(scale = "medium", returnclass = "sf")
Europe <- world[which(world$continent == "Europe"),]
e2=st_transform(Europe,  crs=3035)
w2=st_transform(world,  crs=3035)
tsize=12
osize=12
limi=c(-100,100)
basemap=w2
ggplot(basemap) +
  geom_sf(fill="gray95", color=NA) +
  geom_sf(data=parpl,aes(geometry=geometry,fill=dseason,size=UpA),color="transparent",alpha=.9,shape=21,stroke=0)+ 
  geom_sf(fill=NA, color="grey20") +
  geom_sf(data=parpl,aes(geometry=geometry,size=UpA),col="black",alpha=1,stroke=0.1,shape=1)+ 
  scale_x_continuous(breaks=seq(-30,40, by=5)) +
  scale_size(range = c(1, 3), trans="sqrt")+
  scale_fill_gradientn(
    colors=palet,n.breaks=5,oob = scales::squish, limits=limi, name="Qpeak season difference") +
  coord_sf(xlim = c(min(nco[,1]),max(nco[,1])), ylim = c(min(nco[,2]),max(nco[,2])))+
  labs(x="Longitude", y = "Latitude") +
  guides(fill = guide_colourbar(barwidth = 0.5, barheight = 12,reverse=F),size="none")+
  theme(axis.title=element_text(size=tsize),
        panel.background = element_rect(fill = "aliceblue", colour = "grey1"),
        panel.border = element_rect(linetype = "solid", fill = NA, colour="black"),
        legend.title = element_text(size=tsize),
        legend.text = element_text(size=osize),
        legend.position = "right",
        panel.grid.major = element_line(colour = "grey85",linetype="dashed"),
        panel.grid.minor = element_line(colour = "grey90"),
        legend.key = element_rect(fill = "transparent", colour = "transparent"),
        legend.key.size = unit(.8, "cm"))
ggsave("Plots/Spatial_SMindiff.jpg", width=20, height=15, units=c("cm"),dpi=1500)

names(hit_out)=names(hit_outl)
Perfres=rbind(hit_out,hit_outl)
Perfres$hit_outl=Perfres$hit_outl*100
Perfres$group=1
Perfres$group[c(2902:5802)]=2
ggplot(Perfres, aes(x=factor(group), y=hit_outl)) +
  # labs(title=tit1,x="Warmimg levels (?C)", y = y1)+
  geom_boxplot(notch=F,position=position_dodge(.9),alpha=.8,aes(fill=factor(group)),linewidth=0.8,outlier.alpha = 0.4)+
  scale_y_continuous(limits = c(0,100),name="Hit rate (%)",breaks = seq(0,100,by=10),minor_breaks = seq(-100,100,5))+
  scale_x_discrete(labels=c("1" = "Annual maxs", "2" = "Annual mins"),name="")+
  scale_fill_manual(values = c("1"="royalblue","2" ="tomato"), name = "")+
  theme(axis.title=element_text(size=16, face="bold"),
        axis.text = element_text(size=12),
        panel.background = element_rect(fill = "white", colour = "white"),
        panel.grid = element_blank(),
        panel.border = element_rect(linetype = "solid", fill = NA, colour="black"),
        legend.title = element_text(size=14),
        legend.text = element_text(size=12),
        legend.position = "none",
        panel.grid.major = element_line(colour = "grey60"),
        panel.grid.minor.y = element_line(colour = "grey80",linetype="dashed"),
        legend.key = element_rect(fill = "transparent", colour = "transparent"),
        legend.key.size = unit(.8, "cm"))
ggsave("Plots/Bxplot_Hitrate.jpg", width=20, height=15, units=c("cm"),dpi=1500)



Perfres=rbind(seasonMax_obs[,c(1,4)],seasonMin_obs[,c(1,4)])
Perfres$group=1
Perfres$group[c(2902:5802)]=2

ggplot(Perfres, aes(x=factor(group), y=dseason)) +
  # labs(title=tit1,x="Warmimg levels (?C)", y = y1)+
  geom_violin(aes(fill=factor(group)),position=position_dodge(.9),alpha=0.8,size=0.9)+
  geom_boxplot(width=0.05,notch=F,size=0.8)+
  coord_flip()+
  #geom_boxplot(notch=F,position=position_dodge(.9),alpha=.8,aes(fill=factor(group)),linewidth=0.8,outlier.alpha = 0.4)+
  scale_y_continuous(limits = c(-60,60),name="Error (days)",breaks = seq(-60,60,by=10),minor_breaks = seq(-100,100,5))+
  scale_x_discrete(labels=c("1" = "Annual maxs", "2" = "Annual mins"),name="")+
  scale_fill_manual(values = c("1"="#9B98F6","2" ="#FF8888"), name = "")+
  theme(axis.title=element_text(size=16),
        axis.text = element_text(size=14),
        panel.background = element_rect(fill = "white", colour = "white"),
        panel.grid = element_blank(),
        panel.border = element_rect(linetype = "solid", fill = NA, colour="black"),
        legend.title = element_text(size=14),
        legend.text = element_text(size=12),
        legend.position = "none",
        panel.grid.major = element_line(colour = "grey60"),
        panel.grid.minor = element_line(colour = "grey80",linetype="dashed"),
        legend.key = element_rect(fill = "transparent", colour = "transparent"),
        legend.key.size = unit(.8, "cm")) 

ggsave("Plots/Bxplot_seasonF.jpg", width=20, height=7, units=c("cm"),dpi=1500)



Perfres=rbind(hit_out,Cataggl_f[,c(1,10)])
Perfres$group=1
Perfres$group[c(2902:5802)]=2
ggplot(Perfres, aes(x=factor(group), y=dPeaks)) +
  # labs(title=tit1,x="Warmimg levels (?C)", y = y1)+
  #geom_violin(scale="width",draw_quantiles = c(0.25, 0.5, 0.75),position=position_dodge(.9),alpha=0.7)
  geom_boxplot(notch=F,position=position_dodge(.9),alpha=.8,aes(fill=factor(group)),linewidth=0.8,outlier.alpha = 0.4)+
  scale_y_continuous(limits = c(-100,100),name="Error (%)",breaks = seq(-100,100,by=20),minor_breaks = seq(-100,100,10))+
  scale_x_discrete(labels=c("1" = "Flood peaks", "2" = "Drought lows"),name="")+
  scale_fill_manual(values = c("1"="royalblue","2" ="tomato"), name = "")+
  theme(axis.title=element_text(size=16, face="bold"),
        axis.text = element_text(size=12),
        panel.background = element_rect(fill = "white", colour = "white"),
        panel.grid = element_blank(),
        panel.border = element_rect(linetype = "solid", fill = NA, colour="black"),
        legend.title = element_text(size=14),
        legend.text = element_text(size=12),
        legend.position = "none",
        panel.grid.major = element_line(colour = "grey60"),
        panel.grid.minor.y = element_line(colour = "grey80",linetype="dashed"),
        legend.key = element_rect(fill = "transparent", colour = "transparent"),
        legend.key.size = unit(.8, "cm"))
ggsave("Plots/Bxplot_peakdiff.jpg", width=20, height=15, units=c("cm"),dpi=1500)

#Now I just need to redu the scatterplot for different quantiles and goodbye

obsinr=rbind(Cataggl_f[,c(1,2,3)])
siminr=rbind(Cataggl_f[,c(1,6,7)])
names(obsinr)=names(siminr)=names(obsin)
obsinr=obsinr[which(obsinr$value>0),]
siminr=siminr[which(siminr$value>0),]
nobs="AnMin_obs"
nsim="AnMin_sim"
scplot="scatterplot_AnMin"
Myplot2=ScatterValid(obsinr,nobs,siminr,nsim,scplot,valid_path,ValidSY)

Myplot2
ggsave(paste0("Plots/",scplot,".jpg"), width=20, height=15, units=c("cm"),dpi=1500)



outf=c()
for( Nsq in 1:88){
  print(Nsq)
  rspace= read.csv(paste0(hydroDir,"/subspace_efas.csv"))
  rspace=rspace[,-1]
  nrspace=rspace[Nsq,]
  outletname="efas_rnet_100km_01min"
  
  outhybas=outletopen(hydroDir,outletname,nrspace)
  Idstart=as.numeric(Nsq)*100000
  if (length(outhybas$outlets)>0){
    outhybas$outlets=seq((Idstart+1),(Idstart+length(outhybas$outlets)))
    outhybas$latlong=paste(round(outhybas$Var1,4),round(outhybas$Var2,4),sep=" ")
    # outcut=which(!is.na(match(outhybas$outlets,Paramsfl$catchment)))
    # zebi=unique(Paramsfl$catchment)
    # outhloc=outhybas[outcut,]
    outf=rbind(outf,outhybas)
  }
  
  
}


### Hybas07 ----
Catchmentrivers7=read.csv(paste0(hydroDir,"/Catchments/from_hybas_eu_onlyid.csv"),encoding = "UTF-8", header = T, stringsAsFactors = F)
hybas07 <- read_sf(dsn = paste0(hydroDir,"/Catchments/hydrosheds/hybas_eu_lev07_v1c.shp"))
hybasf7=fortify(hybas07) 
Catamere07=inner_join(hybasf7,Catchmentrivers7,by= "HYBAS_ID")
Catamere07$llcoord=paste(round(Catamere07$POINT_X,4),round(Catamere07$POINT_Y,4),sep=" ") 
#maybe useless
Catf7=inner_join(Catamere07,outf,by= c("llcoord"="latlong"))
cst7=st_transform(Catf7,  crs=3035)


data2=inner_join(ValidStar,cst7,by=c("Var1","Var2"))
# points <- st_as_sf(data2, coords = c("Var1", "Var2"), crs = 4326)
# points <- st_transform(points, crs = 3035)

#Aggregate by catchment
pointagg=aggregate(list(KGE=data2$kge.mean),
                   by = list(HYBAS_ID=data2$HYBAS_ID),
                   FUN = function(x) c(med=median(x,na.rm=T),mean=mean(x,na.rm=T),dev=sd(x,na.rm=T),len=length(x)))
pointagg <- do.call(data.frame, pointagg)


colorz=c("red","orange","hotpink","darkorchid1","royalblue","darkblue","black")


pointagg$kge=pointagg$KGE.med
pointagg$kge[which(is.infinite(pointagg$kge))]=NA
pointagg$kgecode=0
pointagg$kgecode[which(pointagg$kge>-0.41 & pointagg$kge<=0.2)]=1
pointagg$kgecode[which(pointagg$kge>0.2 & pointagg$kge<=0.5)]=2
pointagg$kgecode[which(pointagg$kge>0.5 & pointagg$kge<=0.7)]=3
pointagg$kgecode[which(pointagg$kge>0.7 & pointagg$kge<=0.8)]=4
pointagg$kgecode[which(pointagg$kge>0.8 & pointagg$kge<=0.9)]=5
pointagg$kgecode[which(pointagg$kge>0.9)]=6

pointplot=inner_join(hybasf7,pointagg,by= "HYBAS_ID")

tsize=16
osize=16
pl2=ggplot(basemap) +
  geom_sf(fill="gray95")+
  geom_sf(data=point,aes(fill=factor(kgecode),geometry=geometry),color="transparent")+ 
  geom_sf(fill=NA, color="grey") +
  geom_sf(data=parpl,aes(geometry=geometry),color="grey60",alpha=.7,size=1.4,stroke=0,shape=16)+
  geom_sf(data=parpef,aes(geometry=geometry),col="black",alpha=.7,size=1.4,stroke=0.5,shape=1)+
  coord_sf(xlim = c(min(nco[,1]),max(nco[,1])), ylim = c(min(nco[,2]),max(nco[,2])))+
  # coord_sf(xlim = c(min(nco[,1]),max(nco[,1])), ylim = c(min(nco[,2]),max(nco[,2])))+
  # geom_point(data=event,aes(x=lon,y=lat,col=dr))+
  # scale_color_distiller(palette = "Spectral",
  #                       direction = -1, limit=c(0,max(eventp$dr)),oob = scales::squish)+
  scale_fill_manual(values = colorz, labels= c("< -0.41","-0.41 - 0.2", "0,2 - 0.5","0.5 - 0.7","0.7 - 0.8","0.8 - 0.9", ">0.9"), name="KGE")   +
  # scale_fill_gradientn(
  #   colors=palet,
#     breaks=br,limits=limi,trans=trans,
#     oob = scales::squish,na.value=colNA, name=legend)   +
  labs(x="Longitude", y = "Latitude")+
  guides(colour = guide_legend(override.aes = list(size = 4)))+
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
pl2
ggsave("Validation_KGEmedcatch.jpg", width=20, height=20, units=c("cm"),dpi=1500)
#Histogram of all considered station upstream area




matp100=matpix[which(matpix$DrainingArea.km2.LDD>100),]
df=matp100
p<-ggplot(df, aes(x=DrainingArea.km2.LDD)) + 
  geom_histogram(color="steelblue", fill="slategray1",bins=15,alpha=0.9,lwd=1)+
  scale_y_continuous(name="Number of stations")+
  scale_x_log10(name="Upstream area (km2)",
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
  annotate("label", x=250000, y=140, label= paste0("n = ",length(df$upAc)),size=6)

p
ggsave("histo_stations.jpg", p, width=20, height=15, units=c("cm"),dpi=1500)





#################OLD shitty station detection########################
shittystations=dat_stat[which(dat_stat$rat1>5 | dat_stat$rat2>5),]

#Remove all strange 
shitty_ori=inner_join(shittystations,dat,by=c("station"="Station_ID"))



shitty_rmv=shitty_ori[which(shitty_ori$qm>10),]
shitty_rmv$uniki=paste0(shitty_rmv$station,shitty_rmv$year)
dat$uniki=paste0(dat$Station_ID,dat$year)
rmv=(match(shitty_rmv$uniki,dat$uniki))
do=dat[rmv,]
dat=dat[-rmv,]



#load stations from EFAS to check if there are faulty matches

efas_stations=read.csv(paste0(valid_path,"hres_calib_stations.csv"))
matefas=match(flefas$StationID,shittystations$station)
#First I join with my matching stations to have all fieds
ptn=inner_join(flefas,efas_stations,by=c("National_Station_Identifier"))

shittystations_dat=inner_join(ptn,shittystations, by=c("StationID"="station"))
shittystations_keep=shittystations_dat[which(!is.na(shittystations_dat$Cal_area)),]

ggplot(shittystations_dat) + geom_point(aes(x=obs, y=sim,col=factor(Station_ID)),alpha=0.3) +
  scale_x_log10(name="observed",
                breaks=c(1,10,100,1000), minor_breaks = log10_minor_break(),
                labels=c("1","10","100","1000"), limits=c(0.001,tmax)) +
  scale_y_log10(name="simulated",
                breaks=c(1,10,100,1000), minor_breaks = log10_minor_break(),
                labels=c("1","10","100","1000"), limits=c(0.001,tmax)) +
  geom_abline(slope=1, intercept=0, lwd =1, alpha=0.5)+
  # scale_x_continuous(name=nobs,breaks=seq(0,tmax, by=bi))+
  # scale_y_continuous(name=nsim,breaks=seq(0,tmax, by=bi))+
  annotate("label", x=800, y=max(dat$sim), label= paste0("R2 = ",r2),size=5)+
  scale_color_discrete()+
  theme(axis.title=element_text(size=16, face="bold"),
        axis.text = element_text(size=16),
        panel.background = element_rect(fill = "white", colour = "white"),
        panel.grid = element_blank(),
        panel.border = element_rect(linetype = "solid", fill = NA, colour="black"),
        legend.title = element_text(size=14),
        legend.text = element_text(size=12),
        legend.position = "none",
        panel.grid.major = element_line(colour = "grey80"),
        panel.grid.minor = element_line(colour = "grey90",linetype="dashed"),
        legend.key = element_rect(fill = "transparent", colour = "transparent"),
        legend.key.size = unit(.8, "cm"))


#Remove all strange 
shitty_ori=inner_join(shittystations,dat,by=c("station"="Station_ID"))
shitty_rmv=shitty_ori[which(shitty_ori$qm>10),]
shitty_rmv$uniki=paste0(shitty_rmv$station,shitty_rmv$year)
dat$uniki=paste0(dat$Station_ID,dat$year)
rmv=(match(shitty_rmv$uniki,dat$uniki))
do=dat[rmv,]
dat=dat[-rmv,]




