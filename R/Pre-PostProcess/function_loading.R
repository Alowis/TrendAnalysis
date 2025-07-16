#Trend functions script

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
disNcopen=function(fname,dir,outloc){
  ncdis=paste0(dir,"/",fname,".nc")
  ncd=nc_open(ncdis)
  name.vb=names(ncd[['var']])
  namev=name.vb[1]
  time <- ncvar_get(ncd,"time")
  lt=length(time)
  
  #timestamp corretion
  name.lon="lon"
  name.lat="lat"
  londat = ncvar_get(ncd,name.lon) 
  llo=length(londat)
  latdat = ncvar_get(ncd,name.lat) 
  lla=length(latdat)
  
  outllplus=matrix(-9999, nrow = lt*length(outloc[,1]), ncol = 5)
  outllplus=as.data.frame(outllplus)
  for ( idc in 1:length(outloc[,1])){
    print(idc)
    idm=1+(idc-1)*lt
    start=c(outloc$idlo[idc],outloc$idla[idc],1)
    count=c(1,1,lt)
    outlets = ncvar_get(ncd,namev,start = start, count= count)
    outlets=as.vector(outlets)
    outid=rep(outloc[idc,1],length(time))
    #lonlatloop=expand.grid(c(1:lla),c(1:lt))
    lon=rep(londat[idc],length(time))
    lat=rep(latdat[idc],length(time))
    outll=data.frame(outlets,outid,lon,lat,time)
    names(outllplus)=names(outll)
    fck=length(c(idm:(idm+lt-1)))
    outllplus[c(idm:(idm+lt-1)),]= outll
  }
  
  
  # outllplus=c()
  # for ( idc in 1:llo){
  #   print(idc)
  #   start=c(idc,1,1)
  #   count=c(1,lla,1)
  #   outlets = ncvar_get(ncd,namev,start = start, count= count)
  #   lna=length(which(!is.na(outlets)))
  #   if (lna>0){
  #     print("data")
  #     exact=(which(!is.na(outlets)))
  #     for (st in exact){
  #       start=c(idc,st,1)
  #       count=c(1,1,lt)
  #       outlets = ncvar_get(ncd,namev,start = start, count= count)
  #       outlets=as.vector(outlets)
  #       #lonlatloop=expand.grid(c(1:lla),c(1:lt))
  #       lon=rep(londat[idc],length(time))
  #       lat=rep(latdat[st],length(time))
  #       outll=data.frame(outlets,lon,lat,time)
  #     }
  #     outllplus=rbind(outllplus,outll)
  #   }
  
  #}
  
  return (outllplus)
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


ComputeReturnLevels<-function(nonStationaryEvaParams, RPgoal, timeIndex,trans=NA){
  
  #GEV
  epsilonGEV <- nonStationaryEvaParams[[1]]$parameters$epsilon
  sigmaGEV <- mean(nonStationaryEvaParams[[1]]$parameters$sigma[timeIndex])
  muGEV <- mean(nonStationaryEvaParams[[1]]$parameters$mu[timeIndex])
  dtSampleYears <- nonStationaryEvaParams[[1]]$parameters$timeDeltaYears
  
  #GPD
  epsilonGPD <- nonStationaryEvaParams[[2]]$parameters$epsilon
  sigmaGPD <- mean(nonStationaryEvaParams[[2]]$parameters$sigma[timeIndex])
  thresholdGPD <- mean(nonStationaryEvaParams[[2]]$parameters$threshold[timeIndex])
  nPeaks <- nonStationaryEvaParams[[2]]$parameters$nPeaks
  thStart <- nonStationaryEvaParams[[2]]$parameters$timeHorizonStart
  thEnd <- nonStationaryEvaParams[[2]]$parameters$timeHorizonEnd
  sampleTimeHorizon <- as.numeric((thEnd - thStart)/365.2425)
  
  if (nonStationaryEvaParams[[1]]$method=="No fit"){
    print("could not fit EVD to this pixel")
    ParamGEV=c(epsilonGEV,sigmaGEV,muGEV,NA, NA, NA)
    names(ParamGEV)=c("epsilonGEV","sigmaGEV","muGEV","epsilonStdErrGEV","sigmaStdErrGEV","muStdErrGEV")
    
    ParamGPD=c(epsilonGPD,sigmaGPD,thresholdGPD,NA,NA, NA,nPeaks,sampleTimeHorizon)
    names(ParamGPD)=c("epsilonGPD","sigmaGPD","thresholdGPD","epsilonStdErrGPD","sigmaStdErrGPD","thresholdStdErrGPD","nPeaks","SampleTimeHorizon")
    return(list(Fit="No fit",Params=c(ParamGEV,ParamGPD)))
  }else{
    #GEV
    # epsilonGEV <- nonStationaryEvaParams[[1]]$parameters$epsilon
    # sigmaGEV <- mean(nonStationaryEvaParams[[1]]$parameters$sigma[timeIndex])
    # muGEV <- mean(nonStationaryEvaParams[[1]]$parameters$mu[timeIndex])
    # dtSampleYears <- nonStationaryEvaParams[[1]]$parameters$timeDeltaYears
    epsilonStdErrGEV <- nonStationaryEvaParams[[1]]$paramErr$epsilonErr
    sigmaStdErrGEV <- mean(nonStationaryEvaParams[[1]]$paramErr$sigmaErr[timeIndex])
    muStdErrGEV <- mean(nonStationaryEvaParams[[1]]$paramErr$muErr[timeIndex])
    
    #GPD
    # epsilonGPD <- nonStationaryEvaParams[[2]]$parameters$epsilon
    # sigmaGPD <- mean(nonStationaryEvaParams[[2]]$parameters$sigma[timeIndex])
    # thresholdGPD <- mean(nonStationaryEvaParams[[2]]$parameters$threshold[timeIndex])
    # nPeaks <- nonStationaryEvaParams[[2]]$parameters$nPeaks
    epsilonStdErrGPD <- nonStationaryEvaParams[[2]]$paramErr$epsilonErr
    sigmaStdErrGPD <- mean(nonStationaryEvaParams[[2]]$paramErr$sigmaErr[timeIndex])
    thresholdStdErrGPD <- mean(nonStationaryEvaParams[[2]]$paramErr$thresholdErr[timeIndex])
    # thStart <- nonStationaryEvaParams[[2]]$parameters$timeHorizonStart
    # thEnd <- nonStationaryEvaParams[[2]]$parameters$timeHorizonEnd
    # sampleTimeHorizon <- as.numeric((thEnd - thStart)/365.2425)
    if (trans=="rev"){
      sigmaGEV=-sigmaGEV
      muGEV=-muGEV
      sigmaGPD=-sigmaGPD
      thresholdGPD=-thresholdGPD
    }
    returnLevelsGEV <- tsEvaComputeReturnLevelsGEV(epsilonGEV, sigmaGEV, muGEV, epsilonStdErrGEV, sigmaStdErrGEV, muStdErrGEV, RPgoal)
    
    returnLevelsGPD <- tsEvaComputeReturnLevelsGPD(epsilonGPD, sigmaGPD, thresholdGPD, epsilonStdErrGPD, sigmaStdErrGPD, thresholdStdErrGPD,
                                                   nPeaks, sampleTimeHorizon, RPgoal)
    rlevGEV=returnLevelsGEV$returnLevels
    rlevGPD=returnLevelsGPD$returnLevels
    
    errGEV=returnLevelsGEV$returnLevelsErr
    errGPD=returnLevelsGPD$returnLevelsErr
    
    ParamGEV=c(epsilonGEV,sigmaGEV,muGEV,epsilonStdErrGEV, sigmaStdErrGEV, muStdErrGEV)
    names(ParamGEV)=c("epsilonGEV","sigmaGEV","muGEV","epsilonStdErrGEV","sigmaStdErrGEV","muStdErrGEV")
    
    ParamGPD=c(epsilonGPD,sigmaGPD,thresholdGPD,epsilonStdErrGPD,sigmaStdErrGPD, thresholdStdErrGPD,nPeaks,sampleTimeHorizon)
    names(ParamGPD)=c("epsilonGPD","sigmaGPD","thresholdGPD","epsilonStdErrGPD","sigmaStdErrGPD","thresholdStdErrGPD","nPeaks","SampleTimeHorizon")
    return(list(Fit="Fitted",ReturnLevels=c(ReturnPeriod=RPgoal, GEV=as.numeric(rlevGEV),GPD=as.numeric(rlevGPD),errGEV=as.numeric(errGEV),errGPD=as.numeric(errGPD)),Params=c(ParamGEV,ParamGPD)))
  }  
  
  
}
RPcalc<- function(params,RPiGEV,RPiGPD){
  paramx=data.frame(t(params))
  qxV=1-exp(-(1+paramx$epsilonGEV*(RPiGEV-paramx$muGEV)/paramx$sigmaGEV)^(-1/paramx$epsilonGEV))
  returnPeriodGEV=1/qxV
  X0 <- paramx$nPeaks/paramx$SampleTimeHorizon
  qxD=(((1+paramx$epsilonGPD*(RPiGPD-paramx$thresholdGPD)/paramx$sigmaGPD)^(-1/paramx$epsilonGPD)))
  returnPeriodGPD=1/(X0*qxD)
  
  return(c(GEV=returnPeriodGEV,GPD=returnPeriodGPD))
}
plotchangemaps=function(basemap,catmap,datar,law="GPD",type,period=c(1950,1990),parlist,valuenames){
  
  
  data=full_join(catmap,datar,by = c("outlets"="unikout"))
  datacol=names(data)
  valcol=match(valuenames,datacol)
  datatwin=data
  valcol2=valcol-1
  st_geometry(datatwin) <- NULL
  datatwin=as.data.frame(datatwin)
  class(datatwin)
  dtc=names(datatwin)
  mcor=unique(match(datar$unikout,parlist$catchment))
  if (length(which(is.na(mcor)))>0) datar=datar[-which(is.na(mcor)),]
  
  cref=paste0("Y",period[1])
  crefloc=match(cref,dtc)
  finalperiod=paste0("Y",period[2])
  colsel=match(finalperiod,datacol)
  
  if (type=="RLchange"){
    palet=c(hcl.colors(11, palette = "BrBG", alpha = NULL, rev = F, fixup = TRUE))
    title=paste0("Relative change in 100-years flood Return Level between ", period[1], " and ", period[2])
    legend="Relative Change (%)"
    tmpval=(datatwin[,valcol2])/(datatwin[,crefloc])*100-100
    for (it in 1:length(tmpval[1,])){
      tmpval[,it][which(tmpval[,it]>80)]=NA
    }
    data[,valcol]=tmpval
    br=c(-100,-80,-60,-40,-20,0,20,40,60,80,100)
    limi=c(-100,100)
    trans=scales::modulus_trans(.8)
    colNA="darkgrey"
  }
  if (type=="RPchange"){
    palet=c(hcl.colors(9, palette = "BrBG", alpha = NULL, rev = T, fixup = TRUE))
    title=paste0(period[2]," return period of a 100-years Return Level flood in ", period[1])
    legend="Return Period (Years)"
    br=c(1,10,100,1000,10000)
    limi=c(1,10000)
    trans="log"
    newRP=RPchangeCal(parlist, yi=period[1], yf=period[2], RetLev=datar,law, valuenames)
    ow=match(data$outlets,newRP$catchment)
    data[,colsel]=newRP$newRP[ow]
    colNA="black"
    
  }
  
  
  names(data)[colsel]="fillplot"
  
  ggplot(basemap) +
    geom_sf(fill="white")+
    geom_sf(data=data,aes(fill=fillplot,geometry=geometry),color="transparent")+ 
    geom_sf(fill=NA, color="grey") +
    library(sf)
  library(ggplot2)
  
  coord_sf(xlim = c(min(nco[,1]),max(nco[,1])), ylim = c(min(nco[,2]),max(nco[,2]))) +
    # coord_sf(xlim = c(min(nco[,1]),max(nco[,1])), ylim = c(min(nco[,2]),max(nco[,2])))+
    # geom_point(data=event,aes(x=lon,y=lat,col=dr))+
    # scale_color_distiller(palette = "Spectral",
    #                       direction = -1, limit=c(0,max(eventp$dr)),oob = scales::squish)+
    scale_fill_gradientn(
      colors=palet,
      breaks=br,limits=limi,trans=trans,
      oob = scales::squish,na.value=colNA, name=legend)   +
    labs(x="Longitude", y = "Latitude") +
    guides(fill = guide_colourbar(barwidth = 15, barheight = .8)) +
    theme(axis.title=element_text(size=tsize),
          panel.background = element_rect(fill = "aliceblue", colour = "grey1"),
          panel.border = element_rect(linetype = "solid", fill = NA, colour="black"),
          legend.title = element_text(size=tsize),
          legend.text = element_text(size=osize),
          legend.position = "bottom",
          panel.grid.major = element_line(colour = "grey70"),
          panel.grid.minor = element_line(colour = "grey90"),
          legend.key = element_rect(fill = "transparent", colour = "transparent"),
          legend.key.size = unit(.8, "cm")) +
    ggtitle(title)
  
}
RPchangeCal=function(parlist, yi, yf, RetLev,law,valuenames){
  #Use TSEVA function for return level computation
  #1 selection of the year I want
  #check if the computation can be done on all the points at the same time
  #ay ay ay
  cname=names(RetLev)
  cna=match(valuenames,cname)
  parlistinit=parlist[which(parlist$Year==yi),]
  parlistf=parlist[which(parlist$Year==yf),]
  cref=paste0("Y",yi)
  crefloc=match(cref,cname)
  RLi=RetLev[,crefloc]
  
  
  if (law=="GEV"){
    paramx=data.frame((parlistf))
    qxV=1-exp(-(1+paramx$epsilonGEV*(RLi-paramx$muGEV)/paramx$sigmaGEV)^(-1/paramx$epsilonGEV))
    returnPeriods=1/qxV
  } else if (law=="GPD"){
    paramx=data.frame((parlistf))
    X0 <- paramx$nPeaks/paramx$SampleTimeHorizon
    qxD=(((1+paramx$epsilonGPD*(RLi-paramx$thresholdGPD)/paramx$sigmaGPD)^(-1/paramx$epsilonGPD)))
    returnPeriods=1/(X0*qxD)
  }
  min(returnPeriods)
  
  return(data.frame(catchment=paramx$catchment,newRP=returnPeriods))
}
plotRegime <- function(data, catch, seuil=NULL){
  #~ data <- as.simu(data)
  names(data)=c("date","Q")
  if( is.null(seuil) ) seuil <- -Inf
  ## Suppression des lignes sans couples Qobs/Qsim
  data <- subset( data, !is.na(Q) & Q > seuil)
  ## Preparation de la sequence de mois
  deb=data$date[1]
  mois.deb <-  seq(as.Date(deb), by="month", length=12)
  
  period=paste0(format(range(data$date)[1],"%b %Y"),"-",format(range(data$date)[2],"%b %Y"))
  main = bquote(.(catch)~ " Regime | Period: "~ .(period))
  jours <- as.numeric(format(data$date,"%j"))
  ## Iddinces des lignes de meme categorie
  ind.j <- tapply(seq(length(jours)), jours, c)
  ind.j <- ind.j[-366]
  ## Calcul
  Qc <- data.frame(date=as.numeric(names(ind.j)),
                   mean=sapply(ind.j, function(x) mean(data$Q[x], na.rm=TRUE)),
                   q10=sapply(ind.j, function(x) quantile(data$Q[x],0.1, na.rm=TRUE)),
                   q90=sapply(ind.j, function(x) quantile(data$Q[x],.9, na.rm=TRUE)))
  
  qlim=c(0,1.2*max(Qc[,4]))
  plot(Qc$date, Qc$mean, type="n", axes=FALSE, ylim=qlim,
       xlab = NA, ylab = expression(paste("Debit (",m^3/s,")")))
  mtext(main,3,font = 2,line = 0.5,cex = .75)
  lines(Qc$date, Qc$mean, col="darkblue",lwd=2)
  lines(Qc$date, Qc$q10, col="blue",lwd=1,lty=2)
  lines(Qc$date, Qc$q90, col="blue",lwd=1,lty=2)
  abline(v = format(mois.deb,"%j"), col="lightgrey", lty=3)
  abline(h = quantile(data$Q,.95),col="red",lwd=2)
  polygon(c(Qc$date,rev(Qc$date)),c(Qc$q10,rev(Qc$q90)),
          col=alpha("lightblue",.4),border="transparent")
  axis(2)
  axis(1, format(mois.deb,"%j"), label=format(mois.deb,"%b"),cex.axis=.8)
  box()
  legend("topleft", leg=c("mean","Q90","Q10"),
         lwd=c(2,2,1,2), col=c("darkblue","blue","blue"),
         cex=.8, lty=c(1,2,2),bty="n")
  return(Qc)
}
#the transformation has to happen later
interid= function(data,WindowSize) {
  
  dt1=min(diff(data$date),na.rm=T)
  dt=as.numeric(dt1)
  tdim=attributes(dt1)$units
  if (tdim=="hours") dt=dt/24
  nRunMn = ceiling(WindowSize/dt);
  colnames(data)
  data$Q7=rollmean(data$Qs,nRunMn, align = "right", fill=NA)
  
  data$Qrev=-(data$Q7)
  data$Qinv=1/data$Q7
  data$Qlninv=-log(data$Q7)
  if (length(which(is.na(data$Q7)))==length(data$Q7)){
    print("no data in this pixel")
    list0=NA
    dis07=data
    l0=NA
    fl=NA
    mindis=NA
    
  }else{
    if (length(which(is.na(data$Qs)))>0){
      print("Na alert")
      #seriefill=tsEvaFillSeries(data$date,data$Qs)
      #data$Qs=seriefill
    }
    #tinversing the discharge for EVA
    #Identification of the intermittent river
    mindis=min(data$Q7,na.rm=T)
    m0=length(which(data$Qs7==mindis))
    l0=length(which(data$Q7<=1e-4))
    #30y running average on this
    dayysbelow=tsEvaNanRunnigBlowTh(data$Q7,1e-4,4*365*30)
    dayysbelow$time=as.Date(data$date[dayysbelow$time])
    yrtot=length(unique(year(data$date)))
    list0=NA
    fl=0 #flag for intermitent river
    if (l0>=1){
      print("intermittent river")
      fl=1
      data$Qinv=NA
      data$Qlninv
      #keep days with 0 discharge
      list0=data$date[which(data$Qs==0)]
    }else if(mindis>0 & m0>=yrtot){ 
      print(paste0("river with floor low flow ",mindis))
      fl=3
      list0=data$date[which(data$Q7==mindis)]
    }
    dis07=data[,c(2,3,4,5,6)]
  }
  #The objective here is to return also the dicharge in a better format for next step of the analysis
  return(list(zerodate=list0,trdis=dis07,DaysBlow=dayysbelow,flags=c(n0d=l0,intertype=fl,mindischarge=mindis)))
}



max_daily_value <- function(timeseries) {
  # Convert the 6-hourly time series to daily time series
  daily_timeseries <- apply.daily(timeseries, max)
  tday=unique(as.Date(timeseries$date))
  return(data.frame(date=tday,Qmd=daily_timeseries))
}
min_daily_value <- function(timeseries) {
  # Convert the 6-hourly time series to daily time series
  daily_timeseries <- apply.daily(timeseries, min)
  tday=unique(as.Date(timeseries[,1]))
  return(data.frame(date=tday,Qmd=daily_timeseries))
}