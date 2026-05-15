
#Compre changes in 10 and 100yr RLs
hydroDir<-("D:/tilloal/Documents/LFRuns_utils/data")

haz="Flood"
load(file=paste0(hydroDir,"/TSEVA/output_plots/",haz,"_pixChange_RL100_v1.Rdata"))
Data100=DataSave
load(file=paste0(hydroDir,"/TSEVA/output_plots/",haz,"_pixChange_v3.Rdata"))
Data10=DataSave

names(Data100)
names(Data10) 
  for (dri in 2:5){
  driver=names(Data100)[[dri]]
  DataC100=Data100[[dri]]
  DataC10=Data10[[dri]]
  
  hist(DataC10$Y2015, xlim=c(-5,5),breaks=1000)
  
  ###5.5.1 Aggregation of relative changes ----
  pointagg <- aggregate(list(Rchange_rel = DataC100$Y2015),
                        by = list(HydroR = DataR$HER),
                        FUN = function(x) c(mean = mean(x, na.rm = TRUE), 
                                            dev = sd(x, na.rm = TRUE), 
                                            len = length(x), 
                                            med = median(x, na.rm = TRUE), 
                                            q1 = quantile(x, 0.025, na.rm = TRUE), 
                                            q3 = quantile(x, 0.975, na.rm = TRUE)))
  point100 <- do.call(data.frame, pointagg)
  
  
  pointagg <- aggregate(list(Rchange_rel = DataC10$Y2015),
                        by = list(HydroR = DataR$HER),
                        FUN = function(x) c(mean = mean(x, na.rm = TRUE), 
                                            dev = sd(x, na.rm = TRUE), 
                                            len = length(x), 
                                            med = median(x, na.rm = TRUE), 
                                            q1 = quantile(x, 0.025, na.rm = TRUE), 
                                            q3 = quantile(x, 0.975, na.rm = TRUE)))
  point10 <- do.call(data.frame, pointagg)
  
  point1000=data.frame(rl10=point10$Rchange_rel.mean,rl100=point100$Rchange_rel.mean)
  data=data.frame(rl100=DataC100$Y2015,rl10=DataC10$Y2015)
  data=data[-which(is.na(data$rl100)),]
  
  cor(point1000)
  tmin=(quantile(data$rl100,.0001))
  tmax=(quantile(data$rl100,.9999))
  #do a nice scatterplot
  
  #data$density <- get_density((data$rl100), (data$rl10), n = 500)
  #dat$density <- get_density(log(dat$obsj), log(dat$simj), n = 100)
  scp=ggplot() + 
    # geom_point(aes(x=obs, y=sim,col=density),stroke=0,size=3,alpha=0.2,shape=16) +
    geom_point(data=data,aes(x=rl10, y=rl100),col="grey",stroke=0,size=3,alpha=0.1,shape=16) +
    geom_point(data=point1000,aes(x=rl10, y=rl100),col="black",stroke=0,size=4,alpha=.7,shape=16) +
    geom_smooth(data=data, aes(x=rl10, y=rl100), method = "lm", color = "red", se = FALSE, lwd = 1)+
    #geom_point(data=data,aes(x=rl10, y=rl100),col="grey",stroke=0,size=3,alpha=0.1,shape=16) +
    # geom_jitter(aes(x=obs, y=sim,col=density),stroke=0,size=3,alpha=0.5,shape=16, height=.1, width=.0)+
    geom_abline(slope=1, intercept=0, lwd =1, alpha=1,col="darkgreen",linetype="dashed")+
    scale_x_continuous(limits=c(tmin,tmax), trans=scales::modulus_trans(0.5)) +
    scale_y_continuous(limits=c(tmin,tmax),trans=scales::modulus_trans(0.5)) +
    # scale_y_log10(name=nsim,
    #               breaks=c(0.1,1,10,100,1000,10000), minor_breaks = log10_minor_break(),
    #               labels=c("0.1","1","10","100","1000","10000"), limits=c(tmin,tmax)) +
    # scale_x_continuous(name=nobs,breaks=seq(0,tmax, by=bi))+
    # scale_y_continuous(name=nsim,breaks=seq(0,tmax, by=bi))+
    #annotate("label", x=3*tmin, y=max(tmax), label= paste0("R2 = ",r2),size=5)+
    labs(x="changes in RL10 (%)", y = "Changes in RL100 (%)")+
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
  
  ggsave(paste0("D:/tilloal/Documents/01_Projects/FlooDrough/Figures/Revisions/",haz,"_scatter_rls_",driver,".jpg"),scp, width=20, height=15, units=c("cm"),dpi=1000)
  
  mloc=which(!is.na(match(DataC100$llcoord,DataC10$llcoord)))
  DataC10=DataC10[mloc,]
  epsilon=as.numeric(sign(DataC10$Y2015)*0.1)
  DataC10$diff=DataC100$Y2015/(DataC10$Y2015+epsilon)
  
  weirdos<-cbind(DataC10$Y2015[which(DataC10$diff<0)],DataC100$Y2015[which(DataC10$diff<0)])
  weirdos=data.frame(weirdos)
  
  
  pointagg <- aggregate(list(Rchange_rel = DataC10$diff),
                        by = list(HydroR = DataC10$HER),
                        FUN = function(x) c(mean = mean(x, na.rm = TRUE), 
                                            dev = sd(x, na.rm = TRUE), 
                                            len = length(x), 
                                            med = median(x, na.rm = TRUE), 
                                            q1 = quantile(x, 0.025, na.rm = TRUE), 
                                            q3 = quantile(x, 0.975, na.rm = TRUE)))
  pointClim <- do.call(data.frame, pointagg)
  
  
  
  DataC10=DataC10[-which(is.na(DataC10$Var1)),]
  points <- st_as_sf(DataC10, coords = c("Var1", "Var2"), crs = 4326)
  points <- st_transform(points, crs = 3035)
  
  Regio=HydroRsf
  
  if (haz=="Flood"){
    br=seq(-2,2,by=0.5)
    labels=br
    limi=c(-2,2)
    tsize=16
    osize=12
    legend2="Change(rl100)/Change(rl10)"
    palet=c(hcl.colors(11, palette = "RdYlBu", alpha = NULL, rev = F, fixup = TRUE))
    paletf=c(hcl.colors(11, palette = "RdBu", alpha = NULL, rev = F, fixup = TRUE))
    # Merge with NUTS3 data to get spatial points
    pag <- inner_join(Regio, pointClim, by = c("CODEB" = "HydroR"))
  
  
      titleX=paste0("Ratio between 10-year and 100-year RL changes for ",haz," attributed \nto ",driver," -  1955-2015")
      fmap<-ggplot(basemap) +
        geom_sf(fill="white",color="darkgrey",size=0.5)+
        # geom_sf(data=pag,aes(fill=Rchange_rel.mean,geometry=geometry),alpha=0.2,color="transparent")+
        geom_sf(data=points,aes(col=diff,geometry=geometry,size=upa),alpha=.9,stroke=0,shape=15)+ 
        
        scale_size(range = c(0.08, 0.4), trans="sqrt",name= expression(paste("Upstream area ", (km^2),
                                                                             sep = " ")),
                   breaks=c(101,1000,10000,100000,500000), labels=c("100","1000", "10 000", "100 000", "500 000"),
                   guide = "none")+
        # scale_fill_gradientn(
        #   colors=paletf,
        #   breaks=br,limits=limi,
        #   oob = scales::squish,na.value=colNA, name=legend2)   +
        guides(fill = "none")+
        coord_sf(xlim = c(min(nco[,1]),max(nco[,1])), ylim = c(min(nco[,2]),max(nco[,2])))+
        scale_color_gradientn(
          colors=palet,
          breaks=br,limits=limi,
          oob = scales::squish, name=legend2)   +
        guides(fill = "none")+
        labs(x="Longitude", y = "Latitude")+
        guides(colour = guide_colourbar(barwidth = 22, barheight = 1))+
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
      
      
      ggsave(paste0("D:/tilloal/Documents/01_Projects/FlooDrough/Figures/Revisions/mapF_",
                    driver,"_",haz,"RLvs1.jpg"), fmap,
             width=22, height=20, units=c("cm"),dpi=800) 
      
      
      DataC10$col="no"
      DataC10$col[which(DataC10$diff<0)]="yes"
      n1=round(length(which(DataC10$col=="yes"))/length(DataC10$col),3)*100
      p<-ggplot(DataC10, aes(x=diff, fill=col)) + 
        geom_histogram(color="gray",breaks=seq(-5,5,by=0.2),alpha=0.9,lwd=1)+
        scale_y_continuous(breaks=seq(0,200000, by=10000),name="Number of pixels")+
        scale_x_continuous(breaks=seq(-5,5, by=1),name="ratio")+
        scale_fill_manual(values=c("no"="slategray1","yes"="lightcoral"),name="change of sign")+
        
        coord_cartesian(xlim = c(-5,5))+
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
              legend.key.size = unit(.8, "cm"))+
        annotate("label", x=-2, y=100000, label= paste0("sign change = ",n1,"%"),size=6)
      
      p
      ggsave(paste0("D:/tilloal/Documents/01_Projects/FlooDrough/Figures/Revisions/,",haz,"_histo_ratio_",driver,".jpg"), p, width=20, height=15, units=c("cm"),dpi=1500)
      
    
      
  }
  }
