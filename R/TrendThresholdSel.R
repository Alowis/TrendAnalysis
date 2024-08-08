#Code to identify trend threshold for each catchment, 
#if threshold is different, I need a flag and I check manually what is the best





main_path = "D:/tilloal/Documents/06_Floodrivers/"
valid_path = paste0(main_path)

dis_path<-paste0(main_path,'HERA_CFS/')
load(file=paste0(dis_path,"HERA_CFS6h_H07_19502020.Rdata"))
Station_data_IDs <- as.vector(t(Q_sim[1, ]))[-1]
tqx=as.POSIXct((Q_sim[,1]*3600*24), origin="1979-01-01 00:00:00")
txx=tqx[-1]
Q_simCFS=Q_sim

dis_path<-paste0(main_path,'HERA/')
load(file=paste0(dis_path,"HERA6h_H07_19502020.Rdata"))
Q_simH=Q_sim

dis_path<-paste0(main_path,'HERA_Rstat/')
load(file=paste0(dis_path,"HERA_CFR6h_H07_19502020.Rdata"))
Q_simCFR=Q_sim

#extract values for 1 location

#remove locations without discharge
wolala=which(is.na(Q_simCFR[2,]))-1
Q_simH=Q_simH[,-which(is.na(Q_simH[2,]))]
Q_simCFR=Q_simCFR[,-which(is.na(Q_simCFR[2,]))]
Q_simCFS=Q_simCFS[,-which(is.na(Q_simCFS[2,]))]
Station_data_IDs=Station_data_IDs[-wolala]

#match HYBAS_id with station_IDs

hybmatch=match(Station_data_IDs,cst7$outlets)

hybas_IDs=cst7$HYBAS_ID[hybmatch]
txx=tqx[-1]
rmv=which(year(txx)==1950)
txx=txx[-rmv]

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

series=timeAndSeriesH$data
timeStamps=timeAndSeriesH$timestamp

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
    # varianceSeries <- tsEvaNanRunningVariance(serieb, rs@nRunMn)
    # varianceSeries <- tsEvaNanRunningMean(varianceSeries, 
    #                                       ceiling(rs@nRunMn/2))
    # norm_trend <- rs@trendSeries/mean(rs@trendSeries, na.rm = TRUE)
    dtr1 = serieb - rs@trendSeries
    lneg = length(which(dtr1 < 0))
    # stab <- cor(nr, norm_trend, use = "pairwise.complete.obs")
    if (iter == 1) 
      nr <- norm_trend
    if (lneg >= 1) 
      lnegs = c(lnegs, lneg)
      # sts <- c(sts, stab)
      pctd = c(pctd, pcts[iter])
    }

  rval = pctd[length(pctd)]
  if (sum(lnegs) > 1) {
    rval = pctd[which.min(lnegs)]
  }
  
  return(rval)
}

Trth_H_list=c()
Trth_R_list=c()
Trth_S_list=c()

for (id in 1:length(Station_data_IDs)){
  
  print(id)
  stid=Station_data_IDs[id]
  catch=stid
  Q_H=Q_simH[,id+1][-1]
  Q_SCF=Q_simCFS[,id+1][-1]
  Q_RCF=Q_simCFR[,id+1][-1]

  Q_H=Q_H[-rmv]
  Q_SCF=Q_SCF[-rmv]
  
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
  
  TrendTh_H1=tsEvaFindTrendThreshold(series=timeAndSeriesH$data, timeStamps=timeAndSeriesH$timestamp, timeWindow)
  TrendTh_R1=tsEvaFindTrendThreshold(series=timeAndSeriesR$data, timeStamps=timeAndSeriesR$timestamp, timeWindow)
  TrendTh_S1=tsEvaFindTrendThreshold(series=timeAndSeriesS$data, timeStamps=timeAndSeriesS$timestamp, timeWindow)
  
  TrendTh_H2=tsEvaFindTrendThreshold2(series=timeAndSeriesH$data, timeStamps=timeAndSeriesH$timestamp, timeWindow)
  TrendTh_R2=tsEvaFindTrendThreshold2(series=timeAndSeriesR$data, timeStamps=timeAndSeriesR$timestamp, timeWindow)
  TrendTh_S2=tsEvaFindTrendThreshold2(series=timeAndSeriesS$data, timeStamps=timeAndSeriesS$timestamp, timeWindow)
  
  Trth_H=c(catch,TrendTh_H1,TrendTh_H2)
  Trth_R=c(catch,TrendTh_R1,TrendTh_R2)
  Trth_S=c(catch,TrendTh_S1,TrendTh_S2)
  
  Trth_H_list=rbind(Trth_H_list,Trth_H)
  Trth_R_list=rbind(Trth_R_list,Trth_R)
  Trth_S_list=rbind(Trth_S_list,Trth_S)
  
}

Trth_H_list-Trth_S_list
