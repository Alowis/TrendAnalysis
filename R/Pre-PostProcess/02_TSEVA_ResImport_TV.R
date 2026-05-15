


setwd(dirname(rstudioapi::getSourceEditorContext()$path))
getwd()
setwd("D:/tilloal/Documents/LFRuns_utils/TrendAnalysis/R/")
source("functions_trends.R")
hydroDir<-("D:/tilloal/Documents/LFRuns_utils/data")


# Load inputs from HPC computation ----------------------------------------



hazard="Drought"

rspace= read.csv(paste0(hydroDir,"/subspace_efas.csv"))
rspace=rspace[,-1]

dataDir=paste0("D:/tilloal/Documents/LFRuns_utils/data/",hazard,"/HPC/Calibrated/revision/TrendVar/")

lsce=c("HistoX")
#lsce=c("RWStatX")

for (sce in lsce){
  lf=list.files(path = paste0(dataDir,"/",sce,"/"), full.names = TRUE, recursive = TRUE)
  lf
  nfiles=length(lf)
  cal=T
  if (nfiles<175){
    print("files missing")
  }
  fi1=lf[1]
  load(fi1)
  time<-Results$Trend[,1]
  #wow=Results$Trend[,2]
  
  Yearlist=unique(year(time))
  #I will retain only one value per year which is the last day of the year
  datex=yday(time)
  dtect=c(diff(datex),-1)
  last_days <- time[which(dtect<0)]
  tindexes=match(last_days,time)
  
  # 1. Setup metadata and pre-allocate
  n_files <- length(lf)
  rspace  <- read.csv(file.path(hydroDir, "subspace_efas.csv"))[,-1]
  
  # Pre-allocate empty matrices (filled with NAs)
  # This reserves the memory block once and for all
  TrendSave       <- Results$Trend[tindexes,1]
  VariabilitySave <- Results$Variability[tindexes,1]
  
  # 2. Loop and fill by index
  for (i in seq_along(lf)) {
    load(lf[i])
    print(i)
    
    # Filename parsing
    fils <- basename(lf[i])
    tt   <- unlist(strsplit(fils, "_"))
    Nsq_raw <- as.numeric(tt[3])
    
    # Logic for Nsq
    Nsq <- if (cal || Nsq_raw > 88) floor(Nsq_raw / 10) else Nsq_raw
    
    # Fill the pre-allocated column directly
    TrendSave      <- cbind(TrendSave,Results$Trend[tindexes,-1])
    VariabilitySave <- cbind(VariabilitySave,Results$Variability[tindexes,-1])
    
    # Explicitly remove 'Results' to free up space for the next iteration
    rm(Results)
  }
  
  # 3. Clean up the global environment
  gc()
  
  ResTV=list(TrendSave,VariabilitySave)
  
  
  
  #print(paste0(hydroDir,"/",hazard,"/RL100.",hazard,".",season,".",scenario,mmx,".Rdata"))
  # Saving outputs of the loop
  saveout=T
  if (saveout==T){
    save(ResTV,file=paste0(dataDir,"/",sce,"/Var_TrendX_agg.Rdata"))
  }
  
  gc()
  
}
