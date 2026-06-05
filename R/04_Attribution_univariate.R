## ============================================================
##  04_Attribution_univariate.R
##  Reorganized: main figures first, supplementary below a flag
## ============================================================

setwd(dirname(rstudioapi::getSourceEditorContext()$path))
source("functions_trends.R")
library(exactextractr)
# ── Toggle this flag to also produce supplementary figures ───────────────────
RUN_SUPP <- FALSE   # set TRUE to generate all supplementary plots
# ─────────────────────────────────────────────────────────────────────────────


# =============================================================
# 0  PATHS & GLOBAL SETTINGS
# =============================================================

hydroDir  <- "D:/tilloal/Documents/LFRuns_utils/ChangingHydroExtremes/data"
# hydroDir2 <- "D:/tilloal/Documents/LFRuns_utils/data"
plotDir   <- "D:/tilloal/Documents/LFRuns_utils/TrendAnalysis/plots"

haz    <- "Flood"          # "Flood" or "Drought"
rl100  <- TRUE             # use 100-year RL matrices
mmx    <- "RL100"          # label used in filenames
it     <- 30               # iteration ID for plot filenames

# Shared aesthetics
tsize  <- 12
osize  <- 12
colorz <- c("Clim" = "dodgerblue4", "LUC" = "gold4",
             "Res"  = "firebrick4",  "WU"  = "olivedrab")
colorn <- c("WU"   = "limegreen",   "Res" = "tomato4",
             "LUC"  = "orange",      "Clim"= "royalblue")
clabels <- c("Climate", "Land use", "Reservoirs", "Water demand")

# Helper: shared map theme
map_theme <- function(tsize = 16, osize = 12, pos="right") {
  theme(
    axis.title  = element_text(size = tsize),
    axis.text   = element_text(size = osize),
    panel.background = element_rect(fill = "aliceblue", colour = "grey1"),
    panel.border     = element_rect(linetype = "solid", fill = NA, colour = "black"),
    legend.title     = element_text(size = osize, margin = margin(t=2,r=2,b=6,l=0)),
    legend.text      = element_text(size = osize),
    legend.spacing.x = unit(0.2, "cm"),
    legend.position  = pos,
    legend.box       = "horizontal",
    panel.grid.major = element_line(colour = "grey70"),
    panel.grid.minor = element_line(colour = "grey90"),
    legend.key       = element_rect(fill = "transparent", colour = "transparent"),
    legend.key.size  = unit(1, "cm")
  )
}

# Helper: shared boxplot/timeline theme
bxp_theme <- function() {
  theme(
    axis.title   = element_text(size = 18, face = "bold"),
    title        = element_text(size = 22, face = "bold"),
    axis.text    = element_text(size = 16),
    axis.text.x  = element_text(size = 16, face = "bold"),
    panel.background = element_rect(fill = "white", colour = "white"),
    panel.grid       = element_blank(),
    panel.border     = element_rect(linetype = "solid", fill = NA, colour = "black"),
    legend.title     = element_text(size = 20, face = "bold"),
    legend.text      = element_text(size = 16),
    axis.ticks.y     = element_blank(),
    panel.grid.major.y = element_line(color = "lightgray", linetype = "dashed"),
    panel.grid.minor.x = element_line(colour = "grey23",   linetype = "dashed"),
    legend.position  = "right",
    legend.key       = element_rect(fill = "transparent", colour = "transparent"),
    legend.key.size  = unit(0.8, "cm")
  )
}

# Helper: aggregate change by HydroRegion
aggregate_change <- function(Data, yrname, y_col = "Y2015") {
  
  pointagg <- aggregate(
    list(Rchange_rel = Data[[y_col]]),
    by  = list(HydroR = Data$HER),
    FUN = function(x) c(mean = mean(x, na.rm=TRUE), dev = sd(x, na.rm=TRUE),
                        len  = length(x),
                        med  = median(x, na.rm=TRUE),
                        q1   = quantile(x, 0.025, na.rm=TRUE),
                        q3   = quantile(x, 0.975, na.rm=TRUE))
  )
  point_out <- do.call(data.frame, pointagg)

  yrange <- match(yrname, colnames(Data))
  trendagg <- aggregate(
    list(Rchange = Data[, yrange]),
    by  = list(HydroR = Data$HER),
    FUN = function(x) c(mean = mean(x, na.rm=TRUE))
  )
  trend_out <- do.call(data.frame, trendagg)

  list(point = point_out, trend = trend_out)
}

# Helper: load one RL100 run
load_rl_run <- function(hydroDir, haz, namefile, rl100, Paramsfl_yr) {
  load(file = paste0(hydroDir, "/", haz, "/params.", namefile, ".Rdata"))
  Paramsfl <- Paramsfl[, -c(4:9)]

  if (rl100) {
    load(file = paste0(hydroDir, "/", haz, "/RL100x.", namefile, ".Rdata"))
    #RL100mat <- RLGPDfl
    # colnames(RL100mat)[1] <- "1951"
    # RL100mat <- cbind(RL100mat, Paramsfl$catch[which(Paramsfl$Year == "1955")])
    # colnames(RL100mat) <- paste0("Y", colnames(RL100mat))
    # colnames(RL100mat)[71] <- "unikout"
    # RLmat <- data.frame(RL100mat)
  } else {
    load(file = paste0(hydroDir, "/", haz, "/RL100.", namefile, ".Rdata"))
    RLmat <- RLGPDfl
  }
  load(file = paste0(hydroDir, "/", haz, "/peaks.", namefile, ".Rdata"))
  list(RLmat = RLmat, Params = data.table(Paramsfl), Peaks = Peaksave)
}

load_rl_only <- function(hydroDir, haz, namefile, rl100 = T) {
  if (rl100) {
    load(file = paste0(hydroDir, "/", haz, "/RL100x.", namefile, ".Rdata"))
    # colnames(RL100mat)[1] <- "1951"
    # RL100mat <- cbind(RL100mat, cats)
    # colnames(RL100mat) <- paste0("Y", colnames(RL100mat))
    # colnames(RL100mat)[71] <- "unikout"
    # RLmat <- data.frame(RL100mat)
    return(RLmat)
  } else {
    load(file = paste0(hydroDir, "/", haz, "/RL100.", namefile, ".Rdata"))
    return(RLGPDfl)
  }
}

load_params_peaks <- function(hydroDir, haz, namefile) {
  load(file = paste0(hydroDir, "/", haz, "/params.", namefile, ".Rdata"))
  Paramsfl <- Paramsfl[, -c(4:9,)]
  Params <- data.table(Paramsfl)
  load(file = paste0(hydroDir, "/", haz, "/peaks.", namefile, ".Rdata"))
  return(list(Params = Params, Peaks = Peaksave))
}



# =============================================================
# 1  SPATIAL DATA
# =============================================================

## 1.1  Outlets ------------------------------------------------
if (!exists("outf")) {
  rspace <- read.csv(paste0(hydroDir, "/subspace_efas.csv"))[, -1]
  outf   <- c()
  for (Nsq in 1:88) {
    print(Nsq)
    nrspace    <- rspace[Nsq, ]
    outletname <- "GeoData/efas_rnet_100km_01min"
    outhybas   <- outletopen(hydroDir, outletname, nrspace)
    Idstart    <- as.numeric(Nsq) * 10000
    Idstart2   <- as.numeric(Nsq) * 100000
    if (length(outhybas$outlets) > 0) {
      outhybas$outlets <- seq(Idstart + 1, Idstart + length(outhybas$outlets))
      outhybas$outl2   <- seq(Idstart2 + 1, Idstart2 + length(outhybas$outlets))
      outhybas$latlong <- paste(round(outhybas$Var1, 4), round(outhybas$Var2, 4), sep = " ")
      outf <- rbind(outf, outhybas)
    }
  }
}

## 1.2  HyBas07 ------------------------------------------------
Catchmentrivers7 <- read.csv(
  paste0(hydroDir, "/GeoData/HYBAS07/from_hybas_eu_onlyid.csv"),
  encoding = "UTF-8", header = TRUE, stringsAsFactors = FALSE
)
hybas07  <- read_sf(dsn = paste0(hydroDir, "/GeoData/HYBAS07/hybas_eu_lev07_v1c.shp"))
hybasf7  <- fortify(hybas07)
Catamere07 <- inner_join(hybasf7, Catchmentrivers7, by = "HYBAS_ID")
Catamere07$llcoord <- paste(round(Catamere07$POINT_X, 4),
                             round(Catamere07$POINT_Y, 4), sep = " ")
cst7 <- right_join(Catamere07, outf, by = c("llcoord" = "latlong"))
GNF  <- cst7
st_geometry(GNF) <- NULL
rm(Catamere07)

outhybas07 <- outletopen(hydroDir, "/GeoData/HYBAS07/outletsv8_hybas07_01min")
outhybas07$latlong <- paste(round(outhybas07$Var1, 4),
                             round(outhybas07$Var2, 4), sep = " ")
mhy <- match(outhybas07$latlong, outf$latlong)
outhybas07$outID <- outf$outl2[mhy]

## 1.3  Biogeographic regions ----------------------------------
biogeo    <- read_sf(dsn = paste0(hydroDir, "/GeoData/eea_3035_biogeo-regions_2016/BiogeoRegions2016_wag84.shp"))
biogeof   <- fortify(biogeo)
st_geometry(biogeof) <- NULL

biogeoregions <- raster(paste0(hydroDir, "/GeoData/eea_3035_biogeo-regions_2016/Biogeo_rasterized_wsg84.tif"))
Gbiogeoregions <- as.data.frame(biogeoregions, xy = TRUE)
biogeomatch <- inner_join(biogeof, Gbiogeoregions, by = c("PK_UID" = "Biogeo_rasterized_wsg84"))
biogeomatch$latlong <- paste(round(biogeomatch$x, 4), round(biogeomatch$y, 4), sep = " ")
biogeo_rivers <- right_join(biogeomatch, outf, by = "latlong")

## 1.4  HydroRegions -------------------------------------------
GridHR  <- raster(paste0(hydroDir, "/GeoData/HER/HydroRegions_raster_WGS84.tif"))
GHR     <- as.data.frame(GridHR, xy = TRUE)
GHR     <- GHR[!is.na(GHR[, 3]), ]
GHR$llcoord <- paste(round(GHR$x, 4), round(GHR$y, 4), sep = " ")
GHR_riv <- inner_join(GHR, outf, by = c("llcoord" = "latlong"))
GHshpp  <- read_sf(dsn = paste0(hydroDir, "/GeoData/HER/her_all_adjusted.shp"))
HydroRsf <- fortify(GHshpp)

## 1.5  Basemap & coords ---------------------------------------
outletname <- "/GeoData/efas_rnet_100km_01min"
outll      <- outletopen(hydroDir, outletname)
cord.dec   <- SpatialPoints(outll[, c(2, 3)], proj4string = CRS("+proj=longlat"))
cord.UTM   <- spTransform(cord.dec, CRS("+init=epsg:3035"))
nco        <- cord.UTM@coords
world   <- ne_countries(scale = "medium", returnclass = "sf")
Europe  <- world[world$continent == "Europe", ]
w2      <- st_transform(world, crs = 3035)
basemap <- w2

Impdates   <- seq(1950, 2020, by = 10)
valuenames <- paste0("Y", Impdates)
catmap     <- cst7
rm(cst7)

## 1.6  Reservoirs ---------------------------------------------
res2020 <- resOpen(hydroDir2, "/reservoirs/reservoirs_volumes_2020_Domain2.nc")
res2020$idla  <- 2970 - res2020$idla + 1
res2020$idlalo <- paste(res2020$idlo, res2020$idla, sep = " ")
res1951 <- resOpen(hydroDir2, "/reservoirs/reservoirs_volumes_1951.nc")

matres  <- na.omit(match(res1951$idlalo, res2020$idlalo))
res_old <- res2020[ matres, ];  res_old$status <- "old"
res_new <- res2020[-matres, ];  res_new$status <- "new"
res_comp <- left_join(res_old, res1951, by = "idlalo")
res_f    <- rbind(res_new, res_old)

pointout <- st_as_sf(res_new, coords = c("Var1", "Var2"), crs = 4326)
pointout <- st_transform(pointout, crs = 3035)

## 1.7  Upstream area ------------------------------------------
outf$idlalo <- paste(outf$idlo, outf$idla, sep = " ")
UpArea      <- UpAopen(hydroDir, "/GeoData/upArea_European_01min.nc", outf)

#range of flood separation
flood_sep<-function(x){
  sep=5+log(x/2.59)
  return(sep)
}

range=c(min(flood_sep(UpArea$upa)),max(flood_sep(UpArea$upa)))

max(UpArea$upa)
# =============================================================
# 2  LOAD MODEL RUNS
# =============================================================
haz="Drought"
# Filenames by hazard type
run_names <- if (haz == "Flood") {
  list(H    = "flood.year.Histo",
       SCF  = "Flood.year.socCF",
       RWCF = "flood.year.RWStat",
       WCF  = "flood.year.WStat")
} else {
  list(H    = "Drought.nonfrost.Histo",
       SCF  = "Drought.nonfrost.SocCF",
       RWCF = "Drought.nonfrost.RWStat",
       WCF  = "Drought.nonfrost.WStat")
}

#just load the SCF scenario for parameters

# runH    <- load_rl_run(hydroDir, haz, run_names$H,    rl100)
runSCF  <- load_rl_run(hydroDir, haz, namefile=run_names$SCF,  rl100=T)
RLGPDflSCF  <- as.data.table(runSCF$RLmat);  ParamsflSCF  <- runSCF$Params;  PeakSCF  <- runSCF$Peaks

#RLGPDflSCF=as.data.table(RLGPDflSCF)
# runRWCF <- load_rl_run(hydroDir, haz, run_names$RWCF, rl100)
# runWCF  <- load_rl_run(hydroDir, haz, run_names$WCF,  rl100)

cats=ParamsflSCF$catch[which(ParamsflSCF$Year == "1955")]
RLGPDflH <- load_rl_only(hydroDir, haz, run_names$H, rl100=T)
RLGPDflRWCF<- load_rl_only(hydroDir, haz, run_names$RWCF,rl100=T)
RLGPDflWCF <- load_rl_only(hydroDir, haz, run_names$WCF,rl100=T)

RLGPDflH=data.table(RLGPDflH)
RLGPDflRWCF=data.table(RLGPDflRWCF)
RLGPDflWCF=data.table(RLGPDflWCF)
gc()


# =============================================================
# 3  DATA CLEANING
# =============================================================

## 3.1  Fill missing last-year RL values -----------------------
fix_nan_rl <- function(RLcheck) {
  pbloc <- RLcheck$unikout[which(is.nan(RLcheck$Y2020))]
  if (length(pbloc) > 1) {
    for (fi in pbloc) {
      row_vec <- as.numeric(RLcheck[RLcheck$unikout == fi, ])
      last <- which(is.nan(row_vec))[1] - 1
      row_vec[which(is.nan(row_vec))] <- row_vec[last]
      
      # Assign back as a list
      RLcheck[RLcheck$unikout == fi, (names(RLcheck)) := as.list(row_vec)]
    }
  }
  RLcheck
}
RLGPDflSCF  <- fix_nan_rl(RLGPDflSCF)
RLGPDflRWCF <- fix_nan_rl(RLGPDflRWCF)
RLGPDflWCF  <- fix_nan_rl(RLGPDflWCF)
RLGPDflH    <- fix_nan_rl(RLGPDflH)

## 3.2  IRES status (drought only) -----------------------------
if (haz == "Drought") {
  load_ires <- function(hydroDir2, season, suffix) {
    load(file = paste0(hydroDir2, "/Drought/IRES.", season, ".", suffix, ".Rdata"))
    IRES_save$IRES[IRES_save$IRES == 2] <- 1
    IRES_save$IRES[IRES_save$IRES == 3] <- 2
    IRES_save
  }
  IRES_Histo <- load_ires(hydroDir, "nonfrost", "Histo")
  IRES_WCF   <- load_ires(hydroDir, "nonfrost", "WCF")
  IRES_RWCF  <- load_ires(hydroDir, "nonfrost", "RWCF")
  IRES_SocCF <- load_ires(hydroDir, "nonfrost", "SocCF")

  IRES_comb <- Reduce(function(a, b) full_join(a, b, by = "catlist"),
                      list(IRES_Histo, IRES_SocCF, IRES_WCF, IRES_RWCF))
  for (col in c("IRES.y", "IRES.x.x", "IRES.y.y"))
    IRES_comb[[col]][is.na(IRES_comb[[col]])] <- 0
  IRES_comb$gen_IR <- ceiling(rowSums(IRES_comb[, grep("IRES", names(IRES_comb))]) / 4)
  names(IRES_comb)[2:5] <- c("Histo", "SCF", "WCF", "RWCF")
  IR_locs <- which(IRES_comb$SCF == 1)
  IRpoints <- UpArea
} else {
  IRpoints <- UpArea
}

## 3.3  Drought RL sign correction (drought only) --------------
correct_drought_rl <- function(RL, IR_locs) {
  RL[IR_locs, 1:70] <- NA
  RL[, 1:70] <- -RL[, 1:70]
  yr_cols <- names(RL)[1:70]
  
  for (col in yr_cols) {
    RL[get(col) < 0  & !is.na(get(col)), (col) := 0 ]
    RL[is.infinite(get(col)),            (col) := NA]
  }
  RL
}
if (haz == "Drought") {
  RLGPDflH    <- correct_drought_rl(RLGPDflH,    IR_locs)
  RLGPDflWCF  <- correct_drought_rl(RLGPDflWCF,  IR_locs)
  RLGPDflSCF  <- correct_drought_rl(RLGPDflSCF,  IR_locs)
  RLGPDflRWCF <- correct_drought_rl(RLGPDflRWCF, IR_locs)
}

## 3.4  Remove pixels with unrealistic shape parameter ---------
shp_bnd <- if (haz == "Flood") c(1, -0.5) else c(0, -1.5)

get_shape <- function(Params) Params[, c(1, 2, 5)]

ShapeparSCF <- get_shape(ParamsflSCF)

Shapeparf <- ShapeparSCF[ShapeparSCF$Year == 2015, ]

out_of_bounds <- function(sp) sp > shp_bnd[1] | sp < shp_bnd[2]
rmpixs <- unique(c(
  which(out_of_bounds(Shapeparf$epsilonGPD))
))

rmp2 <- na.omit(unique(Shapeparf$catchment[rmpixs]))
length(rmp2)/length(Shapeparf$catchment)*100
Shapeparf <- Shapeparf[-match(rmp2, Shapeparf$catchment), ]
# Shapeparf$mean <- rowMeans(Shapeparf[, c("epsilonGPD", "V2", "V3", "V4")])
# Shapeparf$sd   <- apply(Shapeparf[, c("epsilonGPD","V2","V3","V4")], 1,
#                         function(x) sqrt(mean((x - mean(x))^2)))

## 3.5  Remove bad pixels from all RL matrices -----------------
if (length(rmpixs)>1){
for (obj in c("RLGPDflSCF", "RLGPDflWCF", "RLGPDflRWCF", "RLGPDflH")) {
  tmp <- get(obj)
  assign(obj, tmp[-match(rmp2, tmp$unikout), ])
}
}
#rm(ParamsflSCF); gc()


# =============================================================
# 4  CHANGE ATTRIBUTION
# =============================================================

cd <- as.numeric(which(colnames(RLGPDflSCF) == "Y1955"))

# Climtrend <- RLGPDflSCF[, 1:70] - RLGPDflSCF[, cd]          # Climate
# Soctrend  <- RLGPDflRWCF[, -71] - RLGPDflSCF[, -71]         # Land use
# Restrend  <- RLGPDflWCF[,  -71] - RLGPDflRWCF[, -71]        # Reservoirs
# Wutrend   <- RLGPDflH[,   -71]  - RLGPDflWCF[, -71]         # Water demand
# Totaltrend <- RLGPDflH[, 1:70]  - RLGPDflSCF[, cd]          # Total

Climtrend <- as.data.frame(RLGPDflSCF[, 1:70] - RLGPDflSCF[[cd]])
Soctrend  <- as.data.frame(RLGPDflRWCF[, -71] - RLGPDflSCF[, -71])
Restrend  <- as.data.frame(RLGPDflWCF[, -71] - RLGPDflRWCF[, -71])
Wutrend   <- as.data.frame(RLGPDflH[, -71] - RLGPDflWCF[, -71])
Totaltrend <- as.data.frame(RLGPDflH[, 1:70] - RLGPDflSCF[[cd]])

UpAvec  <- UpArea[, c(12, 3)]
yrname  <- colnames(Climtrend)
unikout <- RLGPDflSCF[, 71]
Regio   <- HydroRsf

# Reference RL (symmetric average between SCF-1955 and Historical)
DataI        <- data.frame(RLGPDflH)
c1           <- as.numeric(which(colnames(DataI) == "Y1955"))
irange       <- match(yrname, colnames(DataI))
DataI$Init   <- RLGPDflSCF[[c1]]
DataI[, irange] <- (DataI$Init + DataI[, irange]) / 2

# Biogeo – HyBas matching for regional aggregation
bio_hybas   <- inner_join(biogeo_rivers, GNF, by = "outl2")
HRM         <- na.omit(unique(match(bio_hybas$HYBAS_ID, hybasf7$HYBAS_ID)))
hybasf7_dom <- hybasf7[HRM, ]
bhp_m       <- na.omit(unique(match(hybasf7_dom$HYBAS_ID, bio_hybas$HYBAS_ID)))
hybasf7_dom$biogeoreg <- bio_hybas$code[bhp_m]

herziz      <- na.omit(match(GHR_riv$HydroRegions_raster_WGS84, GHshpp$Id))
GHR_riv$HER <- GHshpp$CODEB[herziz]
rmpixels    <- rmpixs


# =============================================================
# 5  DRIVER DATA PROCESSING
# =============================================================

## 5.1  Land use -----------------------------------------------
DataL <- ComputeChange(Soctrend, unikout, DataI, outhybas07,
                       rmpixels, UpAvec, GHR_riv, HydroRsf, yrname, "socio", eps=0.1)
yrange <- match(yrname, colnames(DataL))

res_L <- aggregate_change(DataL, yrname); pointSoc <- res_L$point; trendSoc <- res_L$trend

## 5.2  Water demand -------------------------------------------
DataW <- ComputeChange(Wutrend, unikout, DataI, outhybas07,
                       rmpixels, UpAvec, GHR_riv, HydroRsf, yrname, "socio", eps=0.1)

res_W <- aggregate_change(DataW, yrname); pointWu <- res_W$point; trendWu <- res_W$trend

## 5.3  Total change -------------------------------------------
DataT <- ComputeChange(Totaltrend, unikout, DataI, outhybas07,
                       rmpixels, UpAvec, GHR_riv, HydroRsf, yrname, "total", eps=0.1)

res_T <- aggregate_change(DataT, yrname); pointTot <- res_T$point; trendTot <- res_T$trend

## 5.4  Climate ------------------------------------------------
DataC <- ComputeChange(Climtrend, unikout, DataI, outhybas07,
                       rmpixels, UpAvec, GHR_riv, HydroRsf, yrname, "clim", eps=0.1)
res_C <- aggregate_change(DataC, yrname); pointClim <- res_C$point; trendClim <- res_C$trend

## 5.5  Reservoirs ---------------------------------------------
DataR <- ComputeChange(Restrend, unikout, DataI, outhybas07,
                       rmpixels, UpAvec, GHR_riv, HydroRsf, yrname, "socio", eps=0.1)


res_R <- aggregate_change(DataR, yrname); pointRes <- res_R$point; trendRes <- res_R$trend

## 5.6  Trend data for biogeographic summaries -----------------
TdataClim <- processTrendData(trendClim, DataC, id_var = "HydroR")
TdataLuse <- processTrendData(trendSoc,  DataL, id_var = "HydroR")
TdataRes  <- processTrendData(trendRes,  DataR, id_var = "HydroR")
TdataWuse <- processTrendData(trendWu,   DataW, id_var = "HydroR")

## 5.7  Driver dominance at pixel & HER level ------------------
# HER level
PointAI <- data.frame(
  clim = abs(pointClim$Rchange_rel.mean), lu  = abs(pointSoc$Rchange_rel.mean),
  res  = abs(pointRes$Rchange_rel.mean),  wu  = abs(pointWu$Rchange_rel.mean)
)
max_col_HER <- as.numeric(apply(PointAI, 1, which.max))

pointap <- full_join(GHshpp, pointClim, by = c("CODEB" = "HydroR"))
st_geometry(pointap) <- NULL
cmat    <- match(Regio$Id, pointClim$HydroR)
pointap <- pointap[!is.na(cmat), ]; cmat <- cmat[!is.na(cmat)]
pointap$maxcol <- max_col_HER[cmat]
pointap <- pointap[!is.na(pointap$maxcol), ]
pointplot <- left_join(HydroRsf, pointap, by = c("CODEB" = "Id"))
pointplot <- pointplot[!is.na(pointplot$maxcol), ]
nutplot   <- st_transform(pointplot, crs = 3035)

# Pixel level
idc   <- match(c("Y2020","upa"), colnames(DataC))
idl   <- c(3, idc)
DataAI <- data.frame(
  clim = abs(DataC[, idl]$Y2020), lu   = abs(DataL[, idl]$Y2020),
  resw = abs(DataR[, idl]$Y2020), wu   = abs(DataW[, idl]$Y2020)
)
max_col_pix <- as.numeric(apply(DataAI, 1, which.max))

datap <- DataC; datap$maxcol <- max_col_pix
datap <- datap[!is.na(datap$Var1), ]
points_pix <- st_as_sf(datap, coords = c("Var1","Var2"), crs = 4326)
points_pix <- st_transform(points_pix, crs = 3035)
points_pix <- points_pix[!is.na(points_pix$maxcol), ]

# Summary tables for trend loop
pointsAD_list <- list(
  Dchangelist = list(),   # filled in section 6
  paggC = NULL, paggL = NULL, paggR = NULL, paggW = NULL
)


# =============================================================
# 6  MAIN FIGURES
# =============================================================

## ── Figure 1 ─────────────────────────────────────────────────
## Boxplot in time – pan-European driver contributions
trtF_EU <- rbind(TdataClim$tGlobal, TdataLuse$tGlobal,
                 TdataRes$tGlobal,   TdataWuse$tGlobal)
trtF_EU$year   <- rep(seq(1950, 2010, 10), 4)
trtF_EU$driver <- rep(c("Clim","LUC","Res","WU"), each = 7)
names(trtF_EU)[c(2,4,5,6,7,8)] <- c("changeC","med","cq1","cq2","w1","w2")

xlabs <- seq(1950, 2010, 10)
br_bxp <- c(seq(-100,-10,10), seq(-5,5,5), seq(10,100,10))

fig1 <- ggplot(trtF_EU) +
  geom_linerange(aes(x=year, ymin=w1, ymax=w2, color=factor(driver), group=factor(driver)),
                 position = position_dodge2(width = 9), lwd = 1, alpha = 0.8) +
  scale_color_manual(values = colorn, name = "Drivers", labels = clabels) +
  new_scale_color() +
  geom_rect(aes(xmin=year-4.5, xmax=year+4.5, ymin=cq1, ymax=cq2,
                fill=factor(driver), group=factor(driver)),
            alpha = 0.5, position = position_dodge(width = 9)) +
  geom_point(aes(x=year, y=changeC, color=factor(driver), group=factor(driver)),
             position = position_dodge(width = 9), size = 3) +
  geom_rect(aes(xmin=year-4.5, xmax=year+4.5, ymin=med-1e-1, ymax=med+1e-1,
                fill=factor(driver), group=factor(driver)),
            alpha = 1, position = position_dodge(width = 9)) +
  scale_y_continuous(name = "Mean change (% 10yRL)", breaks = br_bxp,
                     trans = scales::modulus_trans(.6)) +
  scale_x_continuous(breaks = xlabs, labels = xlabs, name = "Decades",
                     minor_breaks = seq(1955,2005,10), expand = c(.01,.01)) +
  scale_fill_manual(values = colorn, name = "Drivers", labels = clabels) +
  scale_color_manual(values = colorz, name = "Drivers", labels = clabels) +
  guides(color = guide_legend(override.aes = list(color = colorz))) +
  bxp_theme() +
  ggtitle("Europe")

fig1
# ggsave(paste0(plotDir, "/Fig1_bxp_EU_", haz, "_FINAL.jpg"),
#        fig1, width=30, height=20, units="cm", dpi=1000)
# 
# write.csv(trtF_EU,file=paste0(plotDir,"/",haz,"_agchanges.csv"))

## ── Figures 2 & 3 ────────────────────────────────────────────
## Spatial maps: change per driver (loop)
yrlist     <- 1951:2020
driverlist <- c("climate","all")
Dchangelist <- list()

colorz_drv <- c("4"="limegreen","3"="tomato4","2"="orange","1"="royalblue")
lab_drv    <- c("Climate","Land use","Reservoirs","Water demand")

for (driver in driverlist) {
  trendPlot <- switch(driver,
    climate    = trendClim,  all = trendTot)
  datap_d   <- switch(driver,
    climate    = DataC,  all = DataT)
  pointagg_d <- switch(driver,
    climate = pointClim, all = pointTot)

  Pplot <- calculatePoints(trendPlot, yrlist, pointagg_d, Regio, GHshpp, datap_d)
  Dchangelist <- c(Dchangelist, list(Pplot))

  colNA   <- "transparent"
  points  <- Pplot$points
  pag     <- Pplot$PagD
  pointsInside <- Pplot$psp
  length(unique(pointsInside$CODEB[which(pointsInside$change>0)]))

  if (haz == "Flood") {
    br    <- c(-50,-20,-10,-5,0,5,10,20,50); limi <- c(-50,50)
    palet  <- hcl.colors(11, "RdYlBu", rev=FALSE)
    paletf <- hcl.colors(11, "RdBu",   rev=FALSE)
    tsize_m <- 16; osize_m <- 12
    legend2 <- "Change (%)    "

    if (driver %in% c("climate","all")) {
      titleX <- paste0("Change in 10-year ",haz," attributed\nto ",driver,
                       " changes (% of 100y flood) — 1955–2015")
      fig_map <- ggplot(basemap) +
        geom_sf(fill="white", color="darkgrey", size=0.5) +
        geom_sf(data=pag, aes(fill=Rchange_rel.mean, geometry=geometry),
                alpha=0.2, color="transparent") +
        geom_sf(data=points, aes(col=Y2015, geometry=geometry, size=upa),
                alpha=.9, stroke=0, shape=15) +
        scale_size(range=c(0.08,0.4), trans="sqrt", guide="none") +
        scale_fill_gradientn(colors=palet, breaks=br, limits=limi,
                             trans=scales::modulus_trans(.3),
                             oob=scales::squish, na.value=colNA, name=legend2) +
        guides(fill="none") + new_scale_fill() +
        geom_sf(data=pointsInside,aes(geometry=geometry, fill=sign),
                alpha=.6,size=.5,stroke=0,shape=21,color="black")+
        scale_fill_manual(values=c("tomato4","steelblue4"), name="Significant trends") +
        coord_sf(xlim=c(min(nco[,1]),max(nco[,1])), ylim=c(min(nco[,2]),max(nco[,2]))) +
        scale_color_gradientn(colors=palet, breaks=br, limits=limi,
                              trans=scales::modulus_trans(.3),
                              oob=scales::squish, na.value=colNA, name=legend2) +
        labs(x="Longitude", y="Latitude") +
        guides(colour=guide_colourbar(barwidth=22, barheight=1),
               fill=guide_legend(override.aes=list(size=10))) +
        map_theme(tsize_m, osize_m,pos="bottom") +
       ggtitle(titleX)
    }
    # } else {
    #   titleX <- paste0("Change in 10-year ",haz," attributed\nto ",driver,
    #                    " changes (% of 100y flood) — 1955–2015")
    #   fig_map <- ggplot(basemap) +
    #     geom_sf(fill="white", color="darkgrey", size=0.5) +
    #     geom_sf(data=pag, aes(fill=Rchange_rel.mean, geometry=geometry),
    #             alpha=0.2, color="transparent") +
    #     geom_sf(data=points, aes(col=Y2015, geometry=geometry, size=upa),
    #             alpha=.9, stroke=0, shape=15) +
    #     scale_size(range=c(0.1,0.5), trans="sqrt", guide="none") +
    #     scale_fill_gradientn(colors=paletf, breaks=br, limits=limi,
    #                          trans=scales::modulus_trans(.3), labels=br,
    #                          oob=scales::squish, na.value=colNA, name="Change (%)") +
    #     coord_sf(xlim=c(min(nco[,1]),max(nco[,1])), ylim=c(min(nco[,2]),max(nco[,2]))) +
    #     scale_color_gradientn(colors=paletf, breaks=br, limits=limi, labels=br,
    #                           trans=scales::modulus_trans(.3),
    #                           oob=scales::squish, na.value=colNA, name="Change (%)") +
    #     labs(x="Longitude", y="Latitude") +
    #     guides(colour=guide_colourbar(barwidth=1.5, barheight=14),
    #            fill=guide_colourbar(barwidth=1.5, barheight=14)) +
    #     map_theme(tsize_m, osize_m) + ggtitle(titleX)
    # }

  } else {
    # Drought palette scheme
    br <- c(-50,-20,-10,-5,0,5,10,20,50); limi <- c(-50,50)
    palet  <- hcl.colors(11, "RdYlBu", rev=FALSE)
    tsize_m <- 16; osize_m <- 12
    legend2 <- "Change (%)    "

    titleX <- paste0("Change in 10-year ",haz," attributed\nto ",driver,
                     " changes (% of 10y drought) — 1955–2015")
    fig_map <- ggplot(basemap) +
      geom_sf(fill="white", color="darkgrey", size=0.5) +
      geom_sf(data=pag, aes(fill=Rchange_rel.mean, geometry=geometry),
              alpha=0.2, color="transparent") +
      geom_sf(data=points, aes(col=Y2015, geometry=geometry, size=upa),
              alpha=.9, stroke=0, shape=15) +
      scale_size(range=c(0.08,0.4), trans="sqrt", guide="none") +
      scale_fill_gradientn(colors=palet, breaks=br, limits=limi,
                           trans=scales::modulus_trans(.3),
                           oob=scales::squish, na.value=colNA, name=legend2) +
      guides(fill="none") + new_scale_fill() +
      geom_sf(data=pointsInside,aes(geometry=geometry, fill=sign),
              alpha=.6,size=.5,stroke=0,shape=21,color="black")+
      scale_fill_manual(values=c("tomato4","steelblue4"), name="Significant trends") +
      coord_sf(xlim=c(min(nco[,1]),max(nco[,1])), ylim=c(min(nco[,2]),max(nco[,2]))) +
      scale_color_gradientn(colors=palet, breaks=br, limits=limi,
                            trans=scales::modulus_trans(.3),
                            oob=scales::squish, na.value=colNA, name=legend2) +
      labs(x="Longitude", y="Latitude") +
      guides(colour=guide_colourbar(barwidth=22, barheight=1),
             fill=guide_legend(override.aes=list(size=10))) +
      map_theme(tsize_m, osize_m,pos="bottom") + ggtitle(titleX)
  }

  fname <- if (haz=="Flood" && driver %in% c("climate","all")) "Fig2" else "Fig3"
  ggsave(paste0(plotDir, "/", fname, "_mapF_", driver, "_", haz, "_FINAL.jpg"),
         fig_map, width=23, height=20, units="cm", dpi=1000)
}



## ============================================================
##  Figures 2 & 3 – Driver-specific spatial plots
##  Replaces the generic driver loop in 03_CHEX_Plot_univariateF_rev.R
## ============================================================
##
##  Assumed pre-loaded objects (from the main script):
##    basemap, nco, haz, it, plotDir
##    trendClim/Soc/Res/Wu/Tot, DataC/L/R/W/T
##    pointClim/Soc/Res/Wu/Tot, Regio, GHshpp
##    Dchangelist (filled below), pointsAD
##    map_theme()  (helper defined in main script)
##  Paths:
##    workDir  <- "D:/tilloal/Documents/06_Floodrivers"
##    hydroDir <- "D:/tilloal/Documents/LFRuns_utils/data"
## ============================================================


# =============================================================
# 0  SHARED AESTHETICS
# =============================================================

colNA    <- "transparent"
yrlist   <- 1951:2020

# Colour scale (same for climate & total – kept identical to original)
palet_clim  <- hcl.colors(11, "RdYlBu", rev = FALSE)

if (haz == "Flood") {
  br_map <- c(-50, -20, -10, -5, 0, 5, 10, 20, 50)
  limi   <- c(-50, 50)
  palet  <- palet_clim
  pct_label <- "% of 10-y flood"
} else {
  br_map <- c(-50, -20, -10, -5, 0, 5, 10, 20, 50)
  limi   <- c(-50, 50)
  palet <- palet_clim
  pct_label <- "% of 10-y drought"
}

tsize_m <- 16
osize_m <- 12

# =============================================================
# 2  RESERVOIRS  –  newly built dams overlay
# =============================================================

## 2.1  Load data ----------------------------------------------

# New reservoirs (post-1951, not present in 1951 file)
res2020 <- resOpen(hydroDir, "/reservoirs/reservoirs_volumes_2020_Domain2.nc")
res2020$idla  <- 2970 - res2020$idla + 1
res2020$idlalo <- paste(res2020$idlo, res2020$idla, sep = " ")
res1951 <- resOpen(hydroDir, "/reservoirs/reservoirs_volumes_1951.nc")

matres  <- na.omit(match(res1951$idlalo, res2020$idlalo))
res_new <- res2020[-matres, ]   # dams built AFTER 1951

# Convert to sf, project to EPSG:3035
res_new_sf <- st_as_sf(res_new, coords = c("Var1", "Var2"), crs = 4326)
res_new_sf <- st_transform(res_new_sf, crs = 3035)

res_ratio=outletopen(hydroDir,"/reservoirs/res_ratio_diff_2020-1951")
#match res_ratio with reservoirs
res_ratio$latlong=paste(round(res_ratio$Var1,4), round(res_ratio$Var2,4), sep = " ")
res_new$latlong=paste(round(res_new$Var1,4), round(res_new$Var2,4), sep = " ")
res_ratioX=res_ratio[which(res_ratio$outlets>5),]


direct_match <- match(res_new$latlong, res_ratio$latlong)
res_matched   <- res_new[ !is.na(direct_match), ]
res_matched$res_ratio_change <- res_ratio$outlets[direct_match[!is.na(direct_match)]]

# Reservoirs that need downstream snapping
res_unmatched <- res_new[is.na(direct_match), ]

# ── 2. Build sf objects in EPSG:3035 ────────────────────────────────────────
res_unmatched_sf <- st_transform(
  st_as_sf(res_unmatched, coords = c("Var1", "Var2"), crs = 4326),
  crs = 3035
)
res_ratioX_sf <- st_transform(
  st_as_sf(res_ratio, coords = c("Var1", "Var2"), crs = 4326),
  crs = 3035
)

# ── 3. Snap each unmatched reservoir to its nearest downstream pixel ─────────
#
#  Strategy: among all res_ratioX pixels, the "nearest downstream" one is
#  simply the spatially closest — because at 1-min resolution, moving
#  downstream always means moving to an adjacent cell, so the nearest large
#  pixel in res_ratioX is overwhelmingly likely to be downstream.
#  If you have a flow-accumulation raster you can add a filter (see note below).

nearest_idx <- st_nearest_feature(res_unmatched_sf, res_ratioX_sf)

res_unmatched$res_ratio_change <- res_ratio$outlets[nearest_idx]
res_unmatched$matched_latlong  <- res_ratio$latlong[nearest_idx]   # for QC

# Distance to matched pixel — flag anything suspiciously far (> ~50 km)
snap_dist <- as.numeric(st_distance(
  res_unmatched_sf,
  res_ratioX_sf[nearest_idx, ],
  by_element = TRUE
))
res_unmatched$snap_dist_km <- snap_dist / 1000
res_unmatched <- res_unmatched[res_unmatched$snap_dist_km < 50, ]   # drop outliers

# ── 4. Combine and flag large ratio increases ────────────────────────────────
res_ratioR <- bind_rows(
  res_matched,
  res_unmatched
)

# Threshold: reservoirs where the ratio change is in the top quartile
# (adjust the percentile or use a fixed threshold like > 1 depending on your units)
thresh <- quantile(res_ratioR$res_ratio_change, 1, na.rm = TRUE)
res_ratioR$large_increase <- res_ratioR$res_ratio_change > 3
res_ratioRX=res_ratioR[which(res_ratioR$large_increase==T),]
# sf version for plotting
res_ratioR_sf <- st_transform(
  st_as_sf(res_ratioRX, coords = c("Var1", "Var2"), crs = 4326),
  crs = 3035
)

## 2.2  Background: reservoir RL change map --------------------

datap_res        <- DataR
Pplot_res        <- calculatePoints(trendRes, yrlist, pointRes, Regio, GHshpp, datap_res)
Dchangelist      <- c(Dchangelist, list(Pplot_res))

pts_res  <- Pplot_res$points
pag_res  <- Pplot_res$PagD
ptsIn_res <- Pplot_res$psp

titleX_res <- paste0("Change in 10-year ", haz, " attributed to reservoir changes\n(",
                     pct_label, ") — 1955–2015  |  diamonds = new reservoirs (post-1951)")

## 2.3  Plot ---------------------------------------------------
fig_res <- ggplot(basemap) +
  geom_sf(fill = "white", color = "darkgrey", size = 0.5) +
  # HER-level background shading
  geom_sf(data = pag_res, aes(fill = Rchange_rel.mean, geometry = geometry),
          alpha = 0.6, color = "transparent") +
  # River pixel colours
  geom_sf(data = pts_res, aes(col = Y2015, geometry = geometry, size = upa),
          alpha = .9, stroke = 0, shape = 15) +
  scale_size(range = c(0.08, 0.4), trans = "sqrt", guide = "none") +
  scale_fill_gradientn(colors = palet, breaks = br_map, limits = limi,
                       trans = scales::modulus_trans(.3),
                       oob = scales::squish, na.value = colNA, name = "Change (%)") +
  scale_color_gradientn(colors = palet, breaks = br_map, limits = limi,
                        trans = scales::modulus_trans(.3),
                        oob = scales::squish, na.value = colNA, name = "Change (%)") +
  guides(fill = "none") +
  new_scale_fill() +
  #new_scale_size() +
  new_scale("size")+
  # New reservoirs: diamonds sized by volume
  geom_sf(data = res_new_sf, aes( size=res,geometry = geometry),
          color = "black", fill = "orange", shape = 23, alpha = 0.55, stroke = 0.4) +
  scale_size(range = c(0.1, 1), trans = "sqrt",
             name = expression(paste("Reservoir volume (m"^3, ")")),
             breaks = c(1e5, 1e6, 1e7, 1e8, 1e9),
             labels = c("100 k", "1 M", "10 M", "100 M", "1 B"),
             guide = guide_legend(direction = "horizontal",
                                  title.position = "top", 
                                  label.position = "top",
                                  nrow = 1)) +
  coord_sf(xlim = c(min(nco[, 1]), max(nco[, 1])),
           ylim = c(min(nco[, 2]), max(nco[, 2]))) +

  labs(x = "Longitude", y = "Latitude") +
  guides(colour = guide_colourbar(barwidth = 1.5, barheight = 14)) +
  map_theme(tsize_m, osize_m) +
  ggtitle(titleX_res)

ggsave(paste0(plotDir, "/Fig3_mapF_reservoirs_", haz, "_FINAL.jpg"),
       fig_res, width = 23, height = 20, units = "cm", dpi = 800)


# =============================================================
# 3  LAND USE  –  sealed-area bubbles + forest-gain contour
# =============================================================

## 3.1  Load land use rasters ----------------------------------
workDir <- "D:/tilloal/Documents/06_Floodrivers"

rast_forest <- raster(paste0(workDir, "/landuse/fracforest_ch20201951.tif"))
rast_sealed <- raster(paste0(workDir, "/landuse/fracsealed_ch20201951.tif"))
rast_water <- raster(paste0(workDir, "/landuse/fracwater_ch20201951.tif"))

# water_df= as.data.frame(rast_water, xy = TRUE, na.rm = TRUE)
# df_lakes_new <- water_df %>% 
#   filter(fracwater_ch20201951 >=1)
# 
# df_lakes_new_laea <- df_lakes_new %>%
#   st_as_sf(coords = c("x", "y"), crs = 4326) %>%
#   st_transform(3035) %>%
#   mutate(x = st_coordinates(.)[,1],
#          y = st_coordinates(.)[,2]) %>%
#   st_drop_geometry()

water_pts <- as.data.frame(rast_water, xy = TRUE, na.rm = TRUE) %>%
  filter(fracwater_ch20201951 >= 1) %>%
  st_as_sf(coords = c("x", "y"), crs = 4326) %>%
  st_transform(3035)   # now water_pts is an sf object with geometry

## 3.2  Aggregate to HER level ---------------------------------
# Project HER shapefile to raster CRS for extraction
HydroRsf_wgs <- st_transform(HydroRsf, crs = st_crs(rast_sealed))

# Extract mean sealed and forest change per HER polygon
HydroRsf_wgs$sealed_change <- exact_extract(rast_sealed, HydroRsf_wgs, "mean")*100
HydroRsf_wgs$forest_change <- exact_extract(rast_forest, HydroRsf_wgs, "mean")*100

# Compute centroids in WGS84, then project to 3035 for plotting
her_centroids_wgs  <- st_point_on_surface(HydroRsf_wgs)
her_centroids_3035 <- st_transform(her_centroids_wgs, crs = 3035)

her_centroid_high <- her_centroids_3035[!is.na(her_centroids_3035$sealed_change) &
                                          her_centroids_3035$sealed_change > 1, ]

# HER polygons in 3035 for the forest-gain contour
HydroRsf_3035 <- st_transform(HydroRsf_wgs, crs = 3035)

# Subset HERs where mean forest fraction increased > 0.2
her_forest_high <- HydroRsf_3035[!is.na(HydroRsf_3035$forest_change) &
                                   HydroRsf_3035$forest_change > 10, ]

res_new_sf_large <- res_new_sf[which(res_new_sf$res>1000000000),]
## 3.3  Background: land use RL change map ---------------------
Pplot_lu   <- calculatePoints(trendSoc, yrlist, pointSoc, Regio, GHshpp, DataL)
Dchangelist <- c(Dchangelist, list(Pplot_lu))

pts_lu    <- Pplot_lu$points
pag_lu    <- Pplot_lu$PagD

titleX_lu <- paste0("Change in 10-year ", haz,
                    " attributed to land use changes\n(",
                    pct_label,
                    ") — 1955–2015\n",
                    "Circles = Sealed area increas > 1%  |  ",
                    "Dark green outline = forest increase > 10%")

## 3.4  Sealed-area bubble colour scale  -----------------------
# Centre colour scale around 0 (increase = red, decrease = blue)
sealed_lim   <- max(abs(her_centroids_3035$sealed_change), na.rm = TRUE)
sealed_lim   <- ceiling(sealed_lim * 10) / 10   # round up to 1 dp


## 3.5  Plot ---------------------------------------------------
fig_lu <- ggplot(basemap) +
  geom_sf(fill = "white", color = "darkgrey", size = 0.5) +
  # HER background shading (RL change)
  geom_sf(data = pag_lu, aes(fill = Rchange_rel.mean, geometry = geometry),
          alpha = 0.6, color = "transparent") +
  # River pixel colours
  geom_sf(data = pts_lu, aes(col = Y2015, geometry = geometry, size = upa),
          alpha = .9, stroke = 0, shape = 15) +
  scale_size(range = c(0.08, 0.4), trans = "sqrt", guide = "none") +
  scale_fill_gradientn(colors = palet_clim, breaks = br_map, limits = limi,
                       trans = scales::modulus_trans(.3),
                       oob = scales::squish, na.value = colNA, name = "RL change (%)") +
  scale_color_gradientn(colors = palet_clim, breaks = br_map, limits = limi,
                        trans = scales::modulus_trans(.3),
                        oob = scales::squish, na.value = colNA, name = "Change (%)") +
  guides(fill = "none") +
  new_scale_fill() +
  new_scale_color() +
  new_scale("size") +
# 
#   geom_sf(data = res_ratioR_sf, aes(geometry = geometry),
#           color = "black", fill = "purple",size=0.5, shape = 23, alpha = 0.55, stroke = 0.1) +
#   
  #Forest-gain contour (dark green outline on HERs with high forest increase)
  geom_sf(data = her_forest_high, fill = "transparent",
          color = "darkgreen", linewidth = .2, linetype = "solid") +
  # Sealed-area bubbles at HER centroids
  geom_sf(data = her_centroid_high,
          aes(size = sealed_change, geometry = geometry),
          shape = 21, color = "grey20", fill="grey20" ,stroke = 0.1, alpha = 0.6) +
  # scale_fill_gradient(low = "grey70", high = "grey30",
  #                     
  #                     name = "Mean sealed\narea change\n(%)",
  #                      guide = guide_colorbar(
  #                        direction = "horizontal", 
  #                        title.position = "top", 
  #                        barwidth = 10,         # You can control the width of the bar here
  #                        order = 4
  #                      ))+
  # scale_fill_gradient2(low  = "royalblue3", mid = "white", high = "firebrick3",
  #                      midpoint = 0,
  #                      limits = c(-sealed_lim, sealed_lim),
  #                      oob = scales::squish,
  #                      name = "Mean sealed\narea change\n(fraction)") +
  scale_size_continuous(range = c(.6, 4),
                        name  = "sealed change (%)",
                        breaks = c(1, 2, 3,4,5),
                        labels = c("1","2", "3", "4","5"),
                        guide = guide_legend(direction = "horizontal",
                                             title.position = "top", 
                                             label.position = "top",
                                             nrow = 1)) +
  coord_sf(xlim = c(min(nco[, 1]), max(nco[, 1])),
           ylim = c(min(nco[, 2]), max(nco[, 2]))) +
  labs(x = "Longitude", y = "Latitude") +
  guides(fill = guide_colourbar(barwidth = 1.5, barheight = 12)) +
  map_theme(tsize_m, osize_m) +
  theme(legend.box = "vertical")+ # Stacks the multiple legends)
  ggtitle(titleX_lu)

ggsave(paste0(plotDir, "/Fig3_mapF_landuse_", haz, "_FINAL.jpg"),
       fig_lu, width = 23, height = 20, units = "cm", dpi = 1000)


# =============================================================
# 4  WATER DEMAND  –  demand-change bubbles at NUTS3 level
# =============================================================

## 4.1  Load water demand rasters & NUTS3 ----------------------
nuts3_shp     <- read_sf(paste0(hydroDir, "/GeoData/NUTS3/NUTS3_modified.shp"))
rast_wd_total <- raster(paste0(hydroDir, "/wateruse/all_ysum_ch20201951.tif"))
rast_wd_2020  <- raster(paste0(hydroDir, "/wateruse/all_demands_2020.tif"))
rast_wd_1951  <- raster(paste0(hydroDir, "/wateruse/all_demands_1951.tif"))

# Project NUTS3 to raster CRS for extraction, then reproject to 3035 for plotting
Wd_wgs <- st_transform(nuts3_shp, crs = st_crs(rast_wd_total))

## 4.2  Aggregate water demand to NUTS3 ------------------------
Wd_wgs$wd_change    <- exact_extract(rast_wd_total, Wd_wgs, "sum")
Wd_wgs$wd_2020      <- exact_extract(rast_wd_2020,  Wd_wgs, "sum")
Wd_wgs$wd_1951      <- exact_extract(rast_wd_1951,  Wd_wgs, "sum")

Wd_wgs$wd_pct_change=Wd_wgs$wd_change/Wd_wgs$wd_1951*100
# Remove regions with negligible demand change
length(which(abs(Wd_wgs$wd_pct_change) > 50))/length(Wd_wgs$wd_change)
Wd_wgs <- Wd_wgs[ abs(Wd_wgs$wd_pct_change) > 50, ]
Wd_wgs$sign="1"
Wd_wgs$sign[which(Wd_wgs$wd_change<0)]="-1"
nuts3_3035 <- st_transform(Wd_wgs, crs = 3035)
# NUTS3 centroids in 3035
nuts3_centroids_3035 <- st_transform(st_point_on_surface(Wd_wgs), crs = 3035)

## 4.3  Background: water demand RL change map -----------------
Pplot_wu    <- calculatePoints(trendWu, yrlist, pointWu, Regio, GHshpp, DataW)
Dchangelist <- c(Dchangelist, list(Pplot_wu))

pts_wu  <- Pplot_wu$points
pag_wu  <- Pplot_wu$PagD

titleX_wu <- paste0("Change in 10-year ", haz,
                    " attributed to water demand changes\n(",
                    pct_label,
                    ") — 1955–2015\n",
                    "Circles = NUTS3 total water demand change (sum, m³/yr) | ",
                    "colour = sign of change")

## 4.4  Demand-change colour scale  ----------------------------
wd_lim   <- quantile(abs(nuts3_centroids_3035$wd_change), 0.95, na.rm = TRUE)
wd_lim   <- max(wd_lim, 1e4)   # floor to avoid degenerate scale

## 4.5  Plot ---------------------------------------------------
fig_wu <- ggplot(basemap) +
  geom_sf(fill = "white", color = "darkgrey", size = 0.5) +
  # HER background shading (RL change)
  geom_sf(data = pag_wu, aes(fill = Rchange_rel.mean, geometry = geometry),
          alpha = 0.6, color = "transparent") +
  # River pixel colours
  geom_sf(data = pts_wu, aes(col = Y2015, geometry = geometry, size = upa),
          alpha = .9, stroke = 0, shape = 15) +
  scale_size(range = c(0.08, 0.4), trans = "sqrt", guide = "none") +
  scale_fill_gradientn(colors = palet_clim, breaks = br_map, limits = limi,
                       trans = scales::modulus_trans(.3),
                       oob = scales::squish, na.value = colNA, name = "RL change (%)") +
  scale_color_gradientn(colors = palet_clim, breaks = br_map, limits = limi,
                        trans = scales::modulus_trans(.3),
                        oob = scales::squish, na.value = colNA, name = "Change (%)") +
  guides(fill = "none") +
  new_scale_fill() +
  new_scale_color() +
  new_scale("size")+
  # Water demand bubbles at NUTS3 centroids
  # Size = absolute change (so both increases & decreases show as large circles)
  # Fill = signed change (blue = decrease, red = increase)
  # Option A: Bicolour without mid‑range white
  geom_sf(data = nuts3_3035, aes(color=sign), fill = NA, linewidth = .1,
          linetype = "solid",alpha=0.2) +
  # geom_sf(data = nuts3_centroids_3035,
  #         aes(size = abs(wd_change), color = sign, geometry = geometry),
  #         shape = 21, fill = NA, stroke = 0.5, alpha = 0.5) +
  
  scale_size_continuous(range = c(.3, 2),
                        trans  = "sqrt",
                        name   = "|Demand change|\n(m³/yr)",
                        breaks = c(1e5, 1e6, 1e7, 1e8),
                        labels = c("100 k", "1 M", "10 M", "100 M"),
                        guide = guide_legend(direction = "horizontal",
                                             title.position = "top", 
                                             label.position = "top",
                                             nrow = 1)) +
  
  scale_color_manual(values = c("-1" = "darkblue","1" = "darkred"),
                    name = "Demand change\n(positive / negative)")+
  
  coord_sf(xlim = c(min(nco[, 1]), max(nco[, 1])),
           ylim = c(min(nco[, 2]), max(nco[, 2]))) +
  labs(x = "Longitude", y = "Latitude") +
  map_theme(tsize_m, osize_m) +
  ggtitle(titleX_wu)

ggsave(paste0(plotDir, "/Fig3_mapF_wateruse_", haz, "_FINAL.jpg"),
       fig_wu, width = 23, height = 20, units = "cm", dpi = 1000)



#Suorva Dam example
suorva=5701511
SuorvaC=DataC[which(DataC$outl2==suorva),]
SuorvaL=DataL[which(DataL$outl2==suorva),]
SuorvaR=DataR[which(DataR$outl2==suorva),]
SuorvaC=DataC[which(DataC$outl2==suorva),]
plot(as.numeric(SuorvaR[,c(16:84)]))
# =============================================================
# 5  COLLECT AGGREGATED RESULTS  (feeds Section 8 in main script)
# =============================================================

paggC <- Dchangelist[[1]]$PagD; paggC$driver <- "Clim"   # climate
paggL <- Dchangelist[[4]]$PagD; paggL$driver <- "Lu"     # land use
paggR <- Dchangelist[[3]]$PagD; paggR$driver <- "Res"    # reservoirs
paggW <- Dchangelist[[5]]$PagD; paggW$driver <- "Wu"     # water use
pointsAD <- rbind(paggC, paggL, paggR, paggW)







  ## Collect aggregated spatial results
paggC <- Dchangelist[[1]]$PagD; paggC$driver <- "Clim"
paggL <- Dchangelist[[2]]$PagD; paggL$driver <- "Lu"
paggR <- Dchangelist[[3]]$PagD; paggR$driver <- "Res"
paggW <- Dchangelist[[4]]$PagD; paggW$driver <- "Wu"
pointsAD <- rbind(paggC, paggL, paggR, paggW)


## ── Figure 4 ─────────────────────────────────────────────────
## Dominant driver map – pixel level
fig4 <- ggplot(basemap) +
  geom_sf(fill="gray95", color="gray10", size=0.5) +
  geom_sf(data=nutplot,   aes(fill=factor(maxcol), geometry=geometry),
          color="transparent", alpha=.6, size=0.25, stroke=0, shape=15) +
  geom_sf(data=points_pix, aes(col=factor(maxcol), geometry=geometry, size=upa),
          alpha=.9, stroke=0, shape=15) +
  coord_sf(xlim=c(min(nco[,1]),max(nco[,1])), ylim=c(min(nco[,2]),max(nco[,2]))) +
  scale_colour_manual(values=colorz_drv, name="Largest change driver", labels=lab_drv) +
  scale_fill_manual(values=colorz_drv,   name="Largest change driver", labels=lab_drv) +
  scale_size(range=c(0.08,0.4), trans="sqrt", guide="none") +
  labs(x="Longitude", y="Latitude") +
  guides(colour=guide_legend(override.aes=list(size=10))) +
  map_theme(14, 12)

ggsave(paste0(plotDir, "/Fig4_dominantDriver_pix_", haz, ".jpg"),
       fig4, width=23, height=20, units="cm", dpi=1000)


# =============================================================
# 7  SUPPLEMENTARY FIGURES  (only run if RUN_SUPP == TRUE)
# =============================================================

if (RUN_SUPP) {
  
  # =============================================================
  #  Barplot – % of pixels with |ΔRL| > 5%, ranked by driver
  # =============================================================
  
  # ── 1. Total number of valid pixels (use DataC as reference) ─
  n_total <- sum(!is.na(DataC$Y2015))
  
  # ── 2. Count pixels per driver & sign ────────────────────────
  count_pixels <- function(Data, driver_name) {
    Data %>%
      filter(!is.na(Y2015)) %>%
      mutate(category = case_when(
        Y2015 >  5 ~ "increase",
        Y2015 < -5 ~ "decrease",
        TRUE       ~ NA_character_
      )) %>%
      filter(!is.na(category)) %>%
      group_by(category) %>%
      summarise(n = n(), .groups = "drop") %>%
      mutate(
        driver = driver_name,
        pct    = n / n_total * 100      # % of total pixel count
      )
  }
  
  pixel_counts <- bind_rows(
    count_pixels(DataC, "Climate"),
    count_pixels(DataL, "Land use"),
    count_pixels(DataR, "Reservoirs"),
    count_pixels(DataW, "Water demand")
  )
  
  # ── 3. Rank drivers by total % affected (top of plot = highest) ──
  driver_order <- pixel_counts %>%
    group_by(driver) %>%
    summarise(total_pct = sum(pct)) %>%
    arrange(total_pct) %>%          # ascending so coord_flip puts highest at top
    pull(driver)
  
  pixel_counts <- pixel_counts %>%
    mutate(
      driver   = factor(driver, levels = driver_order),
      category = factor(category, levels = c("decrease", "increase"))  # stack order
    )
  
  # ── 4. Compute label positions (midpoint of each segment) ────
  pixel_counts <- pixel_counts %>%
    group_by(driver) %>%
    arrange(driver, category) %>%
    mutate(
      cumulative = cumsum(pct),
      label_pos  = cumulative - pct / 2
    ) %>%
    ungroup()
  
  # ── 5. Total % label (end of full bar) ───────────────────────
  bar_totals <- pixel_counts %>%
    group_by(driver) %>%
    summarise(total_pct = sum(pct))
  
  # ── 6. Plot ───────────────────────────────────────────────────
  fill_cols <- c("decrease" = "tomato3", "increase" = "steelblue3")
  sign_labs  <- c("increase" = "ΔRL > +5%", "decrease" = "ΔRL < −5%")
  
  fig_bar <- ggplot(pixel_counts,
                    aes(x = driver, y = pct, fill = category)) +
    geom_col(width = 0.6, color = "grey20", linewidth = 0.3) +
    # Segment % labels inside bars
    # geom_text(aes(y = label_pos,
    #               label = paste0(round(pct, 1), "%")),
    #           size = 4, fontface = "bold", color = "white") +
    # Total % label at end of bar
    geom_text(data = bar_totals,
              aes(x = driver, y = total_pct,
                  label = paste0(round(total_pct, 1), "%")),
              inherit.aes = FALSE,
              hjust = -0.15, size = 4.5, fontface = "bold", color = "grey20") +
    scale_fill_manual(values = fill_cols, labels = sign_labs, name = NULL) +
    scale_y_continuous(
      name   = paste0("Share of river pixels with |Δ", haz, " RL| > 5%  (%)"),
      limits = c(0, max(bar_totals$total_pct) * 1.15),
      expand = c(0, 0)
    ) +
    scale_x_discrete(name = NULL) +
    coord_flip() +
    theme(
      axis.title.x       = element_text(size = 14, face = "bold"),
      axis.text          = element_text(size = 13),
      axis.text.y        = element_text(face = "bold"),
      axis.ticks.y       = element_blank(),
      panel.background   = element_rect(fill = "white"),
      panel.grid.major.x = element_line(color = "grey85", linetype = "dashed"),
      panel.grid.major.y = element_blank(),
      panel.border       = element_rect(fill = NA, color = "black", linewidth = 1),
      legend.position    = "bottom",
      legend.text        = element_text(size = 12),
      legend.key.size    = unit(0.7, "cm")
    ) +
    ggtitle(paste0("Share of river pixels with |Δ", haz,
                   " return level| > 5%  (1955–2015)"))
  fig_bar
  
  ggsave(paste0(plotDir, "/FigX_pixel_pct_barplot_", haz, "_FINAL.jpg"),
         fig_bar, width = 20, height = 14, units = "cm", dpi = 1000)

  ## ── S1 · IRES map (drought only) ─────────────────────────
  if (haz == "Drought") {
    IRpoints_sf <- inner_join(IRES_comb, UpArea, by = c("catlist" = "outl2"))
    length(which(IRES_comb$gen_IR==1))/length(IRES_comb$gen_IR)
    colIR   <- c("0"="royalblue","1"="lightblue","2"="orangered","3"="tomato","4"="purple")
    pts_IR  <- st_transform(st_as_sf(IRpoints_sf, coords=c("Var1.x","Var2.x"), crs=4326), crs=3035)

    figS1 <- ggplot(basemap) +
      geom_sf(fill="gray95", color="gray10", size=0.5) +
      geom_sf(data=pts_IR, aes(col=factor(gen_IR), geometry=geometry, size=upa),
              alpha=.9, stroke=0, shape=15) +
      coord_sf(xlim=c(min(nco[,1]),max(nco[,1])), ylim=c(min(nco[,2]),max(nco[,2]))) +
      scale_colour_manual(values=colIR, name="IR",
                          labels=c("0"="perennial","1"="casi-perennial","2"="IRES")) +
      scale_size(range=c(0.08,0.4), trans="sqrt", guide="none") +
      labs(x="Longitude", y="Latitude") +
      guides(colour=guide_legend(override.aes=list(size=10))) +
      map_theme()
    ggsave(paste0(plotDir, "/FigS1_IRES_", haz, ".jpg"), figS1,
           width=20, height=20, units="cm", dpi=1000)
  }

  ## ── S2 · Mean change in hazard intensity with CI ─────────
  RlevErrtH <- c()
  for (yr in 1952:2020) {
    ParamSpec <- as.data.frame(ParamsflSCF[ParamsflSCF$Year == yr, ])
    ParamSpec <- ParamSpec[-match(rmp2, ParamSpec$catchment), ]
    Rlev <- calculate_return_levels(ParamSpec, ci = 2)
    mapu <- na.omit(match(ParamSpec$catchment, UpArea$outl2))
    ParamSpec$upa <- UpArea$upa[mapu]
    merd=which.max(Rlev$returnLevels)
    ParamSpec[merd,]
    RlevErrtH <- rbind(c(mean(Rlev$returnLevels, na.rm=TRUE),
                         mean(Rlev$returnLevelErr, na.rm=TRUE), yr), RlevErrtH)
  }
  RlevErrtH <- data.frame(RlevErrtH)
  if (haz == "Drought") RlevErrtH$X1 <- -RlevErrtH$X1
  idb <- which(RlevErrtH$X3 == 1955)
  RlevErrtH$X4 <- RlevErrtH$X1 / RlevErrtH$X1[idb] * 100
  RlevErrtH$X5 <- RlevErrtH$X2 / RlevErrtH$X1[idb] * 100

  figS2 <- ggplot(RlevErrtH) +
    geom_line(aes(x=X3, y=X4), lwd=2) +
    geom_ribbon(aes(x=X3, ymin=X4-X5, ymax=X4+X5), fill="blue", alpha=0.2) +
    scale_y_continuous(name="Change in 10yRL", breaks=seq(80,120,2),
                       trans=scales::modulus_trans(.6)) +
    scale_x_continuous(breaks=seq(1950,2020,5), labels=seq(1950,2020,5),
                       name="Years", limits=c(1955,2020),
                       minor_breaks=seq(1955,2005,10), expand=c(.001,.001)) +
    bxp_theme() + ggtitle("Europe")

  figS2
  ggsave(paste0(plotDir, "/FigS2_meanChangeCI_", haz, "FINAL.jpg"),
         figS2, width=30, height=20, units="cm", dpi=1000)

  ## ── S3 · Relative RL error in 1955 ───────────────────────
  paraU <- as.data.frame(ParamsflSCF[ParamsflSCF$Year == 1955, ])
  paraU <- paraU[-match(rmp2, paraU$catchment), ]
  Paramf <- ParamsflSCF[ParamsflSCF$Year == 2015, ]
  Paramf <- Paramf[-match(rmp2, Paramf$catchment), ]

  X0    <- Paramf$nPeaks / 70
  RlevErrf <- calculate_return_levels(Paramf, X0 * 100, ci = 1); RlevErrf$year <- 2015
  RlevErri <- calculate_return_levels(paraU, ci = 1); RlevErri$year <- 1955

  if (haz == "Drought") {
    for (obj in c("RlevErri","RlevErrf")) {
      tmp <- get(obj)
      tmp$returnLevels <- -tmp$returnLevels
      tmp$returnLevels[tmp$returnLevels < 0 & !is.na(tmp$returnLevels)] <- 
        tmp$returnLevelErr[tmp$returnLevels < 0 & !is.na(tmp$returnLevels)]
      #RlevErri$returnLevels[which(RlevErri$returnLevels<0)]=RlevErri$returnLevelErr[which(RlevErri$returnLevels<0)]
      
      assign(obj, tmp)
    }
  }
  RlevErri$catchment <- paraU$catchment
  RlevErri$relErr    <- (RlevErri$returnLevelErr / RlevErri$returnLevels + 1e-4) * 100
  RlevErri$rlf       <- RlevErrf$returnLevels
  RlevErri$ErrVsCh   <- RlevErri$rlf - RlevErri$returnLevelErr
  RlevErri$S2n       <- ifelse(is.nan(RlevErri$ErrVsCh), "Unstable",
                          ifelse(RlevErri$ErrVsCh >= 0, "Change > Err", "Err > Change"))
  
  length(which(RlevErri$S2n=="Unstable"))/length(RlevErri$S2n)*100
  ParamPlot <- inner_join(RlevErri, UpArea, by = c("catchment" = "outl2"))
  Paraplot  <- st_transform(st_as_sf(ParamPlot, coords=c("Var1.x","Var2.x"), crs=4326), crs=3035)

  paletS <- hcl.colors(11, "YlGnBu", rev=TRUE)
  figS3 <- ggplot(basemap) +
    geom_sf(fill="gray90", color="darkgrey", size=0.5) +
    geom_sf(data=Paraplot, aes(col=relErr, geometry=geometry, size=upa),
            alpha=1, stroke=0, shape=15) +
    scale_size(range=c(0.1,0.5), trans="sqrt", guide="none") +
    coord_sf(xlim=c(min(nco[,1]),max(nco[,1])), ylim=c(min(nco[,2]),max(nco[,2]))) +
    scale_color_gradientn(colors=paletS, breaks=seq(0,200,20), limits=c(0,200),
                          oob=scales::squish, na.value="white",
                          name="RL(1955) RelErr (%)") +
    labs(x="Longitude", y="Latitude") +
    guides(colour=guide_colourbar(barheight=16, barwidth=1.5)) +
    map_theme()
  ggsave(paste0(plotDir, "/FigS3_RLerror_", haz, "FINAL.jpg"),
         figS3, width=23, height=20, units="cm", dpi=1000)

  ## ── S4 · Error class map ─────────────────────────────────
  figS4 <- ggplot(basemap) +
    geom_sf(fill="gray90", color="darkgrey", size=0.5) +
    geom_sf(data=Paraplot, aes(col=S2n, geometry=geometry, size=upa),
            alpha=1, stroke=0, shape=15) +
    scale_size(range=c(0.1,0.5), trans="sqrt", guide="none") +
    coord_sf(xlim=c(min(nco[,1]),max(nco[,1])), ylim=c(min(nco[,2]),max(nco[,2]))) +
    scale_color_manual(values=c("Unstable"="purple","Err > Change"="darkred",
                                "Change > Err"="royalblue"), name="Error classes",guide = guide_colorbar(
                                  direction = "horizontal", 
                                  title.position = "top", 
                                  barwidth = 10         # You can control the width of the bar here
                                )) +
    labs(x="Longitude", y="Latitude") +
    guides(colour=guide_legend(override.aes=list(size=10))) +
    map_theme(pos="bottom")
  ggsave(paste0(plotDir, "/FigS4_errorClass_", haz, "FINAL.jpg"),
         figS4, width=23, height=20, units="cm", dpi=1000)

  ## ── S5 · Ordered change at HER level (climate only) ──────
  pointP <- Dchangelist[[1]]$PagD
  pointP <- pointP[match(unique(pointP$CODEB), pointP$CODEB), ]
  pointP <- pointP[order(pointP$Rchange_rel.mean), ]
  pointP$id <- seq_len(nrow(pointP))
  manualcol <- c("-2"="#A51122","-1"="#F1C363","1"="#ACD2BB","2"="#324DA0")
  manualab  <- c("sig. decrease","decrease","increase","sig. increase")

  figS5 <- ggplot() +
    coord_cartesian(ylim = c(-200, 200)) +
    geom_hline(yintercept = 0, lwd = 1, col = "black") +
    geom_segment(data=pointP, aes(x=id, xend=id,
                                  y=Rchange_rel.q1.2.5., yend=Rchange_rel.q3.97.5.,
                                  color=factor(change), size=Rchange_rel.len), alpha=.99) +
    scale_color_manual(values=manualcol, breaks=manualab, labels=manualab, name="") +
    geom_point(data=pointP, aes(x=id, y=Rchange_rel.mean),
               pch=21, fill="white", colour="gray3", size=3, stroke=1) +
    geom_text(data=pointP, aes(x=id, y=Rchange_rel.mean, label=CODEB),
              size=1.5, color="black", fontface="bold") +
    scale_size(range=c(0.8,4), trans="sqrt", guide="none") +
    scale_y_continuous(limits=c(-5000,5000),
                       breaks=c(-200,-100,-50,-10,0,10,50,100,200),
                       name="Change (%)", trans=scales::modulus_trans(.3)) +
    scale_x_continuous(name="HER", expand=c(.01,.01), breaks=c(-100,200)) +
    guides(colour=guide_legend(override.aes=list(size=10))) +
    theme(axis.title=element_text(size=20, face="bold", color="black"),
          axis.text=element_text(size=18, color="black"),
          panel.background=element_rect(fill="white", colour="white"),
          panel.grid=element_blank(),
          panel.border=element_rect(linetype="solid", fill=NA, colour="black", linewidth=2),
          legend.title=element_text(size=20), legend.text=element_text(size=18, color="black"),
          legend.position="none",
          legend.key=element_rect(fill="transparent", colour="transparent"),
          legend.key.size=unit(.8,"cm"))
  ggsave(paste0(plotDir, "/FigS5_orderedHR_", haz, ".jpg"),
         figS5, width=40, height=8, units="cm", dpi=400)

  ## ── S6 · Boxplot by biogeographic region ─────────────────
  bio_names <- unique(biogeo$code)[c(1,3,4,6,7,9,11)]
  br_bgr    <- c(seq(-200,-10,10), seq(-5,5,5), seq(10,200,10))
  xlabs_bgr <- seq(1950, 2010, 10)

  for (bn in bio_names) {
    trtF_bn <- rbind(
      TdataClim$BgData[TdataClim$BgData$loc == bn, ],
      TdataLuse$BgData[TdataLuse$BgData$loc == bn, ],
      TdataRes$BgData[ TdataRes$BgData$loc  == bn, ],
      TdataWuse$BgData[TdataWuse$BgData$loc == bn, ]
    )
    trtF_bn$year   <- rep(seq(1950,2010,10), 4)
    trtF_bn$driver <- rep(c("Clim","LUC","Res","WU"), each = 7)
    names(trtF_bn)[c(3,5,6,7,8,9)] <- c("changeC","med","cq1","cq2","w1","w2")

    figS6_bn <- ggplot(trtF_bn) +
      geom_linerange(aes(x=year, ymin=w1, ymax=w2, color=factor(driver), group=factor(driver)),
                     position=position_dodge2(width=9), lwd=1, alpha=0.8) +
      scale_color_manual(values=colorn, name="Drivers", labels=clabels) +
      new_scale_color() +
      geom_rect(aes(xmin=year-4.5, xmax=year+4.5, ymin=cq1, ymax=cq2,
                    fill=factor(driver), group=factor(driver)),
                alpha=0.5, position=position_dodge(width=9)) +
      geom_point(aes(x=year, y=changeC, color=factor(driver), group=factor(driver)),
                 position=position_dodge(width=9), size=3) +
      geom_rect(aes(xmin=year-4.5, xmax=year+4.5, ymin=med-1e-1, ymax=med+1e-1,
                    fill=factor(driver), group=factor(driver)),
                alpha=1, position=position_dodge(width=9)) +
      scale_y_continuous(name="Mean change (% of 10Y RL)", breaks=br_bgr,
                         trans=scales::modulus_trans(.6)) +
      scale_x_continuous(breaks=xlabs_bgr, labels=xlabs_bgr, name="Decades",
                         minor_breaks=seq(1955,2005,10), expand=c(.01,.01)) +
      scale_fill_manual(values=colorn, name="Drivers", labels=clabels) +
      scale_color_manual(values=colorz, name="Drivers", labels=clabels) +
      guides(color=guide_legend(override.aes=list(color=colorz))) +
      bxp_theme() + ggtitle(bn)

    ggsave(paste0(plotDir, "/FigS6_bxp_BGR_", bn, "_", haz, "_23.jpg"),
           figS6_bn, width=30, height=20, units="cm", dpi=1000)
  }

  ## ── S7 · Boxplot by catchment size class ─────────────────
  ClimUpag <- UpATrendData(DataC)
  LuseUpag <- UpATrendData(DataL)
  ResUpag  <- UpATrendData(DataR)
  WdemUpag <- UpATrendData(DataW)
  UpaF     <- rbind(ClimUpag, LuseUpag, ResUpag, WdemUpag)
  UpaF$driver <- rep(c("Clim","LUC","Res","WU"), each=5)
  UpaF <- UpaF[, -1]
  names(UpaF)[c(2,4,5,6,7,8)] <- c("changeC","med","cq1","cq2","w1","w2")
  br_upa <- c(seq(-300,-20,20), seq(-20,20,10), seq(20,300,20))

  figS7 <- ggplot(UpaF) +
    geom_linerange(aes(x=loc, ymin=w1, ymax=w2, color=factor(driver), group=factor(driver)),
                   position=position_dodge2(width=.9), lwd=1, alpha=0.8) +
    scale_color_manual(values=colorn, name="Drivers", labels=clabels) +
    new_scale_color() +
    geom_rect(aes(xmin=loc-.45, xmax=loc+.45, ymin=cq1, ymax=cq2,
                  fill=factor(driver), group=factor(driver)),
              alpha=0.5, position=position_dodge(width=.9)) +
    geom_point(aes(x=loc, y=changeC, color=factor(driver), group=factor(driver)),
               position=position_dodge(width=.9), size=3) +
    geom_rect(aes(xmin=loc-.45, xmax=loc+.45, ymin=med-1e-1, ymax=med+1e-1,
                  fill=factor(driver), group=factor(driver)),
              alpha=1, position=position_dodge(width=.9)) +
    scale_y_continuous(name="Mean change (% of mean 10yRL)", breaks=br_upa,
                       trans=scales::modulus_trans(.6)) +
    scale_x_continuous(breaks=1:5,
                       labels=c("100–200","200–500","500–1000","1000–10 000",">10 000"),
                       name="Catchment Area (km²)") +
    scale_fill_manual(values=colorn, name="Drivers", labels=clabels) +
    scale_color_manual(values=colorz, name="Drivers", labels=clabels) +
    guides(color=guide_legend(override.aes=list(color=colorz))) +
    bxp_theme() + ggtitle("Europe")

  ggsave(paste0(plotDir, "/FigS7_bxp_UpA_", haz, "_23.jpg"),
         figS7, width=30, height=20, units="cm", dpi=1000)

  ## ── S8 · Dominant driver at HER level (choropleth) ───────
  figS8 <- ggplot(basemap) +
    geom_sf(fill="gray95") + geom_sf(fill=NA, color="grey") +
    geom_sf(data=nutplot, aes(fill=factor(maxcol), geometry=geometry),
            color="transparent", alpha=1, size=0.25, stroke=0, shape=15) +
    coord_sf(xlim=c(min(nco[,1]),max(nco[,1])), ylim=c(min(nco[,2]),max(nco[,2]))) +
    scale_fill_manual(values=colorz_drv, name=" ") +
    labs(x="Longitude", y="Latitude") +
    map_theme()
  ggsave(paste0(plotDir, "/FigS8_dominantDriver_HER_", haz, ".jpg"),
         figS8, width=23, height=20, units="cm", dpi=1000)


} 


# =============================================================
# 8  SAVE OUTPUTS
# =============================================================

DataT$driver <- "Total"; DataC$driver <- "Clim"
DataR$driver <- "Reservoirs"; DataL$driver <- "Landuse"; DataW$driver <- "Wateruse"
Alltrend <- rbind(DataC, DataL, DataR, DataW, DataT)

trendTot$driver <- "Total"; trendClim$driver <- "Clim"
trendRes$driver <- "Reservoirs"; trendSoc$driver <- "Landuse"; trendWu$driver <- "Wateruse"
trendRegio <- rbind(trendClim, trendSoc, trendRes, trendWu, trendTot)

outfile <- if (haz == "Flood") {
  Output_fl_year <- list(TrendPix=Alltrend, TrendRegio=trendRegio,
                          Out2020=pointsAD, DataI=DataI)
  save(Output_fl_year, file=paste0( paste0(hydroDir, "/Flood/outputs_flood_year_relxHR_FINAL.Rdata")))
} else {
  Output_dr_nonfrost <- list(TrendPix=Alltrend, TrendRegio=trendRegio,
                              Out2020=pointsAD, DataI=DataI)
  save(Output_dr_nonfrost, file=paste0paste0(hydroDir, "/Drought/outputs_drought_nonfrost_relxHR_FINAL.Rdata"))
}
