## ============================================================
##  04_CHEX_Plot_bivariateF.R  –  reorganized
##  Main figures first, supplementary below RUN_SUPP flag
## ============================================================

setwd(dirname(rstudioapi::getSourceEditorContext()$path))
source("functions_trends.R")
library(ggsankey)
RUN_SUPP <- FALSE   # set TRUE to produce all supplementary figures


# =============================================================
# 0  PATHS & GLOBAL SETTINGS
# =============================================================

hydroDir <- "D:/tilloal/Documents/LFRuns_utils/ChangingHydroExtremes/data"
plotDir  <- "D:/tilloal/Documents/LFRuns_utils/ChangingHydroExtremes/plots"

tsize <- 16; osize <- 12

# Bivariate colour palettes
colorp <- c(
  "1-1" = "#dd6a40", "2-1" = "#d9926a", "3-1" = "#DEB887", "4-1" = "#FFD39B",
  "1-2" = "#a36229", "2-2" = "#FFFFFF", "3-2" = "#FFFFFF", "4-2" = "#9cc4d2",
  "1-3" = "#635929", "2-3" = "#FFFFFF", "3-3" = "#FFFFFF", "4-3" = "#5fb2d1",
  "1-4" = "#174f28", "2-4" = "#166d68", "3-4" = "#16869e", "4-4" = "#169dd0"
)
colors <- c(
  "1-1" = "#dd6a29", "2-1" = "#d9926a", "3-1" = "#d6b3a0", "4-1" = "#d3d3d3",
  "1-2" = "#a36229", "2-2" = "#a08769", "3-2" = "#9ea69f", "4-2" = "#9cc4d2",
  "1-3" = "#635929", "2-3" = "#617b69", "3-3" = "#60979f", "4-3" = "#5fb2d1",
  "1-4" = "#174f28", "2-4" = "#166d68", "3-4" = "#16869e", "4-4" = "#169dd0"
)
colord <- c(
  "1-1"=12,"2-1"=11,"3-1"=10,"4-1"=9,
  "1-2"=13,"2-2"=4, "3-2"=3, "4-2"=8,
  "1-3"=14,"2-3"=1, "3-3"=2, "4-3"=7,
  "1-4"=15,"2-4"=16,"3-4"=5, "4-4"=6
)
loscolors <- c("Accelerating"="#174f28","Drying"="#dd6a29",
               "Stable"="gray60","Wetting"="#169dd0","Decelerating"="burlywood")
colorn    <- c("WaterDemand"="limegreen","Reservoirs"="tomato4",
               "Landuse"="orange","Climate"="royalblue")
colorz_haz <- c("Flood"="darkblue","Drought"="darkorange")
reg_labels <- c("Alpine"="ALP","Atlantic"="ATL","Boreal"="BOR",
                "Continental"="CON","Mediterranean"="MED")
colx_bgr   <- c("Alpine"="deeppink","Atlantic"="forestgreen","Boreal"="darkviolet",
                "Continental"="orange3","Mediterranean"="gold")

# Classification thresholds
breaker1 <- 0; breaker2 <- 5

# Quadrant annotations (reused across many plots)
quadrant_annots <- list(
  annotate("text", x= 35, y= 45, label="Wetting",      color="#169dd0", size=6, fontface="bold"),
  annotate("text", x=-35, y= 45, label="Accelerating",  color="#174f28", size=6, fontface="bold"),
  annotate("text", x=-35, y=-45, label="Drying",        color="#dd6a29", size=6, fontface="bold"),
  annotate("text", x= 35, y=-45, label="Decelerating",  color="burlywood",size=6,fontface="bold")
)

# Shared map theme
map_theme_bi <- function(ts = tsize, os = osize) {
  theme(
    axis.title       = element_text(size = ts),
    panel.background = element_rect(fill = "aliceblue", colour = "grey1"),
    panel.border     = element_rect(linetype = "solid", fill = NA, colour = "black"),
    legend.title     = element_text(size = ts),
    legend.text      = element_text(size = os),
    legend.position  = "bottom",
    panel.grid.major = element_line(colour = "grey70"),
    panel.grid.minor = element_line(colour = "grey90"),
    legend.key       = element_rect(fill = "transparent", colour = "transparent"),
    legend.key.size  = unit(.8, "cm")
  )
}

bv_scatter_theme <- function() {
  theme(
    axis.title   = element_text(size = 18, face = "bold"),
    title        = element_text(size = 22, face = "bold"),
    axis.text    = element_text(size = 16, face = "bold"),
    panel.background = element_rect(fill = "white", colour = "white"),
    panel.grid   = element_blank(),
    panel.border = element_rect(linetype = "solid", fill = NA, colour = "black"),
    legend.title = element_text(size = 20, face = "bold"),
    legend.text  = element_text(size = 16),
    legend.position = "bottom",
    panel.grid.major = element_line(colour = "grey60"),
    panel.grid.minor = element_line(colour = "grey70", linetype = "dashed")
  )
}


# =============================================================
# 1  DATA LOADING
# =============================================================

## 1.1  Outlets ------------------------------------------------
if (!exists("outf")) {
  rspace <- read.csv(paste0(hydroDir, "/subspace_efas.csv"))[, -1]
  outf   <- c()
  for (Nsq in 1:88) {
    nrspace  <- rspace[Nsq, ]
    outhybas <- outletopen(hydroDir, "GeoData/efas_rnet_100km_01min", nrspace)
    Idstart  <- as.numeric(Nsq) * 10000
    Idstart2 <- as.numeric(Nsq) * 100000
    if (length(outhybas$outlets) > 0) {
      outhybas$outlets <- seq(Idstart + 1, Idstart + length(outhybas$outlets))
      outhybas$outl2   <- seq(Idstart2 + 1, Idstart2 + length(outhybas$outlets))
      outhybas$latlong <- paste(round(outhybas$Var1, 4), round(outhybas$Var2, 4), sep = " ")
      outf <- rbind(outf, outhybas)
    }
  }
}

## 1.2  Flood & drought outputs --------------------------------
load(file = paste0(hydroDir, "/Flood/outputs_flood_year_relxHR_FINAL.Rdata"))
load(file = paste0(hydroDir, "/Drought/outputs_drought_nonfrost_relxHR_FINAL.Rdata"))

FloodTrends   <- Output_fl_year$TrendRegio
DroughtTrends <- Output_dr_nonfrost$TrendRegio
FloodTrendsP  <- Output_fl_year$TrendPix
DroughtTrendsP <- Output_dr_nonfrost$TrendPix

driver <- unique(DroughtTrends$driver)   # [1]Clim [2]Lu [3]Res [4]Wu [5]Total

## 1.3  Spatial layers -----------------------------------------
biogeo        <- read_sf(dsn = paste0(hydroDir, "/GeoData/eea_3035_biogeo-regions_2016/BiogeoRegions2016_wag84.shp"))
biogeof       <- fortify(biogeo); st_geometry(biogeof) <- NULL
biogeoregions <- raster(paste0(hydroDir, "/GeoData/eea_3035_biogeo-regions_2016/Biogeo_rasterized_wsg84.tif"))
Gbiogeoregions <- as.data.frame(biogeoregions, xy = TRUE)
biogeomatch   <- inner_join(biogeof, Gbiogeoregions,
                            by = c("PK_UID" = "Biogeo_rasterized_wsg84"))
biogeomatch$latlong  <- paste(round(biogeomatch$x, 4), round(biogeomatch$y, 4), sep = " ")
biogeo_rivers <- right_join(biogeomatch, outf, by = "latlong")

hybas07  <- read_sf(dsn = paste0(hydroDir, "/GeoData/HYBAS07/hybas_eu_lev07_v1c.shp"))
hybasf7  <- fortify(hybas07)
Catamere07 <- inner_join(hybasf7,
                          read.csv(paste0(hydroDir, "/GeoData/HYBAS07/from_hybas_eu_onlyid.csv"),
                                   encoding = "UTF-8", stringsAsFactors = FALSE),
                          by = "HYBAS_ID")
Catamere07$llcoord <- paste(round(Catamere07$POINT_X, 4),
                             round(Catamere07$POINT_Y, 4), sep = " ")
cst7 <- right_join(Catamere07, outf, by = c("llcoord" = "latlong"))
GNF  <- cst7; st_geometry(GNF) <- NULL

GHshpp   <- read_sf(dsn = "Z:/ClimateRun4/nahaUsers/tilloal/HydroRegions/her_all_adjusted_cor.shp")
HydroRsf <- fortify(GHshpp)

GridHR   <- raster(paste0(hydroDir, "/GeoData/HER/HydroRegions_raster_WGS84.tif"))
GHR      <- as.data.frame(GridHR, xy = TRUE)
GHR      <- GHR[!is.na(GHR[, 3]), ]
GHR$llcoord <- paste(round(GHR$x, 4), round(GHR$y, 4), sep = " ")
GHR_riv  <- inner_join(GHR, outf, by = c("llcoord" = "latlong"))

outll    <- outletopen(hydroDir, "GeoData/efas_rnet_100km_01min")
cord.dec <- SpatialPoints(outll[, c(2, 3)], proj4string = CRS("+proj=longlat"))
cord.UTM <- spTransform(cord.dec, CRS("+init=epsg:3035"))
nco      <- cord.UTM@coords
world    <- ne_countries(scale = "medium", returnclass = "sf")
basemap  <- st_transform(world, crs = 3035)
biobase  <- st_transform(biogeo, crs = 3035)


# =============================================================
# 2  HELPER FUNCTIONS
# =============================================================

## Aggregate trend data to decade-level box stats
agg_decade <- function(Trends, drv) {
  Totrend <- Trends[Trends$driver == drv, ]
  colnames(Totrend)[2:71] <- 1951:2020
  tr <- suppressWarnings(
    melt(Totrend, id.vars = "HydroR", variable.name = "variable", value.name = "value")
  )
  tr$yr    <- as.numeric(as.character(tr$variable))
  tr$value <- as.numeric(tr$value)

  # Continuous time series (for ribbon plot)
  cts <- do.call(data.frame,
    aggregate(list(value = tr$value), by = list(yr = tr$yr),
              FUN = function(x) c(mean = mean(x, na.rm=TRUE),
                                  med  = median(x, na.rm=TRUE),
                                  w1   = quantile(x, 0.025, na.rm=TRUE),
                                  w2   = quantile(x, 0.975, na.rm=TRUE))))

  # Decade summaries
  tr_dec  <- tr[tr$yr %in% seq(1955, 2015, 10), ]
  dec <- do.call(data.frame,
    aggregate(list(value = tr_dec$value), by = list(yr = tr_dec$yr),
              FUN = function(x) c(mean = mean(x, na.rm=TRUE),
                                  med  = median(x, na.rm=TRUE),
                                  ql   = quantile(x, 0.25, na.rm=TRUE),
                                  qh   = quantile(x, 0.75, na.rm=TRUE),
                                  w1   = quantile(x, 0.025, na.rm=TRUE),
                                  w2   = quantile(x, 0.975, na.rm=TRUE))))
  names(dec)[c(2,3,4,5,6,7)] <- c("changeC","med","cq1","cq2","w1","w2")
  list(cts = cts, dec = dec)
}

## Assign bivariate class from two numeric vectors
make_bi_class <- function(x, y, b1 = breaker1, b2 = breaker2) {
  classify <- function(v) {
    cl <- rep(NA_integer_, length(v))
    cl[v <= -b2]              <- 1L
    cl[v > -b2 & v < b1]     <- 2L
    cl[v >= b1 & v < b2]     <- 3L
    cl[v >= b2]               <- 4L
    cl
  }
  paste(classify(y), classify(x), sep = "-")
}


## Assign trajectory category from bi_class string
assign_trcat <- function(bi) {
  tr <- rep(NA, length(bi))
  tr[bi %in% c("2-2","3-2","2-3","3-3")]                <- "Stable"
  tr[bi %in% c("1-1","1-2","2-1")]                       <- "Drying"
  tr[bi %in% c("1-4","2-4","1-3")]                       <- "Accelerating"
  tr[bi %in% c("4-1","4-2","3-1")]                       <- "Decelerating"
  tr[bi %in% c("4-4","4-3","3-4")]                       <- "Wetting"
  tr[is.na(bi) | bi == "NA-NA"]                          <- NA
  tr
}

## HER-level spatial join for one driver
her_bidata <- function(FloodTrends, DroughtTrends, drv, GHshpp, HydroRsf) {
  flt <- FloodTrends[FloodTrends$driver == drv, ]
  drt <- DroughtTrends[DroughtTrends$driver == drv, ]
  flt$HydroR <- FloodTrends[FloodTrends$driver == driver[1], 1]
  drt$HydroR <- DroughtTrends[DroughtTrends$driver == driver[1], 1]

  FlP <- full_join(GHshpp, flt, by = c("CODEB" = "HydroR"))
  st_geometry(FlP) <- NULL
  Flpl <- st_transform(inner_join(HydroRsf, FlP, by = "Id"), crs = 3035)

  DrP <- full_join(GHshpp, drt, by = c("CODEB" = "HydroR"))
  st_geometry(DrP) <- NULL
  Drpl <- st_transform(inner_join(HydroRsf, DrP, by = "Id"), crs = 3035)

  Flpl$x <- Flpl$Rchange.Y2015
  Flpl$y <- Drpl$Rchange.Y2015
  Flpl$bi_class          <- make_bi_class(Flpl$x, Flpl$y)
  Flpl$combined_category <- Flpl$bi_class
  Flpl
}

## Pixel-level bivariate data for one driver
pix_bidata <- function(FloodPix, DroughtPix, base_sf) {
  db       <- base_sf
  db$x     <- FloodPix$Y2015
  db$y     <- DroughtPix$Y2015
  db$bi_class <- make_bi_class(db$x, db$y)
  db$bi_class[is.na(db$x) | is.na(db$y)] <- NA
  db
}

## Stacked barplot by biogeographic region
biogeo_barplot <- function(databipi2, colorp, tsize = 22, osize = 28) {
  databipi2$Biogeo_id[databipi2$Biogeo_id == "Pannonian"] <- "Continental"
  databipi3 <- databipi2[!is.na(databipi2$bi_class), ]

  RegioSize <- aggregate(list(val = databipi2$upa),
                         by = list(region = databipi2$Biogeo_id),
                         FUN = length)

  databipi2$bi_class[which(is.na(databipi2$bi_class))]="nocal"
  magg <- do.call(data.frame,
    aggregate(list(val = databipi2$upa),
              by = list(reg = databipi2$Biogeo_id, traj = databipi2$bi_class),
              FUN = length))
  magg <- magg[!magg$reg %in% c("Steppic","BlackSea","Arctic"), ]
  magg$trcat <- assign_trcat(magg$traj)

  # Ordering by colord
  colorp[which(colorp=="#999999")]="white"
  magg$ordf <- colord[match(magg$traj, names(colord))]
  magg      <- magg %>% mutate(traj = reorder(traj, ordf, FUN = mean))
  
  magg$rm=0
  magg$rm[magg$traj %in% c("2-2","3-2","2-3","3-3")] <- 1
  magg$rm[is.na(magg$trcat)] <-1
  # ── 5. Label positions (midpoint of each stacked segment) ─────
  magg <- magg %>%
    group_by(reg) %>%
    arrange(reg, traj, rm) %>%
    mutate(
      cum_pct   = sum(val),
    ) %>%
    ungroup()
  magg=magg[-which(magg$rm==1),]
  magg <- magg %>%
    group_by(reg) %>%
    arrange(reg, traj) %>%
    mutate(
      pct   = (val)/cum_pct,
    ) %>%
    ungroup()
  
  
  # bar_totals <- magg %>%
  #   group_by(reg) %>%
  #   summarise(total_pct =pct)
  
  # Total % at end of each bar
  bar_totals <- magg %>%
    group_by(reg) %>%
    summarise(total_pct = sum(val)/mean(cum_pct))
  
  # magg$rm=0
  # magg$rm[magg$traj %in% c("2-2","3-2","2-3","3-3")] <- 1
  # magg=magg[-which(magg$rm==1),]

  ggplot(magg, aes(x = reg, y = pct, fill = traj)) +
    geom_col(width = 0.7, color = "transparent")+
    # geom_bar(stat = "identity", position = "fill") +
    scale_fill_manual(values = colorp) +
    scale_x_discrete(name = "Biogeographic regions", labels = reg_labels, expand = c(0, 0)) +
    # scale_y_continuous(breaks = c(0,.25,.50,.75,1),
    #                    labels = c(0, 25, 50, 75, 100),
    #                    name = "Share of river pixels (%)", expand = c(0.005, 0.005)) +
    scale_y_continuous(
      name   = "Share of river network (%)",
      breaks = c(0,.25,.50,.75,1),
      labels = c(0, 25, 50, 75, 100),
      limits = c(0, 1.2),
      expand = c(0, 0)
    ) +
    geom_text(data = bar_totals,
              aes(x = reg, y = total_pct+.05,
                  label = paste0(round(total_pct*100, 0), "%")),
              inherit.aes = FALSE,
              hjust = 0.5, size = 4.5, fontface = "bold", color = "grey20") +
    geom_vline(xintercept = 0.5,
               color = "black", linewidth = 1)+
    geom_hline(yintercept = 0,
               color = "black", linewidth = 1)+
    theme(
      axis.title      = element_text(size = 16, face = "bold"),
      axis.text          = element_text(size = 15),
      axis.text.y        = element_text(face = "bold"),
      # axis.text.x        = element_blank(),
      axis.ticks         = element_blank(),
      panel.background   = element_rect(fill = "white"),
      # panel.grid.major.x = element_line(color = "grey85", linetype = "dashed"),
      
      panel.grid.major = element_blank(),
      panel.border       = element_blank(), 
      legend.position    = "none",
      # legend.text        = element_text(size = 11),
      # legend.title       = element_text(size = 12, face = "bold"),
      legend.key.size    = unit(0.7, "cm")
    ) 
}
## Sankey diagram
make_sankey <- function(class_list, stage_names, loscolors, ltot) {
  d  <- data.frame(do.call(cbind, class_list))
  d<-d[-which(is.na(d$X1)),]
  names(d) <- stage_names
  df <- d %>% make_long(!!!syms(stage_names))

  brk <- na.omit(unique(df$node))
  df$node      <- factor(df$node,      levels = brk)
  df$next_node <- factor(df$next_node, levels = brk)

  reagg <- df %>% group_by(node, x) %>% tally()
  df2   <- merge(df, reagg, by.x = "node", by.y = "node", all.x = FALSE)
  df2   <- df2[df2$x.x == df2$x.y, ]
  df2$np <- round(df2$n / ltot * 100)

  ggplot(df2, aes(x = x.x, next_x = next_x, node = node, next_node = next_node,
                  fill = factor(node),
                  label = paste0(node, "\n", np, "%"))) +
    geom_sankey(flow.alpha = .5, node.color = "black", show.legend = TRUE) +
    geom_sankey_label(size = 3, color = "black", fill = "white") +
    theme_sankey(base_size = 18) +
    theme(legend.position = "none",
          axis.title = element_blank(), axis.text.y = element_blank(),
          axis.ticks = element_blank(), panel.grid = element_blank()) +
    scale_fill_manual(values = loscolors)
}

## HER-level aggregation (drought & flood separately, then join)
agg_her <- function(data_hr, xy_col, group_col = "HER") {
  do.call(data.frame,
    aggregate(list(val = data_hr[[xy_col]]),
              by = list(HR = data_hr[[group_col]]),
              FUN = function(x) c(mean = mean(x,na.rm=T),
                                  med  = median(x,na.rm=T),
                                  q1   = quantile(x,.05,na.rm=T),
                                  q2   = quantile(x,.95,na.rm=T),
                                  l    = length(x),
                                  sd   = sd(x,na.rm=T))))
}


# =============================================================
# 3  PRE-PROCESSING
# =============================================================

## 3.1  Trend significance at HER level ------------------------
trendFlood  <- calculateTrendSig(
  trendPlot = FloodTrends[FloodTrends$driver == "Clim", ],
  pointagg  = FloodTrends[FloodTrends$driver == "Clim", c(1, 71, 72)])
trendDrought <- calculateTrendSig(
  trendPlot = DroughtTrends[DroughtTrends$driver == "Clim", ],
  pointagg  = DroughtTrends[DroughtTrends$driver == "Clim", c(1, 71, 72)])

## 3.2  Pixel-level trend data per driver ----------------------
clean_pix <- function(tp) tp[!is.na(tp$Var1), ]
TotFPix  <- clean_pix(FloodTrendsP[FloodTrendsP$driver   == driver[5], ])
TotDPix  <- clean_pix(DroughtTrendsP[DroughtTrendsP$driver == driver[5], ])
CliFPix  <- clean_pix(FloodTrendsP[FloodTrendsP$driver   == driver[1], ])
CliDPix  <- clean_pix(DroughtTrendsP[DroughtTrendsP$driver == driver[1], ])
LuFPix   <- clean_pix(FloodTrendsP[FloodTrendsP$driver   == driver[2], ])
LuDPix   <- clean_pix(DroughtTrendsP[DroughtTrendsP$driver == driver[2], ])
ResFPix  <- clean_pix(FloodTrendsP[FloodTrendsP$driver   == driver[3], ])
ResDPix  <- clean_pix(DroughtTrendsP[DroughtTrendsP$driver == driver[3], ])
WuFPix   <- clean_pix(FloodTrendsP[FloodTrendsP$driver   == driver[4], ])
WuDPix   <- clean_pix(DroughtTrendsP[DroughtTrendsP$driver == driver[4], ])

## Base pixel sf object (geometry from Total Flood pixels)
base_pix_sf <- st_transform(
  st_as_sf(TotFPix, coords = c("Var1","Var2"), crs = 4326), crs = 3035)

## 3.3  Pixel bivariate data per driver ------------------------
databipi    <- pix_bidata(TotFPix,  TotDPix,  base_pix_sf)
databipic   <- pix_bidata(CliFPix,  CliDPix,  base_pix_sf)
databipilu  <- pix_bidata(LuFPix,   LuDPix,   base_pix_sf)
databipire  <- pix_bidata(ResFPix,  ResDPix,  base_pix_sf)
databipiwd  <- pix_bidata(WuFPix,   WuDPix,   base_pix_sf)


#river network with reduced drought magnitude
length(which(databipic$y>0))/length(databipic$x[which(!is.na(databipic$y))])

#river network with increased flood magnitude
length(which(databipic$x>0))/length(databipic$x[which(!is.na(databipic$x))])

#Stats for Section 3.4
traj="Decelerating"
length(which(databipic$trcat==traj))/length(databipic$trcat)

bg="Mediterranean"
length(which(databipic$trcat==traj & databipic$Biogeo_id==bg))/length(databipic$trcat[which(databipic$Biogeo_id==bg)])

length(which(is.na(databipic$bi_class)))
# Add trajectory category & biogeo matching
for (obj in c("databipi","databipic","databipilu","databipire","databipiwd")) {
  d <- get(obj)
  d$trcat     <- assign_trcat(d$bi_class)
  d$Biogeo_id <- biogeo_rivers$code[match(d$outl2, biogeo_rivers$outl2)]
  assign(obj, d)
}

## 3.4  HER-level bivariate data per driver --------------------
databitot  <- her_bidata(FloodTrends, DroughtTrends, driver[5], GHshpp, HydroRsf)
databiclim <- her_bidata(FloodTrends, DroughtTrends, driver[1], GHshpp, HydroRsf)
databilu   <- her_bidata(FloodTrends, DroughtTrends, driver[2], GHshpp, HydroRsf)
databire   <- her_bidata(FloodTrends, DroughtTrends, driver[3], GHshpp, HydroRsf)
databiwd   <- her_bidata(FloodTrends, DroughtTrends, driver[4], GHshpp, HydroRsf)

#HER with reduced drought magnitude
length(unique(databiclim$CODEB.x[which(databiclim$y>0)]))

## 3.5  Cumulative driver combinations (HER & pixel) -----------
# HER level
databicr   <- databiclim; databicr$x  <- databiclim$x + databilu$x;                          databicr$y  <- databiclim$y + databilu$y
databicrl  <- databiclim; databicrl$x <- databiclim$x + databilu$x + databire$x;             databicrl$y <- databiclim$y + databilu$y + databire$y
databicrlw <- databiclim; databicrlw$x <- databiclim$x + databilu$x + databire$x + databiwd$x; databicrlw$y <- databiclim$y + databilu$y + databire$y + databiwd$y



for (obj in c("databicr","databicrl","databicrlw")) {
  d <- get(obj); d$bi_class <- make_bi_class(d$x, d$y)
  d$combined_category <- d$bi_class; d$maxicat <- assign_trcat(d$combined_category)
  assign(obj, d)
}


databiclim$maxicat <- assign_trcat(databiclim$bi_class)

# Pixel level
databipicr   <- databipic; databipicr$x  <- databipic$x + databipilu$x;                                  databipicr$y  <- databipic$y + databipilu$y
databipicrl  <- databipic; databipicrl$x <- databipic$x + databipilu$x + databipire$x;                   databipicrl$y <- databipic$y + databipilu$y + databipire$y
databipicrlw <- databipic; databipicrlw$x <- databipic$x + databipilu$x + databipire$x + databipiwd$x;   databipicrlw$y <- databipic$y + databipilu$y + databipire$y + databipiwd$y


#only socioeconomic drivers
databise        <- databire
databise$x      <- databire$x + databilu$x + databiwd$x
databise$y      <- databire$y + databilu$y + databiwd$y
databise$bi_class <- make_bi_class(databise$x, databise$y)
databise$combined_category <- databise$bi_class


traj="Accelerating"
bg="Atlantic"

length(which(databipire$trcat==traj & databipire$Biogeo_id==bg))/length(databipire$trcat[which(databipire$Biogeo_id==bg)])
length(which(databipise$trcat==traj & databipise$Biogeo_id==bg))/length(databipise$trcat[which(databipise$Biogeo_id==bg)])
length(which(databipise$trcat=="Wetting"))

for (obj in c("databipicr","databipicrl","databipicrlw","databipise")) {
  print(obj)
  d <- get(obj); d$bi_class <- make_bi_class(d$x, d$y)
 #d$bi_class[is.na(d$x) | is.na(d$y)] <- NA
  d$combined_category <- d$bi_class; d$trcat <- assign_trcat(d$bi_class)
  assign(obj, d)
}
length(which(is.na(databipic$trcat)))

for (obj in c("databipi","databipic","databipilu","databipire","databipiwd","databipise")) {
  print(obj)
  d <- get(obj); 
  head(d)
  d$maxicat <- assign_trcat(d$bi_class); 
  assign(obj, d)
}


## 3.8  Biogeo map for Extended Figure -------------------------


### match biogeoregions with HRs ----
databitotxHR=databipi
RegioRLi=aggregate(list(val=databitotxHR$y),
                   by = list(HydroR=databitotxHR$HydroRegions_raster_WGS84),
                   FUN = function(x) c(mean=mean(x,na.rm=T),dev=sd(x,na.rm=T),len=length(x),med=median(x,na.rm=T),q1=quantile(x, 0.05, na.rm=T),q3=quantile(x, 0.95, na.rm=T)))
RegioRLi <- do.call(data.frame, RegioRLi)

HRM=na.omit((match(RegioRLi$HydroR,HydroRsf$Id)))
HydroRsf_dom=HydroRsf[HRM,]


layerMatch=match(GHR_riv$outl2,biogeo_rivers$outl2)

GHR_riv$biogeoR=biogeo_rivers$code[layerMatch]
unique(GHR_riv$HydroRegions_raster_WGS84)
verif=(match(GHR_riv$HydroRegions_raster_WGS84,HydroRsf_dom$Id))
unique(verif)
unique(HydroRsf_dom$Id[verif])
GHR_riv$HydrRName=HydroRsf_dom$IRST_NAMEB[verif]


bhp_m=na.omit(unique(match(HydroRsf_dom$Id,GHR_riv$HydroRegions_raster_WGS84)))
HydroRsf_dom$biogeoreg=GHR_riv$biogeoR[bhp_m]

sample=GHR_riv[bhp_m,]

domain_union   <- st_union(st_make_valid(HydroRsf_dom))

biogeo<- st_transform(biogeo, st_crs(HydroRsf_dom))
biogeof_merged <- biogeo %>%
  mutate(name    = if_else(name    %in% c("Pannonian","Steppic"), "Continental", name),
         pre_2012 = if_else(pre_2012 %in% c("PAN","STE"),         "CON",         pre_2012)) %>%
  group_by(name, pre_2012) %>%
  summarize(geometry = st_union(geometry), .groups = "drop") %>%
  st_make_valid()

biogeof_clipped <- st_intersection(biogeof_merged, domain_union)

biogeof_segments <- st_cast(st_cast(biogeof_clipped, "MULTIPOLYGON"), "POLYGON")
biogeof_segments$area <- st_area(biogeof_segments)
label_points <- biogeof_segments %>%
  filter((pre_2012 != "ALP" & pre_2012 != "MED" & area > units::set_units(250000, km^2)) |
           (pre_2012 == "ALP" & area > units::set_units(16000, km^2)) |
           (pre_2012 == "MED" & area > units::set_units(50000, km^2))) %>%
  st_point_on_surface()
biogeof_clipped <- st_transform(biogeof_clipped, st_crs(HydroRsf_dom))

layerMatch=match(GHR_riv$outl2,biogeo_rivers$outl2)

GHR_riv$biogeoR=biogeo_rivers$code[layerMatch]

bhp_m=na.omit(unique(match(HydroRsf_dom$Id,GHR_riv$HydroRegions_raster_WGS84)))
HydroRsf_dom$biogeoreg=GHR_riv$biogeoR[bhp_m]


RegioAg = aggregate(list(oc=GHR_riv$HydrRName),
                    by = list(BR=GHR_riv$biogeoR,HR=GHR_riv$HydroRegions_raster_WGS84),
                    FUN = function(x) c(len=length(x)))
RegioAg <- do.call(data.frame, RegioAg)

ur=unique(RegioAg$HR)
Hydroplot=c()
for (id in 1:length(ur)){
  myr=ur[id]
  maty=which(!is.na((match(RegioAg$HR,myr))))
  lar=RegioAg[maty,]
  if(length(lar$oc)>1){
    sumi=sum(lar$oc)
    lar=lar[which.max(lar$oc),]
    lar$oc=lar$oc/sumi*100
  }else{  lar$oc=100}
  Hydroplot=rbind(Hydroplot,lar)
} 

bioplot <- inner_join(HydroRsf_dom,Hydroplot, by=c("Id"="HR"))
bioplot$BR[which(bioplot$BR=="Pannonian")]="Continental"
bioplot$BR[which(bioplot$BR=="Steppic")]="Continental"

## 3.6  HER-level aggregation of pixel data -------------------
# Match pixel data to HER and biogeo regions
add_hr  <- function(df) inner_join(GHR_riv, df, by = "outl2")
add_bgr <- function(df) inner_join(biogeo_rivers, df, by = "outl2")

clixHR  <- add_hr(databipic);   lxHR  <- add_hr(databipilu)
rexHR   <- add_hr(databipire);  wxHR  <- add_hr(databipiwd); totxHR <- add_hr(databipi)

# Aggregate drought (y) and flood (x) means per HR
agg_both <- function(data_hr) {
  d_agg <- agg_her(data_hr, "y.y")
  f_agg <- agg_her(data_hr, "x.y")
  names(d_agg)[-1] <- paste0("d_", names(d_agg)[-1])
  names(f_agg)[-1] <- paste0("f_", names(f_agg)[-1])
  inner_join(d_agg, f_agg, by = "HR")
}

mbcli <- agg_both(clixHR);  mblu  <- agg_both(lxHR)
mbre  <- agg_both(rexHR);   mbwd  <- agg_both(wxHR); mball <- agg_both(totxHR)

l1 <- nrow(mbcli)
mbfX <- bind_rows(
  mutate(mbcli, names = "Climate"),
  mutate(mbre,  names = "Reservoirs"),
  mutate(mblu,  names = "Landuse"),
  mutate(mbwd,  names = "WaterDemand")
)
mbfH <- mutate(mball, names = "Historical")

colnames(mbfX)[c(2,4,5,8,9,10)] <- c("x","xq1","xq2","y","yq1","yq2")
colnames(mbfH)[c(2,4,5,8,9,10)] <- c("x","xq1","xq2","y","yq1","yq2")

# Attach biogeo and significance to Historical
RegionName <- inner_join(mbcli, HydroRsf_dom, by = c("HR" = "CODEB"))
bioplot    <- inner_join(bioplot, mbfH, by = c("CODEB" = "HR"))
bioplot$BR[bioplot$BR == "Pannonian"] <- "Continental"
bioplot$BR[bioplot$BR == "Steppic"]   <- "Continental"

h2p         <- match(mbfH$HR, bioplot$CODEB)
mbfH$biogeo <- bioplot$BR[h2p]
mbfH$CODEB  <- bioplot$CODEB[h2p]

n2t           <- match(trendFlood$HydroR,  bioplot$CODEB)
trendFlood$Hname  <- bioplot$IRST_NAMEB[n2t]
n2t           <- match(trendDrought$HydroR, bioplot$CODEB)
trendDrought$Hname <- bioplot$IRST_NAMEB[n2t]

mbfH$fsig  <- trendFlood$change[match(mbfH$HR,  trendFlood$HydroR)]
mbfH$dsig  <- trendDrought$change[match(mbfH$HR, trendDrought$HydroR)]
mbfH$fdsig <- abs(mbfH$fsig) + abs(mbfH$dsig)
mbfHsig    <- mbfH[!is.na(mbfH$fdsig) & mbfH$fdsig == 4, ]

## 3.7  Pan-European temporal aggregations --------------------
FloodAgg   <- agg_decade(FloodTrends,   "Total")
DroughtAgg <- agg_decade(DroughtTrends, "Total")

FDcts <- rbind(cbind(FloodAgg$cts,   haz = "Flood"),
               cbind(DroughtAgg$cts, haz = "Drought"))
FDdec <- rbind(cbind(FloodAgg$dec,   hazard = "Flood", year = seq(1955,2015,10)),
               cbind(DroughtAgg$dec, hazard = "Drought", year = seq(1955,2015,10)))


# =============================================================
# 4  MAIN FIGURES
# =============================================================

lalim <- 70

## ── Figure 1 ─────────────────────────────────────────────────
## Temporal change: flood vs drought (ribbon + line)
br_line <- c(seq(-100,-20,20), c(-10,-5,-2,0,2,5,10), seq(20,100,20))

fig1 <- ggplot(FDcts, aes(x = yr, y = value.mean, color = haz, fill = haz)) +
  geom_hline(yintercept = 0, color = "black", linewidth = 1) +
  geom_ribbon(aes(ymin = value.w1.2.5., ymax = value.w2.97.5.), alpha = 0.2, color = NA) +
  geom_line(linewidth = 2) +
  scale_color_manual(values = colorz_haz, name = "Hazard") +
  scale_fill_manual(values  = colorz_haz, name = "Hazard") +
  scale_y_continuous(name = "Change (% 10yRL)", breaks = br_line,
                     trans = scales::modulus_trans(.4)) +
  scale_x_continuous(breaks = seq(1950,2020,5), labels = seq(1950,2020,5),
                     name = "Years", limits = c(1955,2020), expand = c(.001,.001)) +
  theme(axis.title = element_text(size=18,face="bold"),
        title      = element_text(size=22,face="bold"),
        axis.text  = element_text(size=16),
        axis.text.x = element_text(size=16,face="bold"),
        panel.background = element_rect(fill="white",colour="white"),
        panel.grid = element_blank(),
        panel.border = element_rect(linetype="solid",fill=NA,colour="black"),
        legend.title = element_text(size=20,face="bold"),
        legend.text  = element_text(size=16),
        axis.ticks.y = element_blank(),
        panel.grid.major.y = element_line(color="lightgray",linetype="dashed"),
        legend.position = "right",
        legend.key = element_rect(fill="transparent",colour="transparent"),
        legend.key.size = unit(.8,"cm")) +
  ggtitle("Europe – total change in flood & drought return levels")

fig1
ggsave(paste0(plotDir, "/Fig1_temporalCh_fldr_FINAL.jpg"),
       fig1, width=30, height=10, units="cm", dpi=1000)


## ── Figure 2 ─────────────────────────────────────────────────
## Bivariate map: total change at HER (background) + pixel level
databipi$bi_class[databipi$bi_class == "NA-NA" | is.na(databipi$Y2015) |
                    is.na(databipi$y)] <- NA

map_tot <- ggplot(basemap) +
  geom_sf(fill = "white") +
  geom_sf(data = databitot, aes(fill = combined_category, geometry = geometry),
          alpha = 0.7, color = "transparent", size = 0.01, show.legend = FALSE) +
  geom_sf(data = databipi, aes(col = bi_class, geometry = geometry, size = upa),
          alpha = 1, stroke = 0, shape = 15, show.legend = FALSE) +
  geom_sf(fill = NA, color = "gray42") +
  scale_fill_manual(values = colorp) +
  scale_color_manual(values = colorp, na.value = NA) +
  scale_size(range = c(0.08,0.4), trans = "sqrt", guide = "none") +
  coord_sf(xlim = c(min(nco[,1]),max(nco[,1])), ylim = c(min(nco[,2]),max(nco[,2]))) +
  map_theme_bi()

legend_tot <- bi_legend(pal = colorp, dim = 4,
                        xlab = "  +  Drought intensity  -  ",
                        ylab = "  -  Flood intensity  +  ",
                        size = 16, arrows = FALSE)

fig2 <- ggarrange(map_tot, legend_tot, ncol = 2, nrow = 1,
                  widths = c(2,1), heights = c(1,1), vjust = -1)
ggsave(paste0(plotDir, "/Fig2_totchange_bvPIX_HR_FINAL.jpg"),
       fig2, width=20, height=20, units="cm", dpi=500)


## ── Figure 3 ─────────────────────────────────────────────────
## Bivariate map: climate driver only
map_cli <- ggplot(basemap) +
  geom_sf(fill = "white") +
  geom_sf(data = databiclim, aes(fill = combined_category, geometry = geometry),
          alpha = 0.7, color = "transparent", show.legend = FALSE) +
  geom_sf(data = databipic, aes(col = bi_class, geometry = geometry, size = upa),
          alpha = 1, stroke = 0, shape = 15, show.legend = FALSE) +
  geom_sf(fill = NA, color = "gray42") +
  scale_fill_manual(values = colorp) +
  scale_color_manual(values = colorp, na.value = NA) +
  scale_size(range = c(0.08,0.4), trans = "sqrt", guide = "none") +
  coord_sf(xlim = c(min(nco[,1]),max(nco[,1])), ylim = c(min(nco[,2]),max(nco[,2]))) +
  map_theme_bi()

fig3 <- ggarrange(map_cli, legend_tot, ncol = 2, nrow = 1,
                  widths = c(2,1), heights = c(1,1), vjust = -1)
ggsave(paste0(plotDir, "/Fig3_climchange_bvPIX_HR_FINAL.jpg"),
       fig3, width=20, height=20, units="cm", dpi=500)




## ── Figure 5 ─────────────────────────────────────────────────
## Bivariate scatter: driver contribution per HER (coloured by driver)

fig5 <- ggplot() +
  geom_vline(xintercept = 0, color = "black", linewidth = 1.5) +
  geom_hline(yintercept = 0, color = "black", linewidth = 1.5) +
  geom_point(data = mbfX,
             aes(x = x, y = y, fill = names, size = d_val.l), shape = 21,
             alpha = .5, stroke = 0) +
  geom_point(data = mbfX,
             aes(x = x, y = y, color = names, size = d_val.l), shape = 1,
             stroke = .8, alpha = .6, show.legend = FALSE) +
  scale_fill_manual(values = colorn, name = "Drivers",
                    labels = c("Climate","Land use","Reservoirs","Water demand")) +
  scale_color_manual(values = colorn, guide = "none") +
  scale_size(range = c(2,10), guide = "none") +
  coord_cartesian(xlim = c(-lalim,lalim), ylim = c(-lalim,lalim)) +
  quadrant_annots +
  labs(x = "Change in drought intensity (%)", y = "Change in flood intensity (%)") +
  scale_x_continuous(trans = scales::modulus_trans(.5),
                     breaks = c(-200,-100,-50,-20,-10,-5,0,5,10,20,50,100,200),
                     labels=rev(c(-200,-100,-50,-20,-10,-5,0,5,10,20,50,100,200)),
                     minor_breaks = seq(-100,100,10)) +
  scale_y_continuous(trans = scales::modulus_trans(.5),
                     breaks = c(-200,-100,-50,-20,-10,-5,0,5,10 ,20,50,100,200),
                     minor_breaks = seq(-100,100,10), expand = c(0,0)) +
  guides(fill = guide_legend(override.aes = list(size = 6))) +
  bv_scatter_theme()
fig5

ggsave(paste0(plotDir, "/Fig5_Contribution_HR_FINAL.jpg"),
       fig5, width=20, height=20, units="cm", dpi=500)


## ── Figure 4 ─────────────────────────────────────────────────
## Stacked barplot of trajectory by biogeoregion (total change, pixel level)
databipi2      <- databipi
databipi2$trcat <- assign_trcat(databipi2$bi_class)
fig4 <- biogeo_barplot(databipi2, colorp)
fig4
ggsave(paste0(plotDir, "/Fig4_Regional_trajectories_tot_FINAL.jpg"),
       fig4, width=20, height=20, units="cm", dpi=800)

## ── Figure 6 ─────────────────────────────────────────────────
## Climate stacked barplot by biogeoregion
databipic2      <- databipic
length(which(is.na(databipic$bi_class)))/length(databipic$bi_class)
databipic2$trcat <- assign_trcat(databipic2$bi_class)
fig6 <- biogeo_barplot(databipic2, colorp)
fig6
ggsave(paste0(plotDir, "/Fig6_Regional_trajectories_clim_FINAL.jpg"),
       fig6, width=20, height=20, units="cm", dpi=800)

## ── S8 · Trajectory barplot – socioeconomic combined signal
databipi2se       <- databipise
databipi2se$trcat <- assign_trcat(databipi2se$bi_class)
length(which(databipi2se$trcat=="Wetting"))
figS8 <- biogeo_barplot(databipi2se, colorp)
figS8
ggsave(paste0(plotDir, "/Fig5_Regional_trajectories_SE_FINAL.jpg"),
       figS8, width=20, height=20, units="cm", dpi=800)

## ── Figure 7 ─────────────────────────────────────────────────
## Sankey at HER level (cumulative driver addition)
databiclim$maxicat <- assign_trcat(databiclim$bi_class)
databicr$maxicat   <- assign_trcat(databicr$bi_class)
databicrl$maxicat  <- assign_trcat(databicrl$bi_class)
databicrlw$maxicat <- assign_trcat(databicrlw$bi_class)

fig7 <- make_sankey(
  list(databiclim$maxicat, databicr$maxicat, databicrl$maxicat, databicrlw$maxicat),
  c("OnlyClimate","ClimateLanduse","ClimateLanduseReservoir","AllDrivers"),
  loscolors, ltot = nrow(databiclim)
) + ggtitle("Change in Europe's hydrological cycle (1950s–2010s) – HER level")

ggsave(paste0(plotDir, "/Fig7_Sankey_HER_FINAL.jpg"),
       fig7, width=30, height=20, units="cm", dpi=300)



# =============================================================
#  Barplot – % of pixels per bivariate trajectory, ranked by driver
# =============================================================

# ── 1. Total pixel count (reference denominator) ─────────────
n_total <- length((databipic$bi_class))

# ── 2. Count pixels per driver & trajectory ──────────────────
count_traj <- function(data, driver_name) {
  data %>%
    filter(!is.na(bi_class), !is.na(trcat)) %>%
    group_by(bi_class) %>%
    summarise(n = n(), .groups = "drop") %>%
    mutate(
      driver = driver_name,
      pct    = n / n_total * 100
    )
}


traj_counts <- bind_rows(
  count_traj(databipise,   "All socioeconomics"),
  count_traj(databipic,  "Climate"),
  count_traj(databipilu, "Land use"),
  count_traj(databipire, "Reservoirs"),
  count_traj(databipiwd, "Water demand")
)


# ── 3. Rank drivers by total % of classified pixels ──────────
driver_order <- traj_counts %>%
  group_by(driver) %>%
  summarise(total_pct = sum(pct)) %>%
  arrange(total_pct) %>%          # ascending → coord_flip puts highest at top
  pull(driver)

# ── 4. Fix stack order: meaningful trajectory sequence ────────
traj_order <- c("Drying", "Decelerating", "Stable", "Wetting", "Accelerating")

# traj_counts <- traj_counts %>%
#   mutate(
#     driver = factor(driver, levels = driver_order),
#     trcat  = factor(trcat,  levels = traj_order)
#   )
head(traj_counts)
traj_counts$rm=0
traj_counts$rm[traj_counts$bi_class %in% c("2-2","3-2","2-3","3-3")] <- 1
traj_counts=traj_counts[-which(traj_counts$rm==1),]
# ── 5. Label positions (midpoint of each stacked segment) ─────
traj_counts <- traj_counts %>%
  group_by(driver) %>%
  arrange(driver, bi_class) %>%
  mutate(
    cum_pct   = cumsum(pct),
    label_pos = cum_pct - pct / 2
  ) %>%
  ungroup()

# ── 3. Rank drivers by total % of classified pixels ──────────
driver_order <- traj_counts %>%
  group_by(driver) %>%
  summarise(total_pct = sum(pct)) %>%
  arrange(total_pct) %>%          # ascending → coord_flip puts highest at top
  pull(driver)

driver_order

# Total % at end of each bar
bar_totals <- traj_counts %>%
  group_by(driver) %>%
  summarise(total_pct = sum(pct))

# Ordering by colord
traj_counts$ordf <- colord[match(traj_counts$bi_class, names(colord))]
traj_counts      <- traj_counts %>% mutate(bi_class = reorder(bi_class, ordf, FUN = mean))

# ── 6. Plot ───────────────────────────────────────────────────
fig_bv_bar <- ggplot(traj_counts,
                     aes(x = driver, y = pct, fill = bi_class)) +
  geom_col(width = 0.65, color = "transparent", linewidth = 0.3) +
  # Segment % labels — only show if segment wide enough to read
  # geom_text(data = filter(traj_counts, pct > 2),
  #           aes(y = label_pos, label = paste0(round(pct, 1), "%")),
  #           size = 3.5, fontface = "bold", color = "white") +
  # Total % at bar end
  geom_text(data = bar_totals,
            aes(x = driver, y = total_pct,
                label = paste0(round(total_pct, 0), "%")),
            inherit.aes = FALSE,
            hjust = -0.15, size = 4.5, fontface = "bold", color = "grey20") +
  
  scale_fill_manual(values = colorp) +
  scale_y_continuous(
    name   = "Share of river network (%)",
    limits = c(0, max(bar_totals$total_pct) * 1.2),
    expand = c(0, 0)
  ) +
  scale_x_discrete(limits = driver_order, name = NULL) +
  geom_hline(yintercept =0,
             color = "black", linewidth = 1) +
  # scale_x_discrete(name = NULL) +
  coord_flip() +
  theme(
    axis.title.x       = element_text(size = 14, face = "bold"),
    axis.text          = element_text(size = 13),
    axis.text.y        = element_text(face = "bold"),
    axis.text.x        = element_blank(),
    axis.ticks         = element_blank(),
    panel.background   = element_rect(fill = "white"),
    # panel.grid.major.x = element_line(color = "grey85", linetype = "dashed"),
    
    panel.grid.major = element_blank(),
    panel.border       = element_blank(), 
    legend.position    = "none",
    legend.text        = element_text(size = 11),
    legend.title       = element_text(size = 12, face = "bold"),
    legend.key.size    = unit(0.7, "cm")
  ) 
  #ggtitle("Share of river pixels by joint flood–drought trajectory (1955–2015)")

fig_bv_bar

ggsave(paste0(plotDir, "/FigX_bv_traj_barplotv2_FINAL.jpg"),
       fig_bv_bar, width = 16, height = 22, units = "cm", dpi = 1000)

## ── Figure Supp ─────────────────────────────────────────────────
## Number of pixels where trajectory could not be computed

biogeo_barplot_na <- function(databipi2, tsize = 22, osize = 28) {
  
  # Recode Pannonian -> Continental
  databipi2$Biogeo_id[databipi2$Biogeo_id == "Pannonian"] <- "Continental"
  
  # Exclude unwanted regions
  databipi2 <- databipi2[!databipi2$Biogeo_id %in% c("Steppic", "BlackSea", "Arctic"), ]
  databipi2 <- databipi2[!databipi2$Biogeo_id %in% c(NA), ]
  
  # Per-region: total pixels and NA pixels
  region_stats <- databipi2 %>%
    group_by(reg = Biogeo_id) %>%
    summarise(
      total = n(),
      na_count = sum(!is.na(bi_class)),
      .groups = "drop"
    ) %>%
    mutate(pct_na = na_count / total)
  
  # Whole-domain horizontal line value
  domain_pct_na <- sum(!is.na(databipi2$bi_class)) / nrow(databipi2)
  
  ggplot(region_stats, aes(x = reg, y = pct_na)) +
    geom_col(width = 0.7, fill = "black", color = "transparent") +
    geom_hline(yintercept = domain_pct_na,
               color = "red", linewidth = 1, linetype = "dashed") +
    annotate("text",
             x = Inf, y = domain_pct_na,
             label = paste0("Domain: ", round(domain_pct_na * 100, 1), "%"),
             hjust = 1.05, vjust = -0.5,
             size = 4.5, fontface = "bold", color = "red") +
    geom_text(aes(y = pct_na + 0.01,
                  label = paste0(round(pct_na * 100, 1), "%")),
              hjust = 0.5, vjust = 0,
              size = 4.5, fontface = "bold", color = "grey20") +
    scale_x_discrete(name = "Biogeographic regions",
                     labels = reg_labels,
                     expand = c(0, 0)) +
    scale_y_continuous(
      name   = "Analysed river network (%)",
      labels = function(x) paste0(x * 100),
      expand = c(0, 0),
      limits = c(0, max(region_stats$pct_na) * 1.25)
    ) +
    geom_vline(xintercept = 0.5, color = "black", linewidth = 1) +
    geom_hline(yintercept = 0,   color = "black", linewidth = 1) +
    theme(
      axis.title       = element_text(size = 16, face = "bold"),
      axis.text        = element_text(size = 15),
      axis.text.y      = element_text(face = "bold"),
      axis.ticks       = element_blank(),
      panel.background = element_rect(fill = "white"),
      panel.grid.major = element_blank(),
      panel.border     = element_blank(),
      legend.position  = "none"
    )
}

figNA<-biogeo_barplot_na(databipic)
figNA
ggsave(paste0(plotDir, "/FigS_nona_FINAL.jpg"),
       figNA, width = 16, height = 16, units = "cm", dpi = 1000)
## ── Figure 8 ─────────────────────────────────────────────────
## Extended Figure: biogeoregion + HER map
fig8 <- ggplot(basemap) +
  geom_sf(fill = "gray95") +
  geom_sf(fill = NA, color = "grey") +
  geom_sf(data = bioplot, aes(fill = factor(BR), geometry = geometry),
          color = "black", alpha = .7, size = 0.0001) +
  coord_sf(xlim = c(min(nco[,1]),max(nco[,1])), ylim = c(min(nco[,2]),max(nco[,2]))) +
  scale_fill_manual(values = colx_bgr, name = "BioGeoRegions") +
  labs(x = "Longitude", y = "Latitude") +
  map_theme_bi(tsize, osize) +
  theme(legend.position = "right")

ggsave(paste0(plotDir, "/Fig8_Map_HRxBG_FINAL.jpg"),
       fig8, width=25, height=20, units="cm", dpi=300)


## ── Figure 9 ─────────────────────────────────────────────────
## Historical HER scatter coloured by biogeoregion
fig9 <- ggplot() +
  geom_vline(xintercept = 0, color = "black", linewidth = 1.5) +
  geom_hline(yintercept = 0, color = "black", linewidth = 1.5) +
  geom_point(data = mbfH, aes(x = d_val.mean, y = f_val.mean, color = biogeo, size = d_val.l),
             stroke = 0, alpha = .7) +
  geom_text(data = mbfH, aes(x = d_val.mean, y = f_val.mean, label = CODEB),
            size = 3, color = "black", fontface = "bold") +
  geom_point(data = mbfHsig, aes(x = d_val.mean, y = f_val.mean, size = d_val.l),
             fill = NA, shape = 1, stroke = 1, col = "black", alpha = .8,
             show.legend = FALSE) +
  scale_color_manual(values = colx_bgr, name = "Regions") +
  scale_size(range = c(6,14), guide = "none") +
  coord_cartesian(xlim = c(-lalim, lalim+10), ylim = c(-lalim, lalim)) +
  quadrant_annots +
  labs(x = "Change in drought intensity (%)", y = "Change in flood intensity (%)") +
  scale_x_continuous(trans = scales::modulus_trans(.5),
                     breaks = c(-200,-100,-50,-20,-10,-5,0,5,10,20,50,100,200),
                     minor_breaks = seq(-100,100,10)) +
  scale_y_continuous(trans = scales::modulus_trans(.5),
                     breaks = c(-200,-100,-50,-10,-5,0,5,10,50,100,200),
                     minor_breaks = seq(-100,100,10), expand = c(0,0)) +
  guides(color = guide_legend(override.aes = list(size = 6))) +
  bv_scatter_theme() + theme(legend.position = "right")
fig9
ggsave(paste0(plotDir, "/Fig9_Contribution_HRxBG_FINAL.jpg"),
       fig9, width=25, height=20, units="cm", dpi=300)


# =============================================================
# 5  SUPPLEMENTARY FIGURES  (only if RUN_SUPP == TRUE)
# =============================================================

if (RUN_SUPP) {

  ## ── S1 · Decade boxplot – flood vs drought ────────────────
  br_bxp <- c(seq(-100,-10,10), seq(-5,5,5), seq(10,100,10))
  xlabs  <- seq(1950,2010,10)

  figS1 <- ggplot(FDdec) +
    geom_linerange(aes(x=year, ymin=w1, ymax=w2, color=hazard, group=hazard),
                   position=position_dodge2(width=9), lwd=1, alpha=.8) +
    scale_color_manual(values=colorz_haz, name="Hazard") +
    new_scale_color() +
    geom_rect(aes(xmin=year-4.5, xmax=year+4.5, ymin=cq1, ymax=cq2,
                  fill=hazard, group=hazard),
              alpha=.5, position=position_dodge(width=9)) +
    geom_point(aes(x=year, y=changeC, color=hazard, group=hazard),
               position=position_dodge(width=9), size=3) +
    geom_rect(aes(xmin=year-4.5, xmax=year+4.5, ymin=med-1e-1, ymax=med+1e-1,
                  fill=hazard, group=hazard),
              alpha=1, position=position_dodge(width=9)) +
    scale_y_continuous(name="Mean change (% 10yRL)", breaks=br_bxp,
                       trans=scales::modulus_trans(.6)) +
    scale_x_continuous(breaks=xlabs, labels=xlabs, name="Decades",
                       minor_breaks=seq(1955,2005,10), expand=c(.01,.01)) +
    scale_fill_manual(values=colorz_haz, name="Hazard") +
    scale_color_manual(values=c("Flood"="darkblue","Drought"="darkorange"),
                       name="Hazard") +
    theme(axis.title=element_text(size=18,face="bold"),
          title=element_text(size=22,face="bold"),
          axis.text=element_text(size=16),
          axis.text.x=element_text(size=16,face="bold"),
          panel.background=element_rect(fill="white",colour="white"),
          panel.grid=element_blank(),
          panel.border=element_rect(linetype="solid",fill=NA,colour="black"),
          legend.title=element_text(size=20,face="bold"),
          legend.text=element_text(size=16),
          axis.ticks.y=element_blank(),
          panel.grid.major.y=element_line(color="lightgray",linetype="dashed"),
          panel.grid.minor.x=element_line(colour="grey23",linetype="dashed"),
          legend.position="right",
          legend.key=element_rect(fill="transparent",colour="transparent"),
          legend.key.size=unit(.8,"cm")) +
    ggtitle("Europe")

  ggsave(paste0(plotDir, "/FigS1_bxplot_fldr_FINAL.jpg"),
         figS1, width=30, height=20, units="cm", dpi=1000)


  ## ── S2 · Bivariate aggregate scatter per driver (EU mean) ─
  mbf <- data.frame(
    x    = c(mean(databiclim$y,na.rm=T), mean(databire$y,na.rm=T),
             mean(databilu$y,na.rm=T),   mean(databiwd$y,na.rm=T)),
    y    = c(mean(databiclim$x,na.rm=T), mean(databire$x,na.rm=T),
             mean(databilu$x,na.rm=T),   mean(databiwd$x,na.rm=T)),
    xq1  = c(quantile(databiclim$y,.05,na.rm=T), quantile(databire$y,.05,na.rm=T),
             quantile(databilu$y,.05,na.rm=T),   quantile(databiwd$y,.05,na.rm=T)),
    xq2  = c(quantile(databiclim$y,.95,na.rm=T), quantile(databire$y,.95,na.rm=T),
             quantile(databilu$y,.95,na.rm=T),   quantile(databiwd$y,.95,na.rm=T)),
    yq1  = c(quantile(databiclim$x,.05,na.rm=T), quantile(databire$x,.05,na.rm=T),
             quantile(databilu$x,.05,na.rm=T),   quantile(databiwd$x,.05,na.rm=T)),
    yq2  = c(quantile(databiclim$x,.95,na.rm=T), quantile(databire$x,.95,na.rm=T),
             quantile(databilu$x,.95,na.rm=T),   quantile(databiwd$x,.95,na.rm=T)),
    names = c("Climate","Reservoirs","Landuse","WaterDemand")
  )

  figS2 <- ggplot(mbf, aes(x=x, y=y, fill=names)) +
    geom_linerange(aes(ymin=yq1, ymax=yq2, color=names), lwd=1, alpha=.8) +
    geom_linerange(aes(xmin=xq1, xmax=xq2, color=names), lwd=1, alpha=.8) +
    geom_point(shape=21, color="black", size=5) +
    scale_fill_manual(values=colorn, name="Drivers",
                      labels=c("Climate","Land use","Reservoirs","Water demand")) +
    scale_color_manual(values=colorn, name="Drivers",
                       labels=c("Climate","Land use","Reservoirs","Water demand")) +
    coord_cartesian(xlim=c(-100,100), ylim=c(-50,50)) +
    annotate("text",x= 50,y= 20,label="Wetting",     color="#169dd0",size=7,fontface="bold") +
    annotate("text",x=-50,y= 20,label="Accelerating",color="#174f28",size=7,fontface="bold") +
    annotate("text",x=-50,y=-20,label="Drying",       color="#dd6a29",size=7,fontface="bold") +
    annotate("text",x= 50,y=-20,label="Decelerating", color="burlywood",size=7,fontface="bold") +
    labs(x="Change in drought (%)", y="Change in flood (%)") +
    scale_x_continuous(trans=scales::modulus_trans(.5),
                       breaks=seq(-100,100,50), minor_breaks=seq(-100,100,10)) +
    scale_y_continuous(trans=scales::modulus_trans(.5),
                       breaks=seq(-100,100,25), minor_breaks=seq(-100,100,10)) +
    guides(fill=guide_legend(override.aes=list(size=6))) +
    bv_scatter_theme()

  ggsave(paste0(plotDir, "/FigS2_Contribution_bv_FINAL.jpg"),
         figS2, width=25, height=20, units="cm", dpi=300)


  ## ── S3 · Sankey at pixel level ───────────────────────────
  databipic$maxicat    <- assign_trcat(databipic$bi_class)
  databipicr$maxicat   <- assign_trcat(databipicr$bi_class)
  databipicrl$maxicat  <- assign_trcat(databipicrl$bi_class)
  databipicrlw$maxicat <- assign_trcat(databipicrlw$bi_class)

  figS3 <- make_sankey(
    class_list=list(databipic$maxicat, databipicr$maxicat,
         databipicrl$maxicat, databipicrlw$maxicat),
    c("OnlyClimate","ClimateLanduse","ClimateLanduseReservoir","AllDrivers"),
    loscolors, ltot = nrow(databipic)
  ) + ggtitle("Change in Europe's hydrological cycle (1950s–2010s) – network level")

  ggsave(paste0(plotDir, "/FigS3_Sankey_pixels_FINAL.jpg"),
         figS3, width=30, height=20, units="cm", dpi=500)


  ## ── S4 · Trajectory proportion barplot by driver ─────────
  s1    <- length(unique(mbfH$HR))
  mbfX$traj <- case_when(
    mbfX$y >= 0 & mbfX$x >  0 ~ "Wetting",
    mbfX$y <= 0 & mbfX$x >= 0 ~ "Decelerating",
    mbfX$y <= 0 & mbfX$x <= 0 ~ "Drying",
    mbfX$y >= 0 & mbfX$x <= 0 ~ "Accelerating",
    TRUE ~ "none"
  )
  magg2 <- do.call(data.frame,
    aggregate(list(val = mbfX$HR),
              by  = list(reg = mbfX$names, traj = mbfX$traj),
              FUN = length))
  magg2$perc <- magg2$val / s1 * 100
  atraj      <- c("Wetting","Decelerating","Drying","Accelerating")
  clabels_d  <- c("Climate","Landuse","Reservoirs","WaterDemand")

  for (driv in atraj) {
    magg3 <- magg2[magg2$traj == driv, ]
    miss  <- setdiff(clabels_d, magg3$reg)
    if (length(miss) > 0)
      magg3 <- rbind(magg3, data.frame(reg=miss, traj=driv, val=0, perc=0))
    magg3$perc <- as.numeric(magg3$val) / s1 * 100

    pS4 <- ggplot(magg3, aes(x=reg, y=perc, fill=reg)) +
      geom_bar(stat="identity") +
      scale_fill_manual(values=colorn) +
      scale_x_discrete(name=" ") +
      scale_y_continuous("HydroRegion (%)", limits=c(0,100)) +
      theme(axis.title=element_text(size=16), axis.text=element_text(size=12),
            panel.background=element_rect(fill="white",colour="grey1"),
            panel.border=element_rect(linetype="solid",fill=NA,colour="black"),
            legend.position="none", plot.title=element_text(size=20),
            legend.key=element_rect(fill="transparent",colour="transparent"),
            legend.key.size=unit(.8,"cm")) +
      ggtitle(driv)

    ggsave(paste0(plotDir, "/FigS4_driverxtraj_", driv, ".jpg"),
           pS4, width=12, height=8, units="cm", dpi=300)
  }


  ## ── 5 · biogeo barplot – climate only ─────────
  figS5 <- biogeo_barplot(databipic2, colorp)
  ggsave(paste0(plotDir, "/Fig5_Regional_trajectories_cli_FINAL.jpg"),
         figS5, width=20, height=20, units="cm", dpi=800)


  ## ── S6 · HER map coloured by HR identity ─────────────────
  databic2 <- inner_join(mbfH, GHshpp, by = c("HR" = "IRST_NAMEB"))
  coHR     <- sample(grDevices::colors()[grep("gr(a|e)y|white|black",
                       grDevices::colors(), invert=TRUE)],
                     length(unique(databiclim$HR)))
  figS6 <- ggplot(basemap) +
    geom_sf(fill="white") +
    geom_sf(data=databic2, aes(fill=factor(HR),geometry=geometry),
            alpha=.7, color="black", size=.001, show.legend=FALSE) +
    geom_sf(fill=NA, color="gray42") +
    scale_fill_manual(values=coHR) +
    coord_sf(xlim=c(min(nco[,1]),max(nco[,1])), ylim=c(min(nco[,2]),max(nco[,2]))) +
    map_theme_bi()

  ggsave(paste0(plotDir, "/FigS6_mapreg_HR.jpg"),
         figS6, width=30, height=20, units="cm", dpi=300)


  ## ── S7 · Bivariate map: socio-economic combined signal ────
  databipi$bi_class[databipi$bi_class == "NA-NA" | is.na(databipi$Y2015) |
                      is.na(databipi$y)] <- NA
  databipise        <- databipire
  databipise$x      <- databipire$x + databipilu$x + databipiwd$x
  databipise$y      <- databipire$y + databipilu$y + databipiwd$y
  databipise$bi_class <- make_bi_class(databipise$x, databipise$y)
  databipise$bi_class[is.na(databipise$x) | is.na(databipise$y)] <- NA

  databise        <- databire
  databise$x      <- databire$x + databilu$x + databiwd$x
  databise$y      <- databire$y + databilu$y + databiwd$y
  databise$bi_class <- make_bi_class(databise$x, databise$y)
  databise$combined_category <- databise$bi_class

  map_se <- ggplot(basemap) +
    geom_sf(fill="white") +
    geom_sf(data=databise, aes(fill=combined_category,geometry=geometry),
            alpha=.7, color="transparent", show.legend=FALSE) +
    geom_sf(data=databipise, aes(col=bi_class,geometry=geometry,size=upa),
            alpha=1, stroke=0, shape=15, show.legend=FALSE) +
    geom_sf(fill=NA, color="gray42") +
    scale_fill_manual(values=colorp) +
    scale_color_manual(values=colorp, na.value=NA) +
    scale_size(range=c(0.08,0.4), trans="sqrt", guide="none") +
    coord_sf(xlim=c(min(nco[,1]),max(nco[,1])), ylim=c(min(nco[,2]),max(nco[,2]))) +
    map_theme_bi()

  figS7 <- ggarrange(map_se, legend_tot, ncol=2, nrow=1,
                     widths=c(2,1), heights=c(1,1), vjust=-1)
  ggsave(paste0(plotDir, "/FigS7_sechange_bvPIX_HR_FINAL.jpg"),
         figS7, width=20, height=20, units="cm", dpi=500)



} # end if (RUN_SUPP)


# =============================================================
# 6  SAVE OUTPUTS
# =============================================================

write.csv(mbfH, file = paste0(hydroDir, "/Trajectories/Histo_res_bvHR_F.csv"),  row.names = FALSE)
write.csv(mbfX, file = paste0(hydroDir, "/Trajectories/Drivers_res_bvHR_F.csv"), row.names = FALSE)
save(databipise, file = paste0(hydroDir, "/Trajectories/SEchanges_bivariate_F.Rdata"))
save(databipic, file = paste0(hydroDir, "/Trajectories/CLchanges_bivariate_F.Rdata"))
