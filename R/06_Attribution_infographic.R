## ============================================================
##  06_Attribution_infographic.R
##  Simplified public-facing figure derived from
##  05_Attribution_bivariate.R
##
##  "What is changing Europe's floods and droughts?"
##   - Left  : map of climate-only effect per hydro region
##             (Wetter / Drier / More extreme / Calmer)
##   - Right : summary table of flood & drought arrows for
##             climate categories and human activities
## ============================================================

get_script_dir <- function() {
  a <- commandArgs(FALSE)
  fa <- grep("^--file=", a, value = TRUE)
  if (length(fa) > 0) {
    return(dirname(normalizePath(sub("^--file=", "", fa[1]))))
  }
  if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
    return(dirname(rstudioapi::getSourceEditorContext()$path))
  }
  getwd()
}
setwd(get_script_dir())
source("functions_trends.R")
source("config_paths.R") # hydroDir, plotDir, geoDir, ID_MULT, ...

library(patchwork)

# =============================================================
# 0  PALETTE & LABELS (matches target infographic)
# =============================================================

# Climate-effect categories -> colours (blue / orange / purple / green)
cat_colors <- c(
  "Wetting" = "#2f78c0", # NW & Central Europe
  "Drying" = "#e8871e", # Iberia, Finland, Baltics, Bulgaria, Greece
  "Accelerating" = "#8e44ad", # Norway, Italy, Scotland
  "Decelerating" = "#2e9e6f", # Sweden, the Alps, Romania
  "Stable" = "grey"
)



# Category description shown in the right-hand table
cat_desc <- c(
  "Wetting" = "North-west and Central Europe",
  "Drying" = "Iberia, the Baltics",
  "Accelerating" = "Norway, Italy, Scotland",
  "Decelerating" = "Sweden, Poland, Romania",
  "Stable" = "Cevennes"
)

# Arrow direction per category: more severe (up), less severe (down)
#   value: "up" = more severe, "down" = less severe, "flat" = little effect
cat_arrows <- data.frame(
  category = c("Wetting", "Drying", "Accelerating", "Decelerating", "Stable"),
  flood = c("up", "down", "up", "down", "flat"),
  drought = c("down", "up", "up", "down", "flat"),
  stringsAsFactors = FALSE
)

# Human-activity arrows (qualitative summary, same direction across Europe)
human_arrows <- data.frame(
  activity = c("Land use", "Dams and reservoirs", "Water use"),
  desc = c(
    "Changes in land cover, such as growing cities",
    "Store water and even out river flows",
    "More water taken from rivers"
  ),
  flood = c("up", "down", "flat"),
  drought = c("up", "down", "up"),
  stringsAsFactors = FALSE
)


# =============================================================
# 1  DATA LOADING  (same sources as the main script)
# =============================================================

## 1.1  Outlets ------------------------------------------------
if (!exists("outf")) {
  rspace <- read.csv(file.path(hydroDir, "subspace_efas.csv"))[, -1]
  outf <- c()
  for (Nsq in 1:88) {
    nrspace <- rspace[Nsq, ]
    outhybas <- outletopen(hydroDir, "GeoData/efas_rnet_100km_01min", nrspace)
    Idstart <- as.numeric(Nsq) * 10000
    Idstart2 <- as.numeric(Nsq) * ID_MULT
    if (length(outhybas$outlets) > 0) {
      outhybas$outlets <- seq(Idstart + 1, Idstart + length(outhybas$outlets))
      outhybas$outl2 <- seq(Idstart2 + 1, Idstart2 + length(outhybas$outlets))
      outhybas$latlong <- paste(round(outhybas$Var1, 4), round(outhybas$Var2, 4), sep = " ")
      outf <- rbind(outf, outhybas)
    }
  }
}

## 1.2  Flood & drought outputs --------------------------------
load(file = file.path(floodDir, "outputs_flood_year_relxHR_FINAL.Rdata"))
load(file = file.path(droughtDir, "outputs_drought_nonfrost_relxHR_FINAL.Rdata"))

FloodTrends <- Output_fl_year$TrendRegio
DroughtTrends <- Output_dr_nonfrost$TrendRegio
FloodTrendsP <- Output_fl_year$TrendPix
DroughtTrendsP <- Output_dr_nonfrost$TrendPix
driver <- unique(DroughtTrends$driver) # [1]Clim [2]Lu [3]Res [4]Wu [5]Total

## 1.3  Spatial layers -----------------------------------------
GHshpp <- read_sf(dsn = file.path(geoDir, "HER/her_all_adjusted.shp"))
HydroRsf <- fortify(GHshpp)

# Pixel -> hydroregion mapping (river pixels tagged with their HER id)
GridHR <- raster(file.path(geoDir, "HER/HydroRegions_raster_WGS84.tif"))
GHR <- as.data.frame(GridHR, xy = TRUE)
GHR <- GHR[!is.na(GHR[, 3]), ]
GHR$llcoord <- paste(round(GHR$x, 4), round(GHR$y, 4), sep = " ")
GHR_riv <- inner_join(GHR, outf, by = c("llcoord" = "latlong"))

outll <- outletopen(hydroDir, "GeoData/efas_rnet_100km_01min")
cord.dec <- SpatialPoints(outll[, c(2, 3)], proj4string = CRS("+proj=longlat"))
cord.UTM <- spTransform(cord.dec, CRS("+init=epsg:3035"))
nco <- cord.UTM@coords
world <- ne_countries(scale = "medium", returnclass = "sf")
basemap <- st_transform(world, crs = 3035)

# Biogeographic regions (for the clipped border overlay)
biogeo <- read_sf(dsn = file.path(geoDir, "eea_3035_biogeo-regions_2016/BiogeoRegions2016_wag84.shp"))


# =============================================================
# 2  HELPERS (bivariate classification -> trajectory)
# =============================================================
breaker1 <- 0
breaker2 <- 5

make_bi_class <- function(x, y, b1 = breaker1, b2 = breaker2) {
  classify <- function(v) {
    cl <- rep(NA_integer_, length(v))
    cl[v <= -b2] <- 1L
    cl[v > -b2 & v < b1] <- 2L
    cl[v >= b1 & v < b2] <- 3L
    cl[v >= b2] <- 4L
    cl
  }
  paste(classify(y), classify(x), sep = "-")
}

assign_trcat <- function(bi) {
  tr <- rep(NA, length(bi))
  tr[bi %in% c("2-2", "3-2", "2-3", "3-3")] <- "Stable"
  tr[bi %in% c("1-1", "1-2", "2-1")] <- "Drying"
  tr[bi %in% c("1-4", "2-4", "1-3")] <- "Accelerating"
  tr[bi %in% c("4-1", "4-2", "3-1")] <- "Decelerating"
  tr[bi %in% c("4-4", "4-3", "3-4")] <- "Wetting"
  tr[is.na(bi) | bi == "NA-NA"] <- NA
  tr
}

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
  Flpl$bi_class <- make_bi_class(Flpl$x, Flpl$y)
  Flpl
}

clean_pix <- function(tp) tp[!is.na(tp$Var1), ]

## Pixel-level bivariate data for one driver
pix_bidata <- function(FloodPix, DroughtPix, base_sf) {
  db <- base_sf
  db$x <- FloodPix$Y2015
  db$y <- DroughtPix$Y2015
  db$bi_class <- make_bi_class(db$x, db$y)
  db$bi_class[is.na(db$x) | is.na(db$y)] <- NA
  db
}


# =============================================================
# 3  CLIMATE-ONLY EFFECT  (pixel-based share per hydroregion)
# =============================================================

## 3.1  Pixel-level climate bivariate classification -----------
CliFPix <- clean_pix(FloodTrendsP[FloodTrendsP$driver == driver[1], ])
CliDPix <- clean_pix(DroughtTrendsP[DroughtTrendsP$driver == driver[1], ])

# Base pixel sf (geometry from the flood pixels, matches script 05)
base_pix_sf <- st_transform(
  st_as_sf(CliFPix, coords = c("Var1", "Var2"), crs = 4326),
  crs = 3035
)

databipic <- pix_bidata(CliFPix, CliDPix, base_pix_sf)
databipic$trcat <- assign_trcat(databipic$bi_class)
databipic$category <- traj_to_cat[databipic$trcat]

# Attach the hydroregion id to each river pixel (via outl2)
databipic$HER <- GHR_riv$HydroRegions_raster_WGS84[match(databipic$outl2, GHR_riv$outl2)]

## 3.1b  Pixel-level bivariate data for the other drivers ------
## (needed for the "share of river network by driver" barplot)
LuFPix <- clean_pix(FloodTrendsP[FloodTrendsP$driver == driver[2], ])
LuDPix <- clean_pix(DroughtTrendsP[DroughtTrendsP$driver == driver[2], ])
ResFPix <- clean_pix(FloodTrendsP[FloodTrendsP$driver == driver[3], ])
ResDPix <- clean_pix(DroughtTrendsP[DroughtTrendsP$driver == driver[3], ])
WuFPix <- clean_pix(FloodTrendsP[FloodTrendsP$driver == driver[4], ])
WuDPix <- clean_pix(DroughtTrendsP[DroughtTrendsP$driver == driver[4], ])

databipilu <- pix_bidata(LuFPix, LuDPix, base_pix_sf)
databipire <- pix_bidata(ResFPix, ResDPix, base_pix_sf)
databipiwd <- pix_bidata(WuFPix, WuDPix, base_pix_sf)

# Combined socioeconomic signal (land use + reservoirs + water demand)
databipise <- databipire
databipise$x <- databipire$x + databipilu$x + databipiwd$x
databipise$y <- databipire$y + databipilu$y + databipiwd$y
databipise$bi_class <- make_bi_class(databipise$x, databipise$y)
databipise$bi_class[is.na(databipise$x) | is.na(databipise$y)] <- NA

for (obj in c("databipilu", "databipire", "databipiwd", "databipise")) {
  d <- get(obj)
  d$trcat <- assign_trcat(d$bi_class)
  assign(obj, d)
}

## 3.1c  Share of river network with a non-stable trajectory ----
## (per driver; same values shown on the barplot and in the table)
n_total_bar <- length(databipic$bi_class)

count_traj <- function(data, driver_name) {
  d <- data
  if (inherits(d, "sf")) st_geometry(d) <- NULL
  d %>%
    filter(!is.na(trcat), trcat != "Stable") %>%
    group_by(trcat) %>%
    summarise(n = dplyr::n(), .groups = "drop") %>%
    mutate(driver = driver_name, pct = n / n_total_bar * 100)
}

traj_counts <- bind_rows(
  count_traj(databipise, "All\nsocioeconomics"),
  count_traj(databipic, "Climate"),
  count_traj(databipilu, "Land use"),
  count_traj(databipire, "Reservoirs"),
  count_traj(databipiwd, "Water demand")
)

# Total % of non-stable river network per driver
bar_totals <- traj_counts %>%
  group_by(driver) %>%
  summarise(total_pct = sum(pct), .groups = "drop")

# Lookup: table row label -> driver name used in bar_totals
driver_pct <- function(driver_name) {
  v <- bar_totals$total_pct[bar_totals$driver == driver_name]
  if (length(v) == 0) NA_real_ else round(v[1])
}
pct_climate <- driver_pct("Climate")
pct_landuse <- driver_pct("Land use")
pct_reservoir <- driver_pct("Reservoirs")
pct_water <- driver_pct("Water demand")
pct_human <- driver_pct("All\nsocioeconomics") # combined human activities

## 3.2  Share of river network per category, within each HER ---
pix_df <- databipic
st_geometry(pix_df) <- NULL
pix_df <- pix_df[!is.na(pix_df$HER) & !is.na(pix_df$trcat), ]


# Count pixels per HER x category, then convert to a within-HER share
her_cat <- pix_df %>%
  group_by(HER, trcat) %>%
  summarise(n = dplyr::n(), .groups = "drop") %>%
  group_by(HER) %>%
  mutate(share = n / sum(n)) %>%
  ungroup()

# Dominant category per HER (plurality) and its share (for map fill + intensity)
her_dom <- her_cat %>%
  group_by(HER) %>%
  slice_max(order_by = share, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  dplyr::rename(dom_category = trcat, dom_share = share)

## 3.3  Join dominant category back to the HER polygons --------
databiclim <- st_transform(HydroRsf, crs = 3035)
databiclim$dom_category <- her_dom$dom_category[match(databiclim$Id, her_dom$HER)]
databiclim$dom_share <- her_dom$dom_share[match(databiclim$Id, her_dom$HER)]
# databiclim$trcat <- names(traj_to_cat)[match(databiclim$dom_category, traj_to_cat)]
# databiclim$category <- factor(databiclim$dom_category, levels = names(cat_colors))


# =============================================================
# 3b  BIOGEO-REGION BORDERS CLIPPED TO THE HYDROREGION DOMAIN
# =============================================================
# Mirrors section 3.8 of 05_Attribution_bivariate.R: merge the biogeo
# polygons (folding Pannonian & Steppic into Continental), then clip them
# to the union of the hydroregions we actually plot (databiclim).
domain_union <- st_union(st_make_valid(databiclim))

biogeo <- st_transform(biogeo, st_crs(databiclim))
biogeof_merged <- biogeo %>%
  mutate(
    name = if_else(name %in% c("Pannonian", "Steppic"), "Continental", name),
    pre_2012 = if_else(pre_2012 %in% c("PAN", "STE"), "CON", pre_2012)
  ) %>%
  group_by(name, pre_2012) %>%
  dplyr::summarize(geometry = st_union(geometry), .groups = "drop") %>%
  st_make_valid()

biogeof_clipped <- st_intersection(biogeof_merged, domain_union)
biogeof_clipped <- st_transform(biogeof_clipped, st_crs(databiclim))


# =============================================================
# 4  LEFT PANEL : MAP
# =============================================================

tsize <- 16
osize <- 12

# Shared map theme
map_theme_bi <- function(ts = tsize, os = osize) {
  theme(
    axis.title       = element_text(size = ts),
    panel.background = element_rect(fill = "aliceblue", colour = NA),
    # panel.border     = element_rect(linetype = "solid", fill = NA, colour = NA),
    legend.title     = element_text(size = ts),
    legend.text      = element_text(size = os),
    legend.position  = "bottom",
    # panel.grid.major = element_line(colour = "grey70"),
    # panel.grid.minor = element_line(colour = "grey90"),
    legend.key       = element_rect(fill = "transparent", colour = "transparent"),
    legend.key.size  = unit(.8, "cm")
  )
}


map_panel <- ggplot(basemap) +
  geom_sf(fill = "white") +
  geom_sf(
    data = databiclim,
    aes(fill = dom_category, geometry = geometry),
    color = NA, linewidth = 0.1
  ) +
  geom_sf(fill = NA, color = "gray42") +
  # Biogeo-region borders clipped to the hydroregion domain
  # geom_sf(
  #   data = biogeof_clipped,
  #   fill = NA, color = "grey25", linewidth = 0.35
  # ) +
  scale_fill_manual(values = cat_colors, na.value = "#e9ecef", guide = "none") +
  coord_sf(
    crs = 3035,
    xlim = c(min(nco[, 1]), max(nco[, 1])),
    ylim = c(min(nco[, 2]), max(nco[, 2]))
  ) +
  labs(
    title = "Effect of climate alone",
    subtitle = "Human activities are not included in this map"
  ) +
  theme_void() +
  theme(
    plot.title = element_text(size = 15, face = "bold", color = "#1f3b57"),
    plot.subtitle = element_text(size = 10, color = "#6c757d"),
    # plot.background = element_rect(fill = "#f4f7fb", color = NA),
    # panel.background = element_rect(fill = "aliceblue", color = NA)
  ) +
  map_theme_bi()

map_panel
# =============================================================
# 5  RIGHT PANEL : ARROW SUMMARY TABLE
# =============================================================
# Build a tidy table of rows; each row has a y position, label, colour box,
# and two arrow glyphs (flood / drought). We draw it with ggplot primitives
# so it integrates cleanly next to the map.

arrow_char <- function(dir) {
  switch(dir,
    up   = "\u2191", # <U+2191>  more severe
    down = "\u2193", # <U+2193>  less severe
    flat = "\u2014", # —  little effect
    ""
  )
}

# X positions for the columns
# Columns 1 (label/description) and 2 (% changed) are widened by pushing
# the % and arrow columns further right and extending the x-scale.
x_box <- 0.15
x_label <- 0.55
x_pct <- 8.6 # column 2: "River network changed" %
x_flood <- 10.6 # column 3: flood arrows
x_drought <- 12.1 # column 4: drought arrows

rows <- list()
yy <- 10

add_row <- function(rows, y, type, title = "", desc = "",
                    flood = NA, drought = NA, color = NA, pct = NA) {
  rows[[length(rows) + 1]] <- data.frame(
    y = y, type = type, title = title, desc = desc,
    flood = ifelse(is.na(flood), "", flood),
    drought = ifelse(is.na(drought), "", drought),
    color = ifelse(is.na(color), NA, color),
    pct = ifelse(is.na(pct), NA_real_, pct),
    stringsAsFactors = FALSE
  )
  rows
}

# Section header: Climate (with % of river network showing a non-stable change)
rows <- add_row(
  rows, yy, "header", "Climate",
  "The largest driver. Shown on the map",
  pct = pct_climate
)
yy <- yy - 1

for (i in seq_len(nrow(cat_arrows))) {
  cat <- cat_arrows$category[i]
  rows <- add_row(
    rows, yy, "cat", cat, cat_desc[[cat]],
    cat_arrows$flood[i], cat_arrows$drought[i], cat_colors[[cat]]
  )
  yy <- yy - 1
}

yy <- yy - 0.3
# Section header: Human activities (with combined % of river network changed)
rows <- add_row(
  rows, yy, "header", "Human activities",
  "Driving consistent local trajectories",
  pct = pct_human
)
yy <- yy - 1

# Map each human-activity row to its driver percentage
human_pct <- c(
  "Land use"            = pct_landuse,
  "Dams and reservoirs" = pct_reservoir,
  "Water use"           = pct_water
)
for (i in seq_len(nrow(human_arrows))) {
  act <- human_arrows$activity[i]
  rows <- add_row(
    rows, yy, "human", act,
    human_arrows$desc[i], human_arrows$flood[i], human_arrows$drought[i], NA,
    pct = unname(human_pct[act])
  )
  yy <- yy - 1
}

tbl <- do.call(rbind, rows)
tbl$flood_glyph <- vapply(tbl$flood, function(d) if (d == "") "" else arrow_char(d), "")
tbl$drought_glyph <- vapply(tbl$drought, function(d) if (d == "") "" else arrow_char(d), "")

headers <- tbl[tbl$type == "header", ]
cats <- tbl[tbl$type == "cat", ]
humans <- tbl[tbl$type == "human", ]
arrows_all <- tbl[tbl$type %in% c("cat", "human"), ]

# Rows that carry a driver percentage (Climate header + human activities)
pct_rows <- tbl[!is.na(tbl$pct), ]
pct_rows$pct_label <- paste0(round(pct_rows$pct), "%")

table_panel <- ggplot() +
  # Legend row (arrow meanings) - near the title, larger
  annotate("text",
    x = x_label, y = 11.9, hjust = 0,
    label = "\u2191 more severe      \u2193 less severe      \u2014 little effect",
    size = 5.2, color = "#495057"
  ) +
  # Column headers
  annotate("text",
    x = x_pct, y = 11, label = "River\nnetwork\nchanged",
    fontface = "bold", size = 4, color = "#2b2b2b", lineheight = 0.9
  ) +
  annotate("text",
    x = x_flood, y = 11, label = "Floods",
    fontface = "bold", size = 4, color = "#2b2b2b"
  ) +
  annotate("text",
    x = x_drought, y = 11, label = "Droughts",
    fontface = "bold", size = 4, color = "#2b2b2b"
  ) +
  # Section headers
  geom_text(
    data = headers, aes(x = x_box, y = y, label = title),
    hjust = 0, fontface = "bold", size = 5.8, color = "#2b2b2b"
  ) +
  geom_text(
    data = headers, aes(x = x_box, y = y - 0.33, label = desc),
    hjust = 0, size = 4.1, color = "#6c757d"
  ) +
  # Climate colour chips
  geom_tile(
    data = cats, aes(x = x_box + 0.1, y = y, fill = color),
    width = 0.3, height = 0.45
  ) +
  scale_fill_identity() +
  # Category / activity titles + descriptions
  geom_text(
    data = arrows_all, aes(x = x_label, y = y, label = title),
    hjust = 0, fontface = "bold", size = 5, color = "#212529"
  ) +
  geom_text(
    data = arrows_all, aes(x = x_label, y = y - 0.33, label = desc),
    hjust = 0, size = 3.9, color = "#6c757d"
  ) +
  # Arrow glyphs
  geom_text(
    data = arrows_all, aes(x = x_flood, y = y, label = flood_glyph),
    size = 8, fontface = "bold", color = "#2b2b2b"
  ) +
  geom_text(
    data = arrows_all, aes(x = x_drought, y = y, label = drought_glyph),
    size = 8, fontface = "bold", color = "#2b2b2b"
  ) +
  # % of river network with a non-stable trajectory, next to each driver
  geom_text(
    data = pct_rows, aes(x = x_pct, y = y, label = pct_label),
    fontface = "bold", size = 4.2, color = "#2b2b2b"
  ) +
  scale_x_continuous(limits = c(0, 12.8)) +
  scale_y_continuous(limits = c(min(tbl$y) - 0.6, 12.3)) +
  labs(
    title = "Effect of each driver",
    subtitle = "Arrows show the direction of change in flood and drought severity"
  ) +
  theme_void() +
  theme(
    plot.title = element_text(size = 15, face = "bold", color = "#1f3b57"),
    plot.subtitle = element_text(size = 10, color = "#6c757d"),
    plot.background = element_rect(fill = NA, color = NA),
    panel.background = element_rect(fill = NA, color = NA)
  )

table_panel

# =============================================================
# 6  BARPLOT – % of river network per trajectory, ranked by driver
# =============================================================
# Adapted from 05_Attribution_bivariate.R. The per-driver counts
# (traj_counts, bar_totals) are computed once in section 3.1c and reused
# both here and in the table panel.

# Stack order (meaningful trajectory sequence) + colours from cat_colors
traj_order <- c("Drying", "Decelerating", "Wetting", "Accelerating")
traj_counts$trcat <- factor(traj_counts$trcat, levels = traj_order)

# Rank drivers by total classified share (ascending -> highest on top after flip)
driver_order <- traj_counts %>%
  group_by(driver) %>%
  summarise(total_pct = sum(pct), .groups = "drop") %>%
  arrange(total_pct) %>%
  pull(driver)

fig_bv_bar <- ggplot(
  traj_counts,
  aes(x = driver, y = pct, fill = trcat)
) +
  geom_col(width = 0.65, color = "transparent") +
  geom_text(
    data = bar_totals,
    aes(x = driver, y = total_pct, label = paste0(round(total_pct, 0), "%")),
    inherit.aes = FALSE,
    hjust = -0.15, size = 5, fontface = "bold", color = "grey20"
  ) +
  scale_fill_manual(values = cat_colors, name = NULL, drop = FALSE) +
  scale_y_continuous(
    name   = "Share of river network (%)",
    limits = c(0, max(bar_totals$total_pct) * 1.2),
    expand = c(0, 0)
  ) +
  scale_x_discrete(limits = driver_order, name = NULL) +
  geom_hline(yintercept = 0, color = "black", linewidth = 1) +
  coord_flip() +
  theme(
    axis.title.x = element_text(size = 16, face = "bold"),
    axis.text = element_text(size = 14),
    axis.text.y = element_text(face = "bold"),
    axis.text.x = element_blank(),
    axis.ticks = element_blank(),
    panel.background = element_rect(fill = "white"),
    panel.grid.major = element_blank(),
    panel.border = element_blank(),
    legend.position = "none"
  )

fig_bv_bar

ggsave(file.path(plotDir, "Fig_bv_bar_by_driver.jpg"),
  fig_bv_bar,
  width = 18, height = 14, units = "cm", dpi = 600
)


# =============================================================
# 7  ASSEMBLE & SAVE
# =============================================================
caption_panel <- ggplot() +
  annotate("text",
    x = 0, y = 0.6, hjust = 0,
    label = "Climate sets the big picture. Human choices on land, dams and water use shape what happens locally.",
    fontface = "bold", size = 4, color = "white"
  ) +
  annotate("text",
    x = 0, y = 0.2, hjust = 0,
    label = "Source: Tilloy et al. Map shows climate-driven change only, simplified from the study's results.",
    size = 2.9, color = "#c7d2de"
  ) +
  scale_x_continuous(limits = c(0, 10)) +
  scale_y_continuous(limits = c(0, 0.9)) +
  theme_void() +
  theme(
    plot.background = element_rect(fill = "#1f3b57", color = NA),
    panel.background = element_rect(fill = "#1f3b57", color = NA)
  )

infographic <- (map_panel | table_panel) +
  plot_layout(widths = c(.9, 1)) &
  theme(
    plot.background = element_rect(fill = NA, color = NA),
    plot.margin = margin(t = 8, r = 8, b = 8, l = 8)
  )

infographic

ggsave(file.path(plotDir, "Fig_infographic_floods_droughts.jpg"),
  infographic,
  width = 34, height = 20, units = "cm", dpi = 400
)
ggsave(file.path(plotDir, "Fig_infographic_floods_droughts.pdf"),
  infographic,
  width = 34, height = 18, units = "cm"
)

message("Infographic saved to: ", file.path(plotDir, "Fig_infographic_floods_droughts.jpg"))
