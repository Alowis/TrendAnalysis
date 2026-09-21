
<!-- README.md is generated from README.Rmd. Please edit that file -->

# Codes of “Attribution of changes in European flood and drought magnitudes to climatic and socioeconomic drivers”

<!-- badges: start -->
<!-- badges: end -->

This repository contains the codes used to analyse discharge data from
the HERA hydrological reanalysis and to reproduce the different steps
and figures in the article: “Attribution of changes in European flood
and drought magnitudes to climatic and socioeconomic drivers”.

## Introduction

The HERA high-resolution pan-European hydrological reanalysis
(1951-2020) dataset
[10.2905/a605a675-9444-4017-8b34-d66be5b18c95](10.2905/a605a675-9444-4017-8b34-d66be5b18c95)
is the result of a joint effort between the JRC and PIK to produce a
long term hydrological reanalysis with downscaled and bias-corrected
climate reanalysis (ERA5-land) and dynamic socioeconomic inputs. It
includes maps of climate variables (evaporation, evapotranspiration,
precipitation, temperature), dynamic socioeconomic inputs (land use,
water demand, reservoir maps) required for hydrological modelling with
LISFLOOD (<https://github.com/ec-jrc/lisflood-code>) and river discharge
with European extent at 1 arc minute (~1.5 km) grid resolution and
6-hourly time step.

The analysis buids on the HERA dataset and three counterfactual scenarii
to attribute changes in hydrological extremes to four drivers:

- Climate (change and variability)
- Reservoir construction
- Land use changes (six land use categories)
- Water demand changes

## Content

The repository is composed of R scripts under the `/R` folder. Inside
the `/R` folder, the numbered scripts represent the successive steps of
the analysis:

1.  [Threshold definition for non-stationary
    analysis](#thresold-definition) — `01_TSEVA_TrendThresholdSel.R`
2.  [Non-stationary EVA (TSEVA)](#tseva) — `02a_TSEVA_TrendxVar.R` and
    `02b_TSEVA_Run_SCF.R`
3.  [Trend and variability analysis](#trendvar) —
    `03_TrendVar_analysis.R`
4.  [Changing Hydrological Extremes - univariate attribution](#che-u) —
    `04_Attribution_univariate.R`
5.  [Changing Hydrological Extremes - bivariate attribution](#che-b) —
    `05_Attribution_bivariate.R`
6.  [Intra-catchment trajectories](#intracatch) —
    `06_IntraCatchment_trajs.R`
7.  [Functions](#functions) — `functions_trends.R`

Two additional helper files support all scripts:

- `functions_trends.R`: the shared library of analysis functions,
  sourced at the top of every script.
- `config_paths.R`: a single source of truth for all data and output
  paths (see [Paths and configuration](#paths)).

A set of supplementary scripts is also provided under
`/Pre-PostProcess`, notably scripts to determine the frost season of
European catchments and to map reservoirs used in HERA. Superseded
versions of the scripts are kept under `/R/old`.

## Paths and configuration <a id="paths"></a>

All scripts resolve their own location automatically (they work both
interactively in RStudio and via `Rscript` on an HPC) and then source
`config_paths.R`, which defines every data and output directory from a
single project root.

The canonical data directory is
`<PROJECT_ROOT>/ChangingHydroExtremes/data`. The project root defaults
to a local path but can be overridden with the `TREND_PROJECT_ROOT`
environment variable so the same code runs unchanged locally and on the
HPC:

    # HPC example
    export TREND_PROJECT_ROOT=/scratch/.../LFRuns_utils

`config_paths.R` also exposes the derived sub-directories (`geoDir`,
`threshDir`, `trendVarDir`, `droughtDir`, `floodDir`, `riverDir`,
`resDir`, `plotDir`) and the pixel-ID multiplier `ID_MULT`, which is
shared across scripts so that identifiers produced in step 1 line up
with those used in the later steps.

The HPC-oriented scripts (`01`, `02a`, `02b`) also accept command-line
arguments, e.g.:

    Rscript 01_TSEVA_TrendThresholdSel.R <Nsq> <tail> <sce>
    Rscript 02a_TSEVA_TrendxVar.R <Nsq> <haz> <sce>
    Rscript 02b_TSEVA_Run_SCF.R "<Nsq> <haz> <sce> <startid> <endid>"

### 1. Threshold definition <a id="thresold-definition"></a>

#### Script: *01_TSEVA_TrendThresholdSel.R*

#### Script Description

This script identifies thresholds for extreme trend assessment using the
TSEVA (Time Series Extreme Value Analysis) method. The method is applied
to a subsample of the 282,521 river pixels analyzed in a research
article. The script processes discharge data from the HERA (Hydrological
Ensemble Prediction System) and performs threshold identification for
trend analysis on selected scenarios.

#### Script Overview

- **Set Working Directory and Load Functions**:
  - Resolves the script location (works interactively and under
    `Rscript`).
  - Sources `functions_trends.R` and `config_paths.R`.
- **Define Parameters**:
  - Reads the square number (`Nsq`), tail type (`tail`) and scenario
    (`sce`) from the command line (HPC) or from local defaults.
  - Loads metadata about the spatial division of the domain into chunks.
- **Load River Pixel Locations**:
  - Loads river pixel locations and identifiers using the `outletopen`
    function.
- **Load Discharge Data**:
  - Loads discharge data for two scenarios: historical run (`Histo`) and
    socioeconomic counterfactual (`SCF`).
  - Extracts discharge data for the first river pixel.
- **Prepare Data for Analysis**:
  - Prepares the discharge data for trend analysis, including time
    stamps and series extraction.
  - Registers a parallel backend for parallel processing.
- **Parallel Processing**:
  - Iterates over the river pixels using a parallel loop to find the
    trend threshold for each pixel.
  - Saves the results to CSV files (`trenTH_x_*` and `xtrempoints_*`)
    under `threshDir`.

#### Key Functions Used

- **`outletopen`**: Opens netCDF outlet files and extracts relevant
  data.
- **`disNcopenloc`**: Extracts discharge time series data for a specific
  location.
- **`tsEvaFindTrendThreshold2`**: Finds the optimal extreme threshold
  for trend estimation.
- **`foreach` and `doParallel`**: Used for parallel processing of the
  data.

#### Output

The script outputs CSV files containing the trend thresholds and the
extreme points for the specified scenarios and tail types.

### 2. Non-stationary EVA (TSEVA) <a id="tseva"></a>

#### Scripts: *02a_TSEVA_TrendxVar.R* and *02b_TSEVA_Run_SCF.R*

This step is split into two scripts. Both are designed to run one square
of the HERA domain on a High-Performance Computing (HPC) cluster and
perform extreme trend assessment for either drought or flood hazards
using the TSEVA (Time Series Extreme Value Analysis) method.

#### 2a. Trend and variability per variable — *02a_TSEVA_TrendxVar.R*

This script extracts, for each river pixel of a square, the long-term
trend and variability of the extremes using the thresholds produced in
step 1. It reuses the extreme points identified in the socioeconomic
counterfactual so that all scenarios share a common set of events.

**Overview**

1.  **Setup**: resolves the script location, sources
    `functions_trends.R` and `config_paths.R`, and reads `Nsq`, `haz`
    and `sce` from the command line or local defaults.
2.  **Load river pixel locations** with `outletopen` and upstream areas
    with `UpAopen`.
3.  **Scenario differentiation** (Historical, SCF, Water demand,
    Reservoir+Water demand).
4.  **Load discharge data** for the scenario and square.
5.  **Load frost data** (for drought) and read the thresholds/extreme
    points from step 1.
6.  **Loop over pixels**: compute the trend and standard-deviation
    series on the extremes.
7.  **Save** the trend and variability series to `trendVarDir`
    (`TrendVarX_*.Rdata`).

#### 2b. Full non-stationary EVA run — *02b_TSEVA_Run_SCF.R*

This script performs the complete non-stationary EVA fit for one square,
producing GEV/GPD parameters, return levels and return periods per pixel
and per year.

**Overview**

1.  **Setup**: sources `config_paths.R` (the data root) and reads
    arguments (`Nsq`, `haz`, `sce`, `startid`, `endid`).
2.  **Load river pixel locations** with `outletopen`.
3.  **Scenario differentiation**.
4.  **Load discharge data** for the scenario and square.
5.  **Load frost data** (for drought analysis).
6.  **Threshold identification**: reads the SCF and Historical threshold
    files from step 1 and retains the SCF threshold unless it is NA.
7.  **Loop over pixels**: detects intermittent rivers, runs `TsEvaNs`,
    and computes return periods and levels.
8.  **Save** the results (parameters, return levels, return periods,
    peaks) under the hazard folder.

#### Key Functions Used

- **`outletopen`**: Opens netCDF outlet files and extracts relevant
  data.
- **`disNcopenloc`**: Extracts discharge time series data for a specific
  location.
- **`interid`**: Detects intermittent rivers.
- **`TsEvaNs`**: Performs non-stationary extreme value analysis.
- **`ComputeReturnLevels`**: Computes return levels associated with a
  return period.
- **`RPcalc`**: Calculates return periods associated with a return
  level.

#### Output

`02a` outputs the trend and variability series per square. `02b` outputs
a list containing parameters, return levels, return periods, peaks, and
catchment rest data for further analysis.

### 3. Trend and variability analysis <a id="trendvar"></a>

#### Script: *03_TrendVar_analysis.R*

#### Script Description

This script assembles the per-square trend and variability outputs from
step 2 across all scenarios, attributes changes in the location, scale
and return-level parameters to the four drivers, and produces the
associated pan-European maps and aggregated trend plots.

#### Key Steps

1.  **Setup**: resolves the script location and sources
    `functions_trends.R` and `config_paths.R`.
2.  **Spatial data**: reads outlets for all squares, HydroBASINS level 7
    catchments, European biogeographic regions, HydroEcoRegions, the
    base map and upstream areas.
3.  **Load scenario data** with `process_hazard_data()` for SCF,
    Historical, WStat and RWStat.
4.  **Merge scenarios** into a single `df_Main` with
    `attach_scenario_cols()`.
5.  **Driver analysis** with `analyze_hazard_symmetric()` for the
    location parameter (threshold), scale parameter (Sigma) and return
    level (RP).
6.  **Aggregate global trends** by HydroRegion with
    `aggregate_global_trends()`.
7.  **Save** driver maps and aggregation plots to `plotDir`, and the
    return-level matrices under the hazard folder.

#### Key Functions Used

- **`outletopen`**, **`UpAopen`**, **`ReservoirOpen`**: spatial data
  extraction from netCDF files.
- **`process_hazard_data`**: loads and reshapes the per-scenario
  trend/variability and parameter files.
- **`attach_scenario_cols`**: joins scenario-specific parameters onto
  the main table.
- **`analyze_hazard_symmetric`**: computes symmetric relative changes
  per driver and variable.
- **`aggregate_global_trends`**: weighted aggregation of trends by
  HydroRegion.
- **`calcGPDReturnLevel_Single`**: vectorised GPD/Gumbel return-level
  computation.

#### Output

Driver maps (threshold, Sigma, RP) per driver and hazard, global trend
aggregation plots, and saved return-level matrices for each scenario.

### 4. Changing Hydrological Extremes - univariate attribution <a id="che-u"></a>

#### Script: *04_Attribution_univariate.R*

#### Script Description

This script is designed to load and process pre-loaded results for
analyzing river discharge data across various scenarios and hazards
(flood and drought). It performs several key steps, including loading
spatial data, fitting results, and computing changes in return levels
due to different drivers such as climate, land use, reservoirs, and
water demand. The script also handles data cleaning, spatial smoothing,
and aggregation for further analysis and visualization.

#### Key Steps

1.  **Set Working Directory and Load Functions**:
    - Resolves the script location and sources `functions_trends.R` and
      `config_paths.R`.
2.  **Load Outlets Data**:
    - Reads and processes outlet data for all squares in the HERA
      domain.
    - Loads river pixel locations and identifiers.
3.  **Load Spatial Data for Catchments**:
    - Loads metadata about the domain’s spatial division into chunks.
    - Loads spatial data for catchments, including Hybas07 and European
      Biogeo regions.
4.  **Load Fitting Results**:
    - Loads fitting results from various runs (Historical, Socioeconomic
      Counterfactual, Water Demand Counterfactual, Reservoir+Water
      Demand Counterfactual) for all river pixels in the domain.
5.  **Data Cleaning**:
    - Cleans the data by identifying and handling missing values,
      especially for the last years’ return levels.
6.  **Load IRES Status**:
    - Loads the Intermittent River Status (IRES) for drought analysis.
    - Combines IRES data from different runs to detect any river pixel
      that is intermittent in at least one run.
7.  **Drought Return Level Corrections**:
    - Applies corrections to drought return levels, including reversing
      values and setting negative levels to zero.
8.  **Remove Irrealistic Shape Parameters**:
    - Identifies and removes pixels with unrealistic shape parameters to
      ensure data quality.
9.  **Compute Large-Scale Errors**:
    - Computes large-scale errors in return levels for different years
      and scenarios.
10. **Generate Plots**:
    - Generates various plots, including maps of return level errors,
      changes in hazard intensity, and shape parameter instability.
11. **Change Attribution**:
    - Computes changes in return levels attributed to different drivers
      (climate, land use, reservoirs, water demand).
    - Aggregates changes by regions and performs spatial smoothing to
      remove noise from unstable Generalized Pareto Distribution (GPD)
      fits.
12. **Save Outputs**:
    - Saves the processed data and results for further analysis and
      visualization.

#### Key Functions Used

- **`outletopen`**: Opens netCDF outlet files and extracts relevant data
  based on specified parameters.
- **`UpAopen`**: Opens upstream area files for selected pixels and
  returns a data frame with upstream area information.
- **`disNcopenloc`**: Extracts discharge time series data for a specific
  location from netCDF files.
- **`ReservoirOpen`**: Opens reservoir location files and returns a data
  frame with reservoir location information.
- **`RPcalc`**: Calculates return periods associated with a return level
  using Generalized Extreme Value (GEV) and Generalized Pareto
  Distribution (GPD) parameters.
- **`RPchangeCal`**: Computes changes in return periods associated with
  a return level for different years and distribution laws.
- **`interid`**: Detects intermittent rivers by identifying periods of
  low or zero discharge and classifies them based on specified criteria.
- **`check_timeserie2`**: Checks if there are not gaps bigger than two
  years in extreme value time series.
- **`tsEvaFindTrendThreshold2`**: Finds the optimal extreme threshold
  for trend estimation by iterating over threshold values and assessing
  the stability of the trend.
- **`ComputeReturnLevels`**: Computes return levels associated with a
  return period using non-stationary extreme value analysis parameters.
- **`GPDLargeRLs`**: Computes return levels for any return period using
  Generalized Pareto Distribution (GPD) parameters.
- **`calculate_return_levels`**: Calculates return levels and their
  errors based on Generalized Pareto Distribution (GPD) parameters.
- **`weighted_average`**: Calculates the weighted average of points
  within a maximum distance from a specific point.
- **`neighbour_finder`**: Finds neighboring points within a maximum
  distance from a specific point.
- **`ComputeChange`**: Computes changes in driver trends and normalizes
  data by area for specified regions and years.
- **`calculatePoints`**: Calculates significant points and their trends,
  and returns a list with significant points data and grid points.
- **`processTrendData`**: Processes trend data and aggregates it by
  decade and location for further analysis.
- **`UpATrendData`**: Processes trend data for upstream areas and
  aggregates it by upstream area groups.
- **`calculateTrendSig`**: Calculates trend significance using the
  Mann-Kendall test and categorizes changes based on significance
  levels.
- **`get_density`**: Calculates the density of points using kernel
  density estimation for specified grid points.

#### Output

The script performs a comprehensive analysis of river discharge data and
generates several outputs, which can be categorized into intermediate
data processing results and final visualizations. Here are the key
outputs of the script:

##### Intermediate Data Processing Results

1.  **Processed Outlet Data**:
    - `outf`: A data frame containing outlet information for all squares
      in the HERA domain.
2.  **Spatial Data for Catchments**:
    - `Catchmentrivers7`: A data frame containing metadata about
      catchment rivers.
    - `hybas07`: A spatial data frame containing hydrological basin
      data.
    - `Catamere07`: A data frame combining catchment metadata with
      hydrological basin data.
    - `GNF`: A data frame with joined catchment and outlet data.
    - `GHR_riv`: A data frame with HydroRegions data.
3.  **Fitting Results**:
    - `ParamsflH`, `ParamsflSCF`, `ParamsflRWCF`, `ParamsflWCF`: Data
      frames containing fitting parameters for different scenarios
      (Historical, Socioeconomic Counterfactual, Reservoir+Water Demand
      Counterfactual, Water Demand Counterfactual).
    - `PeakH`, `PeakSCF`, `PeakRWCF`, `PeakWCF`: Data frames containing
      peak data for different scenarios.
    - `RLGPDflH`, `RLGPDflSCF`, `RLGPDflRWCF`, `RLGPDflWCF`: Data frames
      containing return level data for different scenarios.
4.  **IRES Status**:
    - `IRES_Histo`, `IRES_SocCF`, `IRES_WCF`, `IRES_RWCF`: Data frames
      containing Intermittent River Status (IRES) for different
      scenarios.
    - `IRES_comb`: A combined data frame of IRES status from all
      scenarios.
5.  **Corrected Return Levels**:
    - Corrected return level data frames for drought analysis
      (`RLGPDflH`, `RLGPDflSCF`, `RLGPDflRWCF`, `RLGPDflWCF`).
6.  **Data filtering**:
    - `rmpixs`: A vector of pixel indices with unrealistic shape
      parameters.
    - `Shapeparf`: A data frame with filtered shape parameters.
7.  **Large-Scale Error Computation**:
    - `RlevErrtH`: A data frame containing large-scale error
      computations for return levels.
8.  **Change Attribution**:
    - `Climtrend`, `Soctrend`, `Restrend`, `Wutrend`, `Totaltrend`: Data
      frames containing changes attributed to different drivers
      (climate, land use, reservoirs, water demand, total).
9.  **Spatial Aggregation**:
    - `DataL`, `DataW`, `DataT`, `DataC`, `DataR`: Data frames
      containing spatially aggregated changes for different drivers.
10. **Aggregated Changes**:
    - `pointSoc`, `pointWu`, `pointTot`, `pointClim`, `pointRes`: Data
      frames containing aggregated changes at the HydroRegion level.

##### Final Visualizations

1.  **Plots**:
    - Maps of return level errors.
    - Changes in hazard intensity with confidence intervals.
    - Relative error of 10-year return levels in 1955.
    - Comparison of error and 2015-1955 changes.
    - Intermittent rivers plot.
    - Mean 10-year return level in specific discharge.
    - Lower bound of fitted GPD (for low flows).
    - Shape parameter instability check at the regional level.
    - Shape parameter plot.
    - Ordered change aggregated at the HydroRegion level.
    - Maps of changes in 10-year return levels driven by different
      drivers.
    - Boxplots of changes in time by biogeographical regions and
      catchment size.
2.  **Saved Outputs**:
    - `Output_fl_year` and `Output_dr_nonfrost`: Saved lists containing
      trend data at the pixel and regional levels, output for 2020, and
      initial data for flood and drought scenarios, respectively.

### 5. Changing Hydrological Extremes - bivariate attribution <a id="che-b"></a>

#### Script: *05_Attribution_bivariate.R*

#### Script Description

This script performs bivariate analysis on river discharge data for
flood and drought hazards across Europe. It loads pre-loaded results
from previous analyses, categorizes changes at both the catchment and
pixel levels, and generates various visualizations to illustrate the
changes in hydrological extremes. The script also performs bivariate
categorization for different drivers (climate, land use, reservoirs,
water demand) and combines these drivers to analyze their joint effects
on hydrological changes.

### Key Steps

1.  **Library Calling and Data Loading**:
    - Resolves the script location and sources `functions_trends.R` and
      `config_paths.R`.
    - Loads pre-loaded results for flood and drought hazards from
      previous analyses.
2.  **Bivariate Results Plot**:
    - Aggregates and plots total changes in time for flood and drought
      hazards.
    - Generates a boxplot of changes in flood and drought hazards by
      decade.
3.  **Spatial Data Loading**:
    - Loads biogeographic regions and hydrological basin data.
    - Loads HydroRegions data and prepares plot parameters.
4.  **Trend Significance Recomputation**:
    - Recomputes trend significance for flood and drought data.
5.  **Bivariate Categorization**:
    - Categorizes changes at the catchment and pixel levels for total,
      climate, land use, reservoir, and water demand drivers.
    - Generates bivariate plots showing the combined effects of
      different drivers on hydrological changes.
6.  **Aggregation by Biogeoregion**:
    - Aggregates changes by biogeographic regions and generates stacked
      barplots of change trajectories.
7.  **Driver Contribution Analysis**:
    - Analyzes the contribution of different drivers to hydrological
      changes.
    - Generates bivariate plots of driver contributions to hydro-extreme
      changes.
8.  **Temporal Evolution and Significance**:
    - Analyzes the temporal evolution of changes per region.
    - Generates plots showing the significance of trends and the
      proportion of trajectories by drivers.
9.  **Output Saving**:
    - Saves the processed data and results for further analysis and
      visualization.

### Key Functions Used

- **`outletopen`**: Opens netCDF outlet files and extracts relevant data
  based on specified parameters.
- **`calculateTrendSig`**: Calculates trend significance using the
  Mann-Kendall test and categorizes changes based on significance
  levels.
- **`bi_legend`**: Generates a legend for bivariate plots.
- **`bi_pal`**: Generates a color palette for bivariate plots.
- **`bi_scale_fill`**: Scales fill colors for bivariate plots.
- **`bi_scale_color`**: Scales colors for bivariate plots.
- **`ggarrange`**: Arranges multiple ggplot objects in a grid.
- **`ggplot`**: Generates plots using the ggplot2 package.
- **`geom_sf`**: Adds spatial features to ggplot objects.
- **`scale_fill_manual`**: Manually scales fill colors in ggplot.
- **`scale_color_manual`**: Manually scales colors in ggplot.
- **`coord_sf`**: Sets coordinate reference systems for spatial data in
  ggplot.
- **`scale_size`**: Scales the size of points in ggplot.
- **`theme`**: Customizes the appearance of ggplot objects.

### Outputs

1.  **Plots**:
    - Line plot showing the total changes in time for flood and drought
      hazards.
    - Boxplot showing changes in flood and drought hazards by decade.
    - Bivariate plots showing the combined effects of different drivers
      on hydrological changes at both the catchment and pixel levels.
    - Stacked barplots of change trajectories by biogeographic regions.
    - Bivariate plots showing the contribution of different drivers to
      hydro-extreme changes.
    - Plots showing the temporal evolution of changes per region.
    - Plots showing the significance of trends and the proportion of
      trajectories by drivers.
2.  **Saved outputs**:
    - `mbfX`: Data frame containing aggregated changes by HydroRegions
      for different drivers and `mbfH`: Data frame containing historical
      changes by HydroRegions.
    - Bivariate change objects (`CLchanges_bivariate_F.Rdata`,
      `SEchanges_bivariate_F.Rdata`) saved under the `Trajectories`
      folder for use in step 6.

### 6. Intra-catchment trajectories <a id="intracatch"></a>

#### Script: *06_IntraCatchment_trajs.R*

#### Script Description

This script produces a detailed intra-catchment case study (the Rhône
basin by default) that overlays the climatic and socioeconomic drivers
of change on a single catchment. It combines a digital elevation model,
an aridity-index change layer, land-use and water-demand changes, new
reservoirs, and the bivariate hydrological trajectories from step 5 into
two side-by-side maps.

#### Key Steps

1.  **Setup**: resolves the script location and sources
    `functions_trends.R` and `config_paths.R`.
2.  **Extract the target catchment** and its sub-catchments from the
    HydroBASINS layers.
3.  **Prepare the DEM** (crop, reproject to LAEA, hillshade).
4.  **Load / compute the aridity-index change** raster.
5.  **Load land use, water demand (NUTS3) and reservoir layers** and
    clip them to the catchment.
6.  **Load the bivariate trajectories** (`CLchanges_bivariate_F.Rdata`,
    `SEchanges_bivariate_F.Rdata`) and crop them to the catchment.
7.  **Build the maps**: a climatic-changes panel and a
    socioeconomic-changes panel, combined with a bivariate legend.
8.  **Save** the combined figure to `plotDir`.

#### Key Functions Used

- **`outletopen`**, **`resOpen`**: spatial data extraction from netCDF
  files.
- **`sens.slope`** / **`mmkh`** (from `modifiedmk`/`trend`): pixel-wise
  trend estimation.
- **`exact_extract`** (from `exactextractr`): zonal statistics for water
  demand.
- Standard `terra`/`raster`/`sf`/`ggplot2` spatial and plotting
  functions.

#### Output

A combined multi-panel JPEG (climatic and socioeconomic changes for the
case-study catchment) saved under `plotDir`.

### 7. Functions <a id="functions"></a>

#### 7.1 Miscellaneous Functions

##### 7.1.1 Functions for TSEVA

- **`check_timeserie2`**: Checks if there are not gaps bigger than two
  years in extreme value time series.

- **`tsEvaFindTrendThreshold2`**: Finds the optimal extreme threshold
  for trend estimation by iterating over threshold values and assessing
  the stability of the trend.

- **`outletopen`**: Opens netCDF outlet files and extracts relevant data
  based on specified parameters.

- **`ComputeReturnLevels`**: Computes return levels associated with a
  return period using non-stationary extreme value analysis parameters.

- **`GPDLargeRLs`**: Computes return levels for any return period using
  Generalized Pareto Distribution (GPD) parameters.

##### 7.1.2 Data Extraction Functions

- **`UpAopen`**: Opens upstream area files for selected pixels and
  returns a data frame with upstream area information.

- **`ReservoirOpen`**: Opens reservoir location files and returns a data
  frame with reservoir location information.

- **`disNcopen`**: Extracts discharge time series data for selected
  locations from netCDF files.

- **`disNcopenloc`**: Extracts discharge time series data for a specific
  location from netCDF files.

##### 7.1.3 Miscellaneous Calculations

- **`RPcalc`**: Calculates return periods associated with a return level
  using Generalized Extreme Value (GEV) and Generalized Pareto
  Distribution (GPD) parameters.

- **`RPchangeCal`**: Computes changes in return periods associated with
  a return level for different years and distribution laws.

- **`interid`**: Detects intermittent rivers by identifying periods of
  low or zero discharge and classifies them based on specified criteria.

#### 7.2 Univariate Analysis

- **`calculate_return_levels`**: Calculates return levels and their
  errors based on Generalized Pareto Distribution (GPD) parameters.

- **`weighted_average`**: Calculates the weighted average of points
  within a maximum distance from a specific point.

- **`neighbour_finder`**: Finds neighboring points within a maximum
  distance from a specific point.

- **`ComputeChange`**: Computes changes in driver trends and normalizes
  data by area for specified regions and years.

- **`calculatePoints`**: Calculates significant points and their trends,
  and returns a list with significant points data and grid points.

- **`processTrendData`**: Processes trend data and aggregates it by
  decade and location for further analysis.

- **`UpATrendData`**: Processes trend data for upstream areas and
  aggregates it by upstream area groups.

#### 7.3 Bivariate Analysis

- **`calculateTrendSig`**: Calculates trend significance using the
  Mann-Kendall test and categorizes changes based on significance
  levels.

- **`get_density`**: Calculates the density of points using kernel
  density estimation for specified grid points.
