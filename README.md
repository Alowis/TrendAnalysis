
<!-- README.md is generated from README.Rmd. Please edit that file -->

# Codes of “Climatic and socioeconomic drivers of changing hydrological extremes in Europe”

<!-- badges: start -->
<!-- badges: end -->

This repository contains the codes used to analyse discharge data from
the HERA hydrological reanalysis and to reproduce the different steps
and figures in the article: “Climatic and socioeconomic drivers of
changing hydrological extremes in Europe”.

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

The repository is composed of R scripts under the /R folder. Inside the
/R folder, one can find 5 scripts representing the 4 main steps of the
analysis:

1.  [Thresold definition for non-stationnary
    analysis](#thresold-definition)
2.  [Non-stationnary EVA (TSEVA)](#tseva)
3.  [Changing Hydrological Extremes - univariate analysis](#che-u)
4.  [Changing Hydrological Extremes - bivariate analysis](#che-b)
5.  [Functions](#functions)

A set of supplementary functions is also provided under
/Pre-PostProcess, notably functions to determine the frost season of
European catchments and to maps reservoirs used in HERA.

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
  - Sets the working directory to the location of the script.
  - Sources the `functions_trends.R` file to load necessary functions.
- **Set Data Directory**:
  - Specifies the directory containing the hydrodynamic data.
- **Define Parameters**:
  - Specifies the square number (`Nsq`), tail type (`tail`), scenario
    (`sce`), outlets, and variable (`var`).
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
  - Saves the results to a CSV file.

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

The script outputs a CSV file containing the trend thresholds for the
specified scenarios and tail types.

### 2. Non-stationnary EVA (TSEVA) <a id="tseva"></a>

#### Script: *02_TSEVA_Run.R*

#### Script Description

This script is designed to run one square of the HERA (Hydrological
Ensemble Prediction System) domain on a High-Performance Computing (HPC)
cluster. It performs extreme trend assessment for either drought or
flood hazards using the TSEVA (Time Series Extreme Value Analysis)
method.

#### Script Overview

1.  **Set Working Directory and Load Functions**:
    - Sets the working directory to the location of the script.
    - Sources the `functions_trends.R` file to load necessary functions.
2.  **Set Data Directory**:
    - Specifies the directory containing the hydrodynamic data.
3.  **Define Parameters**:
    - Specifies the tail type (`tail`), hazard type (`haz`), variable
      (`var`), outlets, outlet name, season, square number (`Nsq`), and
      scenario (`sce`).
    - Loads metadata about the spatial division of the domain into
      chunks.
4.  **Load River Pixel Locations**:
    - Loads river pixel locations and identifiers using the `outletopen`
      function.
5.  **Scenario Differentiation**:
    - Differentiates between scenarios (Historical, Socioeconomic
      counterfactual, Water demand counterfactual, Reservoir+Water
      demand counterfactual).
6.  **Load Discharge Data**:
    - Loads discharge data for the specified scenario and square.
    - Prepares the discharge data for trend analysis, including time
      stamps and series extraction.
7.  **Load Frost Data (for drought analysis)**:
    - Loads frost data for drought analysis and prepares it for further
      processing.
8.  **Threshold Identification**:
    - Reads threshold files for historical and socioeconomic
      counterfactual scenarios.
    - Retains thresholds from the historical run unless they are NA.
9.  **Parallel Processing**:
    - Iterates over the river pixels to compute return periods and
      levels using the TSEVA method.
    - Saves the results, including parameters, return levels, and peaks.

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

The script outputs a list containing parameters, return levels, return
periods, peaks, and catchment rest data for further analysis.

### 3. Changing Hydrological Extremes - univariate analysis <a id="che-u"></a>

#### Script: *03_CHEX_Plot_univariateF.R*

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
    - Sets the working directory to the location of the script.
    - Sources the `functions_trends.R` file to load necessary functions.
2.  **Set Data Directory**:
    - Specifies the directory containing the hydrodynamic data.
3.  **Load Outlets Data**:
    - Reads and processes outlet data for all squares in the HERA
      domain.
    - Loads river pixel locations and identifiers.
4.  **Load Spatial Data for Catchments**:
    - Loads metadata about the domain’s spatial division into chunks.
    - Loads spatial data for catchments, including Hybas07 and European
      Biogeo regions.
5.  **Load Fitting Results**:
    - Loads fitting results from various runs (Historical, Socioeconomic
      Counterfactual, Water Demand Counterfactual, Reservoir+Water
      Demand Counterfactual) for all river pixels in the domain.
6.  **Data Cleaning**:
    - Cleans the data by identifying and handling missing values,
      especially for the last years’ return levels.
7.  **Load IRES Status**:
    - Loads the Intermittent River Status (IRES) for drought analysis.
    - Combines IRES data from different runs to detect any river pixel
      that is intermittent in at least one run.
8.  **Drought Return Level Corrections**:
    - Applies corrections to drought return levels, including reversing
      values and setting negative levels to zero.
9.  **Remove Irrealistic Shape Parameters**:
    - Identifies and removes pixels with unrealistic shape parameters to
      ensure data quality.
10. **Compute Large-Scale Errors**:
    - Computes large-scale errors in return levels for different years
      and scenarios.
11. **Generate Plots**:
    - Generates various plots, including maps of return level errors,
      changes in hazard intensity, and shape parameter instability.
12. **Change Attribution**:
    - Computes changes in return levels attributed to different drivers
      (climate, land use, reservoirs, water demand).
    - Aggregates changes by regions and performs spatial smoothing to
      remove noise from unstable Generalized Pareto Distribution (GPD)
      fits.
13. **Save Outputs**:
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

### 4. Changing Hydrological Extremes - bivariate analysis <a id="che-b"></a>

#### Script: *04_CHEX_Plot_bivariateF.R*

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
    - Sets the working directory and sources necessary functions.
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

### 5. Functions <a id="functions"></a>

#### 5.1 Miscellaneous Functions

##### 5.1.1 Functions for TSEVA

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

##### 5.1.2 Data Extraction Functions

- **`UpAopen`**: Opens upstream area files for selected pixels and
  returns a data frame with upstream area information.

- **`ReservoirOpen`**: Opens reservoir location files and returns a data
  frame with reservoir location information.

- **`disNcopen`**: Extracts discharge time series data for selected
  locations from netCDF files.

- **`disNcopenloc`**: Extracts discharge time series data for a specific
  location from netCDF files.

##### 5.1.3 Miscellaneous Calculations

- **`RPcalc`**: Calculates return periods associated with a return level
  using Generalized Extreme Value (GEV) and Generalized Pareto
  Distribution (GPD) parameters.

- **`RPchangeCal`**: Computes changes in return periods associated with
  a return level for different years and distribution laws.

- **`interid`**: Detects intermittent rivers by identifying periods of
  low or zero discharge and classifies them based on specified criteria.

#### 5.2 Univariate Analysis

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

#### 5.3 Bivariate Analysis

- **`calculateTrendSig`**: Calculates trend significance using the
  Mann-Kendall test and categorizes changes based on significance
  levels.

- **`get_density`**: Calculates the density of points using kernel
  density estimation for specified grid points.
