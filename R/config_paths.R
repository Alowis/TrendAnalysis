# =============================================================================
#  config_paths.R  --  single source of truth for all data / output paths
# =============================================================================
#
#  Sourced at the top of every analysis script (after functions_trends.R).
#
#  All data lives under the canonical data directory:
#     <PROJECT_ROOT>/ChangingHydroExtremes/data
#
#  The project root can be overridden with an environment variable so the
#  same code runs unchanged on a local machine and on the HPC:
#
#     Local (default) : D:/tilloal/Documents/LFRuns_utils
#     HPC             : export TREND_PROJECT_ROOT=/scratch/.../LFRuns_utils
#                       (or set Sys.setenv(TREND_PROJECT_ROOT=...) before source)
#
#  Paths are built with file.path() so there are never leading-slash /
#  double-slash ambiguities.  Do NOT prefix sub-paths with "/".
# =============================================================================

# --- Project root (env-var override for HPC / portability) -------------------
project_root <- Sys.getenv(
    "TREND_PROJECT_ROOT",
    unset = "D:/tilloal/Documents/LFRuns_utils"
)

# --- Canonical directories ---------------------------------------------------
hydroDir <- file.path(project_root, "ChangingHydroExtremes", "data") # all data
plotDir <- file.path(project_root, "ChangingHydroExtremes", "plots") # all plots

# --- Common sub-directories (derived, canonical) -----------------------------
geoDir <- file.path(hydroDir, "GeoData")
threshDir <- file.path(hydroDir, "Thresholds")
trendVarDir <- file.path(hydroDir, "TrendVar")
droughtDir <- file.path(hydroDir, "Drought")
floodDir <- file.path(hydroDir, "Flood")
riverDir <- file.path(hydroDir, "RiverData")
resDir <- file.path(hydroDir, "reservoirs")

# --- Domain ID scheme --------------------------------------------------------
# ONE definition used across every script that builds pixel identifiers.
# Pixel id = Nsq * ID_MULT + local_index
ID_MULT <- 100000

# --- Ensure output directories exist -----------------------------------------
for (d in c(plotDir, threshDir, trendVarDir)) {
    if (!dir.exists(d)) dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

message("config_paths.R loaded. Data root: ", hydroDir)
