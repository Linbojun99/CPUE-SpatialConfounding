# ==============================================================================
# 00_settings.R
#
# Settings shared by all scripts: packages, folders, the study period, the
# scenarios of the fishery simulation (Table 2 and Supplementary Table S2)
# and the estimators (Table 3).
#
# Every script starts with source("00_settings.R"). Run all scripts with the
# working directory set to the folder that contains this file.
# ==============================================================================


# 1. Packages ------------------------------------------------------------------

library(data.table)
library(Matrix)
library(fmesher)


# 2. Folders -------------------------------------------------------------------

# Input files (see data/README.md). The set-level logbook data are
# confidential and are not distributed with the code.
dir_data <- "data"

# Everything the scripts write goes below output/.
dir_output      <- "output"
dir_meshes      <- file.path(dir_output, "meshes")
dir_logbook     <- file.path(dir_output, "logbook_model")
dir_simulation  <- file.path(dir_output, "simulation")
dir_fits        <- file.path(dir_output, "fits")
dir_rsr         <- file.path(dir_output, "rsr")
dir_performance <- file.path(dir_output, "performance")
dir_x1_both     <- file.path(dir_output, "x1_density_and_catchability")
dir_mesh_runs   <- file.path(dir_output, "mesh_comparison")
dir_gaussian    <- file.path(dir_output, "gaussian")
dir_tables      <- file.path(dir_output, "tables")
dir_figures     <- file.path(dir_output, "figures")

for (d in c(dir_meshes, dir_logbook, dir_simulation, dir_fits, dir_rsr,
            dir_performance, dir_x1_both, dir_mesh_runs, dir_gaussian,
            dir_tables, dir_figures)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}


# 3. Study period and projection -----------------------------------------------

years   <- 2015:2023
n_years <- length(years)

# Lambert azimuthal equal-area projection centred on the fishing area (km).
crs_km <- paste(
  "+proj=laea +lat_0=0 +lon_0=-25 +x_0=0 +y_0=0",
  "+datum=WGS84 +units=km +no_defs"
)


# 4. Design of the fishery simulation ------------------------------------------

n_replicates <- 50       # replicates per scenario
n_sets       <- 10000    # sets drawn in each replicate
base_seed    <- 20261014

# True relative abundance: a 30% decline between the first and last year.
abundance_trend <- exp(seq(0, log(0.7), length.out = n_years))

# Scenarios (Table 2). Coefficients are multiples of the base effect size.
#   rho_omega    correlation between X1 and the spatial field omega
#   rho_12       correlation between X1 and X2
#   range_ratio  range of the unconfounded component of X1 relative to the
#                range of omega (c_r)
#   beta*_d      density coefficients of X1 and X2
#   beta*_q      catchability coefficients of X1 and X2
#   x1_trend     trend added to X1 between the first and last year, in units
#                of the SD of X1 without the trend
#   role         whether the covariates enter the index as density covariates
#                or are held at reference values as catchability covariates
scenarios <- data.table(
  scenario    = c("A1", "A2", "A3", "A4", "A5", "A6",
                  "B1", "B2",
                  "C1", "C2",
                  "D1", "D2", "D3", "D4",
                  "E1", "E2", "E3"),
  rho_omega   = c(0, 0, 0.6, 0.6, 0.9, 0.6,
                  0.6, 0,
                  0, 0.6,
                  0.6, 0.6, 0, 0,
                  0.6, 0.6, 0.6),
  rho_12      = c(0, 0.8, 0, 0.8, 0.8, 0.95,
                  0.8, 0.8,
                  0.8, 0.8,
                  0.8, 0.8, 0.8, 0.8,
                  0.8, 0.8, 0.8),
  range_ratio = c(1, 1, 1, 1, 1, 1,
                  1, 1,
                  1, 1,
                  0.3, 3, 0.3, 3,
                  1, 1, 1),
  beta1_d     = c(1, 1, 1, 1, 1, 1,
                  1, 1,
                  0, 0,
                  1, 1, 1, 1,
                  1, 0, 0),
  beta2_d     = c(0, 0, 0, 0, 0, 0,
                  1, 1,
                  0, 0,
                  0, 0, 0, 0,
                  0, 0, 0),
  beta1_q     = c(0, 0, 0, 0, 0, 0,
                  0, 0,
                  0, 0,
                  0, 0, 0, 0,
                  0, 1, 1),
  beta2_q     = 0,
  x1_trend    = c(0, 0, 0, 0, 0, 0,
                  0, 0,
                  0, 0,
                  0, 0, 0, 0,
                  1, 1, 0),
  role        = c(rep("density", 15), "catchability", "catchability")
)

# Two further scenarios, B3 and B4, are reported in Supplementary Table S2
# only. They are simulated, fitted and summarized by the same scripts as the
# 17 scenarios of Table 2 and are marked in the column `supplementary`.
supplementary_scenarios <- data.table(
  scenario    = c("B3", "B4"),
  rho_omega   = c(0.6, 0),
  rho_12      = 0,
  range_ratio = 1,
  beta1_d     = 1,
  beta2_d     = 1,
  beta1_q     = 0,
  beta2_q     = 0,
  x1_trend    = 0,
  role        = "density"
)
scenarios[, supplementary := FALSE]
supplementary_scenarios[, supplementary := TRUE]
scenarios <- rbind(scenarios, supplementary_scenarios)
scenarios[, group := substr(scenario, 1, 1)]

# The 17 scenarios of Table 2, in the order used in tables and figures.
main_scenarios <- scenarios[supplementary == FALSE, scenario]


# 5. Estimators ----------------------------------------------------------------

# Six models fitted to every simulated data set (Table 3). The seventh
# estimator, restricted spatial regression (RSR), is computed afterwards from
# the fits of ST (X1, X2) (05_rsr.R).
estimators <- data.table(
  model        = c("M1", "M2", "M3", "M4", "M5", "M6"),
  label        = c("GLMM (X1, X2)", "ST (none)", "ST (X1)",
                   "ST (X1, X2)", "ST spatial+", "ST (X2)"),
  covariates   = list(c("X1", "X2"), character(0), "X1",
                      c("X1", "X2"), c("X1", "X2"), "X2"),
  spatial      = c(FALSE, TRUE, TRUE, TRUE, TRUE, TRUE),
  spatial_plus = c(FALSE, FALSE, FALSE, FALSE, TRUE, FALSE)
)

# Settings of the estimation step.
k_first_stage   <- 100     # basis dimension of the spatial+ first-stage spline
n_draws_index   <- 300     # joint draws used for index intervals
gradient_limit  <- 0.001   # largest absolute gradient of a converged fit
se_limit        <- 100     # largest standard error of a converged fit
field_sd_limit  <- 0.01    # a field SD below this value is at the boundary
