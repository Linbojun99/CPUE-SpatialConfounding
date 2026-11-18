# ==============================================================================
# 09_tables.R
#
# Assembles the numerical tables of the paper from the outputs of the other
# scripts: Table 1 (values of the operating-model parameters) and
# Supplementary Tables S1, S2 and S4 to S9. Each table is written as a CSV
# file with the rows, the columns, the order and the rounding of the printed
# table. Supplementary Table S3 is written by
# 07_x1_density_and_catchability.R (table_s3_x1_density_and_catchability.csv
# in output/tables/) and is not repeated here.
#
# A cell with one number is written as a number, rounded as in the printed
# table. A cell that combines several numbers (the values of Table 1, the
# pairs of Table S1, the cells of Tables S7 and S8) is written as text in
# the printed form. An empty cell is a "-" of the printed table: a
# coefficient that is not in the model or a quantity that was not calculated.
#
# Input, by table:
#   Table 1    output/logbook_model/logbook_estimates.csv         (script 02)
#              output/simulation/operating_model_parameters.rds   (script 03)
#              output/simulation/<scenario>/rep<NNNN>.rds         (script 03)
#   Table S1   output/meshes/mesh_settings.csv                    (script 01)
#   Table S2   output/performance/performance_scenario_pairs.csv  (script 06)
#   Table S4   output/performance/performance_coefficient_summary.csv    (06)
#              output/performance/rsr_coefficient_summary.csv            (05)
#   Table S5   output/gaussian/all_summaries.csv                 (script 08b)
#              output/gaussian/paired_summary.csv                (script 08b)
#              output/gaussian/calibration_selection.csv         (script 08b)
#              output/gaussian/g0_pooled.csv                     (script 08b)
#   Table S6   output/simulation/<scenario>/rep<NNNN>.rds         (script 03)
#   Table S7   output/performance/performance_coefficient_summary.csv    (06)
#              output/performance/rsr_coefficient_summary.csv            (05)
#   Table S8   output/performance/performance_summary.csv         (script 06)
#              output/performance/rsr_index_summary.csv           (script 05)
#   Table S9   output/gaussian/all_summaries.csv                 (script 08b)
#              output/gaussian/calibration_selection.csv         (script 08b)
# Output (in output/tables/):
#   table_1_operating_model.csv             Table 1, with the two values of
#                                           its note in the last two rows
#   table_s1_mesh_construction.csv          Table S1
#   table_s2_second_density_effect.csv      Table S2
#   table_s4_false_positive_rates.csv       Table S4
#   table_s5_gaussian_sensitivity.csv       Table S5
#   table_s6_simulated_data_checks.csv      Table S6
#   table_s7_coefficient_bias_coverage.csv  Table S7
#   table_s8_index_rmse.csv                 Table S8
#   table_s9_gaussian_simulation.csv        Table S9
#
# Run time: from about ten seconds to a few minutes depending on the speed of
# file access, mostly the reading of the 850 simulated data sets (section 2).
# Sets are read for their catches and numbers of hooks only.
# ==============================================================================

source("00_settings.R")


# 1. Estimators and rows shared by the tables ----------------------------------

# Labels of the seven estimators of the fishery simulation. RSR has the model
# identifier M7 in the tables written by 05_rsr.R and 06_performance.R.
estimator_labels        <- c(estimators$label, "RSR")
names(estimator_labels) <- c(estimators$model, "M7")

# Order of the estimators in the rows of Table S2 and in the columns of
# Table S8. ST (none) has no coefficient and is left out of Tables S4 and S7.
estimator_order        <- c("M1", "M2", "M3", "M4", "M6", "M5", "M7")
estimators_coefficient <- c("M1", "M3", "M4", "M6", "M5", "M7")

# Name, in the tables written by 05_rsr.R, of the RSR index used in the
# paper (Section 2.6).
rsr_index_method <- "inherited_M4_bias_ratio"

# Estimators of the Gaussian simulation: names in the tables written by
# 08b_gaussian_summaries.R and labels in Tables S5 and S9.
gaussian_models <- c("GLM", "SGLMM", "RSR", "spatial+")
gaussian_labels <- c("GLM", "Spatial", "RSR", "spatial+")

# One row per scenario of Table 2 and coefficient. The true value of a
# coefficient is the sum of its density and catchability effects (Eq. 11).
# Table S4 has the rows with a true value of zero, Table S7 the other rows.
design <- scenarios[supplementary == FALSE]
coefficient_rows <- rbind(
  design[, .(scenario, term = "X1", coefficient = "beta_1",
             true_zero = beta1_d + beta1_q == 0)],
  design[, .(scenario, term = "X2", coefficient = "beta_2",
             true_zero = beta2_d + beta2_q == 0)]
)
coefficient_rows <- coefficient_rows[
  order(match(scenario, main_scenarios), term)
]


# 2. Values read from the simulated data sets ----------------------------------

# For each replicate of the 17 scenarios of Table 2: s_X, the SD of X1
# without the trend over all cells and years (Table 1), the realized
# correlations stored with the simulated data (gen_diag), the share of sets
# with a catch of zero and the median catch per 1,000 hooks (Table S6).
# "grid_equal" marks correlations over all cells of the simulation grid and
# all years with equal weights, "sets" correlations at the simulated sets.
rows_simulated <- list()

for (scenario_id in main_scenarios) {
  for (replicate_id in seq_len(n_replicates)) {

    simulated <- readRDS(file.path(dir_simulation, scenario_id,
                                   sprintf("rep%04d.rds", replicate_id)))
    stored     <- as.data.table(simulated$gen_diag)[stat == "cor"]
    catch      <- simulated$sets$bet_n
    catch_rate <- catch / simulated$sets$hooks * 1000

    rows_simulated[[length(rows_simulated) + 1]] <- data.table(
      scenario            = scenario_id,
      rep                 = replicate_id,
      s_x                 = simulated$params$s_X,
      cor_x1sp_omega_grid = stored[
        version == "grid_equal" & var == "X1_spatial-omega", value
      ],
      cor_x1_omega_grid   = stored[
        version == "grid_equal" & var == "X1-omega", value
      ],
      cor_x1_x2_grid      = stored[
        version == "grid_equal" & var == "X1-X2", value
      ],
      cor_x1_omega_sets   = stored[
        version == "sets" & var == "X1-omega", value
      ],
      cor_x1_x2_sets      = stored[
        version == "sets" & var == "X1-X2", value
      ],
      zero_share          = mean(catch == 0),
      median_catch_rate   = median(catch_rate)
    )
  }
  cat("Read scenario", scenario_id, "\n")
}

simulated_values <- rbindlist(rows_simulated)


# 3. Table 1: parameters and settings of the operating model -------------------

# Estimates of the model fitted to the logbook data (Eq. 2) and the
# parameter values of the operating model derived from them.
logbook_table     <- fread(file.path(dir_logbook, "logbook_estimates.csv"))
logbook_estimates <- setNames(logbook_table$value, logbook_table$parameter)
operating_model   <- readRDS(
  file.path(dir_simulation, "operating_model_parameters.rds")
)

# SD of delta, the year-to-year part of X1 (Eq. 6). It is a setting of
# 03_simulate_fishery_data.R (sd_delta) and is not stored in its output.
sd_delta <- 0.5

# s_X is calculated for each replicate. The table gives its mean over all
# replicates of the 17 scenarios.
s_x_mean <- mean(simulated_values$s_x)

# Values of the scenario factors in the 17 scenarios (Table 2).
values_rho_omega <- sort(unique(design$rho_omega))
values_rho_12    <- sort(unique(design$rho_12))
values_c_r       <- sort(unique(design$range_ratio))

# Symbols, descriptions and equation numbers as in the printed table. The
# values are text because several cells hold more than one number. The two
# HBF effects are separated by a comma here and by a semicolon in the
# printed table.
table_1 <- data.table(
  symbol      = c("beta_0", "tau_t", "gamma_h", "r", "sigma_O", "sigma_omega",
                  "beta_1^d, beta_1^q, beta_2^d", "sigma_E", "sigma_V", "phi",
                  "rho_omega", "rho_12", "c_r", "SD(delta)", "s_X"),
  description = c(
    "Intercept (log fish per 1,000 hooks)",
    "Underlying trend in density, 2015-2023",
    "HBF effect (15 or fewer, 18 or more)",
    "Range of omega and epsilon (km)",
    "SD of the spatial field in Eq. 2",
    "SD of omega in the operating model",
    "Non-zero covariate effects",
    "SD of epsilon",
    "SD of vessel effects",
    "Dispersion of the negative binomial",
    "Correlation between X_1^sp and omega",
    "Correlation between X_1^0 and X_2",
    "Range of Z, delta and W relative to r",
    "SD of the year-to-year part of X_1",
    "SD of X_1^0 over all cells and years"
  ),
  equation    = c("3", "3", "4", "2, 3", "2", "3", "3, 4", "2, 3", "2, 4",
                  "1", "5", "7", "5-7", "6", "7, 8"),
  value       = c(
    sprintf("%.2f", operating_model$b0),
    paste(abundance_trend[1], "to",
          sprintf("%.2f", abundance_trend[n_years])),
    paste(sprintf("%.3f", operating_model$gamma_hbf[["le15"]]),
          sprintf("%.3f", operating_model$gamma_hbf[["ge18"]]), sep = ", "),
    formatC(round(operating_model$range), format = "d", big.mark = ","),
    sprintf("%.2f", logbook_estimates[["sigma_O"]]),
    sprintf("%.2f", operating_model$sigma_omega),
    sprintf("%.2f", operating_model$beta1),
    sprintf("%.2f", operating_model$sigma_E),
    sprintf("%.2f", operating_model$sigma_V),
    sprintf("%.2f", operating_model$phi),
    paste(values_rho_omega, collapse = ", "),
    paste(values_rho_12, collapse = ", "),
    paste(values_c_r, collapse = ", "),
    paste(sd_delta),
    sprintf("%.2f", s_x_mean)
  )
)

# The two parts of the intercept given in the note below the table: the mean
# of the estimated year effects and the mean of the fitted omega + epsilon
# at the logbook sets.
mean_year_effect   <- logbook_estimates[["mean of alpha_t"]]
mean_field_at_sets <-
  logbook_estimates[["mean of omega + epsilon at the logbook sets"]]

table_1_note <- data.table(
  symbol      = c("mean of alpha_t", "mean of omega + epsilon"),
  description = c(
    "Note: mean of the estimated year effects in Eq. 2",
    "Note: mean of the fitted omega + epsilon at the logbook sets"
  ),
  equation    = "2",
  value       = sprintf("%.2f", c(mean_year_effect, mean_field_at_sets))
)
table_1 <- rbind(table_1, table_1_note)

fwrite(table_1, file.path(dir_tables, "table_1_operating_model.csv"))
print(table_1)


# 4. Table S1: mesh construction settings --------------------------------------

# Numbers of vertices and settings of the meshes, written by
# 01_build_meshes.R. The mesh used to fit Eq. 2 has no maximum edge lengths
# in that table. It is described in the text above Table S1 and is not a row
# of it.
mesh_settings <- fread(file.path(dir_meshes, "mesh_settings.csv"))
mesh_settings <- mesh_settings[!is.na(max_edge_inner)]

# The values for the interior and the outer region are written as a pair.
table_s1 <- mesh_settings[, .(
  use,
  vertices,
  cutoff,
  max_edge_lengths = paste(
    formatC(max_edge_inner, format = "d", big.mark = ","),
    formatC(max_edge_outer, format = "d", big.mark = ","), sep = ", "
  ),
  extensions       = paste(
    formatC(extension_inner, format = "d", big.mark = ","),
    formatC(extension_outer, format = "d", big.mark = ","), sep = ", "
  )
)]

fwrite(table_s1, file.path(dir_tables, "table_s1_mesh_construction.csv"))
print(table_s1)


# 5. Table S2: adding a second density effect (B3, B4) -------------------------

# Mean relative errors of the coefficients in B3 and B4 and means of the
# within-replicate differences from A3 and A1. Relative errors are stored as
# fractions and written here in percent, their differences in percentage
# points.
scenario_pairs <- fread(
  file.path(dir_performance, "performance_scenario_pairs.csv")
)
scenario_pairs <- scenario_pairs[
  order(match(comparison, unique(comparison)),
        match(model_id, estimator_order))
]

table_s2 <- scenario_pairs[, .(
  pair               = sub(" minus ", " - ", comparison),
  estimator          = unname(estimator_labels[model_id]),
  x1_error_percent   = round(100 * x1_rel_err, 1),
  delta_x1_error_pp  = round(100 * delta_x1_rel_err, 1),
  x2_error_percent   = round(100 * x2_rel_err, 1),
  delta_rmse         = round(delta_rmse_cen, 4),
  delta_change_error = round(delta_change_err, 4)
)]

fwrite(table_s2, file.path(dir_tables, "table_s2_second_density_effect.csv"))
print(table_s2)


# 6. Table S4: false positive rates --------------------------------------------

# Coefficient summaries of the six fitted models and of RSR, in one table:
# mean relative error, coverage of the 95% confidence intervals and share of
# the fits in which a zero coefficient was rejected, which is the false
# positive rate when the true value is zero.
coefficient_summary <- fread(
  file.path(dir_performance, "performance_coefficient_summary.csv")
)
rsr_coefficient_summary <- fread(
  file.path(dir_performance, "rsr_coefficient_summary.csv")
)
coefficient_values <- rbind(
  coefficient_summary[, .(
    scenario = scenario_id, model_id, term,
    relative_error = mean_rel_err,
    coverage       = cover95,
    false_positive = reject_zero
  )],
  rsr_coefficient_summary[, .(
    scenario = scenario_id, model_id = "M7", term,
    relative_error = relative_bias,
    coverage,
    false_positive
  )]
)

# Rows: coefficients with a true value of zero. Columns: the estimators.
# Rates are percentages.
rows_s4  <- coefficient_rows[true_zero == TRUE]
table_s4 <- rows_s4[, .(scenario, coefficient)]

for (model in estimators_coefficient) {
  values  <- coefficient_values[model_id == model]
  matched <- match(paste(rows_s4$scenario, rows_s4$term),
                   paste(values$scenario, values$term))
  set(table_s4, j = estimator_labels[[model]],
      value = round(100 * values$false_positive[matched]))
}

fwrite(table_s4, file.path(dir_tables, "table_s4_false_positive_rates.csv"))
print(table_s4)


# 7. Table S5: Gaussian sensitivity to basis dimension and annual variation ----

# Summaries of the Gaussian simulation by run, scenario and estimator
# (bias, coverage and false positive rates in percent), the mean paired
# G2-minus-G1t change error, the basis dimension k selected by the
# calibration and the pooled G0 counts.
all_summaries  <- fread(file.path(dir_gaussian, "all_summaries.csv"))
paired_summary <- fread(file.path(dir_gaussian, "paired_summary.csv"))
calibration    <- fread(file.path(dir_gaussian, "calibration_selection.csv"))
g0_pooled      <- fread(file.path(dir_gaussian, "g0_pooled.csv"))

selected_k <- calibration[selected == TRUE, k]
other_k    <- calibration[selected == FALSE, k]

# Runs in the table: the reporting runs at the two values of k that were not
# selected, and the two runs at the selected k with an annual component in
# X1 only and in neither covariate. G1t and G2 were simulated in the
# reporting runs only.
sensitivity_runs <- data.table(
  setting = c(rep("reporting", length(other_k)), "x1_annual",
              "time_invariant"),
  k       = c(other_k, selected_k, selected_k),
  label   = c(rep("Both annual", length(other_k)), "X1 annual",
              "Time-invariant")
)

rows_s5 <- list()

for (i in seq_len(nrow(sensitivity_runs))) {
  run_setting <- sensitivity_runs$setting[i]
  run_k       <- sensitivity_runs$k[i]

  for (j in seq_along(gaussian_models)) {
    run_values <- all_summaries[
      setting == run_setting & k == run_k & model == gaussian_models[j]
    ]
    values_g0 <- run_values[scenario == "G0"]
    values_g1 <- run_values[scenario == "G1"]

    paired_difference <- NA_real_
    if (run_setting == "reporting") {
      paired_difference <- paired_summary[
        k == run_k & model == gaussian_models[j], paired
      ]
    }

    rows_s5[[length(rows_s5) + 1]] <- data.table(
      setting                    = paste0(sensitivity_runs$label[i],
                                          ", k = ", run_k),
      estimator                  = gaussian_labels[j],
      g0_bias_percent            = round(values_g0$bias, 1),
      g0_coverage_percent        = round(values_g0$coverage),
      g0_false_positives_percent = round(values_g0$fpr),
      g1_bias_percent            = round(values_g1$bias, 1),
      g1_sd                      = round(values_g1$sd, 3),
      g2_minus_g1t               = round(paired_difference, 4)
    )
  }
}

# Pooled G0 rows: coverage and false positive rates over the calibration and
# reporting replicates at the selected k. The other quantities were not
# pooled.
pooled <- g0_pooled[batch == "pooled"]

for (j in seq_along(gaussian_models)) {
  values_pooled <- pooled[model == gaussian_models[j]]

  rows_s5[[length(rows_s5) + 1]] <- data.table(
    setting                    = paste0("G0 pooled, k = ", selected_k,
                                        ", n = ", values_pooled$n),
    estimator                  = gaussian_labels[j],
    g0_bias_percent            = NA_real_,
    g0_coverage_percent        = round(
      100 * values_pooled$covered_X1 / values_pooled$n
    ),
    g0_false_positives_percent = round(
      100 * values_pooled$false_positive_X2 / values_pooled$n
    ),
    g1_bias_percent            = NA_real_,
    g1_sd                      = NA_real_,
    g2_minus_g1t               = NA_real_
  )
}

table_s5 <- rbindlist(rows_s5)

fwrite(table_s5, file.path(dir_tables, "table_s5_gaussian_sensitivity.csv"))
print(table_s5)


# 8. Table S6: checks of the simulated data ------------------------------------

# Means over the 50 replicates of each scenario of the values read in
# section 2.
table_s6 <- simulated_values[, .(
  cor_x1sp_omega_grid  = round(mean(cor_x1sp_omega_grid), 2),
  cor_x1_omega_grid    = round(mean(cor_x1_omega_grid), 2),
  cor_x1_x2_grid       = round(mean(cor_x1_x2_grid), 2),
  cor_x1_omega_sets    = round(mean(cor_x1_omega_sets), 2),
  cor_x1_x2_sets       = round(mean(cor_x1_x2_sets), 2),
  zero_catches_percent = round(100 * mean(zero_share), 1),
  median_catch_rate    = round(mean(median_catch_rate), 1)
), by = scenario]

fwrite(table_s6, file.path(dir_tables, "table_s6_simulated_data_checks.csv"))
print(table_s6)


# 9. Table S7: bias and interval coverage of estimated covariate effects -------

# Rows: coefficients with a non-zero true value. Columns: the estimators.
# Each cell gives the mean relative error (%) and, in parentheses, the
# coverage of the 95% confidence intervals (%), from the coefficient
# summaries read in section 6.
rows_s7  <- coefficient_rows[true_zero == FALSE]
table_s7 <- rows_s7[, .(scenario, coefficient)]

for (model in estimators_coefficient) {
  values  <- coefficient_values[model_id == model]
  matched <- match(paste(rows_s7$scenario, rows_s7$term),
                   paste(values$scenario, values$term))
  cells   <- sprintf(
    "%d (%d)",
    as.integer(round(100 * values$relative_error[matched])),
    as.integer(round(100 * values$coverage[matched]))
  )
  cells[is.na(matched)] <- NA
  set(table_s7, j = estimator_labels[[model]], value = cells)
}

fwrite(table_s7,
       file.path(dir_tables, "table_s7_coefficient_bias_coverage.csv"))
print(table_s7)


# 10. Table S8: accuracy of the index ------------------------------------------

# Median and 25% and 75% quantiles of the log-index RMSE (Eq. 20) over the
# converged fits of each scenario and model, and the same for the RSR index.
performance_summary <- fread(
  file.path(dir_performance, "performance_summary.csv")
)
rsr_index_summary <- fread(
  file.path(dir_performance, "rsr_index_summary.csv")
)
rmse_values <- rbind(
  performance_summary[, .(
    scenario = scenario_id, model_id,
    rmse_median = rmse_cen_median,
    rmse_q25    = rmse_cen_q25,
    rmse_q75    = rmse_cen_q75
  )],
  rsr_index_summary[index_method == rsr_index_method, .(
    scenario = scenario_id, model_id = "M7",
    rmse_median = median_rmse,
    rmse_q25    = q25_rmse,
    rmse_q75    = q75_rmse
  )]
)

# Rows: the 17 scenarios of Table 2. Columns: the estimators. Each cell
# gives the median and, in parentheses, the two quantiles.
table_s8 <- data.table(scenario = main_scenarios)

for (model in estimator_order) {
  values  <- rmse_values[model_id == model]
  matched <- match(table_s8$scenario, values$scenario)
  set(table_s8, j = estimator_labels[[model]], value = sprintf(
    "%.2f (%.2f-%.2f)", values$rmse_median[matched],
    values$rmse_q25[matched], values$rmse_q75[matched]
  ))
}

fwrite(table_s8, file.path(dir_tables, "table_s8_index_rmse.csv"))
print(table_s8)


# 11. Table S9: Gaussian simulation with annual variation in both covariates ---

# Reporting run at the selected k, from the summaries read in section 7.
# The spatial model, RSR and spatial+ use a thin plate regression spline
# (TPRS) with basis dimension k, which is part of their labels.
gaussian_scenarios <- c("G0", "G1", "G1t", "G2")
reporting_values   <- all_summaries[setting == "reporting" & k == selected_k]
reporting_values   <- reporting_values[
  order(match(scenario, gaussian_scenarios), match(model, gaussian_models))
]

labels_s9 <- gaussian_labels[match(reporting_values$model, gaussian_models)]
with_k    <- reporting_values$model != "GLM"
labels_s9[with_k] <- paste0(labels_s9[with_k], " (TPRS, k = ", selected_k, ")")

table_s9 <- reporting_values[, .(
  scenario,
  estimator                  = labels_s9,
  x1_bias_percent            = round(bias, 1),
  x1_sd                      = round(sd, 3),
  coverage_percent           = round(coverage),
  x2_false_positives_percent = round(fpr),
  index_rmse                 = round(rmse, 3),
  change_error               = round(change, 3)
)]

fwrite(table_s9, file.path(dir_tables, "table_s9_gaussian_simulation.csv"))
print(table_s9)
