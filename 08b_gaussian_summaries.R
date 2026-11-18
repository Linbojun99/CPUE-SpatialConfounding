# ==============================================================================
# 08b_gaussian_summaries.R
#
# Summaries of the Gaussian simulation (Section 3.6, Fig. 9, Tables S5 and
# S9): coefficient and index performance of every run, the calibration rule
# that selects the basis dimension k, the paired G2-minus-G1t change error
# with its first-order approximation, the theoretical bias, and the pooled
# G0 coverage and false positive rates.
#
# Input (in each run folder output/gaussian/<setting>/k<k>/, written by
# 08a_gaussian_simulation.R):
#   coefficients.csv, performance.csv, diagnostics.csv, approximation.csv
# Output:
#   output/gaussian/<setting>/k<k>/coefficient_summary.csv
#   output/gaussian/<setting>/k<k>/performance_summary.csv
#   output/gaussian/calibration_selection.csv   calibration of k (Table S5)
#   output/gaussian/all_summaries.csv           values of Tables S5 and S9
#   output/gaussian/paired_replicates.csv       G2-minus-G1t by replicate
#   output/gaussian/paired_summary.csv          means of the paired values
#   output/gaussian/diagnostic_summary.csv      theoretical bias and
#                                               spatial+ residual variance
#   output/gaussian/g0_pooled.csv               pooled G0 rows of Table S5
# ==============================================================================

source("00_settings.R")


# 1. Runs and estimators -------------------------------------------------------

# The eight runs of 08a_gaussian_simulation.R (setting and basis dimension k)
# and the estimators in the order of the tables.
runs <- data.table(
  setting = c("calibration", "calibration", "calibration",
              "reporting", "reporting", "reporting",
              "x1_annual", "time_invariant"),
  k       = c(25, 50, 100, 25, 50, 100, 50, 50)
)
runs[, folder := file.path(dir_gaussian, setting, paste0("k", k))]

gaussian_estimators <- c("GLM", "SGLMM", "RSR", "spatial+")
beta1_true          <- 0.5    # true coefficient of X1


# 2. Coefficient and index summaries of each run -------------------------------

for (run in seq_len(nrow(runs))) {

  coefficients <- fread(file.path(runs$folder[run], "coefficients.csv"))
  performance  <- fread(file.path(runs$folder[run], "performance.csv"))

  # Coefficients: mean error, empirical SD, mean SE and coverage of the 95%
  # intervals and, for a coefficient with a true value of zero, the false
  # positive rate of the two-sided 5% Wald test. Columns ending in
  # "_conditional" use the covariance without smoothing-parameter
  # uncertainty.
  coefficient_summary <- coefficients[, .(
    n                    = .N,
    truth                = truth[1],
    mean_estimate        = mean(estimate),
    mean_error           = mean(error),
    mean_relative_error  = mean(relative_error),
    empirical_sd         = sd(estimate),
    mean_se              = mean(se),
    mean_se_conditional  = mean(se_conditional),
    coverage_conditional = mean(covered_conditional),
    false_positive_rate_conditional = mean(rejected_zero_conditional),
    coverage             = mean(covered),
    false_positive_rate  = mean(rejected_zero)
  ), by = .(scenario, model, term)]
  coefficient_summary[truth != 0, false_positive_rate_conditional := NA_real_]
  coefficient_summary[truth != 0, false_positive_rate := NA_real_]
  coefficient_summary <- coefficient_summary[
    order(scenario, match(model, gaussian_estimators), term)
  ]

  # Index: mean, median and 10% and 90% quantiles of the log-index RMSE
  # (Eq. 20), mean change error (Eq. 22) and mean RMSE of the predicted log
  # density.
  performance_summary <- performance[, .(
    n                        = .N,
    mean_centered_log_rmse   = mean(centered_log_rmse),
    median_centered_log_rmse = median(centered_log_rmse),
    rmse_p10                 = quantile(centered_log_rmse, 0.1),
    rmse_p90                 = quantile(centered_log_rmse, 0.9),
    mean_first_last_error    = mean(first_last_error),
    mean_grid_log_rmse       = mean(grid_log_rmse)
  ), by = .(scenario, model)]
  performance_summary <- performance_summary[
    order(scenario, match(model, gaussian_estimators))
  ]

  fwrite(coefficient_summary,
         file.path(runs$folder[run], "coefficient_summary.csv"))
  fwrite(performance_summary,
         file.path(runs$folder[run], "performance_summary.csv"))
}


# 3. Calibration of the basis dimension k --------------------------------------

# G0 replicates independent of the reporting replicates, for k = 25, 50 and
# 100: relative bias of the spatial+ coefficient of X1, coverage of the
# spatial model for the coefficient of X1, and the false positive rates of
# both for X2 (reported, not used in the selection).
calibration <- data.table(
  k                  = runs[setting == "calibration", k],
  plus_relative_bias = NA_real_,
  spatial_coverage   = NA_real_,
  spatial_X2_fpr     = NA_real_,
  plus_X2_fpr        = NA_real_
)
for (i in seq_len(nrow(calibration))) {
  coefficient_summary <- fread(file.path(
    dir_gaussian, "calibration", paste0("k", calibration$k[i]),
    "coefficient_summary.csv"
  ))
  stopifnot(coefficient_summary$n == n_replicates)
  calibration$plus_relative_bias[i] <- coefficient_summary[
    model == "spatial+" & term == "X1", mean_relative_error]
  calibration$spatial_coverage[i] <- coefficient_summary[
    model == "SGLMM" & term == "X1", coverage]
  calibration$spatial_X2_fpr[i] <- coefficient_summary[
    model == "SGLMM" & term == "X2", false_positive_rate]
  calibration$plus_X2_fpr[i] <- coefficient_summary[
    model == "spatial+" & term == "X2", false_positive_rate]
}

# Selection rule. Candidates with an absolute spatial+ relative bias of at
# most 5% pass the bias criterion. Among them the one whose spatial-model
# coverage is closest to 95% is selected, the smaller k in case of a tie. If
# no candidate passes, the one with the smallest absolute spatial+ bias is
# selected, ties broken by coverage closest to 95% and then by the smaller k.
calibration[, passes_bias := abs(plus_relative_bias) <= 0.05]
calibration[, coverage_distance := abs(spatial_coverage - 0.95)]
if (any(calibration$passes_bias)) {
  candidates <- calibration[passes_bias == TRUE][order(coverage_distance, k)]
} else {
  candidates <- calibration[
    order(abs(plus_relative_bias), coverage_distance, k)
  ]
}
selected_k <- candidates$k[1]
calibration[, selected := k == selected_k]
print(calibration)

# The sensitivity runs of 08a_gaussian_simulation.R were made with k = 50.
stopifnot(selected_k == 50)

setcolorder(calibration, c("k", "plus_relative_bias", "spatial_coverage",
                           "passes_bias", "coverage_distance",
                           "spatial_X2_fpr", "plus_X2_fpr", "selected"))
fwrite(calibration, file.path(dir_gaussian, "calibration_selection.csv"))


# 4. Summary by run, scenario and estimator (Tables S5 and S9) -----------------

# For the reporting and sensitivity runs:
#   bias      mean relative error of the coefficient of X1 (%)
#   sd        empirical SD of the estimates of the coefficient of X1
#   se        mean SE of the coefficient of X1
#   coverage  coverage of the 95% interval for the coefficient of X1 (%)
#   fpr       false positive rate for the coefficient of X2 (%)
#   rmse      median log-index RMSE
#   change    mean change error
summary_rows <- list()
for (run in which(runs$setting != "calibration")) {

  coefficient_summary <- fread(
    file.path(runs$folder[run], "coefficient_summary.csv")
  )
  performance_summary <- fread(
    file.path(runs$folder[run], "performance_summary.csv")
  )

  x1_summary <- coefficient_summary[term == "X1", .(
    scenario, model, n,
    bias     = 100 * mean_relative_error,
    sd       = empirical_sd,
    se       = mean_se,
    coverage = 100 * coverage
  )]
  x2_summary <- coefficient_summary[term == "X2", .(
    scenario, model,
    fpr = 100 * false_positive_rate
  )]
  index_summary <- performance_summary[, .(
    scenario, model,
    rmse   = median_centered_log_rmse,
    change = mean_first_last_error
  )]

  run_summary <- merge(x1_summary, x2_summary, by = c("scenario", "model"))
  run_summary <- merge(run_summary, index_summary,
                       by = c("scenario", "model"))
  run_summary <- run_summary[
    order(scenario, match(model, gaussian_estimators))
  ]
  run_summary[, setting := runs$setting[run]]
  run_summary[, k := runs$k[run]]
  summary_rows[[length(summary_rows) + 1]] <- run_summary
}
all_summaries <- rbindlist(summary_rows)
setcolorder(all_summaries, c("setting", "k"))
fwrite(all_summaries, file.path(dir_gaussian, "all_summaries.csv"))

# Table S9: the reporting run with the selected k.
print(all_summaries[setting == "reporting" & k == selected_k])


# 5. Paired G2-minus-G1t change error and its first-order approximation --------

# G1t and G2 share data and fits. Prediction keeps X1 in G1t and sets it to
# 0 in G2. Within each replicate of the reporting runs:
#   paired_change_error  change error in G2 minus change error in G1t
#   nominal_first_order  -(estimate - true coefficient of X1) * trend of X1
#                        (Section 2.5), the approximation shown in Fig. 9d
#   weighted_first_order the same with the change in the density-weighted
#                        mean of X1 in place of the nominal trend
#   remainder            paired_change_error - nominal_first_order
# For RSR, rsr_minus_source is its change error in G2 minus that of the
# spatial model and rsr_shift_first_order its first-order approximation.
paired_rows <- list()
for (run in which(runs$setting == "reporting")) {

  performance   <- fread(file.path(runs$folder[run], "performance.csv"))
  approximation <- fread(file.path(runs$folder[run], "approximation.csv"))

  paired <- merge(
    performance[scenario == "G2",
                .(model, replicate, first_last_error_G2 = first_last_error)],
    performance[scenario == "G1t",
                .(model, replicate, first_last_error_G1t = first_last_error)],
    by = c("model", "replicate")
  )
  paired <- merge(
    paired,
    approximation[, .(model, replicate, nominal_first_order,
                      weighted_first_order, rsr_minus_source,
                      rsr_shift_first_order)],
    by = c("model", "replicate")
  )
  stopifnot(nrow(paired) == nrow(approximation))

  paired[, paired_change_error := first_last_error_G2 - first_last_error_G1t]
  paired[, remainder := paired_change_error - nominal_first_order]
  paired[, k := runs$k[run]]
  paired_rows[[length(paired_rows) + 1]] <- paired
}
paired_replicates <- rbindlist(paired_rows)
setcolorder(paired_replicates, c("k", "model", "replicate"))
fwrite(paired_replicates, file.path(dir_gaussian, "paired_replicates.csv"))

# Means over replicates and the root mean square of the remainder.
paired_summary <- paired_replicates[, .(
  paired           = mean(paired_change_error),
  nominal          = mean(nominal_first_order),
  weighted         = mean(weighted_first_order),
  remainder        = mean(remainder),
  rms              = sqrt(mean(remainder^2)),
  rsr_minus_source = mean(rsr_minus_source),
  rsr_first_order  = mean(rsr_shift_first_order)
), by = .(k, model)]
paired_summary <- paired_summary[
  order(k, match(model, gaussian_estimators))
]
fwrite(paired_summary, file.path(dir_gaussian, "paired_summary.csv"))
print(paired_summary[k == selected_k])


# 6. Theoretical bias and spatial+ first-stage residual variance ---------------

# For the reporting and sensitivity runs, means over replicates of
#   projection  the theoretical bias (%): the coefficient of X1 in the
#               linear projection of the spatial effect on the year effects,
#               X1 and X2 over the grid, divided by the true coefficient of
#               X1 (Eq. 18, dashed line in Fig. 9b)
#   retained    the share of the variance of X1 at the sampled locations
#               that remains in the spatial+ first-stage residuals
diagnostic_rows <- list()
for (run in which(runs$setting != "calibration")) {

  diagnostics <- fread(file.path(runs$folder[run], "diagnostics.csv"))

  diagnostic_summary <- diagnostics[, .(
    projection = 100 * mean(grid_omitted_field_projection_X1) / beta1_true,
    retained   = mean(X1_first_stage_variance_retained)
  ), by = scenario]
  diagnostic_summary[, setting := runs$setting[run]]
  diagnostic_summary[, k := runs$k[run]]
  diagnostic_rows[[length(diagnostic_rows) + 1]] <- diagnostic_summary
}
diagnostic_summary <- rbindlist(diagnostic_rows)
setcolorder(diagnostic_summary, c("setting", "k", "scenario"))
fwrite(diagnostic_summary, file.path(dir_gaussian, "diagnostic_summary.csv"))
print(diagnostic_summary)


# 7. Pooled G0 coverage and false positive rates (Table S5) --------------------

# The 50 calibration and 50 reporting replicates of G0 at the selected k,
# separately and pooled: number of replicates, number of 95% intervals that
# cover the true coefficient of X1 and number of rejections of a zero
# coefficient of X2. The calibration replicates were used to select k, so
# the pooled values are descriptive.
g0_calibration <- fread(file.path(
  dir_gaussian, "calibration", paste0("k", selected_k), "coefficients.csv"
))[scenario == "G0"]
g0_reporting <- fread(file.path(
  dir_gaussian, "reporting", paste0("k", selected_k), "coefficients.csv"
))[scenario == "G0"]
g0_calibration[, batch := "calibration"]
g0_reporting[, batch := "reporting"]
g0_both   <- rbind(g0_calibration, g0_reporting)
g0_pooled <- copy(g0_both)
g0_pooled[, batch := "pooled"]
g0_both   <- rbind(g0_both, g0_pooled)

g0_summary <- g0_both[, .(
  n                 = sum(term == "X1"),
  covered_X1        = sum(covered[term == "X1"]),
  false_positive_X2 = sum(rejected_zero[term == "X2"])
), by = .(model, batch)]
g0_summary <- g0_summary[
  order(match(model, gaussian_estimators),
        match(batch, c("calibration", "reporting", "pooled")))
]
fwrite(g0_summary, file.path(dir_gaussian, "g0_pooled.csv"))
print(g0_summary)
