# ==============================================================================
# 06_performance.R
#
# Calculates the performance measures of Section 2.6 for the six models
# fitted to the fishery simulation (04_fit_estimators.R). Indices, true
# indices and predicted densities are restricted to the prediction domain of
# 1,484 cells. For each fit the script calculates
#   - for the coefficients of X1 and X2: the error, the relative error,
#     whether the 95% confidence interval covers the true value and whether
#     a two-sided Wald test rejects a zero coefficient at the 5% level (a
#     false positive when the true value is zero),
#   - for the abundance index: the log-index RMSE (Eq. 20), the slope lambda
#     (Eq. 21) and the change error (Eq. 22),
#   - for the predicted log density: the area-weighted mean error and RMSE
#     over the cells of the prediction domain.
# The measures are then summarized by scenario and model (sections 5 and 6).
# Sections 7 to 10 add the measures of RSR, which 05_rsr.R calculated from
# the fits of ST (X1, X2), and compare the seven estimators within
# replicates: log-index RMSE against ST (none) (Fig. 6), the change error in
# E2 against its first-order approximation (Section 3.5) and the scenario
# pairs of Table S2. Run 05_rsr.R before this script.
#
# Input:
#   output/fits/<scenario>/<model>_rep<NNNN>.rds   fits of the six models
#   data/prediction_domain_cells.csv               cells of the simulation
#                                                  grid, in_domain marks the
#                                                  prediction domain
#   output/rsr/rsr_fit.csv                         index measures of RSR
#   output/rsr/rsr_coefficients.csv                coefficients of RSR
#   output/simulation/E2/rep<NNNN>.rds             s_X, the change added to
#                                                  X1 in E2 (section 9)
# Output (in output/performance/):
#   performance_fits.csv                 one row per fit: convergence checks
#                                        and the measures of the index and
#                                        of the predicted density
#   performance_coefficients.csv         one row per fit and coefficient
#   performance_index_years.csv          one row per fit and year: error of
#                                        the log index
#   performance_summary.csv              by scenario and model
#   performance_coefficient_summary.csv  by scenario, model and coefficient
#   performance_paired.csv               mean within-replicate differences
#                                        in RMSE from ST (none), with RSR
#   performance_e2_first_order.csv       change error in E2 minus that in E1
#                                        or E3 and its first-order
#                                        approximation, with RSR
#   performance_scenario_pairs.csv       B3 minus A3 and B4 minus A1
#                                        (Table S2), with RSR
#
# Run time: mostly the reading of the fit files (5,700 files, about 7 GB,
# for the 19 scenarios), from several minutes to more than an hour depending
# on the speed of file access.
# ==============================================================================

source("00_settings.R")


# 1. Settings for this run -----------------------------------------------------

# Scenarios, replicates and models to summarize. The defaults are those of
# the paper: the 17 scenarios of Table 2 and the two of Table S2 (B3, B4),
# replicates 1 to 50 and the six models.
run_scenarios  <- scenarios$scenario
run_replicates <- 1:n_replicates
run_models     <- estimators$model


# 2. Prediction domain and fixed settings --------------------------------------

# Cells of the prediction domain (1,484 of the 1,592 cells of the grid).
domain_table <- fread(file.path(dir_data, "prediction_domain_cells.csv"))
domain_cells <- domain_table[in_domain == TRUE, cell]

# The slope lambda (Eq. 21) is not calculated when the SD of the true
# centred log index over years is below this value.
slope_min_sd <- 0.02

# Model from which the differences in RMSE of section 8 are taken:
# ST (none).
reference_model <- "M2"

# Name, in the tables written by 05_rsr.R, of the RSR index used in the
# paper (the uncorrected index multiplied by the ratio of bias-corrected to
# uncorrected index of the source model, Section 2.6). The other index of
# those tables, "plugin", is the uncorrected index. RSR has the model
# identifier M7 in those tables.
rsr_index_method <- "inherited_M4_bias_ratio"

# Models that do not estimate the coefficient of X1: ST (none) and ST (X2).
models_without_x1 <- c("M2", "M6")

# Scenario pairs of Table S2: B3 and B4 have the fields, covariates and sets
# of A3 and A1 and add a density effect of X2.
scenario_pairs <- data.table(scenario_b = c("B3", "B4"),
                             scenario_a = c("A3", "A1"))


# 3. Measures of each fit ------------------------------------------------------

rows_fits         <- list()
rows_coefficients <- list()
rows_years        <- list()
n_files_missing   <- 0

for (scenario_id in run_scenarios) {
  for (model_id in run_models) {
    for (replicate_id in run_replicates) {

      file_fit <- file.path(
        dir_fits, scenario_id,
        sprintf("%s_rep%04d.rds", model_id, replicate_id)
      )
      if (!file.exists(file_fit)) {
        n_files_missing <- n_files_missing + 1
        next
      }
      fit      <- readRDS(file_fit)
      fit_name <- paste(scenario_id, model_id, replicate_id)

      # Scenario, model, replicate and whether the fit ran without error.
      error_message <- NA
      if (!is.null(fit$error)) {
        error_message <- fit$error
      }
      fit_id <- data.table(
        scenario_id = fit$scenario_id,
        model_id    = fit$model_id,
        rep         = fit$rep,
        ok          = isTRUE(fit$ok),
        stage       = fit$stage,
        error       = error_message
      )
      if (!isTRUE(fit$ok)) {
        rows_fits[[fit_name]] <- fit_id
        next
      }

      # 3.1 Index and true index over the prediction domain. The
      # bias-corrected index was calculated over all 1,592 cells. It is
      # multiplied by the ratio of the uncorrected sum over the prediction
      # domain (Eq. 19) to that over all cells in the same year (Section
      # 2.6). The true index is the sum of area times true density over the
      # prediction domain (Eq. 9).
      grid <- merge(fit$grid_pred, fit$grid_truth, by = c("year", "cell"))
      grid[, in_domain := cell %in% domain_cells]
      index <- grid[, .(
        sum_all    = sum(area_km2 * exp(eta)),
        sum_domain = sum(area_km2 * exp(eta) * in_domain),
        true_index = sum(area_km2[in_domain] * exp(log_d[in_domain]))
      ), by = year]
      index <- merge(fit$index[, .(year, est_all = est)], index, by = "year")
      setorder(index, year)
      index[, est      := sum_domain * est_all / sum_all]
      index[, log_est  := log(est)]
      index[, log_true := log(true_index)]

      # 3.2 Measures of the index. Both series are centred on their means
      # over years (l_hat and l_true). err_cen is their difference and
      # err_unc the error of the log index without centring.
      index[, l_hat   := log_est - mean(log_est)]
      index[, l_true  := log_true - mean(log_true)]
      index[, err_unc := log_est - log_true]
      index[, err_cen := l_hat - l_true]

      # Log-index RMSE (Eq. 20).
      rmse_cen <- sqrt(mean(index$err_cen^2))

      # Slope lambda (Eq. 21). The two series are centred, so the slope of
      # the regression with an intercept equals that of Eq. 21.
      slope <- NA_real_
      if (sd(index$l_true) > slope_min_sd) {
        slope <- unname(coef(lm(l_hat ~ l_true, data = index))[2])
      }

      # Change of the true log index between the first and the last year
      # and the change error (Eq. 22). change_err_plugin is the change
      # error of the uncorrected index (Eq. 19).
      true_change <- index$log_true[n_years] - index$log_true[1]
      change_err  <- (index$log_est[n_years] - index$log_est[1]) -
        true_change
      change_err_plugin <-
        (log(index$sum_domain[n_years]) - log(index$sum_domain[1])) -
        true_change

      rows_years[[fit_name]] <- cbind(
        fit_id[, .(scenario_id, model_id, rep)],
        index[, .(year, err_unc, err_cen)]
      )

      # 3.3 Error of the predicted log density in the cells of the
      # prediction domain: area-weighted mean error (bias) and RMSE over
      # all cells and years, over those with sets in the year (sampled)
      # and over those without sets (unsampled).
      domain_grid <- grid[in_domain == TRUE]
      domain_grid[, err := eta - log_d]
      density_error <- c(
        grid_bias           = NA_real_, grid_rmse           = NA_real_,
        grid_bias_sampled   = NA_real_, grid_rmse_sampled   = NA_real_,
        grid_bias_unsampled = NA_real_, grid_rmse_unsampled = NA_real_
      )
      for (cells_used in c("", "_sampled", "_unsampled")) {
        rows <- rep(TRUE, nrow(domain_grid))
        if (cells_used == "_sampled") {
          rows <- domain_grid$sampled
        }
        if (cells_used == "_unsampled") {
          rows <- !domain_grid$sampled
        }
        if (any(rows)) {
          area <- domain_grid$area_km2[rows]
          err  <- domain_grid$err[rows]
          density_error[[paste0("grid_bias", cells_used)]] <-
            sum(area * err) / sum(area)
          density_error[[paste0("grid_rmse", cells_used)]] <-
            sqrt(sum(area * err^2) / sum(area))
        }
      }

      # 3.4 Measures of the coefficients of the covariates in the model.
      if (nrow(fit$coef) > 0) {
        coefficients_fit <- copy(fit$coef)
        coefficients_fit[, truth := fit$beta_true[term]]
        coefficients_fit[, err := estimate - truth]
        coefficients_fit[, covered := conf.low <= truth & truth <= conf.high]
        coefficients_fit[
          , reject_zero := abs(estimate / std.error) > qnorm(0.975)
        ]
        coefficients_fit[
          , rel_err := ifelse(truth != 0, (estimate - truth) / truth, NA_real_)
        ]
        rows_coefficients[[fit_name]] <- cbind(
          fit_id[, .(scenario_id, model_id, rep)],
          coefficients_fit
        )
      }

      # 3.5 One row for the fit. Coefficients are evaluated when the fit
      # converged and all standard errors are finite (coef_eval), the index
      # and the predicted density when the fit converged and the index is
      # finite (pred_eval).
      converged <- isTRUE(fit$converged_numeric)
      coef_eval <- converged &&
        (nrow(fit$coef) == 0 || all(is.finite(fit$coef$std.error)))
      pred_eval <- converged && all(is.finite(index$log_est))
      refit_spec <- ""
      if (!is.null(fit$refit_spec)) {
        refit_spec <- paste(fit$refit_spec, collapse = "/")
      }
      checks <- fit$diag$sanity

      rows_fits[[fit_name]] <- cbind(fit_id, data.table(
        converged_numeric   = fit$converged_numeric,
        boundary_sigma      = fit$boundary_sigma,
        boundary_refit      = isTRUE(fit$boundary_refit),
        refit_spec          = refit_spec,
        sanity_all_ok       = checks[["all_ok"]],
        hessian_ok          = checks[["hessian_ok"]],
        eigen_values_ok     = checks[["eigen_values_ok"]],
        nlminb_ok           = checks[["nlminb_ok"]],
        gradients_ok        = checks[["gradients_ok"]],
        se_na_ok            = checks[["se_na_ok"]],
        se_magnitude_ok     = checks[["se_magnitude_ok"]],
        sigmas_ok           = checks[["sigmas_ok"]],
        range_ok            = checks[["range_ok"]],
        max_abs_gradient    = fit$diag$max_abs_gradient,
        sigma_O             = fit$diag$sigma_O,
        sigma_E             = fit$diag$sigma_E,
        range               = fit$diag$range,
        coef_eval           = coef_eval,
        pred_eval           = pred_eval,
        rmse_cen            = rmse_cen,
        mean_err_unc        = mean(index$err_unc),
        slope               = slope,
        true_change         = true_change,
        change_err          = change_err,
        change_err_plugin   = change_err_plugin,
        grid_bias           = density_error[["grid_bias"]],
        grid_rmse           = density_error[["grid_rmse"]],
        grid_bias_sampled   = density_error[["grid_bias_sampled"]],
        grid_rmse_sampled   = density_error[["grid_rmse_sampled"]],
        grid_bias_unsampled = density_error[["grid_bias_unsampled"]],
        grid_rmse_unsampled = density_error[["grid_rmse_unsampled"]],
        minutes             = fit$minutes
      ))
    }
  }
  cat("Scenario", scenario_id, "finished\n")
}
cat("Fits read:", length(rows_fits), " files not found:", n_files_missing,
    "\n")


# 4. Tables with one row per fit, coefficient or year --------------------------

fits              <- rbindlist(rows_fits, fill = TRUE)
coefficient_table <- rbindlist(rows_coefficients, fill = TRUE)
index_years       <- rbindlist(rows_years, fill = TRUE)

fwrite(fits, file.path(dir_performance, "performance_fits.csv"))
fwrite(coefficient_table,
       file.path(dir_performance, "performance_coefficients.csv"))
fwrite(index_years,
       file.path(dir_performance, "performance_index_years.csv"))


# 5. Summary by scenario and model ---------------------------------------------

# Rates of fits that did not converge, of fits with a field SD at the
# boundary and of refits, then means, medians and quartiles of the index and
# density measures over the fits that converged (pred_eval).
performance_summary <- fits[, .(
  n_fits                 = .N,
  numeric_fail_rate      = mean(!ok | !converged_numeric %in% TRUE),
  boundary_rate          = mean(boundary_sigma %in% TRUE),
  boundary_refit_rate    = mean(boundary_refit %in% TRUE),
  boundary_refit_ok_rate = mean(boundary_refit %in% TRUE &
                                  converged_numeric %in% TRUE),
  n_pred                 = sum(pred_eval %in% TRUE),
  rmse_cen               = mean(rmse_cen[pred_eval %in% TRUE]),
  mean_err_unc           = mean(mean_err_unc[pred_eval %in% TRUE]),
  grid_rmse              = mean(grid_rmse[pred_eval %in% TRUE]),
  rmse_cen_median        = median(rmse_cen[pred_eval %in% TRUE]),
  rmse_cen_q25           = quantile(rmse_cen[pred_eval %in% TRUE], 0.25),
  rmse_cen_q75           = quantile(rmse_cen[pred_eval %in% TRUE], 0.75),
  slope_median           = median(slope[pred_eval %in% TRUE], na.rm = TRUE),
  change_err_mean        = mean(change_err[pred_eval %in% TRUE]),
  change_err_median      = median(change_err[pred_eval %in% TRUE]),
  change_err_sd          = sd(change_err[pred_eval %in% TRUE]),
  grid_rmse_sampled_median   =
    median(grid_rmse_sampled[pred_eval %in% TRUE]),
  grid_rmse_unsampled_median =
    median(grid_rmse_unsampled[pred_eval %in% TRUE])
), by = .(scenario_id, model_id)]

fwrite(performance_summary,
       file.path(dir_performance, "performance_summary.csv"))


# 6. Summary of the coefficients -----------------------------------------------

# Coefficients of the fits that converged with finite standard errors
# (coef_eval).
fits_evaluated <- fits[coef_eval %in% TRUE, .(scenario_id, model_id, rep)]
coefficients_evaluated <- coefficient_table[
  fits_evaluated, on = .(scenario_id, model_id, rep), nomatch = NULL
]

# Mean error, SD of the estimates across replicates (the empirical SD, as
# the true value is the same in all replicates), root mean square error, mean
# of the estimated standard errors, coverage of the 95% confidence
# intervals and the share of fits in which a zero coefficient was rejected.
# For a coefficient with a true value of zero, reject_zero is the false
# positive rate. The mean relative error is the mean error divided by the
# true value.
coefficient_summary <- coefficients_evaluated[, .(
  n           = .N,
  truth       = truth[1],
  mean_err    = mean(err),
  emp_sd      = sd(err),
  rmse        = sqrt(mean(err^2)),
  mean_se     = mean(std.error),
  cover95     = mean(covered),
  reject_zero = mean(reject_zero)
), by = .(scenario_id, model_id, term)]
coefficient_summary[
  , mean_rel_err := ifelse(truth != 0, mean_err / truth, NA_real_)
]

fwrite(coefficient_summary,
       file.path(dir_performance, "performance_coefficient_summary.csv"))


# 7. Measures of all seven estimators ------------------------------------------

# RSR is not fitted. Its coefficients and indices were calculated by
# 05_rsr.R from the converged fits of ST (X1, X2) and are read here.
rsr_fit <- fread(file.path(dir_rsr, "rsr_fit.csv"))
rsr_fit <- rsr_fit[scenario_id %in% run_scenarios & rep %in% run_replicates]
rsr_coefficients <- fread(file.path(dir_rsr, "rsr_coefficients.csv"))
rsr_coefficients <- rsr_coefficients[
  scenario_id %in% run_scenarios & rep %in% run_replicates
]

# One row per RSR result: the measures of the RSR index of the paper and
# the change error of the uncorrected index.
rsr_measures <- merge(
  rsr_fit[
    index_method == rsr_index_method,
    .(scenario_id, model_id, rep, coef_eval, pred_eval, rmse_cen, grid_rmse,
      change_err = first_last_err)
  ],
  rsr_fit[
    index_method == "plugin",
    .(scenario_id, rep, change_err_plugin = first_last_err)
  ],
  by = c("scenario_id", "rep")
)

# Index measures of the six models (fits that converged, pred_eval) and of
# RSR.
index_measures <- rbind(
  fits[
    pred_eval %in% TRUE,
    .(scenario_id, model_id, rep, rmse_cen, grid_rmse, change_err,
      change_err_plugin)
  ],
  rsr_measures[
    pred_eval == TRUE,
    .(scenario_id, model_id, rep, rmse_cen, grid_rmse, change_err,
      change_err_plugin)
  ]
)

# Coefficient estimates of the six models (coef_eval) and of RSR.
rsr_coefficients <- merge(
  rsr_coefficients,
  rsr_measures[coef_eval == TRUE, .(scenario_id, rep)],
  by = c("scenario_id", "rep")
)
coefficient_measures <- rbind(
  coefficients_evaluated[
    , .(scenario_id, model_id, rep, term, estimate, truth, err, rel_err)
  ],
  rsr_coefficients[
    , .(scenario_id, model_id, rep, term, estimate, truth,
        err = estimate - truth, rel_err)
  ]
)


# 8. Log-index RMSE paired with ST (none) --------------------------------------

# Models are paired within replicates in which both fits converged (Section
# 2.6). For each scenario and each of the other five models and RSR (a), the
# table gives the mean of the within-replicate differences from ST (none)
# (b) in log-index RMSE (Fig. 6, Section 3.4) and in the RMSE of the
# predicted log density.
reference_measures <- index_measures[
  model_id == reference_model,
  .(scenario_id, rep, rmse_reference = rmse_cen, grid_reference = grid_rmse)
]
both <- merge(index_measures[model_id != reference_model],
              reference_measures, by = c("scenario_id", "rep"))
paired <- both[, .(
  b              = reference_model,
  n              = .N,
  diff_rmse_cen  = mean(rmse_cen - rmse_reference),
  diff_grid_rmse = mean(grid_rmse - grid_reference)
), by = .(scenario_id, a = model_id)]
paired <- paired[order(match(scenario_id, scenarios$scenario), a)]

fwrite(paired, file.path(dir_performance, "performance_paired.csv"))


# 9. Change error in E2 and its first-order approximation ----------------------

# E1 and E2 have the same sets and fits. In E2 X1 is a catchability
# covariate with a trend and the index is predicted with X1 = 0. To first
# order the change error (Eq. 22) then differs from that in E1 by
# -(estimate of beta1 - beta1) * s_X for an estimator with X1 and by
# +beta1 * s_X for one without X1, where s_X is the change added to X1
# between the first and the last year (Sections 2.5 and 3.5). The table
# gives, for each estimator, the mean within-replicate difference in change
# error (E2 minus E1, and E2 minus E3, which has no trend in X1) for the
# index of the paper and for the uncorrected index, and the mean of the
# first-order approximation.
if (all(c("E1", "E2") %in% run_scenarios)) {

  # s_X of each replicate of E2.
  trend_e2 <- data.table(rep = run_replicates, s_X = NA_real_)
  for (i in seq_along(run_replicates)) {
    simulated <- readRDS(file.path(
      dir_simulation, "E2", sprintf("rep%04d.rds", run_replicates[i])
    ))
    trend_e2$s_X[i] <- simulated$params$s_X
  }

  # Change errors in E2, the estimate of beta1 (0 for an estimator without
  # X1) and the first-order approximation.
  beta1_true <- coefficient_measures[scenario_id == "E2" & term == "X1",
                                     truth][1]
  change_e2 <- merge(
    index_measures[scenario_id == "E2",
                   .(model_id, rep, change_err, change_err_plugin)],
    coefficient_measures[scenario_id == "E2" & term == "X1",
                         .(model_id, rep, estimate)],
    by = c("model_id", "rep"), all.x = TRUE
  )
  change_e2[model_id %in% models_without_x1, estimate := 0]
  change_e2 <- merge(change_e2, trend_e2, by = "rep")
  change_e2[, first_order := -(estimate - beta1_true) * s_X]

  rows_contrasts <- list()
  for (scenario_reference in intersect(c("E1", "E3"), run_scenarios)) {
    both <- merge(
      change_e2,
      index_measures[scenario_id == scenario_reference,
                     .(model_id, rep, change_err_reference = change_err,
                       change_err_plugin_reference = change_err_plugin)],
      by = c("model_id", "rep")
    )
    both[, difference := change_err - change_err_reference]
    both[, difference_plugin :=
           change_err_plugin - change_err_plugin_reference]
    rows_contrasts[[scenario_reference]] <- both[, .(
      contrast               = paste0("E2_minus_", scenario_reference),
      n                      = .N,
      mean_first_order       = mean(first_order),
      mean_difference        = mean(difference),
      sd_difference          = sd(difference),
      mean_difference_plugin = mean(difference_plugin),
      mean_remainder         = mean(difference - first_order),
      min_remainder          = min(difference - first_order),
      max_remainder          = max(difference - first_order)
    ), by = model_id]
  }
  e2_first_order <- rbindlist(rows_contrasts)
  setcolorder(e2_first_order, c("contrast", "model_id"))

  fwrite(e2_first_order,
         file.path(dir_performance, "performance_e2_first_order.csv"))
}


# 10. Scenario pairs of Table S2 -----------------------------------------------

# For each pair (B3 and A3, B4 and A1) and each estimator: the mean relative
# errors of the coefficients of X1 and X2 in the B scenario, and the means
# of the within-replicate differences (B minus A) in the relative error of
# the X1 coefficient, in log-index RMSE and in change error. Relative errors
# are fractions of the true value (Table S2 gives percentages).
rows_pairs <- list()
for (i in seq_len(nrow(scenario_pairs))) {
  scenario_b <- scenario_pairs$scenario_b[i]
  scenario_a <- scenario_pairs$scenario_a[i]
  if (!all(c(scenario_a, scenario_b) %in% run_scenarios)) {
    next
  }

  # Mean relative errors in the B scenario.
  x1_b <- coefficient_measures[
    scenario_id == scenario_b & term == "X1",
    .(n_coef = .N, x1_rel_err = mean(rel_err)), by = model_id
  ]
  x2_b <- coefficient_measures[
    scenario_id == scenario_b & term == "X2",
    .(x2_rel_err = mean(rel_err)), by = model_id
  ]

  # Within-replicate differences in the relative error of the X1 coefficient.
  both <- merge(
    coefficient_measures[scenario_id == scenario_b & term == "X1",
                         .(model_id, rep, rel_err_b = rel_err)],
    coefficient_measures[scenario_id == scenario_a & term == "X1",
                         .(model_id, rep, rel_err_a = rel_err)],
    by = c("model_id", "rep")
  )
  x1_difference <- both[, .(delta_x1_rel_err = mean(rel_err_b - rel_err_a)),
                        by = model_id]

  # Within-replicate differences in log-index RMSE and in change error.
  both <- merge(
    index_measures[scenario_id == scenario_b,
                   .(model_id, rep, rmse_b = rmse_cen, change_b = change_err)],
    index_measures[scenario_id == scenario_a,
                   .(model_id, rep, rmse_a = rmse_cen, change_a = change_err)],
    by = c("model_id", "rep")
  )
  index_difference <- both[, .(
    n                = .N,
    delta_rmse_cen   = mean(rmse_b - rmse_a),
    delta_change_err = mean(change_b - change_a)
  ), by = model_id]

  pair_table <- merge(index_difference, x1_b, by = "model_id", all.x = TRUE)
  pair_table <- merge(pair_table, x1_difference, by = "model_id",
                      all.x = TRUE)
  pair_table <- merge(pair_table, x2_b, by = "model_id", all.x = TRUE)
  pair_table[, comparison := paste(scenario_b, "minus", scenario_a)]
  rows_pairs[[i]] <- pair_table[, .(
    comparison, model_id, n_coef, x1_rel_err, delta_x1_rel_err, x2_rel_err,
    n, delta_rmse_cen, delta_change_err
  )]
}

# The table is written when at least one pair is among the scenarios of the
# run.
if (length(rows_pairs) > 0) {
  scenario_pair_table <- rbindlist(rows_pairs)
  fwrite(scenario_pair_table,
         file.path(dir_performance, "performance_scenario_pairs.csv"))
}
