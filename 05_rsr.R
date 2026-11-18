# ==============================================================================
# 05_rsr.R
#
# Restricted spatial regression (RSR) for the fishery simulation. RSR is not
# fitted. It is obtained from each converged fit of ST (X1, X2) by projecting
# the fitted random fields at the simulated sets onto the two covariate
# columns and adding the projected coefficients to the fitted coefficients
# (Eq. 17, Section 2.4). The script calculates the RSR coefficients, their
# standard errors by the delta method from the joint precision of the fitted
# fixed and random effects, 95% Wald intervals, the RSR abundance index over
# the prediction domain (Section 2.6) and the performance measures used for
# the other estimators. No model is optimized again.
#
# Input:
#   data/prediction_domain_cells.csv             cells of the simulation grid,
#                                                in_domain marks the
#                                                prediction domain
#   output/meshes/mesh_estimation.rds            estimation mesh
#   output/simulation/<scenario>/rep<NNNN>.rds   simulated data
#   output/fits/<scenario>/M4_rep<NNNN>.rds      fits of ST (X1, X2)
# Output:
#   output/rsr/<scenario>/rsr_rep<NNNN>.rds   results for one source fit
#   output/rsr/run_status.csv                 one row per source fit
#   output/rsr/rsr_coefficients.csv           coefficients, SEs, intervals
#   output/rsr/rsr_years.csv                  annual log indices
#   output/rsr/rsr_fit.csv                    index and density errors
#   output/performance/rsr_coefficient_summary.csv   by scenario and
#                                                    coefficient
#   output/performance/rsr_index_summary.csv         by scenario and index
#                                                    calculation
#   output/performance/rsr_e2_change_error.csv       change error of RSR and
#                                                    of its source model in
#                                                    E2, by replicate
#
# Run time: about one minute per source fit on one core, about 16 hours for
# all 950 fits. A source fit whose result file exists is skipped, so the
# script can be run in several sessions with different replicates.
# ==============================================================================

source("00_settings.R")
library(sdmTMB)
library(TMB)


# 1. Settings for this run -----------------------------------------------------

# Scenarios and replicates to process. The paper uses all 19 scenarios (the
# 17 of Table 2 and B3 and B4 of Table S2) and replicates 1 to 50.
run_scenarios  <- scenarios$scenario
run_replicates <- 1:n_replicates


# 2. Inputs shared by all source fits ------------------------------------------

# Cells of the prediction domain (1,484 of the 1,592 grid cells).
domain_table <- fread(file.path(dir_data, "prediction_domain_cells.csv"))
domain_cells <- domain_table[in_domain == TRUE, cell]
stopifnot(length(domain_cells) == 1484)

# Mesh of the estimation models.
mesh_estimation <- readRDS(file.path(dir_meshes, "mesh_estimation.rds"))
n_vertices      <- mesh_estimation$n

# Normal quantile of the 95% Wald intervals and the two-sided 5% tests.
z_975 <- qnorm(0.975)

# Largest change in the random effects allowed when the joint precision is
# evaluated at the stored fixed-effect estimates.
mode_tolerance <- 1e-4

# RSR is stored as the seventh estimator. The index is calculated twice:
# "plugin" is the uncorrected index (Eq. 19), and "inherited_M4_bias_ratio"
# multiplies it by the source model's ratio of bias-corrected to uncorrected
# index in each year (Section 2.6). The paper reports the second and uses the
# first as a sensitivity comparison (Supplementary material).
rsr_model_id  <- "M7"
index_methods <- c("plugin", "inherited_M4_bias_ratio")


# 3. RSR for each source fit ---------------------------------------------------

for (scenario_id in run_scenarios) {

  dir.create(file.path(dir_rsr, scenario_id), showWarnings = FALSE)
  covariate_role <- scenarios[scenario == scenario_id, role]

  for (replicate_id in run_replicates) {

    file_result <- file.path(dir_rsr, scenario_id,
                             sprintf("rsr_rep%04d.rds", replicate_id))
    if (file.exists(file_result)) {
      next
    }

    # 3.1 Source fit of ST (X1, X2) and the simulated data it was fitted to.
    source_fit <- readRDS(file.path(dir_fits, scenario_id,
                                    sprintf("M4_rep%04d.rds", replicate_id)))
    simulated  <- readRDS(file.path(dir_simulation, scenario_id,
                                    sprintf("rep%04d.rds", replicate_id)))

    # Only converged source fits are transformed (Section 2.4).
    if (!isTRUE(source_fit$ok) || !isTRUE(source_fit$converged_numeric)) {
      status <- data.table(scenario_id = scenario_id, rep = replicate_id,
                           status = "source fit not converged")
      saveRDS(list(status = status), file_result)
      cat(scenario_id, "replicate", replicate_id,
          ": source fit not converged\n")
      next
    }

    # Stored estimates: all fixed and random effects (par_all) and the fixed
    # effects alone (par_fixed). A source fit that was refitted without a
    # field whose SD was at the boundary has no parameters for that field.
    par_all     <- source_fit$par_all
    par_fixed   <- source_fit$par_fixed
    has_omega   <- "omega_s" %in% names(par_all)
    has_epsilon <- "epsilon_st" %in% names(par_all)

    # 3.2 Rebuild the model object of the source fit without optimizing it,
    # and set its parameters to the stored estimates. The object is needed
    # only for the joint precision matrix (step 3.5).
    sets <- as.data.frame(simulated$sets)
    sets$log_hooks <- log(sets$hooks / 1000)
    sets$year_f    <- factor(sets$year, levels = years)
    sets$vessel    <- factor(sets$vessel)
    sets$hbf_cl    <- factor(sets$hbf_cl, levels = c("16-17", "le15", "ge18"))

    model <- sdmTMB(
      bet_n ~ 0 + year_f + hbf_cl + (1 | vessel) + X1 + X2,
      data           = sets,
      mesh           = make_mesh(sets, c("x", "y"), mesh = mesh_estimation),
      time           = "year",
      offset         = "log_hooks",
      family         = nbinom2(),
      spatial        = if (has_omega) "on" else "off",
      spatiotemporal = if (has_epsilon) "iid" else "off",
      share_range    = TRUE,
      control        = sdmTMBcontrol(multiphase = FALSE, getsd = FALSE),
      silent         = TRUE,
      do_fit         = FALSE
    )
    tmb_object <- model$tmb_obj
    stopifnot(identical(names(tmb_object$env$last.par.best), names(par_all)))
    tmb_object$env$last.par.best <- par_all
    tmb_object$env$last.par      <- par_all

    # Positions of the coefficients of X1 and X2 and of the field values at
    # the mesh vertices in the parameter vector. The fixed effects are ordered
    # as the columns of the design matrix, and epsilon is stored year by
    # year.
    design_columns   <- colnames(model.matrix(~ 0 + year_f + hbf_cl + X1 + X2,
                                              data = sets))
    position_beta    <- which(names(par_all) == "b_j")[
      match(c("X1", "X2"), design_columns)]
    position_omega   <- which(names(par_all) == "omega_s")
    position_epsilon <- which(names(par_all) == "epsilon_st")

    # 3.3 Covariate matrix X at the simulated sets and the projection matrix
    # (X'X)^-1 X' of Eq. 17 (2 rows, one column per set).
    X          <- as.matrix(sets[, c("X1", "X2")])
    projection <- solve(crossprod(X), t(X))

    # Sum of the fitted fields at the sets, u = A omega + A epsilon_t, where A
    # interpolates from the mesh vertices to the set positions. Because u is
    # linear in the field values, the Jacobian of the RSR coefficients
    # beta + (X'X)^-1 X' u with respect to all parameters has ones at the two
    # coefficients and (X'X)^-1 X' A at the field values.
    set_positions <- as.matrix(sets[, c("x", "y")])
    A_sets        <- fm_basis(mesh_estimation, loc = set_positions)
    year_index    <- match(sets$year, years)

    field_at_sets <- numeric(nrow(sets))
    jacobian      <- matrix(0, nrow = 2, ncol = length(par_all))
    jacobian[cbind(1:2, position_beta)] <- 1

    if (has_omega) {
      field_at_sets <- field_at_sets +
        as.numeric(A_sets %*% par_all[position_omega])
      jacobian[, position_omega] <- as.matrix(projection %*% A_sets)
    }
    if (has_epsilon) {
      for (i_year in seq_len(n_years)) {
        in_year       <- year_index == i_year
        position_year <- position_epsilon[
          (i_year - 1) * n_vertices + seq_len(n_vertices)]
        A_year        <- A_sets[in_year, , drop = FALSE]
        field_at_sets[in_year] <- field_at_sets[in_year] +
          as.numeric(A_year %*% par_all[position_year])
        jacobian[, position_year] <- as.matrix(
          projection[, in_year, drop = FALSE] %*% A_year)
      }
    }

    # 3.4 RSR transformation (Eq. 17). The year, HBF and vessel effects are
    # not adjusted.
    beta_hat   <- source_fit$coef$estimate[
      match(c("X1", "X2"), source_fit$coef$term)]
    delta_beta <- as.numeric(projection %*% field_at_sets)
    beta_rsr   <- beta_hat + delta_beta
    field_rsr  <- field_at_sets - as.numeric(X %*% delta_beta)

    # The transformation leaves the fitted predictor at the sets unchanged
    # and makes the adjusted field orthogonal to the covariates. Both are
    # recorded (Supplementary material).
    max_fitted_discrepancy <- max(abs(
      as.numeric(X %*% beta_hat) + field_at_sets -
        as.numeric(X %*% beta_rsr) - field_rsr))
    max_projection_residual <- max(abs(crossprod(X, field_rsr))) / nrow(X)

    # 3.5 Standard errors by the delta method. The joint precision matrix Q
    # of all fixed and random effects is evaluated at the stored estimates,
    # and the covariance of the RSR coefficients is G Q^-1 G', with G the
    # Jacobian of step 3.3. If Q is not available the SEs are missing.
    se_rsr                       <- c(NA_real_, NA_real_)
    se_status                    <- "joint_delta_approximation"
    conditional_mode_discrepancy <- NA_real_

    sd_report <- tryCatch(
      sdreport(tmb_object, par.fixed = par_fixed, getJointPrecision = TRUE,
               getReportCovariance = FALSE, skip.delta.method = TRUE),
      error = function(e) e
    )
    if (inherits(sd_report, "error")) {
      se_status <- paste0("unavailable: ", conditionMessage(sd_report))
    } else if (!isTRUE(sd_report$pdHess) ||
               is.null(sd_report$jointPrecision) ||
               !identical(rownames(sd_report$jointPrecision),
                          names(par_all))) {
      se_status <- "unavailable: joint precision not available"
    } else {
      # The random effects must not move away from the stored estimates.
      conditional_mode_discrepancy <- max(abs(
        sd_report$par.random - par_all[tmb_object$env$random]))
      if (!is.finite(conditional_mode_discrepancy) ||
          conditional_mode_discrepancy > mode_tolerance) {
        se_status <- "unavailable: random effects moved"
      } else {
        covariance_rsr <- as.matrix(
          jacobian %*% solve(sd_report$jointPrecision, t(jacobian)))
        se_rsr <- sqrt(diag(covariance_rsr))
      }
    }

    # 3.6 Coefficient table with 95% Wald intervals. The true coefficient is
    # the effect of the covariate on the expected catch (density plus
    # catchability). "shift" is the coefficient adjustment of Eq. 17.
    beta_true <- unname(source_fit$beta_true[c("X1", "X2")])

    coefficients <- data.table(
      scenario_id = scenario_id,
      model_id    = rsr_model_id,
      rep         = replicate_id,
      term        = c("X1", "X2"),
      estimate    = beta_rsr,
      std.error   = se_rsr,
      truth       = beta_true,
      shift       = delta_beta,
      rel_err     = ifelse(beta_true != 0,
                           (beta_rsr - beta_true) / beta_true, NA_real_),
      conf.low    = beta_rsr - z_975 * se_rsr,
      conf.high   = beta_rsr + z_975 * se_rsr,
      reject_zero = ifelse(is.finite(se_rsr) & se_rsr > 0,
                           abs(beta_rsr / se_rsr) > z_975, NA),
      se_status   = se_status
    )
    coefficients[, covered := conf.low <= truth & truth <= conf.high]

    # 3.7 Predictor on the prediction grid. eta is the linear predictor of
    # the source fit in each cell and year, calculated with the reference
    # covariate values: the simulated values for density covariates and zero
    # for catchability covariates (Section 2.6). The adjusted field uses the
    # simulated covariate values, u_RSR = u - X delta_beta, so the RSR
    # predictor differs from eta by (X_ref - X) delta_beta. The difference is
    # zero for density covariates.
    grid <- merge(
      source_fit$grid_pred,
      simulated$grid[, .(year, cell, X1, X2, area_km2, log_d, sampled)],
      by = c("year", "cell")
    )
    stopifnot(nrow(grid) == nrow(source_fit$grid_pred))

    X_grid <- as.matrix(grid[, .(X1, X2)])
    if (covariate_role == "catchability") {
      X_reference <- X_grid * 0
    } else {
      X_reference <- X_grid
    }
    grid[, eta_rsr := eta + as.numeric((X_reference - X_grid) %*% delta_beta)]

    # 3.8 Annual indices. The uncorrected index is the area-weighted sum of
    # the exponentiated predictor (Eq. 19). The source model's ratio of
    # bias-corrected to uncorrected index is calculated over all 1,592 cells.
    # The indices of RSR and of the source model and the true index are sums
    # over the prediction domain.
    index_all <- grid[, .(log_source_all = log(sum(area_km2 * exp(eta)))),
                      keyby = year]
    source_index   <- source_fit$index
    bias_corrected <- source_index$est[match(index_all$year, source_index$year)]
    index_all[, log_ratio := log(bias_corrected) - log_source_all]

    grid_domain  <- grid[cell %in% domain_cells]
    stopifnot(uniqueN(grid_domain$cell) == length(domain_cells))
    index_domain <- grid_domain[, .(
      log_rsr    = log(sum(area_km2 * exp(eta_rsr))),
      log_source = log(sum(area_km2 * exp(eta))),
      log_true   = log(sum(area_km2 * exp(log_d)))
    ), keyby = year]
    stopifnot(identical(index_all$year, index_domain$year))

    # 3.9 Area-weighted RMSE of log density over the prediction domain, in
    # all cell-years, in those with simulated sets and in those without.
    grid_domain[, squared_error := (eta_rsr - log_d)^2]
    grid_rmse <- grid_domain[
      , sqrt(sum(area_km2 * squared_error) / sum(area_km2))]
    grid_rmse_sampled <- grid_domain[
      sampled == TRUE, sqrt(sum(area_km2 * squared_error) / sum(area_km2))]
    grid_rmse_unsampled <- grid_domain[
      sampled == FALSE, sqrt(sum(area_km2 * squared_error) / sum(area_km2))]

    # 3.10 Index errors for the two index calculations. Log indices are
    # centred on their means over years before they are compared with the
    # centred true log index (Eq. 20). The change error is the error in the
    # change between the first and last year (Eq. 22). Columns ending in
    # "_st" refer to the source model, ST (X1, X2).
    years_by_method <- list()
    fit_by_method   <- list()

    for (index_method in index_methods) {

      if (index_method == "plugin") {
        log_correction <- rep(0, n_years)
      } else {
        log_correction <- index_all$log_ratio
      }
      log_est    <- index_domain$log_rsr + log_correction
      log_est_st <- index_domain$log_source + log_correction
      log_true   <- index_domain$log_true

      err_cen    <- (log_est - mean(log_est)) - (log_true - mean(log_true))
      err_cen_st <- (log_est_st - mean(log_est_st)) -
        (log_true - mean(log_true))

      years_by_method[[index_method]] <- data.table(
        scenario_id  = scenario_id,
        model_id     = rsr_model_id,
        rep          = replicate_id,
        year         = index_domain$year,
        index_method = index_method,
        log_est      = log_est,
        log_est_st   = log_est_st,
        log_true     = log_true,
        err_cen      = err_cen,
        err_cen_st   = err_cen_st
      )

      fit_by_method[[index_method]] <- data.table(
        scenario_id                  = scenario_id,
        model_id                     = rsr_model_id,
        rep                          = replicate_id,
        index_method                 = index_method,
        rmse_cen                     = sqrt(mean(err_cen^2)),
        rmse_cen_st                  = sqrt(mean(err_cen_st^2)),
        first_last_err               = err_cen[n_years] - err_cen[1],
        first_last_err_st            = err_cen_st[n_years] - err_cen_st[1],
        max_log_index_difference     = max(abs(log_est - log_est_st)),
        max_fitted_discrepancy       = max_fitted_discrepancy,
        max_projection_residual      = max_projection_residual,
        grid_rmse                    = grid_rmse,
        grid_rmse_sampled            = grid_rmse_sampled,
        grid_rmse_unsampled          = grid_rmse_unsampled,
        coef_eval                    = all(is.finite(se_rsr) & se_rsr > 0),
        pred_eval                    = all(is.finite(log_est)),
        se_status                    = se_status,
        conditional_mode_discrepancy = conditional_mode_discrepancy,
        s_X                          = simulated$params$s_X
      )
    }

    # 3.11 Save the results of this source fit. Only summaries are stored.
    status <- data.table(scenario_id = scenario_id, rep = replicate_id,
                         status = "computed")
    saveRDS(
      list(status       = status,
           coefficients = coefficients,
           years        = rbindlist(years_by_method),
           fit          = rbindlist(fit_by_method)),
      file_result
    )
    cat(scenario_id, "replicate", replicate_id, ": computed\n")

    # Release the model object before the next source fit.
    rm(model, tmb_object, sd_report)
    gc()
  }
}


# 4. Collect the results of all source fits ------------------------------------

status_list      <- list()
coefficient_list <- list()
years_list       <- list()
fit_list         <- list()

for (scenario_id in run_scenarios) {
  for (replicate_id in run_replicates) {
    result <- readRDS(file.path(dir_rsr, scenario_id,
                                sprintf("rsr_rep%04d.rds", replicate_id)))
    i_result <- length(status_list) + 1
    status_list[[i_result]] <- result$status
    if (result$status$status == "computed") {
      coefficient_list[[i_result]] <- result$coefficients
      years_list[[i_result]]       <- result$years
      fit_list[[i_result]]         <- result$fit
    }
  }
}

run_status       <- rbindlist(status_list)
rsr_coefficients <- rbindlist(coefficient_list)
rsr_years        <- rbindlist(years_list)
rsr_fit          <- rbindlist(fit_list)

fwrite(run_status,       file.path(dir_rsr, "run_status.csv"))
fwrite(rsr_coefficients, file.path(dir_rsr, "rsr_coefficients.csv"))
fwrite(rsr_years,        file.path(dir_rsr, "rsr_years.csv"))
fwrite(rsr_fit,          file.path(dir_rsr, "rsr_fit.csv"))


# 5. Summaries over replicates -------------------------------------------------

# Coefficients: mean error, mean relative error (for true values other than
# zero), empirical SD, mean SE, coverage of the 95% intervals and, for true
# values of zero, the false positive rate (Section 2.6). Coverage and false
# positive rates use the estimates with a finite SE.
rsr_coefficients[, err := estimate - truth]
rsr_coefficients[, valid_se := is.finite(std.error) & std.error > 0 &
                   se_status == "joint_delta_approximation"]

coefficient_summary <- rsr_coefficients[, .(
  n              = .N,
  n_se           = sum(valid_se),
  truth          = truth[1],
  mean_estimate  = mean(estimate),
  mean_error     = mean(err),
  relative_bias  = if (truth[1] != 0) mean(rel_err) else NA_real_,
  q10_relative   = if (truth[1] != 0) quantile(rel_err, 0.1) else NA_real_,
  q90_relative   = if (truth[1] != 0) quantile(rel_err, 0.9) else NA_real_,
  empirical_sd   = sd(estimate),
  mean_se        = if (any(valid_se)) mean(std.error[valid_se]) else NA_real_,
  coverage       = if (any(valid_se)) mean(covered[valid_se]) else NA_real_,
  false_positive = if (truth[1] == 0 && any(valid_se)) {
    mean(reject_zero[valid_se])
  } else {
    NA_real_
  }
), by = .(scenario_id, term)]

# Indices: log-index RMSE (Eq. 20), change error (Eq. 22), the difference in
# change error from the source model, and the largest values of the three
# checks recorded for each source fit.
index_summary <- rsr_fit[, .(
  n                               = .N,
  median_rmse                     = median(rmse_cen),
  q25_rmse                        = quantile(rmse_cen, 0.25),
  q75_rmse                        = quantile(rmse_cen, 0.75),
  mean_rmse                       = mean(rmse_cen),
  mean_change_error               = mean(first_last_err),
  q10_change_error                = quantile(first_last_err, 0.1),
  q90_change_error                = quantile(first_last_err, 0.9),
  sd_change_error                 = sd(first_last_err),
  mean_change_difference_M4       = mean(first_last_err - first_last_err_st),
  max_log_difference_M4           = max(max_log_index_difference),
  max_fitted_difference           = max(max_fitted_discrepancy),
  max_projection_residual         = max(max_projection_residual),
  max_conditional_mode_difference = max(conditional_mode_discrepancy,
                                        na.rm = TRUE)
), by = .(scenario_id, index_method)]

fwrite(coefficient_summary,
       file.path(dir_performance, "rsr_coefficient_summary.csv"))
fwrite(index_summary,
       file.path(dir_performance, "rsr_index_summary.csv"))


# 6. Change error of RSR relative to its source model in E2 --------------------

# In E2 the catchability covariate X1 has a trend, and the change error of RSR
# differs from that of its source model. To first order the difference is
# -delta_beta1 * s_X, where delta_beta1 is the adjustment of the X1
# coefficient (Eq. 17) and s_X, the SD of X1 without the trend, is the change
# in X1 added between the first and last year (Eq. 8, Section 3.4 and
# Supplementary material). The difference does not depend on the index
# calculation, because both estimators use the same ratio.
if ("E2" %in% run_scenarios) {

  e2_change_error <- merge(
    rsr_fit[scenario_id == "E2",
            .(rep, index_method, first_last_err, first_last_err_st, s_X)],
    rsr_coefficients[scenario_id == "E2" & term == "X1",
                     .(rep, d1 = shift)],
    by = "rep"
  )
  e2_change_error[, paired_change := first_last_err - first_last_err_st]
  e2_change_error[, trend_only    := -d1 * s_X]
  e2_change_error[, residual      := paired_change - trend_only]

  fwrite(e2_change_error[, .(rep, index_method, d1, s_X, first_last_err,
                             first_last_err_st, paired_change, trend_only,
                             residual)],
         file.path(dir_performance, "rsr_e2_change_error.csv"))

  # Means over replicates for each index calculation.
  print(e2_change_error[, .(
    mean_change_error_rsr    = mean(first_last_err),
    mean_change_error_source = mean(first_last_err_st),
    mean_paired_change       = mean(paired_change),
    mean_delta_beta1         = mean(d1),
    mean_s_X                 = mean(s_X),
    mean_first_order         = mean(trend_only),
    mean_residual            = mean(residual),
    min_residual             = min(residual),
    max_residual             = max(residual)
  ), by = index_method])
}
