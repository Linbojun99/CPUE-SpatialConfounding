# ==============================================================================
# 04_fit_estimators.R
#
# Fits the six estimation models of Table 3 to every simulated data set of
# the fishery simulation: GLMM (X1, X2) (Eq. 11), ST (none) (Eq. 12),
# ST (X1) (Eq. 13), ST (X1, X2) (Eq. 14), ST (X2) (Eq. 15) and ST spatial+
# (Eq. 16). For every fit the script saves the coefficients of X1 and X2 with
# their standard errors, the convergence checks, the refit made when the SD
# of a random field is at the boundary (Section 2.4), the predicted log
# density in the 1-degree cells, and the bias-corrected abundance index with
# joint draws of the log index. The same script fits the models on the
# coarser and the finer mesh for the mesh comparison of Table S1 and
# Fig. S8 (see section 1).
#
# Random numbers: only the joint draws of the log index use them. The seed
# of the draws of one fit is
#   base_seed + 100000 * i + 10 * r + m
# where i is the row of the scenario in `scenarios`, r the replicate and m
# the row of the model in `estimators` (00_settings.R). The seeds are
# distinct for all scenarios, models and replicates up to 9999.
#
# Input:
#   output/simulation/<scenario>/rep<NNNN>.rds   simulated sets, cell values
#                                                and true index of one
#                                                replicate
#   data/simulation_grid.rds                     1,592 cells of the simulation
#                                                grid (columns cell, x, y)
#   output/meshes/mesh_estimation.rds            estimation mesh with 1,676
#                                                vertices (for the mesh
#                                                comparison also the files
#                                                mesh_estimation_coarser.rds
#                                                and mesh_estimation_finer.rds)
# Output:
#   output/fits/<scenario>/<model>_rep<NNNN>.rds   one file per scenario,
#     model (M1 to M6, see 00_settings.R) and replicate. The elements of
#     the saved list are described in step 3.6.
#   output/mesh_comparison/<mesh>/<scenario>/<model>_rep<NNNN>.rds   the
#     same for the fits of the mesh comparison
#
# Run time: on one core less than 1 minute for a GLMM fit, about 8 minutes
# for an ST fit and 15 minutes for an ST spatial+ fit, about 710 core hours
# for the 5,700 fits of the 19 scenarios. To run in parallel, start several
# R sessions with different scenarios in section 1.
# ==============================================================================

source("00_settings.R")
library(sdmTMB)
library(mgcv)


# 1. Settings for this run -----------------------------------------------------

# Scenarios, replicates and models to fit. The defaults are those of the
# paper: the 17 scenarios of Table 2 and the two of Table S2 (B3, B4),
# replicates 1 to 50 and the six models. A fit whose output file exists is
# skipped, so an interrupted run can be started again. To run in parallel,
# start several R sessions and give each session different scenarios here.
# Keep A4 with E3 and E1 with E2 in the same session (see section 2).
run_scenarios  <- scenarios$scenario
run_replicates <- 1:n_replicates
run_models     <- estimators$model

# Estimation mesh and folder of the fits. The defaults are the mesh used in
# the paper (1,676 vertices) and output/fits.
mesh_file    <- file.path(dir_meshes, "mesh_estimation.rds")
dir_fits_run <- dir_fits

# Mesh comparison (Table S1, Fig. S8). It uses the preliminary replicates
# 9001 to 9010 of A2 and D1, which are separate from replicates 1 to 50
# (simulate them with 03_simulate_fishery_data.R first). Run this script
# three times with run_scenarios <- c("A2", "D1") and
#   selected mesh (1,676 vertices): replicates 9001 to 9010, ST (X1) and
#   ST (X1, X2)
#     run_replicates <- 9001:9010
#     run_models     <- c("M3", "M4")
#     mesh_file      <- file.path(dir_meshes, "mesh_estimation.rds")
#     dir_fits_run   <- file.path(dir_mesh_runs, "selected")
#   coarser mesh (620 vertices): the same replicates and models
#     mesh_file      <- file.path(dir_meshes, "mesh_estimation_coarser.rds")
#     dir_fits_run   <- file.path(dir_mesh_runs, "coarser")
#   finer mesh (3,653 vertices): replicates 9001 to 9003, ST (X1, X2)
#     run_replicates <- 9001:9003
#     run_models     <- "M4"
#     mesh_file      <- file.path(dir_meshes, "mesh_estimation_finer.rds")
#     dir_fits_run   <- file.path(dir_mesh_runs, "finer")
# figures/figS08_mesh_sensitivity.R reads the three folders.


# 2. Inputs and fixed settings -------------------------------------------------

# Estimation mesh and the centres of the 1-degree cells (projected, km).
mesh_estimation <- readRDS(mesh_file)
simulation_grid <- readRDS(file.path(dir_data, "simulation_grid.rds"))
cell_centres    <- as.data.table(simulation_grid)[, .(cell, x, y)]

# E3 has the same simulated sets as A4, and E2 the same as E1. Within each
# pair the scenarios differ only in whether X1 enters the index as a density
# covariate or is held at its reference value as a catchability covariate
# (Table 2). When both scenarios of a pair are in the run, each model is
# fitted once to the shared sets and the index is predicted for both
# scenarios.
same_sets_as <- c(E3 = "A4", E2 = "E1")

# Columns of the simulated sets used by the models.
columns_fitted <- c("year", "vessel", "x", "y", "hooks", "hbf_cl",
                    "X1", "X2", "bet_n")

# Classes of hooks between floats (HBF). The first class is the reference
# level at which the index is predicted (Eq. 3).
hbf_levels <- c("16-17", "le15", "ge18")

# Checks of sdmTMB::sanity() that define a converged fit (Section 2.4): the
# optimizer reported convergence, the Hessian is positive definite without
# extreme eigenvalues, and all absolute gradients are below gradient_limit.
convergence_checks <- c("nlminb_ok", "hessian_ok", "eigen_values_ok",
                        "gradients_ok")

# Steps taken to obtain a converged fit of ST (X1, X2) for the weights of
# the spatial+ smooths (step 3.3).
weights_steps <- c("original", "same_model_extra_optimization",
                   "boundary_refit", "boundary_extra_optimization")


# 3. Fits ----------------------------------------------------------------------

for (scenario_id in run_scenarios) {

  # A scenario that shares its sets with another scenario of this run is
  # predicted from the fits to that scenario (section 2).
  if (scenario_id %in% names(same_sets_as) &&
      same_sets_as[[scenario_id]] %in% run_scenarios) {
    next
  }
  scenarios_predicted <- c(
    scenario_id,
    intersect(names(same_sets_as)[same_sets_as == scenario_id], run_scenarios)
  )

  for (replicate_id in run_replicates) {

    # 3.1 Simulated data of every scenario predicted from these fits. The
    # columns used by the models must be identical among the scenarios.
    simulated <- list()
    for (scenario_predicted in scenarios_predicted) {
      simulated[[scenario_predicted]] <- readRDS(file.path(
        dir_simulation, scenario_predicted,
        sprintf("rep%04d.rds", replicate_id)
      ))
      stopifnot(isTRUE(all.equal(
        simulated[[scenario_predicted]]$sets[, ..columns_fitted],
        simulated[[scenario_id]]$sets[, ..columns_fitted],
        tolerance = 0
      )))
    }

    for (model_id in run_models) {

      # Skip the fit when its output files exist.
      files_output <- file.path(
        dir_fits_run, scenarios_predicted,
        sprintf("%s_rep%04d.rds", model_id, replicate_id)
      )
      if (all(file.exists(files_output))) {
        next
      }

      time_start <- Sys.time()
      estimator  <- estimators[model == model_id]
      covariates <- estimator$covariates[[1]]
      fit_failed <- FALSE
      result     <- list(model_id = model_id, label = estimator$label,
                         covs = covariates, stage = "start")

      # 3.2 Sets, model formula and mesh. All models have the offset
      # log(hooks / 1000), year and HBF effects and a random vessel effect
      # (Eq. 10). The covariates of the model are added to these terms.
      sets <- copy(simulated[[scenario_id]]$sets)
      sets[, log_hooks := log(hooks / 1000)]
      sets[, year_f    := factor(year, levels = years)]
      sets[, vessel    := factor(vessel)]
      sets[, hbf_cl    := factor(hbf_cl, levels = hbf_levels)]
      model_formula <- reformulate(
        c("0", "year_f", "hbf_cl", "(1 | vessel)", covariates),
        response = "bet_n"
      )
      mesh <- make_mesh(sets, c("x", "y"), mesh = mesh_estimation)

      # 3.3 ST spatial+ only: each covariate is replaced by its residual
      # from a smooth of position (Eq. 16). The smooths are weighted by the
      # iterative weights of a separate fit of ST (X1, X2) to the same sets,
      # which must pass the convergence checks. If the fit of step 1 does
      # not pass them, its optimization is continued (step 2). If it still
      # does not pass and the SD of a field was at the boundary in step 1,
      # the model is refitted without that field (step 3) and, if needed,
      # the optimization of the refit is continued (step 4).
      smooths     <- list()
      first_stage <- NULL
      if (estimator$spatial_plus) {
        result$stage      <- "first-stage weights"
        steps_done        <- character(0)
        weights_converged <- FALSE
        spatial_on        <- TRUE
        spatiotemporal_on <- TRUE
        for (step in 1:4) {
          if (step == 3) {
            if (!any(boundary_first)) {
              break
            }
            spatial_on        <- !boundary_first[["sigma_O"]]
            spatiotemporal_on <- !boundary_first[["sigma_E"]]
          }
          if (step == 1 || step == 3) {
            weights_fit <- tryCatch(
              sdmTMB(
                model_formula,
                data           = sets,
                mesh           = mesh,
                time           = "year",
                offset         = "log_hooks",
                family         = nbinom2(),
                spatial        = if (spatial_on) "on" else "off",
                spatiotemporal = if (spatiotemporal_on) "iid" else "off",
                share_range    = TRUE,
                control        = sdmTMBcontrol(),
                silent         = TRUE
              ),
              error = function(e) e
            )
          } else {
            weights_fit <- tryCatch(
              run_extra_optimization(weights_fit_start, nlminb_loops = 2,
                                     newton_loops = 1),
              error = function(e) e
            )
          }
          steps_done <- c(steps_done, weights_steps[step])
          if (inherits(weights_fit, "error")) {
            if (step == 2) {
              next
            }
            break
          }
          if (step == 2 || step == 4) {
            # run_extra_optimization() does not update this element, which
            # sanity() reads.
            weights_fit$pos_def_hessian <-
              isTRUE(weights_fit$sd_report$pdHess)
          }

          # Convergence checks and, in step 1, whether a field SD is at
          # the boundary.
          checks <- unlist(suppressMessages(sanity(
            weights_fit,
            big_sd_log10    = log10(se_limit),
            gradient_thresh = gradient_limit,
            silent          = TRUE
          )))
          weights_converged <- isTRUE(all(checks[convergence_checks]))
          if (step == 1) {
            random_parameters <- as.data.table(tidy(weights_fit, "ran_pars"))
            sigma_O <- random_parameters[term == "sigma_O", estimate][1]
            sigma_E <- random_parameters[term == "sigma_E", estimate][1]
            boundary_first <- c(sigma_O = isTRUE(sigma_O < field_sd_limit),
                                sigma_E = isTRUE(sigma_E < field_sd_limit))
          }
          if (step == 1 || step == 3) {
            weights_fit_start <- weights_fit
          }
          if (weights_converged) {
            break
          }
        }

        if (!weights_converged) {
          fit_failed   <- TRUE
          result$error <- paste(
            "weights of the spatial+ smooths unavailable:",
            "the fit of ST (X1, X2) did not converge"
          )
          first_stage <- list(weights = list(converged_numeric = FALSE,
                                             attempts = steps_done))
        } else {
          # Weights w = mu / (1 + mu / phi) with the fitted mean mu (offset
          # included) and the dispersion phi, scaled to a mean of 1.
          fitted_mean    <- exp(predict(weights_fit)$est)
          phi            <- exp(weights_fit$model$par[["ln_phi"]])
          smooth_weights <- fitted_mean / (1 + fitted_mean / phi)
          smooth_weights <- smooth_weights / mean(smooth_weights)
          first_stage <- list(weights = list(
            phi               = phi,
            summary           = quantile(smooth_weights,
                                         c(0, 0.01, 0.25, 0.5, 0.75, 0.99, 1)),
            converged_numeric = TRUE,
            attempts          = steps_done,
            repaired          = length(steps_done) > 1
          ))

          # Thin plate regression spline of position fitted to each
          # covariate by REML. The covariate is replaced by its residual.
          for (covariate in covariates) {
            smooths[[covariate]] <- gam(
              reformulate(sprintf("s(x, y, k = %d)", k_first_stage),
                          response = covariate),
              data    = sets,
              weights = smooth_weights,
              method  = "REML"
            )
            set(sets, j = covariate,
                value = sets[[covariate]] - fitted(smooths[[covariate]]))
          }
        }
      }

      # 3.4 Fit, convergence checks and prediction. When the SD of a random
      # field is below field_sd_limit (at the boundary) and the fit is not
      # converged or the joint draws fail, the model is refitted without
      # that field (attempt 2) and the refit replaces the first fit
      # (Section 2.4). If the refit stops with an error, the first fit is
      # kept.
      predictions <- list()
      if (!fit_failed) {
        spatial_on        <- estimator$spatial
        spatiotemporal_on <- estimator$spatial
        for (attempt in 1:2) {
          result$stage <- if (attempt == 1) "fit" else "boundary refit"
          fit_attempt <- tryCatch(
            sdmTMB(
              model_formula,
              data           = sets,
              mesh           = mesh,
              time           = "year",
              offset         = "log_hooks",
              family         = nbinom2(),
              spatial        = if (spatial_on) "on" else "off",
              spatiotemporal = if (spatiotemporal_on) "iid" else "off",
              share_range    = TRUE,
              control        = sdmTMBcontrol(),
              silent         = TRUE
            ),
            error = function(e) e
          )
          if (inherits(fit_attempt, "error") && attempt == 1) {
            fit_failed   <- TRUE
            result$error <- conditionMessage(fit_attempt)
            break
          }

          if (inherits(fit_attempt, "error")) {
            result$refit_error <- conditionMessage(fit_attempt)
          } else {
            fit <- fit_attempt

            # Convergence checks, estimated SDs and range of the fields,
            # coefficients of the covariates and year effects.
            result$stage <- "diagnostics"
            checks <- unlist(suppressMessages(sanity(
              fit,
              big_sd_log10    = log10(se_limit),
              gradient_thresh = gradient_limit,
              silent          = TRUE
            )))
            converged <- isTRUE(all(checks[convergence_checks]))
            random_parameters <- as.data.table(tidy(fit, "ran_pars"))
            fixed_effects <- as.data.table(
              tidy(fit, "fixed", conf.int = TRUE)
            )
            sigma_O <- random_parameters[term == "sigma_O", estimate][1]
            sigma_E <- random_parameters[term == "sigma_E", estimate][1]
            sigma_V <- NA_real_
            if (any(random_parameters$group_name %in% "vessel")) {
              sigma_V <- random_parameters[group_name %in% "vessel",
                                           estimate][1]
            }
            boundary <- c(sigma_O = isTRUE(sigma_O < field_sd_limit),
                          sigma_E = isTRUE(sigma_E < field_sd_limit))

            result$diag <- list(
              sanity             = checks,
              nlminb_convergence = fit$model$convergence,
              nlminb_message     = fit$model$message,
              max_abs_gradient   = max(abs(fit$gradients)),
              pos_def_hessian    = isTRUE(fit$pos_def_hessian),
              sigma_O            = sigma_O,
              sigma_E            = sigma_E,
              range = random_parameters[term == "range", estimate][1],
              phi   = random_parameters[term == "phi", estimate][1],
              sigma_V            = sigma_V,
              control            = fit$control
            )
            result$converged_numeric <- converged
            result$par_fixed         <- fit$model$par
            result$par_all           <- fit$tmb_obj$env$last.par.best
            result$coef <- fixed_effects[
              term %in% covariates,
              .(term, estimate, std.error, conf.low, conf.high)
            ]
            result$year_eff <- fixed_effects[
              grepl("^year_f", term),
              .(year = as.integer(sub("^year_f", "", term)),
                year_eff = estimate)
            ]
            if (attempt == 1) {
              result$boundary_sigma <- any(boundary)
            }
            result$boundary_refit <- attempt == 2
          }

          # A first fit that is at the boundary and not converged is
          # refitted without being predicted. If the refit stopped with an
          # error, the first fit is predicted unless that was done before.
          refit_needed <- attempt == 1 && any(boundary) && !converged
          predict_now  <- !refit_needed &&
            (length(predictions) == 0 || !inherits(fit_attempt, "error"))

          if (predict_now) {
            result$stage <- "prediction"
            for (scenario_predicted in scenarios_predicted) {

              # 3.5 Prediction for one scenario. Prediction grid: all cells
              # and years at the reference levels of Eq. 3 (1,000 hooks,
              # HBF 16-17, no vessel effect). A density covariate keeps its
              # value in the cell and year. A catchability covariate is
              # set to 0 (Section 2.6). For ST spatial+ the inputs are
              # these values minus the smooth at the centre of the cell.
              role <- scenarios[scenario == scenario_predicted, role]
              grid <- merge(
                simulated[[scenario_predicted]]$grid[
                  , .(year, cell, X1, X2, area_km2)
                ],
                cell_centres,
                by = "cell"
              )
              setorder(grid, year, cell)
              grid[, year_f    := factor(year, levels = years)]
              grid[, log_hooks := 0]
              grid[, hbf_cl    := factor(hbf_levels[1], levels = hbf_levels)]
              grid[, vessel    := NA]
              if (role == "catchability") {
                grid[, X1 := 0]
                grid[, X2 := 0]
              }
              if (estimator$spatial_plus) {
                for (covariate in covariates) {
                  smooth_cells <- as.vector(
                    predict(smooths[[covariate]], newdata = grid)
                  )
                  set(grid, j = covariate,
                      value = grid[[covariate]] - smooth_cells)
                }
              }

              # Predicted log density in every cell and year (eta) and the
              # epsilon bias-corrected index over all cells (Section 2.6).
              prediction <- predict(
                fit,
                newdata           = grid,
                re_form_iid       = NA,
                offset            = grid$log_hooks,
                return_tmb_object = TRUE
              )
              index <- as.data.table(get_index(
                prediction,
                area         = grid$area_km2,
                bias_correct = TRUE
              ))
              index <- index[, .(year, est, lwr, upr, log_est, se)]
              grid_pred <- data.table(year = grid$year, cell = grid$cell,
                                      eta = prediction$data$est)

              # Seed of the joint draws, from the scenario, the replicate
              # and the model, so that the draws do not depend on the order
              # of the fits (see the header).
              seed <- as.integer(
                base_seed +
                  100000 * match(scenario_predicted, scenarios$scenario) +
                  10 * replicate_id +
                  match(model_id, estimators$model)
              )

              # Joint draws of the predicted log density from the joint
              # precision of fixed and random effects. For each draw the
              # log index of every year is log(sum(area * exp(eta))),
              # computed after subtracting the column maximum. Intervals
              # are given for the log index (unc) and for the log index
              # centred on its mean over years within each draw (cen).
              set.seed(seed)
              draws <- tryCatch(
                predict(fit, newdata = grid, nsim = n_draws_index,
                        re_form_iid = NA, offset = grid$log_hooks),
                error = function(e) e
              )
              logI_draws  <- NULL
              draws_error <- NULL
              if (inherits(draws, "error")) {
                draws_error <- conditionMessage(draws)
              } else {
                stopifnot(!anyNA(draws), !anyNA(grid$area_km2))
                logI_draws <- matrix(NA_real_, n_years, ncol(draws),
                                     dimnames = list(years, NULL))
                for (i_year in seq_len(n_years)) {
                  rows      <- which(grid$year == years[i_year])
                  log_terms <- draws[rows, , drop = FALSE] +
                    log(grid$area_km2[rows])
                  log_max   <- apply(log_terms, 2, max)
                  logI_draws[i_year, ] <- log_max +
                    log(colSums(exp(sweep(log_terms, 2, log_max))))
                }
                centred_draws <- sweep(logI_draws, 2, colMeans(logI_draws))
                limits_unc <- apply(logI_draws, 1, quantile,
                                    probs = c(0.025, 0.975))
                limits_cen <- apply(centred_draws, 1, quantile,
                                    probs = c(0.025, 0.975))
                index[, unc_sd  := apply(logI_draws, 1, sd)]
                index[, unc_lwr := limits_unc[1, ]]
                index[, unc_upr := limits_unc[2, ]]
                index[, cen_sd  := apply(centred_draws, 1, sd)]
                index[, cen_lwr := limits_cen[1, ]]
                index[, cen_upr := limits_cen[2, ]]
              }

              predictions[[scenario_predicted]] <- list(
                role = role, seed = seed, index = index,
                grid_pred = grid_pred
              )
              predictions[[scenario_predicted]]$logI_draws  <- logI_draws
              predictions[[scenario_predicted]]$draws_error <- draws_error
            }

            # The joint draws use the same fit for every scenario, so the
            # draws of the fitted scenario show whether they failed.
            draws_failed <- !is.null(predictions[[scenario_id]]$draws_error)
            refit_needed <- attempt == 1 && any(boundary) && draws_failed
          }

          if (!refit_needed) {
            break
          }

          # Keep the checks of the fit at the boundary and switch off the
          # field whose SD is at the boundary.
          result$boundary_original <- list(
            diag        = result$diag,
            coef        = result$coef,
            par_fixed   = result$par_fixed,
            draws_error = predictions[[scenario_id]]$draws_error
          )
          spatial_on        <- spatial_on && !boundary[["sigma_O"]]
          spatiotemporal_on <- spatiotemporal_on && !boundary[["sigma_E"]]
          result$refit_spec <- c(
            spatial        = if (spatial_on) "on" else "off",
            spatiotemporal = if (spatiotemporal_on) "iid" else "off"
          )
        }
      }

      # 3.6 Save one file per predicted scenario. Elements of the list:
      #   model_id, label, covs   model, its label and its covariates
      #   scenario_id, rep        scenario and replicate
      #   ok, stage, error        FALSE, the step reached and the message
      #                           when the fit stopped with an error
      #   converged_numeric       TRUE when the fit passed the checks
      #   diag                    checks of sanity(), gradient, SDs and
      #                           range of the fields, optimizer settings
      #   boundary_sigma          TRUE when a field SD of the first fit was
      #                           at the boundary
      #   boundary_refit          TRUE when the result is the refit
      #   refit_spec              fields of the refit (on or off)
      #   boundary_original       checks and coefficients of the first fit
      #                           when a refit was made
      #   coef                    coefficients of X1 and X2 with standard
      #                           errors and 95% confidence limits
      #   year_eff                estimated year effects
      #   par_fixed, par_all      fixed parameters, and all parameters
      #                           with the random effects (used by 05_rsr.R)
      #   first_stage             ST spatial+: dispersion and quantiles of
      #                           the weights, steps taken for the weights
      #   role, seed              role of the covariates in the index and
      #                           seed of the joint draws
      #   index                   by year: bias-corrected index over all
      #                           1,592 cells (est, lwr, upr, log_est, se)
      #                           and SD and 95% limits of the joint draws
      #   grid_pred               predicted log density (eta) by year and
      #                           cell
      #   logI_draws              joint draws of the log index (years by
      #                           draws), or draws_error when they failed
      #   truth, beta_true,       true index, true coefficients and true
      #   grid_truth              log density by cell (from the simulation)
      #   shared_fit              scenarios predicted from the same fit
      #   minutes                 run time
      result$ok <- !fit_failed
      if (!fit_failed) {
        result$stage <- "done"
      }
      result$first_stage <- first_stage
      result$shared_fit  <- NA_character_
      if (length(scenarios_predicted) > 1) {
        result$shared_fit <- paste(scenarios_predicted, collapse = "/")
      }
      result$minutes <- as.numeric(
        difftime(Sys.time(), time_start, units = "mins")
      )
      for (scenario_predicted in scenarios_predicted) {
        output <- c(result, predictions[[scenario_predicted]])
        output$scenario_id <- scenario_predicted
        output$rep         <- replicate_id
        output$truth       <- simulated[[scenario_predicted]]$truth
        true_coefficients  <- simulated[[scenario_predicted]]$params$coef
        output$beta_true   <- c(
          X1 = true_coefficients[["b1d"]] + true_coefficients[["b1q"]],
          X2 = true_coefficients[["b2d"]] + true_coefficients[["b2q"]]
        )
        output$grid_truth <- simulated[[scenario_predicted]]$grid[
          , .(year, cell, log_d, area_km2, sampled)
        ]

        # Write to a temporary file first, so that an interrupted run
        # leaves no incomplete output file.
        file_output <- file.path(
          dir_fits_run, scenario_predicted,
          sprintf("%s_rep%04d.rds", model_id, replicate_id)
        )
        dir.create(dirname(file_output), recursive = TRUE,
                   showWarnings = FALSE)
        saveRDS(output, paste0(file_output, ".tmp"))
        file.rename(paste0(file_output, ".tmp"), file_output)
      }

      cat(sprintf(
        "%s %s replicate %d %s: ok %s, converged %s, refit %s, %.1f min\n",
        format(Sys.time(), "%H:%M"),
        paste(scenarios_predicted, collapse = "/"), replicate_id, model_id,
        result$ok, isTRUE(result$converged_numeric),
        isTRUE(result$boundary_refit), result$minutes
      ))
    }
  }
}
