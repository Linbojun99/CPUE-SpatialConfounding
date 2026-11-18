# ==============================================================================
# 07_x1_density_and_catchability.R
#
# Sensitivity analysis of E1 and E2 in which X1 affects both density and
# catchability (Section 2.3; Supplementary material, "Paired comparison of
# density and catchability indices"; Supplementary Table S3). E1 and E2 have
# the same simulated sets, catches and fits. They differ only in the role of
# X1: its effect is on density in E1 and on catchability in E2. Here half of
# the true effect of X1 is assigned to density and half to catchability. The
# effect on the catch is unchanged, so nothing is simulated or fitted again.
# For each replicate and estimator the script
#   - calculates the true index from a true log density that includes half
#     of the effect of X1 (Eq. 9),
#   - starts from the predicted log density of E1, which contains the whole
#     fitted effect of X1, and removes the part of that effect assigned to
#     catchability, so that the prediction retains 0%, 50% or 100% of the
#     fitted effect,
#   - calculates the uncorrected index over the prediction domain (Eq. 19),
#     the log-index RMSE (Eq. 20) and the change error (Eq. 22).
# The median RMSE and the mean change error over replicates are the entries
# of Table S3.
#
# Input:
#   data/prediction_domain_cells.csv              cells of the simulation
#                                                 grid, in_domain marks the
#                                                 prediction domain
#   output/simulation/E1/rep<NNNN>.rds            simulated data of E1
#   output/simulation/E2/rep<NNNN>.rds            simulated data of E2
#   output/fits/E1/<model>_rep<NNNN>.rds          fits of the six models
#   output/rsr/rsr_coefficients.csv               RSR coefficients (05_rsr.R)
# Output (in output/x1_density_and_catchability/):
#   x1_density_and_catchability_index_years.csv   one row per replicate,
#                                                 estimator, retained
#                                                 fraction and year: true
#                                                 and estimated log index
#   x1_density_and_catchability_replicates.csv    one row per replicate,
#                                                 estimator and retained
#                                                 fraction: RMSE and change
#                                                 error
#   x1_density_and_catchability_summary.csv       by estimator and retained
#                                                 fraction
# Output (in output/tables/):
#   table_s3_x1_density_and_catchability.csv      Supplementary Table S3
#
# Run time: less than a minute to a few minutes, mostly the reading of the
# 350 fit files.
# ==============================================================================

source("00_settings.R")


# 1. Settings for this run -----------------------------------------------------

# Replicates to use. The paper uses replicates 1 to 50.
run_replicates <- 1:n_replicates


# 2. Fixed settings and shared inputs ------------------------------------------

# The two scenarios with the same sets, catches and fits. The effect of X1
# is on density in the first and on catchability in the second.
scenario_density      <- "E1"
scenario_catchability <- "E2"

# Share of the true effect of X1 assigned to density. The rest is assigned
# to catchability.
density_share <- 0.5

# Fractions of the fitted effect of X1 retained in the prediction (columns
# of Table S3). The fraction 0.5 equals the share in the operating model.
retained_fractions <- c(0, 0.5, 1)

# Estimators: the six fitted models and RSR, which is obtained from the fit
# of ST (X1, X2) (05_rsr.R).
rsr_source_model <- "M4"
estimator_ids    <- c(estimators$model, "RSR")
estimator_labels <- c(estimators$label, "RSR")

# Order of the rows of Table S3.
table_order <- c("M1", "M2", "M3", "M4", "M6", "M5", "RSR")

# Columns of the simulated sets that must be the same in the two scenarios.
columns_shared <- c("year", "set_src", "vessel", "hooks", "hbf", "X1", "X2",
                    "bet_n")

# Cells of the prediction domain (1,484 of the 1,592 cells of the grid).
domain_table <- fread(file.path(dir_data, "prediction_domain_cells.csv"))
domain_cells <- domain_table[in_domain == TRUE, cell]

# RSR coefficients of the density scenario. The RSR coefficients of the two
# scenarios are the same because the source fits are the same.
rsr_coefficients <- fread(file.path(dir_rsr, "rsr_coefficients.csv"))
rsr_coefficients <- rsr_coefficients[scenario_id == scenario_density]


# 3. Indices of each replicate and estimator -----------------------------------

rows_years <- list()

for (replicate_id in run_replicates) {

  # 3.1 Simulated data of the two scenarios, with the cells in the same
  # order (year, then cell).
  simulated_density <- readRDS(file.path(
    dir_simulation, scenario_density, sprintf("rep%04d.rds", replicate_id)
  ))
  simulated_catchability <- readRDS(file.path(
    dir_simulation, scenario_catchability,
    sprintf("rep%04d.rds", replicate_id)
  ))
  grid_density      <- simulated_density$grid[order(year, cell)]
  grid_catchability <- simulated_catchability$grid[order(year, cell)]

  # The two scenarios have the same sets, catches and covariates.
  for (column in columns_shared) {
    stopifnot(identical(simulated_density$sets[[column]],
                        simulated_catchability$sets[[column]]))
  }
  stopifnot(
    identical(grid_density$year, grid_catchability$year),
    identical(grid_density$cell, grid_catchability$cell),
    identical(grid_density$X1, grid_catchability$X1),
    identical(grid_density$X2, grid_catchability$X2),
    identical(grid_density$area_km2, grid_catchability$area_km2)
  )

  # 3.2 True log density when X1 affects both density and catchability. The
  # total effect of X1 is its density coefficient in E1 (Eq. 3), which
  # equals its catchability coefficient in E2 (Eq. 4). The true log density
  # of E2 has no effect of X1, so the share assigned to density is added to
  # it.
  beta1_total <- simulated_density$params$coef[["b1d"]]
  stopifnot(
    identical(beta1_total, simulated_catchability$params$coef[["b1q"]])
  )
  log_density_true <- grid_catchability$log_d +
    density_share * beta1_total * grid_density$X1

  # 3.3 True index over the prediction domain in each year (Eq. 9).
  in_domain <- grid_density$cell %in% domain_cells
  area      <- grid_density$area_km2
  log_true  <- numeric(n_years)
  for (i_year in seq_len(n_years)) {
    rows <- which(grid_density$year == years[i_year] & in_domain)
    log_true[i_year] <-
      log(sum(area[rows] * exp(log_density_true[rows])))
  }

  for (estimator_id in estimator_ids) {

    # 3.4 Fit of the density scenario. Its predicted log density uses the
    # covariates at their values in each cell and year, so it contains the
    # whole fitted effect of X1. RSR has the same predicted log density as
    # its source model ST (X1, X2) at these covariate values (Section 2.5).
    model_id <- estimator_id
    if (estimator_id == "RSR") {
      model_id <- rsr_source_model
    }
    fit <- readRDS(file.path(
      dir_fits, scenario_density,
      sprintf("%s_rep%04d.rds", model_id, replicate_id)
    ))
    if (!isTRUE(fit$ok) || !isTRUE(fit$converged_numeric)) {
      next
    }
    predicted <- fit$grid_pred[order(year, cell)]
    stopifnot(
      identical(predicted$year, grid_density$year),
      identical(predicted$cell, grid_density$cell)
    )

    # 3.5 Coefficient of X1 used to remove the catchability part of its
    # effect: the fitted coefficient, the RSR coefficient (Eq. 17) for RSR,
    # and zero for a model without X1, whose prediction then does not
    # depend on the retained fraction.
    beta1_hat <- 0
    if ("X1" %in% fit$coef$term) {
      beta1_hat <- fit$coef$estimate[fit$coef$term == "X1"]
    }
    if (estimator_id == "RSR") {
      beta1_hat <- rsr_coefficients[
        rep == replicate_id & term == "X1", estimate
      ]
      if (length(beta1_hat) != 1) {
        next
      }
    }

    for (retained in retained_fractions) {

      # 3.6 Predicted log density that retains the given fraction of the
      # fitted effect of X1. X2 and the fitted spatial terms are unchanged.
      log_density_hat <- predicted$eta +
        (retained - 1) * beta1_hat * grid_density$X1

      # 3.7 Uncorrected index over the prediction domain (Eq. 19).
      log_est <- numeric(n_years)
      for (i_year in seq_len(n_years)) {
        rows <- which(grid_density$year == years[i_year] & in_domain)
        log_est[i_year] <-
          log(sum(area[rows] * exp(log_density_hat[rows])))
      }

      rows_years[[paste(replicate_id, estimator_id, retained)]] <-
        data.table(
          rep      = replicate_id,
          model    = estimator_id,
          retained = retained,
          year     = years,
          log_true = log_true,
          log_est  = log_est
        )
    }
  }
  cat("Replicate", replicate_id, "finished\n")
}

index_years <- rbindlist(rows_years)
setorder(index_years, rep, model, retained, year)

fwrite(index_years, file.path(
  dir_x1_both, "x1_density_and_catchability_index_years.csv"
))


# 4. Measures of each replicate ------------------------------------------------

# The error of the log index is centred on its mean over years, which is the
# difference between the centred estimated and true series. rmse is the
# log-index RMSE (Eq. 20). change_error is the error in the last year minus
# that in the first year (Eq. 22).
index_years[, error := log_est - log_true]
replicate_measures <- index_years[, .(
  rmse         = sqrt(mean((error - mean(error))^2)),
  change_error = error[.N] - error[1]
), by = .(rep, model, retained)]

fwrite(replicate_measures, file.path(
  dir_x1_both, "x1_density_and_catchability_replicates.csv"
))


# 5. Summary over replicates and Table S3 --------------------------------------

# Median RMSE and mean change error by estimator and retained fraction.
measure_summary <- replicate_measures[, .(
  n                 = .N,
  median_rmse       = median(rmse),
  mean_change_error = mean(change_error)
), by = .(model, retained)]
setorder(measure_summary, model, retained)

fwrite(measure_summary, file.path(
  dir_x1_both, "x1_density_and_catchability_summary.csv"
))

# Table S3: one row per estimator, the RMSE and then the change error for
# the three retained fractions, rounded to three decimal places.
table_s3 <- data.table(
  estimator = estimator_labels[match(table_order, estimator_ids)]
)
for (measure in c("median_rmse", "mean_change_error")) {
  column_start <- "rmse"
  if (measure == "mean_change_error") {
    column_start <- "change_error"
  }
  for (retained_fraction in retained_fractions) {
    values <- measure_summary[retained == retained_fraction]
    values <- values[match(table_order, model)]
    set(
      table_s3,
      j     = paste0(column_start, "_", 100 * retained_fraction),
      value = round(values[[measure]], 3)
    )
  }
}

fwrite(table_s3, file.path(
  dir_tables, "table_s3_x1_density_and_catchability.csv"
))
print(table_s3)
