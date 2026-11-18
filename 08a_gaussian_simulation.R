# ==============================================================================
# 08a_gaussian_simulation.R
#
# The Gaussian simulation (Section 2.7 and Supplementary material, "Gaussian
# simulation on a unit square"). For every replicate the script draws the
# locations, the spatial effect omega and the covariates X1 and X2 on a unit
# square, simulates Gaussian observations on the log scale for the scenarios
# G0, G1, G1t and G2, fits a Gaussian GLM, a spatial model with a thin plate
# regression spline and spatial+, derives RSR from the spatial model
# (Eq. 17), and calculates the index and the performance measures of
# Section 2.6 for the four estimators.
#
# Eight runs are made: the calibration of the basis dimension k (G0 only,
# k = 25, 50, 100, replicates independent of all other runs), the reporting
# runs (four scenarios, k = 25, 50, 100, of which k = 50 gives Fig. S16 and
# Table S9) and two sensitivity runs at k = 50 without the annual component
# of X2 or of both covariates (Table S5).
#
# Input:
#   none (all data are simulated in this script)
# Output (one folder per run, output/gaussian/<setting>/k<k>/):
#   coefficients.csv    estimates, SEs, coverage and Wald tests of the
#                       coefficients of X1 and X2
#   indices.csv         estimated and true log index by year
#   performance.csv     log-index RMSE (Eq. 20), change error (Eq. 22) and
#                       RMSE of the predicted log density
#   approximation.csv   first-order approximation of the change error in G2
#                       (reporting runs only)
#   diagnostics.csv     theoretical bias on the grid, share of the variance
#                       of X1 left in the spatial+ first-stage residuals and
#                       checks of the RSR projection
#
# Run time: about 10 minutes on one core.
# ==============================================================================

source("00_settings.R")
library(mgcv)


# 1. Settings for this run -----------------------------------------------------

# Runs to carry out (row numbers of the table `runs` in section 2) and
# replicates of each run. The paper uses all eight runs and replicates 1 to
# 50.
run_rows       <- 1:8
run_replicates <- 1:n_replicates


# 2. Design of the Gaussian simulation -----------------------------------------

n_years_gaussian <- 5      # years
n_locations_year <- 60     # locations sampled in each year
n_locations      <- n_years_gaussian * n_locations_year
grid_side        <- 20     # the index is integrated over a 20 x 20 grid
n_cells          <- grid_side^2
n_points         <- n_cells + n_locations

beta_true      <- c(X1 = 0.5, X2 = 0)   # true coefficients of X1 and X2
sd_omega       <- 0.7                   # SD of the spatial effect omega
sd_observation <- 0.5                   # SD of the Gaussian response error

# Year effects alpha_t: a 30% decline between the first and last year.
time_fraction <- seq(0, 1, length.out = n_years_gaussian)
alpha_true    <- log(0.7) * time_fraction

# Matern correlation with smoothness 3/2, C(d) = (1 + z) exp(-z) with
# z = matern_scale * d / matern_range. The range is the distance at which the
# correlation is 0.1, so matern_scale is the solution of
# (1 + z) exp(-z) = 0.1, written here to full double precision.
matern_range  <- 0.25
matern_scale  <- 3.8897201698674273
matern_jitter <- 1e-10    # added to the diagonal before the factorization
stopifnot(abs((1 + matern_scale) * exp(-matern_scale) - 0.1) < 1e-12)

# Scenarios.
#   rho_omega  correlation between the time-invariant part of X1 and omega
#   x1_trend   increase in X1 between the first and last year
#   role       X1 is a density covariate (kept at its values in prediction)
#              or a catchability covariate (set to 0 in prediction)
gaussian_scenarios <- data.table(
  scenario  = c("G0", "G1", "G1t", "G2"),
  rho_omega = c(0, 0.6, 0.6, 0.6),
  x1_trend  = c(0, 0, 1, 1),
  role      = c("density", "density", "density", "catchability")
)

# Runs.
#   setting       "calibration"     selection of k (Table S5)
#                 "reporting"       results of the paper (Fig. S16, Tables S5
#                                   and S9)
#                 "x1_annual"       sensitivity: no annual component in X2
#                 "time_invariant"  sensitivity: no annual component in X1
#                                   and X2
#   k             basis dimension of all thin plate regression splines
#   seed          the seed of replicate r is seed + r
#   sd_annual_x1  SD over the grid of the annual component of X1
#   sd_annual_x2  SD over the grid of the annual component of X2
#   scenarios     scenarios simulated in the run
# The two sensitivity runs use k = 50, the value selected by the calibration
# (08b_gaussian_summaries.R, section 3).
runs <- data.table(
  setting      = c("calibration", "calibration", "calibration",
                   "reporting", "reporting", "reporting",
                   "x1_annual", "time_invariant"),
  k            = c(25, 50, 100, 25, 50, 100, 50, 50),
  seed         = c(1370001, 1370001, 1370001,
                   270001, 270001, 270001, 270001, 270001),
  sd_annual_x1 = c(0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0),
  sd_annual_x2 = c(0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0, 0),
  scenarios    = list("G0", "G0", "G0",
                      c("G0", "G1", "G1t", "G2"),
                      c("G0", "G1", "G1t", "G2"),
                      c("G0", "G1", "G1t", "G2"),
                      c("G0", "G1"), c("G0", "G1"))
)

# Estimators, in the order in which they are written to the output files.
#   "GLM"       Gaussian generalized linear model without a spatial term
#   "SGLMM"     the spatial model: thin plate regression spline, REML
#   "RSR"       restricted spatial regression derived from the spatial model
#   "spatial+"  spatial+ with first-stage smooths of X1 and X2
gaussian_estimators <- c("GLM", "SGLMM", "RSR", "spatial+")
covariate_names     <- c("X1", "X2")
residual_names      <- c("X1_residual", "X2_residual")

# Cell centres of the grid on the unit square. Rows 1 to n_cells of every
# object indexed by "points" below are the cells, the remaining rows are the
# sampled locations.
cell_centres  <- (seq_len(grid_side) - 0.5) / grid_side
grid_cells    <- expand.grid(x = cell_centres, y = cell_centres)
cell_rows     <- seq_len(n_cells)
location_rows <- n_cells + seq_len(n_locations)

# Year of every sampled location and of every row of the prediction grid
# (the cells repeated for each year).
location_year <- rep(seq_len(n_years_gaussian), each = n_locations_year)
grid_year     <- rep(seq_len(n_years_gaussian), each = n_cells)
grid_cell     <- rep(cell_rows, times = n_years_gaussian)


# 3. Simulation, fits and performance measures ---------------------------------

for (run in run_rows) {

  setting_name     <- runs$setting[run]
  k_basis          <- runs$k[run]
  run_seed         <- runs$seed[run]
  sd_annual        <- c(runs$sd_annual_x1[run], runs$sd_annual_x2[run])
  scenarios_in_run <- runs$scenarios[[run]]

  dir_run <- file.path(dir_gaussian, setting_name, paste0("k", k_basis))
  dir.create(dir_run, recursive = TRUE, showWarnings = FALSE)
  cat("Run", run, ":", setting_name, "k =", k_basis, "\n")

  coefficient_rows   <- list()
  index_rows         <- list()
  performance_rows   <- list()
  approximation_rows <- list()
  diagnostic_rows    <- list()

  for (replicate_id in run_replicates) {

    # 3.1 Sampled locations. The same seed is used for all scenarios of a
    # replicate, so the scenarios share locations, fields and response
    # errors.
    set.seed(run_seed + replicate_id)
    observations <- data.frame(
      x    = runif(n_locations),
      y    = runif(n_locations),
      year = factor(location_year, levels = seq_len(n_years_gaussian))
    )

    # 3.2 Three Matern fields drawn jointly at the cell centres and the
    # sampled locations.
    all_points  <- rbind(grid_cells, observations[, c("x", "y")])
    distance    <- as.matrix(dist(all_points))
    scaled      <- matern_scale * distance / matern_range
    correlation <- (1 + scaled) * exp(-scaled)
    diag(correlation) <- diag(correlation) + matern_jitter
    cholesky    <- chol(correlation)
    raw_fields  <- t(cholesky) %*% matrix(rnorm(n_points * 3), ncol = 3)

    # 3.3 Centre, orthogonalize and scale the three fields over the grid
    # (Gram-Schmidt). The result is the standardized spatial effect ("omega")
    # and the components "Z" and "W" of the time-invariant parts of X1 and
    # X2, which have mean 0, mean square 1 and are mutually orthogonal over
    # the grid. The same transformation is applied at the locations.
    fields <- matrix(0, n_points, 3,
                     dimnames = list(NULL, c("omega", "Z", "W")))
    for (j in 1:3) {
      field_j <- raw_fields[, j] - mean(raw_fields[cell_rows, j])
      if (j > 1) {
        for (h in seq_len(j - 1)) {
          projection_h <- sum(field_j[cell_rows] * fields[cell_rows, h]) /
            sum(fields[cell_rows, h]^2)
          field_j <- field_j - projection_h * fields[, h]
        }
      }
      fields[, j] <- field_j / sqrt(mean(field_j[cell_rows]^2))
    }

    # 3.4 Response errors. They are drawn here because the annual components
    # below use their own seeds. The random number stream of this seed is
    # locations, fields, response errors.
    response_error <- rnorm(n_locations, sd = sd_observation)

    # 3.5 Annual components of X1 (v = 1) and X2 (v = 2): one independent
    # Matern field per year with the same correlation function, centred and
    # scaled to the given SD over the grid in each year. Each covariate has
    # its own seed, so that switching one component off leaves everything
    # else unchanged.
    annual <- array(0, dim = c(n_points, n_years_gaussian, 2))
    for (v in 1:2) {
      if (sd_annual[v] > 0) {
        set.seed(run_seed + replicate_id + v * 1000000)
        annual_raw <- t(cholesky) %*%
          matrix(rnorm(n_points * n_years_gaussian), ncol = n_years_gaussian)
        for (i_year in seq_len(n_years_gaussian)) {
          annual[, i_year, v] <- sd_annual[v] *
            (annual_raw[, i_year] - mean(annual_raw[cell_rows, i_year])) /
            sd(annual_raw[cell_rows, i_year])
        }
      }
    }
    annual_location_index <- cbind(location_rows, location_year)
    annual_x1_locations   <- annual[, , 1][annual_location_index]
    annual_x2_locations   <- annual[, , 2][annual_location_index]
    annual_x1_grid        <- as.vector(annual[cell_rows, , 1])
    annual_x2_grid        <- as.vector(annual[cell_rows, , 2])

    for (scenario_id in scenarios_in_run) {

      scenario_row <- match(scenario_id, gaussian_scenarios$scenario)
      rho_omega    <- gaussian_scenarios$rho_omega[scenario_row]
      x1_trend     <- gaussian_scenarios$x1_trend[scenario_row]
      role         <- gaussian_scenarios$role[scenario_row]

      # 3.6 Covariates, spatial effect and response at the sampled locations.
      # X1 is the sum of its time-invariant part (correlation rho_omega with
      # omega over the grid), its annual component and the trend. The
      # response includes the effect of X1 in every scenario.
      observations$omega <- sd_omega * fields[location_rows, "omega"]
      observations$X1 <- rho_omega * fields[location_rows, "omega"] +
        sqrt(1 - rho_omega^2) * fields[location_rows, "Z"] +
        annual_x1_locations + x1_trend * time_fraction[location_year]
      observations$X2 <- fields[location_rows, "W"] + annual_x2_locations
      observations$response <- alpha_true[location_year] +
        observations$omega +
        beta_true[["X1"]] * observations$X1 +
        beta_true[["X2"]] * observations$X2 +
        response_error

      # 3.7 The same quantities in every cell and year, and the true log
      # density. X1 contributes to density when it is a density covariate
      # and not when it is a catchability covariate (G2).
      grid <- data.frame(
        x    = grid_cells$x[grid_cell],
        y    = grid_cells$y[grid_cell],
        year = factor(grid_year, levels = seq_len(n_years_gaussian)),
        area = 1 / n_cells
      )
      grid$omega <- sd_omega * fields[grid_cell, "omega"]
      grid$X1 <- rho_omega * fields[grid_cell, "omega"] +
        sqrt(1 - rho_omega^2) * fields[grid_cell, "Z"] +
        annual_x1_grid + x1_trend * time_fraction[grid_year]
      grid$X2 <- fields[grid_cell, "W"] + annual_x2_grid
      if (role == "density") {
        grid$log_density <- alpha_true[grid_year] + grid$omega +
          beta_true[["X1"]] * grid$X1
      } else {
        grid$log_density <- alpha_true[grid_year] + grid$omega
      }

      # G2 has the same observations as G1t and uses the fits of G1t. Only
      # the prediction (3.11) and the true density differ, so steps 3.8 to
      # 3.10 are skipped for G2.
      if (scenario_id != "G2") {

        # 3.8 Gaussian GLM and the spatial model.
        glm_fit <- glm(response ~ year + X1 + X2, data = observations,
                       family = gaussian())
        spatial_fit <- gam(
          response ~ year + X1 + X2 + s(x, y, bs = "tp", k = k_basis),
          data = observations, family = gaussian(), method = "REML"
        )

        # 3.9 Spatial+ (Dupont et al., 2022). First stage: each covariate is
        # regressed on a thin plate regression spline of the coordinates.
        # Second stage: the spatial model with the covariates replaced by
        # their first-stage residuals.
        first_stage_x1 <- gam(X1 ~ s(x, y, bs = "tp", k = k_basis),
                              data = observations, method = "REML")
        first_stage_x2 <- gam(X2 ~ s(x, y, bs = "tp", k = k_basis),
                              data = observations, method = "REML")
        observations$X1_residual <- observations$X1 -
          predict(first_stage_x1, newdata = observations)
        observations$X2_residual <- observations$X2 -
          predict(first_stage_x2, newdata = observations)
        spatial_plus_fit <- gam(
          response ~ year + X1_residual + X2_residual +
            s(x, y, bs = "tp", k = k_basis),
          data = observations, family = gaussian(), method = "REML"
        )
        stopifnot(glm_fit$converged, spatial_fit$converged,
                  first_stage_x1$converged, first_stage_x2$converged,
                  spatial_plus_fit$converged)

        # 3.10 RSR (Eq. 17). The fitted spatial effect at the locations is
        # projected onto the two covariate columns. The projection is added
        # to the coefficients and removed from the spatial effect, which
        # leaves the fitted values unchanged.
        x_locations    <- as.matrix(observations[, covariate_names])
        design_matrix  <- predict(spatial_fit, newdata = observations,
                                  type = "lpmatrix")
        spline_columns <- seq.int(spatial_fit$smooth[[1]]$first.para,
                                  spatial_fit$smooth[[1]]$last.para)
        spline_basis   <- design_matrix[, spline_columns, drop = FALSE]
        spatial_effect <- drop(
          spline_basis %*% coef(spatial_fit)[spline_columns]
        )

        # Projection matrix (X'X)^-1 X', from the QR decomposition of X.
        qr_x       <- qr(x_locations, tol = 1e-9)
        stopifnot(qr_x$rank == 2)
        projection <- backsolve(qr.R(qr_x), t(qr.Q(qr_x)))
        projection <- projection[order(qr_x$pivot), , drop = FALSE]

        delta_beta <- drop(projection %*% spatial_effect)
        beta_rsr   <- as.numeric(coef(spatial_fit)[covariate_names]) +
          delta_beta
        spatial_effect_rsr <- spatial_effect -
          drop(x_locations %*% delta_beta)
        names(delta_beta) <- covariate_names
        names(beta_rsr)   <- covariate_names

        # Check: fitted values before and after the adjustment.
        fitted_spatial <- drop(design_matrix %*% coef(spatial_fit))
        other_terms    <- fitted_spatial -
          drop(x_locations %*% coef(spatial_fit)[covariate_names]) -
          spatial_effect
        fitted_rsr     <- other_terms + drop(x_locations %*% beta_rsr) +
          spatial_effect_rsr
        rsr_fitted_max_abs <- max(abs(fitted_spatial - fitted_rsr))

        # Covariance of the RSR coefficients by the delta method. The RSR
        # coefficients are a linear function of all coefficients of the
        # spatial model. The Jacobian has 1 for the coefficient itself and
        # the projected spline basis for the spline coefficients. It is
        # applied to the covariance that includes smoothing-parameter
        # uncertainty (unconditional) and to the one that does not.
        jacobian <- matrix(
          0, 2, length(coef(spatial_fit)),
          dimnames = list(covariate_names, names(coef(spatial_fit)))
        )
        jacobian[cbind(1:2, match(covariate_names, colnames(jacobian)))] <- 1
        jacobian[, spline_columns] <- projection %*% spline_basis

        covariance_rsr <- jacobian %*%
          vcov(spatial_fit, unconditional = TRUE) %*% t(jacobian)
        covariance_rsr <- (covariance_rsr + t(covariance_rsr)) / 2
        diag(covariance_rsr) <- pmax(0, diag(covariance_rsr))

        covariance_rsr_conditional <- jacobian %*%
          vcov(spatial_fit, unconditional = FALSE) %*% t(jacobian)
        covariance_rsr_conditional <-
          (covariance_rsr_conditional + t(covariance_rsr_conditional)) / 2
        diag(covariance_rsr_conditional) <-
          pmax(0, diag(covariance_rsr_conditional))

        # Coefficients of X1 and X2 and their covariance for each estimator.
        # For the two GAM fits the covariance includes smoothing-parameter
        # uncertainty, and the conditional version is kept for comparison.
        # The spatial+ covariance is conditional on the first-stage fits.
        beta_hat <- list(
          GLM        = coef(glm_fit)[covariate_names],
          SGLMM      = coef(spatial_fit)[covariate_names],
          RSR        = beta_rsr,
          `spatial+` = coef(spatial_plus_fit)[residual_names]
        )
        covariance_hat <- list(
          GLM        = vcov(glm_fit)[covariate_names, covariate_names],
          SGLMM      = vcov(spatial_fit, unconditional = TRUE)[
            covariate_names, covariate_names],
          RSR        = covariance_rsr,
          `spatial+` = vcov(spatial_plus_fit, unconditional = TRUE)[
            residual_names, residual_names]
        )
        covariance_conditional <- list(
          GLM        = vcov(glm_fit)[covariate_names, covariate_names],
          SGLMM      = vcov(spatial_fit, unconditional = FALSE)[
            covariate_names, covariate_names],
          RSR        = covariance_rsr_conditional,
          `spatial+` = vcov(spatial_plus_fit, unconditional = FALSE)[
            residual_names, residual_names]
        )
      }

      # 3.11 Predicted log density in every cell and year. A density
      # covariate enters at its value in the cell and year, and a
      # catchability covariate (X1 in G2) is set to 0. For RSR the adjusted
      # spatial effect is defined with the original covariates, so the RSR
      # predictor differs from that of the spatial model by
      # (X_reference - X) delta_beta (Section 2.5).
      reference <- grid
      if (role == "catchability") {
        reference$X1 <- 0
      }
      x_grid      <- as.matrix(grid[, covariate_names])
      x_reference <- as.matrix(reference[, covariate_names])
      reference$X1_residual <- reference$X1 -
        predict(first_stage_x1, newdata = grid)
      reference$X2_residual <- reference$X2 -
        predict(first_stage_x2, newdata = grid)

      log_density_hat <- matrix(
        NA_real_, nrow(grid), 4,
        dimnames = list(NULL, gaussian_estimators)
      )
      log_density_hat[, "GLM"] <- as.numeric(
        predict(glm_fit, newdata = reference, type = "link")
      )
      log_density_hat[, "SGLMM"] <- as.numeric(
        predict(spatial_fit, newdata = reference, type = "link")
      )
      log_density_hat[, "RSR"] <- log_density_hat[, "SGLMM"] +
        drop((x_reference - x_grid) %*% delta_beta)
      log_density_hat[, "spatial+"] <- as.numeric(
        predict(spatial_plus_fit, newdata = reference, type = "link")
      )

      # 3.12 True and estimated log index by year: the log of the
      # area-weighted sum of the density over the cells (Eq. 19).
      log_index_true <- numeric(n_years_gaussian)
      log_index_hat  <- matrix(
        NA_real_, n_years_gaussian, 4,
        dimnames = list(NULL, gaussian_estimators)
      )
      for (i_year in seq_len(n_years_gaussian)) {
        year_rows <- which(grid_year == i_year)
        largest   <- max(grid$log_density[year_rows])
        log_index_true[i_year] <- largest + log(sum(
          grid$area[year_rows] * exp(grid$log_density[year_rows] - largest)
        ))
        for (estimator in gaussian_estimators) {
          largest <- max(log_density_hat[year_rows, estimator])
          log_index_hat[i_year, estimator] <- largest + log(sum(
            grid$area[year_rows] *
              exp(log_density_hat[year_rows, estimator] - largest)
          ))
        }
      }

      # 3.13 Coefficient errors, interval coverage, Wald tests and index
      # errors of each estimator (Section 2.6).
      change_error <- numeric(4)
      names(change_error) <- gaussian_estimators

      for (estimator in gaussian_estimators) {

        estimate       <- unname(beta_hat[[estimator]])
        se             <- unname(sqrt(diag(covariance_hat[[estimator]])))
        se_conditional <- unname(
          sqrt(diag(covariance_conditional[[estimator]]))
        )
        truth          <- unname(beta_true)
        stopifnot(is.finite(se), se > 0)

        relative_error <- (estimate - truth) / truth
        relative_error[truth == 0] <- NA_real_

        coefficient_rows[[length(coefficient_rows) + 1]] <- data.table(
          scenario                  = scenario_id,
          replicate                 = replicate_id,
          model                     = estimator,
          term                      = covariate_names,
          truth                     = truth,
          estimate                  = estimate,
          se                        = se,
          se_conditional            = se_conditional,
          covered_conditional       =
            abs(estimate - truth) <= qnorm(0.975) * se_conditional,
          rejected_zero_conditional =
            abs(estimate / se_conditional) > qnorm(0.975),
          error                     = estimate - truth,
          relative_error            = relative_error,
          covered                   =
            abs(estimate - truth) <= qnorm(0.975) * se,
          rejected_zero             = abs(estimate / se) > qnorm(0.975)
        )

        # Error of the log index after subtracting the mean over the years
        # from the estimated and from the true log index.
        centered_log_error <-
          (log_index_hat[, estimator] - mean(log_index_hat[, estimator])) -
          (log_index_true - mean(log_index_true))

        index_rows[[length(index_rows) + 1]] <- data.table(
          scenario           = scenario_id,
          replicate          = replicate_id,
          model              = estimator,
          year               = seq_len(n_years_gaussian),
          log_index          = log_index_hat[, estimator],
          true_log_index     = log_index_true,
          centered_log_error = centered_log_error
        )

        # Change error (Eq. 22): estimated minus true change in the log
        # index between the first and last year.
        change_error[estimator] <-
          (log_index_hat[n_years_gaussian, estimator] -
             log_index_hat[1, estimator]) -
          (log_index_true[n_years_gaussian] - log_index_true[1])

        performance_rows[[length(performance_rows) + 1]] <- data.table(
          scenario          = scenario_id,
          replicate         = replicate_id,
          model             = estimator,
          centered_log_rmse = sqrt(mean(centered_log_error^2)),
          first_last_error  = change_error[[estimator]],
          grid_log_rmse     = sqrt(mean(
            (log_density_hat[, estimator] - grid$log_density)^2
          ))
        )
      }

      # 3.14 First-order approximation of the change error when X1 is a
      # catchability covariate (G2, Section 2.5). The change in X1 between
      # the first and last year is the nominal trend, or the change in the
      # density-weighted mean of X1 over the grid, with the true density or
      # with the density predicted by the spatial model as weights.
      if (role == "catchability") {

        x1_mean_true    <- numeric(n_years_gaussian)
        x1_mean_spatial <- numeric(n_years_gaussian)
        for (i_year in seq_len(n_years_gaussian)) {
          year_rows <- which(grid_year == i_year)
          weight <- grid$area[year_rows] * exp(
            grid$log_density[year_rows] - max(grid$log_density[year_rows])
          )
          x1_mean_true[i_year] <-
            sum(weight * grid$X1[year_rows]) / sum(weight)
          weight <- grid$area[year_rows] * exp(
            log_density_hat[year_rows, "SGLMM"] -
              max(log_density_hat[year_rows, "SGLMM"])
          )
          x1_mean_spatial[i_year] <-
            sum(weight * grid$X1[year_rows]) / sum(weight)
        }
        x1_change_true    <- x1_mean_true[n_years_gaussian] - x1_mean_true[1]
        x1_change_spatial <- x1_mean_spatial[n_years_gaussian] -
          x1_mean_spatial[1]

        for (estimator in gaussian_estimators) {

          coefficient_error <- unname(beta_hat[[estimator]][1]) -
            beta_true[["X1"]]

          # For RSR: change error of RSR minus that of the spatial model
          # from which it is derived, and its approximation from the
          # coefficient adjustment delta_beta.
          rsr_minus_source      <- NA_real_
          rsr_shift_first_order <- NA_real_
          if (estimator == "RSR") {
            rsr_minus_source <- change_error[["RSR"]] -
              change_error[["SGLMM"]]
            rsr_shift_first_order <- -delta_beta[["X1"]] * x1_change_spatial
          }

          approximation_rows[[length(approximation_rows) + 1]] <- data.table(
            scenario                              = scenario_id,
            replicate                             = replicate_id,
            model                                 = estimator,
            X1_trend                              = x1_trend,
            true_density_weighted_X1_change       = x1_change_true,
            source_fit_density_weighted_X1_change = x1_change_spatial,
            coefficient_error                     = coefficient_error,
            nominal_first_order                   =
              -coefficient_error * x1_trend,
            weighted_first_order                  =
              -coefficient_error * x1_change_true,
            first_last_error                      = change_error[[estimator]],
            rsr_minus_source                      = rsr_minus_source,
            rsr_shift_first_order                 = rsr_shift_first_order
          )
        }
      }

      # 3.15 Diagnostics of the data set.
      # Theoretical bias (Eq. 18): coefficient of X1 in the linear projection
      # of omega on the year effects, X1 and X2 over the grid. Divided by the
      # true coefficient of X1 it is the relative bias expected when the
      # spatial effect is omitted.
      omega_projection <- lm(omega ~ year + X1 + X2, data = grid)

      # With a density covariate the RSR prediction equals that of the
      # spatial model in every cell.
      rsr_state_grid_max_abs <- NA_real_
      if (role == "density") {
        rsr_state_grid_max_abs <- max(abs(
          log_density_hat[, "RSR"] - log_density_hat[, "SGLMM"]
        ))
      }

      diagnostic_rows[[length(diagnostic_rows) + 1]] <- data.table(
        scenario                         = scenario_id,
        replicate                        = replicate_id,
        fitted_basis_k                   = k_basis,
        yearly_sd                        = sd_annual[1],
        yearly_x2_sd                     = sd_annual[2],
        grid_omitted_field_projection_X1 =
          unname(coef(omega_projection)["X1"]),
        X1_first_stage_variance_retained =
          var(observations$X1_residual) / var(observations$X1),
        rsr_fitted_max_abs               = rsr_fitted_max_abs,
        rsr_state_grid_max_abs           = rsr_state_grid_max_abs
      )
    }
  }

  # 3.16 Write the results of the run.
  fwrite(rbindlist(coefficient_rows),
         file.path(dir_run, "coefficients.csv"))
  fwrite(rbindlist(index_rows), file.path(dir_run, "indices.csv"))
  fwrite(rbindlist(performance_rows),
         file.path(dir_run, "performance.csv"))
  fwrite(rbindlist(diagnostic_rows),
         file.path(dir_run, "diagnostics.csv"))
  if (length(approximation_rows) > 0) {
    fwrite(rbindlist(approximation_rows),
           file.path(dir_run, "approximation.csv"))
  }
}
