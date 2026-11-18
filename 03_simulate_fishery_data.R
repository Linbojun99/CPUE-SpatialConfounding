# ==============================================================================
# 03_simulate_fishery_data.R
#
# Generates the data sets of the fishery simulation: the 17 scenarios of
# Table 2 and the two of Table S2 (B3, B4), with 50 replicates each. The
# preliminary replicates 9001 to 9010 of A2 and D1, used for the mesh
# comparison of Table S1 and Fig. S8, are made by setting run_scenarios and
# run_replicates in section 1. For every replicate the script draws the spatial
# field omega and the spatiotemporal field epsilon, builds the covariates X1
# and X2 (Eqs. 5-8), calculates the true density in every cell and year
# (Eq. 3) and the true index, draws 10,000 of the 39,440 sets retained for
# resampling, and simulates their catches (Eqs. 1 and 4).
#
# Random numbers: each group of random draws of replicate r starts from its
# own seed, base_seed + 100 r + j, where j is 1 for omega, 2 for epsilon,
# 3 for Z, 4 for delta, 5 for W, 6 for the sets, 7 for the vessel effects and
# 8 for the catches. The seeds do not depend on the scenario, so all
# scenarios share omega, epsilon, the sets and the vessel effects of a
# replicate.
#
# Input:
#   data/logbook_sets_resampling.rds             39,440 sets retained for
#                                                resampling
#   data/simulation_grid.rds                     1,592 cells of the simulation
#                                                grid
#   output/meshes/mesh_simulation.rds            from 01_build_meshes.R
#   output/logbook_model/logbook_estimates.rds   from 02_fit_logbook_model.R
# Output:
#   output/simulation/operating_model_parameters.rds   parameter values of
#   output/simulation/operating_model_parameters.csv   the operating model
#   output/simulation/<scenario>/rep<NNNN>.rds   one file per replicate, a
#     list with the elements
#       sets      simulated sets: year, position, cell, vessel, number of
#                 hooks, HBF class, X1 and X2 at the set, true log expected
#                 catch (log_mu_true) and simulated catch (bet_n)
#       grid      one row per cell and year: X1, X2, X1 without the trend
#                 (X1_0), omega, epsilon (eps), true log density (log_d),
#                 cell area and whether the cell had a set in that year
#       truth     true index by year, summed over all 1,592 cells of the
#                 simulation grid. The true index over the prediction domain
#                 of 1,484 cells (Eq. 9) is calculated from grid$log_d when
#                 performance is summarized.
#       gen_diag  checks of the simulated data: means, SDs and correlations
#                 of the covariates and fields over the grid and at the sets
#       params    scenario, replicate, seed and parameter values
#
# Run time: about 5 minutes for all 19 scenarios and 50 replicates.
# ==============================================================================

source("00_settings.R")


# 1. Settings for this run -----------------------------------------------------

# Scenarios and replicates to simulate. The defaults are those of the paper:
# all 19 scenarios and replicates 1 to 50. For the mesh comparison (Table S1,
# Fig. S8) also run c("A2", "D1") with replicates 9001:9010.
run_scenarios  <- scenarios$scenario
run_replicates <- 1:n_replicates


# 2. Fixed settings of the operating model -------------------------------------

# Correlation between X1 and omega at which the SD of the spatial field in
# Eq. 2 (sigma_O) is divided between omega and X1 (Section 2.2).
reference_rho_omega <- 0.6

# SD of the year-to-year part of X1, delta, within each year (Eq. 6). The
# spatiotemporal component of W (Eq. 7) is multiplied by the same value.
sd_delta <- 0.5


# 3. Input ---------------------------------------------------------------------

resampling_sets <- readRDS(file.path(dir_data, "logbook_sets_resampling.rds"))
simulation_grid <- readRDS(file.path(dir_data, "simulation_grid.rds"))
simulation_mesh <- readRDS(file.path(dir_meshes, "mesh_simulation.rds"))
logbook_estimates <- readRDS(file.path(dir_logbook, "logbook_estimates.rds"))

# Mesh and interpolation matrices: A_grid gives values at the cell centres
# and A_sets at the positions of the sets retained for resampling, from
# values at the mesh vertices.
mesh       <- simulation_mesh$mesh
A_grid     <- simulation_mesh$A_grid
A_sets     <- simulation_mesh$A_sets
n_vertices <- mesh$n
n_cells    <- nrow(simulation_grid)
stopifnot(nrow(A_grid) == n_cells, nrow(A_sets) == nrow(resampling_sets))


# 4. Parameter values of the operating model (Table 1) -------------------------

# Estimates from the model fitted to the logbook data (Eq. 2).
beta_0      <- logbook_estimates$b0
field_range <- logbook_estimates$range
sigma_E     <- logbook_estimates$sigma_E
sigma_V     <- logbook_estimates$sigma_V
phi         <- logbook_estimates$phi
gamma_hbf   <- logbook_estimates$gamma_hbf

# SD of omega: sigma_O divided between omega and X1 at the reference
# correlation, sigma_omega = sigma_O / sqrt(3.2). Every non-zero true
# coefficient equals sigma_omega (the base effect size).
sigma_omega <- logbook_estimates$sigma_O / sqrt(2 + 2 * reference_rho_omega)
base_effect <- sigma_omega

operating_model <- list(
  b0          = beta_0,
  range       = field_range,
  sigma_omega = sigma_omega,
  beta1       = base_effect,
  sigma_E     = sigma_E,
  sigma_V     = sigma_V,
  phi         = phi,
  gamma_hbf   = gamma_hbf
)
saveRDS(
  operating_model,
  file.path(dir_simulation, "operating_model_parameters.rds")
)

# The same values as a table, written with 17 significant digits.
fwrite(
  data.table(
    parameter = c("b0", "range", "sigma_omega", "beta1", "sigma_E",
                  "sigma_V", "phi",
                  paste0("gamma_hbf[", names(gamma_hbf), "]")),
    value     = format(c(beta_0, field_range, sigma_omega, base_effect,
                         sigma_E, sigma_V, phi, gamma_hbf), digits = 17)
  ),
  file.path(dir_simulation, "operating_model_parameters.csv")
)


# 5. Number of sets drawn in each year -----------------------------------------

# The 10,000 sets of a replicate are allocated among years in proportion to
# the annual counts of the sets retained for resampling (Section 2.2).
sets_per_year <- resampling_sets[, .N, by = year][match(years, year), N]
n_sets_year   <- round(n_sets * sets_per_year / sum(sets_per_year))

# Vessels, in the order in which their effects are drawn.
vessels <- sort(unique(resampling_sets$vessel))


# 6. Precision matrix of omega and epsilon -------------------------------------

# A Gaussian field with a Matern covariance (smoothness 1) is represented by
# its values at the mesh vertices, which have the sparse precision matrix Q
# of the stochastic partial differential equation approach. With the
# Cholesky factorization Q = P' L L' P and z a vector of independent standard
# normal values, x = P' solve(L', z) has precision Q. The fields are drawn
# with unit variance and scaled afterwards.
precision_field <- fm_matern_precision(
  mesh,
  alpha = 2,
  rho   = field_range,
  sigma = 1
)
precision_field <- forceSymmetric(as(precision_field, "CsparseMatrix"))
chol_field      <- Cholesky(precision_field, LDL = FALSE, perm = TRUE)


# 7. Simulate the data sets ----------------------------------------------------

# Pairs of variables whose correlation is recorded in the checks.
pair_first  <- c("X1", "X1", "X2", "X1", "X2", "X1_0", "X1_0")
pair_second <- c("X2", "omega", "omega", "eps", "eps", "X2", "omega")

for (scenario_id in run_scenarios) {

  design <- scenarios[scenario == scenario_id]
  dir.create(file.path(dir_simulation, scenario_id), showWarnings = FALSE)

  # True coefficients of X1 and X2 on density (d) and catchability (q).
  true_coefficients <- c(
    b1d = design$beta1_d,
    b2d = design$beta2_d,
    b1q = design$beta1_q,
    b2q = design$beta2_q
  ) * base_effect

  # Precision matrix of the fields Z, delta and W, which have the range
  # c_r r (c_r is the range ratio of the scenario).
  precision_covariate <- fm_matern_precision(
    mesh,
    alpha = 2,
    rho   = field_range * design$range_ratio,
    sigma = 1
  )
  precision_covariate <- forceSymmetric(
    as(precision_covariate, "CsparseMatrix")
  )
  chol_covariate <- Cholesky(precision_covariate, LDL = FALSE, perm = TRUE)

  for (replicate_id in run_replicates) {

    seed_replicate <- base_seed + replicate_id * 100

    # Spatial field omega (Eq. 3) at the mesh vertices, and omega centred and
    # scaled to unit variance over the grid.
    set.seed(seed_replicate + 1)
    white_noise <- matrix(rnorm(n_vertices), n_vertices, 1)
    half_solved <- solve(chol_field, white_noise, system = "Lt")
    omega_unit  <- as.matrix(solve(chol_field, half_solved, system = "Pt"))
    omega_unit  <- omega_unit[, 1]
    omega       <- omega_unit * sigma_omega

    omega_unit_cells <- as.vector(A_grid %*% omega_unit)
    omega_std        <- (omega_unit - mean(omega_unit_cells)) /
      sd(omega_unit_cells)

    # Spatiotemporal field epsilon (Eq. 3), independent among years: one
    # column per year.
    set.seed(seed_replicate + 2)
    white_noise <- matrix(rnorm(n_vertices * n_years), n_vertices, n_years)
    half_solved <- solve(chol_field, white_noise, system = "Lt")
    epsilon     <- as.matrix(solve(chol_field, half_solved, system = "Pt"))
    epsilon     <- epsilon * sigma_E

    # Z (Eq. 5): a spatial field of range c_r r is regressed on omega over
    # the grid, and the residuals are centred and scaled to unit variance
    # over the grid, so that Z is uncorrelated with omega over the grid.
    set.seed(seed_replicate + 3)
    white_noise <- matrix(rnorm(n_vertices), n_vertices, 1)
    half_solved <- solve(chol_covariate, white_noise, system = "Lt")
    z_raw       <- as.matrix(solve(chol_covariate, half_solved, system = "Pt"))
    z_raw       <- z_raw[, 1]

    z_raw_cells      <- as.vector(A_grid %*% z_raw)
    omega_std_cells  <- as.vector(A_grid %*% omega_std)
    regression       <- coef(lm(z_raw_cells ~ omega_std_cells))
    z_residual       <- z_raw - regression[[1]] - regression[[2]] * omega_std
    z_residual_cells <- z_raw_cells - regression[[1]] -
      regression[[2]] * omega_std_cells
    z_std            <- (z_residual - mean(z_residual_cells)) /
      sd(z_residual_cells)

    # Time-invariant part of X1 (Eq. 5): unit variance and correlation
    # rho_omega with omega over the grid.
    x1_spatial <- design$rho_omega * omega_std +
      sqrt(1 - design$rho_omega^2) * z_std

    # delta (Eq. 6): a spatiotemporal field of range c_r r, independent among
    # years. Each year is centred and scaled over the grid to an SD of 0.5.
    set.seed(seed_replicate + 4)
    white_noise <- matrix(rnorm(n_vertices * n_years), n_vertices, n_years)
    half_solved <- solve(chol_covariate, white_noise, system = "Lt")
    delta       <- as.matrix(solve(chol_covariate, half_solved, system = "Pt"))
    delta_cells <- as.matrix(A_grid %*% delta)
    for (i_year in seq_len(n_years)) {
      delta[, i_year] <- (delta[, i_year] - mean(delta_cells[, i_year])) /
        sd(delta_cells[, i_year]) * sd_delta
    }

    # X1 without a trend (Eq. 6), one row per vertex and one column per
    # year, its SD over all cells and years (s_X), and the same covariate
    # centred and scaled to unit variance over all cells and years.
    X1_0       <- x1_spatial + delta
    X1_0_cells <- as.vector(as.matrix(A_grid %*% X1_0))
    s_X        <- sd(X1_0_cells)
    X1_0_std   <- (X1_0 - mean(X1_0_cells)) / s_X

    # W (Eq. 7): a spatial component plus a spatiotemporal component (both
    # of range c_r r) is regressed on X1 without a trend over all cells and
    # years, and the residuals are centred and scaled to unit variance.
    set.seed(seed_replicate + 5)
    white_noise <- matrix(rnorm(n_vertices), n_vertices, 1)
    half_solved <- solve(chol_covariate, white_noise, system = "Lt")
    w_spatial   <- as.matrix(solve(chol_covariate, half_solved, system = "Pt"))
    w_spatial   <- w_spatial[, 1]
    white_noise <- matrix(rnorm(n_vertices * n_years), n_vertices, n_years)
    half_solved <- solve(chol_covariate, white_noise, system = "Lt")
    w_annual    <- as.matrix(solve(chol_covariate, half_solved, system = "Pt"))
    w_raw       <- w_spatial + w_annual * sd_delta

    w_raw_cells      <- as.vector(as.matrix(A_grid %*% w_raw))
    X1_0_std_cells   <- as.vector(as.matrix(A_grid %*% X1_0_std))
    regression       <- coef(lm(w_raw_cells ~ X1_0_std_cells))
    w_residual       <- w_raw - regression[[1]] - regression[[2]] * X1_0_std
    w_residual_cells <- w_raw_cells - regression[[1]] -
      regression[[2]] * X1_0_std_cells
    w_std            <- (w_residual - mean(w_residual_cells)) /
      sd(w_residual_cells)

    # X2 (Eq. 7): correlation rho_12 with X1 without a trend.
    X2 <- design$rho_12 * X1_0_std + sqrt(1 - design$rho_12^2) * w_std

    # X1 (Eq. 8): a shift g_t common to all cells is added in the scenarios
    # with a trend (E1, E2), so that X1 increases by s_X between the first
    # and last year. The shift is zero in the other scenarios.
    x1_shift <- design$x1_trend * s_X *
      ((seq_len(n_years) - 1) / (n_years - 1) - 0.5)
    X1       <- X1_0 + rep(x1_shift, each = n_vertices)

    # True log density at the mesh vertices (Eq. 3).
    log_density <- matrix(0, n_vertices, n_years)
    for (i_year in seq_len(n_years)) {
      log_density[, i_year] <- beta_0 + log(abundance_trend[i_year]) + omega +
        epsilon[, i_year] + true_coefficients[["b1d"]] * X1[, i_year] +
        true_coefficients[["b2d"]] * X2[, i_year]
    }

    # Covariates, fields and true log density at the cell centres, one row
    # per cell and year.
    cell_values <- data.table(
      year     = rep(years, each = n_cells),
      cell     = rep(simulation_grid$cell, n_years),
      X1       = as.vector(as.matrix(A_grid %*% X1)),
      X2       = as.vector(as.matrix(A_grid %*% X2)),
      X1_0     = X1_0_cells,
      omega    = rep(as.vector(A_grid %*% omega), n_years),
      eps      = as.vector(as.matrix(A_grid %*% epsilon)),
      log_d    = as.vector(as.matrix(A_grid %*% log_density)),
      area_km2 = rep(simulation_grid$area_km2, n_years)
    )

    # True index: density summed over the cells, weighted by cell area, and
    # the centred log index.
    truth <- cell_values[, .(true_index = sum(area_km2 * exp(log_d))),
                         by = year]
    truth[, log_true := log(true_index)]
    truth[, l_true := log_true - mean(log_true)]

    # Sets: drawn without replacement within each year from the sets
    # retained for resampling.
    set.seed(seed_replicate + 6)
    sampled_rows <- integer(0)
    sampled_year <- integer(0)
    for (i_year in seq_len(n_years)) {
      rows_year    <- which(resampling_sets$year == years[i_year])
      sampled_rows <- c(sampled_rows, sample(rows_year, n_sets_year[i_year]))
      sampled_year <- c(sampled_year, rep(years[i_year], n_sets_year[i_year]))
    }
    sampled_sets <- resampling_sets[sampled_rows]

    # Covariates, fields and true log density at the set positions, in the
    # year of each set.
    A_sampled        <- A_sets[sampled_rows, , drop = FALSE]
    set_and_year     <- cbind(seq_along(sampled_rows),
                              match(sampled_year, years))
    x1_sets          <- as.matrix(A_sampled %*% X1)[set_and_year]
    x2_sets          <- as.matrix(A_sampled %*% X2)[set_and_year]
    x1_0_sets        <- as.matrix(A_sampled %*% X1_0)[set_and_year]
    log_density_sets <- as.matrix(A_sampled %*% log_density)[set_and_year]
    omega_sets       <- as.vector(A_sampled %*% omega)
    epsilon_sets     <- as.matrix(A_sampled %*% epsilon)[set_and_year]

    # Vessel effects, drawn anew in each replicate.
    set.seed(seed_replicate + 7)
    vessel_effects <- setNames(rnorm(length(vessels), 0, sigma_V), vessels)

    # Expected catch (Eq. 4).
    log_mu <- log(sampled_sets$hooks / 1000) +
      gamma_hbf[as.character(sampled_sets$hbf_cl)] +
      vessel_effects[sampled_sets$vessel] +
      log_density_sets +
      true_coefficients[["b1q"]] * x1_sets +
      true_coefficients[["b2q"]] * x2_sets
    log_mu <- unname(log_mu)

    # Catch in numbers (Eq. 1).
    set.seed(seed_replicate + 8)
    catch <- as.integer(rnbinom(length(log_mu), size = phi, mu = exp(log_mu)))

    simulated_sets <- data.table(
      year        = sampled_year,
      month       = month(sampled_sets$date),
      set_src     = sampled_sets$set_id,
      vessel      = sampled_sets$vessel,
      lon         = sampled_sets$lon,
      lat         = sampled_sets$lat,
      x           = sampled_sets$x,
      y           = sampled_sets$y,
      cell        = sampled_sets$cell,
      hooks       = sampled_sets$hooks,
      hbf         = sampled_sets$hbf,
      hbf_cl      = sampled_sets$hbf_cl,
      X1          = x1_sets,
      X2          = x2_sets,
      log_mu_true = log_mu,
      bet_n       = catch
    )

    # Mark the cells with at least one set in the year.
    cell_values[, sampled := paste(year, cell) %in%
                  paste(simulated_sets$year, simulated_sets$cell)]
    setcolorder(cell_values, c("year", "cell", "sampled"))

    # Checks of the simulated data (Supplementary Table S6): weighted means
    # and SDs of the covariates and weighted correlations among covariates
    # and fields, over the grid with equal weights, over the grid weighted
    # by cell area, and at the simulated sets.
    set_values <- data.table(
      X1    = x1_sets,
      X2    = x2_sets,
      X1_0  = x1_0_sets,
      omega = omega_sets,
      eps   = epsilon_sets
    )
    checks <- list()
    for (weighting in c("grid_equal", "grid_area", "sets")) {
      if (weighting == "grid_equal") {
        values  <- cell_values
        weights <- rep(1, nrow(cell_values))
      }
      if (weighting == "grid_area") {
        values  <- cell_values
        weights <- cell_values$area_km2
      }
      if (weighting == "sets") {
        values  <- set_values
        weights <- rep(1, nrow(set_values))
      }
      for (variable in c("X1", "X2", "X1_0")) {
        x      <- values[[variable]]
        mean_x <- sum(weights * x) / sum(weights)
        sd_x   <- sqrt(sum(weights * (x - mean_x)^2) / sum(weights))
        checks[[length(checks) + 1]] <- data.table(
          version = weighting,
          stat    = c("mean", "sd"),
          var     = variable,
          value   = c(mean_x, sd_x)
        )
      }
      for (p in seq_along(pair_first)) {
        x          <- values[[pair_first[p]]]
        y          <- values[[pair_second[p]]]
        mean_x     <- sum(weights * x) / sum(weights)
        mean_y     <- sum(weights * y) / sum(weights)
        sd_x       <- sqrt(sum(weights * (x - mean_x)^2) / sum(weights))
        sd_y       <- sqrt(sum(weights * (y - mean_y)^2) / sum(weights))
        covariance <- sum(weights * ((x - mean_x) * (y - mean_y))) /
          sum(weights)
        checks[[length(checks) + 1]] <- data.table(
          version = weighting,
          stat    = "cor",
          var     = paste(pair_first[p], pair_second[p], sep = "-"),
          value   = covariance / (sd_x * sd_y)
        )
      }
    }

    # Further checks: correlation between the time-invariant part of X1 and
    # omega over the grid, variances of the covariate term and of the linear
    # predictors, mean expected catch, share of sets with zero catch and the
    # shift added to X1 between the first and last year.
    covariate_term_cells <-
      (true_coefficients[["b1d"]] + true_coefficients[["b1q"]]) *
      cell_values$X1 +
      (true_coefficients[["b2d"]] + true_coefficients[["b2q"]]) *
      cell_values$X2
    covariate_term_sets <-
      (true_coefficients[["b1d"]] + true_coefficients[["b1q"]]) * x1_sets +
      (true_coefficients[["b2d"]] + true_coefficients[["b2q"]]) * x2_sets
    x1_spatial_cells <- as.vector(A_grid %*% x1_spatial)
    checks[[length(checks) + 1]] <- data.table(
      version = "grid_equal",
      stat    = c("cor", "var", "var", "var", "var", "mean", "mean", "trend"),
      var     = c("X1_spatial-omega", "covariate_term_grid",
                  "covariate_term_sets", "log_d_grid", "log_mu_sets",
                  "mu_sets", "zero_share", "X1_first_last_added"),
      value   = c(cor(x1_spatial_cells, omega_std_cells),
                  var(covariate_term_cells), var(covariate_term_sets),
                  var(cell_values$log_d), var(log_mu), mean(exp(log_mu)),
                  mean(catch == 0), x1_shift[n_years] - x1_shift[1])
    )
    gen_diag <- rbindlist(checks)

    # Scenario, replicate, seed and parameter values.
    params <- list(
      scenario    = as.list(design),
      rep         = replicate_id,
      seed0       = seed_replicate,
      b           = base_effect,
      coef        = true_coefficients,
      b0          = beta_0,
      range       = field_range,
      sigma_omega = sigma_omega,
      sigma_E     = sigma_E,
      sigma_V     = sigma_V,
      phi         = phi,
      gamma_hbf   = gamma_hbf,
      tau         = abundance_trend,
      s_X         = s_X,
      m_t         = x1_shift,
      n_sets_year = simulated_sets[, .N, by = year][order(year)],
      n_sets      = nrow(simulated_sets)
    )

    saveRDS(
      list(sets = simulated_sets, grid = cell_values, truth = truth,
           gen_diag = gen_diag, params = params),
      file.path(dir_simulation, scenario_id,
                sprintf("rep%04d.rds", replicate_id))
    )
  }
  cat("Scenario", scenario_id, "finished\n")
}
