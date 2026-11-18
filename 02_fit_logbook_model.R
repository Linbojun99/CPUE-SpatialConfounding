# ==============================================================================
# 02_fit_logbook_model.R
#
# Fits the model in Eq. 2 to the logbook catches. The catch follows a
# negative binomial distribution (Eq. 1), with year and HBF class as fixed
# effects, a random vessel effect, a spatial field and a spatiotemporal field
# that is independent among years. The script then extracts the estimates
# that give the parameter values of the operating model (Table 1). The model
# is also fitted with first-order autoregressive (AR1) annual fields for the
# comparison by AIC reported in the Supplementary material (Model
# diagnostics).
#
# Input:
#   data/logbook_sets.rds                 41,035 sets used to fit Eq. 2
#   output/meshes/mesh_logbook.rds        mesh from 01_build_meshes.R
# Output:
#   output/logbook_model/logbook_fit.rds            fitted model (Eq. 2)
#   output/logbook_model/logbook_fit_ar1.rds        fitted model, AR1 fields
#   output/logbook_model/logbook_estimates.rds      estimates used by the
#                                                   operating model
#   output/logbook_model/logbook_estimates.csv      the same, as a table
#   output/logbook_model/comparison_iid_ar1.csv     AIC and variance
#                                                   parameters of both fits
#
# Run time: about 3 minutes for Eq. 2 and about 10 minutes for the AR1 fit.
# ==============================================================================

source("00_settings.R")
library(sdmTMB)


# 1. Data and mesh -------------------------------------------------------------

logbook_sets <- readRDS(file.path(dir_data, "logbook_sets.rds"))

# Offset log(H / 1000) with H the number of hooks, year and vessel as factors.
# The HBF class is a factor with the class 16-17 as the reference level.
logbook_sets[, log_hooks := log(hooks / 1000)]
logbook_sets[, year_f := factor(year)]
logbook_sets[, vessel := factor(vessel)]

logbook_mesh <- make_mesh(
  logbook_sets,
  xy_cols = c("x", "y"),
  mesh    = readRDS(file.path(dir_meshes, "mesh_logbook.rds"))
)
cat("Vertices of the mesh:", logbook_mesh$mesh$n, "\n")


# 2. Fit the model in Eq. 2 ----------------------------------------------------

# log mu = log(H / 1000) + alpha_t + gamma_h + v_k + omega(s) + epsilon(s, t)
logbook_fit <- sdmTMB(
  bet_n ~ 0 + year_f + hbf_cl + (1 | vessel),
  data           = logbook_sets,
  mesh           = logbook_mesh,
  time           = "year",
  offset         = "log_hooks",
  family         = nbinom2(),
  spatial        = "on",
  spatiotemporal = "iid",
  silent         = TRUE
)

# Convergence checks: optimizer, gradients and Hessian.
sanity(logbook_fit)

saveRDS(logbook_fit, file.path(dir_logbook, "logbook_fit.rds"))


# 3. Estimates used by the operating model (Table 1) ---------------------------

fixed_effects  <- as.data.table(tidy(logbook_fit, "fixed"))
random_effects <- as.data.table(tidy(logbook_fit, "ran_pars"))
print(fixed_effects)
print(random_effects)

year_effects <- fixed_effects[grepl("^year_f", term), estimate]

# The logbook sets are concentrated where the fitted fields are high, whereas
# the operating model draws new fields that are independent of the set
# positions. The intercept beta_0 is therefore the mean year effect plus the
# mean of the fitted omega + epsilon at the logbook sets (note of Table 1).
fitted_at_sets     <- as.data.table(predict(logbook_fit))
mean_field_at_sets <- mean(fitted_at_sets$omega_s + fitted_at_sets$epsilon_st)

logbook_estimates <- list(
  b0        = mean(year_effects) + mean_field_at_sets,
  year_eff  = setNames(year_effects, years),
  gamma_hbf = c(
    "16-17" = 0,
    le15    = fixed_effects[term == "hbf_clle15", estimate],
    ge18    = fixed_effects[term == "hbf_clge18", estimate]
  ),
  range     = random_effects$estimate[random_effects$term == "range"],
  sigma_O   = random_effects$estimate[random_effects$term == "sigma_O"],
  sigma_E   = random_effects$estimate[random_effects$term == "sigma_E"],
  sigma_V   = random_effects$estimate[random_effects$group_name %in% "vessel"],
  phi       = random_effects$estimate[random_effects$term == "phi"]
)
saveRDS(logbook_estimates, file.path(dir_logbook, "logbook_estimates.rds"))

# The same estimates as a table, with the two parts of the intercept.
estimates_table <- data.table(
  parameter = c("beta_0", "mean of alpha_t",
                "mean of omega + epsilon at the logbook sets",
                "gamma_h (15 or fewer)", "gamma_h (18 or more)",
                "r (km)", "sigma_O", "sigma_E", "sigma_V", "phi"),
  value     = c(logbook_estimates$b0, mean(year_effects), mean_field_at_sets,
                logbook_estimates$gamma_hbf[["le15"]],
                logbook_estimates$gamma_hbf[["ge18"]],
                logbook_estimates$range, logbook_estimates$sigma_O,
                logbook_estimates$sigma_E, logbook_estimates$sigma_V,
                logbook_estimates$phi)
)
print(estimates_table)

fwrite(estimates_table, file.path(dir_logbook, "logbook_estimates.csv"))


# 4. Comparison with AR1 annual fields -----------------------------------------

# Same model and mesh, with a first-order autoregressive spatiotemporal field.
logbook_fit_ar1 <- sdmTMB(
  bet_n ~ 0 + year_f + hbf_cl + (1 | vessel),
  data           = logbook_sets,
  mesh           = logbook_mesh,
  time           = "year",
  offset         = "log_hooks",
  family         = nbinom2(),
  spatial        = "on",
  spatiotemporal = "ar1",
  silent         = TRUE
)
sanity(logbook_fit_ar1)

saveRDS(logbook_fit_ar1, file.path(dir_logbook, "logbook_fit_ar1.rds"))

# Variance parameters with 95% confidence intervals. The SD of the vessel
# effects is the only term with a group name and is labelled sigma_V. In the
# AR1 fit the term rho is the correlation of the fields between years.
variance_iid <- as.data.table(
  tidy(logbook_fit, "ran_pars", conf.int = TRUE)
)
variance_iid[!is.na(group_name), term := "sigma_V"]
variance_iid <- variance_iid[, .(
  model     = "IID",
  term,
  estimate  = signif(estimate, 3),
  conf.low  = signif(conf.low, 3),
  conf.high = signif(conf.high, 3)
)]

variance_ar1 <- as.data.table(
  tidy(logbook_fit_ar1, "ran_pars", conf.int = TRUE)
)
variance_ar1[!is.na(group_name), term := "sigma_V"]
variance_ar1 <- variance_ar1[, .(
  model     = "AR1",
  term,
  estimate  = signif(estimate, 3),
  conf.low  = signif(conf.low, 3),
  conf.high = signif(conf.high, 3)
)]

aic_table <- data.table(
  model = c("IID", "AR1"),
  AIC   = round(c(AIC(logbook_fit), AIC(logbook_fit_ar1)), 1)
)
aic_table[, dAIC := AIC - min(AIC)]
print(aic_table)

comparison <- merge(rbind(variance_iid, variance_ar1), aic_table, by = "model")
print(comparison)

fwrite(comparison, file.path(dir_logbook, "comparison_iid_ar1.csv"))
