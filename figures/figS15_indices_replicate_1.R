# ==============================================================================
# figures/figS15_indices_replicate_1.R
#
# Figure S15: true relative abundance and the estimates of the six fitted
# models and RSR in replicate 1 of four scenarios.
#   (a)  A4, X1 affects density, without a trend in X1
#   (b)  E1, X1 affects density, with a trend
#   (c)  E3, X1 affects catchability, without a trend
#   (d)  E2, X1 affects catchability, with a trend
# True and estimated indices over the prediction domain of 1,484 cells are
# each divided by their own 2015 value. The black dashed line is the true
# index.
#
# Input:
#   output/fits/<scenario>/<model>_rep0001.rds   04_fit_estimators.R (A4, E1,
#                                                E3, E2; six models)
#   data/prediction_domain_cells.csv             cells of the simulation
#                                                grid, in_domain marks the
#                                                prediction domain
#   output/rsr/rsr_years.csv                     05_rsr.R
# Output:
#   output/figures/figS15_indices_replicate_1.png
#   output/figures/figS15_indices_replicate_1.pdf
# ==============================================================================

source("00_settings.R")
source("figures/00_figure_style.R")


# 1. Settings of the figure ----------------------------------------------------

replicate_shown <- 1

panel_scenarios <- c("A4", "E1", "E3", "E2")
panel_titles <- c(
  A4 = "(a) A4 (density, no trend)",
  E1 = "(b) E1 (density, trend)",
  E3 = "(c) E3 (catchability, no trend)",
  E2 = "(d) E2 (catchability, trend)"
)


# 2. Indices of the six fitted models over the prediction domain ---------------

domain_table <- fread(file.path(dir_data, "prediction_domain_cells.csv"))
domain_cells <- domain_table[in_domain == TRUE, cell]

index_rows <- list()

for (scenario_id in panel_scenarios) {
  for (model_id in estimators$model) {

    fit <- readRDS(file.path(
      dir_fits, scenario_id,
      sprintf("%s_rep%04d.rds", model_id, replicate_shown)
    ))

    # As in 06_performance.R: the bias-corrected index over all 1,592 cells
    # is multiplied by the ratio of the uncorrected sum over the prediction
    # domain (Eq. 19) to that over all cells in the same year (Section 2.6).
    # The true index is the sum of area times true density over the
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

    index_rows[[length(index_rows) + 1]] <- index[, .(
      scenario_id = scenario_id,
      model_id    = model_id,
      year,
      log_est     = log(sum_domain * est_all / sum_all),
      log_true    = log(true_index)
    )]
  }
}

indices_fitted <- rbindlist(index_rows)


# 3. Relative indices ----------------------------------------------------------

# RSR: the uncorrected index times the ratio of the bias-corrected to the
# uncorrected index of ST (X1, X2), as calculated by 05_rsr.R.
rsr_years   <- fread(file.path(dir_rsr, "rsr_years.csv"))
indices_rsr <- rsr_years[
  index_method == "inherited_M4_bias_ratio" & rep == replicate_shown &
    scenario_id %in% panel_scenarios,
  .(scenario_id, model_id, year, log_est, log_true)
]

indices <- rbind(indices_fitted, indices_rsr)
setorder(indices, scenario_id, model_id, year)

# Each index relative to its value in the first year.
indices[, estimate_relative := exp(log_est - log_est[1]),
        by = .(scenario_id, model_id)]
indices[, true_relative := exp(log_true - log_true[1]),
        by = .(scenario_id, model_id)]

# One series per estimator and one for the true index, which is the same
# for all estimators and is taken from the rows of the first model.
index_series <- rbind(
  indices[, .(scenario_id, year, series = model_id,
              value = estimate_relative)],
  indices[model_id == estimators$model[1],
          .(scenario_id, year, series = "true", value = true_relative)]
)
series_order <- c(estimator_order, "true")
index_series[, series := factor(series, levels = series_order)]
index_series[, panel := factor(panel_titles[scenario_id],
                               levels = panel_titles)]


# 4. Figure --------------------------------------------------------------------

# Labels, colours, point shapes and line types of the eight series. The
# true index is a black dashed line without points.
series_labels <- parse(text = c(
  unname(estimator_labels[estimator_order]),
  "paste('True index',''[phantom(1)])"
))
series_colours <- c(unname(estimator_colours[estimator_order]), "black")
series_fills   <- c(unname(estimator_colours[estimator_order]), NA)
series_shapes  <- c(unname(estimator_shapes[estimator_order]), NA)
series_lines   <- c(rep("solid", length(estimator_order)), linetype_dashed)
names(series_colours) <- series_order
names(series_fills)   <- series_order
names(series_shapes)  <- series_order
names(series_lines)   <- series_order

figure <- ggplot(
  index_series,
  aes(year, value, colour = series, shape = series, linetype = series,
      fill = series)
) +
  geom_line(data = index_series[series != "true"], linewidth = 0.35) +
  geom_point(data = index_series[series != "true"], size = 1.25,
             stroke = 0.35) +
  geom_line(data = index_series[series == "true"], linewidth = 0.4) +
  facet_wrap(~panel, ncol = 2) +
  scale_colour_manual(values = series_colours, breaks = series_order,
                      labels = series_labels, limits = series_order,
                      drop = FALSE) +
  scale_fill_manual(values = series_fills, breaks = series_order,
                    labels = series_labels, limits = series_order,
                    drop = FALSE, guide = "none") +
  scale_shape_manual(values = series_shapes, breaks = series_order,
                     labels = series_labels, limits = series_order,
                     drop = FALSE) +
  scale_linetype_manual(values = series_lines, breaks = series_order,
                        labels = series_labels, limits = series_order,
                        drop = FALSE) +
  guides(
    colour   = guide_legend(nrow = 2, byrow = TRUE,
                            override.aes = list(fill = series_fills)),
    shape    = guide_legend(nrow = 2, byrow = TRUE),
    linetype = guide_legend(nrow = 2, byrow = TRUE)
  ) +
  scale_x_continuous(breaks = year_breaks,
                     expand = expansion(mult = 0.05)) +
  scale_y_continuous(breaks = seq(0.4, 1.6, 0.2)) +
  labs(x = "Year", y = "Relative abundance (2015 = 1)") +
  theme_figure +
  theme(
    legend.key.width = unit(8, "mm"),
    legend.text      = element_text(size = 7.5, hjust = 0),
    panel.spacing.x  = unit(5, "mm")
  )


# 5. Save ----------------------------------------------------------------------

ggsave(
  file.path(dir_figures, "figS15_indices_replicate_1.png"),
  plot   = figure,
  width  = figure_width_double,
  height = 125,
  units  = "mm",
  dpi    = figure_dpi,
  bg     = "white",
  device = ragg::agg_png
)
ggsave(
  file.path(dir_figures, "figS15_indices_replicate_1.pdf"),
  plot   = figure,
  width  = figure_width_double,
  height = 125,
  units  = "mm",
  dpi    = figure_dpi,
  bg     = "white",
  device = cairo_pdf
)
