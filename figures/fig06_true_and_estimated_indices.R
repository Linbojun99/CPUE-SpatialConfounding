# ==============================================================================
# figures/fig06_true_and_estimated_indices.R
#
# Figure 6: true and estimated relative abundance in E1 and E2 for
# ST (X1, X2). The two scenarios share catches and the fitted model. In
# prediction, X1 is retained at its local annual values in E1, where it
# affects density, and set to zero in E2, where it affects catchability.
# Within each replicate, the true index and the estimated index over the
# prediction domain are each divided by their own 2015 value. Lines show
# the annual medians over replicates of the true index (black, dashed) and
# of the estimate (blue, solid); the band shows the 10% and 90% quantiles of
# the estimate.
#
# Input:
#   output/rsr/rsr_years.csv   05_rsr.R (next to the RSR index it holds the
#                              log index of ST (X1, X2), log_est_st, and
#                              the true log index, log_true)
# Output:
#   output/figures/fig06_true_and_estimated_indices.png
#   output/figures/fig06_true_and_estimated_indices.pdf
# ==============================================================================

source("00_settings.R")
source("figures/00_figure_style.R")


# 1. Settings of the figure ----------------------------------------------------

panel_scenarios <- c("E1", "E2")
panel_titles    <- c(
  E1 = "(a) E1   Density covariate",
  E2 = "(b) E2   Catchability covariate"
)

# The two series, in the order of the legend, with their labels, colours,
# line types and point shapes. The estimate takes the colour of ST (X1, X2).
series_order   <- c("true", "estimate")
series_labels  <- c("Simulated truth", "ST estimate")
series_colours <- c(true = "#181818", estimate = estimator_colours[["M4"]])
series_lines   <- c(true = linetype_dashed, estimate = "solid")
series_shapes  <- c(true = 1, estimate = 15)

# Fill of the band between the 10% and 90% quantiles of the estimate.
colour_band <- "#BCD3E6"

figure_height <- 70


# 2. Input ---------------------------------------------------------------------

# Annual log indices over the prediction domain in each replicate. Next to
# the index of RSR, 05_rsr.R stores the index of its source model,
# ST (X1, X2), and the true index. In the rows of the index used in the
# paper, log_est_st is the bias-corrected index of ST (X1, X2) (Section 2.6)
# and log_true the true index (Eq. 9).
rsr_years <- fread(file.path(dir_rsr, "rsr_years.csv"))
indices   <- rsr_years[
  index_method == "inherited_M4_bias_ratio" &
    scenario_id %in% panel_scenarios,
  .(scenario_id, rep, year, log_est = log_est_st, log_true)
]
setorder(indices, scenario_id, rep, year)


# 3. Relative indices and their summaries over replicates ----------------------

# Each index relative to its own value in the first year, within replicates.
indices[, true_relative := exp(log_true - log_true[1]),
        by = .(scenario_id, rep)]
indices[, estimate_relative := exp(log_est - log_est[1]),
        by = .(scenario_id, rep)]

# Medians of the true and the estimated relative index, and the 10% and 90%
# quantiles of the estimate, by scenario and year.
index_summary <- indices[, .(
  true     = median(true_relative),
  estimate = median(estimate_relative),
  low      = quantile(estimate_relative, 0.1, names = FALSE),
  high     = quantile(estimate_relative, 0.9, names = FALSE)
), by = .(scenario_id, year)]
index_summary[, panel := factor(scenario_id, levels = panel_scenarios,
                                labels = panel_titles[panel_scenarios])]

# One row per scenario, year and series for the two lines.
index_series <- melt(
  index_summary,
  id.vars       = c("scenario_id", "panel", "year"),
  measure.vars  = series_order,
  variable.name = "series",
  value.name    = "value"
)


# 4. Figure --------------------------------------------------------------------

figure <- ggplot(
  index_series,
  aes(year, value, colour = series, linetype = series, shape = series)
) +
  geom_ribbon(data = index_summary,
              aes(x = year, ymin = low, ymax = high, group = panel),
              fill = colour_band, colour = NA, inherit.aes = FALSE,
              show.legend = FALSE) +
  geom_line(linewidth = 0.6) +
  geom_point(size = 1.6, stroke = 0.4, fill = "white") +
  facet_wrap(~panel, nrow = 1) +
  scale_colour_manual(values = series_colours, labels = series_labels) +
  scale_linetype_manual(values = series_lines, labels = series_labels) +
  scale_shape_manual(values = series_shapes, labels = series_labels) +
  scale_x_continuous(breaks = seq(2015, 2023, 2),
                     expand = expansion(mult = c(0.055, 0.055))) +
  scale_y_continuous(breaks = scales::breaks_width(0.2),
                     expand = expansion(mult = c(0.055, 0.055))) +
  labs(x = "Year", y = "Relative abundance index (2015 = 1)") +
  theme_figure +
  theme(
    strip.text       = element_text(size = 8.5),
    legend.spacing.x = unit(2, "mm"),
    panel.spacing    = unit(5, "mm")
  )


# 5. Save ----------------------------------------------------------------------

ggsave(
  file.path(dir_figures, "fig06_true_and_estimated_indices.png"),
  plot   = figure,
  width  = figure_width_double,
  height = figure_height,
  units  = "mm",
  dpi    = figure_dpi,
  bg     = "white",
  device = ragg::agg_png
)
ggsave(
  file.path(dir_figures, "fig06_true_and_estimated_indices.pdf"),
  plot   = figure,
  width  = figure_width_double,
  height = figure_height,
  units  = "mm",
  dpi    = figure_dpi,
  bg     = "white",
  device = cairo_pdf
)
