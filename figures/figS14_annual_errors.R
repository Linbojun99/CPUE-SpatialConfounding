# ==============================================================================
# figures/figS14_annual_errors.R
#
# Figure S14: annual errors of the index, log(estimated index / true index),
# in A4, E1 and E3. Columns show the six fitted models and RSR and rows the
# scenarios. Lines show the median and bands the 10% and 90% quantiles over
# replicates; the grey line marks equality of the estimated and the true
# index. The panels for E2 are in Figure 7.
#
# Input:
#   output/performance/performance_index_years.csv   06_performance.R
#   output/performance/performance_fits.csv          06_performance.R
#   output/rsr/rsr_years.csv                         05_rsr.R
# Output:
#   output/figures/figS14_annual_errors.png
#   output/figures/figS14_annual_errors.pdf
# ==============================================================================

source("00_settings.R")
source("figures/00_figure_style.R")


# 1. Input ---------------------------------------------------------------------

index_years <- fread(file.path(dir_performance,
                               "performance_index_years.csv"))
fits        <- fread(file.path(dir_performance, "performance_fits.csv"))
rsr_years   <- fread(file.path(dir_rsr, "rsr_years.csv"))


# 2. Annual errors in A4, E1 and E3 --------------------------------------------

row_scenarios <- c("A4", "E1", "E3")

# Fitted models: error of the log index in the fits that converged.
fits_evaluated <- fits[pred_eval == TRUE, .(scenario_id, model_id, rep)]
errors_fitted  <- merge(index_years, fits_evaluated,
                        by = c("scenario_id", "model_id", "rep"))

# RSR: the index used in the paper, the uncorrected index times the ratio
# of the bias-corrected to the uncorrected index of ST (X1, X2).
errors_rsr <- rsr_years[index_method == "inherited_M4_bias_ratio"]
errors_rsr[, err_unc := log_est - log_true]

annual_errors <- rbind(
  errors_fitted[, .(scenario_id, model_id, rep, year, error = err_unc)],
  errors_rsr[, .(scenario_id, model_id, rep, year, error = err_unc)]
)
annual_errors <- annual_errors[scenario_id %in% row_scenarios]

# Median and 10% and 90% quantiles over replicates.
error_summary <- annual_errors[, .(
  median = median(error),
  low    = quantile(error, 0.1),
  high   = quantile(error, 0.9)
), by = .(scenario_id, model_id, year)]

# One column per estimator, titled with the estimator label, and one row
# per scenario.
error_summary[, scenario := factor(scenario_id, levels = row_scenarios)]
error_summary[, model := factor(model_id, levels = estimator_order)]
error_summary[, panel := factor(estimator_labels[model_id],
                                levels = estimator_labels[estimator_order])]


# 3. Figure --------------------------------------------------------------------

figure <- ggplot(
  error_summary,
  aes(year, median, colour = model, fill = model, shape = model)
) +
  geom_hline(yintercept = 0, colour = colour_reference,
             linewidth = line_width_reference) +
  geom_ribbon(aes(ymin = low, ymax = high), alpha = 0.18, colour = NA) +
  geom_line(linewidth = 0.4) +
  geom_point(size = 1.1, stroke = 0.3) +
  facet_grid(scenario ~ panel, labeller = labeller(panel = label_parsed)) +
  scale_colour_manual(values = estimator_colours,
                      limits = estimator_order) +
  scale_fill_manual(values = estimator_colours, limits = estimator_order) +
  scale_shape_manual(values = estimator_shapes, limits = estimator_order) +
  scale_x_continuous(breaks = year_breaks, labels = year_labels_narrow) +
  scale_y_continuous(
    limits = c(-0.5, 0.5),
    breaks = c(-0.4, -0.2, 0, 0.2, 0.4),
    labels = scales::label_number(style_negative = "minus")
  ) +
  labs(x = "Year", y = "log(estimated index / true index)") +
  theme_figure +
  theme(
    axis.text.x     = suppressWarnings(
      element_text(size = 7.5, hjust = c(0, 0.5, 1))
    ),
    legend.position = "none",
    strip.text.x    = element_text(size = 7.5),
    strip.text.y    = element_text(size = 8.5, angle = 0),
    panel.spacing   = unit(1.6, "mm")
  )


# 4. Save ----------------------------------------------------------------------

ggsave(
  file.path(dir_figures, "figS14_annual_errors.png"),
  plot   = figure,
  width  = figure_width_double,
  height = 108,
  units  = "mm",
  dpi    = figure_dpi,
  bg     = "white",
  device = ragg::agg_png
)
ggsave(
  file.path(dir_figures, "figS14_annual_errors.pdf"),
  plot   = figure,
  width  = figure_width_double,
  height = 108,
  units  = "mm",
  dpi    = figure_dpi,
  bg     = "white",
  device = cairo_pdf
)
