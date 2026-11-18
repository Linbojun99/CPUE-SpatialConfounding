# ==============================================================================
# figures/figS10_density_rmse.R
#
# Figure S10: area-weighted RMSE of the predicted log density over the
# prediction domain of 1,484 cells, for the six fitted models and RSR.
#   (a)  cell x year combinations with at least one simulated set
#   (b)  cell x year combinations without sets
# The converged fits of all 17 scenarios are pooled. Boxes show the median
# and interquartile range and whiskers the most extreme values within 1.5
# times the interquartile range; values beyond the whiskers are not drawn.
#
# Input:
#   output/performance/performance_fits.csv   06_performance.R
#   output/rsr/rsr_fit.csv                    05_rsr.R
# Output:
#   output/figures/figS10_density_rmse.png
#   output/figures/figS10_density_rmse.pdf
# ==============================================================================

source("00_settings.R")
source("figures/00_figure_style.R")


# 1. Input ---------------------------------------------------------------------

fits     <- fread(file.path(dir_performance, "performance_fits.csv"))
rsr_fits <- fread(file.path(dir_rsr, "rsr_fit.csv"))

# Fitted models: fits that converged. RSR: the density error does not depend
# on how the index is calculated, so the rows of one index calculation are
# used.
fits     <- fits[scenario_id %in% scenario_order & pred_eval == TRUE]
rsr_fits <- rsr_fits[scenario_id %in% scenario_order &
                       index_method == "inherited_M4_bias_ratio"]


# 2. RMSE of log density by fit ------------------------------------------------

panel_titles <- c(
  sampled   = "(a) Cell \u00d7 year combinations with sets",
  unsampled = "(b) Cell \u00d7 year combinations without sets"
)

columns_used <- c("scenario_id", "model_id", "rep", "grid_rmse_sampled",
                  "grid_rmse_unsampled")
density_rmse <- rbind(fits[, ..columns_used], rsr_fits[, ..columns_used])

density_rmse <- rbind(
  density_rmse[, .(scenario_id, model_id, rep,
                   panel = panel_titles[["sampled"]],
                   value = grid_rmse_sampled)],
  density_rmse[, .(scenario_id, model_id, rep,
                   panel = panel_titles[["unsampled"]],
                   value = grid_rmse_unsampled)]
)
density_rmse[, panel := factor(panel, levels = panel_titles)]
density_rmse[, model := factor(model_id, levels = estimator_order)]

# Ends of the whiskers of each box; the vertical axis covers them with a
# margin of 0.02.
whiskers <- density_rmse[, .(
  low  = boxplot.stats(value)$stats[1],
  high = boxplot.stats(value)$stats[5]
), by = .(panel, model_id)]
y_limits <- c(min(whiskers$low) - 0.02, max(whiskers$high) + 0.02)

# Median of ST (X1, X2), redrawn in white because the dark box hides a dark
# line.
median_st <- density_rmse[model_id == "M4", .(median = median(value)),
                          by = .(panel, model)]


# 3. Figure --------------------------------------------------------------------

box_width <- 0.62

figure <- ggplot(density_rmse,
                 aes(model, value, fill = model, colour = model)) +
  geom_boxplot(outlier.shape = NA, width = box_width, linewidth = 0.3) +
  geom_errorbar(data = median_st,
                aes(ymin = median, ymax = median, y = NULL, fill = NULL,
                    colour = NULL),
                colour = "white", width = box_width, linewidth = 0.45) +
  facet_wrap(~panel, nrow = 1) +
  scale_fill_manual(values = estimator_colours, limits = estimator_order,
                    drop = FALSE) +
  scale_colour_manual(
    values = rep(colour_box_outline, length(estimator_order)),
    limits = estimator_order,
    drop   = FALSE
  ) +
  scale_x_discrete(
    labels = parse(text = unname(estimator_labels_two_lines[
      estimator_order
    ]))
  ) +
  coord_cartesian(ylim = y_limits) +
  labs(x = NULL, y = "Area-weighted RMSE of log density") +
  theme_figure +
  theme(
    legend.position = "none",
    axis.text.x     = element_text(size = 7),
    panel.spacing.x = unit(4, "mm")
  )


# 4. Save ----------------------------------------------------------------------

ggsave(
  file.path(dir_figures, "figS10_density_rmse.png"),
  plot   = figure,
  width  = figure_width_double,
  height = 92,
  units  = "mm",
  dpi    = figure_dpi,
  bg     = "white",
  device = ragg::agg_png
)
ggsave(
  file.path(dir_figures, "figS10_density_rmse.pdf"),
  plot   = figure,
  width  = figure_width_double,
  height = 92,
  units  = "mm",
  dpi    = figure_dpi,
  bg     = "white",
  device = cairo_pdf
)
