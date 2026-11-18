# ==============================================================================
# figures/fig07_coefficient_bias_change_error.R
#
# Figure 7: mean relative error of the X1 coefficient (horizontal axis)
# against mean change error (Eq. 22, vertical axis) for the five estimators
# that estimate the X1 coefficient.
#   (a)  A4, where X1 is a density covariate without a trend
#   (b)  E2, where X1 is a catchability covariate with a trend
# Each point is the mean of one estimator over the replicates of the
# scenario. The two panels have the same units per mm on the vertical axis
# and each covers the range of its own points. The grey line marks zero; the
# dashed line in (b) is the first-order approximation with slope
# -beta_1 s_X / 100 (Section 2.5), not a fitted regression.
#
# Input:
#   output/performance/performance_coefficient_summary.csv   06_performance.R
#   output/performance/performance_summary.csv               06_performance.R
#   output/performance/rsr_coefficient_summary.csv           05_rsr.R
#   output/performance/rsr_index_summary.csv                 05_rsr.R
# Output:
#   output/figures/fig07_coefficient_bias_change_error.png
#   output/figures/fig07_coefficient_bias_change_error.pdf
# ==============================================================================

source("00_settings.R")
source("figures/00_figure_style.R")


# 1. Settings of the figure ----------------------------------------------------

panel_scenarios  <- c("A4", "E2")
panel_estimators <- c("M1", "M3", "M4", "M5", "M7")

panel_titles <- c(
  A4 = "(a) A4 (density covariate without a trend)",
  E2 = "(b) E2 (catchability covariate with a trend)"
)

# Axis ranges. The heights of the panels are proportional to their vertical
# ranges.
x_limits <- c(24, 74)
y_limits <- list(A4 = c(-0.045, 0.07), E2 = c(-0.215, 0.02))

# Breaks of the vertical axis and their labels with two decimals. In the
# figure of the paper the label of zero in panel (b) ends with a space, which
# moves it slightly away from its tick; the space is kept here so that the
# image is the same.
y_breaks <- seq(-0.2, 0.05, 0.05)
y_labels <- list(
  A4 = c(paste0(minus_sign, c("0.20", "0.15", "0.10", "0.05")),
         "0.00", "0.05"),
  E2 = c(paste0(minus_sign, c("0.20", "0.15", "0.10", "0.05")),
         "0.00 ", "0.05")
)

# Labels next to the points, placed by hand. ST (X1) and ST (X1, X2) nearly
# coincide and share one label.
label_both_st <- "paste('ST (',X[1],') and ST (',X[1],', ',X[2],')')"
point_labels <- rbind(
  data.table(
    scenario_id = "A4",
    x           = c(34.8, 45.25, 58.2, 58.3),
    y           = c(0.0205, 0.0500, 0.0326, -0.0148),
    hjust       = c(0.5, 0.5, 0, 0),
    label       = c(estimator_labels[["M5"]], label_both_st,
                    estimator_labels[["M7"]], estimator_labels[["M1"]])
  ),
  data.table(
    scenario_id = "E2",
    x           = c(32.0, 47.4, 55.4, 54.5),
    y           = c(-0.0415, -0.1070, -0.1345, -0.1775),
    hjust       = c(0.5, 0, 0, 1),
    label       = c(estimator_labels[["M5"]], label_both_st,
                    estimator_labels[["M7"]], estimator_labels[["M1"]])
  )
)


# 2. Input ---------------------------------------------------------------------

coefficient_summary <- fread(
  file.path(dir_performance, "performance_coefficient_summary.csv")
)
performance_summary <- fread(
  file.path(dir_performance, "performance_summary.csv")
)
rsr_coefficient_summary <- fread(
  file.path(dir_performance, "rsr_coefficient_summary.csv")
)
rsr_index_summary <- fread(
  file.path(dir_performance, "rsr_index_summary.csv")
)


# 3. Means by scenario and estimator -------------------------------------------

# Mean relative error of the X1 coefficient (as a fraction) and the true
# coefficient beta_1.
coefficient_means <- rbind(
  coefficient_summary[term == "X1",
                      .(scenario_id, model_id, truth,
                        relative_error = mean_rel_err)],
  rsr_coefficient_summary[term == "X1",
                          .(scenario_id, model_id = "M7", truth,
                            relative_error = relative_bias)]
)

# Mean change error. The RSR index is the uncorrected index times the ratio
# of the bias-corrected to the uncorrected index of ST (X1, X2).
change_means <- rbind(
  performance_summary[, .(scenario_id, model_id,
                          change_error = change_err_mean)],
  rsr_index_summary[index_method == "inherited_M4_bias_ratio",
                    .(scenario_id, model_id = "M7",
                      change_error = mean_change_error)]
)

means <- merge(coefficient_means, change_means,
               by = c("scenario_id", "model_id"))
means <- means[scenario_id %in% panel_scenarios &
                 model_id %in% panel_estimators]
means[, relative_error_percent := 100 * relative_error]
means[, model := factor(model_id, levels = estimator_order)]

# Order of drawing: the square of ST (X1, X2) before the triangle of ST (X1).
means <- means[order(scenario_id,
                     match(model_id, c("M1", "M7", "M5", "M4", "M3")))]

# Slope of the first-order approximation in E2 (Section 2.5): an error of
# 1% in the X1 coefficient changes the change error by -beta_1 s_X / 100.
# s_X, the SD of X1 without the trend, is the change added to X1 between
# the first and the last year (Eq. 8). The line uses its approximate value
# sqrt(1 + 0.5^2), about 1.12 (Section 2.5). beta_1 is the same in all rows.
beta_1 <- means$truth[1]
slope_first_order <- -beta_1 * sqrt(1 + 0.5^2) / 100


# 4. Panels --------------------------------------------------------------------

panels <- list()

for (panel_scenario in panel_scenarios) {

  is_last_panel <- panel_scenario == "E2"

  panel <- ggplot(
    means[scenario_id == panel_scenario],
    aes(relative_error_percent, change_error, colour = model, fill = model,
        shape = model)
  ) +
    geom_hline(yintercept = 0, colour = colour_reference,
               linewidth = line_width_reference)

  # First-order approximation and its label.
  if (panel_scenario == "E2") {
    panel <- panel +
      geom_abline(slope = slope_first_order, intercept = 0,
                  colour = "black", linetype = linetype_dashed,
                  linewidth = line_width_reference) +
      annotate("text", x = 25.5, y = -0.130, hjust = 0, parse = TRUE,
               size = 7 / points_per_mm, family = font_family,
               label = paste0("paste('Slope = ", minus_sign, "',",
                              symbol_beta_1, ",s[X],' / 100')"))
  }

  # The horizontal axis is labelled in the lower panel only.
  x_title     <- NULL
  x_axis_text <- element_blank()
  if (is_last_panel) {
    x_title     <- expression(paste("Mean relative error in the ", X[1],
                                    " coefficient (%)"))
    x_axis_text <- element_text()
  }

  panels[[panel_scenario]] <- panel +
    geom_point(size = 2, stroke = 0.5) +
    geom_text(data = point_labels[scenario_id == panel_scenario],
              aes(x, y, label = label, hjust = hjust), inherit.aes = FALSE,
              parse = TRUE, size = 7 / points_per_mm,
              family = font_family) +
    scale_colour_manual(values = estimator_colours,
                        limits = panel_estimators) +
    scale_fill_manual(values = estimator_colours,
                      limits = panel_estimators) +
    scale_shape_manual(values = estimator_shapes,
                       limits = panel_estimators) +
    scale_x_continuous(breaks = seq(30, 70, 10)) +
    scale_y_continuous(breaks = y_breaks,
                       labels = y_labels[[panel_scenario]]) +
    coord_cartesian(xlim = x_limits, ylim = y_limits[[panel_scenario]],
                    expand = FALSE) +
    labs(x = x_title, y = "Mean change error",
         title = panel_titles[[panel_scenario]]) +
    theme_figure +
    theme(legend.position = "none", axis.text.x = x_axis_text)
}


# 5. Combine and save ----------------------------------------------------------

figure <- panels[["A4"]] / panels[["E2"]] +
  plot_layout(
    heights     = c(diff(y_limits[["A4"]]), diff(y_limits[["E2"]])),
    axis_titles = "collect"
  )

ggsave(
  file.path(dir_figures, "fig07_coefficient_bias_change_error.png"),
  plot   = figure,
  width  = figure_width_single,
  height = 112,
  units  = "mm",
  dpi    = figure_dpi,
  bg     = "white",
  device = ragg::agg_png
)
ggsave(
  file.path(dir_figures, "fig07_coefficient_bias_change_error.pdf"),
  plot   = figure,
  width  = figure_width_single,
  height = 112,
  units  = "mm",
  dpi    = figure_dpi,
  bg     = "white",
  device = cairo_pdf
)
