# ==============================================================================
# figures/fig04_coefficient_error_spatial_scale.R
#
# Figure 4: relative error in the estimated X1 coefficient at three spatial
# scales, for the five estimators that estimate this coefficient. The range
# of the unconfounded component of X1 is 0.3, 1 or 3 times the range of
# omega (c_r) in D1, A4 and D2, respectively; all three scenarios have
# rho_omega = 0.6 and rho_12 = 0.8. One panel per scenario. Boxes show the
# median and interquartile range over replicates and whiskers the most
# extreme values within 1.5 times the interquartile range; values beyond the
# whiskers are not drawn. White circles show the mean, the grey line zero
# error and the dashed line the theoretical bias (Eq. 18).
#
# Input:
#   output/performance/performance_coefficients.csv   06_performance.R
#   output/performance/performance_fits.csv           06_performance.R
#   output/rsr/rsr_coefficients.csv                   05_rsr.R
# Output:
#   output/figures/fig04_coefficient_error_spatial_scale.png
#   output/figures/fig04_coefficient_error_spatial_scale.pdf
# ==============================================================================

source("00_settings.R")
source("figures/00_figure_style.R")


# 1. Settings of the figure ----------------------------------------------------

# Scenarios, one panel each, from the finest to the coarsest scale of the
# unconfounded component of X1, and the titles of the panels.
panel_scenarios <- c("D1", "A4", "D2")
panel_titles    <- c(
  D1 = "(a) D1   Finer scale",
  A4 = "(b) A4   Equal scale",
  D2 = "(c) D2   Coarser scale"
)

# Estimators from left to right in a panel. ST (none) and ST (X2) do not
# estimate the X1 coefficient.
panel_estimators <- c("M1", "M3", "M4", "M5", "M7")

# Width of the boxes and of the white median line of ST (X1, X2).
box_width    <- 0.6
median_width <- 0.43

figure_height <- 72


# 2. Input ---------------------------------------------------------------------

coefficient_table <- fread(
  file.path(dir_performance, "performance_coefficients.csv")
)
fits <- fread(file.path(dir_performance, "performance_fits.csv"))
rsr_coefficients <- fread(file.path(dir_rsr, "rsr_coefficients.csv"))


# 3. Relative errors of the X1 coefficient -------------------------------------

# Fitted models: fits that converged with finite standard errors. RSR: all
# source fits for which it was calculated.
fits_evaluated <- fits[coef_eval == TRUE, .(scenario_id, model_id, rep)]
errors_fitted  <- merge(
  coefficient_table[term == "X1" & truth != 0],
  fits_evaluated,
  by = c("scenario_id", "model_id", "rep")
)

relative_errors <- rbind(
  errors_fitted[, .(scenario_id, model_id, rep, rel_err)],
  rsr_coefficients[term == "X1" & truth != 0,
                   .(scenario_id, model_id, rep, rel_err)]
)
relative_errors <- relative_errors[
  scenario_id %in% panel_scenarios & model_id %in% panel_estimators
]

# Relative error in percent.
relative_errors[, value := 100 * rel_err]


# 4. Box statistics ------------------------------------------------------------

# The boxes are drawn from their statistics, so that the vertical axis
# covers the whiskers and not the values beyond them. Mean, median and
# quartiles of each scenario and estimator.
box_statistics <- relative_errors[, .(
  n      = .N,
  mean   = 100 * mean(rel_err),
  q1     = quantile(value, 0.25, names = FALSE),
  median = median(value),
  q3     = quantile(value, 0.75, names = FALSE)
), keyby = .(scenario_id, model_id)]

# Whiskers: the most extreme values within 1.5 times the interquartile
# range from the quartiles.
relative_errors <- merge(relative_errors, box_statistics,
                         by = c("scenario_id", "model_id"))
relative_errors[, iqr := q3 - q1]
relative_errors[, within_whiskers :=
                  value >= q1 - 1.5 * iqr & value <= q3 + 1.5 * iqr]
whiskers <- relative_errors[, .(
  whisker_low  = min(q1[1], value[within_whiskers]),
  whisker_high = max(q3[1], value[within_whiskers])
), keyby = .(scenario_id, model_id)]
box_statistics <- merge(box_statistics, whiskers,
                        by = c("scenario_id", "model_id"))

box_statistics[, panel := factor(scenario_id, levels = panel_scenarios,
                                 labels = panel_titles[panel_scenarios])]
box_statistics[, model := factor(model_id, levels = panel_estimators)]

# Theoretical bias of the X1 coefficient when the spatial field is omitted,
# 0.8 rho_omega as a percentage (Eq. 18). rho_omega is the same in the three
# scenarios.
scenario_design  <- scenarios[match(panel_scenarios, scenario)]
theoretical_bias <- 100 * 0.8 * unique(scenario_design$rho_omega)


# 5. Figure --------------------------------------------------------------------

# The median line of ST (X1, X2) is redrawn in white, because the dark box
# hides a dark line.
figure <- ggplot(box_statistics, aes(model, fill = model)) +
  geom_hline(yintercept = 0, colour = colour_reference,
             linewidth = line_width_reference) +
  geom_hline(yintercept = theoretical_bias, linewidth = 0.3,
             linetype = linetype_dashed) +
  geom_boxplot(aes(ymin = whisker_low, lower = q1, middle = median,
                   upper = q3, ymax = whisker_high),
               stat = "identity", width = box_width, colour = "#343434",
               linewidth = 0.32) +
  geom_errorbar(data = box_statistics[model_id == "M4"],
                aes(ymin = median, ymax = median), width = median_width,
                colour = "white", linewidth = 0.45) +
  geom_point(aes(y = mean), colour = "black", fill = "white", shape = 21,
             size = 1.55, stroke = 0.35) +
  facet_wrap(~panel, nrow = 1) +
  scale_fill_manual(
    values = estimator_colours[panel_estimators],
    breaks = panel_estimators,
    labels = parse(text = unname(estimator_labels[panel_estimators]))
  ) +
  scale_y_continuous(
    breaks = seq(-50, 200, 50),
    labels = scales::label_number(style_negative = "minus")
  ) +
  labs(
    x = NULL,
    y = expression(paste("Relative error in the estimated ", X[1],
                         " coefficient (%)"))
  ) +
  theme_figure +
  theme(
    axis.text.x   = element_blank(),
    axis.ticks.x  = element_blank(),
    strip.text    = element_text(size = 8.3),
    panel.spacing = unit(3, "mm")
  )


# 6. Save ----------------------------------------------------------------------

ggsave(
  file.path(dir_figures, "fig04_coefficient_error_spatial_scale.png"),
  plot   = figure,
  width  = figure_width_double,
  height = figure_height,
  units  = "mm",
  dpi    = figure_dpi,
  bg     = "white",
  device = ragg::agg_png
)
ggsave(
  file.path(dir_figures, "fig04_coefficient_error_spatial_scale.pdf"),
  plot   = figure,
  width  = figure_width_double,
  height = figure_height,
  units  = "mm",
  dpi    = figure_dpi,
  bg     = "white",
  device = cairo_pdf
)
