# ==============================================================================
# figures/fig03_coefficient_error_groups_AB.R
#
# Figure 3: relative error in the estimated X1 coefficient in the eight
# scenarios of groups A and B, for the five estimators that estimate this
# coefficient. One panel per scenario, with the scenario label and its values
# of rho_omega, rho_12 and c_r at the top left. Boxes show the median and
# interquartile range over replicates and whiskers the most extreme values
# within 1.5 times the interquartile range; values beyond the whiskers are
# not drawn. White circles show the mean, the grey line zero error and the
# dashed line the theoretical bias (Eq. 18).
#
# Input:
#   output/performance/performance_coefficients.csv   06_performance.R
#   output/performance/performance_fits.csv           06_performance.R
#   output/rsr/rsr_coefficients.csv                   05_rsr.R
# Output:
#   output/figures/fig03_coefficient_error_groups_AB.png
#   output/figures/fig03_coefficient_error_groups_AB.pdf
# ==============================================================================

source("00_settings.R")
source("figures/00_figure_style.R")


# 1. Settings of the figure ----------------------------------------------------

# Scenarios, one panel each, and estimators from left to right in a panel.
# ST (none) and ST (X2) do not estimate the X1 coefficient.
panel_scenarios  <- c("A1", "A2", "A3", "A4", "A5", "A6", "B1", "B2")
panel_estimators <- c("M1", "M3", "M4", "M5", "M7")

# Range of the vertical axis, height of the scenario label and width of the
# boxes.
y_limits   <- c(-150, 300)
y_label    <- 262
box_width  <- 0.66

figure_height <- 95


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
relative_errors[, panel := factor(scenario_id, levels = panel_scenarios)]
relative_errors[, model := factor(model_id, levels = panel_estimators)]

# Mean and median of each scenario and estimator.
error_summary <- relative_errors[, .(
  mean   = mean(value),
  median = median(value)
), by = .(scenario_id, model_id, panel, model)]

# Number of values beyond the whiskers, which are not drawn (reported in the
# caption).
relative_errors[, beyond_whiskers :=
                  value < quantile(value, 0.25) - 1.5 * IQR(value) |
                  value > quantile(value, 0.75) + 1.5 * IQR(value),
                by = .(scenario_id, model_id)]
cat("Values:", nrow(relative_errors), " beyond the whiskers:",
    sum(relative_errors$beyond_whiskers), "\n")


# 4. Scenario labels and theoretical bias --------------------------------------

scenario_design <- scenarios[match(panel_scenarios, scenario)]

# Label of each panel: the scenario in bold, then rho_omega, rho_12 and c_r.
panel_labels <- data.table(
  panel = factor(panel_scenarios, levels = panel_scenarios),
  label = sprintf(
    "paste(bold('%s'),'  ',%s,' = %s, ',%s,' = %s, ',c[r],' = %s')",
    scenario_design$scenario,
    symbol_rho_omega, scenario_design$rho_omega,
    symbol_rho_12, scenario_design$rho_12,
    scenario_design$range_ratio
  )
)

# Theoretical bias of the X1 coefficient when the spatial field is omitted,
# 0.8 rho_omega as a percentage (Eq. 18); no line where rho_omega is zero.
theoretical_bias <- data.table(
  panel = factor(panel_scenarios, levels = panel_scenarios),
  value = 100 * 0.8 * scenario_design$rho_omega
)
theoretical_bias <- theoretical_bias[value > 0]


# 5. Panels --------------------------------------------------------------------

# The median line of ST (X1, X2) is redrawn in white, because the dark box
# hides a dark line.
panels <- ggplot(relative_errors, aes(model, value, fill = model)) +
  geom_hline(yintercept = 0, colour = colour_reference,
             linewidth = line_width_reference) +
  geom_hline(data = theoretical_bias, aes(yintercept = value),
             colour = "black", linewidth = 0.35,
             linetype = linetype_dashed) +
  geom_boxplot(outlier.shape = NA, width = box_width, linewidth = 0.3,
               median.linewidth = 0.42, colour = colour_box_outline) +
  geom_errorbar(data = error_summary[model_id == "M4"],
                aes(y = NULL, ymin = median, ymax = median),
                colour = "white", width = box_width * 0.55,
                linewidth = 0.42) +
  geom_point(data = error_summary, aes(y = mean), shape = 21,
             fill = "white", colour = "black", size = 1.15, stroke = 0.35) +
  geom_text(data = panel_labels, aes(x = 0.565, y = y_label, label = label),
            parse = TRUE, hjust = 0, inherit.aes = FALSE,
            size = 7 / points_per_mm, family = font_family,
            colour = "black") +
  scale_fill_manual(values = estimator_colours[panel_estimators],
                    limits = panel_estimators) +
  scale_y_continuous(
    breaks = seq(-100, 200, 100),
    labels = scales::label_number(style_negative = "minus")
  ) +
  coord_cartesian(xlim = c(0.42, 5.58), ylim = y_limits, expand = FALSE) +
  facet_wrap(~panel, ncol = 4) +
  labs(
    x = NULL,
    y = expression(paste("Relative error in the estimated ", X[1],
                         " coefficient (%)"))
  ) +
  theme_figure +
  theme(
    legend.position = "none",
    strip.text      = element_blank(),
    axis.text.x     = element_blank(),
    axis.ticks.x    = element_blank(),
    panel.spacing   = unit(2.5, "mm")
  )


# 6. Key -----------------------------------------------------------------------

# The key is one row below the panels: a box for each estimator, the circle
# of the mean and the dashed line of the theoretical bias. Horizontal
# positions are in mm across the width of the figure.
key_labels <- c(
  sub(",''[phantom(1)]", "", unname(estimator_labels[panel_estimators]),
      fixed = TRUE),
  "'Mean'",
  "'Theoretical bias (Eq. 18)'"
)

# Width of each label (mm) at 7 points, measured on the device that writes
# the PNG file.
ragg::agg_png(tempfile(fileext = ".png"), width = figure_width_double,
              height = 20, units = "mm", res = figure_dpi)
key_label_width <- numeric(length(key_labels))
for (i in seq_along(key_labels)) {
  key_label_width[i] <- grid::convertWidth(
    grid::grobWidth(grid::textGrob(
      parse(text = key_labels[i]),
      gp = grid::gpar(fontsize = 7, fontfamily = font_family)
    )),
    "mm", valueOnly = TRUE
  )
}
invisible(dev.off())

# Each entry is a symbol 4.5 mm wide, a space of 1.2 mm and the label; the
# entries are spread evenly over 165 mm.
key_symbol_width <- 4.5
key_entry_width  <- key_label_width + key_symbol_width + 1.2
key_gap          <- (165 - sum(key_entry_width)) / 6
key_x            <- 2.5 + cumsum(c(0, head(key_entry_width + key_gap, -1)))

# Left edge of the box of each estimator, and the ends and colour of its
# median line (white and shorter for ST (X1, X2), as in the panels).
key_box_x        <- key_x[1:5]
key_median_start <- ifelse(panel_estimators == "M4",
                           key_box_x + key_symbol_width * 0.225, key_box_x)
key_median_end   <- ifelse(panel_estimators == "M4",
                           key_box_x + key_symbol_width * 0.775,
                           key_box_x + key_symbol_width)
key_median_colour <- ifelse(panel_estimators == "M4", "white",
                            colour_box_outline)

key <- ggplot() +
  annotate("rect", xmin = key_box_x, xmax = key_box_x + key_symbol_width,
           ymin = 0.35, ymax = 0.65,
           fill = unname(estimator_colours[panel_estimators]),
           colour = colour_box_outline, linewidth = 0.3) +
  annotate("segment", x = key_median_start, xend = key_median_end,
           y = 0.5, yend = 0.5, colour = key_median_colour,
           linewidth = 0.4) +
  annotate("point", x = key_x[6] + key_symbol_width / 2, y = 0.5,
           shape = 21, fill = "white", colour = "black", size = 1.15,
           stroke = 0.35) +
  annotate("segment", x = key_x[7] - 0.5,
           xend = key_x[7] + key_symbol_width + 0.5, y = 0.5, yend = 0.5,
           linewidth = 0.35, colour = "black", linetype = linetype_dashed) +
  annotate("text", x = key_x + key_symbol_width + 1.2, y = 0.5, hjust = 0,
           label = key_labels, parse = TRUE, size = 7 / points_per_mm,
           family = font_family) +
  coord_cartesian(xlim = c(0, figure_width_double), ylim = c(0, 1),
                  expand = FALSE, clip = "off") +
  theme_void() +
  theme(plot.margin = margin(0, 0, 0, 0))


# 7. Combine and save ----------------------------------------------------------

figure <- panels / key +
  plot_layout(heights = c(figure_height - 7, 5))

ggsave(
  file.path(dir_figures, "fig03_coefficient_error_groups_AB.png"),
  plot   = figure,
  width  = figure_width_double,
  height = figure_height,
  units  = "mm",
  dpi    = figure_dpi,
  bg     = "white",
  device = ragg::agg_png
)
ggsave(
  file.path(dir_figures, "fig03_coefficient_error_groups_AB.pdf"),
  plot   = figure,
  width  = figure_width_double,
  height = figure_height,
  units  = "mm",
  dpi    = figure_dpi,
  bg     = "white",
  device = cairo_pdf
)
