# ==============================================================================
# figures/fig09_gaussian_simulation.R
#
# Figure 9: results of the Gaussian simulation with basis dimension k = 50.
# Points show means and horizontal segments the 10% and 90% quantiles over
# the 50 replicates, for the four estimators.
#   (a), (b)  relative error of the X1 coefficient in G0 and G1, with the
#             theoretical bias (Eq. 18) as a dashed line in (b)
#   (c)       change error (Eq. 22) in G1t and G2, which share data and fits
#   (d)       G2-minus-G1t difference in change error minus its first-order
#             approximation, calculated within each replicate
#
# Input:
#   output/gaussian/reporting/k50/coefficients.csv   08a_gaussian_simulation.R
#   output/gaussian/reporting/k50/performance.csv    08a_gaussian_simulation.R
#   output/gaussian/paired_replicates.csv            08b_gaussian_summaries.R
#   output/gaussian/diagnostic_summary.csv           08b_gaussian_summaries.R
# Output:
#   output/figures/fig09_gaussian_simulation.png
#   output/figures/fig09_gaussian_simulation.pdf
# ==============================================================================

source("00_settings.R")
source("figures/00_figure_style.R")


# 1. Input ---------------------------------------------------------------------

run_folder <- file.path(dir_gaussian, "reporting", "k50")

coefficient_estimates <- fread(file.path(run_folder, "coefficients.csv"))
performance_measures  <- fread(file.path(run_folder, "performance.csv"))

paired_replicates  <- fread(file.path(dir_gaussian, "paired_replicates.csv"))
diagnostic_summary <- fread(file.path(dir_gaussian, "diagnostic_summary.csv"))


# 2. Means and quantiles over replicates ---------------------------------------

# (a), (b) Relative error of the X1 coefficient (%).
coefficient_summary <- coefficient_estimates[
  term == "X1" & scenario %in% c("G0", "G1"),
  .(mean = 100 * mean(relative_error),
    lo   = 100 * quantile(relative_error, 0.1),
    hi   = 100 * quantile(relative_error, 0.9)),
  by = .(scenario, model)
]

# (c) Change error in G1t (X1 retained in prediction) and G2 (X1 fixed at
# zero in prediction).
change_summary <- performance_measures[
  scenario %in% c("G1t", "G2"),
  .(mean = mean(first_last_error),
    lo   = quantile(first_last_error, 0.1),
    hi   = quantile(first_last_error, 0.9)),
  by = .(scenario, model)
]

# (d) Within each replicate, the change error in G2 minus that in G1t, minus
# the first-order approximation (the column "remainder").
remainder_summary <- paired_replicates[
  k == 50,
  .(mean = mean(remainder),
    lo   = quantile(remainder, 0.1),
    hi   = quantile(remainder, 0.9)),
  by = model
]

# Theoretical bias of the X1 coefficient in G1 (%).
theoretical_bias <- diagnostic_summary[
  setting == "reporting" & k == 50 & scenario == "G1", projection
]

# Vertical position of each estimator: the first estimator is at the top.
# In (c), G1t is drawn above G2 within the row of the estimator.
n_estimators <- length(gaussian_estimator_order)

coefficient_summary[, y := n_estimators + 1 -
                      match(model, gaussian_estimator_order)]
change_summary[, y := n_estimators + 1 -
                 match(model, gaussian_estimator_order)]
remainder_summary[, y := n_estimators + 1 -
                    match(model, gaussian_estimator_order)]
change_summary[, y_scenario := y + ifelse(scenario == "G1t", 0.16, -0.16)]


# 3. Scales and theme shared by the panels -------------------------------------

estimator_scales <- list(
  scale_colour_manual(values = gaussian_estimator_colours),
  scale_fill_manual(values = gaussian_estimator_colours),
  scale_shape_manual(values = gaussian_estimator_shapes),
  guides(colour = "none", fill = "none", shape = "none"),
  scale_y_continuous(
    breaks = n_estimators:1,
    labels = unname(gaussian_estimator_labels[gaussian_estimator_order])
  )
)

theme_panel <- theme_figure +
  theme(
    legend.position = "none",
    plot.title      = element_text(size = 9, margin = margin(b = 3)),
    axis.title.x    = element_text(size = 8, margin = margin(t = 4)),
    plot.margin     = margin(2, 3, 2, 2, "mm")
  )


# 4. Panels (a) and (b): relative error of the X1 coefficient ------------------

coefficient_panels <- list()

for (panel_scenario in c("G0", "G1")) {

  panel <- ggplot(
    coefficient_summary[scenario == panel_scenario],
    aes(y = y, colour = model, fill = model, shape = model)
  ) +
    geom_vline(xintercept = 0, colour = colour_reference,
               linewidth = line_width_reference)

  # Theoretical bias in G1, with its rounded value above the line.
  if (panel_scenario == "G1") {
    panel <- panel +
      annotate("segment", x = theoretical_bias, xend = theoretical_bias,
               y = 0.5, yend = 4.28, linetype = linetype_dashed,
               linewidth = 0.35, colour = "black") +
      annotate("text", x = theoretical_bias, y = 4.55,
               label = paste0(round(theoretical_bias), "%"), hjust = 0.5,
               size = 7.5 / points_per_mm, family = font_family)
  }

  panel <- panel +
    geom_segment(aes(x = lo, xend = hi, yend = y), linewidth = 0.45) +
    geom_point(aes(x = mean), size = 2.25, stroke = 0.55) +
    estimator_scales +
    scale_x_continuous(
      breaks = c(-25, 0, 25, 50, 75, 100),
      labels = scales::label_number(style_negative = "minus")
    ) +
    coord_cartesian(xlim = c(-30, 105), ylim = c(0.5, 4.8), expand = FALSE) +
    labs(
      x     = expression(paste("Relative error in the ", X[1],
                               " coefficient (%)")),
      y     = NULL,
      title = ifelse(panel_scenario == "G0", "(a) G0", "(b) G1")
    ) +
    theme_panel

  coefficient_panels[[panel_scenario]] <- panel
}


# 5. Panel (c): change error in G1t and G2 -------------------------------------

panel_change <- ggplot(
  change_summary,
  aes(colour = model, fill = model, shape = model)
) +
  geom_vline(xintercept = 0, colour = colour_reference,
             linewidth = line_width_reference) +
  geom_segment(aes(x = lo, xend = hi, y = y_scenario, yend = y_scenario,
                   linetype = scenario), linewidth = 0.45) +
  geom_point(aes(x = mean, y = y_scenario), size = 2.2, stroke = 0.55) +
  estimator_scales +
  scale_linetype_manual(values = c(G1t = "solid", G2 = linetype_dashed),
                        breaks = c("G1t", "G2")) +
  scale_x_continuous(
    breaks = c(-0.6, -0.4, -0.2, 0, 0.2),
    labels = scales::label_number(style_negative = "minus")
  ) +
  coord_cartesian(xlim = c(-0.6, 0.3), ylim = c(0.5, 4.6), expand = FALSE) +
  labs(x = "Change error", y = NULL, title = "(c) G1t and G2") +
  theme_panel +
  theme(
    legend.position  = "bottom",
    legend.text      = element_text(size = font_size_legend),
    legend.key.width = unit(8, "mm")
  ) +
  guides(
    colour   = "none",
    fill     = "none",
    shape    = "none",
    linetype = guide_legend(
      nrow         = 1,
      override.aes = list(colour = "black", linewidth = 0.45)
    )
  )


# 6. Panel (d): deviation from the first-order approximation -------------------

panel_remainder <- ggplot(
  remainder_summary,
  aes(colour = model, fill = model, shape = model)
) +
  geom_vline(xintercept = 0, colour = colour_reference,
             linewidth = line_width_reference) +
  geom_segment(aes(x = lo, xend = hi, y = y, yend = y), linewidth = 0.45) +
  geom_point(aes(x = mean, y = y), size = 2.25, stroke = 0.55) +
  estimator_scales +
  scale_x_continuous(
    breaks = c(-0.10, -0.05, 0, 0.05, 0.10),
    labels = scales::label_number(style_negative = "minus")
  ) +
  coord_cartesian(xlim = c(-0.105, 0.105), ylim = c(0.5, 4.6),
                  expand = FALSE) +
  labs(
    x     = "Deviation from first-order approximation",
    y     = NULL,
    title = paste("(d) G2", minus_sign, "G1t")
  ) +
  theme_panel


# 7. Combine and save ----------------------------------------------------------

figure <- (coefficient_panels[["G0"]] | coefficient_panels[["G1"]]) /
  (panel_change | panel_remainder) +
  plot_layout(heights = c(1, 1.05), guides = "keep") &
  theme(legend.margin = margin(t = 0, b = 0),
        plot.margin = margin(2, 2, 2, 2, "mm"))

ggsave(
  file.path(dir_figures, "fig09_gaussian_simulation.png"),
  plot   = figure,
  width  = figure_width_double,
  height = 106,
  units  = "mm",
  dpi    = figure_dpi,
  bg     = "white",
  device = ragg::agg_png
)
ggsave(
  file.path(dir_figures, "fig09_gaussian_simulation.pdf"),
  plot   = figure,
  width  = figure_width_double,
  height = 106,
  units  = "mm",
  dpi    = figure_dpi,
  bg     = "white",
  device = cairo_pdf
)
