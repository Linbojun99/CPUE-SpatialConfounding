# ==============================================================================
# figures/figS09_convergence.R
#
# Figure S9: numerical diagnostics of the six fitted models in the 17
# scenarios, with 50 replicates per scenario and model.
#   (a)  largest absolute gradient over the replicates on a log10 scale;
#        outlined entries are at or above the limit of a converged fit
#        (0.001)
#   (b)  number of fits in which an estimated field SD was at the boundary
#        (below 0.01), with the number of refits without that field in
#        parentheses (Section 2.4); empty entries are zero
# RSR has no optimization of its own and is not shown.
#
# Input:
#   output/performance/performance_fits.csv   06_performance.R
# Output:
#   output/figures/figS09_convergence.png
#   output/figures/figS09_convergence.pdf
# ==============================================================================

source("00_settings.R")
source("figures/00_figure_style.R")


# 1. Input ---------------------------------------------------------------------

fits <- fread(file.path(dir_performance, "performance_fits.csv"))
fits <- fits[scenario_id %in% scenario_order]


# 2. Diagnostics by scenario and model -----------------------------------------

# Fitted models in the order of the columns (RSR, M7, is left out).
column_estimators <- estimator_order[estimator_order != "M7"]

diagnostics <- fits[, .(
  boundary     = sum(boundary_sigma),
  refitted     = sum(boundary_refit),
  max_gradient = max(max_abs_gradient, na.rm = TRUE)
), by = .(scenario_id, model_id)]

# Scenarios from top to bottom, models from left to right.
diagnostics[, scenario := factor(scenario_id, levels = rev(scenario_order))]
diagnostics[, model := factor(model_id, levels = column_estimators)]

# Entry of panel (b): the number of flagged fits, followed by the number of
# refits in parentheses when there are any.
diagnostics[, label := as.character(boundary)]
diagnostics[refitted > 0, label := sprintf("%d (%d)", boundary, refitted)]
diagnostics[boundary == 0, label := ""]


# 3. Scales and theme shared by the panels -------------------------------------

# Estimator labels on two lines above the columns.
scale_columns <- scale_x_discrete(
  labels   = parse(text = unname(estimator_labels_two_lines[
    column_estimators
  ])),
  position = "top",
  expand   = c(0, 0)
)

theme_table <- theme_figure +
  theme(
    panel.border = element_blank(),
    axis.ticks   = element_blank(),
    axis.text.x  = element_text(size = 7),
    legend.title = element_text(size = 7.5),
    legend.text  = element_text(size = 7),
    plot.title   = element_text(margin = margin(b = 2, unit = "mm"))
  )


# 4. Panel (a): largest absolute gradient --------------------------------------

panel_gradient <- ggplot(
  diagnostics,
  aes(model, scenario, fill = log10(max_gradient))
) +
  geom_tile(colour = "white", linewidth = 0.3) +
  geom_tile(data = diagnostics[max_gradient >= gradient_limit], fill = NA,
            colour = "black", linewidth = 0.5, width = 0.94, height = 0.9) +
  scale_columns +
  scale_y_discrete(expand = c(0, 0)) +
  scale_fill_gradient(
    low    = "#EEF3F8",
    high   = "#08306B",
    name   = expression(paste("log", ""[10],
                              " of the maximum absolute gradient")),
    labels = scales::label_number(style_negative = "minus"),
    guide  = guide_colourbar(
      title.position  = "top",
      title.hjust     = 0.5,
      barwidth        = unit(42, "mm"),
      barheight       = unit(2, "mm"),
      ticks.colour    = "#333333",
      frame.colour    = "#777777",
      frame.linewidth = 0.2
    )
  ) +
  labs(x = NULL, y = NULL, title = "(a) Maximum absolute gradient") +
  theme_table


# 5. Panel (b): field SD at the boundary ---------------------------------------

panel_boundary <- ggplot(diagnostics, aes(model, scenario, fill = boundary)) +
  geom_tile(colour = "white", linewidth = 0.3) +
  geom_text(aes(label = label), size = 7.5 / points_per_mm,
            family = font_family, colour = "#111111") +
  scale_columns +
  scale_y_discrete(expand = c(0, 0)) +
  scale_fill_gradient(
    low    = "#F7F9FB",
    high   = "#9ECAE1",
    limits = c(0, max(diagnostics$boundary)),
    breaks = c(0, 5, 10),
    guide  = "none"
  ) +
  labs(x = NULL, y = NULL, title = "(b) Field SD boundary flags (refits)") +
  theme_table


# 6. Combine and save ----------------------------------------------------------

figure <- panel_gradient + panel_boundary +
  plot_layout(ncol = 2)

ggsave(
  file.path(dir_figures, "figS09_convergence.png"),
  plot   = figure,
  width  = figure_width_double,
  height = 112,
  units  = "mm",
  dpi    = figure_dpi,
  bg     = "white",
  device = ragg::agg_png
)
ggsave(
  file.path(dir_figures, "figS09_convergence.pdf"),
  plot   = figure,
  width  = figure_width_double,
  height = 112,
  units  = "mm",
  dpi    = figure_dpi,
  bg     = "white",
  device = cairo_pdf
)
