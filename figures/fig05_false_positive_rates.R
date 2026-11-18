# ==============================================================================
# figures/fig05_false_positive_rates.R
#
# Figure 5: false positive rates (%) for coefficients with a true value of
# zero, as a shaded table. Rows are estimators and columns are scenarios.
#   (a)  the X2 coefficient in the 15 scenarios in which X2 has no effect
#        (groups A, C, D and E)
#   (b)  the X1 coefficient in C1 and C2, where neither covariate has an
#        effect
# Above each column are the group name, the scenario label and the values of
# rho_omega, rho_12 and c_r (Table 2). Each entry is the percentage of fits
# in which a two-sided Wald test rejected a zero coefficient at the 5% level;
# an entry is empty where the coefficient is not in the model. The same
# values are listed in Table S4.
#
# Input:
#   output/performance/performance_coefficient_summary.csv   06_performance.R
#   output/performance/rsr_coefficient_summary.csv           05_rsr.R
# Output:
#   output/figures/fig05_false_positive_rates.png
#   output/figures/fig05_false_positive_rates.pdf
# ==============================================================================

source("00_settings.R")
source("figures/00_figure_style.R")


# 1. Settings of the figure ----------------------------------------------------

# Estimators from the top row to the bottom row. ST (none) contains neither
# coefficient.
row_estimators <- c("M1", "M3", "M4", "M6", "M5", "M7")

# Scenarios of block (a) and of block (b).
scenarios_x2 <- scenario_order[substr(scenario_order, 1, 1) %in%
                                 c("A", "C", "D", "E")]
scenarios_x1 <- c("C1", "C2")

# Names of the scenario groups (Section 2.3).
group_names <- c(
  A = "Collinearity and confounding",
  C = "No effects",
  D = "Spatial scale",
  E = "Trend and\ncatchability"
)

# Space between scenario groups, in column widths.
group_gap <- 0.35

# Vertical positions (rows are at 6 to 1): scenario label, the three design
# values and the rule below the group name.
y_scenario  <- 8.92
y_rho_omega <- 8.24
y_rho_12    <- 7.62
y_c_r       <- 7.0
y_rule      <- 9.36
y_title     <- 11.3

# Horizontal position of the right end of the row labels.
x_row_label <- 0.36


# 2. Input ---------------------------------------------------------------------

coefficient_summary <- fread(
  file.path(dir_performance, "performance_coefficient_summary.csv")
)
rsr_coefficient_summary <- fread(
  file.path(dir_performance, "rsr_coefficient_summary.csv")
)


# 3. False positive rates ------------------------------------------------------

# For a coefficient with a true value of zero, the share of fits in which a
# zero coefficient was rejected is the false positive rate.
rates <- rbind(
  coefficient_summary[truth == 0,
                      .(scenario_id, model_id, term, rate = reject_zero)],
  rsr_coefficient_summary[truth == 0,
                          .(scenario_id, model_id = "M7", term,
                            rate = false_positive)]
)
rates <- rates[scenario_id %in% scenario_order & !is.na(rate)]
rates[, percent := 100 * rate]
rates[, label := sprintf("%.0f", percent)]

rates_x2 <- rates[term == "X2" & scenario_id %in% scenarios_x2]
rates_x1 <- rates[term == "X1" & scenario_id %in% scenarios_x1]


# 4. Positions of the cells ----------------------------------------------------

# Columns of block (a): one per scenario, with a gap between groups.
# Columns of block (b): to the right of block (a).
group_x2  <- substr(scenarios_x2, 1, 1)
column_x2 <- seq_along(scenarios_x2) +
  group_gap * (match(group_x2, c("A", "C", "D", "E")) - 1)
column_x1 <- max(column_x2) + c(2.3, 3.3)

column_scenarios <- c(scenarios_x2, scenarios_x1)
column_x         <- c(column_x2, column_x1)

rates_x2[, x := column_x2[match(scenario_id, scenarios_x2)]]
rates_x1[, x := column_x1[match(scenario_id, scenarios_x1)]]
cells <- rbind(rates_x2, rates_x1)
cells[, y := length(row_estimators) + 1 - match(model_id, row_estimators)]

# Outline of every cell, including the empty ones.
cell_outlines <- CJ(x = column_x, y = seq_along(row_estimators))

# Design values above each column.
column_design <- scenarios[match(column_scenarios, scenario)]
design_values <- rbind(
  data.table(x = column_x, y = y_rho_omega,
             label = as.character(column_design$rho_omega)),
  data.table(x = column_x, y = y_rho_12,
             label = as.character(column_design$rho_12)),
  data.table(x = column_x, y = y_c_r,
             label = as.character(column_design$range_ratio))
)

# Rule and name of each scenario group above its columns.
group_spans <- rbind(
  data.table(x = column_x2, group = group_x2)[
    , .(x_start = min(x) - 0.47, x_end = max(x) + 0.47), by = group
  ],
  data.table(group = "C", x_start = column_x1[1] - 0.47,
             x_end = column_x1[2] + 0.47)
)
group_spans[, x_middle := (x_start + x_end) / 2]
group_spans[, name := group_names[group]]

# Row labels: the estimator labels without the empty subscript.
row_labels <- sub(",''[phantom(1)]", "",
                  unname(estimator_labels[row_estimators]), fixed = TRUE)


# 5. Figure --------------------------------------------------------------------

figure <- ggplot(cells, aes(x, y)) +
  geom_tile(data = cell_outlines, fill = NA, colour = "#DDDDDD",
            linewidth = 0.15, width = 0.99, height = 0.96) +
  geom_tile(aes(fill = percent), colour = "#DDDDDD", linewidth = 0.15,
            width = 0.99, height = 0.96) +
  geom_text(aes(label = label), size = 7.5 / points_per_mm,
            colour = "#111111", family = font_family) +
  # Row labels and the symbols of the three design values.
  annotate("text", x = x_row_label, y = rev(seq_along(row_estimators)),
           label = row_labels, parse = TRUE, hjust = 1,
           size = 7.5 / points_per_mm, family = font_family) +
  annotate("text", x = x_row_label, y = c(y_rho_omega, y_rho_12, y_c_r),
           label = c(symbol_rho_omega, symbol_rho_12, "c[r]"), parse = TRUE,
           hjust = 1, size = 7.5 / points_per_mm, family = font_family) +
  # Column headers.
  annotate("text", x = column_x, y = y_scenario, label = column_scenarios,
           fontface = "bold", size = 7.5 / points_per_mm,
           family = font_family) +
  geom_text(data = design_values, aes(label = label),
            size = 7 / points_per_mm, colour = "#111111",
            family = font_family) +
  geom_segment(data = group_spans,
               aes(x = x_start, xend = x_end, y = y_rule, yend = y_rule),
               inherit.aes = FALSE, linewidth = 0.3, colour = "black") +
  geom_text(data = group_spans,
            aes(x = x_middle, y = y_rule + 0.14, label = name),
            inherit.aes = FALSE, vjust = 0, lineheight = 0.9,
            size = 7 / points_per_mm, family = font_family) +
  # Titles of the two blocks.
  annotate("text", x = min(column_x2) - 0.5, y = y_title,
           label = "paste('(a) ',X[2],' coefficient (true value zero)')",
           parse = TRUE, hjust = 0, size = 8.5 / points_per_mm,
           family = font_family) +
  annotate("text", x = column_x1[1] - 0.5, y = y_title,
           label = "paste('(b) ',X[1],' coefficient')",
           parse = TRUE, hjust = 0, size = 8.5 / points_per_mm,
           family = font_family) +
  scale_fill_gradientn(
    colours = c("#FFFFFF", "#FDDCCB", "#F4A582", "#E5846B"),
    values  = c(0, 0.25, 0.65, 1),
    limits  = c(0, 100),
    breaks  = seq(0, 100, 25),
    name    = "False positive rate (%), nominal level 5%"
  ) +
  coord_cartesian(xlim = c(0.5 - 2.78, max(column_x1) + 0.5 + 1.55),
                  ylim = c(0.5, 11.75), expand = FALSE, clip = "off") +
  theme_figure +
  theme(
    panel.border       = element_blank(),
    axis.text          = element_blank(),
    axis.ticks         = element_blank(),
    axis.title         = element_blank(),
    legend.position    = "bottom",
    legend.title       = element_text(size = 7.5),
    legend.text        = element_text(size = 7),
    legend.margin      = margin(0.5, 0, 0, 0, "mm"),
    legend.box.spacing = unit(1.5, "mm"),
    plot.margin        = margin(1, 1, 1, 1, "mm")
  ) +
  guides(fill = guide_colourbar(
    title.position  = "top",
    title.hjust     = 0.5,
    barwidth        = unit(60, "mm"),
    barheight       = unit(2, "mm"),
    display         = "rectangles",
    ticks.colour    = "#333333",
    frame.colour    = "#777777",
    frame.linewidth = 0.2
  ))


# 6. Save ----------------------------------------------------------------------

ggsave(
  file.path(dir_figures, "fig05_false_positive_rates.png"),
  plot   = figure,
  width  = figure_width_double,
  height = 70,
  units  = "mm",
  dpi    = figure_dpi,
  bg     = "white",
  device = ragg::agg_png
)
ggsave(
  file.path(dir_figures, "fig05_false_positive_rates.pdf"),
  plot   = figure,
  width  = figure_width_double,
  height = 70,
  units  = "mm",
  dpi    = figure_dpi,
  bg     = "white",
  device = cairo_pdf
)
