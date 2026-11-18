# ==============================================================================
# figures/fig06_index_rmse_differences.R
#
# Figure 6: mean within-replicate differences in log-index RMSE (Eq. 20)
# from ST (none) in the 17 scenarios, as a shaded table. Rows are estimators
# and columns are scenarios. Above each column are the group name, the
# scenario label and the values of rho_omega, rho_12 and c_r (Table 2).
# Printed values are the mean differences multiplied by 1,000 and rounded to
# the nearest integer; negative values indicate a smaller error than
# ST (none). The colour scale covers the four other ST models and RSR, from
# -25 to +25; the entries of the GLMM lie above this range and are grey.
#
# Input:
#   output/performance/performance_paired.csv   06_performance.R
# Output:
#   output/figures/fig06_index_rmse_differences.png
#   output/figures/fig06_index_rmse_differences.pdf
# ==============================================================================

source("00_settings.R")
source("figures/00_figure_style.R")


# 1. Settings of the figure ----------------------------------------------------

# Estimators from the top row to the bottom row. ST (none), M2, is the
# reference and has no row.
row_estimators <- c("M1", "M3", "M4", "M6", "M5", "M7")

# Names of the scenario groups (Section 2.3).
group_names <- c(
  A = "Collinearity and confounding",
  B = "Two effects",
  C = "No effects",
  D = "Spatial scale",
  E = "Trend and\ncatchability"
)

# Space between scenario groups, in column widths.
group_gap <- 0.35

# Limit of the colour scale, in units of 1000 x the RMSE difference.
colour_limit <- 25

# Vertical positions (rows are at 6 to 1): scenario label, the three design
# values and the rule below the group name.
y_scenario  <- 8.92
y_rho_omega <- 8.24
y_rho_12    <- 7.62
y_c_r       <- 7.0
y_rule      <- 9.36

# Horizontal position of the right end of the row labels.
x_row_label <- 0.36


# 2. Input ---------------------------------------------------------------------

# Mean differences in log-index RMSE of model a from model b, ST (none),
# within replicates.
paired <- fread(file.path(dir_performance, "performance_paired.csv"))
paired <- paired[scenario_id %in% scenario_order & b == "M2" &
                   a %in% row_estimators]


# 3. Positions and labels of the cells -----------------------------------------

# One column per scenario, with a gap between groups.
scenario_group <- substr(scenario_order, 1, 1)
column_x       <- seq_along(scenario_order) +
  group_gap * (match(scenario_group, names(group_names)) - 1)

cells <- paired[, .(scenario_id, model_id = a,
                    difference_1000 = 1000 * diff_rmse_cen)]
cells[, x := column_x[match(scenario_id, scenario_order)]]
cells[, y := length(row_estimators) + 1 - match(model_id, row_estimators)]

# Printed value: the rounded integer with its sign.
cells[, rounded := as.integer(round(difference_1000))]
cells[, label := ifelse(
  rounded == 0, "0",
  ifelse(rounded > 0, paste0("+", rounded),
         paste0(minus_sign, abs(rounded)))
)]

# Design values above each column.
column_design <- scenarios[match(scenario_order, scenario)]
design_values <- rbind(
  data.table(x = column_x, y = y_rho_omega,
             label = as.character(column_design$rho_omega)),
  data.table(x = column_x, y = y_rho_12,
             label = as.character(column_design$rho_12)),
  data.table(x = column_x, y = y_c_r,
             label = as.character(column_design$range_ratio))
)

# Rule and name of each scenario group above its columns.
group_spans <- data.table(x = column_x, group = scenario_group)[
  , .(x_start = min(x) - 0.47, x_end = max(x) + 0.47), by = group
]
group_spans[, x_middle := (x_start + x_end) / 2]
group_spans[, name := group_names[group]]

# Row labels: the estimator labels without the empty subscript.
row_labels <- sub(",''[phantom(1)]", "",
                  unname(estimator_labels[row_estimators]), fixed = TRUE)


# 4. Figure --------------------------------------------------------------------

figure <- ggplot(cells[model_id != "M1"], aes(x, y)) +
  # Entries of the GLMM in grey, the others shaded by their value.
  geom_tile(data = cells[model_id == "M1"], fill = "#EBEBEB",
            colour = "#DDDDDD", linewidth = 0.15, width = 0.99,
            height = 0.96) +
  geom_tile(aes(fill = difference_1000), colour = "#DDDDDD",
            linewidth = 0.15, width = 0.99, height = 0.96) +
  geom_text(data = cells, aes(label = label), size = 7.5 / points_per_mm,
            colour = "#111111", family = font_family) +
  # Row labels and the symbols of the three design values.
  annotate("text", x = x_row_label, y = rev(seq_along(row_estimators)),
           label = row_labels, parse = TRUE, hjust = 1,
           size = 7.5 / points_per_mm, family = font_family) +
  annotate("text", x = x_row_label, y = c(y_rho_omega, y_rho_12, y_c_r),
           label = c(symbol_rho_omega, symbol_rho_12, "c[r]"), parse = TRUE,
           hjust = 1, size = 7.5 / points_per_mm, family = font_family) +
  # Column headers.
  annotate("text", x = column_x, y = y_scenario, label = scenario_order,
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
  scale_fill_gradient2(
    low      = "#4292C6",
    mid      = "white",
    high     = "#F4A582",
    midpoint = 0,
    limits   = c(-colour_limit, colour_limit),
    breaks   = c(-20, -10, 0, 10, 20),
    labels   = c(paste0(minus_sign, "20"), paste0(minus_sign, "10"), "0",
                 "+10", "+20"),
    name     = paste("1000 \u00d7 paired RMSE difference of the spatial",
                     "estimators (negative = smaller error)")
  ) +
  coord_cartesian(xlim = c(0.5 - 2.78, max(column_x) + 0.5 + 0.15),
                  ylim = c(0.5, 10.75), expand = FALSE, clip = "off") +
  labs(
    title = "1000 \u00d7 mean paired RMSE difference from ST (none)",
    x     = NULL,
    y     = NULL
  ) +
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
    plot.margin        = margin(1, 1, 2, 1, "mm")
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


# 5. Save ----------------------------------------------------------------------

ggsave(
  file.path(dir_figures, "fig06_index_rmse_differences.png"),
  plot   = figure,
  width  = figure_width_double,
  height = 70,
  units  = "mm",
  dpi    = figure_dpi,
  bg     = "white",
  device = ragg::agg_png
)
ggsave(
  file.path(dir_figures, "fig06_index_rmse_differences.pdf"),
  plot   = figure,
  width  = figure_width_double,
  height = 70,
  units  = "mm",
  dpi    = figure_dpi,
  bg     = "white",
  device = cairo_pdf
)
