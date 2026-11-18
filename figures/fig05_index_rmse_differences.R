# ==============================================================================
# figures/fig05_index_rmse_differences.R
#
# Figure 5: mean within-replicate differences in log-index RMSE (Eq. 20)
# from ST (none) in the 17 scenarios (Table 2), as a shaded table. Rows are
# estimators and columns are scenarios, with the letter of the scenario
# group above each block of columns. Printed values are the mean differences
# multiplied by 1,000 and rounded to the nearest integer; negative values
# indicate a smaller RMSE than ST (none). The colour scale covers the four
# other ST models and RSR, from -25 to +25; the entries of the GLMM lie
# above this range and are grey.
#
# Input:
#   output/performance/performance_paired.csv   06_performance.R
# Output:
#   output/figures/fig05_index_rmse_differences.png
#   output/figures/fig05_index_rmse_differences.pdf
# ==============================================================================

source("00_settings.R")
source("figures/00_figure_style.R")


# 1. Settings of the figure ----------------------------------------------------

# Estimators from the top row to the bottom row. ST (none), M2, is the
# reference and has no row.
row_estimators <- c("M1", "M3", "M4", "M6", "M5", "M7")

# Space between scenario groups, in column widths.
group_gap <- 0.42

# Limit of the colour scale, in units of 1,000 x the RMSE difference.
colour_limit <- 25

# Vertical positions (rows are at 6 to 1): scenario labels, the rule below
# the group letters and the group letters.
y_scenario <- 6.8
y_rule     <- 7.14
y_group    <- 7.6

figure_height <- 60


# 2. Input ---------------------------------------------------------------------

# Mean differences in log-index RMSE of model a from model b, ST (none),
# within replicates.
paired <- fread(file.path(dir_performance, "performance_paired.csv"))
paired <- paired[scenario_id %in% scenario_order & b == "M2" &
                   a %in% row_estimators]


# 3. Positions and labels of the cells -----------------------------------------

# One column per scenario, with a gap between groups.
scenario_group <- substr(scenario_order, 1, 1)
group_letters  <- unique(scenario_group)
column_x       <- seq_along(scenario_order) +
  group_gap * (match(scenario_group, group_letters) - 1)

cells <- paired[, .(scenario_id, model_id = a,
                    difference_1000 = 1000 * diff_rmse_cen)]
cells[, x := column_x[match(scenario_id, scenario_order)]]
cells[, y := length(row_estimators) + 1 - match(model_id, row_estimators)]

# Cells are drawn row by row from the top, and from left to right in a row.
setorder(cells, -y, x)

# Printed value: the rounded integer with its sign.
cells[, rounded := as.integer(round(difference_1000))]
cells[, label := ifelse(
  rounded == 0, "0",
  ifelse(rounded > 0, paste0("+", rounded),
         paste0(minus_sign, abs(rounded)))
)]

# Rule and letter of each scenario group above its columns.
group_spans <- data.table(x = column_x, group = scenario_group)[
  , .(x_start = min(x) - 0.46, x_end = max(x) + 0.46), by = group
]


# 4. Figure --------------------------------------------------------------------

figure <- ggplot(cells, aes(x, y)) +
  # Entries of the GLMM in grey, the others shaded by their value.
  geom_tile(data = cells[model_id == "M1"], width = 1, height = 0.98,
            fill = "#EAEAEA", colour = "white", linewidth = 0.25) +
  geom_tile(data = cells[model_id != "M1"], aes(fill = difference_1000),
            width = 1, height = 0.98, colour = "#DDDDDD",
            linewidth = 0.16) +
  geom_text(aes(label = label), size = 7.7 / points_per_mm,
            family = font_family) +
  # Column headers.
  annotate("text", x = column_x, y = y_scenario, label = scenario_order,
           size = 7.6 / points_per_mm, family = font_family) +
  geom_segment(data = group_spans,
               aes(x = x_start, xend = x_end, y = y_rule, yend = y_rule),
               inherit.aes = FALSE, linewidth = 0.28) +
  geom_text(data = group_spans,
            aes(x = (x_start + x_end) / 2, y = y_group, label = group),
            inherit.aes = FALSE, family = font_family, fontface = "bold",
            size = 8 / points_per_mm) +
  # Row labels: the estimators on the vertical axis.
  scale_y_continuous(
    breaks = seq_along(row_estimators),
    labels = parse(text = unname(estimator_labels[rev(row_estimators)]))
  ) +
  scale_fill_gradient2(
    low      = "#4292C6",
    mid      = "white",
    high     = "#F4A582",
    midpoint = 0,
    limits   = c(-colour_limit, colour_limit),
    breaks   = c(-20, -10, 0, 10, 20),
    labels   = c(paste0(minus_sign, "20"), paste0(minus_sign, "10"), "0",
                 "10", "20"),
    name     = "Paired RMSE difference \u00d7 1,000"
  ) +
  guides(fill = guide_colourbar(
    title.position  = "top",
    title.hjust     = 0.5,
    barwidth        = unit(58, "mm"),
    barheight       = unit(2, "mm"),
    frame.colour    = "#888888",
    frame.linewidth = 0.2
  )) +
  coord_cartesian(xlim = range(column_x) + c(-0.5, 0.5),
                  ylim = c(0.48, 8.05), expand = FALSE) +
  labs(x = NULL, y = NULL) +
  theme_figure +
  theme(
    panel.border       = element_blank(),
    axis.text.x        = element_blank(),
    axis.ticks         = element_blank(),
    axis.text.y        = element_text(size = 7.7),
    legend.box.spacing = unit(1, "mm")
  )


# 5. Save ----------------------------------------------------------------------

ggsave(
  file.path(dir_figures, "fig05_index_rmse_differences.png"),
  plot   = figure,
  width  = figure_width_double,
  height = figure_height,
  units  = "mm",
  dpi    = figure_dpi,
  bg     = "white",
  device = ragg::agg_png
)
ggsave(
  file.path(dir_figures, "fig05_index_rmse_differences.pdf"),
  plot   = figure,
  width  = figure_width_double,
  height = figure_height,
  units  = "mm",
  dpi    = figure_dpi,
  bg     = "white",
  device = cairo_pdf
)
