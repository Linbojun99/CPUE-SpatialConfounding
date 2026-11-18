# ==============================================================================
# figures/figS13_index_rmse.R
#
# Figure S13: log-index RMSE (Eq. 20) of the six fitted models and RSR in
# the 17 scenarios. Panels (a) to (e) show the scenario groups A to E in two
# rows. Boxes show the median and interquartile range over replicates and
# whiskers the most extreme values within 1.5 times the interquartile range;
# values beyond the whiskers are not drawn.
#
# Input:
#   output/performance/performance_fits.csv   06_performance.R
#   output/rsr/rsr_fit.csv                    05_rsr.R
# Output:
#   output/figures/figS13_index_rmse.png
#   output/figures/figS13_index_rmse.pdf
# ==============================================================================

source("00_settings.R")
source("figures/00_figure_style.R")


# 1. Input ---------------------------------------------------------------------

fits     <- fread(file.path(dir_performance, "performance_fits.csv"))
rsr_fits <- fread(file.path(dir_rsr, "rsr_fit.csv"))

# Fitted models: fits that converged. RSR: the index used in the paper, the
# uncorrected index times the ratio of the bias-corrected to the uncorrected
# index of ST (X1, X2).
fits     <- fits[scenario_id %in% scenario_order & pred_eval == TRUE]
rsr_fits <- rsr_fits[scenario_id %in% scenario_order &
                       index_method == "inherited_M4_bias_ratio"]


# 2. Log-index RMSE by fit -----------------------------------------------------

index_rmse <- rbind(
  fits[, .(scenario_id, model_id, rep, value = rmse_cen)],
  rsr_fits[, .(scenario_id, model_id, rep, value = rmse_cen)]
)
index_rmse[, model := factor(model_id, levels = estimator_order)]
index_rmse[, group := substr(scenario_id, 1, 1)]

# Panel titles: the names of the scenario groups (Section 2.3).
panel_titles <- c(
  A = "(a) Collinearity and confounding",
  B = "(b) Two effects",
  C = "(c) No effects",
  D = "(d) Spatial scale",
  E = "(e) Trend and catchability"
)

# The vertical axis starts at zero and ends 0.01 above the highest whisker.
whiskers <- index_rmse[, .(high = boxplot.stats(value)$stats[5]),
                       by = .(scenario_id, model_id)]
y_limits <- c(0, max(whiskers$high) + 0.01)


# 3. Two rows of panels --------------------------------------------------------

# Groups A and B in the first row, groups C, D and E in the second. Panel
# widths are proportional to the number of scenarios. The median of
# ST (X1, X2), the middle box of each scenario, is redrawn in white.
row_groups <- list(c("A", "B"), c("C", "D", "E"))
rows       <- list()

for (i in seq_along(row_groups)) {

  groups_in_row <- row_groups[[i]]
  row_data      <- index_rmse[group %in% groups_in_row]
  row_data[, scenario := factor(
    scenario_id,
    levels = scenario_order[substr(scenario_order, 1, 1) %in% groups_in_row]
  )]
  row_data[, panel := factor(panel_titles[group],
                             levels = panel_titles[groups_in_row])]

  median_st <- row_data[model_id == "M4", .(median = median(value)),
                        by = .(panel, scenario)]

  rows[[i]] <- ggplot(row_data,
                      aes(scenario, value, fill = model, colour = model)) +
    geom_boxplot(outlier.shape = NA, position = position_dodge(0.8),
                 width = 0.7, linewidth = 0.25, key_glyph = "polygon") +
    geom_errorbar(data = median_st,
                  aes(ymin = median, ymax = median, y = NULL, fill = NULL,
                      colour = NULL),
                  colour = "white", width = 0.1, linewidth = 0.45) +
    facet_grid(. ~ panel, scales = "free_x", space = "free_x") +
    scale_fill_manual(
      values = estimator_colours,
      breaks = estimator_order,
      labels = parse(text = unname(estimator_labels[estimator_order])),
      limits = estimator_order,
      drop   = FALSE
    ) +
    scale_colour_manual(
      values = rep(colour_box_outline, length(estimator_order)),
      breaks = estimator_order,
      labels = parse(text = unname(estimator_labels[estimator_order])),
      limits = estimator_order,
      drop   = FALSE
    ) +
    scale_x_discrete(expand = expansion(add = 0.5)) +
    scale_y_continuous(expand = c(0, 0)) +
    coord_cartesian(ylim = y_limits) +
    labs(x = "Scenario", y = "Log-index RMSE") +
    theme_figure +
    guides(fill = guide_legend(nrow = 1), colour = guide_legend(nrow = 1)) +
    theme(
      legend.text          = element_text(size = 7.5, hjust = 0),
      legend.key.spacing.x = unit(2.5, "mm")
    )
}


# 4. Combine and save ----------------------------------------------------------

# The first row has 8 scenarios and the second 9, so the first row takes 8
# of 9 columns of the layout.
figure <- rows[[1]] + rows[[2]] +
  plot_layout(
    design      = c(area(1, 1, 1, 8), area(2, 1, 2, 9)),
    guides      = "collect",
    axis_titles = "collect"
  ) &
  theme(legend.position = "bottom")

ggsave(
  file.path(dir_figures, "figS13_index_rmse.png"),
  plot   = figure,
  width  = figure_width_double,
  height = 118,
  units  = "mm",
  dpi    = figure_dpi,
  bg     = "white",
  device = ragg::agg_png
)
ggsave(
  file.path(dir_figures, "figS13_index_rmse.pdf"),
  plot   = figure,
  width  = figure_width_double,
  height = 118,
  units  = "mm",
  dpi    = figure_dpi,
  bg     = "white",
  device = cairo_pdf
)
