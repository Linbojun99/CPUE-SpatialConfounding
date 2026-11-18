# ==============================================================================
# figures/figS11_realized_correlations.R
#
# Figure S11: realized correlations in the simulated data of each scenario,
# over the 50 replicates.
#   (a)  between X1 and omega over the cells of the simulation grid
#   (b)  between X1 and X2 at the simulated sets
#   (c)  between X1 and omega at the simulated sets
# Dashed lines show the design values rho_omega (a, c) and rho_12 (b) of
# Table 2. The correlations are those stored with each simulated data set.
#
# Input:
#   output/simulation/<scenario>/rep<NNNN>.rds     03_simulate_fishery_data.R
#                                                  (17 scenarios, 50
#                                                  replicates)
# Output:
#   output/figures/figS11_realized_correlations.png
#   output/figures/figS11_realized_correlations.pdf
# ==============================================================================

source("00_settings.R")
source("figures/00_figure_style.R")


# 1. Panels --------------------------------------------------------------------

# One row per panel: where the correlation is calculated ("grid_equal": over
# all cells and years of the grid with equal weights; "sets": at the
# simulated sets), the pair of variables, the column of the scenario table
# with the design value, and the panel title as a plotmath expression.
panels <- data.table(
  version       = c("grid_equal", "sets", "sets"),
  pair          = c("X1-omega", "X1-X2", "X1-omega"),
  design_column = c("rho_omega", "rho_12", "rho_omega"),
  title         = c(
    paste0("paste('(a) ',X[1],' and ',", symbol_omega,
           ",', simulation grid')"),
    "paste('(b) ',X[1],' and ',X[2],', simulated sets')",
    paste0("paste('(c) ',X[1],' and ',", symbol_omega,
           ",', simulated sets')")
  )
)


# 2. Correlations of each replicate --------------------------------------------

correlation_rows <- list()

for (scenario_id in scenario_order) {
  for (replicate_id in 1:n_replicates) {

    simulated <- readRDS(file.path(dir_simulation, scenario_id,
                                   sprintf("rep%04d.rds", replicate_id)))
    stored    <- as.data.table(simulated$gen_diag)[stat == "cor"]

    for (i in seq_len(nrow(panels))) {
      correlation_rows[[length(correlation_rows) + 1]] <- data.table(
        scenario = scenario_id,
        panel    = panels$title[i],
        value    = stored[version == panels$version[i] &
                            var == panels$pair[i], value]
      )
    }
  }
  cat("Read scenario", scenario_id, "\n")
}

correlations <- rbindlist(correlation_rows)
correlations[, scenario := factor(scenario, levels = scenario_order)]
correlations[, panel := factor(panel, levels = panels$title)]


# 3. Design values -------------------------------------------------------------

# Rows of the scenario table for the 17 scenarios of Table 2.
scenario_design <- scenarios[match(scenario_order, scenario)]

design_rows <- list()

for (i in seq_len(nrow(panels))) {
  design_rows[[i]] <- data.table(
    scenario = scenario_design$scenario,
    panel    = panels$title[i],
    value    = scenario_design[[panels$design_column[i]]]
  )
}

design_values <- rbindlist(design_rows)
design_values[, scenario := factor(scenario, levels = scenario_order)]
design_values[, panel := factor(panel, levels = panels$title)]


# 4. Figure --------------------------------------------------------------------

# Boxes show the median and interquartile range and whiskers the most
# extreme values within 1.5 times the interquartile range; outliers are not
# drawn.
figure <- ggplot(correlations, aes(scenario, value)) +
  geom_boxplot(fill = colour_simulated, alpha = 0.75,
               colour = colour_box_outline, width = 0.65,
               outlier.shape = NA, linewidth = 0.3) +
  geom_errorbar(data = design_values, aes(ymin = value, ymax = value),
                width = 0.8, linetype = linetype_dashed,
                linewidth = line_width_reference, colour = "black") +
  facet_wrap(~panel, ncol = 1, labeller = label_parsed) +
  scale_y_continuous(
    labels = scales::label_number(style_negative = "minus")
  ) +
  coord_cartesian(ylim = c(-0.3, 1)) +
  labs(x = "Scenario", y = "Correlation") +
  theme_figure


# 5. Save ----------------------------------------------------------------------

ggsave(
  file.path(dir_figures, "figS11_realized_correlations.png"),
  plot   = figure,
  width  = figure_width_double,
  height = 110,
  units  = "mm",
  dpi    = figure_dpi,
  bg     = "white",
  device = ragg::agg_png
)
ggsave(
  file.path(dir_figures, "figS11_realized_correlations.pdf"),
  plot   = figure,
  width  = figure_width_double,
  height = 110,
  units  = "mm",
  dpi    = figure_dpi,
  bg     = "white",
  device = cairo_pdf
)
