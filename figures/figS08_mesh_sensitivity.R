# ==============================================================================
# figures/figS08_mesh_sensitivity.R
#
# Figure S8: sensitivity of the coefficient of X1 to the estimation mesh
# (Table S1). The coarser mesh (620 vertices) and the finer mesh (3,653
# vertices) are compared with the selected mesh (1,676 vertices) in the
# preliminary replicates of A2 and D1: 10 replicates for ST (X1) and
# ST (X1, X2) on the coarser mesh and 3 replicates for ST (X1, X2) on the
# finer mesh. The top row shows the estimate on the comparison mesh minus
# that on the selected mesh, the bottom row the ratio of the estimated
# standard errors. Points are means over replicates and vertical segments
# the full range.
#
# Input:
#   output/mesh_comparison/<mesh>/<scenario>/<model>_rep<NNNN>.rds
#     fits of 04_fit_estimators.R with the meshes "selected", "coarser" and
#     "finer" (see section 1 of that script)
#   output/meshes/mesh_estimation_coarser.rds   01_build_meshes.R (number of
#   output/meshes/mesh_estimation_finer.rds     vertices)
# Output:
#   output/figures/figS08_mesh_sensitivity.png
#   output/figures/figS08_mesh_sensitivity.pdf
# ==============================================================================

source("00_settings.R")
source("figures/00_figure_style.R")


# 1. Fits that are compared ----------------------------------------------------

comparison_meshes     <- c("coarser", "finer")
comparison_scenarios  <- c("A2", "D1")
comparison_models     <- c("M3", "M4")
comparison_replicates <- 9001:9010

# Number of vertices of the comparison meshes, written with a thousands
# separator for the axis.
mesh_coarser <- readRDS(file.path(dir_meshes, "mesh_estimation_coarser.rds"))
mesh_finer   <- readRDS(file.path(dir_meshes, "mesh_estimation_finer.rds"))
mesh_labels  <- c(
  coarser = format(mesh_coarser$n, big.mark = ","),
  finer   = format(mesh_finer$n, big.mark = ",")
)


# 2. Coefficient of X1 on the comparison mesh and on the selected mesh ---------

comparison_rows <- list()

for (mesh_name in comparison_meshes) {
  for (scenario_id in comparison_scenarios) {
    for (model_id in comparison_models) {
      for (replicate_id in comparison_replicates) {

        # Not every model and replicate was fitted on the finer mesh.
        file_name <- sprintf("%s_rep%04d.rds", model_id, replicate_id)
        file_comparison <- file.path(dir_mesh_runs, mesh_name,
                                     scenario_id, file_name)
        if (!file.exists(file_comparison)) {
          next
        }
        fit_comparison <- readRDS(file_comparison)
        fit_selected   <- readRDS(file.path(dir_mesh_runs, "selected",
                                            scenario_id, file_name))
        x1_comparison  <- fit_comparison$coef[term == "X1"]
        x1_selected    <- fit_selected$coef[term == "X1"]

        comparison_rows[[length(comparison_rows) + 1]] <- data.table(
          mesh        = mesh_name,
          scenario_id = scenario_id,
          model_id    = model_id,
          replicate   = replicate_id,
          coefficient_difference =
            x1_comparison$estimate - x1_selected$estimate,
          se_ratio    = x1_comparison$std.error / x1_selected$std.error
        )
      }
    }
  }
}

comparisons <- rbindlist(comparison_rows)


# 3. Mean and range over replicates --------------------------------------------

mesh_summary <- rbind(
  comparisons[, .(
    measure = "Coefficient difference",
    mean    = mean(coefficient_difference),
    min     = min(coefficient_difference),
    max     = max(coefficient_difference),
    n       = .N
  ), by = .(scenario_id, model_id, mesh)],
  comparisons[, .(
    measure = "SE ratio",
    mean    = mean(se_ratio),
    min     = min(se_ratio),
    max     = max(se_ratio),
    n       = .N
  ), by = .(scenario_id, model_id, mesh)]
)
mesh_summary[, model := factor(model_id, levels = comparison_models)]
mesh_summary[, vertices := factor(mesh_labels[mesh],
                                  levels = mesh_labels[comparison_meshes])]


# 4. Figure --------------------------------------------------------------------

# One row of panels for each measure. The grey line marks no difference: a
# coefficient difference of 0 and an SE ratio of 1.
measures         <- c("Coefficient difference", "SE ratio")
reference_values <- c(0, 1)
panel_letters    <- list(c("(a)", "(b)"), c("(c)", "(d)"))

panels <- list()
for (i in seq_along(measures)) {
  panel_data <- mesh_summary[measure == measures[i]]
  panel_data[, panel := factor(
    scenario_id,
    levels = comparison_scenarios,
    labels = paste(panel_letters[[i]], comparison_scenarios)
  )]

  panels[[i]] <- ggplot(
    panel_data,
    aes(vertices, mean, colour = model, fill = model, shape = model)
  ) +
    geom_hline(yintercept = reference_values[i], colour = colour_reference,
               linewidth = line_width_reference) +
    geom_linerange(aes(ymin = min, ymax = max),
                   position = position_dodge(0.45), linewidth = 0.4,
                   show.legend = FALSE) +
    geom_point(position = position_dodge(0.45), size = 1.9, stroke = 0.5) +
    facet_wrap(~panel, nrow = 1) +
    scale_colour_manual(
      values = estimator_colours,
      breaks = comparison_models,
      labels = parse(text = unname(estimator_labels[comparison_models])),
      limits = comparison_models,
      drop   = FALSE
    ) +
    scale_fill_manual(
      values = estimator_colours,
      breaks = comparison_models,
      labels = parse(text = unname(estimator_labels[comparison_models])),
      limits = comparison_models,
      drop   = FALSE
    ) +
    scale_shape_manual(
      values = estimator_shapes,
      breaks = comparison_models,
      labels = parse(text = unname(estimator_labels[comparison_models])),
      limits = comparison_models,
      drop   = FALSE
    ) +
    guides(
      colour = guide_legend(
        nrow = 1, byrow = TRUE,
        override.aes = list(
          fill = unname(estimator_colours[comparison_models])
        )
      ),
      shape = guide_legend(nrow = 1, byrow = TRUE),
      fill  = "none"
    ) +
    scale_y_continuous(
      labels = scales::label_number(style_negative = "minus")
    ) +
    labs(x = "Vertices in the comparison mesh", y = measures[i]) +
    theme_figure +
    theme(panel.spacing.x = unit(4, "mm"))
}

figure <- panels[[1]] / panels[[2]] +
  plot_layout(guides = "collect", axis_titles = "collect") &
  theme(legend.position = "bottom")


# 5. Save ----------------------------------------------------------------------

ggsave(
  file.path(dir_figures, "figS08_mesh_sensitivity.png"),
  plot   = figure,
  width  = figure_width_double,
  height = 112,
  units  = "mm",
  dpi    = figure_dpi,
  bg     = "white",
  device = ragg::agg_png
)
ggsave(
  file.path(dir_figures, "figS08_mesh_sensitivity.pdf"),
  plot   = figure,
  width  = figure_width_double,
  height = 112,
  units  = "mm",
  dpi    = figure_dpi,
  bg     = "white",
  device = cairo_pdf
)
