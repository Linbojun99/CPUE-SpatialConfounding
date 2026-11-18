# ==============================================================================
# figures/figS12_true_index.R
#
# Figure S12: true log index over the prediction domain of 1,484 cells
# (Eq. 9) relative to 2015 in A4, E1 and E3: median and 10% and 90%
# quantiles over the 50 replicates. E2 has the same true index as E3. The
# dashed line is the log of the abundance trend tau_t (Eq. 3).
#
# Input:
#   data/prediction_domain_cells.csv               cells of the prediction
#                                                  domain
#   output/simulation/<scenario>/rep<NNNN>.rds     03_simulate_fishery_data.R
#                                                  (A4, E1 and E3, 50
#                                                  replicates)
# Output:
#   output/figures/figS12_true_index.png
#   output/figures/figS12_true_index.pdf
# ==============================================================================

source("00_settings.R")
source("figures/00_figure_style.R")


# 1. Panels --------------------------------------------------------------------

panel_scenarios <- c("A4", "E1", "E3")
panel_titles    <- c("(a) A4 (density, no trend)",
                     "(b) E1 (density, trend)",
                     "(c) E2 and E3 (catchability)")


# 2. True index of each replicate ----------------------------------------------

domain_table <- fread(file.path(dir_data, "prediction_domain_cells.csv"))
domain_cells <- domain_table[in_domain == TRUE, cell]

index_rows <- list()

for (i in seq_along(panel_scenarios)) {
  for (replicate_id in 1:n_replicates) {

    simulated  <- readRDS(file.path(dir_simulation, panel_scenarios[i],
                                    sprintf("rep%04d.rds", replicate_id)))
    cell_years <- as.data.table(simulated$grid)

    # Area-weighted sum of the true density over the cells of the prediction
    # domain (Eq. 9), as the log ratio to 2015.
    true_index <- cell_years[
      cell %in% domain_cells,
      .(log_index = log(sum(area_km2 * exp(log_d)))),
      by = year
    ][order(year)]

    index_rows[[length(index_rows) + 1]] <- data.table(
      panel     = panel_titles[i],
      replicate = replicate_id,
      year      = true_index$year,
      log_ratio = true_index$log_index -
        true_index$log_index[true_index$year == 2015]
    )
  }
}

true_indices <- rbindlist(index_rows)


# 3. Median and quantiles over replicates --------------------------------------

index_summary <- true_indices[, .(
  median = median(log_ratio),
  lo     = quantile(log_ratio, 0.1),
  hi     = quantile(log_ratio, 0.9)
), by = .(panel, year)]
index_summary[, panel := factor(panel, levels = panel_titles)]

# Log of the abundance trend (a 30% decline between 2015 and 2023).
index_summary[, log_trend := log(abundance_trend)[match(year, years)]]


# 4. Figure --------------------------------------------------------------------

figure <- ggplot(index_summary, aes(year, median)) +
  geom_hline(yintercept = 0, colour = colour_reference,
             linewidth = line_width_reference) +
  geom_ribbon(aes(ymin = lo, ymax = hi), fill = colour_simulated,
              alpha = 0.18) +
  geom_line(colour = colour_simulated, linewidth = 0.5) +
  geom_line(aes(y = log_trend), linetype = linetype_dashed,
            linewidth = line_width_reference, colour = "black") +
  facet_wrap(~panel, nrow = 1) +
  scale_x_continuous(breaks = year_breaks, expand = expansion(mult = 0.09)) +
  scale_y_continuous(
    labels = scales::label_number(style_negative = "minus")
  ) +
  labs(x = "Year", y = "Log index relative to 2015") +
  theme_figure +
  theme(panel.spacing.x = unit(4, "mm"))


# 5. Save ----------------------------------------------------------------------

ggsave(
  file.path(dir_figures, "figS12_true_index.png"),
  plot   = figure,
  width  = figure_width_double,
  height = 65,
  units  = "mm",
  dpi    = figure_dpi,
  bg     = "white",
  device = ragg::agg_png
)
ggsave(
  file.path(dir_figures, "figS12_true_index.pdf"),
  plot   = figure,
  width  = figure_width_double,
  height = 65,
  units  = "mm",
  dpi    = figure_dpi,
  bg     = "white",
  device = cairo_pdf
)
