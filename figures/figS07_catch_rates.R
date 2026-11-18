# ==============================================================================
# figures/figS07_catch_rates.R
#
# Figure S7: catch rates in the logbook and in the simulated data.
#   (a)  catch per 1,000 hooks in sets with a catch greater than zero, in
#        the logbook and in the simulated sets of A4 (all replicates)
#   (b)  quantiles of the set-level catch rate in each scenario, pooled over
#        the replicates, with the percentage of sets with zero catch
#   (c)  annual nominal catch rate (total catch divided by total number of
#        hooks) relative to 2015 on the log scale: median and 10% and 90%
#        quantiles over replicates, with the median true index
# The true index in (c) is the one stored with the simulated data, summed
# over all 1,592 cells of the simulation grid. Only catches, numbers of
# hooks and years are read; no set position is used.
#
# Input:
#   data/logbook_sets.rds                          41,035 logbook sets
#   output/simulation/<scenario>/rep<NNNN>.rds     03_simulate_fishery_data.R
#                                                  (17 scenarios, 50
#                                                  replicates)
# Output:
#   output/figures/figS07_catch_rates.png
#   output/figures/figS07_catch_rates.pdf
# ==============================================================================

source("00_settings.R")
source("figures/00_figure_style.R")


# 1. Catch rates in the logbook ------------------------------------------------

logbook_sets  <- readRDS(file.path(dir_data, "logbook_sets.rds"))
logbook_rates <- logbook_sets[, .(catch_rate = bet_n / hooks * 1000)]


# 2. Catch rates in the simulated data -----------------------------------------

scenario_quantiles <- list()   # (b) one row per scenario
annual_rates       <- list()   # (c) one row per scenario, replicate and year
simulated_rates_a4 <- list()   # (a) set-level catch rates of A4

for (scenario_id in scenario_order) {

  scenario_rates <- vector("list", n_replicates)

  for (replicate_id in 1:n_replicates) {

    simulated <- readRDS(file.path(dir_simulation, scenario_id,
                                   sprintf("rep%04d.rds", replicate_id)))

    # Set-level catch per 1,000 hooks.
    sets <- simulated$sets[, .(year, hooks, bet_n)]
    sets[, catch_rate := bet_n / hooks * 1000]
    scenario_rates[[replicate_id]] <- sets[, .(catch_rate)]

    # Annual nominal catch rate and true index, as log ratios to 2015.
    annual <- sets[, .(nominal = sum(bet_n) / sum(hooks) * 1000), by = year]
    annual <- annual[order(year)]
    truth  <- as.data.table(simulated$truth)[order(year)]

    annual_rates[[paste(scenario_id, replicate_id)]] <- data.table(
      scenario   = scenario_id,
      replicate  = replicate_id,
      year       = annual$year,
      nominal    = log(annual$nominal) -
        log(annual$nominal[annual$year == 2015]),
      true_index = log(truth$true_index) -
        log(truth$true_index[truth$year == 2015])
    )
  }

  scenario_rates <- rbindlist(scenario_rates)
  if (scenario_id == "A4") {
    simulated_rates_a4 <- scenario_rates
  }

  # Quantiles of the catch rates greater than zero and the percentage of
  # sets with zero catch, over all replicates.
  positive_quantiles <- unname(quantile(
    scenario_rates$catch_rate[scenario_rates$catch_rate > 0],
    c(0.05, 0.25, 0.5, 0.75, 0.95)
  ))
  scenario_quantiles[[scenario_id]] <- data.table(
    scenario     = scenario_id,
    zero_percent = 100 * mean(scenario_rates$catch_rate == 0),
    q05          = positive_quantiles[1],
    q25          = positive_quantiles[2],
    q50          = positive_quantiles[3],
    q75          = positive_quantiles[4],
    q95          = positive_quantiles[5]
  )
  cat("Read scenario", scenario_id, "\n")
}

scenario_quantiles <- rbindlist(scenario_quantiles)
annual_rates       <- rbindlist(annual_rates)

scenario_quantiles[, scenario := factor(scenario, levels = scenario_order)]


# 3. Panel (a): logbook and simulated sets of A4 -------------------------------

source_names <- c(
  paste0("Logbook (", format(nrow(logbook_rates), big.mark = ","), " sets)"),
  paste0("Simulated (A4, ", n_replicates, " replicates)")
)

logbook_rates[, source := source_names[1]]
simulated_rates_a4[, source := source_names[2]]
set_rates <- rbind(logbook_rates, simulated_rates_a4)
set_rates[, source := factor(source, levels = source_names)]

# Median catch rate and percentage of sets with zero catch, over all sets.
source_summary <- set_rates[, .(
  zero_percent = 100 * mean(catch_rate == 0),
  median       = median(catch_rate)
), by = source][order(source)]

summary_text <- paste(
  sprintf("%s, median %.1f, zero catch %s%%",
          c("Logbook", "Simulated"),
          source_summary$median,
          formatC(source_summary$zero_percent, format = "fg", digits = 2)),
  collapse = "\n"
)

panel_histogram <- ggplot(
  set_rates[catch_rate > 0],
  aes(catch_rate, fill = source)
) +
  geom_histogram(aes(y = after_stat(density)),
                 breaks = 10^seq(-1.5, 2.5, length.out = 41),
                 position = "identity", alpha = 0.6, colour = "black",
                 linewidth = 0.15) +
  annotate("text", x = Inf, y = Inf, label = summary_text, hjust = 1.03,
           vjust = 1.3, size = 7 / points_per_mm, family = font_family,
           lineheight = 1) +
  scale_fill_manual(
    values = setNames(c(colour_neutral, colour_simulated), source_names)
  ) +
  scale_x_log10(labels = scales::label_number(drop0trailing = TRUE)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.25))) +
  labs(x = "Catch per 1,000 hooks (log scale)", y = "Density",
       title = "(a) Sets with a catch greater than zero") +
  theme_figure +
  guides(fill = guide_legend(nrow = 2))


# 4. Panel (b): set-level catch rate by scenario -------------------------------

# Boxes show the 25%, 50% and 75% quantiles and whiskers the 5% and 95%
# quantiles of the catch rates greater than zero. The number above each box
# is the percentage of sets with zero catch.
panel_scenarios <- ggplot(scenario_quantiles, aes(x = scenario)) +
  geom_boxplot(aes(ymin = q05, lower = q25, middle = q50, upper = q75,
                   ymax = q95),
               stat = "identity", fill = colour_simulated, alpha = 0.5,
               colour = colour_box_outline, linewidth = 0.3, width = 0.6) +
  geom_text(aes(y = q95,
                label = formatC(zero_percent, format = "f", digits = 1)),
            vjust = -0.5, size = 6 / points_per_mm, family = font_family) +
  scale_y_log10(labels = scales::label_number(drop0trailing = TRUE),
                expand = expansion(mult = c(0.05, 0.18))) +
  labs(x = NULL, y = "Catch per 1,000 hooks (log scale)",
       title = "(b) Simulated sets by scenario") +
  theme_figure +
  theme(axis.text.x = element_text(size = 7))


# 5. Panel (c): annual nominal catch rate and true index -----------------------

annual_summary <- annual_rates[, .(
  median     = median(nominal),
  lo         = quantile(nominal, 0.1),
  hi         = quantile(nominal, 0.9),
  true_index = median(true_index)
), by = .(scenario, year)]
annual_summary[, scenario := factor(scenario, levels = scenario_order)]

# Colour and line type identify the two lines; the two scales have the same
# (empty) name, so that they share one legend.
line_names <- c("Nominal catch rate", "True index")

panel_annual <- ggplot(annual_summary, aes(year)) +
  geom_hline(yintercept = 0, colour = colour_reference,
             linewidth = line_width_reference) +
  geom_ribbon(aes(ymin = lo, ymax = hi), fill = colour_simulated,
              alpha = 0.25) +
  geom_line(aes(y = median, colour = "Nominal catch rate",
                linetype = "Nominal catch rate"), linewidth = 0.5) +
  geom_line(aes(y = true_index, colour = "True index",
                linetype = "True index"), linewidth = 0.4) +
  facet_wrap(~scenario, ncol = 6) +
  scale_colour_manual(
    name   = NULL,
    values = setNames(c(colour_simulated, "black"), line_names)
  ) +
  scale_linetype_manual(
    name   = NULL,
    values = setNames(c("solid", linetype_dashed), line_names)
  ) +
  scale_x_continuous(breaks = year_breaks, labels = year_labels_narrow,
                     expand = expansion(mult = 0.07)) +
  scale_y_continuous(
    labels = scales::label_number(style_negative = "minus")
  ) +
  labs(x = "Year", y = "Log ratio to 2015",
       title = "(c) Annual nominal catch rate and true index by scenario") +
  theme_figure +
  theme(
    # The labels of the first and last year are aligned inward, so that the
    # labels of adjacent panels do not touch.
    axis.text.x      = suppressWarnings(
      element_text(size = 7, hjust = c(0, 0.5, 1))
    ),
    legend.key.width = unit(8, "mm"),
    panel.spacing.x  = unit(3, "mm")
  )


# 6. Combine and save ----------------------------------------------------------

figure <- (panel_histogram + panel_scenarios +
             plot_layout(widths = c(1, 1.1))) /
  panel_annual +
  plot_layout(heights = c(1, 1.5))

ggsave(
  file.path(dir_figures, "figS07_catch_rates.png"),
  plot   = figure,
  width  = figure_width_double,
  height = 160,
  units  = "mm",
  dpi    = figure_dpi,
  bg     = "white",
  device = ragg::agg_png
)
ggsave(
  file.path(dir_figures, "figS07_catch_rates.pdf"),
  plot   = figure,
  width  = figure_width_double,
  height = 160,
  units  = "mm",
  dpi    = figure_dpi,
  bg     = "white",
  device = cairo_pdf
)
