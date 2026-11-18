# ==============================================================================
# figures/figS01_logbook_residuals.R
#
# Figure S1: residual diagnostics of the two models fitted to the logbook
# data, the model in Eq. 2 with spatiotemporal fields that are independent
# among years (IID) and the same model with first-order autoregressive
# fields (AR1). Randomized quantile residuals are compared with standard
# normal quantiles (top row) and summarized by their median and 10% and 90%
# quantiles in 20 equally populated groups of the fitted linear predictor
# (middle row) and by year (bottom row). The figure shows no set position.
#
# Input:
#   output/logbook_model/logbook_fit.rds       02_fit_logbook_model.R
#   output/logbook_model/logbook_fit_ar1.rds   02_fit_logbook_model.R
# Output:
#   output/figures/figS01_logbook_residuals.png
#   output/figures/figS01_logbook_residuals.pdf
#
# Run time: about 1 minute.
# ==============================================================================

source("00_settings.R")
source("figures/00_figure_style.R")
library(sdmTMB)


# 1. Settings of the figure ----------------------------------------------------

# The randomized quantile residuals depend on random draws; the same seed is
# set before the residuals of each model are calculated.
residual_seed <- 20261023

# Number of groups of the fitted linear predictor and number of points drawn
# in each quantile-quantile panel.
n_groups    <- 20
n_qq_points <- 2000

model_names <- c("IID", "AR1")
model_files <- c(IID = "logbook_fit.rds", AR1 = "logbook_fit_ar1.rds")


# 2. Residuals of the two models -----------------------------------------------

qq_rows    <- list()
group_rows <- list()
year_rows  <- list()

for (model_name in model_names) {

  logbook_fit <- readRDS(file.path(dir_logbook, model_files[[model_name]]))

  # Randomized quantile residuals with the fixed effects at their estimates
  # and the random effects drawn from their approximate distribution.
  set.seed(residual_seed)
  residual <- residuals(logbook_fit, type = "mle-mvn")

  set_residuals <- data.table(
    model    = model_name,
    residual = residual,
    fitted   = predict(logbook_fit)$est,
    year     = logbook_fit$data$year
  )

  # Sorted residuals against standard normal quantiles. Of the 41,035
  # points, 2,000 equally spaced ones are drawn.
  qq_points <- data.table(
    model       = model_name,
    theoretical = qnorm(ppoints(nrow(set_residuals))),
    residual    = sort(set_residuals$residual)
  )
  qq_rows[[model_name]] <- qq_points[
    unique(round(seq(1, nrow(qq_points), length.out = n_qq_points)))
  ]

  # Groups with equal numbers of sets, ordered by the fitted linear predictor.
  set_residuals[, group := pmin(
    n_groups,
    ceiling(frank(fitted, ties.method = "first") / .N * n_groups)
  )]
  group_rows[[model_name]] <- set_residuals[, .(
    fitted = mean(fitted),
    median = median(residual),
    lo     = quantile(residual, 0.1),
    hi     = quantile(residual, 0.9)
  ), by = .(model, group)]

  year_rows[[model_name]] <- set_residuals[, .(
    median = median(residual),
    lo     = quantile(residual, 0.1),
    hi     = quantile(residual, 0.9)
  ), by = .(model, year)]

  rm(logbook_fit, set_residuals)
}

qq_points       <- rbindlist(qq_rows)
group_summaries <- rbindlist(group_rows)
year_summaries  <- rbindlist(year_rows)

qq_points[, model := factor(model, levels = model_names)]
group_summaries[, model := factor(model, levels = model_names)]
year_summaries[, model := factor(model, levels = model_names)]


# 3. Panels --------------------------------------------------------------------

model_colours <- c(IID = colour_scale_low, AR1 = colour_scale_high)

theme_panel <- theme_figure +
  theme(legend.position = "none", panel.spacing.x = unit(5, "mm"))

# (a), (b) Quantile-quantile panels with the one-to-one line.
panel_qq <- ggplot(qq_points, aes(theoretical, residual, colour = model)) +
  geom_abline(slope = 1, intercept = 0, linetype = linetype_dashed,
              linewidth = line_width_reference, colour = "black") +
  geom_point(size = 0.3) +
  facet_wrap(
    ~model,
    nrow     = 1,
    labeller = labeller(model = c(IID = "(a) IID", AR1 = "(b) AR1"))
  ) +
  scale_colour_manual(values = model_colours) +
  scale_x_continuous(
    labels = scales::label_number(style_negative = "minus")
  ) +
  scale_y_continuous(
    labels = scales::label_number(style_negative = "minus")
  ) +
  labs(x = "Theoretical normal quantile", y = "Residual quantile") +
  theme_panel

# (c), (d) Residuals by group of the fitted linear predictor.
panel_fitted <- ggplot(
  group_summaries,
  aes(fitted, median, colour = model, fill = model)
) +
  geom_hline(yintercept = 0, colour = colour_reference,
             linewidth = line_width_reference) +
  geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.18, colour = NA) +
  geom_line(linewidth = 0.4) +
  facet_wrap(
    ~model,
    nrow     = 1,
    labeller = labeller(model = c(IID = "(c) IID", AR1 = "(d) AR1"))
  ) +
  scale_colour_manual(values = model_colours) +
  scale_fill_manual(values = model_colours) +
  scale_y_continuous(
    labels = scales::label_number(style_negative = "minus")
  ) +
  labs(x = "Fitted linear predictor", y = "Residual") +
  theme_panel

# (e), (f) Residuals by year.
panel_year <- ggplot(
  year_summaries,
  aes(year, median, colour = model, fill = model)
) +
  geom_hline(yintercept = 0, colour = colour_reference,
             linewidth = line_width_reference) +
  geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.18, colour = NA) +
  geom_line(linewidth = 0.4) +
  facet_wrap(
    ~model,
    nrow     = 1,
    labeller = labeller(model = c(IID = "(e) IID", AR1 = "(f) AR1"))
  ) +
  scale_colour_manual(values = model_colours) +
  scale_fill_manual(values = model_colours) +
  scale_x_continuous(breaks = year_breaks, expand = expansion(mult = 0.04)) +
  scale_y_continuous(
    labels = scales::label_number(style_negative = "minus")
  ) +
  labs(x = "Year", y = "Residual") +
  theme_panel


# 4. Combine and save ----------------------------------------------------------

figure <- panel_qq / panel_fitted / panel_year +
  plot_layout(heights = c(1.25, 1, 1))

ggsave(
  file.path(dir_figures, "figS01_logbook_residuals.png"),
  plot   = figure,
  width  = figure_width_double,
  height = 190,
  units  = "mm",
  dpi    = figure_dpi,
  bg     = "white",
  device = ragg::agg_png
)
ggsave(
  file.path(dir_figures, "figS01_logbook_residuals.pdf"),
  plot   = figure,
  width  = figure_width_double,
  height = 190,
  units  = "mm",
  dpi    = figure_dpi,
  bg     = "white",
  device = cairo_pdf
)
