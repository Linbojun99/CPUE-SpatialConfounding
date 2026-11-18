# ==============================================================================
# figures/figS02_logbook_fields.R
#
# Figure S2: estimated annual spatiotemporal random fields of the model
# fitted to the logbook data (Eq. 2). The fields are predicted at the centres
# of the cells of the simulation grid, and each of the nine maps shows one
# year in the cells of the prediction domain. The figure shows no set
# position.
#
# Input:
#   data/simulation_grid.rds                cells of the simulation grid
#   data/prediction_domain_cells.csv        cells of the prediction domain
#   output/logbook_model/logbook_fit.rds    02_fit_logbook_model.R
# Output:
#   output/figures/figS02_logbook_fields.png
#   output/figures/figS02_logbook_fields.pdf
# ==============================================================================

source("00_settings.R")
source("figures/00_figure_style.R")
library(sdmTMB)
library(sf)
library(rnaturalearth)


# 1. Input ---------------------------------------------------------------------

logbook_fit <- readRDS(file.path(dir_logbook, "logbook_fit.rds"))

grid_cells <- readRDS(file.path(dir_data, "simulation_grid.rds"))

domain_table <- fread(file.path(dir_data, "prediction_domain_cells.csv"))
domain_cells <- domain_table[in_domain == TRUE, cell]

coastline <- ne_coastline(scale = 50, returnclass = "sf")


# 2. Spatiotemporal fields at the cell centres ---------------------------------

# One row per cell and year. The fields do not depend on the other terms of
# the model; the prediction is for 1,000 hooks, the reference HBF class and
# no vessel effect.
prediction_grid <- grid_cells[rep(seq_len(nrow(grid_cells)), n_years)]
prediction_grid[, year := rep(years, each = nrow(grid_cells))]
prediction_grid[, year_f := factor(year, levels = years)]
prediction_grid[, hbf_cl := factor("16-17",
                                   levels = c("16-17", "le15", "ge18"))]
prediction_grid[, vessel := NA]
prediction_grid[, log_hooks := 0]

prediction <- predict(
  logbook_fit,
  newdata     = prediction_grid,
  re_form_iid = NA,
  offset      = prediction_grid$log_hooks
)
prediction <- as.data.table(prediction)

field_maps <- prediction[
  cell %in% domain_cells,
  .(cell, year, lon = cell_lon, lat = cell_lat, value = epsilon_st)
]


# 3. Colour scale --------------------------------------------------------------

# Symmetric about zero, with the limit rounded up to the next multiple of 0.5.
colour_limit  <- ceiling(max(abs(field_maps$value)) * 2) / 2
colour_step   <- ifelse(colour_limit > 2, 1, 0.5)
colour_breaks <- seq(-colour_limit, colour_limit, by = colour_step)
colour_labels <- sub("-", minus_sign, format(colour_breaks, trim = TRUE),
                     fixed = TRUE)
cat("Largest absolute value:", max(abs(field_maps$value)),
    " colour limit:", colour_limit, "\n")


# 4. Figure --------------------------------------------------------------------

figure <- ggplot() +
  geom_tile(data = field_maps, aes(lon, lat, fill = value),
            width = 1, height = 1) +
  geom_sf(data = coastline, fill = NA, colour = colour_coastline,
          linewidth = line_width_coastline) +
  coord_sf(xlim = map_xlim, ylim = map_ylim, expand = FALSE) +
  scale_x_continuous(breaks = map_breaks_lon) +
  scale_y_continuous(breaks = map_breaks_lat) +
  scale_fill_gradient2(
    low      = colour_scale_low,
    mid      = "white",
    high     = colour_scale_high,
    midpoint = 0,
    limits   = c(-colour_limit, colour_limit),
    oob      = scales::squish,
    breaks   = colour_breaks,
    labels   = colour_labels,
    name     = "Spatiotemporal effect (log scale)",
    guide    = map_colour_bar
  ) +
  facet_wrap(~year, ncol = 3) +
  labs(x = NULL, y = NULL) +
  theme_figure_map


# 5. Save ----------------------------------------------------------------------

ggsave(
  file.path(dir_figures, "figS02_logbook_fields.png"),
  plot   = figure,
  width  = figure_width_double,
  height = 160,
  units  = "mm",
  dpi    = figure_dpi,
  bg     = "white",
  device = ragg::agg_png
)
ggsave(
  file.path(dir_figures, "figS02_logbook_fields.pdf"),
  plot   = figure,
  width  = figure_width_double,
  height = 160,
  units  = "mm",
  dpi    = figure_dpi,
  bg     = "white",
  device = cairo_pdf
)
