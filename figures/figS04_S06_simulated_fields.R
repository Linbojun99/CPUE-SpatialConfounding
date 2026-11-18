# ==============================================================================
# figures/figS04_S06_simulated_fields.R
#
# Figures S4, S5 and S6: the spatial field omega and the covariates X1 and
# X2 in each year of replicate 1 of scenario A4 (rho_omega = 0.6,
# rho_12 = 0.8, c_r = 1). Each figure shows three years (rows) and the three
# fields (columns) in the cells of the prediction domain. Each field is
# standardized with its mean and SD over all cells of the simulation grid
# and all nine years. Only cell values are read; no set position is used.
#
# Input:
#   data/simulation_grid.rds               cells of the simulation grid
#   data/prediction_domain_cells.csv       cells of the prediction domain
#   output/simulation/A4/rep0001.rds       03_simulate_fishery_data.R
# Output:
#   output/figures/figS04_simulated_fields_2015_2017.png and .pdf
#   output/figures/figS05_simulated_fields_2018_2020.png and .pdf
#   output/figures/figS06_simulated_fields_2021_2023.png and .pdf
# ==============================================================================

source("00_settings.R")
source("figures/00_figure_style.R")
library(sf)
library(rnaturalearth)


# 1. Input ---------------------------------------------------------------------

grid_cells <- readRDS(file.path(dir_data, "simulation_grid.rds"))
grid_cells <- grid_cells[, .(cell, lon = cell_lon, lat = cell_lat)]

domain_table <- fread(file.path(dir_data, "prediction_domain_cells.csv"))
domain_cells <- domain_table[in_domain == TRUE, cell]

cell_years <- readRDS(file.path(dir_simulation, "A4", "rep0001.rds"))$grid
cell_years <- as.data.table(cell_years)

coastline <- ne_coastline(scale = 50, returnclass = "sf")


# 2. Standardized fields -------------------------------------------------------

# Column of the simulated data and column title (plotmath) of each field.
field_columns <- c("omega", "X1", "X2")
field_titles  <- c(symbol_omega, "X[1]", "X[2]")

field_maps <- list()

for (i in seq_along(field_columns)) {

  values <- cell_years[[field_columns[i]]]

  field_maps[[i]] <- data.table(
    cell  = cell_years$cell,
    year  = cell_years$year,
    field = field_titles[i],
    value = (values - mean(values)) / sd(values)
  )
}

field_maps <- rbindlist(field_maps)
field_maps <- merge(field_maps[cell %in% domain_cells], grid_cells,
                    by = "cell")
field_maps[, field := factor(field, levels = field_titles)]

# Number of displayed values beyond the limit of the colour scale (stated in
# the caption of Fig. S4).
cat("Displayed values:", nrow(field_maps), "\n")
cat("Values beyond the colour limit:",
    sum(abs(field_maps$value) > map_colour_limit), "\n")


# 3. One figure for each group of three years ----------------------------------

figure_years <- list(2015:2017, 2018:2020, 2021:2023)
figure_names <- c("figS04_simulated_fields_2015_2017",
                  "figS05_simulated_fields_2018_2020",
                  "figS06_simulated_fields_2021_2023")

for (i in seq_along(figure_years)) {

  figure <- ggplot() +
    geom_tile(data = field_maps[year %in% figure_years[[i]]],
              aes(lon, lat, fill = value), width = 1, height = 1) +
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
      limits   = c(-map_colour_limit, map_colour_limit),
      oob      = scales::squish,
      breaks   = map_colour_breaks,
      labels   = map_colour_labels,
      name     = "Standardized value",
      guide    = map_colour_bar
    ) +
    facet_grid(year ~ field, labeller = labeller(field = label_parsed)) +
    labs(x = NULL, y = NULL) +
    theme_figure_map +
    theme(strip.text.y = element_text(angle = 0, hjust = 0))

  ggsave(
    file.path(dir_figures, paste0(figure_names[i], ".png")),
    plot   = figure,
    width  = figure_width_double,
    height = 160,
    units  = "mm",
    dpi    = figure_dpi,
    bg     = "white",
    device = ragg::agg_png
  )
  ggsave(
    file.path(dir_figures, paste0(figure_names[i], ".pdf")),
    plot   = figure,
    width  = figure_width_double,
    height = 160,
    units  = "mm",
    dpi    = figure_dpi,
    bg     = "white",
    device = cairo_pdf
  )
}
