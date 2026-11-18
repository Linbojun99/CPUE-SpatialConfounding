# ==============================================================================
# figures/fig02_simulated_fields.R
#
# Figure 2: simulated spatial fields in replicate 1. Eight maps show the
# spatial field omega, the covariate X1 at three values of rho_omega and two
# values of c_r, and the covariate X2 at two values of rho_12. Each map shows
# the nine-year mean of the field in the cells of the prediction domain,
# centred and scaled with its mean and SD over all cells of the simulation
# grid. Only cell values are read; no set position is used.
#
# Input:
#   data/simulation_grid.rds                   cells of the simulation grid
#   data/prediction_domain_cells.csv           cells of the prediction domain
#   output/simulation/<scenario>/rep0001.rds   03_simulate_fishery_data.R
#                                              (replicate 1 of A1, A3, A4,
#                                              A5, D1 and D2)
# Output:
#   output/figures/fig02_simulated_fields.png
#   output/figures/fig02_simulated_fields.pdf
# ==============================================================================

source("00_settings.R")
source("figures/00_figure_style.R")
library(sf)
library(rnaturalearth)


# 1. Panels --------------------------------------------------------------------

# One row per panel: the scenario and the field that it shows, and the panel
# title as a plotmath expression. Panel (c) (A4: rho_omega = 0.6, c_r = 1) is
# the reference for panels (e) to (h).
panels <- data.table(
  scenario = c("A4", "A1", "A4", "A5", "D1", "D2", "A3", "A4"),
  field    = c("omega", "X1", "X1", "X1", "X1", "X1", "X2", "X2"),
  title    = c(
    paste0("paste('(a) ',", symbol_omega, ")"),
    paste0("paste('(b) ',X[1],'   ',", symbol_rho_omega, ",' = 0')"),
    paste0("paste('(c) ',X[1],'   ',", symbol_rho_omega, ",' = 0.6')"),
    paste0("paste('(d) ',X[1],'   ',", symbol_rho_omega, ",' = 0.9')"),
    "paste('(e) ',X[1],'   ',c[r],' = 0.3')",
    "paste('(f) ',X[1],'   ',c[r],' = 3')",
    paste0("paste('(g) ',X[2],'   ',", symbol_rho_12, ",' = 0')"),
    paste0("paste('(h) ',X[2],'   ',", symbol_rho_12, ",' = 0.8')")
  )
)


# 2. Cells and coastline -------------------------------------------------------

grid_cells <- readRDS(file.path(dir_data, "simulation_grid.rds"))
grid_cells <- grid_cells[, .(cell, lon = cell_lon, lat = cell_lat)]

domain_table <- fread(file.path(dir_data, "prediction_domain_cells.csv"))
domain_cells <- domain_table[in_domain == TRUE, cell]

coastline <- ne_coastline(scale = 50, returnclass = "sf")


# 3. Standardized nine-year mean of each field ---------------------------------

panel_maps <- vector("list", nrow(panels))

for (i in seq_len(nrow(panels))) {

  replicate_file <- file.path(dir_simulation, panels$scenario[i],
                              "rep0001.rds")
  cell_years     <- as.data.table(readRDS(replicate_file)$grid)

  # Mean over the nine years in each cell of the simulation grid.
  cell_means <- cell_years[, .(value = mean(get(panels$field[i]))),
                           by = cell]

  # Centre and scale with the mean and SD over all cells of the grid, then
  # keep the cells of the prediction domain.
  cell_means[, value := (value - mean(value)) / sd(value)]
  cell_means <- merge(cell_means[cell %in% domain_cells], grid_cells,
                      by = "cell")
  cell_means[, panel := panels$title[i]]

  panel_maps[[i]] <- cell_means
}

panel_maps <- rbindlist(panel_maps)
panel_maps[, panel := factor(panel, levels = panels$title)]

# Number of displayed values beyond the limit of the colour scale (stated in
# the caption).
cat("Displayed values:", nrow(panel_maps), "\n")
cat("Values beyond the colour limit:",
    sum(abs(panel_maps$value) > map_colour_limit), "\n")


# 4. Figure --------------------------------------------------------------------

figure <- ggplot() +
  geom_tile(data = panel_maps, aes(lon, lat, fill = value),
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
    limits   = c(-map_colour_limit, map_colour_limit),
    oob      = scales::squish,
    breaks   = map_colour_breaks,
    labels   = map_colour_labels,
    name     = "Standardized nine-year mean",
    guide    = map_colour_bar
  ) +
  facet_wrap(~panel, ncol = 4, labeller = label_parsed) +
  labs(x = NULL, y = NULL) +
  theme_figure_map


# 5. Save ----------------------------------------------------------------------

ggsave(
  file.path(dir_figures, "fig02_simulated_fields.png"),
  plot   = figure,
  width  = figure_width_double,
  height = 100,
  units  = "mm",
  dpi    = figure_dpi,
  bg     = "white",
  device = ragg::agg_png
)
ggsave(
  file.path(dir_figures, "fig02_simulated_fields.pdf"),
  plot   = figure,
  width  = figure_width_double,
  height = 100,
  units  = "mm",
  dpi    = figure_dpi,
  bg     = "white",
  device = cairo_pdf
)
