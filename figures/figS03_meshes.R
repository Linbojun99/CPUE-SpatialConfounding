# ==============================================================================
# figures/figS03_meshes.R
#
# Figure S3: triangular meshes of the operating model and of the estimation
# models (Table S1). Panels (a) and (b) show the whole meshes with the
# centres of the cells of the simulation grid; panels (c) and (d) enlarge
# the area marked by a rectangle in (a) and (b).
#
# Input:
#   data/simulation_grid.rds              cells of the simulation grid
#   output/meshes/mesh_simulation.rds     01_build_meshes.R
#   output/meshes/mesh_estimation.rds     01_build_meshes.R
# Output:
#   output/figures/figS03_meshes.png
#   output/figures/figS03_meshes.pdf
# ==============================================================================

source("00_settings.R")
source("figures/00_figure_style.R")


# 1. Input ---------------------------------------------------------------------

grid_cells   <- readRDS(file.path(dir_data, "simulation_grid.rds"))
cell_centres <- unique(grid_cells[, .(x, y)])

meshes <- list(
  readRDS(file.path(dir_meshes, "mesh_simulation.rds"))$mesh,
  readRDS(file.path(dir_meshes, "mesh_estimation.rds"))
)
mesh_names <- c("Operating model", "Estimation models")


# 2. Edges of the triangles ----------------------------------------------------

# Each triangle has three edges; an edge shared by two triangles is kept
# once. The edges are stored by the coordinates (km) of their two ends.
mesh_edges <- list()

for (i in seq_along(meshes)) {

  triangles <- meshes[[i]]$graph$tv
  vertices  <- meshes[[i]]$loc

  edge_ends <- rbind(triangles[, 1:2], triangles[, 2:3],
                     triangles[, c(3, 1)])
  edge_ends <- unique(cbind(pmin(edge_ends[, 1], edge_ends[, 2]),
                            pmax(edge_ends[, 1], edge_ends[, 2])))

  mesh_edges[[i]] <- data.table(
    x     = vertices[edge_ends[, 1], 1],
    y     = vertices[edge_ends[, 1], 2],
    x_end = vertices[edge_ends[, 2], 1],
    y_end = vertices[edge_ends[, 2], 2]
  )
}


# 3. Extent of the panels and axis labels --------------------------------------

# The overview panels cover both meshes. The enlarged panels show a square
# of 1,400 km centred on the median coordinates of the cell centres.
overview_xlim <- range(c(meshes[[1]]$loc[, 1], meshes[[2]]$loc[, 1]))
overview_ylim <- range(c(meshes[[1]]$loc[, 2], meshes[[2]]$loc[, 2]))

half_width    <- 700
centre_x      <- median(cell_centres$x)
centre_y      <- median(cell_centres$y)
enlarged_xlim <- c(centre_x - half_width, centre_x + half_width)
enlarged_ylim <- c(centre_y - half_width, centre_y + half_width)

# Axis labels with a thousands separator and a minus sign. The outer
# format() pads the labels with spaces to a common width, as in the
# published figure.
overview_breaks <- c(-6000, -3000, 0, 3000, 6000)
overview_labels <- format(sub(
  "-", minus_sign, format(overview_breaks, big.mark = ",", trim = TRUE),
  fixed = TRUE
))
enlarged_breaks <- c(-500, 0, 500)
enlarged_labels <- format(sub(
  "-", minus_sign, format(enlarged_breaks, big.mark = ",", trim = TRUE),
  fixed = TRUE
))


# 4. Panels --------------------------------------------------------------------

panel_titles <- c(
  paste0("(a) ", mesh_names[1], " (",
         format(meshes[[1]]$n, big.mark = ","), " vertices)"),
  paste0("(b) ", mesh_names[2], " (",
         format(meshes[[2]]$n, big.mark = ","), " vertices)"),
  paste0("(c) ", mesh_names[1], ", enlarged area"),
  paste0("(d) ", mesh_names[2], ", enlarged area")
)

panels <- list()

# (a), (b) Whole meshes, with the enlarged area as a black rectangle.
for (i in 1:2) {
  panels[[i]] <- ggplot() +
    geom_segment(data = mesh_edges[[i]],
                 aes(x, y, xend = x_end, yend = y_end),
                 linewidth = 0.06, colour = "grey45") +
    geom_point(data = cell_centres, aes(x, y), colour = colour_simulated,
               size = 0.25) +
    annotate("rect", xmin = enlarged_xlim[1], xmax = enlarged_xlim[2],
             ymin = enlarged_ylim[1], ymax = enlarged_ylim[2],
             fill = NA, colour = "black", linewidth = 0.4) +
    scale_x_continuous(breaks = overview_breaks, labels = overview_labels) +
    scale_y_continuous(breaks = overview_breaks, labels = overview_labels) +
    coord_equal(xlim = overview_xlim, ylim = overview_ylim, expand = FALSE) +
    labs(x = "Easting (km)", y = "Northing (km)", title = panel_titles[i]) +
    theme_figure
}

# (c), (d) Enlarged area.
for (i in 1:2) {
  panels[[i + 2]] <- ggplot() +
    geom_segment(data = mesh_edges[[i]],
                 aes(x, y, xend = x_end, yend = y_end),
                 linewidth = 0.15, colour = "grey45") +
    geom_point(data = cell_centres, aes(x, y), colour = colour_simulated,
               size = 0.4) +
    scale_x_continuous(breaks = enlarged_breaks, labels = enlarged_labels) +
    scale_y_continuous(breaks = enlarged_breaks, labels = enlarged_labels) +
    coord_equal(xlim = enlarged_xlim, ylim = enlarged_ylim, expand = FALSE) +
    labs(x = "Easting (km)", y = "Northing (km)",
         title = panel_titles[i + 2]) +
    theme_figure
}


# 5. Combine and save ----------------------------------------------------------

figure <- wrap_plots(panels, ncol = 2) +
  plot_layout(axis_titles = "collect")

ggsave(
  file.path(dir_figures, "figS03_meshes.png"),
  plot   = figure,
  width  = figure_width_double,
  height = 160,
  units  = "mm",
  dpi    = figure_dpi,
  bg     = "white",
  device = ragg::agg_png
)
ggsave(
  file.path(dir_figures, "figS03_meshes.pdf"),
  plot   = figure,
  width  = figure_width_double,
  height = 160,
  units  = "mm",
  dpi    = figure_dpi,
  bg     = "white",
  device = cairo_pdf
)
