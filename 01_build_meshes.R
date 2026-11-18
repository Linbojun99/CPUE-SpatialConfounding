# ==============================================================================
# 01_build_meshes.R
#
# Builds the triangular meshes on which the spatial and spatiotemporal fields
# are represented (Supplementary Table S1): the mesh used to simulate the
# fields of the operating model, the mesh used by the estimators, the coarser
# and finer meshes of the mesh sensitivity comparison, and the mesh used to
# fit Eq. 2 to the logbook data. Also builds the matrices that interpolate
# values at the vertices of the simulation mesh to the centres of the grid
# cells and to the positions of the sets.
#
# Input:
#   data/simulation_grid.rds            1,592 cells of the simulation grid
#   data/logbook_sets.rds               41,035 sets used to fit Eq. 2
#   data/logbook_sets_resampling.rds    39,440 sets retained for resampling
# Output:
#   output/meshes/mesh_simulation.rds           mesh of the operating model
#                                               and interpolation matrices
#   output/meshes/mesh_estimation.rds           mesh of the estimators
#   output/meshes/mesh_estimation_coarser.rds   coarser mesh (sensitivity)
#   output/meshes/mesh_estimation_finer.rds     finer mesh (sensitivity)
#   output/meshes/mesh_logbook.rds              mesh used to fit Eq. 2
#   output/meshes/mesh_settings.csv             settings and numbers of
#                                               vertices of the meshes
#                                               (Supplementary Table S1)
# ==============================================================================

source("00_settings.R")


# 1. Mesh settings (Supplementary Table S1) ------------------------------------

# Distances are in km in the projection of 00_settings.R. The first and second
# values of a pair refer to the interior and to the outer extension. The
# cutoff is the smallest distance allowed between two input locations.

# Mesh of the operating model. It is finer than the estimation meshes so that
# covariate components with the shortest range (300 km) can be represented.
simulation_cutoff    <- 40
simulation_max_edge  <- c(100, 500)
simulation_extension <- c(300, 2000)

# Mesh of the estimators (Section 2.4).
estimation_cutoff    <- 150
estimation_max_edge  <- c(200, 1500)
estimation_extension <- c(300, 1500)

# Coarser and finer meshes of the mesh sensitivity comparison. They differ
# from the estimation mesh in the cutoff and the interior maximum edge only.
coarser_cutoff   <- 250
coarser_max_edge <- c(350, 1500)
finer_cutoff     <- 100
finer_max_edge   <- c(150, 1500)

# Mesh used to fit Eq. 2 to the logbook data.
logbook_cutoff <- 150


# 2. Locations -----------------------------------------------------------------

# Projected centres (x, y) of the cells of the simulation grid. The
# coordinates are used at full precision. Rounding them changes the meshes.
simulation_grid <- readRDS(file.path(dir_data, "simulation_grid.rds"))
grid_xy         <- as.matrix(simulation_grid[, .(x, y)])

# Projected positions of the logbook sets and of the sets retained for
# resampling.
logbook_sets    <- readRDS(file.path(dir_data, "logbook_sets.rds"))
logbook_xy      <- as.matrix(logbook_sets[, .(x, y)])
resampling_sets <- readRDS(file.path(dir_data, "logbook_sets_resampling.rds"))
resampling_xy   <- as.matrix(resampling_sets[, .(x, y)])


# 3. Mesh of the operating model -----------------------------------------------

simulation_mesh <- fm_mesh_2d_inla(
  loc      = grid_xy,
  max.edge = simulation_max_edge,
  cutoff   = simulation_cutoff,
  offset   = simulation_extension
)
cat("Vertices of the simulation mesh:", simulation_mesh$n, "\n")

# Interpolation matrices (piecewise linear basis functions): one row per grid
# cell or set, one column per mesh vertex.
A_grid <- fm_basis(simulation_mesh, loc = grid_xy)
A_sets <- fm_basis(simulation_mesh, loc = resampling_xy)

saveRDS(
  list(mesh = simulation_mesh, A_grid = A_grid, A_sets = A_sets),
  file.path(dir_meshes, "mesh_simulation.rds")
)


# 4. Mesh of the estimators ----------------------------------------------------

estimation_mesh <- fm_mesh_2d_inla(
  loc      = grid_xy,
  max.edge = estimation_max_edge,
  cutoff   = estimation_cutoff,
  offset   = estimation_extension
)
cat("Vertices of the estimation mesh:", estimation_mesh$n, "\n")

saveRDS(estimation_mesh, file.path(dir_meshes, "mesh_estimation.rds"))


# 5. Coarser and finer meshes of the sensitivity comparison --------------------

coarser_mesh <- fm_mesh_2d_inla(
  loc      = grid_xy,
  max.edge = coarser_max_edge,
  cutoff   = coarser_cutoff,
  offset   = estimation_extension
)
finer_mesh <- fm_mesh_2d_inla(
  loc      = grid_xy,
  max.edge = finer_max_edge,
  cutoff   = finer_cutoff,
  offset   = estimation_extension
)
cat("Vertices of the coarser and finer meshes:",
    coarser_mesh$n, finer_mesh$n, "\n")

saveRDS(coarser_mesh, file.path(dir_meshes, "mesh_estimation_coarser.rds"))
saveRDS(finer_mesh, file.path(dir_meshes, "mesh_estimation_finer.rds"))


# 6. Mesh used to fit Eq. 2 to the logbook data --------------------------------

# This mesh is built from the positions of the logbook sets, without added
# refinement or extension. It is the mesh that sdmTMB::make_mesh() builds from
# the same positions with cutoff = 150.
logbook_mesh <- fm_rcdt_2d_inla(
  loc    = logbook_xy,
  refine = list(),
  extend = list(),
  cutoff = logbook_cutoff
)
cat("Vertices of the mesh used to fit Eq. 2:", logbook_mesh$n, "\n")

saveRDS(logbook_mesh, file.path(dir_meshes, "mesh_logbook.rds"))


# 7. Supplementary Table S1 ----------------------------------------------------

mesh_settings <- data.table(
  use             = c("Operating model",
                      "Coarser estimation mesh",
                      "Selected estimation mesh",
                      "Finer estimation mesh",
                      "Model fitted to the logbook data (Eq. 2)"),
  vertices        = c(simulation_mesh$n, coarser_mesh$n, estimation_mesh$n,
                      finer_mesh$n, logbook_mesh$n),
  cutoff          = c(simulation_cutoff, coarser_cutoff, estimation_cutoff,
                      finer_cutoff, logbook_cutoff),
  max_edge_inner  = c(simulation_max_edge[1], coarser_max_edge[1],
                      estimation_max_edge[1], finer_max_edge[1], NA),
  max_edge_outer  = c(simulation_max_edge[2], coarser_max_edge[2],
                      estimation_max_edge[2], finer_max_edge[2], NA),
  extension_inner = c(simulation_extension[1], estimation_extension[1],
                      estimation_extension[1], estimation_extension[1], NA),
  extension_outer = c(simulation_extension[2], estimation_extension[2],
                      estimation_extension[2], estimation_extension[2], NA)
)
print(mesh_settings)

fwrite(mesh_settings, file.path(dir_meshes, "mesh_settings.csv"))
