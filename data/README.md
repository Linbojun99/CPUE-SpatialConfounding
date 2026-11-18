# Input data

The scripts read the files listed below from this folder. The folder is
distributed without the files.

The set-level logbook data are confidential and are not distributed with the
code (see the Data availability statement of the paper). The same applies to
files that keep the positions of individual sets: the interpolation matrix for
the sets in `output/meshes/mesh_simulation.rds`, the fitted models in
`output/logbook_model/` and the simulated data sets in `output/simulation/`.

The scripts start from the files described here.

## Files

| File | Rows | Used by |
| :-- | :-- | :-- |
| `logbook_sets.rds` | 41,035 sets used to fit Eq. 2 | `01_build_meshes.R`, `02_fit_logbook_model.R`, `figures/figS07_catch_rates.R` |
| `logbook_sets_resampling.rds` | 39,440 sets resampled in the fishery simulation | `01_build_meshes.R`, `03_simulate_fishery_data.R` |
| `simulation_grid.rds` | 1,592 cells of 1 degree of the simulation grid | `01_build_meshes.R`, `03_simulate_fishery_data.R`, `04_fit_estimators.R`, and the figure scripts `fig02`, `figS02`, `figS03` and `figS04_S06` |
| `prediction_domain_cells.csv` | 1,592 cells, with the prediction domain of 1,484 cells marked | `05_rsr.R`, `06_performance.R`, `07_x1_density_and_catchability.R`, and the figure scripts `fig02`, `figS02`, `figS04_S06`, `figS12` and `figS15` |

The `.rds` files hold a `data.table` and are read with `readRDS()`.

## `logbook_sets.rds`

One row per longline set, 2015-2023.

| Column | Meaning |
| :-- | :-- |
| `set_id` | Set identifier (integer) |
| `year` | Year of the set |
| `date` | Date of the set |
| `vessel` | Vessel code |
| `lon`, `lat` | Position of the set in decimal degrees |
| `x`, `y` | Position of the set in km, in the Lambert azimuthal equal-area projection defined in `00_settings.R` |
| `cell` | Index of the cell of the simulation grid that contains the set |
| `hooks` | Number of hooks |
| `hbf` | Hooks between floats (HBF) |
| `hbf_cl` | HBF class, a factor with the levels `16-17` (reference), `le15` (15 or fewer) and `ge18` (18 or more), in this order |
| `bet_n` | Catch of bigeye tuna in numbers |

`01_build_meshes.R` uses `x` and `y`. `02_fit_logbook_model.R` uses `year`,
`vessel`, `x`, `y`, `hooks`, `hbf_cl` and `bet_n`.
`figures/figS07_catch_rates.R` uses `hooks` and `bet_n`.

## `logbook_sets_resampling.rds`

The sets retained for resampling, a subset of the rows of `logbook_sets.rds`
with the same columns. `01_build_meshes.R` uses `x` and `y`.
`03_simulate_fishery_data.R` uses `set_id`, `year`, `date` (for the month),
`vessel`, `lon`, `lat`, `x`, `y`, `cell`, `hooks`, `hbf` and `hbf_cl`. The
observed catch is not used. The order of the rows matters, because sets are
drawn by row number.

## `simulation_grid.rds`

One row per cell of the simulation grid. The rows are ordered by `cell`.

| Column | Meaning |
| :-- | :-- |
| `cell` | Cell index, 1 to 1,592 |
| `cell_id` | Cell label built from the latitude and longitude of the centre |
| `cell_lon`, `cell_lat` | Longitude and latitude of the cell centre in decimal degrees |
| `x`, `y` | Cell centre in km, in the projection defined in `00_settings.R` |
| `area_km2` | Cell area in square kilometres |

`04_fit_estimators.R` uses `cell`, `x` and `y`; the figure scripts use `cell`,
`cell_lon`, `cell_lat`, `x` and `y`. The meshes are built from `x` and `y`.
These coordinates must be stored at full precision (hence the `.rds`
format): coordinates rounded to 15 significant digits give a different mesh.

## `prediction_domain_cells.csv`

One row per cell of the simulation grid.

| Column | Meaning |
| :-- | :-- |
| `cell` | Cell index, as in `simulation_grid.rds` |
| `in_domain` | `TRUE` for the 1,484 cells of the prediction domain, over which the true and estimated indices are calculated (Eq. 9) |
| `has_sets` | `TRUE` if at least one of the sets retained for resampling lies in the cell |

The scripts use `cell` and `in_domain`.
