# Code for "Spatial confounding can bias CPUE-based abundance indices through trending catchability covariates: a simulation study"

R code for the analyses, tables and figures of the manuscript. The scripts
are run in order from this folder. The set-level logbook data are
confidential and are not included.

| Script | What it does |
| :-- | :-- |
| `00_settings.R` | Shared settings: packages, folders, scenarios and estimators |
| `01_build_meshes.R` | Builds the meshes of the operating model and of the estimators |
| `02_fit_logbook_model.R` | Fits the model in Eq. 2 to the logbook data to obtain the operating-model parameters |
| `03_simulate_fishery_data.R` | Simulates the data sets of the fishery simulation |
| `04_fit_estimators.R` | Fits the six models to every simulated data set |
| `05_rsr.R` | Calculates restricted spatial regression (RSR) |
| `06_performance.R` | Calculates the performance measures and their summaries |
| `07_x1_density_and_catchability.R` | Analysis in which X1 affects both density and catchability (Table S3) |
| `08a_gaussian_simulation.R` | Runs the Gaussian simulation |
| `08b_gaussian_summaries.R` | Summarizes the Gaussian simulation |
| `09_tables.R` | Assembles Table 1 and Tables S1, S2 and S4 to S9 |
| `figures/` | One script per figure, named by figure number |
| `data/README.md` | Describes the input files that the scripts expect |
