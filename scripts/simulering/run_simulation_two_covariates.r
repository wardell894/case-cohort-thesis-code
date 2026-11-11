source("scripts/simulering/simulering.r")

dir.create("results", showWarnings = FALSE)

# Kap. 3.3, Samuelsen-settingen: to korrelerte kovariater, K = 4.2, strata paa Z1.
# Z2 = -Z1 + N(0, 0.5^2), cor(Z1, Z2) = -0.5; beta2 = -0.25; felles proxy P for begge.
# 5000 replikasjoner paa aatte kjerner, alle med samme design (m_counts).
beta    = c(Z1 = 1, Z2 = -0.25)
setting = list(cohort_size = 1000, cens_cut = 0.5, Tau1 = 0, Tau2 = 0.3,
               K = 4.2, sd_z2 = 0.5)

design_counts(setting, beta, seed = 4711) |>
  run_split(setting, beta, iterations = 5000, cores = 8, seed = 4711) |>
  save_runs("results/grid_two_covariates_runs.rds") |>
  pool_runs(beta) |>
  write_results("results/grid_two_covariates_results.csv")
