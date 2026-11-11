source("scripts/simulering/simulering.r")

dir.create("results", showWarnings = FALSE)

# Kap. 3, NOPW-settingen: 4 K x 3 nivaaer av strata-stoy, 5000 replikasjoner hver.
# rho = 1/sqrt(1 + 12*sigma^2): sigma 0 / 0.500 / 1.118 -> rho 1.00 / 0.50 / 0.25.
# Forrige runde (backup i results/archive_2026-09-21) brukte rho 1.00 / 0.69 / 0.28.
scenarios = expand.grid(K = c(1, 2, 4.2, 10), strata_noise = c(0, 0.500, 1.118))

scenarios |>
  run_grid(cohort_size = 24411, cens_cut = 0.1865, Tau1 = 0, Tau2 = 0.15,
           iterations = 5000, seed = 4711, cores = 4) |>
  save_runs("results/nopw_grid_runs.rds") |>
  grid_table(scenarios) |>
  write_results("results/nopw_grid_results.csv")
