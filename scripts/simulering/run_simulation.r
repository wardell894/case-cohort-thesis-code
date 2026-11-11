source("scripts/simulering/simulering.r")

dir.create("results", showWarnings = FALSE)

# Kap. 3, Samuelsen-settingen: 4 K x 3 nivaaer av strata-stoy, 5000 replikasjoner hver.
# rho = 1/sqrt(1 + 12*sigma^2), siden Var(Z) = 1/12 for Z ~ U[0,1].
#   sigma 0     -> rho 1.00
#   sigma 0.500 -> rho 0.50
#   sigma 1.118 -> rho 0.25
# Forrige runde brukte c(0, 0.3, 1.0), altsaa rho 1.00 / 0.69 / 0.28.
scenarios = expand.grid(K = c(1, 2, 4.2, 6), strata_noise = c(0, 0.500, 1.118))

scenarios |>
  run_grid(cohort_size = 1000, cens_cut = 0.5, Tau1 = 0, Tau2 = 0.3,
           iterations = 5000, seed = 4711, cores = 4) |>
  save_runs("results/grid_runs.rds") |>
  grid_table(scenarios) |>
  write_results("results/grid_results.csv")

writeLines(capture.output(sessionInfo()), "results/sessionInfo.txt")
