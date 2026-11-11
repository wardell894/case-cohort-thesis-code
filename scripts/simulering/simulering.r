source("scripts/simulering/helpers.R")

# Simuleringen i kap. 3 (en kovariat) og 3.3 (to kovariater) som en pipe.
# Antall kovariater styres av beta:
#   kap. 3:  beta = c(Z = 1)
#   3.3:     beta = c(Z1 = 1, Z2 = -0.25), Z2 = -Z1 + N(0, sd_z2^2)
# Strataene dannes paa forste kovariat; en felles proxy P imputerer alle.

run_casecohort_sim = function(iterations   = 10,
                              cohort_size  = 20000,
                              beta         = c(Z = 1),
                              K            = 1,
                              Tau1         = 0,
                              Tau2         = 0.3,
                              noise_proxy  = 0.1,
                              cens_cut     = 0.5,
                              seed         = NULL,
                              strata       = 2,
                              strata_noise = 0,
                              sd_z2        = 0.5,
                              pilot_reps   = 200,
                              m_counts     = NULL,
                              verbose      = TRUE) {

  if (!is.null(seed)) set.seed(seed)

  covs    = names(beta)
  proxies = setNames(as.list(rep("P", length(covs))), covs)

  # m_counts kan gis utenfra, slik at parallelle kjoringer deler samme design.
  if (is.null(m_counts)) {
    E_cases  = expected_cases_pilot(pilot_reps, strata, covs[1], strata_noise,
                                    cohort_size = cohort_size, beta = beta,
                                    cens_cut = cens_cut, sd_z2 = sd_z2)
    m_counts = eg_subcohort_counts(E_cases, cohort_size, K)
  }

  if (verbose) message("m_counts fra EG-regelen (K = ", K, "): ", paste(m_counts, collapse = ", "))

  estimates = vector("list", iterations)
  ratios    = vector("list", iterations)

  for (i in seq_len(iterations)) {
    dat = make_cohort(cohort_size, beta, cens_cut, noise_proxy, sd_z2) |>
      assign_strata(covs[1], strata, strata_noise) |>
      sample_subcohort(m_counts) |>
      mask_phase2(covs)

    estimates[[i]] = fit_all_methods(dat, covs, proxies, Tau1, Tau2)
    ratios[[i]]    = stratum_ratios(dat)

    if (verbose && i %% 100 == 0)
      cat("K", K, "noise", strata_noise, i, "/", iterations, "\n")
  }

  summarise_results(estimates, beta, do.call(rbind, ratios),
                    fracs = m_counts / (cohort_size / strata))
}

#caffeinate -i Rscript run_simulation.r for aa kjore over natta
