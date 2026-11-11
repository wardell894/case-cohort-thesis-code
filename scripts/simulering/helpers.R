library(CaseCohortCoxSurvival)
library(survival)
library(dplyr)
library(parallel)

# Byggesteiner for simuleringene i kap. 3 (en kovariat) og 3.3 (to kovariater).
# Trekkerekkefolgen er den samme som i den gamle simulering.r, slik at
# tidligere kjoringer reproduserer seg med samme seed (se tests/test_helpers.r).

# --- Kohort --------------------------------------------------------------------

# Weibull-tider (form 2) og uniform sensur, gitt lineaer prediktor lp.
# Trekker U og saa sensur, i samme rekkefolge som for.
simulate_times = function(lp, cens_cut) {
  n  = length(lp)
  U  = runif(n)
  Tt = sqrt(-log(U) / exp(lp))
  Ct = runif(n, 0, cens_cut)
  data.frame(status = as.integer(Tt <= Ct), event_time = pmin(Tt, Ct))
}

# Kohort med kovariater etter navnene i beta. En kovariat: Z ~ U(0,1).
# To: Z2 = -Z1 + N(0, sd_z2^2), som gir cor(Z1, Z2) = -0.5 med sd_z2 = 0.5.
# Proxyen P = Z1 + N(0, noise_proxy^2) dekker alle kovariatene.
# proxy = FALSE hopper over rnorm for P, som i piloten.
make_cohort = function(cohort_size, beta, cens_cut, noise_proxy = 0.1,
                       sd_z2 = 0.5, proxy = TRUE) {
  covs = names(beta)
  Z    = data.frame(runif(cohort_size))
  if (length(beta) == 2) Z[[2]] = -Z[[1]] + rnorm(cohort_size, sd = sd_z2)
  names(Z) = covs

  P  = if (proxy) Z[[1]] + rnorm(cohort_size, sd = noise_proxy)
  lp = drop(as.matrix(Z) %*% beta)

  dat = data.frame(id = seq_len(cohort_size), simulate_times(lp, cens_cut), Z)
  dat[paste0(covs, "_full")] = Z
  if (proxy) dat$P = P
  dat
}

# --- Design ----------------------------------------------------------------------

# Strata som like store grupper paa stratifiseringsvariabelen by, med paafort stoy
# for rho-aksen i kap. 3. noise = 0 hopper over rnorm, slik at RNG-sekvensen er uendret.
assign_strata = function(dat, by, strata, noise = 0) {
  x = dat[[by]]
  dat$Z_strata = if (noise > 0) x + rnorm(nrow(dat), sd = noise) else x
  dat$W = ntile(dat$Z_strata, strata) - 1
  dat
}

# Trekker m_counts[w + 1] individer i hvert stratum w, uavhengig av utfall.
sample_subcohort = function(dat, m_counts) {
  dat$subcohort = 0L
  for (w in 0:(length(m_counts) - 1)) {
    idx = which(dat$W == w)
    dat$subcohort[idx[sample(length(idx), min(m_counts[w + 1], length(idx)))]] = 1L
  }
  dat
}

# Fase-to-kovariatene er bare kjent i subkohorten og for casene.
mask_phase2 = function(dat, covs) {
  in.cc = dat$subcohort == 1 | dat$status == 1
  dat[!in.cc, covs] = NA
  dat
}

# --- EG-allokering ---------------------------------------------------------------

# Forventet antall caser per stratum, fra pilot_reps kohorter uten proxy.
# ... sendes til make_cohort, saa piloten bruker samme generator som simuleringen.
expected_cases_pilot = function(pilot_reps, strata, strata_by, strata_noise = 0, ...) {
  counts = matrix(0, nrow = pilot_reps, ncol = strata)

  for (i in seq_len(pilot_reps)) {
    dat = make_cohort(..., proxy = FALSE) |>
      assign_strata(strata_by, strata, strata_noise)

    for (w in 0:(strata - 1)) {
      counts[i, w + 1] = sum(dat$W == w & dat$status == 1)
    }
  }

  colMeans(counts)
}

# EGs regel: M_l = K * odds_l * N_l, med odds fra forventede caser per stratum.
eg_subcohort_counts = function(E_cases, cohort_size, K) {
  n_j  = cohort_size / length(E_cases)
  p_j  = E_cases / n_j
  odds = p_j / (1 - p_j)
  floor(odds * n_j * K + 0.5)
}

# --- Tilpasning ------------------------------------------------------------------

# Diagonalen av en varians. CaseCohortCoxSurvival gir den allerede som vektor;
# coxph og cch gir matriser.
dg = function(v) unname(if (is.null(dim(v))) v else diag(v))

# Kohortleddet i EG per koeffisient, forste sum i likn. (2.66):
# (n/(n-1)) * sum_i D_ik^2 / w_i, diagonalen av krysseproduktet.
eg_phase1_var = function(coxph, n) {
  D = as.matrix(residuals(coxph, type = "dfbeta", weighted = TRUE))
  unname((n / (n - 1)) * colSums(D^2 / coxph$weights))
}

# Full kohort, EG med design- og kalibrerte vekter, og Borgan II (cch).
# covs: fase-to-kovariatene. proxies: navngitt liste, kovariat -> proxykolonne.
fit_all_methods = function(dat, covs, proxies, Tau1, Tau2) {
  fit_full_cohort = caseCohortCoxSurvival(dat, status = "status",
                                          time = "event_time",
                                          cox.phase1 = paste0(covs, "_full"),
                                          Tau1 = Tau1, Tau2 = Tau2, print = 0)

  fit_uncalibrated = caseCohortCoxSurvival(dat, status = "status",
                                           time = "event_time",
                                           cox.phase2 = covs,
                                           strata = "W",
                                           subcohort = "subcohort",
                                           Tau1 = Tau1, Tau2 = Tau2, print = 0)

  coxph = fit_uncalibrated$coxph.fit

  fit_calibrated = caseCohortCoxSurvival(dat, status = "status",
                                         time = "event_time",
                                         cox.phase2 = covs,
                                         other.covars = unique(unlist(proxies)),
                                         strata = "W", subcohort = "subcohort",
                                         calibrated = TRUE, predict = TRUE,
                                         predictors.cox.phase2 = proxies,
                                         aux.method = "Breslow",
                                         Tau1 = Tau1, Tau2 = Tau2, print = 0)

  cases_subcohort = dat[dat$subcohort == 1 | dat$status == 1, ]

  cch_borgan2 = cch(reformulate(covs, "Surv(event_time, status)"), data = cases_subcohort,
                    subcoh      = ~subcohort,
                    id          = ~id,
                    stratum     = factor(cases_subcohort$W),
                    cohort.size = c(table(dat$W)),
                    method      = "II.Borgan")

  # Dekomponering i kohortledd og samplingledd. Restene regnes i etterkant:
  #   EG fase 2        = var_uncal  - var_uncal_ph1
  #   Borgan II fase 1 = var_borgan - var_borgan_ph2
  # var_uncal_naive er Itilde^-1 fra SAMME fit som var_uncal_ph1.
  first = dat[[paste0(covs[1], "_full")]]

  out = list(
    coef                = covs,
    beta_uncal          = unname(fit_uncalibrated$beta[covs]),
    var_uncal           = dg(fit_uncalibrated$beta.var[covs]),
    beta_cal            = unname(fit_calibrated$beta[covs]),
    var_cal             = dg(fit_calibrated$beta.var[covs]),
    beta_full_cohort    = unname(fit_full_cohort$beta),
    var_full_cohort     = dg(fit_full_cohort$coxph.fit$naive.var),
    var_robust_coxph    = dg(coxph$var),
    var_robust_cccs     = dg(fit_uncalibrated$beta.robustvar[covs]),
    beta_borgan         = unname(cch_borgan2$coefficients),
    var_borgan          = dg(cch_borgan2$var),
    var_uncal_ph1       = eg_phase1_var(coxph, nrow(dat)),
    var_borgan_ph2      = dg(cch_borgan2$phase2var),
    var_uncal_naive     = dg(coxph$naive.var),
    number_of_cases     = sum(dat$status == 1),
    non_cases_subcohort = sum(dat$status == 0 & dat$subcohort == 1),
    correlation_proxy   = cor(first, dat[[proxies[[1]]]]),
    correlation_strata  = cor(first, dat$Z_strata)
  )
  if (length(covs) == 2) out$cor_z1_z2 = cor(first, dat[[paste0(covs[2], "_full")]])
  c(out, as.list(strata_r2(dat, covs)))
}

# --- Diagnostikk og oppsummering -------------------------------------------------

# Ikke-caser i subkohorten per case, per stratum.
stratum_ratios = function(dat) {
  strata = sort(unique(dat$W))
  sapply(strata, function(w) {
    i = dat$W == w
    sum(i & dat$subcohort == 1 & dat$status == 0) / sum(i & dat$status == 1)
  })
}

# R^2 fra regresjon paa stratum, som tab:strata-r2 i kap. 4. Med to kovariater
# ogsaa partiell R^2: residualen av den ene gitt den andre.
strata_r2 = function(dat, covs) {
  W    = factor(dat$W)
  full = dat[paste0(covs, "_full")]
  r2   = function(x) summary(lm(x ~ W))$r.squared

  out = setNames(sapply(full, r2), paste0("r2_", covs))
  if (length(covs) == 2) {
    out[paste0("r2_", covs[1], "_given_", covs[2])] = r2(resid(lm(full[[1]] ~ full[[2]])))
    out[paste0("r2_", covs[2], "_given_", covs[1])] = r2(resid(lm(full[[2]] ~ full[[1]])))
  }
  out
}

coverage_prob = function(beta, variance, true_beta) {
  beta     = as.numeric(beta)
  variance = as.numeric(variance)
  mean(abs(beta - true_beta) <= qnorm(0.975) * sqrt(variance))
}

# Tabell per koeffisient: gjennomsnittlig beta, empirisk og estimert varians, dekning.
# true_beta er samme navngitte vektor som beta i make_cohort.
summarise_results = function(estimates, true_beta, ratios = NULL, fracs = NULL) {
  df = bind_rows(lapply(estimates, as.data.frame))

  out = bind_rows(lapply(names(true_beta), function(k) {
    d  = df[df$coef == k, ]
    tb = true_beta[[k]]
    data.frame(
      coef      = k,
      method    = c("Robust", "Borgan.II", "uncalibrated", "calibrated", "full_cohort"),
      mean_beta = c(mean(d$beta_uncal), mean(d$beta_borgan), mean(d$beta_uncal),
                    mean(d$beta_cal), mean(d$beta_full_cohort)),
      emp_var   = c(var(d$beta_uncal), var(d$beta_borgan), var(d$beta_uncal),
                    var(d$beta_cal), var(d$beta_full_cohort)),
      est_var   = c(mean(d$var_robust_coxph), mean(d$var_borgan), mean(d$var_uncal),
                    mean(d$var_cal), mean(d$var_full_cohort)),
      coverage  = c(coverage_prob(d$beta_uncal, d$var_robust_coxph, tb),
                    coverage_prob(d$beta_borgan, d$var_borgan, tb),
                    coverage_prob(d$beta_uncal, d$var_uncal, tb),
                    coverage_prob(d$beta_cal, d$var_cal, tb),
                    coverage_prob(d$beta_full_cohort, d$var_full_cohort, tb))
    )
  }))
  out[-(1:2)] = round(out[-(1:2)], 5)

  attr(out, "raw") = df

  # Diagnostikken er per datasett, saa den leses fra radene for forste koeffisient.
  d1   = df[df$coef == names(true_beta)[1], ]
  r2   = grep("^r2_|^cor_z1_z2$", names(d1), value = TRUE)
  diag = c(avg_cases              = mean(d1$number_of_cases),
           avg_correlation_proxy  = mean(d1$correlation_proxy),
           avg_correlation_strata = mean(d1$correlation_strata),
           avg_non_cases_cohort   = mean(d1$non_cases_subcohort),
           avg_ratio_total        = mean(d1$non_cases_subcohort / d1$number_of_cases),
           setNames(colMeans(d1[r2]), paste0("avg_", r2)))
  attr(out, "diagnostics") = round(diag, 3)

  if (!is.null(ratios)) {
    r = round(colMeans(ratios), 2)
    names(r) = paste0("ratio_stratum_", seq_along(r) - 1)
    attr(out, "diagnostics") = c(attr(out, "diagnostics"), r)
  }

  if (!is.null(fracs)) {
    f = round(fracs, 3)
    names(f) = paste0("frac_stratum_", seq_along(f) - 1)
    attr(out, "diagnostics") = c(attr(out, "diagnostics"),
                                 f, frac_total = round(mean(fracs), 3))
  }

  out
}

# --- Kjoring ---------------------------------------------------------------------

# Kap. 3: en rad i scenarios (K, strata_noise) per kjerne, alle med samme seed.
run_grid = function(scenarios, ..., iterations, seed, cores) {
  mclapply(seq_len(nrow(scenarios)), function(i) {
    run_casecohort_sim(iterations   = iterations,
                       K            = scenarios$K[i],
                       strata_noise = scenarios$strata_noise[i],
                       seed         = seed,
                       verbose      = TRUE,
                       ...)
  }, mc.cores = cores)
}

# EGs subkohortstorrelser for ett scenario, regnet en gang foer kjoringen.
design_counts = function(setting, beta, seed, pilot_reps = 200, strata = 2, strata_noise = 0) {
  set.seed(seed)
  E_cases  = expected_cases_pilot(pilot_reps, strata, names(beta)[1], strata_noise,
                                  cohort_size = setting$cohort_size, beta = beta,
                                  cens_cut = setting$cens_cut, sd_z2 = setting$sd_z2)
  m_counts = eg_subcohort_counts(E_cases, setting$cohort_size, setting$K)
  message("m_counts fra EG-regelen (K = ", setting$K, "): ", paste(m_counts, collapse = ", "))
  m_counts
}

# 3.3: ett scenario delt paa flere kjerner med hver sin seed og samme design.
run_split = function(m_counts, setting, beta, iterations, cores, seed) {
  seeds = seed + seq_len(cores) - 1
  mclapply(seeds, function(s) {
    do.call(run_casecohort_sim, c(setting, list(iterations = iterations / cores, beta = beta,
                                                seed = s, m_counts = m_counts, verbose = TRUE)))
  }, mc.cores = cores)
}

# Lagrer rds og sender runs videre i pipen.
save_runs = function(runs, file) {
  saveRDS(runs, file)
  runs
}

# Kap. 3: en rad per scenario og metode, med diagnostikken som kolonner.
grid_table = function(runs, scenarios) {
  do.call(rbind, lapply(seq_along(runs), function(i) {
    data.frame(K            = scenarios$K[i],
               strata_noise = scenarios$strata_noise[i],
               runs[[i]],
               as.data.frame(as.list(attr(runs[[i]], "diagnostics"))),
               row.names = NULL)
  }))
}

# 3.3: slaar sammen de raa tallene fra kjernene og oppsummerer paa nytt.
pool_runs = function(runs, beta) {
  raw = do.call(rbind, lapply(runs, attr, "raw"))
  estimates = split(raw, rep(seq_len(nrow(raw) / length(beta)), each = length(beta)))
  summarise_results(estimates, true_beta = beta)
}

# Skriver CSV, med diagnostikken som kolonner hvis den ikke allerede er med.
write_results = function(results, file) {
  d   = attr(results, "diagnostics")
  out = as.data.frame(results)
  if (!is.null(d)) {
    d   = d[setdiff(names(d), names(out))]
    out = data.frame(out, as.list(d), row.names = NULL)
  }
  write.csv(out, file, row.names = FALSE)
  results
}
