# Figurer til kapittel 3. Kjøres fra analysis/:  Rscript scripts/figurer.r
# Leser bare ferdige resultater (results/*.csv, results/*.rds) - ingen simulering.
# Figurene lagres som PDF i figures/ i tekstbredde (14,7 cm).

library(ggplot2)
library(patchwork)
library(scales)
library(ggdag)

out_dir <- "figures"
dir.create(out_dir, showWarnings = FALSE)

text_width <- 5.79   # tommer, = \textwidth i oppgaven (418 pt)

# Farger (validert): rho er ordnet -> én blåtone lys -> mørk; estimatorene kategoriske.
rho_cols <- c("0.25" = "#86b6ef", "0.50" = "#2a78d6", "1.00" = "#104281")
est_cols <- c("Robust" = "#2a78d6", "Borgan II" = "#eb6834", "EG" = "#1baf7a")

theme_thesis <- function() {
  theme_minimal(base_size = 9) +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(colour = "grey90", linewidth = 0.3),
      axis.line = element_line(colour = "grey40", linewidth = 0.3),
      axis.ticks = element_line(colour = "grey40", linewidth = 0.3),
      strip.text = element_text(face = "bold", hjust = 0),
      legend.position = "top",
      legend.justification = "left",
      legend.margin = margin(0, 0, 0, 0),
      legend.box.spacing = unit(2, "pt"),
      plot.margin = margin(4, 6, 4, 4)
    )
}

save_fig <- function(p, name, height, width = text_width) {
  ggsave(file.path(out_dir, paste0(name, ".pdf")), p,
         width = width, height = height, device = pdf)
}

# --- Figur (3.1): hva rho betyr ------------------------------------------------
# Strata: Z* = Z + eps, eps ~ N(0, sigma^2), delt ved medianen 0.5.
# P(høyrisikostratum | Z = z) = Phi((z - 0.5) / sigma); trappetrinn ved sigma = 0.
fig_strata <- function() {
  z <- seq(0, 1, length.out = 401)
  sig <- c("1.00" = 0, "0.50" = 0.5, "0.25" = 1.118)
  d <- do.call(rbind, lapply(names(sig), function(r) {
    s <- sig[[r]]
    p <- if (s == 0) as.numeric(z > 0.5) else pnorm((z - 0.5) / s)
    data.frame(z = z, p = p, rho = r)
  }))
  d$rho <- factor(d$rho, levels = c("1.00", "0.50", "0.25"))
  ends <- d[d$z == 1, ]

  ggplot(d, aes(z, p, colour = rho)) +
    geom_line(linewidth = 0.8) +
    geom_text(data = ends, aes(label = paste0("rho == '", rho, "'")), parse = TRUE,
              hjust = 0, nudge_x = 0.02, size = 3, show.legend = FALSE) +
    scale_colour_manual(values = rho_cols, guide = "none") +
    scale_x_continuous(breaks = seq(0, 1, 0.25), expand = expansion(mult = c(0.01, 0.14))) +
    scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25)) +
    labs(x = expression(paste("Covariate ", italic(Z))),
         y = "P(high-risk stratum | Z)") +
    theme_thesis()
}

save_fig(fig_strata(), "sim_strata", height = 2.4)

# --- Felles: resultater per design fra CSV-ene ----------------------------------
sigma_to_rho <- function(s) factor(c("0" = "1.00", "0.5" = "0.50", "1.118" = "0.25")[as.character(s)],
                                   levels = c("1.00", "0.50", "0.25"))
read_grid <- function(f, setting) {
  d <- read.csv(f)
  d$rho <- sigma_to_rho(d$strata_noise)
  d$setting <- setting
  d
}
grid <- rbind(read_grid("results/grid_results.csv", "Samuelsen"),
              read_grid("results/nopw_grid_results.csv", "NOPW"))
grid$setting <- factor(grid$setting, levels = c("Samuelsen", "NOPW"))

# K-aksen: bare designene som finnes i panelet (6 i Samuelsen, 10 i NOPW).
k_breaks <- function(lim) if (lim[2] > 8) c(1, 2, 4.2, 10) else c(1, 2, 4.2, 6)

# --- Figur (3.2.2): effektivitet (Reilly s. 219) --------------------------------
# Eff = empirisk varians i fullkohorten / empirisk varians i case-kohorten, per design
# (samme 5000 kohorter i teller og nevner). EG og robust deler estimatet (uncalibrated).
re_data <- function() {
  key <- c("setting", "K", "rho")
  fc  <- grid[grid$method == "full_cohort",  c(key, "emp_var")]
  des <- grid[grid$method == "uncalibrated", c(key, "emp_var")]
  cal <- grid[grid$method == "calibrated",   c(key, "emp_var")]
  m <- merge(merge(fc, des, by = key, suffixes = c("_full", "_des")), cal, by = key)
  rbind(data.frame(m[key], weights = "Design weights",     re = m$emp_var_full / m$emp_var_des),
        data.frame(m[key], weights = "Calibrated weights", re = m$emp_var_full / m$emp_var))
}

fig_re <- function(with_calibrated) {
  d <- re_data()
  d$weights <- factor(d$weights, levels = c("Design weights", "Calibrated weights"))
  if (!with_calibrated) d <- d[d$weights == "Design weights", ]
  ends <- do.call(rbind, lapply(split(d[d$weights == "Design weights", ],
                                      ~ setting + rho, drop = TRUE),
                                function(x) x[which.max(x$K), ]))
  p <- ggplot(d, aes(K, re, colour = rho, linetype = weights,
                     group = interaction(rho, weights))) +
    # Tommelfingerregelen K/(K+1) som referanse (teksten i 3.2.2 viser til den)
    geom_function(fun = function(k) k / (k + 1), inherit.aes = FALSE,
                  colour = "grey55", linetype = "dotted", linewidth = 0.5) +
    geom_line(linewidth = 0.7) +
    geom_point(size = 1.8) +
    facet_wrap(~ setting, scales = "free_x") +
    scale_colour_manual(values = rho_cols, name = expression(rho)) +
    scale_linetype_manual(values = c("Design weights" = "solid", "Calibrated weights" = "22"),
                          name = NULL) +
    scale_x_continuous(breaks = k_breaks, labels = function(x) format(x, drop0trailing = TRUE)) +
    scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25),
                       labels = label_percent()) +
    labs(x = expression(paste("Non-cases sampled per case, ", italic(K))),
         y = "Efficiency") +
    theme_thesis() +
    theme(panel.spacing = unit(1.2, "lines"))
  if (!with_calibrated) p <- p + guides(linetype = "none")
  p
}

save_fig(fig_re(FALSE), "sim_re",        height = 2.8)   # 3.2.2: bare designvekter
save_fig(fig_re(TRUE),  "sim_re_both",   height = 2.8)   # 3.2.6: design mot kalibrert

# --- Figur (3.2.3): robust mot Borgan II - treffer de den sanne variansen? ---------
# y = snitt av variansestimatet / empirisk varians av estimatorens eget beta-hat.
# Paneler: setting (rader) x rho (kolonner); facet_wrap for fri K-akse per rad.
fig_robust_bii <- function() {
  d <- grid[grid$method %in% c("Robust", "Borgan.II", "uncalibrated"), ]
  d$estimator <- factor(c(Robust = "Robust", Borgan.II = "Borgan II", uncalibrated = "EG")[d$method],
                        levels = c("Robust", "Borgan II", "EG"))
  d$ratio <- d$est_var / d$emp_var
  d$panel <- factor(paste0(d$setting, "*','~~rho==", "'", d$rho, "'"),
                    levels = as.vector(outer(c("Samuelsen", "NOPW"), c("1.00", "0.50", "0.25"),
                                             function(s, r) paste0(s, "*','~~rho==", "'", r, "'")))[c(1, 3, 5, 2, 4, 6)])
  ggplot(d, aes(K, ratio, colour = estimator, linetype = estimator)) +
    geom_hline(yintercept = 1, colour = "grey40", linewidth = 0.4) +
    geom_line(linewidth = 0.7) +
    geom_point(size = 1.8) +
    facet_wrap(~ panel, nrow = 2, scales = "free_x", labeller = label_parsed) +
    scale_colour_manual(values = est_cols, name = NULL) +
    scale_linetype_manual(values = c("Robust" = "solid", "Borgan II" = "22", "EG" = "11"), name = NULL) +
    scale_x_continuous(breaks = k_breaks, labels = function(x) format(x, drop0trailing = TRUE)) +
    scale_y_continuous(breaks = seq(0.9, 1.4, 0.1)) +
    labs(x = expression(paste("Non-cases sampled per case, ", italic(K))),
         y = "Mean variance estimate / empirical variance") +
    theme_thesis() +
    theme(panel.spacing = unit(1, "lines"))
}

save_fig(fig_robust_bii(), "sim_robust_bii", height = 4.2)

# --- Felles: fasekomponenter per design fra .rds (per replikasjon -> snitt) -------
# EG: kohortdel = var_uncal_ph1, samplingdel = var_uncal - var_uncal_ph1.
# Robust: samplingdel = var_robust_coxph - var_uncal_ph1 (eq:eg-vs-robust).
# Borgan II: samplingdel = var_borgan_ph2, kohortdel = var_borgan - var_borgan_ph2.
# Elementene kobles til K og rho via Robust-raden i CSV-ene (samme est_var, emp_var, dekning).
read_phases <- function(rds, csv, setting) {
  x <- readRDS(rds); g <- read.csv(csv); rob <- g[g$method == "Robust", ]
  do.call(rbind, lapply(x, function(el) {
    er <- el[el$method == "Robust", ]
    i <- which(abs(rob$est_var - er$est_var) < 1e-9 & abs(rob$emp_var - er$emp_var) < 1e-9 &
               abs(rob$coverage - er$coverage) < 1e-9)
    stopifnot(length(i) == 1)
    r <- as.data.frame(attr(el, "raw"))
    data.frame(setting = setting, K = rob$K[i], rho = sigma_to_rho(rob$strata_noise[i]),
               eg_tot = mean(r$var_uncal), eg_ph1 = mean(r$var_uncal_ph1),
               rob_tot = mean(r$var_robust_coxph),
               b2_tot = mean(r$var_borgan), b2_ph2 = mean(r$var_borgan_ph2))
  }))
}
phases <- rbind(read_phases("results/grid_runs.rds", "results/grid_results.csv", "Samuelsen"),
                read_phases("results/nopw_grid_runs.rds", "results/nopw_grid_results.csv", "NOPW"))
phases$setting <- factor(phases$setting, levels = c("Samuelsen", "NOPW"))
phases$eg_ph2  <- phases$eg_tot - phases$eg_ph1
phases$rob_ph2 <- phases$rob_tot - phases$eg_ph1
phases$share   <- phases$eg_ph2 / phases$eg_tot          # samplingandel (EG)
phases$excess  <- phases$rob_ph2 / phases$eg_ph2 - 1     # robusts overestimat av samplingdelen
phases$share_b2 <- phases$b2_ph2 / phases$b2_tot         # samplingandel (Borgan II)

# --- Figur (3.2.4): gapet = samplingandel x overestimat --------------------------
fig_factor <- function() {
  d <- rbind(data.frame(phases[c("setting", "K", "rho")], measure = "Sampling share", value = phases$share),
             data.frame(phases[c("setting", "K", "rho")], measure = "Excess of the robust sampling component",
                        value = phases$excess))
  d$measure <- factor(d$measure, levels = c("Sampling share", "Excess of the robust sampling component"))
  ggplot(d, aes(K, value, colour = rho, linetype = setting,
                group = interaction(rho, setting))) +
    geom_line(linewidth = 0.7) +
    geom_point(size = 1.8) +
    facet_wrap(~ measure, scales = "free_y") +
    scale_colour_manual(values = rho_cols, name = expression(rho)) +
    scale_linetype_manual(values = c("Samuelsen" = "22", "NOPW" = "solid"), name = NULL) +
    scale_x_continuous(breaks = c(1, 2, 4.2, 6, 10),
                       labels = function(x) format(x, drop0trailing = TRUE)) +
    scale_y_continuous(labels = label_percent()) +
    labs(x = expression(paste("Non-cases sampled per case, ", italic(K))), y = NULL) +
    theme_thesis() +
    theme(panel.spacing = unit(1.2, "lines"))
}

save_fig(fig_factor(), "sim_factor", height = 2.8)

# --- Figur (3.2.4): Borgan II mot EG per replikasjon, kohortdel og samplingdel ------
# Ett design: NOPW, K = 4.2, rho = 1.00. Verdier x 10^3. Diagonal = full enighet.
raw_design <- function(rds, csv, K, sigma) {
  x <- readRDS(rds); g <- read.csv(csv)
  rob <- g[g$method == "Robust" & g$K == K & g$strata_noise == sigma, ]
  for (el in x) {
    er <- el[el$method == "Robust", ]
    if (abs(er$est_var - rob$est_var) < 1e-9 && abs(er$emp_var - rob$emp_var) < 1e-9 &&
        abs(er$coverage - rob$coverage) < 1e-9) return(as.data.frame(attr(el, "raw")))
  }
  stop("fant ikke designet")
}

fig_phases_scatter <- function() {
  r <- raw_design("results/nopw_grid_runs.rds", "results/nopw_grid_results.csv", 4.2, 0)
  d <- data.frame(
    eg_ph1 = r$var_uncal_ph1, eg_ph2 = r$var_uncal - r$var_uncal_ph1,
    b2_ph2 = r$var_borgan_ph2, b2_ph1 = r$var_borgan - r$var_borgan_ph2) * 1e3
  panel <- function(x, y, title) {
    lim <- range(c(d[[x]], d[[y]]))
    ggplot(d, aes(.data[[x]], .data[[y]])) +
      geom_abline(slope = 1, intercept = 0, colour = "grey40", linewidth = 0.4) +
      geom_point(size = 0.4, alpha = 0.25, colour = "#2a78d6") +
      coord_equal(xlim = lim, ylim = lim) +
      labs(title = title, x = expression(paste("Borgan II (", 10^-3, ")")),
           y = expression(paste("EG (", 10^-3, ")"))) +
      theme_thesis() +
      theme(plot.title = element_text(size = 9, face = "bold"))
  }
  panel("b2_ph1", "eg_ph1", "Cohort component") + panel("b2_ph2", "eg_ph2", "Sampling component")
}

save_fig(fig_phases_scatter(), "sim_phases_scatter", height = 3.0)

# --- Figur (3.2.5): ett datasett om gangen - robust mot EG og Borgan II per replikasjon --
# y = robust / designbasert - 1 i prosent, per replikasjon. Boks = kvartiler, strek = median,
# værhår = 5- og 95-persentil (ingen enkeltpunkter). NOPW-settingen.
raw_all <- function(rds, csv, setting) {
  x <- readRDS(rds); g <- read.csv(csv); rob <- g[g$method == "Robust", ]
  do.call(rbind, lapply(x, function(el) {
    er <- el[el$method == "Robust", ]
    i <- which(abs(rob$est_var - er$est_var) < 1e-9 & abs(rob$emp_var - er$emp_var) < 1e-9 &
               abs(rob$coverage - er$coverage) < 1e-9)
    stopifnot(length(i) == 1)
    r <- as.data.frame(attr(el, "raw"))
    data.frame(setting = setting, K = rob$K[i], rho = sigma_to_rho(rob$strata_noise[i]),
               vs_eg = 100 * (r$var_robust_coxph / r$var_uncal - 1),
               vs_b2 = 100 * (r$var_robust_coxph / r$var_borgan - 1))
  }))
}

fig_single <- function() {
  r <- raw_all("results/nopw_grid_runs.rds", "results/nopw_grid_results.csv", "NOPW")
  long <- rbind(data.frame(r[c("K", "rho")], estimator = "EG", diff = r$vs_eg),
                data.frame(r[c("K", "rho")], estimator = "Borgan II", diff = r$vs_b2))
  long$estimator <- factor(long$estimator, levels = c("Borgan II", "EG"))
  s <- aggregate(diff ~ K + rho + estimator, long, function(v)
    quantile(v, c(0.05, 0.25, 0.5, 0.75, 0.95)))
  s <- cbind(s[c("K", "rho", "estimator")], as.data.frame(s$diff))
  names(s)[4:8] <- c("q05", "q25", "q50", "q75", "q95")
  s$Kf <- factor(as.character(s$K), levels = c("1", "2", "4.2", "10"))
  s$panel <- factor(paste0("rho==", "'", s$rho, "'"),
                    levels = paste0("rho==", "'", c("1.00", "0.50", "0.25"), "'"))
  ggplot(s, aes(Kf, fill = estimator, group = interaction(Kf, estimator))) +
    geom_hline(yintercept = 0, colour = "grey40", linewidth = 0.4) +
    geom_boxplot(aes(ymin = q05, lower = q25, middle = q50, upper = q75, ymax = q95),
                 stat = "identity", width = 0.7, linewidth = 0.3, colour = "grey25",
                 position = position_dodge(width = 0.8)) +
    facet_wrap(~ panel, labeller = label_parsed) +
    scale_fill_manual(values = est_cols[c("Borgan II", "EG")], name = "Robust relative to") +
    scale_y_continuous(labels = function(x) paste0(x, "%")) +
    labs(x = expression(paste("Non-cases sampled per case, ", italic(K))),
         y = "Robust / design-based estimate - 1") +
    theme_thesis() +
    theme(panel.spacing = unit(1, "lines"), panel.grid.major.x = element_blank())
}

save_fig(fig_single(), "sim_single", height = 2.9)

# --- Figur (3.3): to korrelerte kovariater mot en kovariat, K = 4.2 ----------------
# Overskudd = robusts samplingledd over EGs - 1; kalibrering = emp. varians kal / design - 1.
# Kap. 3-punktene fra K = 4.2-cellene, 3.3-punktene fra *_two_covariates_runs.rds.
two_cov_stats <- function(r, setting, point) {
  data.frame(setting = setting, point = point,
             excess = mean(r$var_robust_coxph - r$var_uncal_ph1) /
                      mean(r$var_uncal - r$var_uncal_ph1) - 1,
             calibration = var(r$beta_cal) / var(r$beta_uncal) - 1)
}

two_cov_data <- function() {
  one <- function(rds, csv, setting) do.call(rbind, lapply(c(0, 0.5, 1.118), function(s)
    two_cov_stats(raw_design(rds, csv, 4.2, s), setting, as.character(sigma_to_rho(s)))))
  two <- function(rds, setting) {
    raw <- do.call(rbind, lapply(readRDS(rds), function(el) as.data.frame(attr(el, "raw"))))
    do.call(rbind, lapply(c("Z1", "Z2"), function(k) two_cov_stats(raw[raw$coef == k, ], setting, k)))
  }
  d <- rbind(one("results/grid_runs.rds", "results/grid_results.csv", "Samuelsen"),
             two("results/grid_two_covariates_runs.rds", "Samuelsen"),
             one("results/nopw_grid_runs.rds", "results/nopw_grid_results.csv", "NOPW"),
             two("results/nopw_grid_two_covariates_runs.rds", "NOPW"))
  long <- rbind(data.frame(d[c("setting", "point")], measure = "Excess of the robust sampling component",
                           value = d$excess),
                data.frame(d[c("setting", "point")], measure = "Change in variance from calibration",
                           value = d$calibration))
  long$measure <- factor(long$measure, levels = c("Excess of the robust sampling component",
                                                  "Change in variance from calibration"))
  long$setting <- factor(long$setting, levels = c("Samuelsen", "NOPW"))
  long
}

# Stiplet boks = en kovariat med samme marginale rho
# (rho = 1 for beta1, 0.5 for beta2), fylt soyle = to kovariater. Kalibrering som reduksjon.
fig_two_cov_bars <- function() {
  d   <- two_cov_data()
  get <- function(setting, point, measure) d$value[d$setting == setting & d$point == point & d$measure == measure]
  ms  <- levels(d$measure)
  rows <- expand.grid(coef = c("Z1", "Z2"), setting = c("Samuelsen", "NOPW"), stringsAsFactors = FALSE)
  rows$ref  <- c("Z1" = "1.00", "Z2" = "0.50")[rows$coef]
  rows$y    <- c(5, 4, 2, 1)
  rows$lab  <- paste0("paste('", rows$setting, ", ', beta[", substr(rows$coef, 2, 2), "])")
  long <- do.call(rbind, lapply(seq_along(ms), function(j) {
    sgn <- if (j == 2) -1 else 1                      # kalibrering vises som reduksjon
    data.frame(rows, measure = ms[j],
               two = sgn * mapply(get, rows$setting, rows$coef, ms[j]),
               one = sgn * mapply(get, rows$setting, rows$ref, ms[j]))
  }))
  long$measure <- factor(long$measure, levels = ms,
                         labels = c("Excess of the robust\nsampling component (vs. EG)",
                                    "Reduction in variance\nfrom calibration"))
  h <- 0.3
  ggplot(long) +
    geom_vline(xintercept = 0, colour = "grey40", linewidth = 0.3) +
    geom_rect(aes(xmin = 0, xmax = two, ymin = y - h + 0.05, ymax = y + h - 0.05,
                  fill = "Two covariates")) +
    geom_rect(aes(xmin = 0, xmax = one, ymin = y - h, ymax = y + h,
                  linetype = "One covariate, same marginal rho"),
              fill = NA, colour = "grey25", linewidth = 0.4) +
    geom_text(aes(x = pmax(two, 0), y = y, label = label_percent(accuracy = 1)(two)),
              hjust = -0.2, size = 2.5) +
    facet_wrap(~ measure, scales = "free_x") +
    scale_fill_manual(values = c("Two covariates" = "grey45"), name = NULL) +
    scale_linetype_manual(values = c("One covariate, same marginal rho" = "22"), name = NULL,
                          labels = expression(paste("One covariate, same marginal ", rho))) +
    scale_x_continuous(labels = label_percent(), expand = expansion(mult = c(0.02, 0.12))) +
    scale_y_continuous(breaks = rows$y, labels = parse(text = rows$lab), expand = expansion(add = 0.5)) +
    labs(x = NULL, y = NULL) +
    theme_thesis() +
    theme(panel.grid.major.y = element_blank(), panel.spacing = unit(1.2, "lines"),
          strip.text = element_text(size = 9, face = "bold"), legend.position = "bottom")
}

save_fig(fig_two_cov_bars(), "sim_two_cov_bars", height = 2.6)

# ================================================================================
# Kapittel 4: NOPW-dataene på kartet fra kapittel 3
# ================================================================================

# --- Figur (4.1): DAG for hovedmodellen ---------------------------------------------
# Tilpasset fra Babigumira et al. (2023, figur S3), benzen-studien i samme kohort: SES og
# arbeidshistorie er latente; SES -> roeyking og arbeidshistorie; arbeidshistorie -> CSD og
# asbest; roeyking, CSD, asbest (IARC) og alder -> lungekreft. Lagt til for vaar modell:
# fodselskohort -> alder, foerste aar, roeyking; foerste aar -> CSD, asbest og lungekreft
# (periodeavhengige forhold som CSD og asbest ikke fanger). Modellens kovariater er da det
# minimale justeringssettet for CSD og asbest (dagitty::adjustmentSets).
nopw_dag <- dagify(
  lc    ~ age + first + smoke + csd + asb,
  csd   ~ first + work,
  asb   ~ first + work,
  smoke ~ cohort + ses,
  work  ~ ses,
  age   ~ cohort,
  first ~ cohort,
  latent = c("ses", "work"), outcome = "lc",
  labels = c(lc = "Lung\ncancer", age = "Age at\nbaseline", first = "First year\noffshore",
             smoke = "Smoking", csd = "Crystalline\nsilica dust", asb = "Asbestos",
             cohort = "Birth\ncohort", ses = "SES", work = "Work\nhistory"),
  coords = list(x = c(cohort = 0, age = 1.6, first = 1.6, ses = 0, work = 1.6, smoke = 3.2,
                      csd = 3.2, asb = 3.2, lc = 4.8),
                y = c(cohort = 1.4, age = 0, first = 1.4, ses = 3.6, work = 3.0, smoke = 3.6,
                      csd = 2.4, asb = 1.4, lc = 1.8)))

fig_nopw_dag <- function() {
  d <- dplyr::mutate(tidy_dagitty(nopw_dag), latent = name %in% c("ses", "work"))
  ggplot(d, aes(x = x, y = y, xend = xend, yend = yend)) +
    geom_dag_edges(aes(edge_linetype = ifelse(latent, "dashed", "solid")),
                   edge_colour = "grey35", edge_width = 0.4,
                   arrow_directed = grid::arrow(length = grid::unit(5, "pt"), type = "closed"),
                   start_cap = ggraph::ellipsis(9, 5, "mm"), end_cap = ggraph::ellipsis(9, 5, "mm")) +
    geom_label(aes(label = label, colour = latent), size = 2.8, lineheight = 0.9,
               linewidth = 0.3, label.r = grid::unit(3, "pt"), fill = "white") +
    scale_colour_manual(values = c("FALSE" = "grey15", "TRUE" = "grey55"), guide = "none") +
    scale_x_continuous(expand = expansion(add = 0.5)) +
    scale_y_continuous(expand = expansion(add = 0.4)) +
    theme_void()
}

save_fig(fig_nopw_dag(), "nopw_dag", height = 2.8)

# --- Figur (4.2/4.3): samplingandel og overskudd per koeffisient ------------------
# Hardkodet fra main.tex: EG-SE og robust/EG - 1 (tab:variance-comparison, regnet på
# urundede tall) og EGs kohortledd (tab:variance-decomp). Samplingandel s = 1 - kohort/EG-varians;
# overskudd = (robust/EG - 1) / s, siden robust - EG = robusts samplingledd - EGs samplingledd.
nopw_coefs <- data.frame(
  coef = c("Age at baseline", "First year offshore",
           "CSD, quartile 1", "CSD, quartile 2", "CSD, quartile 3", "CSD, quartile 4",
           "Current smoker, < median", "Current smoker, > median", "Former smoker, > 15 years",
           "Cumulative asbestos dose"),
  se_eg   = c(0.00644, 0.01153, 0.18224, 0.18837, 0.18333, 0.17591, 0.24040, 0.23442, 0.24285, 0.00339),
  diff_eg = c(11.54, 0.61, -0.24, -0.19, -0.12, -0.26, 0.19, 0.19, 0.08, -0.26) / 100,
  eg_ph1  = c(3.31e-5, 9.95e-5, 2.37e-2, 2.55e-2, 2.33e-2, 2.05e-2, 5.34e-2, 5.02e-2, 5.58e-2, 8.31e-6),
  # Borgan II fra tab:variance-decomp (cch): kohortledd og samplingledd
  b2_ph1  = c(3.49e-5, 8.59e-5, 2.33e-2, 2.53e-2, 2.30e-2, 2.07e-2, 5.38e-2, 5.05e-2, 5.46e-2, 8.58e-6),
  b2_ph2  = c(8.41e-6, 3.32e-5, 9.45e-3, 9.92e-3, 1.026e-2, 1.043e-2, 4.40e-3, 4.79e-3, 3.135e-3, 3.16e-6)
)
nopw_coefs$share  <- 1 - nopw_coefs$eg_ph1 / nopw_coefs$se_eg^2
nopw_coefs$excess <- nopw_coefs$diff_eg / nopw_coefs$share
nopw_coefs$share_b2 <- nopw_coefs$b2_ph2 / (nopw_coefs$b2_ph1 + nopw_coefs$b2_ph2)
nopw_coefs$coef   <- factor(nopw_coefs$coef, levels = rev(nopw_coefs$coef))

fig_nopw_factor <- function() {
  # Venstre: stiplet boks = hele variansen (100 %); to søyler per rad, EG (grønn, øverst) og
  # Borgan II (oransje, nederst) - samme farger som estimatorene i kap. 3.
  # Høyre: overskuddet (robust mot EG) som søyle. Simuleringen (NOPW-settingen, K = 4.2) som
  # tre egne rader nederst; i høyre panel i rho-fargene fra kap. 3.
  ref <- phases[phases$setting == "NOPW" & phases$K == 4.2, c("rho", "share", "share_b2", "excess")]
  ref <- ref[order(ref$rho), ]                       # 1.00, 0.50, 0.25
  d <- rbind(
    data.frame(y = as.numeric(nopw_coefs$coef) + 3.8, share = nopw_coefs$share,
               share_b2 = nopw_coefs$share_b2, excess = nopw_coefs$excess, grp = "NOPW",
               lab = paste0("'", as.character(nopw_coefs$coef), "'")),
    data.frame(y = 3:1, share = ref$share, share_b2 = ref$share_b2, excess = ref$excess,
               grp = as.character(ref$rho),
               lab = paste0("paste('Simulation, ', rho, ' = ", ref$rho, "')")))
  ylab <- scale_y_continuous(breaks = d$y, labels = parse(text = d$lab),
                             expand = expansion(add = 0.5))
  h <- 0.3
  base <- function(title) list(
    labs(title = title, x = NULL, y = NULL),
    theme_thesis(),
    theme(plot.title = element_text(size = 9, face = "bold"), panel.grid.major.y = element_blank()))
  left <- ggplot(d) +
    geom_rect(aes(xmin = 0, xmax = share, ymin = y + 0.02, ymax = y + h - 0.03, fill = "EG")) +
    geom_rect(aes(xmin = 0, xmax = share_b2, ymin = y - h + 0.03, ymax = y - 0.02, fill = "Borgan II")) +
    geom_rect(aes(xmin = 0, xmax = 1, ymin = y - h, ymax = y + h,
                  linetype = "Total variance"), fill = NA, colour = "grey25", linewidth = 0.4) +
    scale_fill_manual(values = est_cols[c("EG", "Borgan II")], breaks = c("EG", "Borgan II"), name = NULL) +
    scale_linetype_manual(values = c("Total variance" = "22"), name = NULL) +
    scale_x_continuous(labels = label_percent(), limits = c(0, 1), breaks = seq(0, 1, 0.25),
                       expand = expansion(mult = c(0.01, 0.03))) +
    ylab + base("Sampling share\nof the variance") +
    guides(fill = guide_legend(order = 1), linetype = guide_legend(order = 2))
  right <- ggplot(d) +
    geom_rect(aes(xmin = 0, xmax = excess, ymin = y - h + 0.05, ymax = y + h - 0.05, fill = grp)) +
    geom_vline(xintercept = 0, colour = "grey40", linewidth = 0.3) +
    scale_fill_manual(values = c("NOPW" = "grey45", rho_cols), guide = "none") +
    scale_x_continuous(labels = label_percent(), expand = expansion(mult = c(0.05, 0.08))) +
    ylab + base("Excess of the robust\nsampling component (vs. EG)") +
    theme(axis.text.y = element_blank())
  (left | right) + plot_layout(guides = "collect") & theme(legend.position = "bottom")
}

save_fig(fig_nopw_factor(), "nopw_factor", height = 4.2)

# --- Figur (4.5): kalibreringsgevinst mot taket ------------------------------------
# Hardkodet fra tab:calibration-decomp (trinnvis kjøring, variansskala): kohortledd og
# samplingledd med design- og kalibrerte vekter. Tak = samplingandelen med designvekter.
# Reduksjon i standardfeil i % av design-SE, delt på de to leddene; stiplet boks = taket.
# (Andel av taket mot simuleringen, cal_ref(), er droppet fra figuren etter brukerens valg.)
nopw_cal <- data.frame(
  coef   = levels(nopw_coefs$coef)[10:1],
  ph1_d  = c(3.30e-5, 9.95e-5, 2.37e-2, 2.55e-2, 2.33e-2, 2.05e-2, 5.34e-2, 5.02e-2, 5.58e-2, 8.30e-6),
  ph1_c  = c(3.24e-5, 1.02e-4, 2.38e-2, 2.54e-2, 2.24e-2, 2.01e-2, 5.34e-2, 5.09e-2, 5.65e-2, 7.58e-6),
  ph2_d  = c(8.51e-6, 3.32e-5, 9.44e-3, 9.91e-3, 1.026e-2, 1.035e-2, 4.40e-3, 4.78e-3, 3.13e-3, 3.15e-6),
  ph2_c  = c(6.26e-6, 3.40e-5, 9.30e-3, 9.20e-3, 9.34e-3, 8.52e-3, 2.54e-3, 2.69e-3, 2.25e-3, 2.32e-6)
)
# SE-skala (som tab:calibration-gain): reduksjon = 1 - sqrt(kal. varians / designvarians),
# tak = 1 - sqrt(1 - samplingandel) = hele samplingleddet fjernet, kohortleddet uendret.
# Samplingdelen regnes på samme måte (kohortleddet holdt på designnivå), så den kan leses
# direkte mot taket; kohortdelen er resten (endring i det *estimerte* kohortleddet).
se_red <- function(a, b, tot) 1 - sqrt((a + b) / tot)
nopw_cal$tot_d   <- nopw_cal$ph1_d + nopw_cal$ph2_d
nopw_cal$ceiling <- with(nopw_cal, se_red(ph1_d, 0, tot_d))
nopw_cal$red_tot <- with(nopw_cal, se_red(ph1_c, ph2_c, tot_d))
nopw_cal$red_ph2 <- with(nopw_cal, se_red(ph1_d, ph2_c, tot_d))
nopw_cal$red_ph1 <- nopw_cal$red_tot - nopw_cal$red_ph2
nopw_cal$frac    <- nopw_cal$red_tot / nopw_cal$ceiling
nopw_cal$coef    <- factor(nopw_cal$coef, levels = levels(nopw_coefs$coef))

cal_ref <- function() {
  key <- c("setting", "K", "rho")
  des <- grid[grid$method == "uncalibrated", c(key, "emp_var")]
  cal <- grid[grid$method == "calibrated",   c(key, "emp_var")]
  m <- merge(merge(des, cal, by = key, suffixes = c("_des", "_cal")), phases[c(key, "share")], by = key)
  m <- m[m$setting == "NOPW" & m$K == 4.2, ]
  data.frame(rho = m$rho, frac = (1 - m$emp_var_cal / m$emp_var_des) / m$share)
}

fig_nopw_calibration <- function() {
  # Grønn søyle = reduksjonen fra samplingleddet (kohortleddet holdt på designnivå), inne i den
  # stiplede boksen = taket. Endringen i det estimerte kohortleddet er utelatt (brukerens valg;
  # tas i teksten og tab:calibration-decomp). Tre rader nederst = simuleringen (NOPW, K = 4.2):
  # SE-reduksjon fra empirisk varians, tak = 1 - sqrt(1 - snittandel).
  d <- data.frame(y = as.numeric(nopw_cal$coef) + 3.8, red = nopw_cal$red_ph2,
                  ceiling = nopw_cal$ceiling, fill = "NOPW",
                  lab = paste0("'", as.character(nopw_cal$coef), "'"))
  key <- c("setting", "K", "rho")
  des <- grid[grid$method == "uncalibrated", c(key, "emp_var")]
  cal <- grid[grid$method == "calibrated",   c(key, "emp_var")]
  m <- merge(merge(des, cal, by = key, suffixes = c("_des", "_cal")), phases[c(key, "share")], by = key)
  m <- m[m$setting == "NOPW" & m$K == 4.2, ]
  m <- m[order(m$rho), ]                                # 1.00, 0.50, 0.25
  sim <- data.frame(y = 3:1, red = 1 - sqrt(m$emp_var_cal / m$emp_var_des),
                    ceiling = 1 - sqrt(1 - m$share), fill = as.character(m$rho),
                    lab = paste0("paste('Simulation, ', rho, ' = ", m$rho, "')"))
  all <- rbind(d, sim)
  h <- 0.3
  ggplot(all) +
    geom_vline(xintercept = 0, colour = "grey40", linewidth = 0.3) +
    geom_rect(aes(xmin = pmin(red, 0), xmax = pmax(red, 0), ymin = y - h, ymax = y + h, fill = fill)) +
    geom_rect(aes(xmin = 0, xmax = ceiling, ymin = y - h, ymax = y + h,
                  linetype = "Ceiling: the whole sampling component removed"),
              fill = NA, colour = "grey25", linewidth = 0.4) +
    scale_fill_manual(values = c("NOPW" = "#1baf7a", rho_cols), guide = "none") +
    scale_linetype_manual(values = c("Ceiling: the whole sampling component removed" = "22"), name = NULL) +
    scale_y_continuous(breaks = all$y, labels = parse(text = all$lab), expand = expansion(add = 0.5)) +
    scale_x_continuous(labels = label_percent(), expand = expansion(mult = c(0.05, 0.08))) +
    labs(x = "Reduction in standard error from the sampling component", y = NULL) +
    theme_thesis() +
    theme(panel.grid.major.y = element_blank(), legend.position = "bottom")
}

save_fig(fig_nopw_calibration(), "nopw_calibration", height = 3.8)

# --- Figur (4.6): ren risiko per profil, ueksponert mot Q4 ---------------------------------
# Hardkodet fra tab:pure-risk (designvekter, [0, 22.5], per 1000, KI på log-skala).
# Teksten til høyre = absolutt økning fra ueksponert til Q4 (samme HR 1.48 i alle profilene).
nopw_pr <- data.frame(
  profile = rep(c("Low", "Medium", "High"), each = 2),
  csd     = rep(c("Unexposed", "Highest quartile (Q4)"), 3),
  pr      = c(6.5, 9.6, 34.8, 51.0, 125.0, 178.7),
  lo      = c(4.1, 5.7, 27.2, 36.1, 101.7, 128.9),
  hi      = c(10.2, 16.1, 44.5, 72.0, 153.6, 247.9)
)
prof_lab <- c(Low = "Low: never smoker or\nformer, < 15 years", Medium = "Medium: former\nsmoker, > 15 years",
              High = "High: current smoker,\n> median pack-years")
nopw_pr$profile <- factor(prof_lab[nopw_pr$profile], levels = prof_lab[c("High", "Medium", "Low")])
nopw_pr$csd     <- factor(nopw_pr$csd, levels = c("Unexposed", "Highest quartile (Q4)"))

fig_nopw_pure_risk <- function() {
  csd_cols <- c("Unexposed" = "#2a78d6", "Highest quartile (Q4)" = "#eb6834")
  dodge <- position_dodge(width = 0.5)
  # Forskjellene fra teksten i 4.6 (urundede tall, nopw/pure_risk/pure_risk_driver.r), ikke fra de avrundede pr.
  gap <- data.frame(profile = factor(levels(nopw_pr$profile), levels = levels(nopw_pr$profile)),
                    diff = c(53.8, 16.1, 3.1))
  gap$x <- tapply(nopw_pr$hi, nopw_pr$profile, max)[as.character(gap$profile)]
  ggplot(nopw_pr, aes(pr, profile, colour = csd)) +
    geom_errorbar(aes(xmin = lo, xmax = hi), width = 0, orientation = "y", linewidth = 0.6, position = dodge) +
    geom_point(size = 2.2, position = dodge) +
    geom_text(data = gap, aes(x = x, y = profile, label = sprintf("+%.1f per 1000", diff)),
              inherit.aes = FALSE, hjust = 0, nudge_x = 6, size = 3, colour = "grey25") +
    scale_colour_manual(values = csd_cols, name = "Crystalline silica dust") +
    scale_x_continuous(limits = c(0, 310), breaks = seq(0, 250, 50)) +
    labs(x = "Pure risk of lung cancer over 22.5 years, per 1000 workers", y = NULL) +
    theme_thesis() +
    theme(panel.grid.major.y = element_blank())
}

save_fig(fig_nopw_pure_risk(), "nopw_pure_risk", height = 2.6)
