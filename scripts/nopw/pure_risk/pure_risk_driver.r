##=============================================================================
## Project: Master's thesis - variance estimation in the NOPW case-cohort study
## Script Purpose: Estimate pure risk of lung cancer over (0, 22.5] for selected
##                 covariate profiles (Section 4.5) with the driver function of
##                 CaseCohortCoxSurvival and the predefined design weights. The
##                 package returns the pure risk with the EG variance and the
##                 robust variance.
## Date: 14.10.2026
## Author: Steinar Jonathan Mpagi
##=============================================================================
## Load datasets and libraries
getwd()
setwd("<DATA PATH>")
load("data/case_cohort_cox_survival_data.rda")
library(CaseCohortCoxSurvival)
library(dplyr)
library(haven)
##-----------------------------------------------------------------------------
# set weights to match new weights
w = full_cohort_and_subcohort$weight
full_cohort_and_subcohort$weight = ifelse(w == 1.838, 1.847,
                                          ifelse(w == 4.351, 4.346, full_cohort_and_subcohort$weight))
##-----------------------------------------------------------------------------
# 1 - data set up, refactor, add columns, change names needed to run the analysis later.
# As factor and set reference level
full_cohort_and_subcohort$Smoke12 = as.factor(full_cohort_and_subcohort$Smoke12)
full_cohort_and_subcohort = within(full_cohort_and_subcohort, Smoke12 <- relevel(Smoke12, ref = "Never/former, dur<15"))

# CSD as factor
full_cohort_and_subcohort$CSD_CUM_quantile = as_factor(full_cohort_and_subcohort$CSD_CUM_quantile)
full_cohort_and_subcohort$CSD_CUM_quantile = relevel(full_cohort_and_subcohort$CSD_CUM_quantile, "Q0")

#set rename levels to some arbitrary name to be able to run model - wont runt with special characters
# The reference is named "0" and not "D": for pure risk the driver turns the
# profile values into a factor again in alphabetical order, so the reference
# must sort first. With "D", a profile with "B" is computed as the reference
# (warning: "columns 'Smoke12D' that are not covariates"). Q0 already sorts first.
# 0 = never or former, duration <= 15 years (reference)
# A = current smoker, < median pack-years   (HR ~ 9.9)
# B = current smoker, > median pack-years   (HR ~ 20.5)
# C = former smoker, duration > 15 years    (HR ~ 5.4)
full_cohort_and_subcohort = full_cohort_and_subcohort |>
  mutate(Smoke12 = if_else(Smoke12 == "Current <median packyrs", "A", Smoke12)) |>
  mutate(Smoke12 = if_else(Smoke12 == "Current >median packyrs", "B", Smoke12)) |>
  mutate(Smoke12 = if_else(Smoke12 == "Former & duration>15", "C", Smoke12)) |>
  mutate(Smoke12 = if_else(Smoke12 == "Never/former, dur<15", "0", Smoke12)) |>
  mutate(Smoke12 = factor(Smoke12)) |>
  mutate(Smoke12 = Smoke12 <- relevel(Smoke12, ref = "0"))

# Set factor levels for stilling_20 + add missing + create indicator for phase2
full_cohort_and_subcohort$phase2 = ifelse(full_cohort_and_subcohort$Subcohort_indicator == 1 | full_cohort_and_subcohort$Case == 1, 1, 0)
full_cohort_and_subcohort = full_cohort_and_subcohort |> mutate(stilling_20 = as_factor(stilling_20))

full_cohort_and_subcohort = within(full_cohort_and_subcohort, Smoke12 <- relevel(Smoke12, ref = "0"))

# the data used for the analysis, here the 81 members without data on Smoke12 are removed.
case_cohort = full_cohort_and_subcohort %>% filter(!((Case == 1 | Subcohort_indicator == 1) & is.na(Smoke12)))

# needed to run uncalibrated driver model. Can't have NA eventough the values set here are not used.
case_cohort$weights = ifelse(is.na(case_cohort$weight), 1, case_cohort$weight)

full_cohort_and_subcohort$weights = ifelse(is.na(full_cohort_and_subcohort$weight), 1, full_cohort_and_subcohort$weight)

################################################################################
# 2 - pure risk using the driver function with design weights.
# Ties in time_to_event are broken at random, so set a seed for reproducibility.
set.seed(1999)

# Covariate profiles. Table profiles (tab:pure-risk): age 55 and first year 1979
# (medians among the 2409 analysed), asbestos 0; smoking at three levels, CSD Q0
# and Q4. Extremes (text in Section 4.5.2): age, first year and asbestos at their
# quartiles. Made with data.frame(), not cbind(), which turns every column into
# text. When several profiles are run together, the reference levels ("0" and
# "Q0") must be among them, otherwise the first level present is taken as reference.
model_reference_worker = data.frame(
  age_at_baseline = 55,
  first_year = 1979,
  Smoke12 = "0",
  CSD_CUM_quantile = "Q0",
  Asb_CumDose = 0)

model_Q4CSD_worker = data.frame(
  age_at_baseline = 55,
  first_year = 1979,
  Smoke12 = "0",
  CSD_CUM_quantile = "Q4",
  Asb_CumDose = 0)

model_BigSmoker_worker = data.frame(
  age_at_baseline = 55,
  first_year = 1979,
  Smoke12 = "B",
  CSD_CUM_quantile = "Q0",
  Asb_CumDose = 0)

model_BigSmoker_Q4CSD_worker = data.frame(
  age_at_baseline = 55,
  first_year = 1979,
  Smoke12 = "B",
  CSD_CUM_quantile = "Q4",
  Asb_CumDose = 0)

model_exSmoker_worker = data.frame(
  age_at_baseline = 55,
  first_year = 1979,
  Smoke12 = "C",
  CSD_CUM_quantile = "Q0",
  Asb_CumDose = 0)

model_exSmoker_Q4CSD_worker = data.frame(
  age_at_baseline = 55,
  first_year = 1979,
  Smoke12 = "C",
  CSD_CUM_quantile = "Q4",
  Asb_CumDose = 0)

model_bestCase_worker = data.frame(
  age_at_baseline = 49,
  first_year = 1985,
  Smoke12 = "0",
  CSD_CUM_quantile = "Q0",
  Asb_CumDose = 0.6)

model_WorstCase_worker = data.frame(
  age_at_baseline = 60,
  first_year = 1977,
  Smoke12 = "B",
  CSD_CUM_quantile = "Q4",
  Asb_CumDose = 17.4)

# data must be a data.frame: with a tibble the pure risk step stops with
# "INTERNAL CODING ERROR 1" and Pi.var is NULL.
cccs_design_weights = caseCohortCoxSurvival(
  data = as.data.frame(case_cohort),
  time = "time_to_event",
  weights.phase2 = "weights",
  status = "Case",
  cox.phase2 = c("age_at_baseline" ,
                 "Smoke12",
                 "first_year",
                 "CSD_CUM_quantile",
                 "Asb_CumDose"
  ),
  Tau1 = 0,
  Tau2 = max(full_cohort_and_subcohort$time_to_event),   # = 22.5
  x = rbind(model_reference_worker,
            model_Q4CSD_worker,
            model_BigSmoker_worker,
            model_BigSmoker_Q4CSD_worker,
            model_exSmoker_worker,
            model_exSmoker_Q4CSD_worker
            ,model_WorstCase_worker,
            model_bestCase_worker),
  subcohort = "Subcohort_indicator",
  strata = "bcohort",
  print = 1
)
cccs_design_weights$Pi.var
## Output (run on the VM, 14.10.2026, set.seed(1999)):
##          risk     variance robust.variance
## 1 0.006497416 2.242987e-06    2.259181e-06    # Low,    Q0
## 2 0.009569461 6.398404e-06    6.426337e-06    # Low,    Q4
## 3 0.124959282 1.726389e-04    1.784504e-04    # High,   Q0
## 4 0.178728980 8.903429e-04    8.991659e-04    # High,   Q4
## 5 0.034834130 1.904805e-05    1.949773e-05    # Medium, Q0
## 6 0.050955652 8.053815e-05    8.131738e-05    # Medium, Q4
## 7 0.287431213 1.821573e-03    1.840149e-03    # Worst case
## 8 0.003032605 4.670661e-07    4.733367e-07    # Best case
##-----------------------------------------------------------------------------
# 3 - per 1000 workers, with 95% CI on the log scale (EG, Section 7.2):
#     exp{log PR +- 1.96 SE / PR}
pr = cccs_design_weights$Pi.var
pr$se = sqrt(pr$variance)
pr$se_robust = sqrt(pr$robust.variance)
pr$lower = pr$risk * exp(-1.96 * pr$se / pr$risk)
pr$upper = pr$risk * exp( 1.96 * pr$se / pr$risk)

profile = c("Low Q0", "Low Q4", "High Q0", "High Q4", "Medium Q0", "Medium Q4", "Worst case", "Best case")
cbind(profile, round(1000 * pr[, c("risk", "se", "lower", "upper", "se_robust")], 1))
## Computed from the Pi.var output above with the code in this step:
##      profile  risk   se lower upper se_robust
##       Low Q0   6.5  1.5   4.1  10.2       1.5
##       Low Q4   9.6  2.5   5.7  16.1       2.5
##      High Q0 125.0 13.1 101.7 153.6      13.4
##      High Q4 178.7 29.8 128.9 247.9      30.0
##    Medium Q0  34.8  4.4  27.2  44.5       4.4
##    Medium Q4  51.0  9.0  36.1  72.0       9.0
##   Worst case 287.4 42.7 214.9 384.5      42.9
##    Best case   3.0  0.7   1.9   4.7       0.7
##-----------------------------------------------------------------------------
