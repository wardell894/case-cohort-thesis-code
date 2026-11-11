##=============================================================================
## Project: Master's thesis - variance estimation in the NOPW case-cohort study
## Script Purpose: Fit the main Cox model for lung cancer (Chapter 4) with the
##                 driver function of CaseCohortCoxSurvival and calibrated
##                 weights (calibrated = TRUE). The package imputes from the
##                 proxies, calibrates the design weights and returns the EG
##                 variance.
## Date:
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
# A = current smoker, < median pack-years   (HR ~ 9.9)
# B = current smoker, > median pack-years   (HR ~ 20.5)
# C = former smoker, duration > 15 years    (HR ~ 5.4)
# D = never or former, duration <= 15 years (reference)
full_cohort_and_subcohort = full_cohort_and_subcohort |>
  mutate(Smoke12 = if_else(Smoke12 == "Current <median packyrs", "A", Smoke12)) |>
  mutate(Smoke12 = if_else(Smoke12 == "Current >median packyrs", "B", Smoke12)) |>
  mutate(Smoke12 = if_else(Smoke12 == "Former & duration>15", "C", Smoke12)) |>
  mutate(Smoke12 = if_else(Smoke12 == "Never/former, dur<15", "D", Smoke12)) |>
  mutate(Smoke12 = factor(Smoke12)) |>
  mutate(Smoke12 = Smoke12 <- relevel(Smoke12, ref = "D"))

# Set factor levels for stilling_20 + add missing + create indicator for phase2
full_cohort_and_subcohort$phase2 = ifelse(full_cohort_and_subcohort$Subcohort_indicator == 1 | full_cohort_and_subcohort$Case == 1, 1, 0)
full_cohort_and_subcohort = full_cohort_and_subcohort |> mutate(stilling_20 = as_factor(stilling_20))
full_cohort_and_subcohort$stilling_20 |> as.factor()
full_cohort_and_subcohort = full_cohort_and_subcohort |> mutate(stilling_20 = if_else(is.na(stilling_20), "Missing", stilling_20))

# strip the different titel to the core to reduce the amount of levels
full_cohort_and_subcohort$stilling_20_trimmed = trimws(
  sub("\\s*(feilstavelse|flere titler|leder|lokalisasjon|lokasjon|samme)\\s*$", "",
      full_cohort_and_subcohort$stilling_20, ignore.case = TRUE))

full_cohort_and_subcohort$weight = ifelse(is.na(full_cohort_and_subcohort$weight), 1, full_cohort_and_subcohort$weight)
full_cohort_and_subcohort$weights = full_cohort_and_subcohort$weight
full_cohort_and_subcohort = full_cohort_and_subcohort %>%
  mutate(phase2 = as.integer(Subcohort_indicator == 1 | Case == 1)) %>%
  mutate(phase2_no_missing = as.integer(phase2 == 1 & !is.na(Smoke12))) %>%
  mutate(Case_driver = ifelse(phase2_no_missing == 1, Case, 0)) %>%
  mutate(Subcohort_indicator_driver = ifelse(phase2_no_missing == 1, Subcohort_indicator, 0)) %>%
  mutate(full_cohort_no_missing = as.integer(!(phase2 == 1 & is.na(Smoke12))))

nopw_cohort = full_cohort_and_subcohort %>% filter(!(phase2 == 1 & is.na(Smoke12)))
nopw_cohort = as.data.frame(nopw_cohort)
for (v in c("typeselsk_1_sp23", "stilling_20_trimmed", "bcohort", "boring_2_sp91", "royk30_sp37b3", "losemdamp_sp2915")) nopw_cohort[[v]] = as_factor(nopw_cohort[[v]])

set.seed(1999)

cccs_calibrated_weights = caseCohortCoxSurvival(
  data = nopw_cohort,
  time = "time_to_event",
  weights.phase2 = "weight",
  status = "Case",
  cox.phase2 = c("age_at_baseline",
                 "Smoke12",
                 "first_year",
                 "CSD_CUM_quantile",
                 "Asb_CumDose"
  ),
  aux.method = "Breslow",
  calibrated = TRUE,
  predictors.cox.phase2 = list(age_at_baseline = c("fodselsaar_fr", "typeselsk_1_sp23"),
                               first_year = c("fodselsaar_fr", "bcohort"),
                               Asb_CumDose = c("stilling_20_trimmed"),
                               CSD_CUM_quantile = c("stilling_20_trimmed", "losemdamp_sp2915", "bcohort", "boring_2_sp91"),
                               Smoke12 = c("royk30_sp37b3")
  ),
  subcohort = "Subcohort_indicator",
  strata = "bcohort",
  print = 1
)

cccs_calibrated_weights$beta
## Output (run on the VM):
## age_at_baseline       Smoke12A        Smoke12B        Smoke12C      first_year
##    0.1033261325   2.3710019629    3.0257367898    1.7300880844   -0.0178535288
## CSD_CUM_quantileQ1 CSD_CUM_quantileQ2 CSD_CUM_quantileQ4 CSD_CUM_quantileQ3   Asb_CumDose
##       0.1071811281      -0.1355182562       0.3753033929       0.0444804839 -0.0006877991

sqrt(cccs_calibrated_weights$beta.var)
## Output (run on the VM):
## age_at_baseline       Smoke12A        Smoke12B        Smoke12C      first_year
##     0.006226837    0.236612806     0.231439832     0.242408663     0.011670663
## CSD_CUM_quantileQ1 CSD_CUM_quantileQ2 CSD_CUM_quantileQ4 CSD_CUM_quantileQ3   Asb_CumDose
##        0.182045159        0.186036125        0.169434740        0.178359785   0.003148739
##-----------------------------------------------------------------------------
