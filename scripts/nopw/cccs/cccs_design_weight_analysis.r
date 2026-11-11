##=============================================================================
## Project: Master's thesis - variance estimation in the NOPW case-cohort study
## Script Purpose: Fit the main Cox model for lung cancer (Chapter 4) with the
##                 driver function of CaseCohortCoxSurvival and the predefined
##                 design weights. The package returns the EG variance and the
##                 robust variance.
## Date: 04.08.2026
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
# 1 - data set up, refactor, add clumns, changa names needed to run the analyssis later.
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
##-----------------------------------------------------------------------------
# Set factor levels for stilling_20 + add missing + create indicator for phase2
full_cohort_and_subcohort$phase2 = ifelse(full_cohort_and_subcohort$Subcohort_indicator == 1 | full_cohort_and_subcohort$Case == 1, 1, 0)
full_cohort_and_subcohort = full_cohort_and_subcohort |> mutate(stilling_20 = as_factor(stilling_20))

full_cohort_and_subcohort = within(full_cohort_and_subcohort, Smoke12 <- relevel(Smoke12, ref = "D"))

# the data used for the analysis, here the 81 membes without NA on smoke12 are removed.
case_cohort = full_cohort_and_subcohort %>% filter(!((Case == 1 | Subcohort_indicator == 1) & is.na(Smoke12)))

# needed to run uncalibrated driver model. Can't have NA eventough the values set here are not used.
case_cohort$weights = ifelse(is.na(case_cohort$weight), 1, case_cohort$weight)

full_cohort_and_subcohort$weights = ifelse(is.na(full_cohort_and_subcohort$weight), 1, full_cohort_and_subcohort$weight)

################################################################################
# 2 - analysis using the driver function with design weights.

max(full_cohort_and_subcohort$time_to_event)
cccs_design_weights = caseCohortCoxSurvival(
  data = case_cohort,
  time = "time_to_event",
  weights.phase2 = "weights",
  status = "Case",
  cox.phase2 = c("age_at_baseline",
                 "Smoke12",
                 "first_year",
                 "CSD_CUM_quantile",
                 "Asb_CumDose"
  ),
  subcohort = "Subcohort_indicator",
  strata = "bcohort",
  print = 3
)

## Note: cccs orders the CSD quartiles as Q1, Q2, Q4, Q3 in its output.
cccs_design_weights$beta
## Output (run on the VM):
## age_at_baseline       Smoke12A        Smoke12B        Smoke12C      first_year
##      0.10817053     2.29265371      3.01931241      1.69363645     -0.01891411
## CSD_CUM_quantileQ1 CSD_CUM_quantileQ2 CSD_CUM_quantileQ4 CSD_CUM_quantileQ3 Asb_CumDose
##         0.09625404        -0.14208926         0.38856303         0.02977941 -0.00205458

## EG standard errors
sqrt(cccs_design_weights$beta.var)
## Output (run on the VM):
## age_at_baseline       Smoke12A        Smoke12B        Smoke12C      first_year
##     0.006442668    0.240405714     0.234420525     0.242850198     0.011526044
## CSD_CUM_quantileQ1 CSD_CUM_quantileQ2 CSD_CUM_quantileQ4 CSD_CUM_quantileQ3 Asb_CumDose
##        0.182241611        0.188382846        0.175916437        0.183337039 0.003388580

## Robust standard errors (from the coxph fit inside the package)
sqrt(cccs_design_weights$beta.robustvar)
## Output (run on the VM):
## age_at_baseline       Smoke12A        Smoke12B        Smoke12C      first_year
##     0.006807985    0.240619767     0.234645760     0.242931659     0.011562003
## CSD_CUM_quantileQ1 CSD_CUM_quantileQ2 CSD_CUM_quantileQ4 CSD_CUM_quantileQ3 Asb_CumDose
##        0.182078668        0.188216409        0.175710430        0.183226553 0.003386676
