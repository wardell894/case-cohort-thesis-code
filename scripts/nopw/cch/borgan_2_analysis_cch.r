##=============================================================================
## Project: Master's thesis - variance estimation in the NOPW case-cohort study
## Script Purpose: Fit the main Cox model for lung cancer (Chapter 4) with cch()
##                 (method "II.Borgan", stratified by five-year birth cohort)
##                 and compute the Borgan II variance and its decomposition
##                 into cohort and sampling components (phase2var).
## Note: cch() computes post-stratified weights itself from cohort.size. The
##       stratum sizes passed in cohort.size are therefore set (w * m + d per
##       stratum) so that cch() reproduces the predefined design weights used
##       by coxph() and CaseCohortCoxSurvival. The two birth cohorts without
##       cases (1915-19, 1975-79) are left out of the fit.
## Date: 25.09.2026
## Author: Steinar Jonathan Mpagi
##=============================================================================
## Load datasets and libraries
getwd()
setwd("<DATA PATH>")
library(survival)
library(dplyr)
load("data/analysis/full_cohort_and_subcohort.RData")
lung_data_with_subcohort_indicator = full_cohort_and_subcohort %>% filter(Case == 1 | Subcohort_indicator == 1)
##-----------------------------------------------------------------------------
# smoke12 as a factor 
lung_data_with_subcohort_indicator =  lung_data_with_subcohort_indicator |> filter(!is.na(Smoke12)) 
lung_data_with_subcohort_indicator$Smoke12 = as.factor(lung_data_with_subcohort_indicator$Smoke12)
lung_data_with_subcohort_indicator = within(lung_data_with_subcohort_indicator, Smoke12 <- relevel(Smoke12, ref = "Never/former, dur<15"))

# manipulated stratasizes to match the weights from the weighting scheme (excluding strata with no cases)
birth_cohort_values = c("2" = 20, "3" = 120, "4" = 372.695, "5" = 886.543, "6" = 1986.818,
                        "7" = 3096.096, "8" = 4017.780, "9" = 4728.210, "10" = 4203.917,
                        "11" = 2789.750, "12" = 1189)

# remove strata with no cases as it is needed to run the cch analysis.
lung_data_with_subcohort_indicator = lung_data_with_subcohort_indicator |> filter(!bcohort %in% c("1", "13"))

# naming in the same order as the dataframe
lev = levels(factor(lung_data_with_subcohort_indicator$bcohort))

# make it a table to use as the stratified command
birth_cohort_strata = as.table(birth_cohort_values[lev])

table(lung_data_with_subcohort_indicator$bcohort)
birth_cohort_strata
##-----------------------------------------------------------------------------
# run borgan model
borgan2_model = cch(Surv(time_to_event, Case)
                    ~ age_at_baseline +
                      Smoke12 +
                      first_year +
                      CSD_CUM_quantile +
                      Asb_CumDose,
                    data = lung_data_with_subcohort_indicator,
                    subcoh = ~Subcohort_indicator,
                    id = ~id,
                    stratum = ~bcohort,
                    cohort.size = birth_cohort_strata,
                    method = "II.Borgan")

summary(borgan2_model)
summary(borgan2_model$weights)
## Output (run on the VM):
## Exposure-stratified case-cohort analysis, II.Borgan method.
##                  10      11   12  2   3       4       5        6        7       8       9
## subcohort   107.000   32.00   11 20 120 216.000 322.000  538.000  584.000  260.00  190.00
## cohort     4203.917 2789.75 1189 20 120 372.695 886.543 1986.818 3096.096 4017.78 4728.21
##
## Coefficients:
##                                  Coef     HR  (95%    CI)     p
## age_at_baseline                 0.108  1.114  1.099  1.128 0.000
## Smoke12Current <median packyrs  2.291  9.888  6.163 15.865 0.000
## Smoke12Current >median packyrs  3.014 20.378 12.853 32.309 0.000
## Smoke12Former & duration>15     1.695  5.449  3.402  8.726 0.000
## first_year                     -0.019  0.981  0.961  1.003 0.085
## CSD_CUM_quantileQ1              0.098  1.103  0.774  1.573 0.587
## CSD_CUM_quantileQ2             -0.140  0.869  0.602  1.256 0.456
## CSD_CUM_quantileQ3              0.032  1.033  0.723  1.476 0.859
## CSD_CUM_quantileQ4              0.394  1.482  1.049  2.095 0.026
## Asb_CumDose                    -0.002  0.998  0.991  1.005 0.554

## Borgan II standard errors (total)
sqrt(diag(borgan2_model$var))
## Output (run on the VM):
## age_at_baseline                0.006583056
## Smoke12Current <median packyrs 0.241197357
## Smoke12Current >median packyrs 0.235146500
## Smoke12Former & duration>15    0.240254086
## first_year                     0.010913471
## CSD_CUM_quantileQ1             0.181065840
## CSD_CUM_quantileQ2             0.187614956
## CSD_CUM_quantileQ3             0.182259961
## CSD_CUM_quantileQ4             0.176532775
## Asb_CumDose                    0.003425894

## Square root of the sampling (phase two) component, same order as above
sqrt(diag(borgan2_model$phase2var))
## Output (run on the VM):
## [1] 0.002900459 0.066331965 0.069226270 0.055989325 0.005764428
## [6] 0.097204200 0.099605244 0.101309083 0.102139356 0.001777189
##-----------------------------------------------------------------------------
