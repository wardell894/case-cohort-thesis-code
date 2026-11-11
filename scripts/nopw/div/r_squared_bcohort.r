##=============================================================================
## Project: Master's thesis - variance estimation in the NOPW case-cohort study
## Script Purpose: Describe how informative the strata (birth cohort, 13 levels)
##                 are about each phase-two covariate (tab:strata-r2, Section 4.1.3).
##                 Each covariate is regressed on birth cohort in the case-cohort
##                 sample, weighted by the design weights so that R^2 refers to
##                 the cohort. Smoking status and the CSD quartile are represented
##                 by the continuous measures behind them, pack-years and
##                 cumulative CSD dose.
## Date:
## Author: Steinar Jonathan Mpagi
##=============================================================================
## Load datasets and libraries
library("dplyr")
library(haven)
getwd()
setwd("<DATA PATH>")
load("data/case_cohort_cox_survival_data.Rda")
##-----------------------------------------------------------------------------
# set weights to match new weights
w = full_cohort_and_subcohort$weight
full_cohort_and_subcohort$weight = ifelse(w == 1.838, 1.847,
                                          ifelse(w == 4.351, 4.346, full_cohort_and_subcohort$weight))
##-----------------------------------------------------------------------------
# 1 - data set up: the case-cohort sample (cases and subcohort)
subcohort_cases = full_cohort_and_subcohort |>  filter(Case == 1 | Subcohort_indicator == 1)
##-----------------------------------------------------------------------------
# 2 - weighted R^2 from regressing each covariate on birth cohort
#asbest
asbest = lm(Asb_CumDose ~ as_factor(bcohort), data = subcohort_cases, weights = weight) |> summary()
#first_year
first_year = lm(first_year ~ as_factor(bcohort), data = subcohort_cases, weights = weight) |> summary()
#age_at_baseline
age_at_baseline = lm(age_at_baseline ~ as_factor(bcohort), data = subcohort_cases, weights = weight) |> summary()
#csd
csd = lm(CSD_CumDose ~ as_factor(bcohort), data = subcohort_cases, weights = weight) |> summary()

#pacyrs -> smoke12
packyrs = lm(packyrs.x ~ as_factor(bcohort), data = subcohort_cases, weights = weight) |> summary()

print(c(asbest$r.squared, first_year$r.squared, age_at_baseline$r.squared, csd$r.squared, packyrs$r.squared))
## Output (run on the VM, 14.10.2026):
## [1] 0.08556769 0.46247023 0.97804257 0.01418465 0.15447867
## Matches tab:strata-r2 (age 0.978, first year 0.462, pack-years 0.154,
## asbestos 0.086, CSD dose 0.014; square roots 0.989, 0.680, 0.393, 0.293, 0.119).
##-----------------------------------------------------------------------------
