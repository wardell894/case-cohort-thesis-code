##=============================================================================
## Project: Master's thesis - variance estimation in the NOPW case-cohort study
## Script Purpose: Rank the candidate phase-one variables (proxies) by gradient
##                 boosting (xgboost), separately for each phase-two covariate,
##                 to choose the predictors used in the imputation models for
##                 calibration (Section 4.4.1, tab:proxy-models).
## Date:
## Author: Steinar Jonathan Mpagi
##=============================================================================
## Load datasets and libraries
library("xgboost")
library("dplyr")
library(haven)
getwd()
setwd("<DATA PATH>")
load("data/case_cohort_cox_survival_data.rda")
##-----------------------------------------------------------------------------
# 1 - data set up: the case-cohort sample, where the phase-two covariates are known
obs = full_cohort_and_subcohort |>
  filter(Case == 1 | Subcohort_indicator == 1)

# phase-one candidates; the full list of columns is long and is left out here
phase1 = full_cohort_and_subcohort |>
  filter(Case == 1 | Subcohort_indicator == 1) |>
  select(...)
ncol(phase1)
##-----------------------------------------------------------------------------
# 2 - xgboost per phase-two covariate
# Shown for pack-years (the continuous measure behind Smoke12); the same steps
# are run with each phase-two covariate as Y.
# phase1 (the select above) holds the phase-one candidates together with the
# phase-two covariates packyrs.x, Asb_CumDose, age_at_baseline, CSD_CumDose and
# first_year.
Y = obs |> filter(!is.na(packyrs.x)) |> select(packyrs.x) |> pull()

X = phase1 |> filter(!is.na(packyrs.x))
nrow(X)
dtrain = xgb.DMatrix(data = X, label = Y)        # X numerisk, NA ok


bst = xgb.train(
  params  = list(objective = "reg:squarederror"),
  data    = dtrain,
  nrounds = 100
)

xgb.importance(model = bst) |> View()

full_cohort_and_subcohort |> filter(is.na(never98_2c)) |> nrow()
##-----------------------------------------------------------------------------
# 3 - linear models with the top-ranked predictors (final models in
#     cccs_calibrated_weight_analysis_manual_variance.r)
#asbest
#lm(Y ~ as_factor(stilling_20)+ as_factor(yrke_idag_sp5), data = phase1) |> summary()

#first_year
#lm(Y ~ fodselsaar_fr + as_factor(typeselsk_1_sp23) , data = phase1) |> summary()

#age_at_baseline
#lm(Y ~ fodselsaar_fr , data = phase1) |> summary()

#csd
#lm(Y ~ as_factor(stilling_20) + as_factor(losemdamp_sp2915) + as_factor(boring_2_sp91) + , data = phase1) |> summary()

#pacyrs -> smoke12

##-----------------------------------------------------------------------------
