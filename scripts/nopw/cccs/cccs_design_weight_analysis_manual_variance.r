##=============================================================================
## Project: Master's thesis - variance estimation in the NOPW case-cohort study
## Script Purpose: Fit the main Cox model for lung cancer (Chapter 4) with the
##                 predefined design weights, step by step with the functions
##                 of CaseCohortCoxSurvival, and compute the EG variance split
##                 into the cohort and sampling components (the design columns
##                 of tab:calibration-gain and tab:calibration-decomp).
##                 Same data set-up as cccs_calibrated_weight_analysis_manual_variance.r,
##                 with design weights in coxph and calibrated = FALSE in
##                 influences.RH and variance.
## Date: 14.10.2026
## Author: Steinar Jonathan Mpagi
##=============================================================================
#load relevant packages
library(CaseCohortCoxSurvival)
library(dplyr)
library(haven)
library(nnet)
################################################################################
# 0 - load data
getwd()
setwd("<DATA PATH>")
load("data/case_cohort_cox_survival_data.rda")

# same weights as in new weight scheme
w = full_cohort_and_subcohort$weight
full_cohort_and_subcohort$weight = ifelse(w == 1.838, 1.847,
                                          ifelse(w == 4.351, 4.346, full_cohort_and_subcohort$weight))
################################################################################
# 1 - data set up

data_generator = function(dataset)
{
  #klargjøring for casecohortcoxsurvival metoder
  dataset$W = dataset$bcohort
  dataset$status = dataset$Case
  dataset$subcohort = dataset$Subcohort_indicator
  sn = table(dataset$bcohort)
  sm = table(dataset$bcohort[dataset$Subcohort_indicator == 1])

  dataset$strata.n = as.vector(sn[as.character(dataset$bcohort)])
  dataset$strata.m = as.vector(sm[as.character(dataset$bcohort)])

  dataset$Smoke12 = as.factor(dataset$Smoke12)
  dataset = within(dataset, Smoke12 <- relevel(Smoke12, ref = "Never/former, dur<15"))

  # NB: levels() renames by position, and after relevel the order is
  # Never/former, Current <median, Current >median, Former >15. The labels are
  # therefore shifted one step: the reference keeps the first position, so the
  # coefficients are right but carry the wrong names:
  #   Smoke12_driverCurrent_greater_than_median_packyrs  = current, < median pack-years
  #   Smoke12_driverFormer_and_duration_greater_than_15  = current, > median pack-years
  #   Smoke12_driverNever_former_duration_less_than_15   = former, duration > 15 years
  dataset$Smoke12_driver = dataset$Smoke12
  levels(dataset$Smoke12_driver) = c("Current_less_than_median_packyrs",
                                     "Current_greater_than_median_packyrs",
                                     "Former_and_duration_greater_than_15",
                                     "Never_former_duration_less_than_15")
  # CSD as factor
  dataset$CSD_CUM_quantile = as_factor(dataset$CSD_CUM_quantile)
  dataset$CSD_CUM_quantile = relevel(dataset$CSD_CUM_quantile, "Q0")

  # Set factor levels for stilling_20 + add missing + create indicator for phase2
  dataset$phase2 = ifelse(dataset$Subcohort_indicator == 1 | dataset$Case == 1, 1, 0)
  dataset = dataset |> mutate(stilling_20 = as_factor(stilling_20))
  dataset$stilling_20 |> as.factor()
  dataset = dataset |> mutate(stilling_20 = if_else(is.na(stilling_20), "Missing", stilling_20))

  # strip the different titel to the core to reduce the amount of levels
  dataset$stilling_20_trimmed = trimws(
    sub("\\s*(feilstavelse|flere titler|leder|lokalisasjon|lokasjon|samme)\\s*$", "",
        dataset$stilling_20, ignore.case = TRUE))
  #set more proxy variables as factor
  dataset$typeselsk_1_sp23 = as_factor(dataset$typeselsk_1_sp23)
  dataset$losemdamp_sp2915 = as_factor(dataset$losemdamp_sp2915)
  dataset$boring_2_sp91 = as_factor(dataset$boring_2_sp91)
  dataset$never98_2c = as_factor(dataset$never98_2c)
  dataset$royk60_sp37b6 = as_factor(dataset$royk60_sp37b6)
  dataset$yrke_idag_sp5 = as_factor(dataset$yrke_idag_sp5)
  #dataset = within(dataset, Smoke12 <- relevel(Smoke12, ref = "D"))

  dataset$weight = ifelse(is.na(dataset$weight), 1, dataset$weight)
  dataset$weights = dataset$weight
  dataset = dataset %>%
    mutate(phase2 = as.integer(Subcohort_indicator == 1 | Case == 1)) %>%
    mutate(phase2_no_missing = as.integer(phase2 == 1 & !is.na(Smoke12))) %>%
    mutate(Case_driver = ifelse(phase2_no_missing == 1, Case, 0)) %>%
    mutate(Subcohort_indicator_driver = ifelse(phase2_no_missing == 1, Subcohort_indicator, 0)) %>%
    mutate(full_cohort_no_missing = as.integer(!(phase2 == 1 & is.na(Smoke12))))

  dataset
}
nopw_cohort = data_generator(full_cohort_and_subcohort) %>% filter(!(phase2 == 1 & is.na(Smoke12)))
nopw_cohort = as.data.frame(nopw_cohort)

################################################################################
# CASE-COHORT SAMPLE
# (the imputation, auxiliary variables and calibration of the calibrated
# version are not needed for design weights and are left out)
casecohort = nopw_cohort[nopw_cohort$phase2 == 1, ]
################################################################################
# RUN COX REGRESSION ON CALIBRATED WEIGHTS
# change to design weights here for regular var
coxph_model = coxph(Surv(time_to_event, Case)
                    ~ age_at_baseline +
                      first_year +
                      Smoke12_driver +
                      CSD_CUM_quantile +
                      Asb_CumDose,
                    data = casecohort,
                    weight = casecohort$weight,
                    robust = TRUE)



################################################################################
# GET INFLUENCES
# calibrated = FALSE for design weights (no auxiliary matrix A needed)
influences_calibrated_breslow = influences.RH(coxph_model, calibrated = FALSE)
################################################################################

# set calibrated to false for design weights (infl2 is then not needed)
var_calib = variance(n = nrow(nopw_cohort),
                     casecohort = casecohort,
                     cohort = nopw_cohort,
                     infl = influences_calibrated_breslow$infl.beta,
                     stratified = TRUE,
                     calibrated = FALSE,
                     variance.phase2 = TRUE)

influences_calibrated_breslow$beta.hat
## Output (run on the VM, 14.10.2026):
##                                     age_at_baseline                                          first_year
##                                         0.107960566                                        -0.018897137
## Smoke12_driverCurrent_greater_than_median_packyrs Smoke12_driverFormer_and_duration_greater_than_15
##                                         2.292083937                                         3.018109517
##  Smoke12_driverNever_former_duration_less_than_15                                  CSD_CUM_quantileQ1
##                                         1.693418953                                         0.096629751
##                                  CSD_CUM_quantileQ2                                  CSD_CUM_quantileQ4
##                                        -0.141651027                                         0.388882813
##                                  CSD_CUM_quantileQ3                                         Asb_CumDose
##                                         0.030465809                                        -0.002052561

var_calib
## Output (run on the VM, 14.10.2026), unnamed, in the same order as beta.hat:
## $variance
##  [1] 4.152505e-05 1.327225e-04 5.779021e-02 5.494445e-02 5.897214e-02 3.314843e-02 3.544859e-02 3.089485e-02 3.357858e-02 1.145308e-05
##
## $variance.phase2
##  [1] 8.508334e-06 3.318171e-05 4.395544e-03 4.782593e-03 3.127875e-03 9.444966e-03 9.914394e-03 1.035448e-02 1.025616e-02 3.149755e-06
##
## Matches the design columns of tab:calibration-gain (beta and SE = sqrt(variance))
## and tab:calibration-decomp (sampling = variance.phase2, cohort = variance - variance.phase2).
################################################################################
