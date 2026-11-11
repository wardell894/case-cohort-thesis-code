##=============================================================================
## Project: Master's thesis - variance estimation in the NOPW case-cohort study
## Script Purpose: Fit the main Cox model for lung cancer (Chapter 4) with
##                 calibrated weights, step by step with the functions of
##                 CaseCohortCoxSurvival: impute the phase-two covariates from
##                 the proxies, build the auxiliary variables from the imputed
##                 model, calibrate the design weights, and compute the EG
##                 variance split into the cohort and sampling components
##                 (tab:calibration-gain and tab:calibration-decomp).
## Date: 27.09.2026
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
#RUN IMPUTATION
for (v in c("typeselsk_1_sp23", "stilling_20_trimmed", "bcohort", "boring_2_sp91", "royk30_sp37b3", "losemdamp_sp2915")) nopw_cohort[[v]] = as_factor(nopw_cohort[[v]])


imputation_age_at_baseline = lm(age_at_baseline ~ fodselsaar_fr + typeselsk_1_sp23, data = nopw_cohort, weights = weight)
imputation_first_year = lm(first_year ~ fodselsaar_fr + bcohort, data = nopw_cohort, weights = weight)
imputation_Asb_CumDose = lm(Asb_CumDose ~ stilling_20_trimmed, data = nopw_cohort, weights = weight)


imputation_CSD_CUM_quantile = multinom(CSD_CUM_quantile ~ stilling_20_trimmed +
                                         losemdamp_sp2915 +
                                         bcohort +
                                         boring_2_sp91,
                                       data = nopw_cohort,
                                       maxit = 1000,
                                       weights = weight)
imputation_Smoke12 = multinom(Smoke12 ~ royk30_sp37b3, data = nopw_cohort,
                              weights = weight)

nopw_cohort = nopw_cohort %>%
  mutate(age_at_baseline_pred = predict(imputation_age_at_baseline, newdata = nopw_cohort)) %>%
  mutate(first_year_pred = predict(imputation_first_year, newdata = nopw_cohort)) %>%
  mutate(Asb_CumDose_pred = predict(imputation_Asb_CumDose, newdata = nopw_cohort)) %>%
  mutate(Smoke12_pred = predict(imputation_Smoke12, newdata = nopw_cohort)) %>%
  mutate(CSD_CUM_quantile_pred = predict(imputation_CSD_CUM_quantile, newdata = nopw_cohort))

nopw_cohort = nopw_cohort %>%
  mutate(age_at_baseline_pred_centered = age_at_baseline_pred - mean(age_at_baseline_pred)) %>%
  mutate(first_year_pred_centered = first_year_pred - mean(first_year_pred))
################################################################################
# RUN COX REGRESSION ON IMPUTED COVARIATES

mod_imputation = coxph(Surv(time_to_event, Case)
                       ~ age_at_baseline_pred +
                         Smoke12_pred +
                         CSD_CUM_quantile_pred +
                         Asb_CumDose_pred +
                         first_year_pred,
                       data = nopw_cohort, robust = T)

summary(mod_imputation)
################################################################################
# CREATE AUXILIARY VARIABLES FROM IMPUTATION MODEL

casecohort = nopw_cohort[nopw_cohort$phase2 == 1, ]
w = casecohort$weight
aux = auxiliary.construction(mod_imputation, Tau1 = 0,
                             Tau2 = max(nopw_cohort$time_to_event),
                             casecohort = casecohort,
                             method = "Breslow",
                             time.on.study = nopw_cohort$time_to_event)

Aux_matrix = as.matrix(cbind(1, aux$A.RH.Breslow))
################################################################################
# WEIGHT CALIBRATION

ph2 = which(nopw_cohort$phase2 == 1)
calibration = calibration(A.phase2 = Aux_matrix[ph2, ],
                          design.weights = w,
                          total = colSums(Aux_matrix),
                          eta0 = rep(0, ncol(Aux_matrix)))
################################################################################
# RUN COX REGRESSION ON CALIBRATED WEIGHTS
coxph_model = coxph(Surv(time_to_event, Case)
                    ~ age_at_baseline +
                      first_year +
                      Smoke12_driver +
                      CSD_CUM_quantile +
                      Asb_CumDose,
                    data = casecohort,
                    weight = c(calibration$calibrated.weights),
                    robust = TRUE)



################################################################################
# GET INFLUENCES

influences_calibrated_breslow = influences.RH(coxph_model, calibrated = TRUE, A = Aux_matrix)
################################################################################

var_calib = variance(n = nrow(nopw_cohort),
                     casecohort = casecohort,
                     cohort = nopw_cohort,
                     infl = influences_calibrated_breslow$infl.beta,
                     infl2 = influences_calibrated_breslow$infl2.beta,
                     stratified = TRUE,
                     calibrated = TRUE,
                     variance.phase2 = TRUE)

influences_calibrated_breslow$beta.hat
## Output (run on the VM, 14.10.2026):
##                                     age_at_baseline                                          first_year
##                                        0.1031658053                                       -0.0178424897
## Smoke12_driverCurrent_greater_than_median_packyrs Smoke12_driverFormer_and_duration_greater_than_15
##                                        2.3711696749                                        3.0250134433
##  Smoke12_driverNever_former_duration_less_than_15                                  CSD_CUM_quantileQ1
##                                        1.7298622837                                        0.1073138242
##                                  CSD_CUM_quantileQ2                                  CSD_CUM_quantileQ4
##                                       -0.1361181612                                        0.3768312045
##                                  CSD_CUM_quantileQ3                                         Asb_CumDose
##                                        0.0466365007                                       -0.0007080985

var_calib
## Output (run on the VM, 14.10.2026):
## $variance
##                                     age_at_baseline                                          first_year
##                                        3.869082e-05                                        1.359883e-04
## Smoke12_driverCurrent_greater_than_median_packyrs Smoke12_driverFormer_and_duration_greater_than_15
##                                        5.597598e-02                                        5.353993e-02
##  Smoke12_driverNever_former_duration_less_than_15                                  CSD_CUM_quantileQ1
##                                        5.876964e-02                                        3.307591e-02
##                                  CSD_CUM_quantileQ2                                  CSD_CUM_quantileQ4
##                                        3.460887e-02                                        2.864168e-02
##                                  CSD_CUM_quantileQ3                                         Asb_CumDose
##                                        3.178185e-02                                        9.900575e-06
##
## $variance.phase2
##                                     age_at_baseline                                          first_year
##                                        6.264637e-06                                        3.400829e-05
## Smoke12_driverCurrent_greater_than_median_packyrs Smoke12_driverFormer_and_duration_greater_than_15
##                                        2.537395e-03                                        2.693605e-03
##  Smoke12_driverNever_former_duration_less_than_15                                  CSD_CUM_quantileQ1
##                                        2.245759e-03                                        9.304600e-03
##                                  CSD_CUM_quantileQ2                                  CSD_CUM_quantileQ4
##                                        9.196832e-03                                        8.518765e-03
##                                  CSD_CUM_quantileQ3                                         Asb_CumDose
##                                        9.343355e-03                                        2.322498e-06
##
## Matches the calibrated columns of tab:calibration-gain (beta and SE = sqrt(variance))
## and tab:calibration-decomp (sampling = variance.phase2, cohort = variance - variance.phase2).
################################################################################
