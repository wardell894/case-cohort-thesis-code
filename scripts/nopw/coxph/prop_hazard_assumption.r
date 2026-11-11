##=============================================================================
## Project: Master's thesis - variance estimation in the NOPW case-cohort study
## Script Purpose: Check the proportional hazards assumption of the main Cox
##                 model (Chapter 4) with scaled Schoenfeld residuals (cox.zph).
## Date: 25.09.2026
## Author: Steinar Jonathan Mpagi
##=============================================================================
## Load datasets and libraries
getwd()
setwd("<DATA PATH>")
load("data/analysis/full_cohort_and_subcohort.RData")
library(survival)
library(dplyr)
Lung_data = full_cohort_and_subcohort %>% filter(Case == 1 | Subcohort_indicator == 1)
##-----------------------------------------------------------------------------
# Set design weights to match from the new weight schema
w = Lung_data$weight
Lung_data$weight =  ifelse(w == 1.838, 1.847,
                            ifelse(w == 4.351, 4.346, Lung_data$weight))
##-----------------------------------------------------------------------------
## The reference level for smoking is set to the lowest  exposure category (optional)
Lung_data$Smoke12 = as.factor(Lung_data$Smoke12)
Lung_data = within(Lung_data, Smoke12 <- relevel(Smoke12, ref = "Never/former, dur<15"))
##-----------------------------------------------------------------------------
## Run the coxph model
coxph_model = coxph(Surv(time_to_event, Case) ~
                     age_at_baseline +
                     Smoke12 +
                     first_year +
                     CSD_CUM_quantile +
                     Asb_CumDose,
                     data = Lung_data,
                     weights = weight,
                     robust = TRUE)
##-----------------------------------------------------------------------------
## Check the proportional hazards assumption using scaled Schoenfeld residuals
ph_test = cox.zph(coxph_model)
ph_test
## Output (run on the VM; default transform = "km"):
##                  chisq df    p
## age_at_baseline  1.127  1 0.29
## Smoke12          2.573  3 0.46
## first_year       0.165  1 0.68
## CSD_CUM_quantile 1.130  4 0.89
## Asb_CumDose      0.110  1 0.74
## GLOBAL           6.658 10 0.76

## Time on study (identity), the version reported in the thesis (tab:ph-test)
ph_test = cox.zph(coxph_model, transform = "identity")
ph_test
## Output (run on the VM):
##                   chisq df    p
## age_at_baseline  1.2320  1 0.27
## Smoke12          2.0106  3 0.57
## first_year       0.0257  1 0.87
## CSD_CUM_quantile 1.3504  4 0.85
## Asb_CumDose      0.4137  1 0.52
## GLOBAL           6.0949 10 0.81
##-----------------------------------------------------------------------------