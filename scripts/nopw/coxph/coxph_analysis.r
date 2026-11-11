##=============================================================================
## Project: Master's thesis - variance estimation in the NOPW case-cohort study
## Script Purpose: Fit the main Cox model for lung cancer (Chapter 4) with the
##                 predefined design weights, and compute the robust variance.
## Date: 25.09.2026
## Author: Steinar Jonathan Mpagi
## Adapted from a script by [NAME] (NOPW group), with some changes and comments added. Aand
##=============================================================================
## Load datasets
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

summary(coxph_model)
## Output (run on the VM):
##   n= 2409, number of events= 463
##    (81 observations deleted due to missingness)
##
##                                   coef exp(coef)  se(coef) robust se      z Pr(>|z|)
## age_at_baseline               0.107961  1.114004  0.005886  0.006804 15.867  < 2e-16 ***
## Smoke12Current <median packyrs 2.292084 9.895538 0.231895  0.240627  9.525  < 2e-16 ***
## Smoke12Current >median packyrs 3.018110 20.452590 0.224678 0.234642 12.863  < 2e-16 ***
## Smoke12Former & duration>15   1.693419  5.438041  0.233651  0.242943  6.970 3.16e-12 ***
## first_year                   -0.018897  0.981280  0.009257  0.011561 -1.635   0.1021
## CSD_CUM_quantileQ1            0.096630  1.101452  0.152685  0.182024  0.531   0.5955
## CSD_CUM_quantileQ2           -0.141651  0.867924  0.158923  0.188191 -0.753   0.4516
## CSD_CUM_quantileQ3            0.030466  1.030935  0.151442  0.183225  0.166   0.8679
## CSD_CUM_quantileQ4            0.388883  1.475332  0.143797  0.175683  2.214   0.0269 *
## Asb_CumDose                  -0.002053  0.997950  0.002928  0.003384 -0.606   0.5442
##
##                               exp(coef) exp(-coef) lower .95 upper .95
## age_at_baseline                  1.1140    0.89766    1.0992     1.129
## Smoke12Current <median packyrs   9.8955    0.10106    6.1747    15.858
## Smoke12Current >median packyrs  20.4526    0.04889   12.9128    32.395
## Smoke12Former & duration>15      5.4380    0.18389    3.3779     8.755
## first_year                       0.9813    1.01908    0.9593     1.004
## CSD_CUM_quantileQ1               1.1015    0.90789    0.7710     1.574
## CSD_CUM_quantileQ2               0.8679    1.15217    0.6002     1.255
## CSD_CUM_quantileQ3               1.0309    0.96999    0.7199     1.476
## CSD_CUM_quantileQ4               1.4753    0.67781    1.0456     2.082
## Asb_CumDose                      0.9979    1.00205    0.9914     1.005
##
## Concordance= 0.853  (se = 0.01 )
## Likelihood ratio test= 904.5  on 10 df,   p=<2e-16
## Wald test            = 609.2  on 10 df,   p=<2e-16
## Score (logrank) test = 872.2  on 10 df,   p=<2e-16,   Robust = 243.4  p=<2e-16
##-----------------------------------------------------------------------------