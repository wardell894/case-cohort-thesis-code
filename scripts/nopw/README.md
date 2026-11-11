# NOPW analysis scripts (Chapter 4)

These scripts produce the results for the NOPW case-cohort study in Chapter 4. They were run in R 4.5.1 on the secure server at the Cancer Registry of Norway, using `CaseCohortCoxSurvival` 0.0.36. They cannot be run outside that server.

- The data path is replaced by `setwd("<DATA PATH>")`.
- The output from each run is pasted into the script as comments. Only aggregate results are included.

## What is not included

- **Data preparation.** Every script starts from a finished analysis file: `case_cohort_cox_survival_data.rda`, or `full_cohort_and_subcohort.RData` for the `coxph` and `cch` scripts. These files were built on the server from the registry data, and that step is not included here.
- **Table 4.2 (`tab:nopw-design`).** The counts by birth cohort come from a spreadsheet, not from a script.
- **The EG columns of `tab:variance-decomp` and the numbers behind `fig:nopw-factor`.** They come from an earlier run. They differ by less than 1% from `cccs/cccs_design_weight_analysis_manual_variance.r`, as the table caption notes.

## Scripts

| Script | Gives |
|---|---|
| `coxph/coxph_analysis.r` | Main Cox model with design weights and robust variance (`tab:baseline-coefs`) |
| `coxph/prop_hazard_assumption.r` | Test of proportional hazards with `cox.zph` (`tab:ph-test`) |
| `cch/borgan_2_analysis_cch.r` | Borgan II variance with `cch`, using adjusted stratum sizes so that the design weights are reproduced (`tab:variance-comparison`, the Borgan II columns of `tab:variance-decomp`) |
| `cccs/cccs_design_weight_analysis.r` | EG and robust variance from the package driver with design weights (`tab:beta-comparison`, `tab:variance-comparison`) |
| `cccs/cccs_design_weight_analysis_manual_variance.r` | EG variance with design weights, computed step by step and split into cohort and sampling components (the design columns of `tab:calibration-gain` and `tab:calibration-decomp`) |
| `calibration/cccs_calibrated_weight_analysis_manual_variance.r` | Imputation, calibration and EG variance step by step, split into cohort and sampling components (the calibrated columns of `tab:calibration-gain` and `tab:calibration-decomp`) |
| `calibration/cccs_calibrated_weight_analysis.r` | The same model with the package driver and `calibrated = TRUE`, as a check (caption of `tab:calibration-gain`) |
| `pure_risk/pure_risk_driver.r` | Pure risk for the covariate profiles, with EG and robust variance (`tab:pure-risk`, `fig:nopw-pure-risk`, Section 4.5) |
| `div/r_squared_bcohort.r` | Weighted R² of each covariate on birth cohort (`tab:strata-r2`) |
| `div/proxy_selection_xgboost.r` | Ranking of candidate proxies with xgboost (`tab:proxy-models`). The list of candidate variables is omitted. |

## Pitfalls in `CaseCohortCoxSurvival` found along the way

- **Tibbles.** Pass `data` as a `data.frame`. With a tibble the calibrated variance is wrong, and pure risk stops with `INTERNAL CODING ERROR 1`.
- **Factor levels in pure risk.** The driver codes the levels in the profiles (`x`) alphabetically. The reference level must therefore come first in alphabetical order, which is why the smoking reference is named `"0"`. When several profiles are run together, the reference level must be among them.
- **Smoking labels in the step-by-step scripts.** In those scripts `levels()` renames `Smoke12_driver` by position, so the smoking labels are shifted one step. The coefficients are right, and the scripts note the correct meaning of each label.
