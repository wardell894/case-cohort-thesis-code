# Code for my masters thesis

R scripts for the master's thesis in Data Science at the University of Oslo (2026) by Steinar Jonathan Mpagi, on variance estimation, calibration and pure risk in the case-cohort design, with the Norwegian Offshore Petroleum Workers (NOPW) cohort as the case study.

- `scripts/simulering/`: the simulation study in Chapter 3. Summary results in `results/`.
- `scripts/nopw/`: the NOPW analysis in Chapter 4. The data are not included; see `scripts/nopw/README.md`.
- `scripts/figurer.r`: the figures.

## Make the figures

From the repository root:

```sh
Rscript scripts/figurer.r
```

Needs `ggplot2`, `patchwork`, `scales`, `ggdag` and `dplyr`. Reads the finished simulation results in `results/` and writes the 12 figures of the thesis as PDF to `figures/` (a few seconds).

## Rerun the simulation (Chapter 3, optional)

Needs `CaseCohortCoxSurvival`, `survival` and `dplyr`. Each of `scripts/simulering/run_*.r` runs one simulation setting, takes hours, and writes its results to `results/` (`*_results.csv` and `*_runs.rds`).

## NOPW analysis (Chapter 4)

Cannot be rerun: the data are registry data on a secure server. See [`scripts/nopw/README.md`](scripts/nopw/README.md).
