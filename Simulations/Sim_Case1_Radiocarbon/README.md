# Case 1 — Radiocarbon

## The dating problem

## How the data is generated

## What the models see

## Result

## Figures

| File | Shows |
|---|---|
| `figures/dataset_anatomy.png` | one dataset per window (Hallstatt plateau vs steep section), before any fitting |
| `figures/single_fit_comparison.png` | one dataset, midpoint vs full distribution, fitted up close |
| `figures/recovery_summary.png` | headline result across all datasets |
| `figures/recovery_by_lab_error.png` | how the result changes with lab error |
| `figures/checks/` | calibrated-date overview, sample-size sweep, deposition sweep |

## Files

Run in order: `01` → `03` → `04`, then `02` and `05` for the illustration figures.

| File | What it does |
|---|---|
| `scripts/simulate.R` | draws true calendar dates, turns each into a 14C determination with lab error |
| `scripts/01_design.R` | builds the list of datasets to simulate, `data/design.csv` |
| `scripts/02_figures.R` | `dataset_anatomy.png` and `checks/overview.png` (the only case with a `02_`) |
| `scripts/03_recovery_study.R` | calibrates each dataset and fits every model to it; slow |
| `scripts/04_recovery_plots.R` | turns the fits into tables (`output/`) and figures |
| `scripts/05_single_fit.R` | `single_fit_comparison.png` |
| `scripts/checks/00_check_rows.R` | calibrated weight rows reproduce `rcarbon`'s own medians |
| `scripts/checks/00_check_fit.R` | simulate → calibrate → fit, on one dataset |
| `scripts/diagnostics/` | side investigations into why the steep section attenuates (`05_`) and the estimated study period (`06_`, `07_`) |

Case 1 fits three models (midpoint, median, marginal). The example figure shows only midpoint
and marginal, the pair compared in every case; median differs from midpoint only for
radiocarbon.

Detailed findings so far: `findings.md`.
