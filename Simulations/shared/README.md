# shared: the models and code all four cases use

Each dating case lives in its own `Sim_Case*/` folder. The models, the weight rows,
the fitting loop and the summary figures live here. The overall design is in
`docs/IMPLEMENTATION_PLAN.md`.

## One idea behind all four cases

Every case gives the model a grid of candidate calendar years and, for each find, a
row of probabilities over that grid. Only the shape of the row changes:

| Case | Dating | Weight row |
|---|---|---|
| 1 | radiocarbon | the calibrated distribution |
| 2 | typochronology | flat over the find's own window |
| 3 | overlapping phases | flat over the assigned phase |
| 4 | merged phases | flat over the run of merged phases |

So one Stan model serves all four cases.

## The models

| Name in figures | In code | What the date becomes | Stan file |
|---|---|---|---|
| Midpoint | `midpoint` | the middle of the range | `midpoint.stan` |
| Calibrated median | `median` | the median of the date distribution | `midpoint.stan` with start = end = median |
| Full distribution | `marginal` | every candidate year, weighted by its probability | `marginal_date.stan` |

For a flat window the median equals the midpoint, so Cases 2-4 fit only the midpoint
and label it "Midpoint / Median" (checked in Case 2's
`checks/00_check_median_identity.R`).

The full-distribution model does not pick one year per find. For each candidate year
it asks how well the trend fits if the find dates to that year, and averages the
answers weighted by that year's probability (the `log_sum_exp` in the Stan code). The
date is integrated out, so the model has only three parameters (intercept, slope,
sigma) whatever the number of finds. This is an errors-in-variables regression with
a known dating distribution per find; sigma is the scatter around the trend, not the
dating error.

`latent_date.stan` samples one date per find inside its window instead. For flat
windows it gives the same answer as the full-distribution model, so it is not
reported, but `checks/00_check_marginal.R` uses it to check that the full-distribution
model is correct. `marginal_date_period.stan` is a prototype that also estimates the
study period, used only by Case 1's `diagnostics/`.

## Layout

```
shared/
  models/    the Stan files above
  scripts/   weight_rows.R          weight rows for all cases
             design.R               ids, seeds and true trends for every 01_
             fit_models.R           Stan data per model, one fit to one dataset
             run_recovery.R         fits every model to every dataset, in parallel
             recovery_summary.R     metrics, tables and figures for every 04_
             single_fit_figure.R    one dataset, two models, for every 05_
             case_dating_figures.R  figures/caseN_dating.png: how each case dates finds
  checks/    00_check_marginal.R    full-distribution vs latent model
```

Each case folder has the same structure: `simulate.R`, then `01_design.R` (writes
`data/design.csv`), `03_recovery_study.R` (fits), `04_recovery_plots.R` (tables and
figures), `05_single_fit.R` (one fit up close). Only Case 1 has a `02_figures.R`. The
`03_` and `04_` scripts read every setting from the design file.

## Running a study

```
Rscript Simulations/Sim_Case1_Radiocarbon/scripts/01_design.R
Rscript Simulations/Sim_Case1_Radiocarbon/scripts/03_recovery_study.R   # hours
Rscript Simulations/Sim_Case1_Radiocarbon/scripts/04_recovery_plots.R
```

`03_` does nothing if `output/recovery_results.csv` already exists; delete it to refit.
A failed fit is recorded in the `error` column instead of stopping the run. For a
quick test:

```
RECOVERY_LIMIT=20 RECOVERY_WORKERS=8 RECOVERY_WARMUP=300 RECOVERY_SAMPLING=300 \
  RECOVERY_OUT=/tmp/test.csv Rscript .../03_recovery_study.R
```

## Metrics

| Metric | Meaning | Target |
|---|---|---|
| bias | mean of posterior median minus the true value | 0 |
| RMSE | root mean squared error of the posterior median | as small as possible |
| empirical SE | SD of the errors across datasets | as small as possible |
| calibration slope | slope of estimated slopes regressed on true slopes | 1 |
| coverage | proportion of 90% equal-tailed credible intervals that contain the true value | 0.90 |
| interval width | mean width of the 90% credible interval | narrower, once coverage is on target |

The calibration slope is reported because bias cannot show attenuation: true slopes
are drawn around zero, so flattening positive and negative slopes cancels out in the
average error. Below 1 the trend is attenuated (flattened), above 1 exaggerated.
Case READMEs call it the "slope ratio".

Note on terms. Earlier versions of these scripts called coverage "accuracy" and
interval width "precision". Both names were dropped because they clash with standard
usage: accuracy is closeness of the estimates to the truth (summarised here by bias
and RMSE), and precision is their spread across datasets (the empirical SE). The
performance measures follow Morris, White and Crowther (2019, Statistics in Medicine
38: 2074-2102). Coverage and width describe the posterior intervals, not the point
estimates.

The error bars in the summary figures show how precisely each number is known from a
limited number of simulated datasets: a Jeffreys interval for coverage (it stays
inside 0-1), an OLS 95% confidence interval for the calibration slope, and +/- 2
Monte Carlo standard errors for bias and interval width. Each fit saves the 50, 80,
90 and 95% intervals, so the reported level can be changed without refitting.

## Checks

Run both after any change to `weight_rows.R`:

- `checks/00_check_marginal.R`: the full-distribution and latent models must agree on
  a Case 2 dataset. Last run: slope 0.01921 vs 0.01920.
- `Sim_Case1_Radiocarbon/scripts/checks/00_check_rows.R`: each calibrated weight row
  must have the same median as rcarbon's calibration.

`Sim_Case1_Radiocarbon/scripts/checks/00_check_fit.R` runs the whole radiocarbon path
on one dataset.

## Things that fail without an error

- **Calendar order of calibrated dates.** rcarbon labels rows in cal BP, oldest first,
  which is already ascending calendar order. The labels need converting
  (year = 1950 - BP) but the rows must not be reversed; reversing them flips the sign of
  the slope. `calmatrix_to_calendar()` checks this.
- **Dates cut at the grid edge.** A calibrated date that loses probability off the grid
  is pulled inward. `calibrated_rows()` stops if a date loses more than 1%, and Case 1
  pads the grid (588 yr) so it never happens.
- **Coarse grids.** On a 5-year grid each cell sums the yearly probabilities around it.
  Taking every fifth year would throw away most of a calibrated distribution.
- **Loading rcarbon draws random numbers.** Case 1's design only reproduces because
  rcarbon is first loaded after `set.seed()`.

## Calendar axis and grid

The calendar range passed to Stan (`time_ref_min`, `time_ref_range`) is fixed per case
in the design, not taken from each dataset, so the prior on the slope means the same
thing in every dataset and model.

The full-distribution model checks every candidate year of every row at every sampler
step, so its cost grows with the grid resolution:

| Grid step | Years per row | Per chain | vs latent model |
|---|---|---|---|
| 1 yr | ~120 | ~37 s | ~25x |
| 5 yr | 26 | ~8.5 s | ~7x |

The 5-year grid is the default. On Case 2 (previous design, windows up to 200 yr) it
gave the same slope as the annual grid (0.01967 vs 0.01970). It has not been checked
separately for Case 1, where multi-peaked calibrated dates could be smoothed by a
coarse grid. Set the step with `MARGINAL_GRID_STEP=5`.
