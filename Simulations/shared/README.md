# shared: the model all four cases use

## One model

`models/linear_dates.stan` fits `value = intercept + slope * date + noise` when the
date of each find is uncertain. Each find comes with a set of candidate years and a
probability for each. The model averages the fit over those years, weighted by
their probability, so the date never has to be fixed to one year.

Every method, in every case, is the same model with a different table of dates:

| Method | Candidate years |
|---|---|
| Midpoint | one year, the middle of the range, probability 1 |
| Calibrated median (Case 1) | one year, the median, probability 1 |
| Full distribution, Case 1 | every grid year, with its calibrated probability |
| Full distribution, Cases 2-4 | every grid year inside the window, equally likely |

For a flat window the median is the midpoint, so Cases 2-4 fit it once and call it
"Midpoint / Median".

Details:

- The grid is every 5 years. In Case 1 each 5-year cell adds up the yearly
  calibrated probabilities around it.
- Each find only uses the run of grid years where it has probability (`first` to
  `last`), which keeps the fits fast.
- The slope is sampled per 100 years and reported per year. Per year it is a very
  small number next to the intercept, which makes sampling about 4 times slower.
- Priors: slope normal(0, 0.05) per year, the same in every case; value at the
  centre of the period normal(0, 10); sigma exponential(1).

## Metrics

From `04_recovery_plots.R` in each case, written to `output/recovery_metrics.csv`:

| Metric | Meaning | Target |
|---|---|---|
| slope bias | mean of posterior median minus true slope | 0 |
| calibration slope | estimated slopes regressed on true slopes | 1 |
| 90% coverage | share of 90% intervals that contain the true slope | 0.90 |
| 90% interval width | mean width of those intervals | narrower, once coverage is right |
| bias in sigma | mean of posterior median minus true sigma | 0 |

True slopes are drawn on both sides of zero, so a method that flattens every slope
can still have zero bias. The calibration slope shows it: below 1 the trend is
flattened, above 1 exaggerated. Case READMEs call it the "slope ratio". Measures
follow Morris, White and Crowther (2019, Statistics in Medicine 38: 2074-2102).

Error bars are 95% intervals for how well each number is known from 100 or so
datasets: +/- 2 standard errors for means, the `lm` confidence interval for the
calibration slope, a Jeffreys interval for coverage.

## Other files

`scripts/` and the other Stan files are the code before the rewrite of October
2026. They are still used by the checks, by Case 1's `diagnostics/`, and by the
`05_single_fit.R` of Cases 2-4. `scripts/case_dating_figures.R` draws
`figures/caseN_dating.png` for the main README.

## Things that fail without an error

- rcarbon labels calibrated years in cal BP, oldest first. That is already
  ascending calendar order: convert the labels (year = 1950 - BP), do not reverse
  the rows, or the slope changes sign.
- A calibrated date cut at the grid edge is pulled inward. Case 1 pads the grid by
  588 years on each side so this never happens.
- Loading rcarbon draws random numbers. Case 1's design only reproduces because
  rcarbon is loaded after `set.seed()`.
