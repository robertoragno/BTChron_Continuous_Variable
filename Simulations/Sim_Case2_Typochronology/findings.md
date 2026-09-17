# Case 2 — finds with mixed dating precision

The method, the model and the comparison models are documented once in
`Simulations/shared/README.md`. This file covers only what is specific to Case 2.
The findings from the previous design are in
`archive/superseded_run_20260916/findings.md`.

## What makes this case different

Each find carries its own window, with no shared phases (those are Cases 3 and
4). Some windows are short (well dated), some long (coarsely dated). The weight
row is flat across the window.

## What to expect, before the full run

The true date is drawn first and the window placed around it. The midpoint then
misses the true date by an amount that does not depend on the date itself, which
flattens the midpoint slope, more so the longer the windows are compared with
the period. Measured with plain least squares on 200,000 finds (no Stan model):

```
coarse window   coarsely dated   midpoint slope ratio
10% of period   50%              0.993
10% of period   100%             0.990
25% of period   50%              0.968
25% of period   100%             0.940
45% of period   50%              0.906
45% of period   100%             0.830
```

These follow `1 / (1 + average(width^2) / period^2)`.

The full-distribution model should do better but is not expected to be exact.
It spreads each find evenly across its window, including any part that runs
past the period, because it is not told when the study period starts and ends
(the no-window-prior decision in `shared/scripts/weight_rows.R`).

Full run started 2026-09-16 (6400 fits), log in
`Simulations/logs/case2_refit_20260916.log`.

## Checks

- `scripts/01_design.R` stops if a window fails to contain its true date, if the
  realised proportion of coarsely dated finds is off, or if deposition growth and the
  precision trend do not show up in a large simulated dataset.
- `scripts/checks/00_check.R` draws the same checks as a figure.
- `scripts/checks/00_check_median_identity.R`: median and midpoint are the same
  estimator for a flat window, so median is not fitted. Passed on the new design.
