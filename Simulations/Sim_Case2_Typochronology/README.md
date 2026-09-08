# Case 2 — typochronological windows

The method, the model and the comparison models are documented once in
`Simulations/shared/README.md`. This file covers only what is specific to typochronology.

## What makes this case different

Windows of the "first half of the 1st c. CE" kind: **each sample carries its own independent
window**, with no shared phase structure — that is Case 3. The weight row is flat across the
window.

The error here is Berkson, not measurement-derived: the window is known to contain the true
date rather than being derived from a noisy measurement of it. So the midpoint model is
expected to be roughly unbiased on the slope, and the interesting result is calibration —
inflated `sigma` and intervals that are too narrow. That makes this case the reference the
others are read against, and the reason it is built first.

Note that "the true date is uniform within its window" is an **assumption**, not something
measured. Sub-case 2b breaks it deliberately.

## Generator

- Widths are drawn as a **share of the studied period** (5–25%), not in absolute years, so
  the same simulation covers a Roman sequence and a prehistoric one. This is Case 2's
  fine↔coarse factor, the counterpart of `H` in the phase-based cases.
- The **true date is drawn first** and the window placed around it. Drawing the window first
  makes wide windows cluster mid-timeline and leaves the true dates triangular rather than
  uniform.
- Position within the window is `Beta(skew_shape, 1)`: 1 is uniform (2a), above 1 leans late
  (2b).
- Window boundaries are multiples of 25 yr, as real datings are.

Windows may overhang the study period. A typological attribution does not know where the
study window stops, and clipping them creates degenerate one-year windows at the edges.

## Checks

```
Rscript Simulations/Sim_Case2_Typochronology/scripts/00_check.R
```

Last run, midpoint slope ratio by dating resolution — the Berkson prediction is that these
stay near 1 even as windows widen:

```
fine    mean width 103 yr   midpoint slope ratio 0.979
even    mean width 125 yr   midpoint slope ratio 0.972
coarse  mean width 147 yr   midpoint slope ratio 0.963
```

The marginal-vs-latent equivalence check in `shared/scripts/00_check_marginal.R` also runs
on Case 2 data, since flat rows are where the two must agree exactly.

## Running the study

```
Rscript Simulations/Sim_Case2_Typochronology/scripts/01_design.R
Rscript Simulations/Sim_Case2_Typochronology/scripts/00_check_median_identity.R
Rscript Simulations/Sim_Case2_Typochronology/scripts/03_recovery_study.R
Rscript Simulations/Sim_Case2_Typochronology/scripts/04_recovery_plots.R
```

2800 fits: 1400 datasets x 2 models.

Cheapest case: no calibration, narrow supports. Develop here and scale up to Case 1.

`00_check_median_identity.R` asserts that median and midpoint are the same estimator for
flat windows, which is why `03_` does not fit median at all here. The figures show one row,
labelled "Midpoint / Median".
