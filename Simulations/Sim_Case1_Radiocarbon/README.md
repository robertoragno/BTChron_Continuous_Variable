# Case 1 — radiocarbon

The method, the model and the comparison models are documented once in
`Simulations/shared/README.md`. This file covers only what is specific to radiocarbon.

## What makes this case different

The date information is a **calibrated posterior**, not a window: peaked, often multimodal,
and derived from a noisy measurement rather than from a range known to contain the truth.
This is the only case where the point summaries disagree with each other on the slope, and
where the calibrated median and the envelope midpoint are different numbers (see the results
below): the calibrated median attenuates the slope, the envelope midpoint over-estimates it
on the plateau, and only the date-marginalised model recovers it.

Dates are calibrated outside Stan with `rcarbon::calibrate(..., calMatrix = TRUE)` and
handed to the model as one probability per candidate year.

## Calendar windows

Two windows, both 400 yr, measured on IntCal20 at ±25:

| Window | mean 95% span | median vs envelope midpoint |
|---|---|---|
| **800–400 BCE** — Hallstatt plateau | 223 yr | 45 yr mean, 98 max |
| **1600–1200 BCE** — control | 147 yr | 9 yr mean, 25 max |

The plateau is where the models should come apart; the control is where they should agree. A
single plateau date can have a 95% HPD in **four disjoint chunks**, which is exactly the
structure a single-interval model must flatten and the date-marginalised model keeps.

## Grid padding

`c14_grid_pad()` measures the widest calibrated support rather than guessing it. Taken as
the maximum over both windows and all three lab errors:

| | ±15 | ±30 | ±50 |
|---|---|---|---|
| plateau | 358 | 505 | **588** |
| control | 343 | 390 | 543 |

`01_design.R` freezes that single 588 yr value across every cell, so `time_ref_range` is
identical (1576 yr) in all 1600 datasets and only `time_ref_min` shifts between windows.
Per-cell padding would make the slope prior differ between cells and break the comparison it
exists to support.

Truncation is an error, not a warning — see the shared README. Verified not to fire even
with every true date pinned to a window edge.

## The prior inside the calibrated posterior

`rcarbon` calibrates against a uniform prior over calendar years, so the row it returns is
proportional to the likelihood of the determination. **That is the intended behaviour and
the prior wanted here**, confirmed with the supervisor. Note in the methods that the prior
is uniform over the whole calibration range, not over the study window: the model does not
use the window as evidence about where a date can fall, which matches what an analyst
actually knows. It is also why the grid is padded rather than clipped — clipping would
impose a window prior by the back door.

## Lab error and curve error

A determination carries two independent uncertainties, and both belong in the generator:

- **lab error** — how precisely the AMS measured this sample. The `lab_error` sweep, ±15 / ±30 / ±50.
- **curve error** — IntCal20 is itself an estimate, worth ±12–16 yr across these windows.
  `rcarbon::uncalibrate()` returns it per year as `ccError`.

`rcarbon::calibrate()` uses both: its kernel is `dnorm(age, mu, sqrt(lab_error^2 + curve_error^2))`.
`simulate_c14()` therefore draws with the same combined sd, so the model being fitted is the
model that produced the data.

Drawing with the lab error alone — what this case did until 2026-09-08 — generated at one sd
and inverted at a wider one. The gap is largest where the lab error is smallest, so it ran
along the sweep rather than across it:

| lab error | sd used to generate | sd `calibrate()` assumes | too wide by |
|---|---|---|---|
| 15 | 15 | 20.5 | 37% |
| 30 | 30 | 33.1 | 10% |
| 50 | 50 | 51.9 | 4% |

That made posteriors broader than the noise warranted, so coverage read better than it was,
and it confounded the lab-error sweep with how misspecified the generator was at each level.

Note this is a correctness fix, not an explanation of the control-window result.
`05_attenuation_diagnostic.R` had already tested folding curve error into the generator as a
candidate explanation for that attenuation and found it does not account for it.

**Remaining limitation.** The curve error is drawn independently per date, whereas there is
one IntCal20 residual and it displaces neighbouring dates together. Correlated curve error
would be a systematic error no per-date model can absorb, including the date-marginalised
one. Out of scope here, stated rather than fixed.

## Deposition

Finds are placed across the window by `rdeposition()`, parameterised by `growth_ratio`: how
many times denser deposition is at the end of the window than at the start. `1` is uniform,
`4` is fourfold denser at the late edge, drawn by inverting the CDF of an exponential
density.

This is not an SPD, so deposition intensity is not the signal — the estimand is a slope of a
measured value on date, and non-uniform density in the dates changes precision, not slope
bias. It is swept here because Case 1 is the one case where dating quality depends on
*calendar position*: growth oversamples the late edge of the window, changing the mix of
well- and badly-determined dates. The flat-window cases have no such interaction and keep
uniform deposition.

The `deposition` block runs `growth_ratio = 4` in both windows at N = 200 and the reference
lab error. Uniform is not re-simulated — the `factor` cells at ±30 are its reference. One
figure per window rather than a pooled one, since the question is whether growth bites on the
plateau and not off it.

## Checks

```
Rscript Simulations/Sim_Case1_Radiocarbon/scripts/00_check_rows.R   # indexing and row order
Rscript Simulations/Sim_Case1_Radiocarbon/scripts/00_check_fit.R    # simulate -> calibrate -> fit
```

Last run of the end-to-end check, 50 dates on the plateau, 5-year grid,
`beta ~ normal(0, 40)`:

```
slope 0.0201 (true 0.020)   sigma 1.99 (true 2.0)   Rhat 1.00   ~6 s/chain
dates correlate 0.850 with truth
OLS on calibrated medians: 0.01843
```

One dataset says nothing about bias — that is the recovery study's job, run by
`03_recovery_study.R`.

## Recovery study results

4800 fits (1600 datasets x 3 models), run 2026-09-08 on the corrected generator: the
curve error folded into the draw, deposition swept, slope prior `beta ~ normal(0, 40)`.
The previous run is kept in `output/superseded_uniform_nocurve/`, and the one before it,
under the tighter prior, in `output/superseded_prior_normal10/`.

The headline table holds the deposition sweep out, so it is the reference condition
(uniform deposition) and stays comparable with the earlier runs. Growth appears only in
its own figures.

Attenuation is the slope of `estimate ~ truth`: 1 is faithful, below 1 flattens, above 1
exaggerates. Sigma error is the signed bias in the noise term (target 0).

**Plateau window (800-400 BCE), 1100 fits per model:**

| model | attenuation | slope coverage (90%) | sigma error |
|---|---|---|---|
| Midpoint (envelope) | 1.11 | 0.83 | +0.28 |
| Median (calibrated) | 0.82 | 0.73 | +0.26 |
| Date-marginalised | 0.96 | 0.84 | -0.03 |

**Control window (1600-1200 BCE), 300 fits per model:**

| model | attenuation | slope coverage (90%) | sigma error |
|---|---|---|---|
| Midpoint (envelope) | 0.80 | 0.43 | +0.14 |
| Median (calibrated) | 0.82 | 0.47 | +0.14 |
| Date-marginalised | 0.83 | 0.49 | +0.01 |

The three separate on the plateau as intended, and the ordering is unchanged from the
earlier runs: the date-marginalised model recovers the slope and the noise term, the
envelope midpoint over-estimates the slope, the calibrated median under-estimates it.
Sigma is where the gap is widest and least ambiguous — -0.03 against +0.26 and +0.28.
Coverage is short of nominal 0.90 for every model, including the date-marginalised one.

### What the curve-error fix changed

Same cells, same prior; only the generator differs.

| | attenuation | coverage | sigma error |
|---|---|---|---|
| Date-marginalised, plateau | 0.983 → 0.961 | 0.859 → 0.842 | -0.048 → -0.029 |
| Median, plateau | 0.863 → 0.822 | 0.774 → 0.726 | +0.275 → +0.258 |
| Midpoint, plateau | 1.172 → 1.113 | 0.772 → 0.828 | +0.309 → +0.283 |
| all three, control | ~0.84 → ~0.82 | ~0.51 → ~0.46 | little change |

Coverage falls almost everywhere, which is the flattering being removed: posteriors are no
longer wider than the noise warrants. The exception is the envelope midpoint on the
plateau, whose coverage *improves* because its slope inflation drops from 17% to 11% —
over-wide posteriors were producing over-wide envelopes whose midpoints spread further
apart than the truth, exaggerating the slope. So the old generator was overstating the
midpoint model's plateau failure. The contrast is a little weaker than previously reported
and it is the honest one.

### Deposition: underpowered, nothing resolved

Uniform against `growth_ratio = 4`, 100 datasets per cell, both windows. **Every
comparison overlaps** — no effect is resolved at this replication, in either window, on
any model.

The one near-miss is the cell where an effect was expected, the date-marginalised model on
the plateau: attenuation 0.986 [0.958, 1.015] uniform against 0.936 [0.912, 0.961] under
growth, intervals overlapping by 0.003. Directionally growth pulls it away from 1, which is
what oversampling the late edge of a plateau should do, but at n = 100 that cannot be
claimed.

This is the same wall Case 2's skew cells hit, and the same fix applies: `N_REP` for the
deposition block should go to 300, as `N_REP_SKEW` did. That is +400 datasets, roughly
45 minutes. Reported as a null in the meantime, not as evidence of no effect.

### Convergence

4800 fits, none errored, **no divergent transitions at all**, worst Rhat 1.0141 with 34
fits above 1.01 and none above 1.05. The `deposition` block is the cleanest of the three
(worst Rhat 1.010), so the growth geometry costs nothing in sampling. `04_` reprints this
on every run.

`04_` also reports that the 95% HPD envelope retains 97.5% of the calibrated mass on
average, and a false-positive rate on the zero-slope cells of 0.115 / 0.130 / 0.140 for
midpoint / median / date-marginalised against a nominal 0.10.

## Running the study

```
Rscript Simulations/Sim_Case1_Radiocarbon/scripts/01_design.R          # once, ~1 min
Rscript Simulations/Sim_Case1_Radiocarbon/scripts/03_recovery_study.R  # 4800 fits, ~3 h on 24 workers
Rscript Simulations/Sim_Case1_Radiocarbon/scripts/04_recovery_plots.R
Rscript Simulations/Sim_Case1_Radiocarbon/scripts/05_attenuation_diagnostic.R  # optional, chases the control-window result
```

4800 fits: 1600 datasets x 3 models — core 800, factor 600, deposition 200. Calibration
happens once per dataset and is shared by all three, so they see identical dates and differ
only in what they do with them.

Case 1 is the expensive case — widest supports, and the plateau widest of all. Develop
against `RECOVERY_LIMIT` and scale up; see the shared README for the env vars.
