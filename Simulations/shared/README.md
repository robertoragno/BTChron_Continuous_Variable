# shared — the model all four cases run through

The paper's four dating cases each live in their own `Sim_Case*/` folder. What they share
lives here: one Stan model, one weight-row builder, one equivalence check. See
`docs/IMPLEMENTATION_PLAN.md` §3–§7 for the design; this file records how the pieces fit.

## The one thing that unifies the four cases

Every case reduces to the same object: a **shared grid of candidate calendar years**, and
**one probability row per sample** over that grid. What differs is the shape of the row.

| Case | Dating | Weight row |
|---|---|---|
| 1 | radiocarbon | the calibrated posterior, peaked or multimodal (`rcarbon` `calMatrix`) |
| 2 | typochronology | flat over the sample's own window |
| 3 | overlapping phases | flat over the assigned phase's window |
| 4 | merged phases | flat over the union of the merged phases |

Because the row is just data, one Stan file covers all four, and a 14C date and a "first
half of the 1st c. CE" attribution differ only in how peaked their row is.

## Layout

```
shared/
  models/   marginal_date.stan   latent_date.stan   midpoint.stan
  scripts/  weight_rows.R          rows and grid indexing
            partition.R            broken-stick periodisation for the phase cases
            fit_models.R           the models; dispatch on MODEL, never on case
            design.R               ids, seeds and nuisance parameters, for every 01_
            run_recovery.R         the study itself: jobs, workers, results file
            recovery_summary.R     metrics, headline plots, sweep_figure() and
                                   dataset_anatomy(), for every 04_ and 00_check
            00_check_marginal.R    marginal vs latent equivalence

Sim_Case1_Radiocarbon/    simulate.R  00_check_rows.R  00_check_fit.R
                          01_design.R  02_figures.R  03_recovery_study.R  04_recovery_plots.R
Sim_Case2_Typochronology/ simulate.R  00_check.R
                          01_design.R  03_recovery_study.R  04_recovery_plots.R
Sim_Case3_OverlappingPhases/   simulate.R  00_check.R  README.md
                          PAUSED - overlap does not bias the slope, see its README
Sim_Case4_MergedPhases/   simulate.R
                          01_design.R  03_recovery_study.R  04_recovery_plots.R
```

Within a case folder: bare names are libraries, `00_` are checks, `01_` builds the design,
`02_` figures, `03_` runs the study, `04_` summarises it. A case's `03_`/`04_` re-declare
nothing — every setting comes from its `data/design.csv`.

`simulate.R`, `01_`, `03_`, `04_` are the spine, and every case that runs has all four and
nothing else structural. The extra numbers are case-specific evidence, not pipeline: Case 1
carries `02_figures.R` (a paper figure) and `05_attenuation_diagnostic.R` (why the steep
window attenuates). `03_` and `04_` are thin by design — each one names the case's own
`prepare_dataset()` and its own sweeps, and the machinery is in `run_recovery.R` and
`recovery_summary.R`.

Every case writes `figures/dataset_anatomy.png` from `dataset_anatomy()` in
`recovery_summary.R`: one representative dataset per level of the case's own factor, drawn
the way the archived Sim_Linear figures drew them — each find at its measured value, its
dating window a faint segment, the true calendar date a dark-red dot, the point-date
estimate a grey square, the true trend dashed. For radiocarbon the segment is the 95% HPD
envelope; for the phase cases it is the window, with the phase boundaries as bands behind
the cloud. It shows what the model is handed before any fitting.

## Running a study

```
Rscript Simulations/Sim_Case1_Radiocarbon/scripts/01_design.R          # once
Rscript Simulations/Sim_Case1_Radiocarbon/scripts/03_recovery_study.R  # hours
Rscript Simulations/Sim_Case1_Radiocarbon/scripts/04_recovery_plots.R
```

`03_` refuses to overwrite an existing `output/recovery_results.csv`; delete it to re-run.
While developing, shrink the job with env vars:

```
RECOVERY_LIMIT=20 RECOVERY_WORKERS=8 RECOVERY_WARMUP=300 RECOVERY_SAMPLING=300 \
  RECOVERY_OUT=/tmp/smoke.csv Rscript .../03_recovery_study.R
```

Two invariants worth knowing, both checked automatically:

- **midpoint and median must be identical for flat windows.** `midpoint.stan` takes
  `(start + end) / 2`, so the median model passes `start = end = median`. One Stan file, two
  point-date models. Case 2's `04_` prints the gap; it should be exactly zero.
- **failed fits are captured, not thrown.** One pathological dataset must not kill a run of
  thousands. `read_recovery()` reports how many were excluded rather than silently
  returning fewer rows.

## Vocabulary

Three words carry the results, and they follow the earlier studies rather than inventing
new terms:

| Term | Meaning | Target |
|---|---|---|
| **bias** | how far the posterior median sits from the truth, on average | 0 |
| **accuracy** | how often the 90% interval contained the truth. What the statistical literature calls *coverage* | 0.90 |
| **precision** | how wide that interval was | narrower — but only once accuracy is on target |

A narrow interval that misses the truth is worse than a wide one that does not, so precision
is never read on its own and carries no target line in the figures.

### Why the summaries carry an interval

A **Monte Carlo simulation** computes something with random numbers that algebra cannot
reach. There is no formula for "what fraction of datasets would this model's interval
cover", so the study generates random datasets and counts. The answer is therefore
approximate, and the **Monte Carlo standard error** measures that approximation: how
precisely a summary is pinned down, given the study stopped at `n_datasets` rather than
running forever. It shrinks as `1/sqrt(n)` and says nothing about the world — only about
whether a gap between two models is real or is simulation noise.

Two Monte Carlo layers are stacked here: the random datasets per cell, which is what the
reported error covers, and the MCMC draws inside each fit, which is folded into that fit's
numbers. MCMC is itself a Monte Carlo method.

**Accuracy is a proportion, so it gets a Jeffreys interval instead** — `Beta(0.5, 0.5)`
prior, hence a `Beta(k + 0.5, n - k + 0.5)` posterior, reported as its equal-tailed
quantiles. The usual normal approximation `p +/- 2 sqrt(p(1-p)/n)` breaks at the
boundaries: at accuracy 0.98 with n = 100 it draws a bar running past 1.0. Accuracy near 1
is the expected result for the date-marginalised model, so that is not a corner case.
Jeffreys stays inside [0, 1] and goes asymmetric near the edges.

The figures therefore carry **two kinds of interval**, and the caption says which is which:
Jeffreys on accuracy, +/- 2 Monte Carlo standard errors elsewhere.

Note that 0.88 is exact as a description of one run of 100 datasets — 88 of them did contain
the truth, and you can count them. The interval is there because the paper's claim is about
the *method*, not about those 100 seeds: run 100 different datasets and the number moves.

## What "date-marginalised" means, and what is marginalised over

**Over the unknown calendar date of each object.** One date per object, each integrated out
separately.

The word is the only jargon here and the idea underneath is not. Every object's true date is
unknown, but its dating evidence says how likely each candidate year is. Rather than picking
one year and pretending it is certain, the model asks *"how well does this trend fit, if the
object dates to year t?"* for every candidate year, and averages the answers weighted by how
likely each year is. The date then disappears from the result — that is what "marginalising
it out" names. In the Stan code this is the `log_sum_exp`, which is just a way of adding
probabilities stored as logarithms without them underflowing to zero.

Two things this is **not**:

- Not an approximation. Adding up over all the possibilities is the definition of the
  quantity wanted, not a shortcut to it. The only approximation is the grid spacing, and the
  5-yr vs annual comparison measures it (negligible: 0.01967 against 0.01970).
- Not a departure from errors-in-variables. It is EIV with a **known, fully specified**
  measurement-error distribution instead of an assumed Gaussian one. `sigma` remains the
  scatter of y around the trend and is never the dating error.

*Terminology.* In distributional terms each object's date carries a **categorical** prior
over candidate years — a multinomial with `n = 1`, what NIMBLE writes as `dcat()`. The
structure is exactly that; the reason no distribution statement appears in the Stan code is
that Stan cannot sample discrete parameters, so the categorical is summed over instead of
drawn from. Do not use `multinomial_lpmf`: that function expects counts, and nothing here is
counted.

*If a clearer name is wanted for the paper*, the honest alternatives are **"date-marginalised"**
(explicit about what is integrated out) or **"full-distribution"** (explicit about what it
consumes, and the most transparent for an archaeological readership). `marginal` is kept as
the short label in code and column values either way.

## The comparison models

Three models are compared, referred to by name and never by number:

| Model | The date becomes | Model |
|---|---|---|
| **midpoint** | one number: the middle of the range | `midpoint.stan` |
| **median** | one number: the median of the date distribution | `midpoint.stan`, with `start = end = median` |
| **marginal** | the full per-year probabilities | `marginal_date.stan` |

For a flat window the median *is* the midpoint — `(start + end) / 2` is the median of a
uniform distribution on `[start, end]` — so those two are the same estimator in Cases 2–4
and only separate for radiocarbon. That separation is the point: it isolates "the naive
summary was the problem" from "any single number is the problem".

**Cases 2–4 therefore do not fit median at all.** Doing so would spend a third of the
compute re-deriving algebra and put a duplicate row in every figure. The identity is
asserted rather than assumed: `00_check_median_identity.R` fits both on a small sample and
requires the date columns to match exactly. The surviving row is labelled
**Midpoint / Median** via `read_recovery(..., merge_median = TRUE)`, passed explicitly so
shared code never guesses what a case meant.

### Why "latent" is not reported

`latent_date.stan` samples one date parameter per observation, uniform inside a single
`[start, end]`. It is not in `MODELS` for two reasons:

- **For radiocarbon it would be a straw man.** Feeding it a calibrated posterior means first
  flattening a four-humped curve into one interval — a reduction nobody performs in practice.
- **For flat windows it is redundant.** It encodes the identical assumption to the marginal
  model by a different route, and `00_check_marginal.R` shows the two agree to five decimal
  places. It would be a duplicate column in Cases 2–4.

**The model file stays in `models/` regardless**, because that equivalence is what proves the
marginal model correct and the check needs something to check against. Adding `"latent"` back
to `MODELS` in `fit_models.R` is all it takes to report it again.

## Files

- `models/marginal_date.stan` — primary model. The date is integrated out over its weight
  row rather than sampled, so the posterior has three parameters (`alpha`, `beta`, `sigma`)
  whatever `N` is.
- `models/latent_date.stan` — **not reported** (see above); kept because the
  marginal-vs-latent equivalence is what validates `marginal_date.stan`. Samples one date
  parameter per observation, uniform in `[start, end]`; cannot represent a multimodal
  calibrated posterior.
- `models/midpoint.stan` — the midpoint and median models; they differ only in the date column.
- `scripts/fit_models.R` — `MODELS`, the Stan-data builder per model, one fit to one row.
- `scripts/recovery_summary.R` — metrics, greyscale plot styling, shared by every `04_`.
- `scripts/weight_rows.R` — builds the ragged weight rows for all cases. Grid indexing and
  renormalisation live here and nowhere else.
- `scripts/00_check_marginal.R` — marginal vs latent equivalence check, below.

## The equivalence check

With flat rows the marginalised and latent models encode the same assumption by different routes,
so they have to agree. `00_check_marginal.R` fits both on one Case 2 dataset and compares.
Last run, `N = 200`, true slope 0.02:

```
slope median  marginal 0.01970  latent 0.01970  diff -0.00001 (mcse ~0.00005)
per-date posterior mean: max abs diff 2.48 yr, correlation 1.0000
```

This catches a wrong grid offset or an unnormalised row **in `uniform_rows()`** — neither of
which errors on its own. **Re-run it whenever the weight construction changes.**

It says nothing about `calibrated_rows()`, which it never calls. That path has its own
check, `Sim_Case1_Radiocarbon/scripts/00_check_rows.R`, which calibrates a plateau date and an off-plateau date and
requires each row's median to reproduce `rcarbon`'s calibrated median:

```
1 yr grid  2450+-25 BP: row spans 370 yr (-771..-402), peak -716, row median -564, rcarbon median -564
1 yr grid  2000+-25 BP: row spans 245 yr (-144..100),  peak   15, row median   12, rcarbon median   13
```

**Run both checks after any change to `weight_rows.R`.**

`Sim_Case1_Radiocarbon/scripts/00_check_fit.R` runs the whole path — simulate, calibrate, fit on the calibrated
posteriors — and is the only check that exercises a peaked, multimodal weight row. It found
a bug the flat-row cases structurally could not: a plateau row holds exact zeros in the gaps
between its humps, and `log(0)` broke the date reconstruction, blanking every reported
quantity while the sampler still looked healthy.

Note the boundary: the `00_*` scripts declare their own settings, deliberately, so a check
stays readable on its own. The recovery study reads every setting from
`Sim_Case1_Radiocarbon/data/design.csv` instead, and must never re-declare them.

## Three things that fail silently (all guarded, none obvious)

- **Row order.** `calMatrix` rownames run cal BP *descending* — but cal BP counts backwards,
  so that is already ascending calendar order. The rows need **relabelling
  (`year = 1950 - BP`), not reordering**; reversing them is the bug. Handled once in
  `calmatrix_to_calendar()`, behind an ascending assertion, because it was got wrong here on
  the first attempt and only that assertion caught it.
- **Truncation at the edge of the grid.** `rcarbon` normalises over its whole `timeRange`; a
  row filtered to the study window no longer sums to 1. The renormalisation itself is
  harmless to the trend — scaling a row adds a constant to the log target — but it is the
  visible sign of the thing that *is* harmful: mass cut off the edge of the grid, which
  asserts the object is definitely inside the window and drags its date inward. Systematic,
  edge-only, and in the same direction as the attenuation the paper measures. Hence the
  padding rule for the grid, and `calibrated_rows()` warning above 1% loss.
- **Coarse grids must bin, not subset.** Taking every fifth calibrated year would discard
  four fifths of the mass and resample a multimodal posterior instead of summarising it.
  This cannot go wrong for Cases 2–4, whose rows are flat — which is why it is easy to miss
  when Case 1 inherits the coarse grid.

## Calendar axis

`time_ref_min` and `time_ref_range` come from the design and are passed to all three models
as data. They are deliberately not derived from each dataset's realised dates: with a
per-dataset range, `normal(0, 10)` on `beta` is a different prior on the calendar-scale
slope in every dataset and every model, which would confound the model comparison with the
centring.

The same argument decides Case 1's grid padding. A calibrated date near the edge of the grid
loses part of its posterior, and renormalising then asserts it is certainly inside the
window — dragging edge dates inward, in the same direction as the attenuation being
measured. `c14_grid_pad()` measures the widest calibrated support rather than guessing, and
`Sim_Case1_Radiocarbon/scripts/01_design.R` freezes **one** value (588 yr) across every window and lab error, so
`time_ref_range` is identical (1576 yr) in all 1400 datasets. Per-cell padding would have
made the slope prior differ between cells and broken the comparison it exists to support.

## Cost

The marginal model scans every candidate year in a row at every leapfrog step, so runtime is
linear in candidate years per row and grid resolution is the main cost lever.

| Grid step | Years per row | Per chain | vs latent model |
|---|---|---|---|
| 1 yr | ~120 | ~37 s | ~25× |
| 5 yr | 26 | ~8.5 s | ~7× |

The 5-year grid costs nothing in accuracy on Case 2 — slope median 0.01967 vs 0.01970,
inside Monte Carlo error — so it is the default for Cases 2–4, with the annual grid kept as
the reference. **Case 1 must be re-checked before inheriting this**: a multimodal plateau
posterior has structure at a scale a flat window does not, and it is the case a coarse grid
could smear.

Set the step with `MARGINAL_GRID_STEP=5 Rscript shared/scripts/00_check_marginal.R`.
