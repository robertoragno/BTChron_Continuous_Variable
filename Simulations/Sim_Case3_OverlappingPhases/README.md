# Case 3 — overlapping phases

The method, the model and the comparison models are documented once in
`Simulations/shared/README.md`. This file covers only what is specific to
overlapping phases, and one measured finding that changes what the case claims.

## What makes this case different

Phases come from the same Dirichlet broken stick as the archived `Sim_Linear`,
but adjacent phases share a strip rather than meeting at a point, as ceramic
phases do. A date in the shared strip is labelled one of the two phases. Two
factors:

- **overlap** — `0`, `0.25`, `0.5`: the fraction of a phase's own width by which
  its recorded window is stretched past each abutting boundary.
- **assign_p** — `0.5`, `0.35`, `0.2`: for a date in the earlier+later overlap,
  the probability it is labelled the *earlier* phase. `0.5` is a coin flip; lower
  values lean assignments toward the later phase, the usual "ambiguous sherd goes
  to the better-attested phase".

The window is flat, so the median IS the midpoint here, exactly as in Case 2.
Only `midpoint` and `marginal` are fitted, and the row is labelled
"Midpoint / Median". Case 2's `00_check_median_identity.R` proves that identity
once for every flat-window case.

## The finding that changes the claim

**The plan (section 5.1) says overlapping phases bias the point-date *slope*. They
do not. Measured, not argued — `00_check.R` step 3, OLS only, no Stan:**

| overlap | assign_p | midpoint slope ratio |
|---|---|---|
| 0.00 | 0.50 | 1.000 |
| 0.25 | 0.50 | 1.000 |
| 0.50 | 0.50 | 0.996 |
| 0.25 | 0.35 | 0.999 |
| 0.50 | 0.35 | 0.994 |
| 0.25 | 0.20 | 1.002 |
| 0.50 | 0.20 | 0.995 |

The ratio is 1 at every level of both factors. The reason is simple once written
down: `assign_p` below 0.5 does thin each phase's early-edge strip less than its
late-edge strip, so the retained true dates are no longer symmetric about the
window midpoint — the depletion figure shows this clearly. But the shift is
roughly the *same number of years for every phase* (it scales with the phase's
own width via the overlap reach, and broken-stick width is uncorrelated with
position on the timeline). A constant shift `d` in every phase's retained mean
gives `E[true | midpoint] = midpoint + d`: it moves the regression **intercept**,
not the slope. So Berkson fails locally and the slope stays unbiased.

### Why not the clamped window

An earlier version stretched the windows but clamped the two end phases at
`[t_min, t_max]`. That produced a midpoint slope ratio of ~1.07 at overlap 0.5 —
apparent slope *inflation* — which looked like a mechanism and was not. Clamping
the first phase on its left while stretching it right moves its midpoint later;
the last phase's moves earlier; the midpoint range compresses against the
true-date range and the OLS slope inflates. This is the same artifact Case 2
found and rejected when it chose overhang over clamped windows
(`Sim_Case2_Typochronology/scripts/01_design.R`). The windows now overhang
`[t_min, t_max]`, and the grid is padded by the widest reach, the way Case 2 pads
by its widest window.

### What Case 3 actually shows

The per-phase Berkson violation is real, and `assign_p` controls it. It surfaces
where Case 2's headline surfaces: **inflated `sigma` and credible intervals that
are too narrow** — Claim A, not Claim B. `assign_p` is the factor that matters,
`overlap` is its amplifier, and `sigma` error and coverage are the headline
metrics. The slope is a null here and is reported as one.

Section 5.1 of the plan needs correcting to match; open question Q3 (symmetric vs
one-sided overlap) is moot — neither orientation breaks the symmetry — and can be
closed.

## Calendar axis

`time_ref_min` / `time_ref_range` come from the design and are passed to all
models as data, same argument as Case 1 and Case 2. Case 3 needs padding because
the stretched windows overhang `[100, 900]`; `01_design.R` freezes one pad value
(the widest reach over every overlap level) so `time_ref_range` is identical
across cells and `normal(0, s)` on beta means the same thing on the calendar
scale everywhere.

Note this pad differs from Case 2's, so `time_ref_range` is not identical between
the two cases. That is fine within the study — section 9 only requires it to be
fixed *within* a case — but a Case 2 vs Case 3 slope number is not exactly
like-for-like.

## Checks

```
Rscript Simulations/Sim_Case3_OverlappingPhases/scripts/00_check.R
```

Asserts that `overlap = 0` reproduces the archived `partition_timeline()` year
for year (the plan's section 12 gate), that every window contains its own true
date with gap-free coverage, prints the slope-ratio table above, and writes the
retained-date depletion figure.

## Running the study

```
Rscript Simulations/Sim_Case3_OverlappingPhases/scripts/01_design.R
Rscript Simulations/Sim_Case3_OverlappingPhases/scripts/03_recovery_study.R
Rscript Simulations/Sim_Case3_OverlappingPhases/scripts/04_recovery_plots.R
```

`partition_timeline()` lives in `Simulations/shared/scripts/partition.R`, extended
from the archived copy with the `overlap` and `assign_p` arguments and carrying
the abutting boundaries as an attribute for Case 4.
