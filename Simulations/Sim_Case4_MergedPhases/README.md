# Case 4 - nested / merged phases

The method, the model and the comparison models are documented once in
`Simulations/shared/README.md`. This file covers only what is specific to nested
phases.

## What makes this case different

A find is dated to a broad period that spans several adjacent fine phases. The
broad period is the union of those phases, so it is still one interval, just
wider than a single phase. The true date stays uniform inside it, so Berkson
holds and the midpoint is unbiased on the slope, as in Case 2.

What the case shows is the **mix of dating resolutions within one dataset**. Each
find draws its own broad-period width from `1..merge_max` fine phases: some keep
a fine-phase label, others sit only in the broad period. That is what an
excavation report looks like, and the usual response is to analyse the finely
dated and broadly dated finds separately, or to drop the coarse ones. One weight
row per find keeps both in the same model, tight for the fine-phase finds and
wide and flat for the broad-period ones.

The expected result is the sigma / coverage story of Cases 2 and 3: the point-date
inflates sigma and its intervals are too narrow on the broadly dated finds, while
the date-marginalised model recovers sigma and still uses the finely dated finds
correctly. The slope is a null here and is reported as one.

## Factors

- **merge_max** - `2`, `3`, `4`: the largest number of fine phases a broad
  period can span. Higher values put coarser finds in the mix. `3` is the
  reference for the core `N` and slope-condition sweep.

The fine phases come from the same Dirichlet broken stick as Case 3
(`shared/scripts/partition.R`, `overlap = 0`). Each find's broad period is a run
of `m` consecutive fine phases positioned randomly around its own phase, so
window width does not track calendar position.

## Median is not fitted

For a flat window `(start + end) / 2` is the median, so it duplicates the
midpoint row. `Sim_Case2_Typochronology/scripts/00_check_median_identity.R`
proves that once for every flat-window case. The surviving row is labelled
"Midpoint / Median".

## Calendar axis

The merged windows never leave `[100, 900]`, so the axis needs no padding.
`time_ref_min` / `time_ref_range` are `100` / `800`, passed to Stan as data and
identical for every dataset and model.

## Running

```
Rscript Simulations/Sim_Case4_MergedPhases/scripts/01_design.R          # writes data/design.csv, runs generator checks
Rscript Simulations/Sim_Case4_MergedPhases/scripts/03_recovery_study.R  # the fits
Rscript Simulations/Sim_Case4_MergedPhases/scripts/04_recovery_plots.R  # tables and figures
```

Same three stages as every other case. `03_` skips the fits if
`output/recovery_results.csv` already exists, so re-running after a crash
resumes rather than discarding hours of work; delete that file to refit.
Shrink a development run with `RECOVERY_LIMIT`, `RECOVERY_WARMUP`,
`RECOVERY_SAMPLING`, `RECOVERY_OUT`.
