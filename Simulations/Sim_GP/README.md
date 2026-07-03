# Simulation 3: Gaussian Process

A continuous measurement is assumed to follow a smooth but wavy temporal trend
`f(t)`, observed with Gaussian noise. Each observation is dated only to a phase
with `[Start_date, End_date]` rather than to a known year. Does treating each
date as a latent parameter within its phase recover the curve better than
collapsing the phase to its midpoint?

Unlike the linear and changepoint simulations, the trend here has no small set of
generating parameters to recover one-for-one: the fitted model is a
non-parametric Gaussian process. The truth is instead a wavy curve — a baseline
plus a few sine waves — and the question is asked of the whole curve.

```
f(t) = baseline + Σ amp_k · sin(2π (t − 100) / period_k + phase_k),   k = 1..3
y_n  ~ Normal(f(true_date_n), sigma)
```

Each dataset draws its own curve: `baseline ~ U(2, 15)`, three amplitudes
`~ U(2, 8)`, three periods of roughly 200–530 years, and random phases. The
periods are deliberately kept on the same scale as the dating phases (a few
hundred years): a curve that bends inside a dating window is where the midpoint
shortcut can go wrong, whereas a curve much flatter than the windows would behave
like the linear case and one much wavier could not be recovered by any model. The
truth is built from sine waves while the model assumes a generic smooth GP, so
recovering it is a fair test rather than one rigged in the model's favour.

## Periodisation and dating resolution

The dating step is identical to the linear and changepoint simulations. The
timeline (100–900 CE) is divided into `K` phases by a Dirichlet broken stick, and
every observation inherits the `[Start_date, End_date]` of the phase its true
date falls in, so observations in the same phase share an identical window and the
dating uncertainty is structured rather than random. Each dataset draws its own
number of phases (`K ~ U{3, …, 10}`) and Dirichlet concentration
(`alpha_conc = 10^(2·Beta(2, 3.5) − 1)`), giving anything from a few even phases
to one long phase dominating.

The resolution of a periodisation is summarised by the Shannon entropy of the
phase weights, `H = −Σ p · log p`. Even phases have high H (the ceiling is about
2.3 for ten equal phases); one dominant coarse phase gives low H. Higher H means
finer, more balanced periods and better-resolved dating, so H is the natural axis
to read accuracy and precision against. A prior-predictive check is in
`00_periodisation_check.R`.

## Scripts

`simulate.R` holds the shared generating process: `simulate_gp()`, the wavy-curve
generator (`draw_curve` / `f_true`), and the `partition_timeline()` helper shared
with the other simulations. It is sourced by both the example and the study so
they use the same process.

| Script | Purpose | Output |
|---|---|---|
| `00_periodisation_check.R` | prior predictive check: draw the priors, simulate datasets, show H, example phases, and the curves and values they imply | `figures/periodisation_H_distribution.png`, `periodisation_examples.png`, `prior_predictive.png` |
| `01_example.R` | one dataset, shown and fit with both models | `figures/exploratory_panel.png`, `model_comparison_single_fit.png`, `individual_date_posteriors.png` |
| `02_recovery_study.R` | many random datasets (each its own curve, noise, sample size, periodisation), both models fit, curve and sigma recorded | `output/recovery_results.csv` |
| `03_recovery_plots.R` | the study figures and the recovery table | `figures/accuracy_vs_entropy.png`, `precision_vs_entropy.png`, `precision_vs_n.png`, `sigma_vs_entropy.png`, `curve_recovery.png`, `output/recovery_table.csv` |

The latent-date GP is the heaviest of the three models to fit (several minutes per
dataset), so the study (`02`) is the expensive step and is run in a tmux session;
plotting (`03`) is kept separate so figures can be restyled without refitting.
Sample size is swept across `N = 50, 200` (low / high) to read precision against N.

## The model

A Hilbert-space approximation to a squared-exponential GP (Solin & Särkkä 2020):
the trend is `f(t) = mu + PHI(t)' · beta`, a low-rank basis expansion that avoids
inverting an N×N covariance matrix. The EIV model treats each true date as a
latent parameter, uniform within its phase, inferred jointly with the GP; the
midpoint model fixes each date at the centre of its phase. Both share the same
priors, so the two are directly comparable. The GP length-scale prior
(`rho ~ inv_gamma(5, 0.9)`, ~80–370 years) and noise prior (`sigma ~ normal(0, 5)`
half-normal, over the generating range) are set for this wavy-curve, varying-truth
design; the latent-date model needs `adapt_delta = 0.99` to sample its harder
geometry.

## Exploratory panel

A single example used in the paper to illustrate the data-generating process.

<p align="center">
<img src="figures/exploratory_panel.png" height="700" text-align="center"/>
</p>

## Model comparison on one dataset

The example dataset fit with both models: recovered curve (median and 50%/90%
band) against the true wavy curve, and the sigma / mu / length-scale posteriors.
This is one dataset for illustration only — the systematic result is the recovery
study below.

<p align="center">
<img src="figures/model_comparison_single_fit.png" height="800" text-align="center"/>
</p>

<p align="center">
<img src="figures/individual_date_posteriors.png" height="900" text-align="center"/>
</p>

## Recovery study

The dataset above is just one example; the real study draws many datasets — each
with its own curve, noise, sample size, and periodisation. Both models are fit to
all of them. Because the GP target is the curve, accuracy is measured as the
**proportion of the true curve inside the credible band**, averaged over a time
grid, at the 50% and 90% levels; precision is the mean band width. `sigma` is kept
as a scalar target, as in the other two simulations, since it carries the
"midpoint mistakes dating spread for noise" story and links the three studies.

Convergence is checked per target, because the two can fail independently: `03`
keeps a fit for the curve figures only if the grid curve (`mu_pred`) mixed
(Rhat ≤ 1.05, no divergences), and for the sigma figure only if sigma mixed, and
reports how many were kept.

### The noise scale is the headline result

As in the linear and changepoint simulations, the robust and consistent finding
is on **sigma**: the midpoint model reads the within-phase date spread as
measurement noise and overestimates it (roughly doubling the true value under
coarse dating), while EIV recovers it. The error grows as phases coarsen (low H).
This is the result the whole simulation section rests on, and it holds across
every fit.

<p align="center">
<img src="figures/sigma_vs_entropy.png" height="340" text-align="center"/>
</p>

### The curve, where it is identifiable

The curve is reported over the subset of datasets where the latent-date GP curve
converged. It is **not** recovered everywhere: with dates left latent, the GP
length-scale is only weakly identified — the data can be explained by a smoother
curve with dates arranged one way or a wavier curve with dates arranged another —
so on coarser-dated datasets the curve posterior is multimodal and does not
converge (chains settle on different curves). This is a property of the
latent-date GP itself, not of the periodisation redesign: the earlier
single-bell study only ever checked sigma's convergence, so the same weak
identifiability was present but unmeasured. Where the curve *does* converge (finer
dating, higher H), EIV matches the midpoint's point accuracy with a narrower,
better-calibrated band; the midpoint's occasional higher coverage is bought with
wider bands and an inflated sigma, not better point accuracy. How deep to take the
GP curve validation (up to full simulation-based calibration) is an open question
for the professor.

<p align="center">
<img src="figures/accuracy_vs_entropy.png" height="340" text-align="center"/>
</p>

<p align="center">
<img src="figures/precision_vs_entropy.png" height="320" text-align="center"/>
<img src="figures/curve_recovery.png" height="320" text-align="center"/>
</p>

### Sample size

Precision against sample size, as a sanity check that more data tightens the
intervals for both models.

<p align="center">
<img src="figures/precision_vs_n.png" height="320" text-align="center"/>
</p>

## Scope

This is a recovery study (simulate → fit → check recovery across many datasets),
not simulation-based calibration from the GP prior. Full SBC remains an open
option for the methods appendix, to be decided with the professor.
