# Simulation 1: Linear Regression

A continuous measurement is assumed to follow a linear temporal trend
$f(t) = \text{baseline} + \text{slope} \cdot t$, observed with Gaussian noise. Each observation is
dated only to a phase with `[Start_date, End_date]` rather than to a known year. Does treating each date as a latent parameter within its phase recover the trend better than collapsing the phase to its midpoint?

As a single example, we could say that: baseline = 5, slope = 0.015,
noise $\sim \mathcal{N}(0, 1.5)$. The timeline is 100–900 CE. The worked example uses a coarse
periodisation, 3 uneven phases with one covering more than half the range. This
is deliberate, as when phases are narrow the two models mostly agree, whereas with a coarse
example we show differences between two models. 

## Periodisation and dating resolution

We assume that, as common in archaeological relative dating, dates vary in the timeline by phase ("The sample X belongs to phase Y, which has a known duration `[Start_date, End_date]`"). To simulate this we divide the timeline into `K` phases by a Dirichlet broken stick, and every
observation inherits the `[Start_date, End_date]` of the phase its true date
falls in. In this way, observations in the same phase share an identical window and the
dating uncertainty is structured rather than random (as would be the case instead for a radiocarbon date). Each dataset draws its own number of phases ($K \sim \mathcal{U}\{3, \ldots, 10\}$) and Dirichlet concentration
($\alpha_{\text{conc}} = 10^{\,2 \cdot \text{Beta}(2,\, 3.5)\, -\, 1}$)[1], which could result either in few even phases or one long phase dominating others in the range. 

The resolution of a periodisation is summarised by the Shannon entropy of the
phase weights, $H = -\sum_k p_k \log p_k$. Even phases have high H (the
ceiling is about 2.3 for ten equal phases), whereas one dominant phase would have instead low H. We also run a
prior-predictive check located in `00_periodisation_check.R`.

*[1] With Beta(2, 3.5) the distribution leans toward low concentrations (coarser, lower-H periodisations). alpha_conc is 10 exponentiated to this Beta so that we get values that are between 0 and 10, but the mass is placed more towards coarser phases.*

<p align="center">
<img src="figures/periodisation_H_distribution.png" height="320" text-align="center"/>
<img src="figures/periodisation_examples.png" height="520" text-align="center"/>
</p>

## Scripts

`simulate.R` contains the shared generating process `simulate_linear()`, the
`partition_timeline()` helper (note to self: might move this to a separate file), and the timeline constants — sourced by the
example and the study so they use the same process.

| Script | Purpose | Output |
|---|---|---|
| `00_periodisation_check.R` | prior predictive check: draw the priors, simulate datasets, show H, example phases, and the values they imply | `figures/periodisation_H_distribution.png`, `periodisation_examples.png`, `prior_predictive.png` |
| `01_example.R` | one dataset, shown and fit with both models | `figures/exploratory_panel.png`, `model_comparison_single_fit.png`, `individual_date_posteriors.png` |
| `02_recovery_study.R` | many random datasets (each its own intercept, slope, noise, sample size, periodisation), both models fit, intervals recorded | `output/recovery_results.csv` |
| `03_recovery_plots.R` | the study figures | `figures/recovery.png`, `accuracy.png`, `accuracy_precision_composite.png`, `precision_vs_entropy.png`, `sigma_vs_entropy.png`, `accuracy_vs_entropy.png`, `precision_boxplots.png`, `precision_vs_n.png` |
| `05_scenarios.R` | four fixed archetype cases (fine dating, coarse dating, skewed deposition, diagnostic sampling), many datasets each, both models fit; also draws one example dataset per case | `output/scenarios.csv`, `figures/scenario_examples.png` |
| `06_scenario_plots.R` | the scenario figures | `figures/scenario_accuracy_precision.png`, `scenario_sigma.png`, `output/scenario_table.csv` |
| `07_berkson_vs_classical.R` | didactic: classical vs Berkson measurement error, side by side, to show why the slope ties between models | `figures/berkson_vs_classical.png` |

The example (`01`) fits in less than a minute; the study (`02`) takes a few minutes (around 10 circa), so
plotting (`03`) is kept separate and figures can be restyled without refitting.
Sample size is swept across $N \in \{50, 100, 200, 400\}$ to read precision against N.

## Exploratory panel
This is a single example that is used as a figure in the paper to illustrate the data generating process.
<p align="center">
<img src="figures/exploratory_panel.png" height="700" text-align="center"/>
</p>

## Errors-in-variables (EIV) inference

The EIV model treats each true date as a latent parameter, uniform within its
phase. The regressor is known only to a phase, not to a year:

$$
\begin{aligned}
\text{trend}(x) &= \alpha + \beta x \\
y_n &\sim \mathcal{N}\!\left(\text{trend}(t_n),\ \sigma\right)
\end{aligned}
$$

where $t_n$ is the latent true date of observation $n$.

Left panels: the recovered trend, with a sample's phase shaded and its observed
value marked on a labelled dashed line. Right panels: each date's posterior, with
the phase window (dashed grey) and the true date (red).

These panels are only valid for the EIV model, as the midpoint model cannot produce them as it assigns each sample the centre of its phase. 
The EIV model uses the observed value together with the
trend to place the object inside its phase: a sample known only to fall somewhere in the widest phase comes back with a date estimate that is more accurate. This estimate is pulled late when its value is high and early when it is low. It does not always work (see sample #195, in which the noise pushed its value up), but the true date stays inside the range.

<p align="center">
<img src="figures/individual_date_posteriors.png" height="900" text-align="center"/>
</p>

## Model comparison on one dataset

The midpoint model fixes each date at the centre of its phase. Both recover the
trend, and here they even agree on how uncertain the trend is. The slope leans on
how far apart the dates sit across the whole 100–900 range, and blurring each
date by a century or two barely moves that spread, so both models pin the slope
about equally. However, the midpoint cannot absorb the within-phase date spread, so it reads that spread as noise and returns a larger sigma. On this coarse
example the midpoint puts sigma near 2.1 against a true 1.5, while the EIV model places sigma closer to 1.5.

<p align="center">
<img src="figures/model_comparison_single_fit.png" height="800" text-align="center"/>
</p>

## Recovery study

The dataset above is just one example for the paper; the real study draws many datasets — each with
its own intercept, slope, noise, sample size, and periodisation. Both models (EIV and midpoint) are
fit to all of them, recording how often the 50% / 90% interval holds the truth
(accuracy) and how wide the interval is (precision).

*Preliminary results*. Both models recover the intercept and slope with about the
right coverage. For a straight line the midpoint is the average of the possible
dates within a symmetric phase, so it does not bias the trend.

Berkson error? Classical ME adds independent noise to a known
true value, and that does bias a regression, flattening the slope. Here however we fix the phase window first and the true date is uniform inside this window, so the midpoint's error are symmetric around a known point and cancel out for a straight line. This might be why the slope correctly survives in the midpoint route.

<p align="center">
<img src="figures/berkson_vs_classical.png" height="340" text-align="center"/>
</p>

Each thin grey line links one observation to itself: the dot is where it
truly was, the square is what you'd actually have to work with for that same
point (same y, shifted x), and the line just shows how far it moved. Same
amount of x-uncertainty in both panels, opposite outcome. Left: a true value
blurred by independent noise flattens the fitted line to about half its true
slope. Right: a fixed window with the true value uniform inside it (exactly
the phase-dating setup used here) leaves the fitted line almost untouched.

The EIV model however gains precision, with its slope interval
about 18% narrower and its intercept interval about 17% narrower, because the
midpoint model throws away the within-phase spread.

Sigma is where the two models differ substantially. The midpoint reads the within-phase date
spread as measurement noise and overestimates it, so its 90% interval holds the
true sigma only about 45% of the time against about 90% for EIV. This behaviour is expected since the midpoint can only place the date uncertainty in the residuals. As phases coarsen, the midpoint gets worse: from
about 26% coverage at the lowest H up to about 66% at the finest. Precision improves with sample size for both models.

The short version: when dating is fine the two models agree on the trend and slope and the midpoint is a fair shortcut. If phases are coarser, EIV can be considered: it provides an honest sigma and narrower trend intervals. Moreover, the per-sample date estimates can be very useful and cannot be produced by the midpoint.

<p align="center">
<img src="figures/recovery.png" height="320" text-align="center"/>
</p>

<p align="center">
<img src="figures/accuracy.png" height="340" text-align="center"/>
</p>

Overall calibration. For each parameter, the number of datasets whose 50% and 90%
interval (dashed line) actually held the true value. The whiskers are 95% Jeffreys intervals.

### Accuracy and precision against dating resolution

These plots show if an interval is accurate (does the 90% interval hold the truth?) and precise (width of the 90% interval). Panel A shows the accuracy and panel B shows the precision, both across the entropy range.


<p align="center">
<img src="figures/accuracy_precision_composite.png" height="620" text-align="center"/>
</p>

Again, no much difference between the two models for the intercept and slope, what changes is sigma. For coarser phases, sigma is better captured by the EIV model: 

<p align="center">
<img src="figures/precision_boxplots.png" height="340" text-align="center"/>
<img src="figures/sigma_vs_entropy.png" height="340" text-align="center"/>
</p>


### Further checks: Sample size

We also check how sample size affects the precision: both models show a steady improvement in precision as sample size increases. 

<p align="center">
<img src="figures/precision_vs_n.png" height="320" text-align="center"/>
</p>

## Scenario study

The recovery study draws fully random datasets and shows the general
relationship. When should we use EIV instead of midpoints? Is there an advantage? We test this with the script 
`05_scenarios.R` with four cases with random intercept, slope and sigma so that each case shows:

- **Fine dating**: many even phases, so the dates are tight.
- **Coarse dating**: few uneven phases, one covering more than half the range.
- **Skewed deposition**: coarse dating with the widest phase late on the
  timeline and dates piling early inside it (deposition follows a Beta(1.5, 5)
  over the timeline). The one configuration the data-generating-process checks
  singled out as able to bend the trend, because the within-phase dates sit
  away from the midpoint in a way that grows with time.
- **Diagnostic sampling**: a periodisation with wide and narrow phases, but the
  narrow ones are over-sampled (as when recognisable ceramics date a short
  phase), so the sample lands in the well dated phases and still spans the
  timeline.

Small samples are not a case of their own; the recovery study already sweeps N.

What the four cases look like as datasets, one draw each with the
worked-example parameters. True dates are dots, midpoints squares. In the
skewed deposition panel the dots pile toward the early edge of the wide late
phase while the midpoints sit at its centre, which is the whole problem:

<p align="center">
<img src="figures/scenario_examples.png" height="560" text-align="center"/>
</p>

The figure reads top over bottom: accuracy (does the 90% interval hold the
truth) above precision (how wide it is), so the two must be read together. A
narrow interval that sits on the wrong value is not a good one, and a narrow
interval that still covers is the win.

<p align="center">
<img src="figures/scenario_accuracy_precision.png" height="620" text-align="center"/>
</p>

When the dates are good (fine dating, and diagnostic sampling once the sample
concentrates in the narrow phases) the two models agree on everything, sigma
included, and the midpoint is a fair shortcut. Under coarse dating they separate:
the EIV slope interval is about twice as narrow (0.0035 against 0.0075) at
comparable coverage, and the EIV sigma sits on the truth while the midpoint
reads it far too high (+1.2 on average), so the midpoint sigma interval holds
the truth only 24% of the time against 92% for EIV.

The skewed deposition case produces worse results for both models. When the dates inside a wide late phase actually pile toward its early edge, the
midpoint places them at the centre (by definition) and flattens the slope to about half its
true value. Similarly, the EIV model does not improve (slope ratio 0.53 against 0.54,
90% coverage 17% for both, intercept coverage similarly collapsed). The reason
is that EIV's within-phase prior is uniform, which is the same wrong
assumption in a softer form, and the observed values are not informative
enough to overrule it. EIV still reads the noise more honestly (77% against
45% sigma coverage), but the trend is lost for both models. In practice, with coarse dating and a deposition process drifting within phases, the trend cannot be easily recovered by both models.

Moreover, under coarse dating some EIV fits do not converge and are
dropped (87 of 100 survive the R-hat and divergence filter). The midpoint sigma interval is a slightly narrower than EIV's, but it is narrow around the
wrong value. The direction of that error,
the midpoint always reading sigma too high, is clearest here:

<p align="center">
<img src="figures/scenario_sigma.png" height="360" text-align="center"/>
</p>
