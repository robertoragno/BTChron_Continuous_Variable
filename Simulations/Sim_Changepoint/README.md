# Simulation 2: Changepoint Regression

This page shows the main findings from the simulated dataset of 300 samples of a continuous variable ranging approximately between ~7 and ~18 (Figure A). Each sample is assigned a time range with a `Start_Date` and an `End_Date`, with a known true date assigned between these boundaries (Figure B). 
The general structure of the dataset assumes an underlying segmented temporal trend (i.e. with a change in its slope) with a single changepoint (Figure C):
`f(t) = baseline + slope_1 * (t - t_min)` for `t <= changepoint`
`f(t) = baseline + slope_1 * (changepoint - t_min) + slope_2 * (t - changepoint)` for `t > changepoint`

The main goal is to verify that the model recovers the known generating parameters that we used in the simulation, listed below:

- **Baseline** = 8
- **Slope 1** = 0.02
- **Slope 2** = -0.01
- **Changepoint** = 500 CE
- **Noise** ~ N(0, 1.5)

## Exploratory panel
<p align="center">
<img src="figures/exploratory_panel.png" height="700" text-align="center"/>
</p>

## Model 1: Changepoint regression with latent dates
This is a Bayesian segmented regression with a latent date, as per the linear model. Each observation's 'true' date is modelled as a uniform draw within its `[Start_date, End_date]` window. The trend is a continuous piecewise-linear function (basically a broken-stick parameterisation):

```
mu(t) = alpha + beta1*t + (beta2 - beta1) * max(0, t - changepoint) // mean
y_n ~ Normal(mu(true_date_n), sigma)
```

where `alpha` is the intercept, `beta1` and `beta2` are the slopes before and after the changepoint, and `sigma` is the observation noise. Before the changepoint, `max(0, t - changepoint) = 0` and the mean reduces to `alpha + beta1*t`. After it, the slope shifts to `beta2`, with the two segments meeting continuously at the changepoint — so `beta2 - beta1` is the change in slope.

The changepoint is a free parameter with a uniform prior over the normalised time range. Stan treats it as unknown and samples a posterior for it alongside all other parameters. For each candidate value, the piecewise mean changes shape and the likelihood scores how well that shape fits the data — so the posterior concentrates where the slope change is most supported by the observations. The true value (500 CE) is only used afterwards to check whether the posterior recovered it correctly.

Panel **A** shows the recovered trend (median and 50%/90% credible intervals) against the true generating function; the rug on the right margin shows the distribution of observed values. Panel **B** is a date recovery plot: each dot is one sample, with the true generating date on the x-axis and the posterior median date on the y-axis — if recovery were perfect, all points would sit on the dashed 1:1 line; vertical bars show the 90% CI. Panel **C** shows posterior distributions for all five generating parameters, with graduated shading by credible interval and dashed lines marking the true values.

<p align="center">
<img src="figures/model_results_panel.png" height="700" text-align="center"/>
</p>

## Model comparison: latent dates vs midpoint dates

A second model uses the midpoint of each sample's `[Start_date, End_date]` window as a fixed date, with no latent date inference. This is the conventional approach — treating the midpoint as the "best guess" and ignoring date uncertainty.

Panels **A–B** compare the trend recovery of both models. Panel **C** shows the posterior distributions for each generating parameter side by side (latent-date model on top, midpoint model on bottom), with graduated shading by credible interval and dashed lines marking the true values.

<p align="center">
<img src="figures/model_comparison.png" height="800" text-align="center"/>
</p>

## Simulation-based calibration (SBC)

To verify that the credible intervals are well calibrated, both models were each run 100 times — each time redrawing true dates and observation noise from the same generating process while keeping the timespans fixed.

Each replication produces (for each parameter) 4,000 posterior draws. The **rank** is where the true value falls among those draws — i.e., how many of the 4,000 draws are below the true value. If the model is well calibrated, the true value is equally likely to land anywhere in the posterior, so across 100 replications the ranks should be roughly uniformly distributed. The histograms below have 20 bins, so each bin should contain about 100 / 20 = 5 replications (dashed line). A flat histogram means the model's uncertainty is honest; a histogram piled up on one side means the model is systematically wrong about that parameter.

The figure has two rows: the **latent-date model** (top) and the **midpoint model** (bottom). Coverage rates for the 50% and 90% credible intervals are annotated in each panel. The latent-date model produces approximately flat rank histograms across all five parameters, with coverage rates close to their nominal levels. The midpoint model is severely miscalibrated: baseline, slope 1, and slope 2 show skewed rank distributions, and sigma collapses entirely to rank 0 in every replication (50% CI: 0%, 90% CI: 0%), indicating that the posterior for sigma is systematically too narrow when date uncertainty is ignored.

<p align="center">
<img src="figures/calibration_coverage.png" height="500" text-align="center"/>
</p>
