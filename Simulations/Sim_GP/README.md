# Simulation 3: Gaussian Process

This page shows the main findings from the simulated dataset of 300 samples of a continuous variable ranging approximately between ~5 and ~25 (Figure A). Each sample is assigned a time range with a `Start_Date` and an `End_Date`, with a known true date assigned between these boundaries (Figure B). 
The general structure of the dataset assumes an underlying temporal trend, with a peak around 450 CE and a decline thereafter (Figure C):
`f(t) = baseline + amplitude * exp(-((t - peak) / width)^2)`

Given the simulated nature of the dataset, the goal is to verify that the model recovers the known generating parameters:

- **Baseline** = 8
- **Amplitude** = 12
- **Peak** = 450 CE
- **Width** = 200
- **Noise** ~ N(0, 2.5)

## Exploratory panel
<p align="center">
<img src="figures/exploratory_panel.png" height="700" text-align="center"/>
</p>

## Prior predictive check

Before fitting any data, 200 parameter sets were drawn from the priors and GP realisations were generated through the HSGP parameterisation. Panel **A** shows the resulting trend curves (trimmed to the 90th percentile by absolute magnitude for legibility); the blue band marks the observed data range. Panel **B** compares the density of prior predictive observations against the actual observed values. The priors are intentionally wide — they should cover the data range without being tightly centred on it.

<p align="center">
<img src="figures/prior_predictive_check.png" height="700" text-align="center"/>
</p>

## Model 1: HSGP regression with latent date inference

A Hilbert-space GP regression with latent date inference. Each observation's true date is modelled as a uniform draw within its `[Start_date, End_date]` window, and the trend is approximated via a low-rank basis expansion of a squared-exponential GP:

```
f(t) = mu + PHI(t)' * beta
y_n  ~ Normal(f(true_date_n), sigma)
```

Panel **A** shows the recovered trend (median and 50%/90% credible intervals) against the true Gaussian bell curve used to generate the data (dashed). Panel **B** is a date recovery plot: each dot is one sample, with the true generating date on the x-axis and the posterior median on the y-axis — vertical bars are 90% CIs. Panel **C** shows posterior distributions for sigma (dashed line marks the true generating value), the GP mean mu, and the GP length-scale rho in original time units; mu and rho do not have simple true values since they are inferred non-parametric properties of the fitted curve.

<p align="center">
<img src="figures/model_results_panel.png" height="700" text-align="center"/>
</p>

## Model comparison: latent dates vs midpoint dates

A second model fixes each observation's date at the midpoint of its `[Start_date, End_date]` window and runs the same HSGP regression on those fixed dates, ignoring temporal uncertainty entirely.

Panels **A–B** compare trend recovery between the two models. Panel **C** shows the posterior for sigma under each model, with the true generating value (2.5) marked by a dashed line. The midpoint model places its sigma posterior noticeably above 2.5 because unmodelled date uncertainty is absorbed into the residual variance — the same pattern seen in the linear and changepoint simulations.

<p align="center">
<img src="figures/model_comparison.png" height="700" text-align="center"/>
</p>

## Simulation-based calibration (SBC)

To check whether the credible intervals are well calibrated, the simulation was repeated 100 times — each time redrawing true dates and observation noise from the same generating process while keeping the date windows fixed.

Sigma is the primary calibration target because it has a known true value (2.5). Under a well-calibrated model the true value is equally likely to land anywhere in the posterior, so posterior ranks across 100 replications should be roughly uniform (dashed line marks the expected bin count of 5). A histogram piled up on the right means the model systematically underestimates sigma — the pattern expected for the midpoint model.

<p align="center">
<img src="figures/calibration_coverage.png" height="500" text-align="center"/>
</p>
