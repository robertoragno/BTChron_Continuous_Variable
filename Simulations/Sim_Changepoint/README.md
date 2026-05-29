# Simulation 2: Changepoint Regression

This page shows the main findings from the simulated dataset of 300 samples of a continuous variable ranging approximately between ~7 and ~18 (Figure A). Each sample is assigned a time range with a `Start_Date` and an `End_Date`, with a known true date assigned between these boundaries (Figure B). 
The general structure of the dataset assumes an underlying piecewise-linear temporal trend with a single changepoint (Figure C):
`f(t) = baseline + slope_1 * (t - t_min)` for `t <= changepoint`
`f(t) = baseline + slope_1 * (changepoint - t_min) + slope_2 * (t - changepoint)` for `t > changepoint`

Given the simulated nature of the dataset, the goal is to verify that the model recovers the known generating parameters:

- **Baseline** = 8
- **Slope 1** = 0.02
- **Slope 2** = -0.01
- **Changepoint** = 500 CE
- **Noise** ~ N(0, 1.5)

## Exploratory panel
<p align="center">
<img src="figures/exploratory_panel.png" height="700" text-align="center"/>
</p>
