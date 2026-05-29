# Simulation 1: Linear Regression

This page shows the main findings from the simulated dataset of 300 samples of a continuous variable ranging approximately between ~5 and ~20 (Figure A). Each sample is assigned a time range with a `Start_Date` and an `End_Date`, with a known true date assigned between these boundaries (Figure B). 
The general structure of the dataset assumes an underlying linear temporal trend (Figure C):
`f(t) = baseline + slope * t`

Given the simulated nature of the dataset, the goal is to verify that the model recovers the known generating parameters:

- **Baseline** = 5
- **Slope** = 0.015
- **Noise** ~ N(0, 1.5)

## Exploratory panel
<p align="center">
<img src="figures/exploratory_panel.png" height="700" text-align="center"/>
</p>
