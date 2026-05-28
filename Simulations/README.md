# Simulations

This page shows the main findings from the simulated dataset of 300 samples of a continuous variable ranging approximately between (approximately) ~3 and ~25 (Figure A). Each sample is assigned a time range with a `Start_Date` and an `End_Date`, with a known true date assigned between these boundaries (Figure B). 
The general structure of the dataset assumes an underlying temporal trend, with a peak around 450 CE and a decline thereafter (Figure C):
`f(t) = baseline + amplitude * exp(-((t - peak) / width)^2)`

Given the simulated nature of the dataset, we want to make sure that our models are able to recover the 'true' trend. In this sense, we assign:

- **Baseline** = 8
- **Amplitude** = 12
- **Peak** = 450 CE
- **Width** = 200
- **Noise** ~ N(0, 2.5)

## Exploratory panel
<p align="center">
<img src="figures/exploratory_panel.png" height="700" text-align="center"/>
</p>
