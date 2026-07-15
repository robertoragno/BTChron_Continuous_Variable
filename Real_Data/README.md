# Real Data

Each `dataset_N` is one case study applying the methods validated in
`Simulations/` to a real archaeological dataset. Raw data and per-dataset fits
live in `data/dataset_N/` and `output/dataset_N/`; the shared `models/`,
`scripts/`, and `figures/` folders hold the code and outputs, with filenames
prefixed by the dataset.

## dataset_1 — Linear case study (GINI database)

The linear model from `Simulations/Sim_Linear`, applied to roofed-area data from
the GINI database (`gini_database_all_records_20240721.csv`). The response is
`log(RoofedArea)`; each house carries a dating window `[BeginHouse, EndHouse]`
rather than a known year. The question is the same as the simulation: does
treating each date as a latent parameter within its window change the inferred
trend, relative to collapsing the window to its midpoint.

Subset: Great Britain, nuclear-family houses, positive roofed area, dated window
(`EndHouse > BeginHouse`), midpoint year in `[-100, 1000]` — **2682 houses**.

> [!IMPORTANT]
>  Alternatively there is this paper, but it is missing a chronological range: https://www.frontiersin.org/journals/ecology-and-evolution/articles/10.3389/fevo.2023.1191274/full

### Scripts

| Script | Purpose | Output |
|---|---|---|
| `01_prepare.R` | filter the GINI database, write the subset, exploratory panel | `data/dataset_1/gini_gb_filtered.csv`, `figures/dataset_1_exploratory_panel.png` |
| `02_fit.R` | fit the latent-date and midpoint models, save fits and diagnostics | `output/dataset_1/fit_latent.rds`, `fit_midpoint.rds`, `diagnostics.txt` |
| `03_compare.R` | comparison figure and per-house date posteriors | `figures/dataset_1_model_comparison.png`, `dataset_1_individual_date_posteriors.png` |

Fitting is separated from plotting so figures can be restyled without refitting.
Both models converged cleanly (0 divergences, max Rhat 1.002, min bulk ESS > 3100).
Dates are mapped to `[-1, 1]` internally (not `[0, 1]`), which keeps alpha and
beta from being strongly correlated in the posterior; on this dataset it
roughly tripled the midpoint model's effective sample size (min bulk ESS 1253
-> 3816) without moving either model's estimates.

### Exploratory panel
<p align="center">
<img src="figures/dataset_1_exploratory_panel.png" height="320" text-align="center"/>
</p>

### Model comparison

The two models agree on the trend within the credible intervals. The latent-date
model returns a slightly steeper slope (0.00186 vs 0.00179) and a lower sigma
(0.963 vs 0.974) than the midpoint model, consistent with the simulation
finding: unmodelled date uncertainty is otherwise absorbed into the residual
variance. With no ground truth on real data, the figures carry no truth
overlay; the contrast between the two models is the point. Per the simulation
study, this gap is expected to be small under uniform-ish, reasonably fine
dating, and to widen under coarser or skewed dating -- it is not, on its own,
evidence that either model is "right" here.

<p align="center">
<img src="figures/dataset_1_model_comparison.png" height="800" text-align="center"/>
</p>

### Per-house date posteriors

For a sample of houses, the left panels show the recovered trend with the house's
dating window (shaded) and observed value (dashed); the right panels show the
posterior for each estimated date, with the window boundaries (dashed). The
posteriors stay close to the uniform window prior but tilt toward the end of the
window that agrees with the trend.

<p align="center">
<img src="figures/dataset_1_individual_date_posteriors.png" height="900" text-align="center"/>
</p>
