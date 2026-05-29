# BTChron Continuous Variable

> [!NOTE]  
> Work in progress for the [first BTChron paper](https://www.overleaf.com/read/jjbchksqjpdr#694367). This readme will be updated soon with more details.

## Repository structure

```
BTChron_Paper_1/
├── Simulations/
│   ├── Sim_Linear/        # Linear regression (baseline + slope)
│   ├── Sim_Changepoint/   # Changepoint regression (two slopes)
│   └── Sim_GP/            # Gaussian Process (bell-curve trend)
│       └── (each contains: data/ figures/ scripts/ models/ output/)
├── Real_Data/
│   ├── data/
│   │   ├── dataset_1/
│   │   ├── dataset_2/
│   │   └── dataset_3/
│   ├── models/            # Stan models (real data)
│   ├── scripts/           # R scripts (real data)
│   ├── output/
│   │   ├── dataset_1/
│   │   ├── dataset_2/
│   │   └── dataset_3/
│   └── figures/
```
