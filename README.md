# BTChron Continuous Variable

> [!NOTE]  
> Work in progress for the [first BTChron paper](https://www.overleaf.com/read/jjbchksqjpdr#694367). This readme will be updated soon with more details.

> [!ALERT]
> Note to self (RR): I need to clean up the code a bit, for instance I can move the partition function to a single helper. I also need to make some comments shorter and more elegant

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
