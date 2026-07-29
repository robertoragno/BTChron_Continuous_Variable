# BTChron Continuous Variable

> [!NOTE]  
> Work in progress for the [first BTChron paper](). This readme will be updated soon with more details.

> [!IMPORTANT]
> Note to self (RR): I need to clean up the code a bit, for instance I can move the partition function to a single helper. I also need to make some comments shorter and more elegant

## Repository structure

```
BTChron_Paper_1/
├── Simulations/
│   ├── Sim_Linear/        # Linear regression (baseline + slope)
│   ├── Sim_Changepoint/   # Changepoint regression (two slopes)
│   ├── archive/Sim_GP/    # Gaussian Process (bell-curve trend), probably remove later
│   └── (each contains: data/ figures/ scripts/ models/ output/)
├── Real_Data/
│   ├── dataset_1/          # Linear case study (GINI database)
│   ├── dataset_2/          # not yet started
│   ├── dataset_3/          # not yet started
│   └── (each contains: data/ figures/ scripts/ models/ output/ archive/)
```
