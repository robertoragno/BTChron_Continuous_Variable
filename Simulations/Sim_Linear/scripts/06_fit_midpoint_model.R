#' Purpose: Fit the midpoint linear model (no latent date inference)

library(here)
library(tidyverse)
library(cmdstanr)

# Load data

sim_data <- read_csv(here("Simulations", "Sim_Linear", "data", "simulated_data.csv"),
                     show_col_types = FALSE)

pred_grid <- seq(min(sim_data$Start_date), max(sim_data$End_date), by = 1)

stan_data <- list(
  N          = nrow(sim_data),
  y          = sim_data$Value,
  start_date = sim_data$Start_date,
  end_date   = sim_data$End_date,
  N_pred     = length(pred_grid),
  x_pred     = pred_grid
)

# Compile and fit

model <- cmdstan_model(here("Simulations", "Sim_Linear", "models",
                            "sim_linear_midpoint.stan"))

fit <- model$sample(
  data            = stan_data,
  chains          = 4,
  parallel_chains = 4,
  iter_warmup     = 1000,
  iter_sampling   = 1000,
  seed            = 42,
  adapt_delta     = 0.95,
  max_treedepth   = 12
)

fit$save_object(here("Simulations", "Sim_Linear", "output",
                     "fit_midpoint.rds"))

fit$cmdstan_diagnose()

cat("\n── Midpoint model parameter recovery ──\n")
print(fit$summary(variables = c("baseline_original", "slope_original", "sigma")))
