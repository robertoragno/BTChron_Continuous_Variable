#' Purpose: Fit the changepoint Stan model with latent date inference
#'          and check recovery of the known generating parameters.

library(here)
library(tidyverse)
library(cmdstanr)

# Load data

sim_data    <- read_csv(here("Simulations", "Sim_Changepoint", "data",
                             "simulated_data.csv"),
                        show_col_types = FALSE)
true_params <- read_csv(here("Simulations", "Sim_Changepoint", "data",
                             "generating_parameters.csv"),
                        show_col_types = FALSE)

# Prediction grid

pred_grid <- seq(min(sim_data$Start_date), max(sim_data$End_date), by = 1)

# Prepare Stan data

stan_data <- list(
  N          = nrow(sim_data),
  y          = sim_data$Value,
  start_date = sim_data$Start_date,
  end_date   = sim_data$End_date,
  N_pred     = length(pred_grid),
  x_pred     = pred_grid
)

# Compile and fit

model <- cmdstan_model(
  here("Simulations", "Sim_Changepoint", "models", "sim_changepoint.stan")
)

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

# Save fit

fit$save_object(
  here("Simulations", "Sim_Changepoint", "output", "fit_changepoint.rds")
)

# Diagnostics

fit$cmdstan_diagnose()
fit$summary(variables = c("alpha", "beta1", "beta2", "cp_norm", "sigma",
                          "baseline_original", "slope1_original",
                          "slope2_original", "cp_actual"))

# Parameter recovery

cat("\n── Generating parameters ──\n")
print(true_params)

cat("\n── Recovered (original scale) ──\n")
fit$summary(variables = c("baseline_original", "slope1_original",
                          "slope2_original", "cp_actual", "sigma"))
