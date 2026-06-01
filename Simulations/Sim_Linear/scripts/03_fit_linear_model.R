#' Purpose: Fit the simple linear Stan model to the simulated data
#'          and check recovery of the known generating parameters.

library(here)
library(tidyverse)
library(cmdstanr)

# Load data

sim_data     <- read_csv(here("Simulations", "Sim_Linear", "data", "simulated_data.csv"),
                         show_col_types = FALSE)
ground_truth <- read_csv(here("Simulations", "Sim_Linear", "data", "ground_truth.csv"),
                         show_col_types = FALSE)
true_params  <- read_csv(here("Simulations", "Sim_Linear", "data", "generating_parameters.csv"),
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

model <- cmdstan_model(here("Simulations", "Sim_Linear", "models", "sim_linear.stan"))

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

fit$save_object(here("Simulations", "Sim_Linear", "output", "fit_linear.rds"))

# Diagnostics

fit$cmdstan_diagnose()
fit$summary(variables = c("alpha", "beta", "sigma",
                          "slope_original", "baseline_original"))

# Parameter recovery

cat("\n── Generating parameters ──\n")
print(true_params)

cat("\n── Recovered (original scale) ──\n")
fit$summary(variables = c("baseline_original", "slope_original", "sigma"))
