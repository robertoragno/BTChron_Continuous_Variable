#' Purpose: Fit the HSGP midpoint model (dates fixed at window centres).

library(here)
library(tidyverse)
library(cmdstanr)

sim_data <- read_csv(here("Simulations", "Sim_GP", "data", "simulated_data.csv"),
                     show_col_types = FALSE)

pred_grid <- seq(min(sim_data$Start_date), max(sim_data$End_date), by = 1)

stan_data <- list(
  N          = nrow(sim_data),
  y          = sim_data$Value,
  start_date = sim_data$Start_date,
  end_date   = sim_data$End_date,
  N_pred     = length(pred_grid),
  x_pred     = pred_grid,
  M          = 20,
  c          = 1.5
)

model <- cmdstan_model(here("Simulations", "Sim_GP", "models", "sim_hsgp_midpoint.stan"))

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

fit$save_object(here("Simulations", "Sim_GP", "output", "fit_hsgp_midpoint.rds"))

fit$cmdstan_diagnose()
fit$summary(variables = c("mu", "alpha", "rho", "sigma", "rho_actual"))
