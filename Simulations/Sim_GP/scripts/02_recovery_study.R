# Repeated calibration study for the GP simulation.
#
# The GP is non-parametric, so this study focuses on a known scalar target:
# sigma. We repeatedly redraw true dates and observation noise, refit both
# models, and save the posterior summaries needed for the recovery plots.

library(here)
library(tidyverse)
library(cmdstanr)

source(here("Simulations", "Sim_GP", "scripts", "simulate.R"))

sim_data <- read_csv(here("Simulations", "Sim_GP", "data", "simulated_data.csv"),
                     show_col_types = FALSE)

n <- nrow(sim_data)
pred_grid <- seq(min(sim_data$Start_date), max(sim_data$End_date), by = 1)
K <- as.integer(Sys.getenv("RECOVERY_K", "20"))
N_CHAINS <- as.integer(Sys.getenv("RECOVERY_CHAINS", "2"))
N_WARMUP <- as.integer(Sys.getenv("RECOVERY_WARMUP", "500"))
N_SAMPLING <- as.integer(Sys.getenv("RECOVERY_SAMPLING", "500"))
ADAPT_DELTA <- as.numeric(Sys.getenv("RECOVERY_ADAPT_DELTA", "0.99"))

extract_results <- function(fit, seed_id, runtime_seconds) {
  draws <- fit$draws(format = "df")
  sg <- draws$sigma

  tibble(
    seed = seed_id,
    runtime_seconds = runtime_seconds,
    sigma_q05 = quantile(sg, 0.05),
    sigma_q25 = quantile(sg, 0.25),
    sigma_q50 = median(sg),
    sigma_q75 = quantile(sg, 0.75),
    sigma_q95 = quantile(sg, 0.95),
    sigma_rank = sum(sg < TRUE_SIGMA)
  )
}

run_one_model <- function(model_path, output_csv, label) {
  if (file.exists(output_csv)) {
    cat(sprintf("Loading existing %s results.\n", label))
    return(read_csv(output_csv, show_col_types = FALSE))
  }

  model <- cmdstan_model(model_path)
  results <- vector("list", K)

  for (k in seq_len(K)) {
    set.seed(k)
    new_true_dates <- runif(n, min = sim_data$Start_date, max = sim_data$End_date)
    new_values <- round(f_true(new_true_dates) + rnorm(n, 0, TRUE_SIGMA), 1)

    stan_data <- list(
      N = n,
      y = new_values,
      start_date = sim_data$Start_date,
      end_date = sim_data$End_date,
      N_pred = length(pred_grid),
      x_pred = pred_grid,
      M = 20,
      c = 1.5
    )

    fit_start <- Sys.time()
    fit <- model$sample(
      data = stan_data,
      chains = N_CHAINS,
      parallel_chains = N_CHAINS,
      iter_warmup = N_WARMUP,
      iter_sampling = N_SAMPLING,
      seed = k,
      adapt_delta = ADAPT_DELTA,
      max_treedepth = 12,
      refresh = 0,
      show_messages = FALSE
    )
    fit_seconds <- as.numeric(difftime(Sys.time(), fit_start, units = "secs"))

    results[[k]] <- extract_results(fit, k, fit_seconds)
    cat(sprintf("[%s] %s seed %3d / %d done (%.1f min)\n",
                format(Sys.time(), "%H:%M:%S"), label, k, K, fit_seconds / 60))
  }

  results <- bind_rows(results)
  write_csv(results, output_csv)
  results
}

latent_csv <- here("Simulations", "Sim_GP", "output", "calibration_results.csv")
midpoint_csv <- here("Simulations", "Sim_GP", "output", "calibration_results_midpoint.csv")
runtime_csv <- here("Simulations", "Sim_GP", "output", "runtime_summary.csv")

latent_results <- run_one_model(
  here("Simulations", "Sim_GP", "models", "sim_hsgp.stan"),
  latent_csv,
  "latent"
)

midpoint_results <- run_one_model(
  here("Simulations", "Sim_GP", "models", "sim_hsgp_midpoint.stan"),
  midpoint_csv,
  "midpoint"
)

runtime_summary <- bind_rows(
  latent_results %>% mutate(model = "latent"),
  midpoint_results %>% mutate(model = "midpoint")
) %>%
  group_by(model) %>%
  summarise(
    fits = n(),
    total_minutes = sum(runtime_seconds) / 60,
    mean_minutes = mean(runtime_seconds) / 60,
    median_minutes = median(runtime_seconds) / 60,
    min_minutes = min(runtime_seconds) / 60,
    max_minutes = max(runtime_seconds) / 60,
    .groups = "drop"
  )

write_csv(runtime_summary, runtime_csv)

cat("Wrote calibration_results.csv, calibration_results_midpoint.csv, and runtime_summary.csv\n")
