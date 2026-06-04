# Recovery study for the changepoint model.
#
# Many random plausible datasets are simulated: each draws its own baseline,
# two slopes, changepoint location, noise sigma, sample size, and (mixed)
# dating-window widths. Both the latent-date and midpoint models are fit to
# each, and we record whether each 50% / 90% credible interval contains the
# true value (recovery) and how wide it is (precision). Nothing is held fixed
# except the structure of the question.
#
# Slope scheme: slope_1 ~ U(-0.03, 0.03); then slope_2 = slope_1 + delta
# where delta = ±U(0.01, 0.05). This guarantees |delta_slope| >= 0.01, so
# every dataset has a real kink and the changepoint is identifiable. Drawing
# both slopes independently would leave ~30% of datasets with no detectable
# changepoint, making cp unidentifiable and washing out the latent vs midpoint
# contrast.

library(here)
library(cmdstanr)
library(posterior)
library(dplyr)
library(tidyr)
library(parallel)

source(here("Simulations", "Sim_Changepoint", "scripts", "simulate.R"))

set.seed(2026)

output_csv <- Sys.getenv("RECOVERY_OUT",
  here("Simulations", "Sim_Changepoint", "output", "recovery_results.csv"))
n_datasets <- as.integer(Sys.getenv("RECOVERY_K", "400"))
n_workers  <- as.integer(Sys.getenv("RECOVERY_WORKERS", "18"))

# Draw one plausible parameter set per dataset

delta_sign <- sample(c(-1, 1), n_datasets, replace = TRUE)

datasets <- tibble(
  dataset_id   = seq_len(n_datasets),
  baseline     = runif(n_datasets, 2, 15),
  slope_1      = runif(n_datasets, -0.03, 0.03),
  delta        = delta_sign * runif(n_datasets, 0.01, 0.05),
  slope_2      = slope_1 + delta,
  cp           = runif(n_datasets, 300, 700),
  sigma        = runif(n_datasets, 0.5, 4),
  mean_width   = runif(n_datasets, 20, 450),
  N            = sample(100:500, n_datasets, replace = TRUE)
) %>% select(-delta)

jobs <- expand_grid(dataset_id = datasets$dataset_id,
                    model = c("latent", "midpoint")) %>%
  left_join(datasets, by = "dataset_id") %>%
  mutate(job_id = row_number())

model_latent   <- cmdstan_model(here("Simulations", "Sim_Changepoint", "models",
                                     "sim_changepoint.stan"))
model_midpoint <- cmdstan_model(here("Simulations", "Sim_Changepoint", "models",
                                     "sim_changepoint_midpoint.stan"))

prediction_grid <- c(TMIN, TMAX)

# Summarise one parameter's posterior: median, whether the 50%/90% interval
# contained the truth, the 90% width, and the signed error.
summarise_parameter <- function(posterior_draws, true_value, name) {
  quantiles <- quantile(posterior_draws, c(0.05, 0.25, 0.5, 0.75, 0.95),
                        names = FALSE)
  tibble(
    !!paste0(name, "_med")     := quantiles[3],
    !!paste0(name, "_cov50")   := true_value >= quantiles[2] & true_value <= quantiles[4],
    !!paste0(name, "_cov90")   := true_value >= quantiles[1] & true_value <= quantiles[5],
    !!paste0(name, "_width90") := quantiles[5] - quantiles[1],
    !!paste0(name, "_err")     := quantiles[3] - true_value
  )
}

# Simulate one dataset, fit one model, return a one-row summary.
run_one_job <- function(job_index) {
  job         <- jobs[job_index, ]
  one_dataset <- simulate_changepoint(
    N = job$N, baseline = job$baseline, slope_1 = job$slope_1,
    slope_2 = job$slope_2, cp = job$cp, sigma = job$sigma,
    mean_width = job$mean_width, seed = job$dataset_id
  )

  stan_data <- list(N = nrow(one_dataset), y = one_dataset$Value,
                    start_date = one_dataset$Start_date,
                    end_date   = one_dataset$End_date,
                    N_pred = length(prediction_grid), x_pred = prediction_grid)

  chosen_model <- if (job$model == "latent") model_latent else model_midpoint

  summary_row <- tryCatch({
    fit <- chosen_model$sample(
      data = stan_data, chains = 4, parallel_chains = 4,
      iter_warmup = 500, iter_sampling = 500,
      adapt_delta = 0.95, max_treedepth = 12,
      refresh = 0, show_messages = FALSE, show_exceptions = FALSE
    )
    draws       <- fit$draws(variables = c("baseline_original", "slope1_original",
                                           "slope2_original", "cp_actual", "sigma"),
                             format = "df")
    diagnostics <- fit$diagnostic_summary(quiet = TRUE)
    rhats       <- fit$summary(c("baseline_original", "slope1_original",
                                 "slope2_original", "cp_actual", "sigma"))$rhat
    bind_cols(
      summarise_parameter(draws$baseline_original, job$baseline,  "baseline"),
      summarise_parameter(draws$slope1_original,   job$slope_1,   "slope1"),
      summarise_parameter(draws$slope2_original,   job$slope_2,   "slope2"),
      summarise_parameter(draws$cp_actual,         job$cp,        "cp"),
      summarise_parameter(draws$sigma,             job$sigma,     "sigma"),
      tibble(n_divergent = sum(diagnostics$num_divergent),
             max_rhat    = max(rhats, na.rm = TRUE), error = NA_character_)
    )
  }, error = function(e) tibble(error = conditionMessage(e)))

  bind_cols(
    job %>% select(dataset_id, model,
                   true_baseline = baseline, true_slope_1 = slope_1,
                   true_slope_2 = slope_2, true_cp = cp, true_sigma = sigma,
                   mean_width, N),
    summary_row
  )
}

if (file.exists(output_csv)) {
  cat("recovery_results.csv already exists; nothing to do.\n")
} else {
  cat(sprintf("Running %d fits (%d datasets x 2 models) on %d workers...\n",
              nrow(jobs), n_datasets, n_workers))
  start_time <- Sys.time()
  results <- bind_rows(mclapply(seq_len(nrow(jobs)), run_one_job,
                                mc.cores = n_workers, mc.preschedule = FALSE))
  readr::write_csv(results, output_csv)
  cat(sprintf("Done in %.1f min. Wrote %s\n",
              as.numeric(difftime(Sys.time(), start_time, units = "mins")),
              output_csv))
  n_errored <- sum(!is.na(results$error))
  if (n_errored > 0)
    cat(sprintf("WARNING: %d/%d fits errored.\n", n_errored, nrow(results)))
}
