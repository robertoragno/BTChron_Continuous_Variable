#' Recovery study for the linear model.
#'
#' Many random plausible datasets are simulated: each draws its own intercept,
#' slope, noise sigma, sample size, and periodisation (a Dirichlet broken stick
#' of K phases with concentration alpha_conc). Each dataset's dating resolution
#' is summarised by the Shannon entropy H of its phase weights. Both the
#' latent-date and midpoint models are fit to each, and we record whether each
#' 50% / 90% credible interval contains the true value (accuracy) and how wide it
# is (precision). Nothing is held fixed except the structure of the question.

# Libraries
library(here)
library(cmdstanr)
library(posterior)
library(dplyr)
library(tidyr)
library(parallel)

# Source simulation helper
source(here("Simulations", "Sim_Linear", "scripts", "simulate.R"))

# Set seed
set.seed(2026)

# Environment variables
output_csv  <- Sys.getenv("RECOVERY_OUT",
                          here("Simulations", "Sim_Linear", "output", "recovery_results.csv"))
n_datasets  <- as.integer(Sys.getenv("RECOVERY_K", "400"))
n_workers   <- as.integer(Sys.getenv("RECOVERY_WORKERS", "18"))  # concurrent fits, change for less performing machines
iter_warmup   <- as.integer(Sys.getenv("RECOVERY_WARMUP", "500"))
iter_sampling <- as.integer(Sys.getenv("RECOVERY_SAMPLING", "500"))

# One parameter set per dataset. K and alpha_conc set the periodisation;
# alpha_conc leans low (Beta over [0.1, 10]) so lumpy phases are common. N sweeps
# four levels to read precision against sample size; everything else is random.

datasets <- tibble(
  dataset_id = seq_len(n_datasets),
  intercept  = runif(n_datasets, 2, 15),        # measurement level
  slope      = runif(n_datasets, -0.03, 0.03),  # trend up or down over ~800 yr
  sigma      = runif(n_datasets, 0.5, 4),       # low to high scatter
  K          = sample(3:10, n_datasets, replace = TRUE),          # number of phases
  #' Alpha still leans towards coarse phases, but by exponentiating 10 to a Beta, 
  #' we get a spread over [0.1, 10] rather than the conventional (Beta) boundaries 
  #' between 0 and 1.
  alpha_conc = 10^(-1 + 2 * rbeta(n_datasets, 2, 3.5)),           
  N          = rep(c(50, 100, 200, 400), length.out = n_datasets) # swept sample size
)

# One job per dataset and model
jobs <- expand_grid(dataset_id = datasets$dataset_id,
                    model = c("latent", "midpoint")) %>%
  left_join(datasets, by = "dataset_id") %>%
  mutate(job_id = row_number())

model_latent   <- cmdstan_model(here("Simulations", "Sim_Linear", "models",
                                     "sim_linear.stan"))
model_midpoint <- cmdstan_model(here("Simulations", "Sim_Linear", "models",
                                     "sim_linear_midpoint.stan"))

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

#' Helper function:
#' Simulate one dataset, fit one model, return a one-row summary.
run_one_job <- function(job_index) {
  job          <- jobs[job_index, ]
  one_dataset  <- simulate_linear(N = job$N, intercept = job$intercept,
                                  slope = job$slope, sigma = job$sigma,
                                  K = job$K, alpha_conc = job$alpha_conc,
                                  seed = job$dataset_id)
  entropy_H    <- attr(one_dataset, "H")
  mean_width   <- mean(one_dataset$End_date - one_dataset$Start_date)

  stan_data <- list(N = nrow(one_dataset), y = one_dataset$Value,
                    start_date = one_dataset$Start_date,
                    end_date = one_dataset$End_date,
                    N_pred = length(prediction_grid), x_pred = prediction_grid)

  chosen_model <- if (job$model == "latent") model_latent else model_midpoint

  summary_row <- tryCatch({
    fit <- chosen_model$sample(
      data = stan_data, chains = 4, parallel_chains = 4,
      iter_warmup = iter_warmup, iter_sampling = iter_sampling,
      adapt_delta = 0.95, max_treedepth = 10,
      refresh = 0, show_messages = FALSE, show_exceptions = FALSE
    )
    draws       <- fit$draws(variables = c("slope_original",
                                           "baseline_original", "sigma"),
                             format = "df")
    diagnostics <- fit$diagnostic_summary(quiet = TRUE)
    rhats       <- fit$summary(c("slope_original", "baseline_original",
                                 "sigma"))$rhat
    bind_cols(
      summarise_parameter(draws$baseline_original, job$intercept, "intercept"),
      summarise_parameter(draws$slope_original,    job$slope,     "slope"),
      summarise_parameter(draws$sigma,             job$sigma,     "sigma"),
      tibble(n_divergent = sum(diagnostics$num_divergent),
             max_rhat    = max(rhats, na.rm = TRUE), error = NA_character_)
    )
  }, error = function(e) tibble(error = conditionMessage(e)))

  bind_cols(
    job %>% select(dataset_id, model, true_intercept = intercept,
                   true_slope = slope, true_sigma = sigma, K, alpha_conc, N),
    tibble(H = entropy_H, mean_width = mean_width),
    summary_row
  )
}

if (file.exists(output_csv)) {
  cat("recovery_results.csv already exists; do not re run.\n")
} else {
  cat(sprintf("Running %d fits (%d datasets x 2 models) on %d workers.\n",
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
    cat(sprintf("WARNING: %d/%d fits with errors\n", n_errored, nrow(results)))
}
