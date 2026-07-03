# Recovery study for the GP model.
#
# Many random plausible datasets are simulated: each draws its own smooth wavy
# curve, noise sigma, sample size, and periodisation (a Dirichlet broken stick of
# K phases with concentration alpha_conc). Each dataset's dating resolution is
# summarised by the Shannon entropy H of its phase weights. Both the EIV
# (latent-date) and midpoint models are fit to each.
#
# The GP is non-parametric, so the curve itself is the recovery target: for each
# fit we record the proportion of the true curve caught by the 50% / 90% band
# over a time grid (accuracy) and the mean band width (precision). sigma is kept
# as a scalar target, as in the linear and changepoint studies, since it carries
# the "midpoint mistakes dating spread for noise" story and links the three sims.
#
# GP fits are heavy, so this is the expensive step. Run it in a tmux session.

library(here)
library(cmdstanr)
library(posterior)
library(dplyr)
library(tidyr)
library(parallel)

source(here("Simulations", "Sim_GP", "scripts", "simulate.R"))

set.seed(2026)

# Environment variables
output_csv  <- Sys.getenv("RECOVERY_OUT",
                          here("Simulations", "Sim_GP", "output", "recovery_results.csv"))
n_datasets    <- as.integer(Sys.getenv("RECOVERY_K", "60"))
n_workers     <- as.integer(Sys.getenv("RECOVERY_WORKERS", "16"))  # concurrent fits (each uses 4 chains)
iter_warmup   <- as.integer(Sys.getenv("RECOVERY_WARMUP", "1000"))
iter_sampling <- as.integer(Sys.getenv("RECOVERY_SAMPLING", "1000"))
M_basis       <- as.integer(Sys.getenv("RECOVERY_M", "20"))
c_boundary    <- as.numeric(Sys.getenv("RECOVERY_C", "1.5"))
# The latent-date GP has a hard posterior geometry, so it needs a high adapt_delta
# (as the earlier fixed-window GP study did) to keep divergences down.
adapt_delta   <- as.numeric(Sys.getenv("RECOVERY_ADAPT_DELTA", "0.99"))

# One parameter set per dataset. K and alpha_conc set the periodisation, matching
# the linear and changepoint studies; alpha_conc leans low (Beta over [0.1, 10])
# so lumpy phases are common. N sweeps two levels (low / high) to read precision
# against sample size. The curve (baseline, amplitudes, periods, phases) is drawn
# fresh per dataset by draw_curve and stored as a list column.
datasets <- tibble(
  dataset_id = seq_len(n_datasets),
  sigma      = runif(n_datasets, 0.5, 4),
  K          = sample(3:10, n_datasets, replace = TRUE),
  alpha_conc = 10^(-1 + 2 * rbeta(n_datasets, 2, 3.5)),
  N          = rep(c(50, 200), length.out = n_datasets),
  curve      = replicate(n_datasets, draw_curve(), simplify = FALSE)
)

# One job per dataset and model.
jobs <- expand_grid(dataset_id = datasets$dataset_id,
                    model = c("latent", "midpoint")) %>%
  left_join(datasets, by = "dataset_id") %>%
  mutate(job_id = row_number())

model_latent   <- cmdstan_model(here("Simulations", "Sim_GP", "models", "sim_hsgp.stan"))
model_midpoint <- cmdstan_model(here("Simulations", "Sim_GP", "models", "sim_hsgp_midpoint.stan"))

# Time grid on which the recovered curve is compared with the truth. Fixed across
# datasets so curve coverage is measured the same way everywhere.
prediction_grid <- seq(TMIN, TMAX, by = 10)

# Scalar-parameter summary (used for sigma): median, whether the 50%/90% interval
# held the truth, the 90% width, and the signed error.
summarise_scalar <- function(posterior_draws, true_value, name) {
  q <- quantile(posterior_draws, c(0.05, 0.25, 0.5, 0.75, 0.95), names = FALSE)
  tibble(
    !!paste0(name, "_med")     := q[3],
    !!paste0(name, "_cov50")   := true_value >= q[2] & true_value <= q[4],
    !!paste0(name, "_cov90")   := true_value >= q[1] & true_value <= q[5],
    !!paste0(name, "_width90") := q[5] - q[1],
    !!paste0(name, "_err")     := q[3] - true_value
  )
}

# Curve summary: over the grid, the proportion of the true curve inside the 50%
# and 90% bands (accuracy) and the mean band width (precision).
summarise_curve <- function(trend_matrix, true_curve) {
  lower90 <- apply(trend_matrix, 2, quantile, 0.05)
  upper90 <- apply(trend_matrix, 2, quantile, 0.95)
  lower50 <- apply(trend_matrix, 2, quantile, 0.25)
  upper50 <- apply(trend_matrix, 2, quantile, 0.75)
  tibble(
    curve_cov50   = mean(true_curve >= lower50 & true_curve <= upper50),
    curve_cov90   = mean(true_curve >= lower90 & true_curve <= upper90),
    curve_width50 = mean(upper50 - lower50),
    curve_width90 = mean(upper90 - lower90),
    curve_rmse    = sqrt(mean((apply(trend_matrix, 2, median) - true_curve)^2))
  )
}

# Simulate one dataset, fit one model, return a one-row summary.
run_one_job <- function(job_index) {
  job         <- jobs[job_index, ]
  curve       <- job$curve[[1]]
  one_dataset <- simulate_gp(N = job$N, sigma = job$sigma, K = job$K,
                             alpha_conc = job$alpha_conc, curve = curve,
                             seed = job$dataset_id)
  entropy_H  <- attr(one_dataset, "H")
  mean_width <- mean(one_dataset$End_date - one_dataset$Start_date)
  true_curve <- f_true(prediction_grid, curve)

  stan_data <- list(N = nrow(one_dataset), y = one_dataset$Value,
                    start_date = one_dataset$Start_date,
                    end_date = one_dataset$End_date,
                    N_pred = length(prediction_grid), x_pred = prediction_grid,
                    M = M_basis, c = c_boundary)

  chosen_model <- if (job$model == "latent") model_latent else model_midpoint

  summary_row <- tryCatch({
    fit <- chosen_model$sample(
      data = stan_data, chains = 4, parallel_chains = 4,
      iter_warmup = iter_warmup, iter_sampling = iter_sampling,
      adapt_delta = adapt_delta, max_treedepth = 12,
      refresh = 0, show_messages = FALSE, show_exceptions = FALSE
    )
    draws        <- fit$draws(format = "df")
    trend_matrix <- as.matrix(draws[, paste0("mu_pred[", seq_along(prediction_grid), "]")])
    diagnostics  <- fit$diagnostic_summary(quiet = TRUE)
    # Convergence is recorded per target so 03 can gate each metric on its own:
    # the curve accuracy/width/RMSE need the grid curve (mu_pred) to have mixed,
    # the sigma metric needs sigma to have mixed. The two can fail independently
    # (a wide-window, low-noise dataset leaves sigma weakly identified even when
    # the curve is fine, and vice versa). The latent dates are left out entirely:
    # they are often legitimately multimodal, so their Rhat should not gate anything.
    curve_rhat <- max(fit$summary("mu_pred")$rhat, na.rm = TRUE)
    sigma_rhat <- fit$summary("sigma")$rhat
    bind_cols(
      summarise_curve(trend_matrix, true_curve),
      summarise_scalar(draws$sigma, job$sigma, "sigma"),
      tibble(n_divergent = sum(diagnostics$num_divergent),
             curve_rhat  = curve_rhat, sigma_rhat = sigma_rhat,
             error = NA_character_)
    )
  }, error = function(e) tibble(error = conditionMessage(e)))

  bind_cols(
    tibble(dataset_id = job$dataset_id, model = job$model,
           true_sigma = job$sigma, K = job$K, alpha_conc = job$alpha_conc,
           N = job$N, H = entropy_H, mean_width = mean_width),
    summary_row
  )
}

if (file.exists(output_csv)) {
  cat("recovery_results.csv already exists; move it aside to re-run.\n")
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
    cat(sprintf("WARNING: %d/%d fits with errors\n", n_errored, nrow(results)))
}
