# Repeated calibration study for the GP simulation.
#
# The GP is non-parametric, so the repeated study tracks a few simple recovery
# targets: sigma plus curve features (peak year, peak value, width at half
# height). We repeatedly redraw true dates and observation noise, refit both
# models, and save the posterior summaries needed for the recovery plots.

library(here)
library(tidyverse)
library(cmdstanr)

source(here("Simulations", "Sim_GP", "scripts", "simulate.R"))

sim_data <- read_csv(here("Simulations", "Sim_GP", "data", "simulated_data.csv"),
                     show_col_types = FALSE)

n <- nrow(sim_data)
mean_width <- mean(sim_data$End_date - sim_data$Start_date)
pred_grid <- seq(min(sim_data$Start_date), max(sim_data$End_date), by = 1)
true_targets <- evaluation_targets() %>% select(target, value)
K <- as.integer(Sys.getenv("RECOVERY_K", "20"))
N_CHAINS <- as.integer(Sys.getenv("RECOVERY_CHAINS", "2"))
N_WARMUP <- as.integer(Sys.getenv("RECOVERY_WARMUP", "500"))
N_SAMPLING <- as.integer(Sys.getenv("RECOVERY_SAMPLING", "500"))
ADAPT_DELTA <- as.numeric(Sys.getenv("RECOVERY_ADAPT_DELTA", "0.99"))
output_dir <- Sys.getenv(
  "GP_RECOVERY_OUTPUT_DIR",
  here("Simulations", "Sim_GP", "output")
)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# These columns are the new curve-feature summaries. If an older cached CSV does
# not have them, we know it came from the earlier sigma-only study and must be
# recomputed before the new recovery plots can be trusted.
feature_columns <- c(
  "peak_year_q05", "peak_year_q25", "peak_year_q50", "peak_year_q75", "peak_year_q95",
  "peak_value_q05", "peak_value_q25", "peak_value_q50", "peak_value_q75", "peak_value_q95",
  "fwhm_years_q05", "fwhm_years_q25", "fwhm_years_q50", "fwhm_years_q75", "fwhm_years_q95"
)

extract_curve_features <- function(curve, years) {
  peak_idx <- which.max(curve)
  peak_year <- years[peak_idx]
  peak_value <- curve[peak_idx]
  # FWHM = full width at half maximum. For this single-peaked bell curve it is
  # a simple way to describe how spread out the recovered peak is.
  half_height <- min(curve) + 0.5 * (peak_value - min(curve))
  above_half <- which(curve >= half_height)
  fwhm_years <- years[max(above_half)] - years[min(above_half)]

  c(
    peak_year = peak_year,
    peak_value = peak_value,
    fwhm_years = fwhm_years
  )
}

read_existing_results <- function(path, expected_rows, label) {
  if (!file.exists(path)) {
    return(NULL)
  }

  existing <- read_csv(path, show_col_types = FALSE)
  if (nrow(existing) == expected_rows) {
    if (!"n_divergent" %in% names(existing)) {
      existing <- existing %>% mutate(n_divergent = 0)
    }
    if (!"max_rhat" %in% names(existing)) {
      existing <- existing %>% mutate(max_rhat = 1)
    }
    if (!"error" %in% names(existing)) {
      existing <- existing %>% mutate(error = NA_character_)
    }
    missing_features <- setdiff(feature_columns, names(existing))
    if (length(missing_features) > 0) {
      cat(sprintf(
        "Existing %s results are missing GP feature summaries (%s). Refitting.\n",
        label,
        paste(missing_features, collapse = ", ")
      ))
      return(NULL)
    }
    cat(sprintf("Loading existing %s results (%d rows).\n", label, nrow(existing)))
    return(existing)
  }

  cat(sprintf(
    "Existing %s results have %d row%s, but RECOVERY_K = %d. Refitting.\n",
    label,
    nrow(existing),
    ifelse(nrow(existing) == 1, "", "s"),
    expected_rows
  ))
  NULL
}

extract_results <- function(fit, seed_id, runtime_seconds) {
  draws <- fit$draws(format = "df")
  sg <- draws$sigma
  trend_cols <- grep("^mu_pred\\[", names(draws), value = TRUE)
  trend_mat <- as.matrix(draws[, trend_cols, drop = FALSE])
  # Each posterior draw gives one whole recovered curve on the prediction grid.
  # We turn each draw into a few interpretable summaries (peak year/value/FWHM),
  # then summarise those across draws exactly like we do for sigma.
  curve_features <- t(vapply(
    seq_len(nrow(trend_mat)),
    function(i) extract_curve_features(trend_mat[i, ], pred_grid),
    numeric(3)
  ))
  colnames(curve_features) <- c("peak_year", "peak_value", "fwhm_years")
  diagnostics <- fit$diagnostic_summary(quiet = TRUE)
  sigma_rhat <- fit$summary("sigma")$rhat
  sigma_q <- quantile(sg, c(0.05, 0.25, 0.5, 0.75, 0.95), names = FALSE)

  feature_q <- function(name) {
    quantile(curve_features[, name], c(0.05, 0.25, 0.5, 0.75, 0.95), names = FALSE)
  }

  peak_year_q <- feature_q("peak_year")
  peak_value_q <- feature_q("peak_value")
  fwhm_years_q <- feature_q("fwhm_years")

  tibble(
    seed = seed_id,
    runtime_seconds = runtime_seconds,
    sigma_q05 = sigma_q[1],
    sigma_q25 = sigma_q[2],
    sigma_q50 = sigma_q[3],
    sigma_q75 = sigma_q[4],
    sigma_q95 = sigma_q[5],
    peak_year_q05 = peak_year_q[1],
    peak_year_q25 = peak_year_q[2],
    peak_year_q50 = peak_year_q[3],
    peak_year_q75 = peak_year_q[4],
    peak_year_q95 = peak_year_q[5],
    peak_value_q05 = peak_value_q[1],
    peak_value_q25 = peak_value_q[2],
    peak_value_q50 = peak_value_q[3],
    peak_value_q75 = peak_value_q[4],
    peak_value_q95 = peak_value_q[5],
    fwhm_years_q05 = fwhm_years_q[1],
    fwhm_years_q25 = fwhm_years_q[2],
    fwhm_years_q50 = fwhm_years_q[3],
    fwhm_years_q75 = fwhm_years_q[4],
    fwhm_years_q95 = fwhm_years_q[5],
    sigma_rank = sum(sg < TRUE_SIGMA),
    n_divergent = sum(diagnostics$num_divergent),
    max_rhat = max(sigma_rhat, na.rm = TRUE),
    error = NA_character_
  )
}

run_one_model <- function(model_path, output_csv, label) {
  existing <- read_existing_results(output_csv, K, label)
  if (!is.null(existing)) {
    return(existing)
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

summarise_feature_recovery <- function(results, model_code) {
  make_target_rows <- function(target_name, q05, q25, q50, q75, q95) {
    truth <- true_targets$value[match(target_name, true_targets$target)]

    tibble(
      dataset_id = results$seed,
      model = model_code,
      target = target_name,
      true_value = truth,
      estimate_med = q50,
      cov50 = truth >= q25 & truth <= q75,
      cov90 = truth >= q05 & truth <= q95,
      width90 = q95 - q05,
      err = q50 - truth,
      runtime_seconds = results$runtime_seconds,
      mean_width = mean_width,
      N = n,
      n_divergent = results$n_divergent,
      max_rhat = results$max_rhat,
      error = results$error
    )
  }

  # Bind the four GP targets into one long table so the plotting script can use
  # the same "target / estimate / coverage / error" structure as the parametric
  # simulations.
  bind_rows(
    make_target_rows("peak_year",
                     results$peak_year_q05, results$peak_year_q25, results$peak_year_q50,
                     results$peak_year_q75, results$peak_year_q95),
    make_target_rows("peak_value",
                     results$peak_value_q05, results$peak_value_q25, results$peak_value_q50,
                     results$peak_value_q75, results$peak_value_q95),
    make_target_rows("fwhm_years",
                     results$fwhm_years_q05, results$fwhm_years_q25, results$fwhm_years_q50,
                     results$fwhm_years_q75, results$fwhm_years_q95),
    make_target_rows("sigma_noise",
                     results$sigma_q05, results$sigma_q25, results$sigma_q50,
                     results$sigma_q75, results$sigma_q95)
  )
}

latent_csv <- file.path(output_dir, "calibration_results.csv")
midpoint_csv <- file.path(output_dir, "calibration_results_midpoint.csv")
recovery_csv <- file.path(output_dir, "recovery_results.csv")
feature_csv <- file.path(output_dir, "feature_recovery_results.csv")
runtime_csv <- file.path(output_dir, "runtime_summary.csv")

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

recovery_results <- bind_rows(
  latent_results %>%
    transmute(
      dataset_id = seed,
      model = "latent",
      true_sigma = TRUE_SIGMA,
      mean_width = mean_width,
      N = n,
      sigma_med = sigma_q50,
      sigma_cov50 = TRUE_SIGMA >= sigma_q25 & TRUE_SIGMA <= sigma_q75,
      sigma_cov90 = TRUE_SIGMA >= sigma_q05 & TRUE_SIGMA <= sigma_q95,
      sigma_width90 = sigma_q95 - sigma_q05,
      sigma_err = sigma_q50 - TRUE_SIGMA,
      n_divergent,
      max_rhat,
      error
    ),
  midpoint_results %>%
    transmute(
      dataset_id = seed,
      model = "midpoint",
      true_sigma = TRUE_SIGMA,
      mean_width = mean_width,
      N = n,
      sigma_med = sigma_q50,
      sigma_cov50 = TRUE_SIGMA >= sigma_q25 & TRUE_SIGMA <= sigma_q75,
      sigma_cov90 = TRUE_SIGMA >= sigma_q05 & TRUE_SIGMA <= sigma_q95,
      sigma_width90 = sigma_q95 - sigma_q05,
      sigma_err = sigma_q50 - TRUE_SIGMA,
      n_divergent,
      max_rhat,
      error
    )
)

write_csv(recovery_results, recovery_csv)

feature_results <- bind_rows(
  summarise_feature_recovery(latent_results, "latent"),
  summarise_feature_recovery(midpoint_results, "midpoint")
)

write_csv(feature_results, feature_csv)

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

cat("Wrote calibration_results.csv, calibration_results_midpoint.csv, recovery_results.csv, feature_recovery_results.csv, and runtime_summary.csv\n")
