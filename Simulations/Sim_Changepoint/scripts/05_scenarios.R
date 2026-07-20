# Scenario study for the changepoint model: fixed archetype cases (vs. 02's
# fully random draws), many reps each, so "when does EIV matter" has clean
# numbers per case. Only within-phase deposition shape can bias the slopes
# or changepoint; phase width and sampling only affect precision. See
# Notes/05_scenarios_notes.md.

library(here)
library(cmdstanr)
library(posterior)
library(dplyr)
library(tidyr)
library(parallel)

source(here("Simulations", "Sim_Changepoint", "scripts", "simulate.R"))

set.seed(2026)

output_csv <- Sys.getenv("SCENARIOS_OUT",
                         here("Simulations", "Sim_Changepoint", "output", "scenarios.csv"))
n_rep      <- as.integer(Sys.getenv("SCENARIOS_REPS", "100"))
n_workers  <- as.integer(Sys.getenv("SCENARIOS_WORKERS", "18"))
iter_warmup   <- as.integer(Sys.getenv("SCENARIOS_WARMUP", "500"))
iter_sampling <- as.integer(Sys.getenv("SCENARIOS_SAMPLING", "500"))

# fine/coarse: uniform dates (no skew), double as skew cases' reference.
# diagnostic: uniform within-window, over-samples narrow phases.
# *_mild/*_strong: skewed (Beta(skew_shape,1)), position_strength = 0.
# coarse_strong_early/late: same, but position_strength > 0 (worked example
# of width correlating with position -- see notes).
scenarios <- tibble(
  case       = c("fine",   "coarse", "diagnostic",
                 "fine_mild", "fine_strong", "coarse_mild", "coarse_strong",
                 "coarse_strong_early", "coarse_strong_late"),
  K          = c(10,3,10,
                 10,10,3,3,
                 3,3),
  alpha_conc = c(5,0.3,0.7,
                 5,5,0.3,0.3,
                 0.3,0.3),
  N          = 200,
  sampling   = c("proportional", "proportional", "diagnostic",
                 "skewed", "skewed", "skewed", "skewed",
                 "skewed", "skewed"),
  skew_shape = c(NA,NA,NA,
                 2,4,2,4,
                 4,4),
  position_strength = c(NA,NA,NA,
                        0,0,0,0,
                        0.7,0.7),
  coarse_early = c(NA,NA,NA,
                   TRUE,TRUE,TRUE,TRUE,
                   TRUE,FALSE))

# Diagnostic sampling: over-samples narrow phases (prob ~ 1/width), uniform
# within-window (see notes).
simulate_diagnostic <- function(N, baseline, slope_1, slope_2, cp, sigma,
                                K, alpha_conc, seed) {
  set.seed(seed)
  weights <- rgamma(K, shape = alpha_conc, rate = 1)
  weights <- weights / sum(weights)
  spans   <- 1 + weights * ((TMAX - TMIN) - K)      # 1 year reserved per phase
  bounds  <- round(TMIN + c(0, cumsum(spans)))
  bounds[1] <- TMIN
  bounds[K + 1] <- TMAX
  width   <- diff(bounds)

  phase     <- sample(K, N, replace = TRUE, prob = 1 / width)
  start     <- bounds[phase]
  end       <- bounds[phase + 1]
  true_date <- round(runif(N, start, end))

  mu <- changepoint_mean(true_date, baseline, slope_1, slope_2, cp)
  out <- data.frame(
    Start_date = start, End_date = end, True_date = true_date,
    Value      = round(mu + rnorm(N, 0, sigma), 1))
  phase_p <- width / (TMAX - TMIN)
  attr(out, "H") <- -sum(phase_p * log(phase_p))
  out
}

# Skewed deposition: fine/coarse widths and sampling, but true dates drawn
# from Beta(skew_shape, 1) within each window instead of uniform, and
# optionally (position_strength > 0) with wide phases nudged toward one end
# (coarse_early). See notes for why only the latter bends the slopes/cp.
simulate_skewed <- function(N, baseline, slope_1, slope_2, cp, sigma,
                            K, alpha_conc, skew_shape, seed,
                            position_strength = 0, coarse_early = TRUE) {
  set.seed(seed)
  weights <- rgamma(K, shape = alpha_conc, rate = 1)
  weights <- weights / sum(weights)
  spans   <- 1 + weights * ((TMAX - TMIN) - K)

  order_key <- position_strength * rank(spans) +
    (1 - position_strength) * rank(runif(K))
  spans <- spans[order(order_key, decreasing = coarse_early)]

  bounds  <- round(TMIN + c(0, cumsum(spans)))
  bounds[1] <- TMIN
  bounds[K + 1] <- TMAX
  width   <- diff(bounds)
  phase_p <- width / (TMAX - TMIN)

  phase <- sample(K, N, replace = TRUE, prob = phase_p)
  start <- bounds[phase]
  end   <- bounds[phase + 1]

  u         <- rbeta(N, skew_shape, 1)
  true_date <- round(start + u * (end - start))

  mu <- changepoint_mean(true_date, baseline, slope_1, slope_2, cp)
  out <- data.frame(
    Start_date = start, End_date = end, True_date = true_date,
    Value      = round(mu + rnorm(N, 0, sigma), 1))
  attr(out, "H") <- -sum(phase_p * log(phase_p))
  out
}

# One dataset per case per replicate. Slopes/cp fixed (not drawn near zero
# like 02) so slope ratios stay meaningful -- see notes.
datasets <- expand_grid(case = scenarios$case, rep = seq_len(n_rep)) %>%
  left_join(scenarios, by = "case") %>%
  mutate(dataset_id = row_number(),
         baseline   = runif(n(), 2, 15),
         slope_1    = 0.015,
         slope_2    = -0.01,
         cp         = 500,
         sigma      = runif(n(), 0.5, 4))

jobs <- expand_grid(dataset_id = datasets$dataset_id,
                    model = c("latent", "midpoint")) %>%
  left_join(datasets, by = "dataset_id") %>%
  mutate(job_id = row_number())

model_latent   <- cmdstan_model(here("Simulations", "Sim_Changepoint", "models",
                                     "sim_changepoint.stan"))
model_midpoint <- cmdstan_model(here("Simulations", "Sim_Changepoint", "models",
                                     "sim_changepoint_midpoint.stan"))

prediction_grid <- c(TMIN, TMAX)

# Median, whether the 50/90 interval held the truth, the 90% width, signed error.
summarise_parameter <- function(posterior_draws, true_value, name) {
  quantiles <- quantile(posterior_draws, c(0.05, 0.25, 0.5, 0.75, 0.95),
                        names = FALSE)
  tibble(
    !!paste0(name, "_med")     := quantiles[3],
    !!paste0(name, "_cov50")   := true_value >= quantiles[2] & true_value <= quantiles[4],
    !!paste0(name, "_cov90")   := true_value >= quantiles[1] & true_value <= quantiles[5],
    !!paste0(name, "_width90") := quantiles[5] - quantiles[1],
    !!paste0(name, "_err")     := quantiles[3] - true_value)
}

# Simulate one dataset (by its case's sampling), fit one model, one-row summary.
run_one_job <- function(job_index) {
  job <- jobs[job_index, ]
  one_dataset <- if (job$sampling == "diagnostic") {
    simulate_diagnostic(job$N, job$baseline, job$slope_1, job$slope_2, job$cp,
                        job$sigma, job$K, job$alpha_conc, seed = job$dataset_id)
  } else if (job$sampling == "skewed") {
    simulate_skewed(job$N, job$baseline, job$slope_1, job$slope_2, job$cp,
                    job$sigma, job$K, job$alpha_conc, job$skew_shape,
                    seed = job$dataset_id, position_strength = job$position_strength,
                    coarse_early = job$coarse_early)
  } else {
    simulate_changepoint(N = job$N, baseline = job$baseline, slope_1 = job$slope_1,
                         slope_2 = job$slope_2, cp = job$cp, sigma = job$sigma,
                         K = job$K, alpha_conc = job$alpha_conc, seed = job$dataset_id)
  }
  entropy_H  <- attr(one_dataset, "H")
  mean_width <- mean(one_dataset$End_date - one_dataset$Start_date)

  stan_data <- list(N = nrow(one_dataset), y = one_dataset$Value,
                    start_date = one_dataset$Start_date,
                    end_date = one_dataset$End_date,
                    N_pred = length(prediction_grid), x_pred = prediction_grid)

  chosen_model <- if (job$model == "latent") model_latent else model_midpoint

  summary_row <- tryCatch({
    fit <- chosen_model$sample(
      data = stan_data, chains = 4, parallel_chains = 4,
      iter_warmup = iter_warmup, iter_sampling = iter_sampling,
      adapt_delta = 0.95, max_treedepth = 12,
      refresh = 0, show_messages = FALSE, show_exceptions = FALSE)
    draws       <- fit$draws(variables = c("baseline_original", "slope1_original",
                                           "slope2_original", "cp_actual", "sigma"),
                             format = "df")
    diagnostics <- fit$diagnostic_summary(quiet = TRUE)
    rhats       <- fit$summary(c("baseline_original", "slope1_original",
                                 "slope2_original", "cp_actual", "sigma"))$rhat
    bind_cols(
      summarise_parameter(draws$baseline_original, job$baseline, "baseline"),
      summarise_parameter(draws$slope1_original,    job$slope_1,  "slope1"),
      summarise_parameter(draws$slope2_original,    job$slope_2,  "slope2"),
      summarise_parameter(draws$cp_actual,           job$cp,       "cp"),
      summarise_parameter(draws$sigma,               job$sigma,    "sigma"),
      tibble(n_divergent = sum(diagnostics$num_divergent),
             max_rhat    = max(rhats, na.rm = TRUE), error = NA_character_))
  }, error = function(e) tibble(error = conditionMessage(e)))

  bind_cols(
    job %>% select(dataset_id, case, model, true_baseline = baseline,
                   true_slope_1 = slope_1, true_slope_2 = slope_2, true_cp = cp,
                   true_sigma = sigma, K, alpha_conc, N),
    tibble(H = entropy_H, mean_width = mean_width),
    summary_row)
}

if (file.exists(output_csv)) {
  cat("scenarios.csv already exists; do not re run.\n")
} else {
  cat(sprintf("Running %d fits (%d cases x %d reps x 2 models) on %d workers...\n",
              nrow(jobs), nrow(scenarios), n_rep, n_workers))
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

# One example dataset per case (fixed nuisance params) so scenarios can be
# seen before any model touches them. Runs whether or not the study above
# was skipped.
library(ggplot2)

example_baseline <- 8
example_slope_1  <- 0.015
example_slope_2  <- -0.01
example_cp       <- 500
example_sigma    <- 1.5
example_seed     <- 1

# fine_mild left out of the gallery (redundant with fine_strong visually);
# still fit/reported in scenario_skew_bias.png -- see notes.
case_labs <- c(fine = "Fine dating", coarse = "Coarse dating",
               diagnostic = "Diagnostic sampling",
               fine_strong = "Fine, strong skew",
               coarse_mild = "Coarse, mild skew",
               coarse_strong = "Coarse, strong skew")

examples <- bind_rows(lapply(which(scenarios$case %in% names(case_labs)), function(i) {
  sc <- scenarios[i, ]
  d <- if (sc$sampling == "diagnostic") {
    simulate_diagnostic(sc$N, example_baseline, example_slope_1, example_slope_2,
                        example_cp, example_sigma, sc$K, sc$alpha_conc,
                        seed = example_seed)
  } else if (sc$sampling == "skewed") {
    simulate_skewed(sc$N, example_baseline, example_slope_1, example_slope_2,
                    example_cp, example_sigma, sc$K, sc$alpha_conc, sc$skew_shape,
                    seed = example_seed)
  } else {
    simulate_changepoint(N = sc$N, baseline = example_baseline,
                         slope_1 = example_slope_1, slope_2 = example_slope_2,
                         cp = example_cp, sigma = example_sigma,
                         K = sc$K, alpha_conc = sc$alpha_conc, seed = example_seed)
  }
  d$case <- sc$case
  d
})) %>%
  mutate(case = factor(case, names(case_labs), case_labs),
         Midpoint = (Start_date + End_date) / 2)

truth_line <- data.frame(
  Year  = c(TMIN, example_cp, TMAX),
  Value = changepoint_mean(c(TMIN, example_cp, TMAX), example_baseline,
                           example_slope_1, example_slope_2, example_cp))

# True dates as dots vs. midpoint squares: shows the skew the windows alone
# can't.
examples_plot <- ggplot(examples) +
  geom_linerange(aes(y = Value, xmin = Start_date, xmax = End_date),
                 linewidth = 0.15, colour = "grey60", alpha = 0.5) +
  geom_point(aes(Midpoint, Value, colour = "Midpoint date"), shape = 15,
             size = 0.6) +
  geom_point(aes(True_date, Value, colour = "True date"), shape = 16,
             size = 0.5) +
  geom_line(data = truth_line, aes(Year, Value), linewidth = 0.6,
            colour = "black", linetype = "dashed") +
  facet_wrap(~ case, nrow = 2) +
  scale_colour_manual(NULL, values = c("Midpoint date" = "grey45",
                                       "True date" = "black")) +
  scale_x_continuous(breaks = seq(100, 900, 200), expand = c(0.01, 0)) +
  guides(colour = guide_legend(override.aes = list(size = 2))) +
  labs(title = "One simulated dataset per scenario",
       x = "Year (CE)", y = "Value") +
  theme_classic(base_size = 11) +
  theme(plot.title = element_text(size = 11, face = "bold", hjust = 0),
        strip.background = element_blank(),
        strip.text = element_text(face = "bold"),
        panel.spacing = unit(1.2, "lines"),
        legend.position = "bottom")

ggsave(here("Simulations", "Sim_Changepoint", "figures", "scenario_examples.png"),
       examples_plot, width = 11, height = 6.5, dpi = 300, bg = "white")
