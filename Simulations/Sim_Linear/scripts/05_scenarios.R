# Scenario study for the linear model.
#
# The recovery study (02) draws fully random datasets and shows the general
# relationship. This script does the complementary thing: it fixes four
# archetype cases and draws many datasets within each, so the "when does EIV
# matter" question has clean, quotable numbers per case. Everything except the
# case-defining factor stays random (intercept, slope, sigma), so each case
# isolates one thing:
#
#   fine        finely dated: many even phases (high H). EIV and midpoint agree.
#   coarse      coarsely dated: few uneven phases (low H). EIV pulls ahead.
#   skewed      coarse dating with the widest phase late on the timeline and a
#               skewed deposition (dates pile early, Beta(1.5, 5)). Inside the
#               wide late phase the dates sit toward its early edge while the
#               midpoint stays at the centre: the one configuration where the
#               midpoint slope is systematically biased (04_dgp_sandbox.R,
#               skew_late). Small samples are not a case of their own; the
#               recovery study sweeps N.
#   diagnostic  mixed-resolution dating (some wide vague phases, some narrow),
#               with the narrow phases over-sampled (as when recognisable
#               ceramics date a short phase), so the sample lands in the well
#               dated phases and still spans the timeline. The EIV gap shrinks.
#
# Both models are fit to every dataset and the 50% / 90% interval coverage
# (accuracy) and width (precision) are recorded, as in 02.

library(here)
library(cmdstanr)
library(posterior)
library(dplyr)
library(tidyr)
library(parallel)

source(here("Simulations", "Sim_Linear", "scripts", "simulate.R"))

set.seed(2026)

output_csv <- Sys.getenv("SCENARIOS_OUT",
                         here("Simulations", "Sim_Linear", "output", "scenarios.csv"))
n_rep      <- as.integer(Sys.getenv("SCENARIOS_REPS", "100"))
n_workers  <- as.integer(Sys.getenv("SCENARIOS_WORKERS", "18"))
iter_warmup   <- as.integer(Sys.getenv("SCENARIOS_WARMUP", "500"))
iter_sampling <- as.integer(Sys.getenv("SCENARIOS_SAMPLING", "500"))

# The four cases turn three separate dials. fine and coarse differ only in the
# PERIODISATION (K, alpha_conc; dates stay uniform over the timeline, drawn
# inside partition_timeline via simulate_linear). diagnostic keeps uniform
# dates within each window but changes the SAMPLING (which phase a sample
# comes from, prob 1/width). skewed changes the DEPOSITION itself (dates
# Beta(1.5, 5) over the timeline, widest phase forced late).
scenarios <- tibble(
  case       = c("fine", "coarse", "skewed", "diagnostic"),
  K          = c(10,      3,        4,        10),
  alpha_conc = c(5,       0.3,      0.3,      0.7),
  N          = c(200,     200,      200,      200),
  sampling   = c("proportional", "proportional", "skewed", "diagnostic"))

# Diagnostic sampling: same Dirichlet broken-stick boundaries as
# partition_timeline, but observations are drawn towards the NARROW phases
# (probability proportional to 1 / phase width) and placed uniformly inside the
# window. Keeps the phase within-window date uniform, so the Stan prior still
# holds; only the sampling intensity across phases changes.
simulate_diagnostic <- function(N, intercept, slope, sigma, K, alpha_conc, seed) {
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

  out <- data.frame(
    Start_date = start, End_date = end, True_date = true_date,
    Value      = round(intercept + slope * true_date + rnorm(N, 0, sigma), 1))
  phase_p <- width / (TMAX - TMIN)
  attr(out, "H") <- -sum(phase_p * log(phase_p))
  out
}

# Skewed deposition: Dirichlet broken-stick boundaries with the phase widths
# sorted so the widest phase falls last, and true dates drawn from Beta(1.5, 5)
# over the timeline, then binned. The within-window date distribution in the
# wide late phase is tilted toward its early edge, which is what the Stan
# uniform prior and the midpoint both get wrong -- deliberately, since this
# case exists to measure the cost of that.
simulate_skewed <- function(N, intercept, slope, sigma, K, alpha_conc, seed) {
  set.seed(seed)
  weights <- sort(rgamma(K, shape = alpha_conc, rate = 1))   # widest phase last
  weights <- weights / sum(weights)
  spans   <- 1 + weights * ((TMAX - TMIN) - K)
  bounds  <- round(TMIN + c(0, cumsum(spans)))
  bounds[1] <- TMIN
  bounds[K + 1] <- TMAX

  true_date <- round(TMIN + rbeta(N, 1.5, 5) * (TMAX - TMIN))
  phase     <- findInterval(true_date, bounds, rightmost.closed = TRUE,
                            all.inside = TRUE)

  out <- data.frame(
    Start_date = bounds[phase], End_date = bounds[phase + 1],
    True_date  = true_date,
    Value      = round(intercept + slope * true_date + rnorm(N, 0, sigma), 1))
  width   <- diff(bounds)
  phase_p <- width / (TMAX - TMIN)
  attr(out, "H") <- -sum(phase_p * log(phase_p))
  out
}

# One dataset per case per replicate; nuisance parameters random, everything
# else set by the case.
datasets <- expand_grid(case = scenarios$case, rep = seq_len(n_rep)) %>%
  left_join(scenarios, by = "case") %>%
  mutate(dataset_id = row_number(),
         intercept  = runif(n(), 2, 15),
         slope      = runif(n(), -0.03, 0.03),
         sigma      = runif(n(), 0.5, 4))

jobs <- expand_grid(dataset_id = datasets$dataset_id,
                    model = c("latent", "midpoint")) %>%
  left_join(datasets, by = "dataset_id") %>%
  mutate(job_id = row_number())

model_latent   <- cmdstan_model(here("Simulations", "Sim_Linear", "models",
                                     "sim_linear.stan"))
model_midpoint <- cmdstan_model(here("Simulations", "Sim_Linear", "models",
                                     "sim_linear_midpoint.stan"))

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
    simulate_diagnostic(job$N, job$intercept, job$slope, job$sigma,
                        job$K, job$alpha_conc, seed = job$dataset_id)
  } else if (job$sampling == "skewed") {
    simulate_skewed(job$N, job$intercept, job$slope, job$sigma,
                    job$K, job$alpha_conc, seed = job$dataset_id)
  } else {
    simulate_linear(N = job$N, intercept = job$intercept, slope = job$slope,
                    sigma = job$sigma, K = job$K, alpha_conc = job$alpha_conc,
                    seed = job$dataset_id)
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
      adapt_delta = 0.95, max_treedepth = 10,
      refresh = 0, show_messages = FALSE, show_exceptions = FALSE)
    draws       <- fit$draws(variables = c("slope_original", "baseline_original",
                                           "sigma"), format = "df")
    diagnostics <- fit$diagnostic_summary(quiet = TRUE)
    rhats       <- fit$summary(c("slope_original", "baseline_original",
                                 "sigma"))$rhat
    bind_cols(
      summarise_parameter(draws$baseline_original, job$intercept, "intercept"),
      summarise_parameter(draws$slope_original,    job$slope,     "slope"),
      summarise_parameter(draws$sigma,             job$sigma,     "sigma"),
      tibble(n_divergent = sum(diagnostics$num_divergent),
             max_rhat    = max(rhats, na.rm = TRUE), error = NA_character_))
  }, error = function(e) tibble(error = conditionMessage(e)))

  bind_cols(
    job %>% select(dataset_id, case, model, true_intercept = intercept,
                   true_slope = slope, true_sigma = sigma, K, alpha_conc, N),
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

# One example dataset per case, with the nuisance parameters fixed at the
# worked-example values, so the four scenarios can be seen as datasets before
# any model touches them. Same look as the values-across-time panel of the
# exploratory figure: one line per sample at its value spanning its phase
# window, midpoint marked, true trend dashed. Runs whether or not the study
# above was skipped.
library(ggplot2)

example_intercept <- 5
example_slope     <- 0.015
example_sigma     <- 1.5
example_seed      <- 1

case_labs <- c(fine = "Fine dating", coarse = "Coarse dating",
               skewed = "Skewed deposition", diagnostic = "Diagnostic sampling")

examples <- bind_rows(lapply(seq_len(nrow(scenarios)), function(i) {
  sc <- scenarios[i, ]
  d <- if (sc$sampling == "diagnostic") {
    simulate_diagnostic(sc$N, example_intercept, example_slope, example_sigma,
                        sc$K, sc$alpha_conc, seed = example_seed)
  } else if (sc$sampling == "skewed") {
    simulate_skewed(sc$N, example_intercept, example_slope, example_sigma,
                    sc$K, sc$alpha_conc, seed = example_seed)
  } else {
    simulate_linear(N = sc$N, intercept = example_intercept,
                    slope = example_slope, sigma = example_sigma,
                    K = sc$K, alpha_conc = sc$alpha_conc, seed = example_seed)
  }
  d$case <- sc$case
  d
})) %>%
  mutate(case = factor(case, names(case_labs), case_labs),
         Midpoint = (Start_date + End_date) / 2)

truth_line <- data.frame(
  Year  = c(TMIN, TMAX),
  Value = example_intercept + example_slope * c(TMIN, TMAX))

# True dates shown as dots: in the skewed case they pile toward the early edge
# of the wide late phase while the midpoints stay at its centre, which is the
# feature that defines the case and the windows alone cannot show.
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

ggsave(here("Simulations", "Sim_Linear", "figures", "scenario_examples.png"),
       examples_plot, width = 9, height = 6.5, dpi = 300, bg = "white")
