# Scenario study for the linear model.
#
# The recovery study (02) draws fully random datasets and shows the general
# relationship. This script does the complementary thing: it fixes archetype
# cases and draws many datasets within each, so the "when does EIV matter"
# question has clean, quotable numbers per case. Everything except the
# case-defining factor stays random (intercept, slope, sigma).
#
# There are three separate things about the dating that could go wrong, and
# only one of them can actually bias a linear slope:
#
#   How finely the timeline is periodised (phase width: K, alpha_conc) only
#   changes how uncertain each date is -- coarser phases give wider credible
#   intervals, but a linear fit under Berkson error is not biased by this.
#   Cases: fine (many even phases, high H), coarse (few uneven phases, low H).
#
#   Which phases the samples come from -- e.g. narrow, well-dated phases
#   being over-sampled because recognisable ceramics happen to date them --
#   also only affects precision (both OLS and Bayesian regression are
#   unbiased under an uneven mix of x-values). Case: diagnostic.
#
#   Where within a phase the true dates actually sit is the thing that can
#   matter. If deposition leans toward one edge of a phase -- e.g. population
#   growth concentrating activity late within it -- the midpoint mislabels
#   that sample's date, but on its own this mostly shifts the intercept
#   rather than bending the slope. Cases: fine/coarse crossed with
#   mild/strong skew, via simulate_skewed()'s skew_shape. skew_shape = 1 (the
#   fine/coarse cases above) is uniform deposition, i.e. no skew, so those
#   cases double as the "no skew" reference points for the skew comparison.
#
#   The slope only bends when phase width also correlates with position on
#   the timeline -- coarse phases clustering at one end. Whether that
#   happens, and which end, is a property of a specific dataset's evidentiary
#   record (which periods happen to be well-typologised), not something safe
#   to assume in general, so it is kept out of the main case grid above and
#   demonstrated separately, as a worked example, in the
#   coarse_strong_early / coarse_strong_late cases further down.
#
# Both models are fit to every dataset and the 50% / 90% interval coverage
# (accuracy), width (precision), and slope bias are recorded, as in 02.

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

# fine and coarse draw dates uniformly over the timeline as before
# (simulate_linear, via partition_timeline) and double as the "no skew"
# reference for the skew cases below. diagnostic keeps uniform dates within
# each window but over-samples the narrow phases. The four *_mild / *_strong
# cases reuse the fine/coarse phase widths but draw true dates from a
# within-phase Beta(skew_shape, 1) instead of uniform, via simulate_skewed(),
# with position_strength = 0 -- phase width uncorrelated with position, the
# general-purpose default with no assumption about which periods are better
# dated. coarse_strong_early / coarse_strong_late are not part of that
# general story: they turn position_strength on in each direction, purely as
# a worked example of what happens if width does correlate with position in
# a real dataset -- see simulate_skewed()'s comment below.
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

# Skewed deposition: same Dirichlet broken-stick widths and between-phase
# sampling (prob proportional to width) as the fine/coarse cases. Two things
# can change relative to those cases.
#
# First, where inside a phase the true date falls: instead of uniform, it is
# drawn from Beta(skew_shape, 1) scaled onto [start, end]. skew_shape = 1 is
# uniform (no skew); > 1 concentrates true dates toward the late edge of each
# phase, whatever the reason (activity intensity, taphonomy, deposition
# practice -- the mechanism doesn't matter here, only the resulting shape).
# This alone mostly shifts the intercept, not the slope.
#
# Second, and off by default (position_strength = 0), WHICH phase gets which
# width can be tied to position on the timeline: with position_strength > 0,
# wide phases are nudged toward one end (coarse_early = TRUE puts them early,
# FALSE puts them late) rather than placed at a random position. This is NOT
# a general archaeological rule -- whether coarse phases cluster early, late,
# or nowhere in particular is a property of a specific dataset's evidentiary
# record (which periods happen to have well-studied typologies, coinage,
# etc.), and can run in either direction depending on the case. It is used
# here only in the coarse_strong_early / coarse_strong_late cases, as a
# worked example of what happens if that correlation exists in a dataset --
# not as an assumption baked into the general-purpose cases above. The nudge
# is soft, not a rigid sort (position_strength blends a width-sorted order
# with a random one), so any one dataset can still come out as something
# other than a clean early-to-late staircase. What matters for the slope is
# the correlation between width and position averaged over many datasets, not
# the exact order in any single one: only when phase width is tied to
# position does that correlation appear and the skew bend the slope; with no
# such correlation, skew still shows up in sigma and the intercept, just not
# in the trend -- which is why the general-purpose cases leave it off.
simulate_skewed <- function(N, intercept, slope, sigma, K, alpha_conc,
                            skew_shape, seed, position_strength = 0,
                            coarse_early = TRUE) {
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

  out <- data.frame(
    Start_date = start, End_date = end, True_date = true_date,
    Value      = round(intercept + slope * true_date + rnorm(N, 0, sigma), 1))
  attr(out, "H") <- -sum(phase_p * log(phase_p))
  out
}

# One dataset per case per replicate; intercept and sigma are random nuisance
# parameters, everything else set by the case. Unlike the recovery study,
# slope is fixed rather than drawn near zero: the scenario table reports
# slope_ratio (estimate / true slope) as its bias measure, and dividing by a
# true slope close to zero makes that ratio blow up or flip sign for reasons
# that have nothing to do with which model is better. A fixed, moderate,
# single-sign slope (the same value as the worked example in 01_example.R)
# keeps that ratio meaningful.
datasets <- expand_grid(case = scenarios$case, rep = seq_len(n_rep)) %>%
  left_join(scenarios, by = "case") %>%
  mutate(dataset_id = row_number(),
         intercept  = runif(n(), 2, 15),
         slope      = 0.015,
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
                    job$K, job$alpha_conc, job$skew_shape, seed = job$dataset_id,
                    position_strength = job$position_strength,
                    coarse_early = job$coarse_early)
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

# fine_mild is left out of the gallery: on narrow phases even strong skew
# barely moves a date off its midpoint (see fine_strong), so the mild version
# adds no visible information and only breaks the grid's symmetry. It is
# still fit and reported in scenario_skew_bias.png, which needs all four
# skew cells for the width x strength design.
case_labs <- c(fine = "Fine dating", coarse = "Coarse dating",
               diagnostic = "Diagnostic sampling",
               fine_strong = "Fine, strong skew",
               coarse_mild = "Coarse, mild skew",
               coarse_strong = "Coarse, strong skew")

examples <- bind_rows(lapply(which(scenarios$case %in% names(case_labs)), function(i) {
  sc <- scenarios[i, ]
  d <- if (sc$sampling == "diagnostic") {
    simulate_diagnostic(sc$N, example_intercept, example_slope, example_sigma,
                        sc$K, sc$alpha_conc, seed = example_seed)
  } else if (sc$sampling == "skewed") {
    simulate_skewed(sc$N, example_intercept, example_slope, example_sigma,
                    sc$K, sc$alpha_conc, sc$skew_shape, seed = example_seed)
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
       examples_plot, width = 11, height = 6.5, dpi = 300, bg = "white")
