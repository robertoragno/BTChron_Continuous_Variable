#' Purpose: Simulation-based calibration (SBC) check for both the latent-date
#'          and midpoint linear models. Re-draw true dates and noise K times,
#'          refit each model, and check credible interval coverage rates.
#'
#' Figure caption suggestion:
#'   "Posterior rank histograms from 100 SBC replications for the latent-date
#'    model (top) and midpoint model (bottom). Dashed line marks the expected
#'    bin count under uniform ranks. Coverage rates for the 50% and 90%
#'    credible intervals are annotated in each panel."

library(here)
library(tidyverse)
library(cmdstanr)
library(patchwork)

# ── Setup ────────────────────────────────────────────────────────────────────

sim_data <- read_csv(here("Simulations", "Sim_Linear", "data", "simulated_data.csv"),
                     show_col_types = FALSE)

true_baseline <- 5
true_slope    <- 0.015
true_sigma    <- 1.5
f_true        <- function(t) true_baseline + true_slope * t

n         <- nrow(sim_data)
pred_grid <- seq(min(sim_data$Start_date), max(sim_data$End_date), by = 1)

K <- 100

# ── Helper: extract quantiles and rank from a fit ────────────────────────────

extract_results <- function(fit, seed_id) {
  draws <- fit$draws(format = "df")
  bl <- draws$baseline_original
  sl <- draws$slope_original
  sg <- draws$sigma

  tibble(
    seed           = seed_id,
    baseline_q05   = quantile(bl, 0.05),
    baseline_q25   = quantile(bl, 0.25),
    baseline_q50   = median(bl),
    baseline_q75   = quantile(bl, 0.75),
    baseline_q95   = quantile(bl, 0.95),
    baseline_rank  = sum(bl < true_baseline),
    slope_q05      = quantile(sl, 0.05),
    slope_q25      = quantile(sl, 0.25),
    slope_q50      = median(sl),
    slope_q75      = quantile(sl, 0.75),
    slope_q95      = quantile(sl, 0.95),
    slope_rank     = sum(sl < true_slope),
    sigma_q05      = quantile(sg, 0.05),
    sigma_q25      = quantile(sg, 0.25),
    sigma_q50      = median(sg),
    sigma_q75      = quantile(sg, 0.75),
    sigma_q95      = quantile(sg, 0.95),
    sigma_rank     = sum(sg < true_sigma)
  )
}

# ── Load existing latent-date results (if available) or rerun ────────────────

latent_csv <- here("Simulations", "Sim_Linear", "output", "calibration_results.csv")

if (file.exists(latent_csv)) {
  cat("Loading existing latent-date SBC results.\n")
  results_latent <- read_csv(latent_csv, show_col_types = FALSE)
} else {
  cat("Running latent-date SBC loop...\n")
  model_latent <- cmdstan_model(here("Simulations", "Sim_Linear", "models", "sim_linear.stan"))
  results_latent <- vector("list", K)

  for (k in seq_len(K)) {
    set.seed(k)
    new_true_dates <- runif(n, min = sim_data$Start_date, max = sim_data$End_date)
    new_values     <- pmax(round(f_true(new_true_dates) + rnorm(n, 0, true_sigma), 1), 0.5)

    stan_data <- list(
      N = n, y = new_values,
      start_date = sim_data$Start_date, end_date = sim_data$End_date,
      N_pred = length(pred_grid), x_pred = pred_grid
    )

    fit <- model_latent$sample(
      data = stan_data, chains = 4, parallel_chains = 4,
      iter_warmup = 1000, iter_sampling = 1000, seed = 42,
      adapt_delta = 0.95, max_treedepth = 12,
      refresh = 0, show_messages = FALSE
    )

    results_latent[[k]] <- extract_results(fit, k)
    cat(sprintf("[%s] latent seed %3d / %d done\n", format(Sys.time(), "%H:%M:%S"), k, K))
  }
  results_latent <- bind_rows(results_latent)
  write_csv(results_latent, latent_csv)
}

# ── Midpoint SBC loop ───────────────────────────────────────────────────────

cat("Running midpoint SBC loop...\n")
model_midpoint <- cmdstan_model(here("Simulations", "Sim_Linear", "models",
                                     "sim_linear_midpoint.stan"))

results_midpoint <- vector("list", K)

for (k in seq_len(K)) {
  set.seed(k)
  new_true_dates <- runif(n, min = sim_data$Start_date, max = sim_data$End_date)
  new_values     <- pmax(round(f_true(new_true_dates) + rnorm(n, 0, true_sigma), 1), 0.5)

  stan_data <- list(
    N = n, y = new_values,
    start_date = sim_data$Start_date, end_date = sim_data$End_date,
    N_pred = length(pred_grid), x_pred = pred_grid
  )

  fit <- model_midpoint$sample(
    data = stan_data, chains = 4, parallel_chains = 4,
    iter_warmup = 1000, iter_sampling = 1000, seed = 42,
    adapt_delta = 0.95, max_treedepth = 12,
    refresh = 0, show_messages = FALSE
  )

  results_midpoint[[k]] <- extract_results(fit, k)
  cat(sprintf("[%s] midpoint seed %3d / %d done\n", format(Sys.time(), "%H:%M:%S"), k, K))
}

results_midpoint <- bind_rows(results_midpoint)
write_csv(results_midpoint,
          here("Simulations", "Sim_Linear", "output", "calibration_results_midpoint.csv"))

# ── Coverage summary ─────────────────────────────────────────────────────────

compute_coverage <- function(res) {
  tibble(
    Parameter = rep(c("Baseline", "Slope", "Sigma"), each = 2),
    CI        = rep(c("50%", "90%"), 3),
    Coverage  = c(
      mean(true_baseline >= res$baseline_q25 & true_baseline <= res$baseline_q75),
      mean(true_baseline >= res$baseline_q05 & true_baseline <= res$baseline_q95),
      mean(true_slope    >= res$slope_q25    & true_slope    <= res$slope_q75),
      mean(true_slope    >= res$slope_q05    & true_slope    <= res$slope_q95),
      mean(true_sigma    >= res$sigma_q25    & true_sigma    <= res$sigma_q75),
      mean(true_sigma    >= res$sigma_q05    & true_sigma    <= res$sigma_q95)
    )
  )
}

cov_latent   <- compute_coverage(results_latent)   %>% mutate(Model = "Latent-date")
cov_midpoint <- compute_coverage(results_midpoint) %>% mutate(Model = "Midpoint")

coverage <- bind_rows(cov_latent, cov_midpoint)

# ── Plotting ─────────────────────────────────────────────────────────────────

theme_panel <- theme_classic(base_size = 11) +
  theme(
    plot.title  = element_text(size = 11, face = "bold", hjust = 0),
    plot.margin = margin(5, 10, 5, 10)
  )

n_draws <- 4000
n_bins  <- 20
expected_count <- K / n_bins

rank_df <- bind_rows(
  tibble(Baseline = results_latent$baseline_rank,
         Slope    = results_latent$slope_rank,
         Sigma    = results_latent$sigma_rank,
         Model    = "Latent-date"),
  tibble(Baseline = results_midpoint$baseline_rank,
         Slope    = results_midpoint$slope_rank,
         Sigma    = results_midpoint$sigma_rank,
         Model    = "Midpoint")
) %>%
  pivot_longer(c(Baseline, Slope, Sigma),
               names_to = "Parameter", values_to = "Rank") %>%
  mutate(
    Parameter = factor(Parameter, levels = c("Baseline", "Slope", "Sigma")),
    Model     = factor(Model, levels = c("Latent-date", "Midpoint"))
  )

cov_labels <- coverage %>%
  mutate(label = sprintf("%s CI: %.0f%%", CI, Coverage * 100)) %>%
  group_by(Model, Parameter) %>%
  summarise(text = paste(label, collapse = "\n"), .groups = "drop") %>%
  mutate(
    Parameter = factor(Parameter, levels = c("Baseline", "Slope", "Sigma")),
    Model     = factor(Model, levels = c("Latent-date", "Midpoint"))
  )

p <- ggplot(rank_df, aes(x = Rank)) +
  geom_histogram(bins = n_bins, fill = "grey60", colour = "black", linewidth = 0.2) +
  geom_hline(yintercept = expected_count, linetype = "dashed", colour = "black") +
  geom_text(data = cov_labels, aes(label = text),
            x = n_draws * 0.95, y = Inf, hjust = 1, vjust = 1.5,
            size = 3, lineheight = 1.1) +
  facet_grid(Model ~ Parameter, scales = "free_x") +
  labs(x = "Rank of true value among posterior draws", y = "Count") +
  theme_panel +
  theme(
    strip.background = element_blank(),
    strip.text       = element_text(face = "bold", size = 10)
  )

ggsave(here("Simulations", "Sim_Linear", "figures", "calibration_coverage.png"),
       p, width = 10, height = 6, dpi = 300, bg = "white")
