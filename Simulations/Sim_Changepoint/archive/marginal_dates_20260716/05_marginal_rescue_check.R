# Does marginalizing the latent dates fix the changepoint model's convergence
# failures? sim_changepoint.stan gives each observation its own free date
# parameter; a bone whose window straddles the changepoint has to pick a side
# (before/after the kink), and different chains can settle on different
# sides -- a multimodality visible as high Rhat / divergences / dropped fits
# in 02_recovery_study.R (about a fifth of fits). sim_changepoint_marginal.stan
# removes the per-observation parameters and integrates each date out over a
# fixed grid instead (quadrature), so there is no per-observation mode for
# chains to disagree about.
#
# This refits ONLY the marginal model, on the exact datasets (same seed, same
# true parameters) that the latent model already failed to converge on in the
# existing recovery study -- a direct rescue-rate check, not a full rerun.
#
# Archived 2026-07-16: experimental, not part of the main recovery/scenario
# pipeline (see Sim_Changepoint/README.md, "Marginalized latent dates"). Paths
# below point out of this archive folder back to the live simulation, so it
# can still be re-run in place if wanted; it was not re-run as part of the
# archiving itself.
#
# Produces: marginal_rescue_check.csv (in this archive folder)

library(here)
library(cmdstanr)
library(posterior)
library(dplyr)
library(parallel)

source(here("Simulations", "Sim_Changepoint", "scripts", "simulate.R"))

output_csv <- here("Simulations", "Sim_Changepoint", "archive",
                   "marginal_dates_20260716", "marginal_rescue_check.csv")
n_workers  <- as.integer(Sys.getenv("RESCUE_WORKERS", "18"))

recovery_results <- readr::read_csv(
  here("Simulations", "Sim_Changepoint", "output", "recovery_results.csv"),
  show_col_types = FALSE)

bad <- recovery_results %>%
  filter(model == "latent", is.na(error), (max_rhat > 1.05 | n_divergent > 0))

model_marginal <- cmdstan_model(here("Simulations", "Sim_Changepoint", "archive",
                                     "marginal_dates_20260716",
                                     "sim_changepoint_marginal.stan"))
prediction_grid <- c(TMIN, TMAX)

refit_one <- function(i) {
  job <- bad[i, ]
  one_dataset <- simulate_changepoint(
    N = job$N, baseline = job$true_baseline, slope_1 = job$true_slope_1,
    slope_2 = job$true_slope_2, cp = job$true_cp, sigma = job$true_sigma,
    K = job$K, alpha_conc = job$alpha_conc, seed = job$dataset_id)

  stan_data <- list(N = nrow(one_dataset), y = one_dataset$Value,
                    start_date = one_dataset$Start_date,
                    end_date   = one_dataset$End_date,
                    N_pred = length(prediction_grid), x_pred = prediction_grid)

  tryCatch({
    fit <- model_marginal$sample(
      data = stan_data, chains = 4, parallel_chains = 4,
      iter_warmup = 500, iter_sampling = 500,
      adapt_delta = 0.95, max_treedepth = 12,
      refresh = 0, show_messages = FALSE, show_exceptions = FALSE)
    diagnostics <- fit$diagnostic_summary(quiet = TRUE)
    rhats <- fit$summary(c("baseline_original", "slope1_original",
                           "slope2_original", "cp_actual", "sigma"))$rhat
    tibble(dataset_id = job$dataset_id, N = job$N, K = job$K, H = job$H,
           mean_width = job$mean_width,
           old_max_rhat = job$max_rhat, old_n_divergent = job$n_divergent,
           new_max_rhat = max(rhats, na.rm = TRUE),
           new_n_divergent = sum(diagnostics$num_divergent), error = NA_character_)
  }, error = function(e) tibble(dataset_id = job$dataset_id, N = job$N, K = job$K,
                                H = job$H, mean_width = job$mean_width,
                                old_max_rhat = job$max_rhat,
                                old_n_divergent = job$n_divergent,
                                new_max_rhat = NA_real_, new_n_divergent = NA_real_,
                                error = conditionMessage(e)))
}

if (file.exists(output_csv)) {
  cat("marginal_rescue_check.csv already exists; nothing to do.\n")
} else {
  cat(sprintf("Refitting the marginal model on all %d previously non-converged latent fits...\n",
              nrow(bad)))
  t0 <- Sys.time()
  results <- bind_rows(mclapply(seq_len(nrow(bad)), refit_one,
                                mc.cores = n_workers, mc.preschedule = FALSE))
  readr::write_csv(results, output_csv)
  cat(sprintf("Done in %.1f min. Wrote %s\n",
              as.numeric(difftime(Sys.time(), t0, units = "mins")), output_csv))

  rescued <- results %>% filter(is.na(error)) %>%
    mutate(rescued = new_max_rhat <= 1.05 & new_n_divergent == 0)
  cat(sprintf("Rescued (Rhat<=1.05, 0 divergences): %d/%d = %.1f%%\n",
              sum(rescued$rescued), nrow(rescued), 100 * mean(rescued$rescued)))
}
