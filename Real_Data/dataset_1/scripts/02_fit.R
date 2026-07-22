# Fit both linear models to the GB subset: latent-date inference vs fixed
# midpoint dates. Slow (one latent date per house, N ~ 2700), so fits are saved
# to output/dataset_1/ and the comparison figures are built separately in
# 03_compare.R. Run this once, in a tmux BTChron session.
#
# Produces: fit_latent.rds, fit_midpoint.rds, diagnostics.txt

library(here)
library(tidyverse)
library(cmdstanr)
library(posterior)

clean_path <- here("Real_Data", "dataset_1", "data", "gini_gb_filtered.csv")
output_dir <- here("Real_Data", "dataset_1", "output")

gb_data <- read_csv(clean_path, show_col_types = FALSE)
n_observations <- nrow(gb_data)

prediction_grid <- seq(min(gb_data$Start_date), max(gb_data$End_date), by = 1)
stan_data <- list(N = n_observations, y = gb_data$Value,
                  start_date = gb_data$Start_date, end_date = gb_data$End_date,
                  N_pred = length(prediction_grid), x_pred = prediction_grid)

fit_args <- list(data = stan_data, chains = 4, parallel_chains = 4,
                 iter_warmup = 1000, iter_sampling = 1000, seed = 42,
                 adapt_delta = 0.95, max_treedepth = 12, refresh = 100)

cat("Fitting latent-date model (N =", n_observations, ")\n")
fit_latent <- do.call(
  cmdstan_model(here("Real_Data", "dataset_1", "models", "linear_latent.stan"))$sample, fit_args)

cat("Fitting midpoint model\n")
fit_midpoint <- do.call(
  cmdstan_model(here("Real_Data", "dataset_1", "models", "linear_midpoint.stan"))$sample, fit_args)

fit_latent$save_object(file.path(output_dir, "fit_latent.rds"))
fit_midpoint$save_object(file.path(output_dir, "fit_midpoint.rds"))

# Diagnostics: gate the figures on these before trusting any comparison ---------

report_diagnostics <- function(fit, label) {
  pars <- c("baseline_original", "slope_original", "sigma")
  summ <- fit$summary(pars)
  sampler <- fit$diagnostic_summary()
  c(sprintf("--- %s ---", label),
    capture.output(print(summ)),
    sprintf("divergences (per chain): %s", paste(sampler$num_divergent, collapse = " ")),
    sprintf("max treedepth hits:      %s", paste(sampler$num_max_treedepth, collapse = " ")),
    sprintf("min bulk ESS: %.0f | min tail ESS: %.0f | max Rhat: %.4f",
            min(summ$ess_bulk), min(summ$ess_tail), max(summ$rhat)),
    "")
}

diagnostics <- c(report_diagnostics(fit_latent, "Latent dates"),
                 report_diagnostics(fit_midpoint, "Midpoint dates"))
writeLines(diagnostics, file.path(output_dir, "diagnostics.txt"))
cat(diagnostics, sep = "\n")
