# Checks marginal_date.stan against the latent-date model on Case 2 data.
#
# With uniform weight rows the two models encode the same assumption - the date is
# uniform in its window - by different routes, so they must agree on the trend
# parameters up to Monte Carlo error. Any disagreement means the weight rows,
# the grid indexing, or the shared normalisation is wrong, and this catches it
# before Case 1 puts a calibrated posterior in the same slot.

library(here)
library(cmdstanr)

source(here("Simulations", "Sim_Case2_Typochronology", "scripts", "simulate.R"))
source(here("Simulations", "shared", "scripts", "weight_rows.R"))

PERIOD_START <- 100
PERIOD_END   <- 900
GRID         <- 25
N            <- 200

# Candidate-year spacing. 1 is the reference; run again at 5 for the
# grid-resolution equivalence check, which is also the main cost lever, since
# runtime is linear in candidate years per row.
STEP <- as.integer(Sys.getenv("MARGINAL_GRID_STEP", "1"))

widths <- typo_widths(PERIOD_START, PERIOD_END, grid = GRID)
time_ref_min   <- PERIOD_START - max(widths)
time_ref_range <- (PERIOD_END + max(widths)) - time_ref_min

sim <- simulate_typo(N, intercept = 8, slope = 0.02, sigma = 2,
                     period_start = PERIOD_START, period_end = PERIOD_END,
                     widths = widths, grid = GRID, seed = 1)

grid <- make_grid(time_ref_min, time_ref_range, step = STEP)
w    <- uniform_rows(sim$Start_date, sim$End_date, grid)
cat(sprintf("grid step %d yr, %d candidate years, mean %.0f per row\n",
            STEP, length(grid), mean(w$row_n_years)))

x_pred <- seq(PERIOD_START, PERIOD_END, length.out = 50)

common <- list(N = N, y = sim$Value,
               time_ref_min = time_ref_min, time_ref_range = time_ref_range,
               N_pred = length(x_pred), x_pred = x_pred)

dat_marginal <- c(common,
                  list(n_years = length(grid), grid_year = grid,
                       n_weights = w$n_weights, log_year_prob_packed = w$log_year_prob_packed,
                       row_first_year = w$row_first_year, row_n_years = w$row_n_years))
dat_latent <- c(common,
                list(start_date = sim$Start_date, end_date = sim$End_date))

models <- here("Simulations", "shared", "models")
fit <- function(file, dat) {
  cmdstan_model(file.path(models, file))$sample(
    data = dat, chains = 4, parallel_chains = 4,
    iter_warmup = 1000, iter_sampling = 1000, seed = 1, refresh = 0)
}

f_marg <- fit("marginal_date.stan", dat_marginal)
f_lat  <- fit("latent_date.stan",   dat_latent)

pars <- c("slope_original", "baseline_original", "sigma")
cat("\n-- marginal --\n"); print(f_marg$summary(pars))
cat("\n-- latent --\n");   print(f_lat$summary(pars))
cat("\ntrue slope 0.02, intercept 8, sigma 2\n")

# Same assumption, two routes: the slope posteriors should overlap closely.
d_marg <- f_marg$draws("slope_original", format = "draws_matrix")[, 1]
d_lat  <- f_lat$draws("slope_original",  format = "draws_matrix")[, 1]
cat(sprintf("slope median  marginal %.5f  latent %.5f  diff %.5f (mcse ~%.5f)\n",
            median(d_marg), median(d_lat), median(d_marg) - median(d_lat),
            sd(d_marg) / sqrt(200)))

# Recovered dates should also match, since both target the same date posterior.
dm <- colMeans(f_marg$draws("date_actual",  format = "draws_matrix"))
dl <- colMeans(f_lat$draws("date_actual",   format = "draws_matrix"))
cat(sprintf("per-date posterior mean: max abs diff %.2f yr, correlation %.4f\n",
            max(abs(dm - dl)), cor(dm, dl)))
