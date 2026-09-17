# Case 4 single-fit comparison figure: one dataset, midpoint vs full
# distribution, both fitted live.

library(here)
library(cmdstanr)

source(here("Simulations", "shared", "scripts", "weight_rows.R"))
source(here("Simulations", "shared", "scripts", "fit_models.R"))
source(here("Simulations", "shared", "scripts", "single_fit_figure.R"))
source(here("Simulations", "Sim_Case4_MergedPhases", "scripts", "simulate.R"))

PERIOD_START <- 100; PERIOD_END <- 900
GRID_STEP <- 5
INTERCEPT <- 8; SLOPE <- 0.02; SIGMA <- 1
N <- 70; MERGE_MAX <- 3

sim <- simulate_merged(N, INTERCEPT, SLOPE, SIGMA, K = 8, alpha_conc = 1,
                       merge_max = MERGE_MAX,
                       period_start = PERIOD_START, period_end = PERIOD_END, seed = 4)

pad <- 0     # merged windows stay inside the period
grid  <- seq(PERIOD_START - pad, PERIOD_END + pad, by = GRID_STEP)
dates <- list(start = sim$Start_date, end = sim$End_date,
              median = (sim$Start_date + sim$End_date) / 2, grid = grid,
              weights = uniform_rows(sim$Start_date, sim$End_date, grid))

x_pred <- seq(PERIOD_START, PERIOD_END, length.out = 60)
truth  <- list(intercept = INTERCEPT, slope = SLOPE, sigma = SIGMA)

fits <- fit_single_example(dates, sim$Value, PERIOD_START - pad,
                           (PERIOD_END + pad) - (PERIOD_START - pad),
                           x_pred, seed = 4)

obs <- data.frame(x = (sim$Start_date + sim$End_date) / 2, y = sim$Value)
out <- here("Simulations", "Sim_Case4_MergedPhases", "figures", "single_fit_comparison.png")
single_fit_comparison_figure(fits, x_pred, truth, out,
                             title = "Case 4: one dataset, broad periods up to 3 phases", obs = obs)
cat("wrote", out, "\n")
