# Case 3 single-fit comparison figure: one dataset, midpoint vs full
# distribution, both fitted live.

library(here)
library(cmdstanr)

source(here("Simulations", "shared", "scripts", "weight_rows.R"))
source(here("Simulations", "shared", "scripts", "fit_models.R"))
source(here("Simulations", "shared", "scripts", "single_fit_figure.R"))
source(here("Simulations", "Sim_Case3_OverlappingPhases", "scripts", "simulate.R"))

PERIOD_START <- 100; PERIOD_END <- 900
GRID_STEP <- 5
INTERCEPT <- 8; SLOPE <- 0.02; SIGMA <- 1
N <- 200
OVERLAP <- 0.5; ASSIGN_P <- 0.5

sim <- simulate_overlap(N, INTERCEPT, SLOPE, SIGMA, K = 6, alpha_conc = 1,
                        overlap = OVERLAP, assign_p = ASSIGN_P,
                        period_start = PERIOD_START, period_end = PERIOD_END, seed = 3)

pad <- 249   # same padding as 01_design.R
grid  <- seq(PERIOD_START - pad, PERIOD_END + pad, by = GRID_STEP)
dates <- list(start = sim$Start_date, end = sim$End_date,
              median = (sim$Start_date + sim$End_date) / 2, grid = grid,
              weights = uniform_rows(sim$Start_date, sim$End_date, grid))

x_pred <- seq(PERIOD_START, PERIOD_END, length.out = 60)
truth  <- list(intercept = INTERCEPT, slope = SLOPE, sigma = SIGMA)

fits <- fit_single_example(dates, sim$Value, PERIOD_START - pad,
                           (PERIOD_END + pad) - (PERIOD_START - pad),
                           x_pred, seed = 3)

# What each model reads: the midpoint of each dating range, and the range itself
obs    <- data.frame(x = (sim$Start_date + sim$End_date) / 2, y = sim$Value)
ranges <- data.frame(xmin = sim$Start_date, xmax = sim$End_date, y = sim$Value)
out <- here("Simulations", "Sim_Case3_OverlappingPhases", "figures", "single_fit_comparison.png")
single_fit_comparison_figure(
  fits, x_pred, truth, out,
  title = "Case 3: one dataset, phases overlapping by half - midpoint vs full distribution",
  obs = obs, ranges = ranges)
cat("wrote", out, "\n")
