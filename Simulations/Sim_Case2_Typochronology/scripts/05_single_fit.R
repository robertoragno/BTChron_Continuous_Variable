# Case 2 single-fit comparison figure: one dataset, midpoint vs full
# distribution, both fitted live. See shared/README.md and
# shared/scripts/single_fit_figure.R for what the figure shows and why.

library(here)
library(cmdstanr)

source(here("Simulations", "shared", "scripts", "weight_rows.R"))
source(here("Simulations", "shared", "scripts", "fit_models.R"))
source(here("Simulations", "shared", "scripts", "single_fit_figure.R"))
source(here("Simulations", "Sim_Case2_Typochronology", "scripts", "simulate.R"))

PERIOD_START <- 100; PERIOD_END <- 900
GRID_STEP <- 5
INTERCEPT <- 8; SLOPE <- 0.02; SIGMA <- 1.5
N <- 60

sim <- simulate_typo(N, INTERCEPT, SLOPE, SIGMA, prop_coarse_samples = 0.5,
                     period_start = PERIOD_START, period_end = PERIOD_END, seed = 7)

pad  <- 240   # longest coarse window at the reference, 1.2 * 25% of 800 yr
grid <- seq(PERIOD_START - pad, PERIOD_END + pad, by = GRID_STEP)
dates <- list(start = sim$Start_date, end = sim$End_date,
             median = (sim$Start_date + sim$End_date) / 2, grid = grid,
             weights = uniform_rows(sim$Start_date, sim$End_date, grid))

x_pred <- seq(PERIOD_START, PERIOD_END, length.out = 60)
truth  <- list(intercept = INTERCEPT, slope = SLOPE, sigma = SIGMA)

fits <- fit_single_example(dates, sim$Value, PERIOD_START - pad,
                           (PERIOD_END + pad) - (PERIOD_START - pad),
                           x_pred, seed = 7)

# What each model reads: the midpoint of each dating range, and the range itself
obs    <- data.frame(x = (sim$Start_date + sim$End_date) / 2, y = sim$Value)
ranges <- data.frame(xmin = sim$Start_date, xmax = sim$End_date, y = sim$Value)
out <- here("Simulations", "Sim_Case2_Typochronology", "figures",
           "single_fit_comparison.png")
single_fit_comparison_figure(
  fits, x_pred, truth, out,
  title = "Case 2: one dataset, half of the finds coarsely dated - midpoint vs full distribution",
  obs = obs, ranges = ranges)
cat("wrote", out, "\n")
