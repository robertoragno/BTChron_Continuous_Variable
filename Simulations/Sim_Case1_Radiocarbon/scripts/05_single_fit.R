# Case 1 single-fit figure: one dataset on the Hallstatt plateau, fitted by the
# midpoint model and the full-distribution model.

library(here)
library(cmdstanr)
suppressMessages(library(rcarbon))

source(here("Simulations", "shared", "scripts", "weight_rows.R"))
source(here("Simulations", "shared", "scripts", "fit_models.R"))
source(here("Simulations", "shared", "scripts", "single_fit_figure.R"))
source(here("Simulations", "Sim_Case1_Radiocarbon", "scripts", "simulate.R"))

WINDOW    <- c(-800, -400)   # the Hallstatt plateau, as in 02_figures.R
LAB_ERROR <- 25
INTERCEPT <- 8; SLOPE <- 0.02; SIGMA <- 2
N <- 50; STEP <- 5

pad  <- c14_grid_pad(WINDOW, LAB_ERROR)
grid <- seq(WINDOW[1] - pad, WINDOW[2] + pad, by = STEP)
sim  <- simulate_c14(N, INTERCEPT, SLOPE, SIGMA, WINDOW, LAB_ERROR, seed = 1)

cal  <- calmatrix_to_calendar(calibrate(sim$CRA, errors = sim$Error,
                                        calCurves = "intcal20", calMatrix = TRUE,
                                        verbose = FALSE))
smry <- calibrated_summaries(cal, grid)
dates <- list(start = smry$start, end = smry$end, median = smry$median,
             grid = grid, weights = calibrated_rows(cal, grid))

x_pred <- seq(WINDOW[1], WINDOW[2], length.out = 60)
truth  <- list(intercept = INTERCEPT, slope = SLOPE, sigma = SIGMA)

fits <- fit_single_example(dates, sim$Value, min(grid), max(grid) - min(grid),
                           x_pred, seed = 1)

# Points sit where the midpoint model sees them: the centre of the 95% range
obs <- data.frame(x = (smry$start + smry$end) / 2, y = sim$Value)

out <- here("Simulations", "Sim_Case1_Radiocarbon", "figures",
           "single_fit_comparison.png")
single_fit_comparison_figure(
  fits, x_pred, truth, out,
  title = "Case 1: one dataset on the Hallstatt plateau, midpoint vs full distribution",
  obs = obs)
cat("wrote", out, "\n")
