# Case 1 single-fit figures: one dataset per window, fitted by the midpoint,
# calibrated-median and full-distribution models.
#   figures/single_fit_comparison.png        Hallstatt plateau
#   figures/single_fit_comparison_steep.png  steep section

library(here)
library(cmdstanr)
suppressMessages(library(rcarbon))

source(here("Simulations", "shared", "scripts", "weight_rows.R"))
source(here("Simulations", "shared", "scripts", "fit_models.R"))
source(here("Simulations", "shared", "scripts", "single_fit_figure.R"))
source(here("Simulations", "Sim_Case1_Radiocarbon", "scripts", "simulate.R"))

LAB_ERROR <- 25
INTERCEPT <- 8; SLOPE <- 0.02; SIGMA <- 2
N <- 50; STEP <- 5

# The two windows, as in 02_figures.R and the recovery study
windows      <- list(c(-800, -400), c(-1600, -1200))
window_names <- c("the Hallstatt plateau", "the steep section")
file_names   <- c("single_fit_comparison.png", "single_fit_comparison_steep.png")

for (k in 1:2) {
  WINDOW <- windows[[k]]

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
                             x_pred, models = c("midpoint", "median", "marginal"),
                             seed = 1)

  # What each model reads from the calibrated dates:
  # midpoint of the 95% range, calibrated median, and the 95% range itself
  midpoints <- data.frame(x = (smry$start + smry$end) / 2, y = sim$Value)
  medians   <- data.frame(x = smry$median, y = sim$Value)
  ranges    <- data.frame(xmin = smry$start, xmax = smry$end, y = sim$Value)

  out <- here("Simulations", "Sim_Case1_Radiocarbon", "figures", file_names[k])
  single_fit_comparison_figure(
    fits, x_pred, truth, out,
    title = paste("Case 1: one dataset on", window_names[k],
                  "- three models"),
    obs = midpoints, median_obs = medians, ranges = ranges,
    obs_label = "midpoint of 95% calibrated range",
    range_label = "95% calibrated range")
  cat("wrote", out, "\n")
}
