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
  # midpoint of the 95% range, calibrated median, and the calibrated shape
  # itself (the same weight rows marginal_date.stan is fitted on)
  midpoints <- data.frame(x = (smry$start + smry$end) / 2, y = sim$Value)
  medians   <- data.frame(x = smry$median, y = sim$Value)

  # Drawing every date's shape floods the plateau panel into a solid wash, so
  # only a subsample is drawn, spread evenly across the dataset's calibrated
  # medians so the whole window is represented. The fit still uses all N.
  BLOB_N   <- 20
  blob_idx <- order(smry$median)[unique(round(seq(1, N, length.out = min(N, BLOB_N))))]
  blobs    <- calibrated_blob_data(dates$weights, grid, sim$Value,
                                   scale = 0.06 * diff(range(sim$Value)),
                                   rows = blob_idx, one_sided = TRUE)

  out <- here("Simulations", "Sim_Case1_Radiocarbon", "figures", file_names[k])
  single_fit_comparison_figure(
    fits, x_pred, truth, out,
    title = paste("Case 1: one dataset on", window_names[k],
                  "- three models"),
    subtitle = if (length(blob_idx) < N)
      sprintf("Panel C shows %d of %d calibrated dates for legibility; every model is fitted on all %d",
              length(blob_idx), N, N),
    obs = midpoints, median_obs = medians, blobs = blobs,
    obs_label = "midpoint of 95% calibrated range")
  cat("wrote", out, "\n")
}
