# Shared data-generating process for the linear simulation.
# Sourced by 01_example.R (one dataset) and 02_recovery_study.R (many datasets),
# so the worked example and the study use the exact same process.

TMIN <- 100
TMAX <- 900

# Simulate one dataset from f(t) = intercept + slope * t with Gaussian noise.
# Each observation gets its own dating window; widths vary across observations
# (0.3x to 1.7x mean_width), so a dataset holds a realistic mix of precise and
# vague dates. Dates are whole years; the true date is uniform within each window.
simulate_linear <- function(N, intercept, slope, sigma, mean_width, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  w          <- runif(N, 0.3 * mean_width, 1.7 * mean_width)
  start_date <- round(runif(N, TMIN, TMAX - w))
  end_date   <- round(start_date + w)
  true_date  <- round(runif(N, start_date, end_date))
  data.frame(
    Start_date = start_date,
    End_date   = end_date,
    True_date  = true_date,
    Value      = round(intercept + slope * true_date + rnorm(N, 0, sigma), 1)
  )
}
