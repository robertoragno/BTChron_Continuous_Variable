# Shared data-generating process for the changepoint simulation.
# Sourced by 01_example.R (one dataset) and 02_recovery_study.R (many datasets),
# so the worked example and the study use the exact same process.

TMIN <- 100
TMAX <- 900

# Simulate one dataset from a two-slope broken-stick model, continuous at cp:
#   f(t) = baseline + slope_1 * (t - TMIN)                              for t <= cp
#   f(t) = baseline + slope_1 * (cp - TMIN) + slope_2 * (t - cp)       for t >  cp
#
# Each observation gets its own dating window; widths vary across observations
# (0.3x to 1.7x mean_width), so a dataset holds a realistic mix of precise and
# vague dates. Dates are whole years; the true date is uniform within each window.
simulate_changepoint <- function(N, baseline, slope_1, slope_2, cp, sigma,
                                 mean_width, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  w          <- runif(N, 0.3 * mean_width, 1.7 * mean_width)
  start_date <- round(runif(N, TMIN, TMAX - w))
  end_date   <- round(start_date + w)
  true_date  <- round(runif(N, start_date, end_date))
  mu <- ifelse(true_date <= cp,
               baseline + slope_1 * (true_date - TMIN),
               baseline + slope_1 * (cp - TMIN) + slope_2 * (true_date - cp))
  data.frame(
    Start_date = start_date,
    End_date   = end_date,
    True_date  = true_date,
    Value      = round(mu + rnorm(N, 0, sigma), 1)
  )
}
