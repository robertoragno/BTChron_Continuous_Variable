# Shared data-generating process for the GP simulation.
# Sourced by 01_example.R (one dataset) and 02_recovery_study.R (many datasets),
# so the worked example and the study use the same process.
#
# The truth is a smooth, gently undulating curve rather than a single bell. A bell
# is locally almost linear away from its peak, so the midpoint shortcut barely
# suffers; the interesting failure for a nonlinear trend is where the curve bends
# inside a dating window. The curve is a sum of two sine waves with periods longer
# than a typical dating window (350-550 years). This keeps two things in tension:
# the curve must bend within the windows for the midpoint shortcut to go wrong,
# but a wave shorter than the window is averaged away by it and cannot be recovered
# by any model (a window spanning a full oscillation carries no information about
# it). Periods a bit longer than the windows sit on the right side of that line.

TMIN <- 100
TMAX <- 900

# Number of sine components summed into the truth.
N_COMPONENTS <- 2

# Dirichlet broken stick (after EC diristick in beyond_aoristic): split
# [TMIN, TMAX] into K phases with Dirichlet(alpha_conc) lengths, then give each
# observation the [start, end] of the phase it falls in. Low alpha_conc makes
# coarser phases, high alpha_conc more even ones. This is the same partition used
# by the linear and changepoint simulations, so dating resolution is comparable
# across all three. Boundaries are whole years; H is the Shannon entropy of the
# phase weights.
partition_timeline <- function(N, K, alpha_conc, min_years = 1) {
  weights <- rgamma(K, shape = alpha_conc, rate = 1)
  weights <- weights / sum(weights)

  # reserve min_years per phase, then Dirichlet-share the rest of the span
  phase_spans <- min_years + weights * ((TMAX - TMIN) - K * min_years)
  phase_p     <- phase_spans / (TMAX - TMIN)

  boundaries        <- round(TMIN + c(0, cumsum(phase_spans)))
  boundaries[1]     <- TMIN
  boundaries[K + 1] <- TMAX

  true_date <- round(runif(N, TMIN, TMAX))
  phase     <- findInterval(true_date, boundaries, rightmost.closed = TRUE,
                            all.inside = TRUE)

  windows <- data.frame(
    Start_date = boundaries[phase],
    End_date   = boundaries[phase + 1],
    True_date  = true_date
  )
  attr(windows, "H")          <- -sum(phase_p * log(phase_p))
  attr(windows, "K")          <- K
  attr(windows, "alpha_conc") <- alpha_conc
  attr(windows, "weights")    <- phase_p
  windows
}

# Draw one random smooth curve: a baseline plus two sine waves. Periods sit in
# 350-550 years, longer than a typical dating window, so the undulation is
# recoverable while still bending inside the windows. Returns the parameters; the
# curve is evaluated by f_true.
draw_curve <- function() {
  list(
    baseline = runif(1, 2, 15),
    amp      = runif(N_COMPONENTS, 3, 8),
    period   = runif(N_COMPONENTS, 350, 550),
    phase    = runif(N_COMPONENTS, 0, 2 * pi)
  )
}

# The true curve value at times t (scalar or vector), for a curve from
# draw_curve. This is the hidden signal we simulate from and later score the fits
# against; the models never see it. They only get noisy values at dates that are
# themselves uncertain (known to fall in a window, not exactly when). The truth
# is built from sine waves on purpose, a different shape from the GP the model
# assumes, so recovering it is a fair test rather than a rigged one.
f_true <- function(t, curve) {
  components <- sapply(seq_along(curve$amp), function(k)
    curve$amp[k] * sin(2 * pi * (t - TMIN) / curve$period[k] + curve$phase[k]))
  if (is.null(dim(components))) components <- matrix(components, nrow = length(t))
  curve$baseline + rowSums(components)
}

# One dataset: dates placed by the broken stick above, values are the true curve
# at each hidden date plus Gaussian noise, rounded like a recorded measurement.
# H is carried along as an attribute.
simulate_gp <- function(N, sigma, K, alpha_conc, curve, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  windows <- partition_timeline(N, K, alpha_conc)
  out <- data.frame(
    Start_date = windows$Start_date,
    End_date   = windows$End_date,
    True_date  = windows$True_date,
    Value      = round(f_true(windows$True_date, curve) + rnorm(N, 0, sigma), 1)
  )
  attr(out, "H")          <- attr(windows, "H")
  attr(out, "K")          <- K
  attr(out, "alpha_conc") <- alpha_conc
  out
}

# Dense version of the true curve on a yearly grid, for plotting and for checking
# how much of the curve the fitted band covers.
truth_grid <- function(curve, by = 1) {
  years <- seq(TMIN, TMAX, by = by)
  data.frame(Year = years, True_value = f_true(years, curve))
}
