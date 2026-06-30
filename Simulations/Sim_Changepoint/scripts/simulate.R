# Shared data-generating process for the changepoint simulation.
# Sourced by 01_example.R (one dataset) and 02_recovery_study.R (many datasets),
# so the worked example and the study use the exact same process.

TMIN <- 100
TMAX <- 900

# Dirichlet broken stick (Crema's diristick): split [TMIN, TMAX] into K phases
# with Dirichlet(alpha_conc) lengths, then give each observation the [start, end]
# of the phase it falls in. Low alpha_conc makes lumpy phases, high alpha_conc
# even ones. Boundaries are whole years. Attaches the phase proportions and entropy H.
partition_timeline <- function(
  N, # How many observations to generate?
  K, # How many phases to split the timeline into?
  alpha_conc, # Dirichlet concentration parameter for the phase lengths
  min_years = 1 # Minimum number of years per phase (default to 1 as a guard against zero-width phases)
  ) {
  weights <- rgamma(K, shape = alpha_conc, rate = 1)
  weights <- weights / sum(weights) # Normalize to sum to 1

  # Reserve min_years per phase, then Dirichlet-share the rest of the span
  # min_years should be small enough that K * min_years < TMAX - TMIN
  phase_spans <- min_years + weights * ((TMAX - TMIN) - K * min_years)
  phase_p     <- phase_spans / (TMAX - TMIN)

  boundaries  <- round(TMIN + c(0, cumsum(phase_spans))) # cumsum so we keep the position of the boundaries
  # Just to make sure that the boundaries are exactly TMIN and TMAX, in case of rounding issues
  boundaries[1]       <- TMIN
  boundaries[K + 1]   <- TMAX

  true_date <- round(runif(N, TMIN, TMAX))
  phase     <- findInterval(true_date, boundaries, rightmost.closed = TRUE,
                            all.inside = TRUE) # Get the phase index for each true_date
                            # Basically asking how many boundaries are to the left of each true_date

  windows <- data.frame(
    Start_date = boundaries[phase],
    End_date   = boundaries[phase + 1],
    True_date  = true_date
  )
  # Add attributes (to avoid meaningless columns in the data frame)
  attr(windows, "H")          <- -sum(phase_p * log(phase_p))
  attr(windows, "K")          <- K
  attr(windows, "alpha_conc") <- alpha_conc
  attr(windows, "weights")    <- phase_p
  windows
}

# One dataset from a two-slope trend that bends at cp (continuous at the join):
#   f(t) = baseline + slope_1 * (t - TMIN)                        for t <= cp
#   f(t) = baseline + slope_1 * (cp - TMIN) + slope_2 * (t - cp)  for t >  cp
# Dates placed by the broken stick above; carries H as an attribute.
simulate_changepoint <- function(N, baseline, slope_1, slope_2, cp, sigma,
                                 K, alpha_conc, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  windows   <- partition_timeline(N, K, alpha_conc)
  true_date <- windows$True_date
  mu <- ifelse(true_date <= cp,
               baseline + slope_1 * (true_date - TMIN),
               baseline + slope_1 * (cp - TMIN) + slope_2 * (true_date - cp))
  out <- data.frame(
    Start_date = windows$Start_date,
    End_date   = windows$End_date,
    True_date  = true_date,
    Value      = round(mu + rnorm(N, 0, sigma), 1)
  )
  attr(out, "H")          <- attr(windows, "H")
  attr(out, "K")          <- attr(windows, "K")
  attr(out, "alpha_conc") <- attr(windows, "alpha_conc")
  out
}
