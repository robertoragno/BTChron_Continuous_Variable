# Shared data-generating process for the linear simulation.
# Sourced by 01_example.R (one dataset) and 02_recovery_study.R (many datasets),
# so the paper example and the study use the same process.

TMIN <- 100
TMAX <- 900

#' Dirichlet broken stick (after EC diristick in beyond_aoristic): split [TMIN, TMAX] 
#' into K phases with Dirichlet(alpha_conc) lengths, then give each observation 
#' the [start, end] of the phase it falls in. 
#' Low alpha_conc makes coarser phases, high alpha_conc even phases. 
#' min_years is a guard against zero-width phases 
#' IMPORTANT: IF WE MAKE MIN_YEARS SOMETHING AS 25, IT WILL FORCE THE H DISTRIBUTION TO AVOID LOW VALUES
#' WHICH IS A BIT UNREALISTIC. This at least for a dataset of 800 years.
#' Boundaries are whole years. Attaches the phase proportions and entropy H.
partition_timeline <- function(N, K, alpha_conc, min_years = 1) {
  weights <- rgamma(K, shape = alpha_conc, rate = 1)
  weights <- weights / sum(weights)

  # reserve min_years per phase, then Dirichlet-share the rest of the span
  phase_spans <- min_years + weights * ((TMAX - TMIN) - K * min_years)
  phase_p     <- phase_spans / (TMAX - TMIN)

  boundaries  <- round(TMIN + c(0, cumsum(phase_spans)))
  boundaries[1]       <- TMIN
  boundaries[K + 1]   <- TMAX

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

# One dataset: y = intercept + slope * true_date + Gaussian noise, with dates
# placed by the broken stick fn above. Includes H as an attribute
simulate_linear <- function(N, intercept, slope, sigma, K, alpha_conc,
                            seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  windows <- partition_timeline(N, K, alpha_conc)
  out <- data.frame(
    Start_date = windows$Start_date,
    End_date   = windows$End_date,
    True_date  = windows$True_date,
    Value      = round(intercept + slope * windows$True_date + rnorm(N, 0, sigma), 1)
  )
  attr(out, "H")          <- attr(windows, "H")
  attr(out, "K")          <- attr(windows, "K")
  attr(out, "alpha_conc") <- attr(windows, "alpha_conc")
  out
}
