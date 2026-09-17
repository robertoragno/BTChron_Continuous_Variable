#' CASE 3: Archaeological finds dated to overlapping phases (example of
#' ceramic phases)
#'
#' The study period is split into K phases of random length. Neighbouring
#' phases overlap: each phase's recorded window is stretched past its
#' boundaries, so a find dated to the shared strip could belong to either
#' phase and is assigned to one of them.
#' Memo for the readme/paper:
#' - Each find is assumed to be independent.
#' - The phase windows are the same for every find in a dataset (unlike Case 2,
#' where each find has its own window).
#' - Inside a stretched window the true dates are not flat: the shared strips
#' hold half as many dates as the phase's own span (a trapezoid, not a rectangle).

#' Function simulate_overlap to simulate a single (Case 3) dataset.
#'
#' @param N                    number of samples within a dataset
#' @param intercept, slope, sigma; the simulated known trend:
#'        Value = intercept + slope * date + noise
#' @param K                    number of phases
#' @param alpha_conc           how even the phase lengths are: low is uneven,
#'                             high is even
#' @param overlap              how far each phase's window is stretched past its
#'                             boundaries, as a share of the phase's own length
#'                             (half on each side). 0 means phases just touch.
#' @param assign_p             for a find in a shared strip, the chance it is
#'                             assigned to the earlier phase. 0.5 is a coin flip.
#' @param period_start, period_end  the study period
#' @param seed                 random seed
simulate_overlap <- function(N, intercept, slope, sigma, K, alpha_conc,
                             overlap = 0.25, assign_p = 0.5,
                             period_start = 100, period_end = 900,
                             seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  span <- period_end - period_start

  # Phase lengths from a Dirichlet "broken stick", at least 1 yr each
  weights <- rgamma(K, shape = alpha_conc, rate = 1)
  weights <- weights / sum(weights)
  lengths <- 1 + weights * (span - K)
  bounds  <- round(period_start + c(0, cumsum(lengths)))
  bounds[c(1, K + 1)] <- c(period_start, period_end)

  true_date <- round(runif(N, period_start, period_end))
  phase     <- findInterval(true_date, bounds, rightmost.closed = TRUE,
                            all.inside = TRUE)

  # Recorded window of each phase, stretched on both sides. The first and last
  # phases run past the study period too.
  reach  <- overlap * diff(bounds) / 2
  win_lo <- round(bounds[-(K + 1)] - reach)
  win_hi <- round(bounds[-1] + reach)

  if (overlap > 0) {
    # A find in the strip shared with the previous or the next phase can be
    # assigned to that phase instead. A find in both strips keeps its own phase.
    in_prev <- phase > 1 & true_date <= win_hi[pmax(phase - 1, 1)]
    in_next <- phase < K & true_date >= win_lo[pmin(phase + 1, K)]
    only_prev <- in_prev & !in_next
    only_next <- in_next & !in_prev

    flip <- runif(N)
    phase[only_prev] <- phase[only_prev] - (flip[only_prev] < assign_p)
    phase[only_next] <- phase[only_next] + (flip[only_next] >= assign_p)
  }

  out <- data.frame(Start_date = win_lo[phase],
                    End_date   = win_hi[phase],
                    True_date  = true_date,
                    Value      = round(intercept + slope * true_date +
                                       rnorm(N, 0, sigma), 1))
  attr(out, "bounds") <- bounds   # phase boundaries, for the figures
  out
}
