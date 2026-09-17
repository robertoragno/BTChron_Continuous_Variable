#' CASE 4: Archaeological finds dated to merged phases
#'
#' The study period is split into K fine phases of random length. Some finds
#' are dated to a single fine phase, others only to a broader period made of
#' 2 to merge_max adjacent fine phases (e.g. "Late Roman" instead of "Late
#' Roman 2").
#' Memo for the readme/paper:
#' - Each find is assumed to be independent.
#' - The broad period is a run of phases drawn per find, not a fixed hierarchy.
#' Whether a fixed (nested) hierarchy is the better design is still open.

#' Function simulate_merged to simulate a single (Case 4) dataset.
#'
#' @param N                    number of samples within a dataset
#' @param intercept, slope, sigma; the simulated known trend:
#'        Value = intercept + slope * date + noise
#' @param K                    number of fine phases
#' @param alpha_conc           how even the phase lengths are: low is uneven,
#'                             high is even
#' @param merge_max            largest number of fine phases a broad period can
#'                             span. Each find draws its own span from 1 to
#'                             merge_max; 1 keeps the find in its fine phase.
#' @param period_start, period_end  the study period
#' @param seed                 random seed
simulate_merged <- function(N, intercept, slope, sigma, K, alpha_conc,
                            merge_max = 3, period_start = 100,
                            period_end = 900, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  if (merge_max > K) stop("merge_max cannot exceed the number of fine phases")
  span <- period_end - period_start

  # Phase lengths from a Dirichlet "broken stick", at least 1 yr each
  weights <- rgamma(K, shape = alpha_conc, rate = 1)
  weights <- weights / sum(weights)
  lengths <- 1 + weights * (span - K)
  bounds  <- round(period_start + c(0, cumsum(lengths)))
  bounds[c(1, K + 1)] <- c(period_start, period_end)

  true_date <- round(runif(N, period_start, period_end))
  fine      <- findInterval(true_date, bounds, rightmost.closed = TRUE,
                            all.inside = TRUE)

  # A find's broad period is a run of m phases that includes its own. Where
  # the run starts is random, so wide windows are not more common early or late.
  m    <- sample.int(merge_max, N, replace = TRUE)
  back <- floor(runif(N) * m)
  lo   <- pmax(fine - back, 1)
  hi   <- pmin(lo + m, K + 1)
  lo   <- pmax(hi - m, 1)

  out <- data.frame(Start_date = bounds[lo],
                    End_date   = bounds[hi],
                    True_date  = true_date,
                    Phases     = hi - lo,   # fine phases in the window
                    Value      = round(intercept + slope * true_date +
                                       rnorm(N, 0, sigma), 1))
  attr(out, "bounds") <- bounds   # phase boundaries, for the figures
  out
}
