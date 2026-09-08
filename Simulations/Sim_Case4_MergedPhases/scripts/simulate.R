# Case 4: nested phases. A find is dated to a broad period that covers several
# adjacent fine phases, so its window is the union of those phases. The union is
# still one interval, just wider than a single phase, and the true date stays
# uniform inside it, so the midpoint is unbiased on the slope here as in Case 2.
# The case is about the mix of resolutions within one dataset: some finds carry
# a fine-phase label, others only sit in the broad period. The usual response is
# to analyse the two groups separately or drop the coarse one. One weight row per
# find keeps them in the same model, tight for the fine-phase finds and wide and
# flat for the broad-period ones.
#
# The fine phases come from the same broken stick as Case 3. partition_timeline()
# is in shared/scripts/partition.R; every script that uses this file sources
# that one first.

#' One nested-phase dataset, in the Start_date / End_date / True_date / Value
#' shape the other cases return.
#'
#' @param merge_max largest number of fine phases a broad period can span. Each
#'   find draws its own span from 1..merge_max, so merge_max sets how coarse the
#'   coarse end of the dataset gets. A span of 1 leaves the find in its own fine
#'   phase.
#' @inheritParams partition_timeline
simulate_merged <- function(N, intercept, slope, sigma, K, alpha_conc,
                            merge_max = 3, t_min = 100, t_max = 900,
                            seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  if (merge_max > K) stop("merge_max cannot exceed the number of fine phases")

  base   <- partition_timeline(N, K, alpha_conc, overlap = 0,
                               t_min = t_min, t_max = t_max)
  bounds <- attr(base, "bounds")
  fine   <- findInterval(base$True_date, bounds, rightmost.closed = TRUE,
                         all.inside = TRUE)

  # A find's broad period is a run of m consecutive fine phases that includes
  # its own. Where the run sits relative to that phase is drawn per find, so
  # window width does not track calendar position. Merging only forwards would
  # give early finds wide windows and clamp late ones to a single phase, and a
  # width that trends along the timeline would bias the slope.
  m    <- sample.int(merge_max, N, replace = TRUE)
  back <- floor(runif(N) * m)
  lo   <- pmax(fine - back, 1L)
  hi   <- pmin(lo + m, K + 1L)
  lo   <- pmax(hi - m, 1L)

  out <- data.frame(
    Start_date = bounds[lo],
    End_date   = bounds[hi],
    True_date  = base$True_date,
    Value      = round(intercept + slope * base$True_date + rnorm(N, 0, sigma), 1)
  )
  attr(out, "H")         <- attr(base, "H")
  attr(out, "K")         <- K
  attr(out, "merge_max") <- merge_max
  attr(out, "bounds")    <- bounds
  attr(out, "span")      <- hi - lo   # fine phases per find, for the checks
  out
}
