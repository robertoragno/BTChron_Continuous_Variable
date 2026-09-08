# Case 3: phase-based dates with overlapping phases, the ceramic-phase case.
#
# Phases come from the same Dirichlet broken stick as the archived Sim_Linear,
# but adjacent phases now share a strip instead of meeting at a point. A date in
# that strip is labelled one of the two phases, and once the labelling is not a
# straight coin flip the retained dates inside a recorded window stop being
# uniform: the window midpoint is no longer their mean, so the point-date slope
# picks up real bias, not just mis-calibration. That is the mechanism the case
# exists to show - see Simulations/shared/README.md and section 5 of the plan.
#
# The window is flat, so the median IS the midpoint here, exactly as in Case 2.
# Only midpoint and marginal are fitted; Case 2's 00_check_median_identity.R
# proves the identity once for every flat-window case.
#
# partition_timeline() lives in shared/scripts/partition.R; every script that
# uses this file sources that one first.

#' One overlapping-phase dataset, in the Start_date / End_date / True_date /
#' Value shape every case returns.
#'
#' @param overlap,assign_p passed straight through to partition_timeline()
#' @inheritParams partition_timeline
simulate_overlap <- function(N, intercept, slope, sigma, K, alpha_conc,
                             overlap = 0.25, assign_p = 0.5,
                             t_min = 100, t_max = 900, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  w <- partition_timeline(N, K, alpha_conc, overlap = overlap,
                          assign_p = assign_p, t_min = t_min, t_max = t_max)
  out <- data.frame(
    Start_date = w$Start_date,
    End_date   = w$End_date,
    True_date  = w$True_date,
    Value      = round(intercept + slope * w$True_date + rnorm(N, 0, sigma), 1)
  )
  attr(out, "H")        <- attr(w, "H")
  attr(out, "K")        <- attr(w, "K")
  attr(out, "overlap")  <- attr(w, "overlap")
  attr(out, "assign_p") <- attr(w, "assign_p")
  attr(out, "bounds")   <- attr(w, "bounds")
  out
}
