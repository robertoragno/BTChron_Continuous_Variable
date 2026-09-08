# Dirichlet broken-stick periodisation, shared by the phase-based cases.
#
# Copied from archive/Sim_Linear/scripts/simulate.R and extended with an
# `overlap` argument (Case 3) and an `assign_p` argument for how ambiguous dates
# in an overlap are labelled. The abutting boundaries ride along as an attribute
# so Case 4's merges do not have to rebuild the broken stick.
#
# The archived copy is left untouched, so its results stay reproducible. At
# overlap = 0 this function draws the same rgamma and runif in the same order as
# the archived one and returns the same windows year for year on a shared seed;
# Sim_Case3_OverlappingPhases/scripts/00_check.R asserts it.
#
# The calendar span is an argument here rather than the globals TMIN / TMAX the
# archived version reads off its own file, so nothing depends on load order.

#' Split [t_min, t_max] into K phases of Dirichlet-random width, then date N
#' observations to the phase each falls in.
#'
#' @param N,K         observations, phases
#' @param alpha_conc  Dirichlet concentration; low is lumpy, high is even
#' @param overlap    fraction of a phase's own width by which its recorded
#'                    window is stretched past each abutting boundary. 0 is the
#'                    archived abutting partition. The stretch is skipped at the
#'                    two ends of [t_min, t_max], so coverage stays gap-free
#'                    and no window runs past the span.
#' @param assign_p   for a date in the overlap of an earlier and a later phase,
#'                    the probability it is labelled the EARLIER one. 0.5 is a
#'                    coin flip and keeps the retained dates symmetric about the
#'                    window midpoint (Berkson holds); below 0.5 leans
#'                    assignments toward the later phase, the usual "ambiguous
#'                    sherd goes to the better-attested phase". Pulls the
#'                    midpoint off the retained mean per phase. Ignored when
#'                    overlap = 0.
#' @param min_years  floor on phase width, a guard against zero-width phases
#' @return data.frame(Start_date, End_date, True_date), with attributes H
#'         (entropy of the phase weights), K, alpha_conc, overlap, assign_p,
#'         and bounds (the K + 1 abutting boundaries, for Case 4)
partition_timeline <- function(N, K, alpha_conc, overlap = 0, assign_p = 0.5,
                               t_min = 100, t_max = 900, min_years = 1) {
  weights <- rgamma(K, shape = alpha_conc, rate = 1)
  weights <- weights / sum(weights)

  # reserve min_years per phase, then Dirichlet-share the rest of the span
  phase_spans <- min_years + weights * ((t_max - t_min) - K * min_years)
  phase_p     <- phase_spans / (t_max - t_min)

  bounds        <- round(t_min + c(0, cumsum(phase_spans)))
  bounds[1]     <- t_min
  bounds[K + 1] <- t_max

  true_date <- round(runif(N, t_min, t_max))

  # Base assignment is the abutting partition: deterministic, and the only place
  # this function would draw an extra random number is the overlap branch below,
  # which is skipped entirely at overlap = 0. So the archived RNG stream is
  # reproduced exactly there.
  phase <- findInterval(true_date, bounds, rightmost.closed = TRUE,
                        all.inside = TRUE)

  # Stretched windows, kept on whole years. The stretch is outward on both
  # sides, including at the two ends of [t_min, t_max]: an analyst's earliest
  # ceramic phase does not stop at the year the study window opens, and clamping
  # the end phases there compresses the midpoint range against the true-date
  # range and inflates the point-date slope, the same artifact Case 2 found and
  # rejected when it chose overhang over clamped windows. The phase-based case
  # that uses this pads its study grid to cover the overhang instead. At
  # overlap = 0 round() is a no-op on the already integer boundaries, which keeps
  # the archived output reproducible.
  reach  <- overlap * diff(bounds) / 2
  win_lo <- round(bounds[-(K + 1)] - reach)
  win_hi <- round(bounds[-1]       + reach)

  if (overlap > 0) {
    # With the windows stretched, a date can now fall inside its own phase and
    # the stretched end of one neighbour. `in_prev` is the strip it shares with
    # the earlier phase, `in_next` the strip it shares with the later one. A date
    # claimed by both (only possible for a phase much narrower than both
    # neighbours) keeps its base phase rather than being pulled twice.
    in_prev <- phase > 1 & true_date <= win_hi[pmax(phase - 1L, 1L)]
    in_next <- phase < K & true_date >= win_lo[pmin(phase + 1L, K)]
    solo_prev <- in_prev & !in_next
    solo_next <- in_next & !in_prev

    flip <- runif(N)
    # earlier phase of the [prev, this] overlap is `phase - 1`
    phase[solo_prev] <- phase[solo_prev] - (flip[solo_prev] < assign_p)
    # earlier phase of the [this, next] overlap is `phase` itself
    phase[solo_next] <- phase[solo_next] + (flip[solo_next] >= assign_p)
  }

  out <- data.frame(Start_date = win_lo[phase],
                    End_date   = win_hi[phase],
                    True_date  = true_date)
  attr(out, "H")          <- -sum(phase_p * log(phase_p))
  attr(out, "K")          <- K
  attr(out, "alpha_conc") <- alpha_conc
  attr(out, "overlap")    <- overlap
  attr(out, "assign_p")   <- assign_p
  attr(out, "bounds")     <- bounds
  out
}
