# Case 2: typochronology-derived dates. Every sample carries its own independent
# window, with no shared phase structure (that is Case 3).
#
# Widths are set as a share of the studied period rather than in years, so the
# same simulation covers a Roman sequence and a prehistoric one. 5-25% of the
# period is 40-200 yr on an 800-yr Roman study and 175-875 yr on a 3500-yr
# Neolithic one.
#
# The true date is drawn first and the window placed around it. Drawing the
# window first also works, but only if the window is not forced to sit inside
# the period: bounding it by (period_end - width) makes wide windows cluster in
# the middle of the timeline and leaves the true dates triangular instead of
# uniform.

#' Width set for a period, from a fraction range. Snapped to the grid, so the
#' boundaries look like real datings rather than arbitrary numbers.
typo_widths <- function(period_start, period_end,
                        min_frac = 0.05, max_frac = 0.25, grid = 25) {
  span <- period_end - period_start
  lo   <- max(grid, round(min_frac * span / grid) * grid)
  hi   <- round(max_frac * span / grid) * grid
  seq(lo, hi, by = grid)
}

#' Relative chance of each width, fine to coarse. This is Case 2's dating
#' resolution factor, the counterpart of H in the phase-based cases: "fine"
#' leans to tight windows, "coarse" to broad ones, "even" is flat.
#' strength is how lopsided the tilt is - 1 is flat whatever the mix, higher
#' values make fine/coarse more extreme.
typo_width_prob <- function(widths, mix = "even", strength = 4) {
  f <- seq(0, 1, length.out = length(widths))
  switch(mix,
         even   = rep(1, length(widths)),
         fine   = strength^(-f),
         coarse = strength^(f),
         stop("mix must be one of: even, fine, coarse"))
}

#' One typochronology-style dataset.
#'
#' @param N            sample size
#' @param intercept, slope, sigma  linear trend, as in simulate_linear()
#' @param period_start, period_end the studied period
#' @param widths       window widths to draw from, see typo_widths()
#' @param width_prob   relative chance of each width, see typo_width_prob()
#' @param skew_shape   Beta(skew_shape, 1) for where in its window the true date
#'                      sits. 1 = uniform (2a), > 1 leans to the late edge (2b).
#' @param grid         window boundaries are multiples of this
#' @param seed         RNG seed
#'
#' Windows may run past period_start / period_end. A typological attribution does
#' not know where the study window stops, and clipping them creates degenerate
#' one-year windows at the edges. At a 25% maximum width this costs about 2% in
#' the midpoint slope, reported as a limitation.
simulate_typo <- function(N, intercept, slope, sigma,
                          period_start, period_end,
                          widths = typo_widths(period_start, period_end),
                          width_prob = NULL, skew_shape = 1, grid = 25,
                          seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  if (period_end <= period_start) stop("period_end must be after period_start")
  if (min(widths) < grid) stop("widths must be at least one grid step")

  true_date <- round(runif(N, period_start, period_end))
  width     <- sample(widths, N, replace = TRUE, prob = width_prob)

  # Pick the window start straight from the grid multiples that already contain
  # the date, so nothing needs repairing afterwards. k indexes those multiples:
  # k_hi puts the date at the start of its window, k_lo at the end, and the Beta
  # draw chooses where in between.
  #
  # floor() over (n_steps + 1) buckets rather than round() over n_steps: a
  # 50-year window on a 25-year grid has only three legal starts, and rounding
  # hands the middle one half the mass instead of a third, which shows up as a
  # hump in the middle of the within-window position and quietly breaks 2a.
  k_lo    <- ceiling((true_date - width) / grid)
  k_hi    <- floor(true_date / grid)
  n_steps <- k_hi - k_lo
  offset  <- pmin(floor(rbeta(N, skew_shape, 1) * (n_steps + 1)), n_steps)
  start_date <- grid * (k_hi - offset)

  data.frame(
    Start_date = start_date,
    End_date   = start_date + width,
    True_date  = true_date,
    Value      = round(intercept + slope * true_date + rnorm(N, 0, sigma), 1)
  )
}
