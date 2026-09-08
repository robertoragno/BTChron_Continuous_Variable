# Case 1: radiocarbon dates. Each sample's date information is a calibrated
# posterior over calendar years, not a window, so the error is classical rather
# than Berkson: the assigned date comes from a noisy measurement of the true
# date, and attenuation is expected in the point-date models.
#
# Calendar years are BCE/CE throughout (negative for BCE), matching the other
# cases. Conversion to cal BP happens only where rcarbon needs it.

#' When finds were deposited, across the study window.
#'
#' growth_ratio is how many times denser deposition is at the end of the window
#' than at the start: 1 is uniform, 4 is fourfold denser at the late edge. A
#' ratio rather than a rate in years, so the same number means the same shape
#' whatever the window length.
#'
#' Drawn by inverting the CDF, which is three lines and exact. The density is
#' f(t) proportional to exp(k * (t - start)) with k = log(ratio) / span, so
#' F(t) = (exp(k * (t - start)) - 1) / (ratio - 1), and solving F(t) = u for t
#' gives the line below. At ratio 1 that k is zero and the formula divides by
#' zero, so uniform is its own branch - and it draws the identical runif() the
#' study used before, which keeps every uniform cell reproducible.
rdeposition <- function(N, window, growth_ratio = 1) {
  if (growth_ratio == 1) return(runif(N, window[1], window[2]))
  k <- log(growth_ratio) / diff(window)
  window[1] + log(1 + runif(N) * (growth_ratio - 1)) / k
}

#' One radiocarbon dataset.
#'
#' The determination is drawn from the IntCal20 mean at the true date, with the
#' lab error and the curve's own error added in quadrature. That second term is
#' what rcarbon puts in the denominator when it calibrates: its kernel is
#' dnorm(age, mu, sqrt(lab_error^2 + curve_error^2)). Drawing with the lab error
#' alone would mean generating the data from one model and inverting a wider
#' one, so the posteriors would be broader than the noise warrants and the
#' mismatch would be largest where the lab error is smallest - precisely along
#' the lab_error sweep. uncalibrate() returns the curve error at each year as
#' ccError, so nothing has to be looked up by hand.
#'
#' Still idealised in one way, and stated as a limitation rather than fixed: the
#' curve error is drawn independently per date, whereas there is only one
#' IntCal20 residual and it displaces nearby dates together.
simulate_c14 <- function(N, intercept, slope, sigma, window, lab_error,
                         growth_ratio = 1, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)

  true_date <- round(rdeposition(N, window, growth_ratio))
  curve     <- rcarbon::uncalibrate(1950 - true_date, verbose = FALSE)
  cra       <- round(rnorm(N, curve$ccCRA,
                           sqrt(lab_error^2 + curve$ccError^2)))

  data.frame(
    CRA       = cra,
    Error     = lab_error,
    True_date = true_date,
    Value     = round(intercept + slope * true_date + rnorm(N, 0, sigma), 1)
  )
}

#' How far the grid must extend beyond the study window.
#'
#' A date near the edge loses part of its calibrated posterior off the end of
#' the grid, and renormalising then asserts it is certainly inside the window,
#' dragging it inward. The padding is measured rather than guessed: calibrate a
#' sweep of determinations across the window and take the widest support seen.
#'
#' Called once per scenario and the result frozen into the design. Recomputing
#' it per dataset would make the same prior mean different things in different
#' datasets, which is the mistake already fixed in the models.
c14_grid_pad <- function(window, lab_error, mass_keep = 0.9999, by = 10) {
  sweep <- seq(window[1], window[2], by = by)
  cra   <- round(rcarbon::uncalibrate(1950 - sweep, verbose = FALSE)$ccCRA)
  x     <- rcarbon::calibrate(cra, errors = rep(lab_error, length(cra)),
                              calCurves = "intcal20", calMatrix = TRUE,
                              verbose = FALSE)
  cal   <- calmatrix_to_calendar(x)

  drop <- (1 - mass_keep) / 2
  span <- apply(cal$prob, 2, function(p) {
    cum <- cumsum(p / sum(p))
    cal$years[which(cum >= 1 - drop)[1]] - cal$years[which(cum >= drop)[1]]
  })
  ceiling(max(span))
}
