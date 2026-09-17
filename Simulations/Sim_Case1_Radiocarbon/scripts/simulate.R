#' CASE 1: radiocarbon dates.
#'
#' Each find has a true calendar date. The lab measures its radiocarbon age with
#' some error, and the archaeologist calibrates that age back into a probability
#' distribution over calendar years.
#'
#' Calendar years are BCE/CE (negative for BCE). They are converted to cal BP
#' only when rcarbon needs it: BP = 1950 - year.

#' True dates of deposition across the study window.
#'
#' growth_ratio is how many times denser finds are at the end of the window than
#' at the start. 1 is even. The draw inverts the cumulative distribution of an
#' exponential density.
rdeposition <- function(N, window, growth_ratio = 1) {
  if (growth_ratio == 1) return(runif(N, window[1], window[2]))
  k <- log(growth_ratio) / diff(window)
  window[1] + log(1 + runif(N) * (growth_ratio - 1)) / k
}

#' One simulated radiocarbon dataset.
#'
#' @param N                        number of finds
#' @param intercept, slope, sigma  the true trend: Value = intercept + slope * date + noise
#' @param window                   study window, c(start, end) in BCE/CE
#' @param lab_error                lab error of each measurement, in 14C years
#' @param growth_ratio             see rdeposition()
#' @param seed                     random seed
simulate_c14 <- function(N, intercept, slope, sigma, window, lab_error,
                         growth_ratio = 1, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)

  # 1. When each find was deposited
  true_date <- round(rdeposition(N, window, growth_ratio))

  # 2. Its radiocarbon age: the IntCal20 value at that date, plus noise. The
  # noise combines the lab error and the curve's own error, which is what
  # rcarbon assumes when it calibrates.
  curve <- rcarbon::uncalibrate(1950 - true_date, verbose = FALSE)
  cra   <- round(rnorm(N, curve$ccCRA, sqrt(lab_error^2 + curve$ccError^2)))

  data.frame(
    CRA       = cra,
    Error     = lab_error,
    True_date = true_date,
    Value     = round(intercept + slope * true_date + rnorm(N, 0, sigma), 1)
  )
}

#' How far the calendar grid must extend past the study window.
#'
#' The model is not told the study window, so every calibrated date has to fit
#' on the grid whole. This calibrates dates across the window and returns the
#' widest range that holds 99.99% of the probability.
c14_grid_pad <- function(window, lab_error) {
  years <- seq(window[1], window[2], by = 10)
  cra   <- round(rcarbon::uncalibrate(1950 - years, verbose = FALSE)$ccCRA)
  cal   <- calmatrix_to_calendar(
    rcarbon::calibrate(cra, errors = rep(lab_error, length(cra)),
                       calCurves = "intcal20", calMatrix = TRUE, verbose = FALSE))
  ranges <- apply(cal$prob, 2, function(p) {
    cum <- cumsum(p / sum(p))
    cal$years[which(cum >= 0.99995)[1]] - cal$years[which(cum >= 0.00005)[1]]
  })
  ceiling(max(ranges))
}
