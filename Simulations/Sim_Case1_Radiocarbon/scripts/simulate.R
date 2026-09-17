#' CASE 1: radiocarbon dates.
#'
#' Simulate radiocarbon dates with a lab error, then calibrate the age back into
#' a probability distribution over calendar years. 
#' Notes on conventions:
#' - Calendar years are BCE/CE (negative for BCE). They are converted to cal BP
#'   only when functions from rcarbon need it: BP = 1950 - year.
#' 

#' True dates of deposition across the study window.
#'
#' growth_ratio: how many times denser finds are at the end of the window compared
#' to the start. 1 is uniform. The draw inverts the cumulative distribution of an
#' exponential density.
rdeposition <- function(N, window, growth_ratio = 1) {
  if (growth_ratio == 1) return(runif(N, window[1], window[2]))
  k <- log(growth_ratio) / diff(window)
  window[1] + log(1 + runif(N) * (growth_ratio - 1)) / k
}

#' One simulated radiocarbon dataset.
#'
#' @param N                        number of samples
#' @param intercept, slope, sigma  simulated trend: Value = intercept + slope * date + noise
#' @param window                   study temporal window, c(start, end) in BCE/CE
#' @param lab_error                lab error of each measurement, in 14C years
#' @param growth_ratio             see rdeposition()
#' @param seed                     random seed
simulate_c14 <- function(N, intercept, slope, sigma, window, lab_error,
                         growth_ratio = 1, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)

  # Simulated deposition date
  true_date <- round(rdeposition(N, window, growth_ratio))

  # Uncalibrate to get the BP age (IntCal20 value at that date + noise)
  # The noise combines the lab error and the curve's own error
  curve <- rcarbon::uncalibrate(1950 - true_date, verbose = FALSE)
  cra   <- round(rnorm(N, curve$ccCRA, sqrt(lab_error^2 + curve$ccError^2)))
  # Memo: ccCRA so we get the uncalibrated age without the calibration curve error 

  data.frame(
    CRA       = cra,
    Error     = lab_error,
    True_date = true_date,
    Value     = round(intercept + slope * true_date + rnorm(N, 0, sigma), 1)
  )
}

#' Grid padding for Case 1
#'
#' The model is not told the study window, so each calibrated date must fit on
#' the grid whole. Calibrates one date every 10 yr across the window and returns
#' the width (yr) of the widest 99.99% range (this should remove the super tiny 
#' probabilities far from the peak - so the grid is smaller), the same share 
#' calibrated_rows() keeps of each date.
c14_grid_pad <- function(window, lab_error) {

  # Calibrate one date every 10 yr across the window
  years <- seq(window[1], window[2], by = 10)
  cra   <- round(rcarbon::uncalibrate(1950 - years, verbose = FALSE)$ccCRA)

  # calmatrix_to_calendar() checks the calMatrix is in ascending calendar order
  cal   <- calmatrix_to_calendar(
    rcarbon::calibrate(cra, errors = rep(lab_error, length(cra)),
                       calCurves = "intcal20", calMatrix = TRUE, verbose = FALSE))
  ranges <- apply(cal$prob, 2, function(p) {
    cum <- cumsum(p / sum(p))
    cal$years[which(cum >= 0.99995)[1]] - cal$years[which(cum >= 0.00005)[1]]
  })
  ceiling(max(ranges))
}
