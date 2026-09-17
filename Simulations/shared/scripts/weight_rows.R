# Weight rows for the full-distribution model (marginal_date.stan).
#
# Every case gives the model the same thing: a shared grid of candidate calendar
# years, and one row of probabilities per find over that grid. A row is flat
# over a window (Cases 2-4) or a calibrated distribution (Case 1).
#
# Rows are stored trimmed to the years they cover, glued end to end into one
# vector, with the first grid position and length of each row. The grid years
# themselves are passed to Stan too, so a row can never be read in the wrong
# calendar order.

#' Pack rows into the form Stan reads.
#'
#' @param rows    list of list(index = first grid position, prob = probabilities)
#' @param T_grid  number of years on the grid
pack_rows <- function(rows, T_grid) {
  row_first_year <- as.integer(vapply(rows, function(r) r$index, numeric(1)))
  probs          <- lapply(rows, function(r) r$prob / sum(r$prob))   # each row sums to 1
  row_n_years    <- vapply(probs, length, integer(1))

  if (any(row_first_year < 1L) || any(row_first_year + row_n_years - 1L > T_grid))
    stop("a weight row falls outside the grid")
  if (any(!vapply(probs, function(p) all(is.finite(p)) && all(p >= 0), logical(1))))
    stop("weight rows must be finite and non-negative")

  list(n_weights            = sum(row_n_years),
       log_year_prob_packed = log(unlist(probs, use.names = FALSE)),
       row_first_year       = row_first_year,
       row_n_years          = row_n_years)
}

#' Cases 2-4: flat over each find's window.
#'
#' A window narrower than one grid step still gets one grid cell.
uniform_rows <- function(start_date, end_date, grid) {
  step <- if (length(grid) > 1) grid[2] - grid[1] else 1
  i_lo <- pmax(1L, ceiling((start_date - grid[1]) / step) + 1L)
  i_hi <- pmin(length(grid), floor((end_date - grid[1]) / step) + 1L)
  i_hi <- pmax(i_hi, i_lo)

  rows <- Map(function(lo, hi) list(index = lo, prob = rep(1, hi - lo + 1L)),
              i_lo, i_hi)
  pack_rows(rows, length(grid))
}

#' rcarbon's calMatrix as probabilities over ascending BCE/CE years.
#'
#' calMatrix rows are labelled in cal BP, oldest first, which is already
#' ascending calendar order: the labels need converting (year = 1950 - BP), the
#' rows must not be reversed. Reversing them would flip the sign of the slope
#' without any error, hence the check.
calmatrix_to_calendar <- function(cal_dates) {
  m <- cal_dates$calmatrix
  if (is.null(m)) stop("calibrate() must be called with calMatrix = TRUE")

  bp <- as.numeric(rownames(m))
  if (anyNA(bp)) stop("calMatrix rownames are not cal BP years")
  years <- 1950 - bp
  if (is.unsorted(years))
    stop("calMatrix rows are not in ascending calendar order after relabelling")

  list(prob = m, years = years)
}

#' Calibrated probabilities summed onto the grid: one column per date, one row
#' per grid year.
#'
#' On a grid coarser than annual, each grid cell collects the yearly
#' probabilities within half a step of it. Taking every fifth year instead would
#' throw away four fifths of each distribution.
#'
#' Stops if a date loses more than max_loss of its probability off the grid.
#' The model is not told the study window (a real analyst does not know it), so
#' calibrated dates are used whole and the grid is padded to hold them; cutting
#' them at the grid edge would pull edge dates inward.
calibrated_on_grid <- function(cal, grid, max_loss = 0.01) {
  step <- if (length(grid) > 1) grid[2] - grid[1] else 1

  # Grid cell each calendar year falls into (1 = first grid year). Several years
  # share a cell when the grid step is more than one year.
  cell   <- round((cal$years - grid[1]) / step) + 1L
  inside <- cell >= 1L & cell <= length(grid)
  summed <- rowsum(cal$prob[inside, , drop = FALSE], cell[inside], reorder = TRUE)

  # rowsum() drops empty cells; put them back so row i is grid year i
  m <- matrix(0, nrow = length(grid), ncol = ncol(cal$prob))
  m[as.integer(rownames(summed)), ] <- summed

  lost <- 1 - colSums(m)
  if (any(lost > max_loss))
    stop(sprintf("date(s) lose up to %.1f%% of their calibrated mass off the grid; widen the padding",
                 100 * max(lost)))
  m
}

#' Case 1: each date's calibrated distribution as a weight row.
#'
#' Each row is trimmed to the shortest continuous span holding mass_keep of its
#' probability. One continuous span, not only the non-zero years, because a
#' plateau date has near-empty gaps between its peaks.
calibrated_rows <- function(cal, grid, mass_keep = 0.9999, max_loss = 0.01) {
  if (nrow(cal$prob) != length(cal$years))
    stop("cal_matrix rows and cal_years must correspond")
  if (is.unsorted(cal$years))
    stop("cal_years must be ascending, matching the grid")
  step <- if (length(grid) > 1) grid[2] - grid[1] else 1
  if (min(cal$years) > min(grid) - step / 2 || max(cal$years) < max(grid) + step / 2)
    stop("the calibration range does not cover the study grid")

  m <- calibrated_on_grid(cal, grid, max_loss)

  rows <- lapply(seq_len(ncol(m)), function(j) {
    p <- m[, j]
    if (sum(p) <= 0) stop("date ", j, " has no probability mass on the grid")
    p    <- p / sum(p)
    drop <- (1 - mass_keep) / 2
    cum  <- cumsum(p)
    lo   <- max(1L, which(cum >= drop)[1])
    hi   <- min(length(p), which(cum >= 1 - drop)[1])
    list(index = lo, prob = p[lo:hi])
  })
  pack_rows(rows, length(grid))
}

#' Shortest single range holding `mass` of a date's probability (95% HPD
#' envelope). A multi-peaked date gets one range spanning all its peaks.
hpd_envelope <- function(p, years, mass = 0.95) {
  p   <- p / sum(p)
  ord <- order(p, decreasing = TRUE)
  cut <- p[ord][which(cumsum(p[ord]) >= mass)[1]]
  yr  <- years[p >= cut]
  c(start = min(yr), end = max(yr))
}

#' Case 1: what the point-date models read from each calibrated date.
#'
#' start, end  the 95% HPD range (the midpoint model uses its centre)
#' median      the calibrated median
#' mass_kept   share of the date's probability inside that range
calibrated_summaries <- function(cal, grid, mass = 0.95, max_loss = 0.01) {
  m <- calibrated_on_grid(cal, grid, max_loss)

  out <- vapply(seq_len(ncol(m)), function(j) {
    p <- m[, j] / sum(m[, j])
    e <- hpd_envelope(p, grid, mass)
    c(e,
      median    = grid[which(cumsum(p) >= 0.5)[1]],
      mass_kept = sum(p[grid >= e["start"] & grid <= e["end"]]))
  }, numeric(4))

  data.frame(start = out[1, ], end = out[2, ], median = out[3, ], mass_kept = out[4, ])
}
