# Builds the ragged weight rows that marginal_date.stan consumes.
#
# Every dating case reduces to the same object: a shared annual grid of candidate
# calendar years, and one probability row per sample over that grid. The cases
# differ only in what the row looks like - flat over a window for typochronology
# and phases, the calibrated posterior for radiocarbon.
#
# Rows are stored trimmed to their support, since a 100-yr window on a 1500-yr
# axis is 93% zeros. Each row is one contiguous span: grid index of its first
# year, its length, and the log weights concatenated into one flat vector.
#
# The grid index to calendar year map is affine and lives in one place:
#   year = grid_min_year + (index - 1) * grid_step
# grid_min_year is the constant that has to travel with the matrix. Losing it,
# or reading a row order backwards, flips the recovered slope with no error
# anywhere, so grid_year is passed to Stan explicitly rather than rebuilt there.

#' Shared candidate-year grid for a case.
#'
#' @param time_ref_min,time_ref_range the calendar axis reference from the design
#' @param step grid resolution in years (1 = annual)
make_grid <- function(time_ref_min, time_ref_range, step = 1) {
  seq(time_ref_min, time_ref_min + time_ref_range, by = step)
}

#' Pack a list of (first grid index, probability vector) into ragged form.
#'
#' Rows are renormalised here, once, after any trimming the caller has done.
#' rcarbon normalises over its whole timeRange, so a row filtered to the study
#' window no longer sums to 1, and the model rejects rows that do not.
pack_rows <- function(rows, T_grid) {
  row_first_year <- as.integer(vapply(rows, function(r) r$index, numeric(1)))
  probs   <- lapply(rows, function(r) r$prob / sum(r$prob))
  row_n_years   <- vapply(probs, length, integer(1))

  if (any(row_first_year < 1L) || any(row_first_year + row_n_years - 1L > T_grid))
    stop("a weight row falls outside the grid")
  if (any(!vapply(probs, function(p) all(is.finite(p)) && all(p >= 0), logical(1))))
    stop("weight rows must be finite and non-negative")

  list(n_weights = sum(row_n_years), log_year_prob_packed = log(unlist(probs, use.names = FALSE)),
       row_first_year = row_first_year, row_n_years = row_n_years)
}

#' Case 2/3/4 rows: flat over each sample's window.
#'
#' Windows are given on the calendar scale and snapped to the grid cells they
#' cover, so a window narrower than one grid step still gets a single cell
#' rather than an empty row.
uniform_rows <- function(start_date, end_date, grid) {
  step <- if (length(grid) > 1) grid[2] - grid[1] else 1
  i_lo <- pmax(1L, ceiling((start_date - grid[1]) / step) + 1L)
  i_hi <- pmin(length(grid), floor((end_date - grid[1]) / step) + 1L)
  i_hi <- pmax(i_hi, i_lo)

  rows <- Map(function(lo, hi) list(index = lo, prob = rep(1, hi - lo + 1L)),
              i_lo, i_hi)
  pack_rows(rows, length(grid))
}

#' Highest-density envelope of a probability row: the smallest single interval
#' spanning every year above the density threshold that holds `mass`.
#'
#' The latent and midpoint models need one [start, end] per date, so a multimodal
#' calibrated posterior has to be collapsed. The gaps between humps are swallowed
#' - that loss is the point of the comparison, and its size is reported.
hpd_envelope <- function(p, years, mass = 0.95) {
  p   <- p / sum(p)
  ord <- order(p, decreasing = TRUE)
  cut <- p[ord][which(cumsum(p[ord]) >= mass)[1]]
  yr  <- years[p >= cut]
  c(start = min(yr), end = max(yr))
}

#' Envelopes and medians for every column of a calibrated matrix on the grid.
#' Returns what the point-date and latent models consume, alongside the weight rows.
calibrated_summaries <- function(cal, grid, mass = 0.95, max_loss = 0.01) {
  step <- if (length(grid) > 1) grid[2] - grid[1] else 1
  cell <- round((cal$years - grid[1]) / step) + 1L
  keep <- cell >= 1L & cell <= length(grid)
  m    <- rowsum(cal$prob[keep, , drop = FALSE], cell[keep], reorder = TRUE)
  full <- matrix(0, length(grid), ncol(cal$prob))
  full[as.integer(rownames(m)), ] <- m

  # Same guard as calibrated_rows(): a truncated date gives a summary that is
  # pulled inward, silently, rather than an obviously wrong one.
  lost <- 1 - colSums(full)
  if (any(lost > max_loss))
    stop(sprintf("date(s) lose up to %.1f%% of their calibrated mass off the grid; widen the padding",
                 100 * max(lost)))

  out <- vapply(seq_len(ncol(full)), function(j) {
    p <- full[, j] / sum(full[, j])
    e <- hpd_envelope(p, grid, mass)
    c(e, median = grid[which(cumsum(p) >= 0.5)[1]],
      mass_kept = sum(p[grid >= e["start"] & grid <= e["end"]]))
  }, numeric(4))

  data.frame(start = out[1, ], end = out[2, ],
             median = out[3, ], mass_kept = out[4, ])
}

#' rcarbon's calMatrix, turned into an ascending BCE/CE matrix.
#'
#' calMatrix is years x dates with its rows labelled in cal BP DESCENDING: row 1
#' is 55000 BP, the oldest year, and the last row is 0 BP. Cal BP counts
#' backwards, so descending BP is already ASCENDING calendar order - the rows
#' need relabelling, not reordering. Reversing them would be the bug, not the
#' fix.
#'
#' The conversion is here and nowhere else on purpose. Relabelling the years
#' without matching the row order passes every assertion downstream and flips
#' the sign of the recovered slope, so the years are taken from the rownames
#' rather than supplied by the caller, and the ascending check below is the
#' guard that catches it.
calmatrix_to_calendar <- function(cal_dates) {
  m <- cal_dates$calmatrix
  if (is.null(m))
    stop("calibrate() must be called with calMatrix = TRUE")

  bp    <- as.numeric(rownames(m))
  if (anyNA(bp)) stop("calMatrix rownames are not cal BP years")
  years <- 1950 - bp

  if (is.unsorted(years))
    stop("calMatrix rows are not in ascending calendar order after relabelling; ",
         "rcarbon's row order has changed and this conversion no longer holds")

  list(prob = m, years = years)
}

#' Case 1 rows: the calibrated posterior on the shared grid.
#'
#' Trimmed twice: to the study grid, then to each date's own support.
#'
#' @param cal        output of calmatrix_to_calendar()
#' @param grid       the shared grid from make_grid(), ascending BCE/CE
#' @param mass_keep  keep the smallest contiguous span holding this much of the
#'                    calibrated mass. Trimming the tails is what makes the
#'                    ragged storage pay off; the threshold is a reported choice.
#' @param max_loss  error if any date loses more than this share of its mass off
#'                    the edge of the grid
calibrated_rows <- function(cal, grid, mass_keep = 0.9999, max_loss = 0.01) {
  cal_matrix <- cal$prob
  cal_years  <- cal$years
  if (nrow(cal_matrix) != length(cal_years))
    stop("cal_matrix rows and cal_years must correspond")
  if (is.unsorted(cal_years))
    stop("cal_years must be ascending, matching the grid")

  #' Shared grid. On a grid that is coarser than annual, the calibrated years
  #' are binned, not subsetted: i.e. taking every fifth year would throw away 
  #' four fifths of the mass and would resample a multimodal posterior rather than
  #' summarise it. Each grid cell collects the annual probabilities within half a
  #' step of it.
  step <- if (length(grid) > 1) grid[2] - grid[1] else 1 # else 1 is just a safeguard
  if (min(cal_years) > min(grid) - step / 2 ||
      max(cal_years) < max(grid) + step / 2)
    stop("the calibration range does not cover the study grid")

  cell <- round((cal_years - grid[1]) / step) + 1L # Creates an index of the grid cell 
  # each year falls into, 1-based. It's a many-to-one mapping, so for a 5 years grid it assigns the same index for each of the 5 yars
  inside <- cell >= 1L & cell <= length(grid)
  m <- rowsum(cal_matrix[inside, , drop = FALSE], cell[inside], reorder = TRUE)
  # rowsum() drops empty cells; put them back so row i is grid year i.
  full <- matrix(0, nrow = length(grid), ncol = ncol(cal_matrix))
  full[as.integer(rownames(m)), ] <- m
  m <- full

  # DECISION: the study window is NOT used as prior information.
  #
  # rcarbon calibrates against a uniform prior over calendar years, and the row
  # is used exactly as it comes back. The model is never told that the true
  # dates were drawn from a 400-yr window, because a real analyst does not know
  # that either - telling it would flatter the method with knowledge nobody has.
  #
  # This check is where that decision is enforced. Clipping a row at the edge of
  # the grid and renormalising would impose a window prior by the back door, and
  # a broken one: asserted by a grid boundary rather than stated as a prior, and
  # dragging edge dates inward, in the same direction as the attenuation the
  # study measures. So truncation is an error rather than something silently
  # absorbed, and the grid is padded (01_design.R) so it never fires. A warning
  # would scroll past unseen in a 1400-dataset loop.
  lost <- 1 - colSums(m)
  if (any(lost > max_loss))
    stop(sprintf("date(s) lose up to %.1f%% of their calibrated mass off the grid; widen the padding",
                 100 * max(lost)))

  rows <- lapply(seq_len(ncol(m)), function(j) {
    p <- m[, j]
    if (sum(p) <= 0) stop("date ", j, " has no probability mass on the grid")
    p <- p / sum(p)
    # Contiguous span rather than the non-zero cells: a plateau posterior is
    # multimodal with near-zero gaps, and one span keeps the indexing simple.
    drop <- (1 - mass_keep) / 2
    cum  <- cumsum(p)
    lo   <- max(1L, which(cum >= drop)[1])
    hi   <- min(length(p), which(cum >= 1 - drop)[1])
    list(index = lo, prob = p[lo:hi])
  })

  pack_rows(rows, length(grid))
}
