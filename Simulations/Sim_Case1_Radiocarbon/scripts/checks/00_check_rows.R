# Checks that a calibrated date keeps its calendar order when it becomes a weight
# row. Getting the order wrong would flip the sign of the slope without any error,
# so the row's median is compared with rcarbon's own calibrated median.

library(here)
suppressMessages(library(rcarbon))

source(here("Simulations", "shared", "scripts", "weight_rows.R"))

# Two dates, deliberately different: one on the Hallstatt plateau, where the
# posterior is wide and multimodal, one off it.
cra    <- c(2450, 2000)
errors <- c(25, 25)

x   <- calibrate(x = cra, errors = errors, calCurves = "intcal20",
                 calMatrix = TRUE, verbose = FALSE)
cal <- calmatrix_to_calendar(x)

cat("calMatrix rows: cal BP", rownames(x$calmatrix)[1], "->",
    tail(rownames(x$calmatrix), 1), "\n")
cat("after conversion: ", cal$years[1], "->", tail(cal$years, 1),
    " BCE/CE, ascending\n\n")

# rcarbon's own answer, converted to BCE/CE for comparison.
median_bcad <- 1950 - summary(x)$MedianBP

for (step in c(1, 5)) {
  # Padded well past both dates: calibrated_rows() errors on truncation, and
  # truncating here would test the padding rather than the indexing.
  grid <- seq(-1400, 600, by = step)
  w    <- calibrated_rows(cal, grid)

  for (j in seq_along(cra)) {
    p    <- exp(w$log_year_prob_packed[seq(sum(w$row_n_years[seq_len(j - 1)]) + 1, length.out = w$row_n_years[j])])
    yrs  <- grid[w$row_first_year[j] + seq_along(p) - 1]
    peak <- yrs[which.max(p)]
    med  <- yrs[which.max(cumsum(p) >= 0.5)]

    cat(sprintf(
      "%d yr grid  %d+-%d BP: row spans %d yr (%d..%d), peak %d, row median %d, rcarbon median %d\n",
      step, cra[j], errors[j], w$row_n_years[j], min(yrs), max(yrs), peak, med,
      round(median_bcad[j])))

    # The peak can sit a little off the median on a skewed or multimodal
    # posterior, so the row median is the strict test; a sign flip would put
    # both hundreds of years out, and on the wrong side.
    if (abs(med - median_bcad[j]) > 3 * step + 5)
      stop("row median disagrees with rcarbon - check the row order")
  }
  cat("\n")
}

cat("row order and grid indexing verified\n")
