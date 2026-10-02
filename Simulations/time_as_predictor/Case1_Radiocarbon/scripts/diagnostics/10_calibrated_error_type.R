# Which kind of error does a calibrated date carry?
#
# A radiocarbon measurement has classical error, but a calibrated date is a
# posterior, and its summaries can behave either way. With e = summary - true:
#   classical: e unrelated to the true date, summaries MORE spread than truth
#   Berkson:   e unrelated to the summary, summaries LESS spread than truth
#
# The last two columns are the slope ratios a regression on that summary gives,
# with the date as the predictor (cov / var of the summary) and as the response
# (cov / var of the true date). No model is fitted: these are the expected
# values, and they reproduce the fitted results in the README and in
# 08_reversed_berkson.R. 1000 dates per window, Case 1's own generator.
#
#   Rscript Simulations/time_as_predictor/Case1_Radiocarbon/scripts/diagnostics/10_calibrated_error_type.R

suppressMessages(library(rcarbon))
library(here)
source(here("Simulations", "shared", "scripts", "weight_rows.R"))
source(here("Simulations", "shared", "dating", "case1_radiocarbon.R"))

windows <- list(steep = c(-1600, -1200), plateau = c(-800, -400))

for (w in names(windows)) {
  win <- windows[[w]]
  sim <- simulate_c14(1000, 8, 0.02, 2, win, 30, seed = 1)

  # Calibrate over the window plus 900 yr each side, which holds every date whole
  bp  <- sort(1950 - (win + c(-900, 900)), decreasing = TRUE)
  cal <- calibrate(sim$CRA, errors = sim$Error, calMatrix = TRUE,
                   timeRange = bp, verbose = FALSE)
  m     <- cal$calmatrix
  years <- 1950 - as.numeric(rownames(m))
  m     <- sweep(m, 2, colSums(m), "/")
  inwin <- years >= win[1] & years <= win[2]

  # The last summary is the posterior mean under the right prior (the window),
  # which is what the model would give if it were told where the dates lie
  summaries <- list(
    "midpoint of 95% range" = apply(m, 2, function(p) mean(hpd_envelope(p, years))),
    "calibrated median"     = apply(m, 2, function(p) years[which(cumsum(p) >= 0.5)[1]]),
    "calibrated mean"       = colSums(m * years),
    "mean, window prior"    = apply(m, 2, function(p) { q <- p * inwin; sum(q * years) / sum(q) })
  )

  t <- sim$True_date
  cat("\n", w, "\n", sep = "")
  cat(sprintf("%-24s %8s %11s %10s %10s %10s\n", "summary", "sd ratio",
              "cor(e,true)", "cor(e,obs)", "as x", "as y"))
  for (s in names(summaries)) {
    o <- summaries[[s]]
    e <- o - t
    cat(sprintf("%-24s %8.2f %11.2f %10.2f %10.2f %10.2f\n", s, sd(o) / sd(t),
                cor(e, t), cor(e, o), cov(o, t) / var(o), cov(o, t) / var(t)))
  }
}
