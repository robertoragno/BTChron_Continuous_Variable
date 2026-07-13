# Base-R sandbox: does the choice of midpoint vs true-date change the recovered
# rate of change, and under which data-generating process? A fast lm-based
# exploration to run before committing any Stan time. Nothing here fits Stan.
#
# Archaeological setup: each sample carries a phase [start, end] but no exact
# date. A quantity (say settlement size) trends with time and we want the slope.
# Two slopes are compared:
#   true-date - lm(Value ~ True_date): uses the real dates, available only in
#              simulation. The best any method could recover; the target.
#   midpoint - lm(Value ~ phase midpoint): the naive baseline, and exactly what
#              the midpoint Stan model does.
# The EIV / latent Stan model cannot be reproduced by lm, but its slope tracks
# the true-date fit whenever the midpoint one does, and its real payoff is honest
# intervals, not a different point estimate. So the true-date-vs-midpoint gap here
# is the quantity the latent model exists to close.
#
# Four data-generating processes for the true dates, all over [TMIN, TMAX]:
#   uniform  - dates flat over the timeline (current simulate.R). Within each
#              window the dates are uniform, hence symmetric about the midpoint.
#   counts   - phases populated unevenly, but dates still uniform *within* each
#              window. Only the per-phase counts change; within-window shape
#              stays symmetric.
#   skew     - dates drawn from a left-skewed timeline distribution, then binned
#              into fixed phases. Within each window the dates pile toward the
#              early edge: non-uniform and asymmetric about the midpoint.
#   skew_late - same skew, but the widest phases sit late on the timeline, so the
#              asymmetric shift grows with time instead of staying constant.
#
# Prediction under Berkson-type error: midpoint recovers the slope whenever the
# within-window dates are symmetric about the midpoint (uniform, counts), and
# only departs when they are skewed (skew), most on the slope when that skew is
# position-dependent (skew_late).

library(here)
source(here("Simulations", "Sim_Linear", "scripts", "simulate.R"))  # TMIN, TMAX

set.seed(2026)

intercept <- 5
slope     <- 0.015   # cities get bigger with time
sigma     <- 1.5

# Fixed broken-stick phase boundaries. coarse_end places the widest phases early,
# late, at both ends, or leaves them in Dirichlet order (even).
make_boundaries <- function(K, alpha_conc, coarse_end = "even") {
  w <- rgamma(K, shape = alpha_conc, rate = 1)
  w <- w / sum(w)
  if (coarse_end == "late")  w <- sort(w)                 # widest phase last
  if (coarse_end == "early") w <- sort(w, decreasing = TRUE)
  if (coarse_end == "both") {                             # widest at both ends
    ord <- order(w, decreasing = TRUE)
    slot <- integer(K); ends <- c(1, K); mid <- 2:(K - 1)
    slot[c(ends, mid)] <- ord                             # 2 widest -> ends
    w <- w[slot]
  }
  spans <- w * (TMAX - TMIN)
  b <- round(TMIN + c(0, cumsum(spans)))
  b[1] <- TMIN
  b[length(b)] <- TMAX
  b
}

# Bin true dates into the fixed phases and attach each phase's window.
attach_windows <- function(true_date, boundaries) {
  ph <- findInterval(true_date, boundaries, rightmost.closed = TRUE,
                     all.inside = TRUE)
  data.frame(Start_date = boundaries[ph],
             End_date   = boundaries[ph + 1],
             True_date  = true_date)
}

# Draw one dataset under a named DGP. Returns Start_date, End_date, True_date,
# Value.
simulate_dgp <- function(dgp, N, K, alpha_conc) {
  if (dgp == "counts") {
    # uneven phase counts, uniform within each window (symmetric)
    boundaries    <- make_boundaries(K, alpha_conc)
    phase_weights <- sort(rgamma(K, shape = 0.5, rate = 1), decreasing = TRUE)
    ph <- sample(K, N, replace = TRUE, prob = phase_weights)
    a  <- boundaries[ph]; b <- boundaries[ph + 1]
    df <- data.frame(Start_date = a, End_date = b,
                     True_date  = round(runif(N, a, b)))
  } else {
    coarse_end <- switch(dgp, skew_late = "late", skew_both = "both", "even")
    boundaries <- make_boundaries(K, alpha_conc, coarse_end)
    u <- switch(dgp,
                uniform   = runif(N),               # symmetric within window
                skew      = rbeta(N, 1.5, 5),        # left-skewed timeline
                skew_late = rbeta(N, 1.5, 5),
                skew_both = rbeta(N, 1.5, 5))
    true_date <- round(TMIN + u * (TMAX - TMIN))
    df <- attach_windows(true_date, boundaries)
  }
  df$Value <- round(intercept + slope * df$True_date + rnorm(N, 0, sigma), 1)
  df
}

# True-date and midpoint slope + intercept for one dataset. A degenerate midpoint
# (all observations in one phase, no spread to regress on) returns NA.
estimates <- function(df) {
  mid <- (df$Start_date + df$End_date) / 2
  truth <- coef(lm(df$Value ~ df$True_date))
  mfit   <- if (length(unique(mid)) > 1) coef(lm(df$Value ~ mid)) else c(NA, NA)
  c(truedate_slope = unname(truth[2]), truedate_int = unname(truth[1]),
    mid_slope    = unname(mfit[2]),   mid_int    = unname(mfit[1]))
}

# Repeat each DGP many times and summarise where midpoint lands relative to the
# true-date fit and the truth. More, more even phases (higher K, higher alpha_conc)
# so the skewed early cluster still gets subdivided.
dgps <- c("uniform", "counts", "skew", "skew_late", "skew_both")
n_rep <- 500
N <- 300; K <- 8; alpha_conc <- 1.5

summary_rows <- lapply(dgps, function(dgp) {
  est <- t(replicate(n_rep, estimates(simulate_dgp(dgp, N, K, alpha_conc))))
  n_ok <- sum(!is.na(est[, "mid_slope"]))
  data.frame(
    DGP            = dgp,
    n_fit          = n_ok,
    midpoint_slope = mean(est[, "mid_slope"], na.rm = TRUE),
    slope_bias     = mean(est[, "mid_slope"], na.rm = TRUE) - slope,
    slope_sd       = sd(est[, "mid_slope"],   na.rm = TRUE),   # precision
    intercept_bias = mean(est[, "mid_int"],   na.rm = TRUE) - intercept
  )
})
results <- do.call(rbind, summary_rows)

cat("\nMidpoint recovery by data-generating process",
    sprintf("(N = %d, K = %d phases, %d replicates; true slope %.3f, intercept %.0f)\n\n",
            N, K, n_rep, slope, intercept))
print(format(results, digits = 3, scientific = FALSE), row.names = FALSE)

cat("\nmidpoint stays on the true slope while within-window dates are symmetric",
    "(uniform, counts);\neven skew mostly shifts the intercept; only",
    "position-dependent skew (skew_late) bends the slope.\n")

# Overconfidence: one periodisation held fixed, only dates + noise resampled.
# Deposition is uniform within each window, so the slope stays unbiased and only
# the honesty of its reported uncertainty is on trial. The same set of phase
# widths is arranged two ways: widest phases in the middle (low regression
# leverage) versus at both ends (high leverage). The midpoint lm folds within-
# window date spread into a single homoskedastic noise term; when the widest,
# noisiest windows sit at the high-leverage ends, that flat SE underweights them
# and understates the true slope spread. The true-date fit is the control.

boundaries_from <- function(weights, where = c("ends", "middle")) {
  where <- match.arg(where)
  K <- length(weights)
  prio <- as.vector(rbind(seq_len(K), rev(seq_len(K))))[seq_len(K)]  # 1,K,2,K-1,..
  if (where == "middle") prio <- rev(prio)                           # widest centred
  slot <- integer(K)
  slot[prio] <- order(weights, decreasing = TRUE)
  spans <- weights[slot] / sum(weights) * (TMAX - TMIN)
  b <- round(TMIN + c(0, cumsum(spans)))
  b[1] <- TMIN
  b[K + 1] <- TMAX
  b
}

# One resample from a fixed periodisation: slope estimate and its reported SE,
# for both the true-date and the midpoint fits.
fixed_fit <- function(boundaries, N) {
  df <- attach_windows(round(runif(N, TMIN, TMAX)), boundaries)
  df$Value <- round(intercept + slope * df$True_date + rnorm(N, 0, sigma), 1)
  mid <- (df$Start_date + df$End_date) / 2
  orc <- summary(lm(df$Value ~ df$True_date))$coefficients[2, 1:2]
  mp  <- if (length(unique(mid)) > 1)
           summary(lm(df$Value ~ mid))$coefficients[2, 1:2] else c(NA, NA)
  c(truedate_est = unname(orc[1]), truedate_se = unname(orc[2]),
    mid_est    = unname(mp[1]),  mid_se    = unname(mp[2]))
}

summarise_oc <- function(fits, prefix, periodisation, model) {
  est <- fits[, paste0(prefix, "_est")]; se <- fits[, paste0(prefix, "_se")]
  data.frame(periodisation = periodisation, model = model,
             true_slope_sd  = sd(est, na.rm = TRUE),     # actual spread
             reported_se    = mean(se, na.rm = TRUE),    # what the fit claims
             overconfidence = sd(est, na.rm = TRUE) / mean(se, na.rm = TRUE),
             ci95_coverage  = mean(abs(est - slope) <= 1.96 * se, na.rm = TRUE))
}

run_oc <- function(boundaries, periodisation, n_rep = 2000, N = 300) {
  fits <- t(replicate(n_rep, fixed_fit(boundaries, N)))
  rbind(summarise_oc(fits, "truedate", periodisation, "true dates"),
        summarise_oc(fits, "mid",    periodisation, "midpoint"))
}

set.seed(7)
phase_weights <- rgamma(8, shape = 0.7, rate = 1)   # spread of narrow to wide

oc <- rbind(run_oc(boundaries_from(phase_weights, "middle"), "coarse in middle"),
            run_oc(boundaries_from(phase_weights, "ends"),   "coarse at both ends"))

cat("\nSlope uncertainty: reported SE vs true spread, one fixed periodisation",
    "\n(same phase widths, arranged two ways; overconfidence > 1 = too sure)\n\n")
print(format(oc, digits = 3), row.names = FALSE)

# Diagnostic sampling: sample intensity coupled to phase width, always uniform
# within the window (Stan prior intact -- this changes only the sampling, not the
# machinery). Real assemblages need not sample phases in proportion to their
# width: a short, diagnostic phase (recognisable ceramics) can be over-sampled.
# When the narrow, well-dated phases carry most of the sample, the within-window
# date scatter that midpoint mistakes for noise shrinks, so the room EIV has to
# improve on midpoint (its edge in sigma and in slope precision) shrinks too.
# The true-date fit is the target EIV would approach.

sample_regime <- function(boundaries, N, regime) {
  K <- length(boundaries) - 1
  width <- diff(boundaries)
  w <- switch(regime,
              proportional = width,        # samples prop to width (uniform global)
              flat         = rep(1, K),     # equal count per phase
              diagnostic   = 1 / width)     # narrow phases over-sampled
  ph <- sample(K, N, replace = TRUE, prob = w)
  a <- boundaries[ph]; b <- boundaries[ph + 1]
  df <- data.frame(Start_date = a, End_date = b,
                   True_date  = round(runif(N, a, b)))
  df$Value <- round(intercept + slope * df$True_date + rnorm(N, 0, sigma), 1)
  df
}

# One resample: mean window width actually sampled, and how far midpoint sits
# from the true-date fit on sigma and on slope-interval width (the room EIV has).
gap_fit <- function(df) {
  mid <- (df$Start_date + df$End_date) / 2
  if (length(unique(mid)) < 2) return(rep(NA, 4))
  om <- summary(lm(df$Value ~ df$True_date))
  mm <- summary(lm(df$Value ~ mid))
  c(sampled_width  = mean(df$End_date - df$Start_date),
    mid_slope      = mm$coefficients[2, 1],
    mid_sigma_bias = mm$sigma - sigma,
    width_ratio    = mm$coefficients[2, 2] / om$coefficients[2, 2])  # mid / true-date
}

set.seed(11)
gap_boundaries <- boundaries_from(rgamma(8, shape = 0.7, rate = 1), "ends")

gap <- do.call(rbind, lapply(c("proportional", "flat", "diagnostic"), function(rg) {
  fits <- t(replicate(1000, gap_fit(sample_regime(gap_boundaries, 300, rg))))
  data.frame(sampling = rg,
             sampled_width  = mean(fits[, "sampled_width"], na.rm = TRUE),
             mid_slope      = mean(fits[, "mid_slope"], na.rm = TRUE),
             mid_sigma_bias = mean(fits[, "mid_sigma_bias"], na.rm = TRUE),
             width_ratio    = mean(fits[, "width_ratio"], na.rm = TRUE))
}))

cat("\nDiagnostic sampling: EIV-vs-midpoint room by sampling intensity",
    "\n(one fixed periodisation; uniform within window; true sigma", sigma, ")\n\n")
print(format(gap, digits = 3), row.names = FALSE)
cat("\nover-sampling the narrow diagnostic phases pulls the mean sampled window",
    "down;\nmidpoint's sigma inflation and its slope-interval excess over the",
    "true-date fit\nboth shrink -- the EIV advantage narrows as the sample gets better",
    "dated.\n")

# Visual: midpoint slope against true-date slope per replicate, one panel per DGP.
png(here("Simulations", "Sim_Linear", "figures", "dgp_sandbox_slopes.png"),
    width = 3000, height = 700, res = 300)
op <- par(mfrow = c(1, 5), mar = c(4, 4, 3, 1))
for (dgp in dgps) {
  est <- t(replicate(300, estimates(simulate_dgp(dgp, N, K, alpha_conc))))
  plot(est[, "truedate_slope"], est[, "mid_slope"], pch = 16, cex = 0.4,
       col = "#78000055", xlab = "true-date slope", ylab = "midpoint slope",
       main = dgp)
  abline(0, 1, col = "grey40", lty = 2)
  abline(v = slope, col = "grey70", lty = 3)
}
par(op)
dev.off()
