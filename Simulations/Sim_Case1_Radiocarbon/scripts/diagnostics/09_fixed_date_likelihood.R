# What happens to the likelihood when one find's date is fixed to a constant
#
# The question: take one find out of the dataset, replace its calibrated
# distribution with a single year, do the same with a second year, and ask
# whether the total likelihood of the model differs between the two.
#
# In marginal_date.stan the date is not a parameter. Each find enters through a
# row of probabilities over the year grid, which is data, so fixing a date means
# handing the model a row one year long instead of the calibrated one. Nothing
# in the Stan file changes: a one-year row still sums to 1 and row_n_years
# allows a length of 1. The parameter space is the same three parameters either
# way, so the two runs are directly comparable and the usual objection about
# comparing lp__ across models of different size does not apply here.
#
# Because a row is normalised, collapsing find n to year t makes its
# contribution to the likelihood exactly
#
#   Normal(y_n | alpha + beta * t, sigma)
#
# and every other find's contribution is untouched. So the difference between
# two constants is that one term evaluated at the two years, and sweeping t
# across the grid traces what the trend alone says about when the find dates to.
# The sweep is computed on the draws of the unmodified fit, so the two years are
# compared at identical alpha, beta and sigma rather than across two samplers.
# The two refits are there to confirm the sweep and to show how little the trend
# moves when one find out of 200 is pinned.
#
# The find is picked as the most multi-peaked calibrated date in the dataset.
# On the Hallstatt plateau a date can have its probability split over several
# humps, and that is where fixing it to one year has something to say.
#
#   Rscript Simulations/Sim_Case1_Radiocarbon/scripts/diagnostics/09_fixed_date_likelihood.R

library(here)
library(cmdstanr)
suppressMessages(library(rcarbon))

source(here("Simulations", "shared", "scripts", "weight_rows.R"))
source(here("Simulations", "shared", "scripts", "fit_models.R"))
source(here("Simulations", "Sim_Case1_Radiocarbon", "scripts", "simulate.R"))

out_dir <- here("Simulations", "Sim_Case1_Radiocarbon", "output", "fixed_date")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# Study settings, the reference cell of the main design ----
WINDOW    <- c(-800, -400)   # Hallstatt plateau
LAB_ERROR <- 30
N         <- 200
GRID_STEP <- 5
INTERCEPT <- 8
SLOPE     <- 0.015
SIGMA     <- 1.5
OFFSET    <- 200             # how far the second constant sits from the first

# One dataset ----
sim  <- simulate_c14(N, INTERCEPT, SLOPE, SIGMA, WINDOW, LAB_ERROR, seed = 9)
PAD  <- c14_grid_pad(WINDOW, LAB_ERROR)
grid <- seq(WINDOW[1] - PAD, WINDOW[2] + PAD, by = GRID_STEP)

cal <- calmatrix_to_calendar(
  calibrate(sim$CRA, errors = sim$Error, calCurves = "intcal20",
            calMatrix = TRUE, verbose = FALSE))
w <- calibrated_rows(cal, grid)

# Unpack the rows so one of them can be replaced, then repack with pack_rows
unpack <- function(w) {
  off <- 0L
  lapply(seq_along(w$row_n_years), function(i) {
    k <- w$row_n_years[i]
    p <- exp(w$log_year_prob_packed[off + seq_len(k)])
    off <<- off + k
    list(index = w$row_first_year[i], prob = p)
  })
}
rows <- unpack(w)

# The find with the most peaks in its calibrated date ----
n_peaks <- vapply(rows, function(r) {
  p <- r$prob / max(r$prob)
  sum(p[-c(1, length(p))] > 0.2 & diff(p)[-length(diff(p))] > 0 & diff(p)[-1] < 0)
}, numeric(1))
target <- which.max(n_peaks)
cat(sprintf("Sample %d has %d peaks, true date %d\n",
            target, n_peaks[target], sim$True_date[target]))

# Fit the model as the study runs it ----
model <- cmdstan_model(here("Simulations", "shared", "models", "marginal_date.stan"))

stan_data <- model_stan_data(
  "marginal", list(grid = grid, weights = w), sim$Value,
  WINDOW[1] - PAD, diff(WINDOW) + 2 * PAD, WINDOW)

fit <- model$sample(data = stan_data, chains = 4, parallel_chains = 4,
                    iter_warmup = 1000, iter_sampling = 1000,
                    adapt_delta = 0.95, refresh = 0)

draws <- fit$draws(format = "df")
alpha <- draws$alpha; beta <- draws$beta; sigma <- draws$sigma
ll_free <- rowSums(as.matrix(draws[, grep("^log_lik\\[", colnames(draws))]))

# Sweeping the constant across the grid ----
# The contribution of the target find if its date were fixed to year t, one
# value per posterior draw. Computed on the rescaled axis the model uses.
grid_norm <- 2 * (grid - stan_data$time_ref_min) / stan_data$time_ref_range - 1
y_t <- sim$Value[target]

sweep <- sapply(grid_norm, function(t)
  dnorm(y_t, alpha + beta * t, sigma, log = TRUE))
sweep_mean <- colMeans(sweep)

# The find's own contribution when its date is left free, for reference
ll_target_free <- draws[[sprintf("log_lik[%d]", target)]]

best  <- grid[which.max(sweep_mean)]
# The true date is included because the claim being tested is whether the year
# the find really dates to is the one the likelihood prefers. It is rounded to
# the grid, like every other constant here.
true_on_grid <- grid[which.min(abs(grid - sim$True_date[target]))]
years <- c(best, true_on_grid, best - OFFSET, best + OFFSET)

cat(sprintf("Best constant on the sweep: %d (true date %d, calibrated median %.0f)\n",
            best, sim$True_date[target],
            calibrated_summaries(cal, grid)$median[target]))

sweep_out <- data.frame(year = grid, mean_log_lik = sweep_mean,
                        lo = apply(sweep, 2, quantile, 0.05),
                        hi = apply(sweep, 2, quantile, 0.95))
write.csv(sweep_out, file.path(out_dir, "likelihood_sweep.csv"), row.names = FALSE)

# Refitting with the date fixed ----
refit <- function(year) {
  rows_fixed <- rows
  rows_fixed[[target]] <- list(index = which.min(abs(grid - year)), prob = 1)
  d <- stan_data
  d[c("n_weights", "log_year_prob_packed", "row_first_year", "row_n_years")] <-
    pack_rows(rows_fixed, length(grid))[
      c("n_weights", "log_year_prob_packed", "row_first_year", "row_n_years")]

  f <- model$sample(data = d, chains = 4, parallel_chains = 4,
                    iter_warmup = 1000, iter_sampling = 1000,
                    adapt_delta = 0.95, refresh = 0)
  dr <- f$draws(format = "df")
  data.frame(
    fixed_year = year,
    total_log_lik = mean(rowSums(as.matrix(dr[, grep("^log_lik\\[", colnames(dr))]))),
    find_log_lik  = mean(dr[[sprintf("log_lik[%d]", target)]]),
    slope = mean(dr$slope_original), slope_sd = sd(dr$slope_original),
    sigma = mean(dr$sigma), sigma_sd = sd(dr$sigma))
}

fixed <- rbind(
  data.frame(fixed_year = NA, total_log_lik = mean(ll_free),
             find_log_lik = mean(ll_target_free),
             slope = mean(draws$slope_original),
             slope_sd = sd(draws$slope_original),
             sigma = mean(draws$sigma), sigma_sd = sd(draws$sigma)),
  do.call(rbind, lapply(years, refit)))
fixed$true_slope <- SLOPE
fixed$true_sigma <- SIGMA
fixed$true_date  <- sim$True_date[target]
fixed$delta_total <- fixed$total_log_lik - fixed$total_log_lik[1]
write.csv(fixed, file.path(out_dir, "fixed_date.csv"), row.names = FALSE)
print(fixed)

# Figure ----
# The calibrated date, what the trend says about it, and what the model ends up
# believing, all scaled to a maximum of 1 so the shapes can be compared.
row_target  <- rows[[target]]
row_years   <- grid[row_target$index + seq_along(row_target$prob) - 1]
sweep_prob  <- exp(sweep_mean - max(sweep_mean))
date_actual <- draws[[sprintf("date_actual[%d]", target)]]
post <- density(date_actual, from = min(grid), to = max(grid), n = 1024)

xlim <- range(row_years) + c(-100, 100)
# the right panel also has to hold the three constants
xlim_fixed <- range(c(xlim, years)) + c(-50, 50)

png(file.path(out_dir, "fixed_date.png"), width = 10, height = 4.5,
    units = "in", res = 300)
par(mfrow = c(1, 2), mar = c(4, 4, 3, 1))

plot(NA, xlim = xlim, ylim = c(0, 1.05), xlab = "calendar year (BCE/CE)",
     ylab = "scaled density", main = sprintf("Sample %d", target))
polygon(c(row_years, rev(row_years)),
        c(row_target$prob / max(row_target$prob), rep(0, length(row_years))),
        border = NA, col = adjustcolor("grey60", 0.5))
lines(grid, sweep_prob, lwd = 2, lty = 2, col = "grey25")
lines(post$x, post$y / max(post$y), lwd = 2, col = "#780000")
abline(v = sim$True_date[target], lty = 3)
legend("topright", bty = "n", cex = 0.8,
       legend = c("Calibrated date", "Likelihood of the sample at each year",
                  "Model's date for the sample", "True date"),
       lwd = c(8, 2, 2, 1), lty = c(1, 2, 1, 3),
       col = c(adjustcolor("grey60", 0.5), "grey25", "#780000", "black"))

# The parabola runs down to about -45 at the edges of the grid, which squashes
# the part that matters against the top of the panel, so the y axis is cut to
# the years actually shown, with headroom for the labels.
sel  <- grid >= xlim_fixed[1] & grid <= xlim_fixed[2]
ylim <- range(sweep_mean[sel]) + c(-0.5, 1.5)

plot(grid, sweep_mean, type = "l", lwd = 2, xlim = xlim_fixed, ylim = ylim,
     xlab = "year the date is fixed to", ylab = "log likelihood of the sample",
     main = "Fixing the date to one year")
points(years, sweep_mean[match(years, grid)], pch = 19, col = "#780000")
text(years, sweep_mean[match(years, grid)], labels = years, pos = 3, cex = 0.8)
abline(h = mean(ll_target_free), lty = 2)
legend("bottomleft", bty = "n", cex = 0.85,
       legend = c("date fixed to that year", "date not fixed (baseline)"),
       lwd = c(2, 1), lty = c(1, 2))
dev.off()

cat("Done. Output in", out_dir, "\n")
