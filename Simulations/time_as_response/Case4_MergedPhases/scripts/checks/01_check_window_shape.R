# Why does the full distribution stay below nominal coverage in Case 4, time as
# response? The full-distribution fit treats every year of a find's window as
# equally likely. In the simulator, the chance that a find gets a given window
# depends on which fine phase its true date falls in: near the first and last
# phases the runs of merged phases are clipped, so windows touching the edges of
# the phase scheme are more likely for some phases than for others.
#
# For each find this computes the exact chance of its window given the fine
# phase of each candidate year, and refits the core datasets with N = 200 and a
# real slope with it. Compared with the flat fits already in recovery_results.csv:
#   exact recovers nominal coverage  -> the window shape (edge clipping) is the cause
#   exact stays below nominal        -> long windows carry too little on sigma
# Writes output/check_window_shape.csv and prints the summary.
library(here)
library(cmdstanr)
library(parallel)

source(here("Simulations", "shared", "dating", "case4_merged_phases.R"))

n.workers <- 40
case.dir <- here("Simulations", "time_as_response", "Case4_MergedPhases")

design <- read.csv(file.path(case.dir, "data", "design.csv"))
design <- design[design$sweep == "core" & design$N == 200 & design$slope_condition == "random", ]
model <- cmdstan_model(here("Simulations", "shared", "models", "linear_dates_response.stan"))

# Chance that a find in fine phase k gets the window made of phases a to b - 1,
# following the merging rule in simulate_merged(): a run of m phases (m from 1
# to merge_max, equally likely) placed at random around the find's phase and
# clipped at the first and last phase
window_probability <- function(a, b, K, merge_max)
{
	prob <- rep(0, K)
	for (k in 1:K)
	{
		for (m in 1:merge_max)
		{
			for (back in 0:(m - 1))
			{
				lo <- max(k - back, 1)
				hi <- min(lo + m, K + 1)
				lo <- max(hi - m, 1)
				if (lo == a && hi == b) prob[k] <- prob[k] + 1 / (merge_max * m)
			}
		}
	}
	prob
}

# Check against the simulator: 200000 finds on one phase scheme, compare how
# often each phase gets each window with window_probability(), allowing for
# sampling noise (4 standard errors; short phases hold few finds)
set.seed(1)
check <- simulate_merged(200000, 0, 0, 1, K = 7, alpha_conc = 1, merge_max = 3,
                         true_date = runif(200000, 100, 900))
bounds <- attr(check, "bounds")
phase <- findInterval(check$True_date, bounds, rightmost.closed = TRUE, all.inside = TRUE)
a <- match(check$Start_date, bounds)
b <- match(check$End_date, bounds)
finds.per.phase <- tabulate(phase, 7)
worst <- 0
for (w in unique(paste(a, b)))
{
	ab <- as.integer(strsplit(w, " ")[[1]])
	observed <- tapply(a == ab[1] & b == ab[2], factor(phase, levels = 1:7), mean)
	expected <- window_probability(ab[1], ab[2], 7, 3)
	se <- sqrt(expected * (1 - expected) / finds.per.phase) + 1e-6
	worst <- max(worst, abs(observed - expected) / se)
}
stopifnot(worst < 4)
cat("window probabilities match the simulator (largest difference", round(worst, 1), "standard errors)\n")

results <- mclapply(1:nrow(design), function(j)
{
	d <- design[j, ]
	# Same data as in 02_recovery_study.R
	set.seed(d$seed)
	x <- runif(d$N, d$x_min, d$x_max)
	true.date <- rnorm(d$N, d$intercept + d$slope * x, d$sigma)
	sim <- simulate_merged(d$N, 0, 0, 1, d$K, d$alpha_conc, merge_max = d$merge_max, true_date = true.date)
	bounds <- attr(sim, "bounds")

	# Candidate years on the grid used in 02_recovery_study.R, weighted by the
	# chance of the find's window given the fine phase of each year
	grid <- seq(floor(min(sim$Start_date) / d$grid_step) * d$grid_step,
	            ceiling(max(sim$End_date) / d$grid_step) * d$grid_step, by = d$grid_step)
	grid.phase <- findInterval(grid, bounds, rightmost.closed = TRUE, all.inside = TRUE)
	year <- matrix(grid, nrow = d$N, ncol = length(grid), byrow = TRUE)
	log_p <- matrix(-Inf, nrow = d$N, ncol = length(grid))
	first <- rep(NA, d$N)
	last <- rep(NA, d$N)
	for (i in 1:d$N)
	{
		inside <- grid >= sim$Start_date[i] & grid <= sim$End_date[i]
		if (!any(inside)) inside[which.min(abs(grid - (sim$Start_date[i] + sim$End_date[i]) / 2))] <- TRUE
		w <- window_probability(match(sim$Start_date[i], bounds), match(sim$End_date[i], bounds),
		                        d$K, d$merge_max)[grid.phase] * inside
		if (sum(w) == 0) w <- as.numeric(inside)
		log_p[i, w > 0] <- log(w[w > 0] / sum(w))
		first[i] <- min(which(w > 0))
		last[i] <- max(which(w > 0))
	}

	dat <- list(n = d$N, n_years = ncol(year), x = x, year = year,
	            log_p = log_p, first = first, last = last,
	            centre = (d$period_start + d$period_end) / 2, x_centre = (d$x_min + d$x_max) / 2,
	            half_cell = d$grid_step / 2)
	inits <- lapply(1:4, function(ch) list(a = runif(1, -0.5, 0.5), b = runif(1, -0.001, 0.001), s = runif(1, 0.5, 1.5)))
	fit <- model$sample(data = dat, chains = 4, parallel_chains = 1, init = inits,
	                    iter_warmup = 500, iter_sampling = 500, adapt_delta = 0.95,
	                    seed = d$seed, refresh = 0, show_messages = FALSE, show_exceptions = FALSE)
	posterior <- as.data.frame(fit$draws(c("slope", "sigma"), format = "draws_df"))
	q <- quantile(posterior$slope, c(0.05, 0.5, 0.95))
	s <- quantile(posterior$sigma, c(0.05, 0.5, 0.95))
	data.frame(dataset_id = d$dataset_id, model = "exact", slope = d$slope,
	           slope_med = q[2], slope_cov90 = d$slope >= q[1] & d$slope <= q[3],
	           sigma_err = s[2] - d$sigma, sigma_cov90 = d$sigma >= s[1] & d$sigma <= s[3],
	           max_rhat = max(fit$summary(c("slope", "sigma"))$rhat))
}, mc.cores = n.workers)

failed <- !sapply(results, is.data.frame)
cat(sum(failed), "of", length(results), "fits failed\n")
results <- do.call(rbind, results[!failed])

# The flat fits are the ones already in the recovery study
flat <- read.csv(file.path(case.dir, "output", "recovery_results.csv"))
flat <- flat[flat$model == "marginal" & flat$dataset_id %in% design$dataset_id, ]
flat <- data.frame(dataset_id = flat$dataset_id, model = "flat",
                   slope = design$slope[match(flat$dataset_id, design$dataset_id)],
                   slope_med = flat$slope_med, slope_cov90 = flat$slope_cov90,
                   sigma_err = flat$sigma_err, sigma_cov90 = flat$sigma_cov90, max_rhat = flat$max_rhat)
results <- rbind(flat, results)
write.csv(results, file.path(case.dir, "output", "check_window_shape.csv"), row.names = FALSE)

for (m in c("flat", "exact"))
{
	s <- results[results$model == m, ]
	cat(sprintf("%-6s slope ratio %.3f | coverage %.2f | sigma bias %+.1f yr | sigma coverage %.2f | R-hat > 1.05: %d\n",
	            m, coef(lm(slope_med ~ slope, data = s))[2], mean(s$slope_cov90),
	            mean(s$sigma_err), mean(s$sigma_cov90), sum(s$max_rhat > 1.05)))
}
