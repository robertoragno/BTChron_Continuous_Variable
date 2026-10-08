# Does the midpoint keep its slope when the overlaps are asymmetric?
# In the main study each phase window is stretched by the same amount on both
# sides, so the thinning of dates at the two edges is balanced and the window
# midpoint stays near the middle of the true dates. Here a share of the same
# total stretch (late_share) goes on the late side of each phase, as when new
# types appear abruptly and old types linger:
#   late_share 0.5   half on each side, as in the main study
#   late_share 0.75  three quarters on the late side
#   late_share 1     all on the late side
# 100 datasets per late_share at overlap 0.5, assign_p 0.5, N = 200, each fitted on the
# midpoint and on the flat window with shared/models/linear_dates.stan.
# The dating step is written out here (not in shared/dating/) so the main
# simulator stays as it is; at late_share 0.5 it gives exactly what simulate_overlap()
# gives, which is checked before fitting.
# Writes output/check_asymmetric_overlap.csv and prints the summary.
library(here)
library(cmdstanr)
library(parallel)

source(here("Simulations", "shared", "dating", "case3_overlapping_phases.R"))
source(here("Simulations", "shared", "scripts", "design.R"))

n.rep <- 100
N <- 200
overlap <- 0.5
assign.p <- 0.5
period <- c(100, 900)
grid.step <- 5
n.workers <- 24
case.dir <- here("Simulations", "time_as_predictor", "Case3_OverlappingPhases")

set.seed(2026)
design <- expand.grid(rep = 1:n.rep, late_share = c(0.5, 0.75, 1), slope_condition = "random",
                      stringsAsFactors = FALSE)
design <- add_true_trend(design)
design$K <- sample(5:10, nrow(design), replace = TRUE)
design$alpha_conc <- 10^(-1 + 2 * rbeta(nrow(design), 2, 3.5))
model <- cmdstan_model(here("Simulations", "shared", "models", "linear_dates.stan"))

# One dataset: the same steps as simulate_overlap(), except that each window is
# stretched by overlap x phase length in total, a share late_share of it on the late side
simulate_overlap_asymmetric <- function(d)
{
	set.seed(d$seed)
	span <- diff(period)
	repeat {
		weights <- rgamma(d$K, shape = d$alpha_conc, rate = 1)
		weights <- weights / sum(weights)
		if (max(weights) <= 0.6) break
	}
	lengths <- 1 + weights * (span - d$K)
	bounds <- round(period[1] + c(0, cumsum(lengths)))
	bounds[c(1, d$K + 1)] <- period
	true.date <- round(runif(N, period[1], period[2]))
	phase <- findInterval(true.date, bounds, rightmost.closed = TRUE, all.inside = TRUE)
	stretch <- overlap * diff(bounds)
	win.lo <- round(bounds[-(d$K + 1)] - stretch * (1 - d$late_share))
	win.hi <- round(bounds[-1] + stretch * d$late_share)
	in.prev <- phase > 1 & true.date <= win.hi[pmax(phase - 1, 1)]
	in.next <- phase < d$K & true.date >= win.lo[pmin(phase + 1, d$K)]
	only.prev <- in.prev & !in.next
	only.next <- in.next & !in.prev
	flip <- runif(N)
	phase[only.prev] <- phase[only.prev] - (flip[only.prev] < assign.p)
	phase[only.next] <- phase[only.next] + (flip[only.next] >= assign.p)
	data.frame(Start_date = win.lo[phase], End_date = win.hi[phase], True_date = true.date,
	           Value = round(d$intercept + d$slope * true.date + rnorm(N, 0, d$sigma), 1))
}

# Check: at late_share 0.5 the written-out steps give exactly what simulate_overlap() gives,
# and at every late_share each window contains its true date
for (r in which(design$late_share == 0.5)[1:10])
{
	d <- design[r, ]
	a <- simulate_overlap_asymmetric(d)
	b <- simulate_overlap(N, d$intercept, d$slope, d$sigma, d$K, d$alpha_conc, overlap = overlap,
	                      assign_p = assign.p, period_start = period[1], period_end = period[2], seed = d$seed)
	attr(b, "bounds") <- NULL
	stopifnot(isTRUE(all.equal(a, b)))
}
for (r in 1:nrow(design))
{
	s <- simulate_overlap_asymmetric(design[r, ])
	stopifnot(all(s$True_date >= s$Start_date & s$True_date <= s$End_date))
}
cat("dating checks ok\n")

jobs <- expand.grid(row = 1:nrow(design), method = c("midpoint", "marginal"), stringsAsFactors = FALSE)
results <- mclapply(1:nrow(jobs), function(j)
{
	d <- design[jobs$row[j], ]
	method <- jobs$method[j]
	sim <- simulate_overlap_asymmetric(d)
	if (method == "midpoint")
	{
		year <- matrix((sim$Start_date + sim$End_date) / 2, ncol = 1)
		log_p <- matrix(0, nrow = N, ncol = 1)
		first <- rep(1, N)
		last <- rep(1, N)
	}
	if (method == "marginal")
	{
		# Flat window, on a grid from the earliest to the latest window of this dataset
		grid <- seq(floor(min(sim$Start_date) / grid.step) * grid.step,
		            ceiling(max(sim$End_date) / grid.step) * grid.step, by = grid.step)
		year <- matrix(grid, nrow = N, ncol = length(grid), byrow = TRUE)
		log_p <- matrix(-Inf, nrow = N, ncol = length(grid))
		first <- rep(NA, N)
		last <- rep(NA, N)
		for (i in 1:N)
		{
			inside <- grid >= sim$Start_date[i] & grid <= sim$End_date[i]
			if (!any(inside)) inside[which.min(abs(grid - (sim$Start_date[i] + sim$End_date[i]) / 2))] <- TRUE
			log_p[i, inside] <- log(1 / sum(inside))
			first[i] <- min(which(inside))
			last[i] <- max(which(inside))
		}
	}
	dat <- list(n = N, n_years = ncol(year), y = sim$Value, year = year,
	            log_p = log_p, first = first, last = last, centre = mean(period))
	fit <- model$sample(data = dat, chains = 4, parallel_chains = 1,
	                    iter_warmup = 500, iter_sampling = 500, adapt_delta = 0.95,
	                    seed = d$seed, refresh = 0, show_messages = FALSE, show_exceptions = FALSE)
	posterior <- as.data.frame(fit$draws(c("slope", "sigma"), format = "draws_df"))
	q <- quantile(posterior$slope, c(0.05, 0.5, 0.95))
	data.frame(dataset_id = d$dataset_id, late_share = d$late_share, model = method, slope = d$slope,
	           slope_med = q[2], slope_cov90 = d$slope >= q[1] & d$slope <= q[3],
	           sigma_err = median(posterior$sigma) - d$sigma,
	           max_rhat = max(fit$summary(c("slope", "sigma"))$rhat))
}, mc.cores = n.workers, mc.preschedule = FALSE)

failed <- !sapply(results, is.data.frame)
cat(sum(failed), "of", length(results), "fits failed\n")
results <- do.call(rbind, results[!failed])
write.csv(results, file.path(case.dir, "output", "check_asymmetric_overlap.csv"), row.names = FALSE)

for (l in c(0.5, 0.75, 1))
{
	for (m in c("midpoint", "marginal"))
	{
		s <- results[results$late_share == l & results$model == m, ]
		cat(sprintf("late share %-4s %-8s calibration slope %.3f | 90%% coverage %.2f | bias in sigma %+.3f | R-hat > 1.01: %d\n",
		            l, m, coef(lm(slope_med ~ slope, data = s))[2], mean(s$slope_cov90),
		            mean(s$sigma_err), sum(s$max_rhat > 1.01)))
	}
}
