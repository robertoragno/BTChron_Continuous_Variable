# Is the full distribution's low coverage in Case 3 caused by the flat window?
# The full-distribution fit treats every year of a phase window as equally
# likely. The simulator does not: a find in a strip shared with a neighbouring
# phase gets this phase's label only with probability assign_p (or 1 - assign_p),
# and the first and last windows run past the period, where no true dates fall.
#
# For each find this computes the exact distribution of its true date given its
# phase label, and refits the core datasets with N = 500 (overlap 0.5,
# assign_p 0.2, random slope) three ways:
#   flat     the window as in 03_recovery_study.R (taken from recovery_results.csv)
#   clipped  the flat window cut to the study period
#   exact    the exact distribution of the true date given the label
# Writes output/check_window_shape.csv and figures/checks/check_window_shape.png
# (the exact and flat distributions of each label for one dataset).
library(here)
library(cmdstanr)
library(parallel)
library(ggplot2)

source(here("Simulations", "shared", "dating", "case3_overlapping_phases.R"))

grid.step <- 5
n.workers <- 48
case.dir <- here("Simulations", "time_as_predictor", "Case3_OverlappingPhases")

design <- read.csv(file.path(case.dir, "data", "design.csv"))
design <- design[design$sweep == "core" & design$N == 500 & design$slope_condition == "random", ]
model <- cmdstan_model(here("Simulations", "shared", "models", "linear_dates.stan"))

# Probability that a find with true date t gets each phase label, for every
# whole year of the period. Follows the assignment rule in simulate_overlap().
label_probability <- function(bounds, overlap, assign_p, period_start, period_end)
{
	K <- length(bounds) - 1
	t <- period_start:period_end
	phase <- findInterval(t, bounds, rightmost.closed = TRUE, all.inside = TRUE)
	reach <- overlap * diff(bounds) / 2
	win_lo <- round(bounds[-(K + 1)] - reach)
	win_hi <- round(bounds[-1] + reach)
	prob <- matrix(0, nrow = length(t), ncol = K)
	prob[cbind(seq_along(t), phase)] <- 1
	if (overlap > 0)
	{
		in_prev <- phase > 1 & t <= win_hi[pmax(phase - 1, 1)]
		in_next <- phase < K & t >= win_lo[pmin(phase + 1, K)]
		only_prev <- which(in_prev & !in_next)
		only_next <- which(in_next & !in_prev)
		prob[cbind(only_prev, phase[only_prev])] <- 1 - assign_p
		prob[cbind(only_prev, phase[only_prev] - 1)] <- assign_p
		prob[cbind(only_next, phase[only_next])] <- assign_p
		prob[cbind(only_next, phase[only_next] + 1)] <- 1 - assign_p
	}
	list(year = t, prob = prob, win_lo = win_lo, win_hi = win_hi)
}

jobs <- expand.grid(row = 1:nrow(design), method = c("clipped", "exact"), stringsAsFactors = FALSE)

results <- mclapply(1:nrow(jobs), function(j)
{
	d <- design[jobs$row[j], ]
	method <- jobs$method[j]
	sim <- simulate_overlap(d$N, d$intercept, d$slope, d$sigma, d$K, d$alpha_conc,
	                        overlap = d$overlap, assign_p = d$assign_p,
	                        period_start = d$period_start, period_end = d$period_end,
	                        seed = d$seed)
	lp <- label_probability(attr(sim, "bounds"), d$overlap, d$assign_p, d$period_start, d$period_end)
	label <- match(paste(sim$Start_date, sim$End_date), paste(lp$win_lo, lp$win_hi))

	grid <- seq(d$time_ref_min, d$time_ref_min + d$time_ref_range, by = grid.step)
	cell <- round((lp$year - grid[1]) / grid.step) + 1
	year <- matrix(grid, nrow = d$N, ncol = length(grid), byrow = TRUE)
	log_p <- matrix(-Inf, nrow = d$N, ncol = length(grid))
	first <- rep(NA, d$N)
	last <- rep(NA, d$N)
	for (i in 1:d$N)
	{
		if (method == "clipped")
		{
			lo <- max(sim$Start_date[i], d$period_start)
			hi <- min(sim$End_date[i], d$period_end)
			w <- as.numeric(grid >= lo & grid <= hi)
			# A clipped window narrower than the grid step gets the grid year nearest its middle
			if (sum(w) == 0) w[which.min(abs(grid - (lo + hi) / 2))] <- 1
		}
		if (method == "exact")
		{
			w <- rep(0, length(grid))
			s <- rowsum(lp$prob[, label[i]], cell)
			w[as.integer(rownames(s))] <- s
		}
		log_p[i, w > 0] <- log(w[w > 0] / sum(w))
		first[i] <- min(which(w > 0))
		last[i] <- max(which(w > 0))
	}

	dat <- list(n = d$N, n_years = ncol(year), y = sim$Value, year = year,
	            log_p = log_p, first = first, last = last, centre = (d$period_start + d$period_end) / 2)
	fit <- model$sample(data = dat, chains = 4, parallel_chains = 1,
	                    iter_warmup = 500, iter_sampling = 500, adapt_delta = 0.95,
	                    seed = d$seed, refresh = 0, show_messages = FALSE,
	                    show_exceptions = FALSE)
	posterior <- as.data.frame(fit$draws(c("slope", "sigma"), format = "draws_df"))
	q <- quantile(posterior$slope, c(0.05, 0.5, 0.95))
	data.frame(dataset_id = d$dataset_id, model = method, slope = d$slope,
	           slope_med = q[2], slope_cov90 = d$slope >= q[1] & d$slope <= q[3],
	           sigma_err = median(posterior$sigma) - d$sigma,
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
                   sigma_err = flat$sigma_err, max_rhat = flat$max_rhat)
results <- rbind(flat, results)
write.csv(results, file.path(case.dir, "output", "check_window_shape.csv"), row.names = FALSE)

for (m in c("flat", "clipped", "exact"))
{
	s <- results[results$model == m, ]
	cat(sprintf("%-8s calibration slope %.3f | 90%% coverage %.2f | bias in sigma %+.3f | R-hat > 1.01: %d\n",
	            m, coef(lm(slope_med ~ slope, data = s))[2], mean(s$slope_cov90),
	            mean(s$sigma_err), sum(s$max_rhat > 1.01)))
}

# Figure: exact and flat distribution of each label, first dataset
d <- design[1, ]
sim <- simulate_overlap(d$N, d$intercept, d$slope, d$sigma, d$K, d$alpha_conc,
                        overlap = d$overlap, assign_p = d$assign_p,
                        period_start = d$period_start, period_end = d$period_end, seed = d$seed)
lp <- label_probability(attr(sim, "bounds"), d$overlap, d$assign_p, d$period_start, d$period_end)
shape <- data.frame()
for (k in 1:ncol(lp$prob))
{
	years <- lp$win_lo[k]:lp$win_hi[k]
	shape <- rbind(shape,
		data.frame(phase = paste("Phase", k), year = lp$year, density = lp$prob[, k] / sum(lp$prob[, k]),
		           distribution = "Exact, given the label"),
		data.frame(phase = paste("Phase", k), year = years, density = 1 / length(years),
		           distribution = "Flat window, as fitted"))
}
shape$phase <- factor(shape$phase, levels = paste("Phase", 1:ncol(lp$prob)))
p <- ggplot(shape, aes(x = year, y = density, colour = distribution)) +
	geom_rect(data = data.frame(xmin = c(-Inf, d$period_end), xmax = c(d$period_start, Inf)),
	          aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf, fill = "Outside the study period"),
	          inherit.aes = FALSE) +
	geom_step() +
	scale_colour_manual(values = c("Exact, given the label" = "#780000", "Flat window, as fitted" = "grey50")) +
	scale_fill_manual(values = c("Outside the study period" = "grey92")) +
	facet_wrap(~phase, scales = "free_y", ncol = 1) +
	labs(title = "Distribution of the true date of a find, given its phase label",
	     subtitle = sprintf("Case 3, one dataset: overlap %.2f, assign_p %.2f", d$overlap, d$assign_p),
	     x = "Year", y = "Probability per year", colour = NULL, fill = NULL) +
	theme_classic() +
	theme(legend.position = "top")
ggsave(file.path(case.dir, "figures", "checks", "check_window_shape.png"), p,
       width = 8, height = 1.2 + 1.3 * ncol(lp$prob), dpi = 300)
