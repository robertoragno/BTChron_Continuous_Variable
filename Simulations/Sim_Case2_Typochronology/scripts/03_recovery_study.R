# Case 2 recovery study: every dataset in data/design.csv is simulated and
# fitted twice with shared/models/linear_dates.stan, once on the midpoint of
# each window and once on the full window. Writes one row per fit.
# For a flat window the median equals the midpoint, so it is not fitted.
library(here)
library(cmdstanr)
library(parallel)

source(here("Simulations", "Sim_Case2_Typochronology", "scripts", "simulate.R"))

grid.step <- 5 # Years between candidate years for the full-distribution fit
n.workers <- 24 # Fits run at the same time
output.file <- here("Simulations", "Sim_Case2_Typochronology", "output", "recovery_results.csv")

design <- read.csv(here("Simulations", "Sim_Case2_Typochronology", "data", "design.csv"))
model <- cmdstan_model(here("Simulations", "shared", "models", "linear_dates.stan"))

# One job per dataset and method ("marginal" is the full distribution)
jobs <- expand.grid(row = 1:nrow(design), method = c("midpoint", "marginal"),
                    stringsAsFactors = FALSE)

results <- mclapply(1:nrow(jobs), function(j)
{
	d <- design[jobs$row[j], ]
	method <- jobs$method[j]
	sim <- simulate_typo(d$N, d$intercept, d$slope, d$sigma,
	                     prop_coarse_samples = d$prop_coarse_samples,
	                     coarse_frac = d$coarse_frac, fine_frac = d$fine_frac,
	                     precision_trend = d$precision_trend, growth_ratio = d$growth_ratio,
	                     period_start = d$period_start, period_end = d$period_end,
	                     seed = d$seed)

	if (method == "midpoint")
	{
		# One candidate year per find, with probability 1
		year <- matrix((sim$Start_date + sim$End_date) / 2, ncol = 1)
		log_p <- matrix(0, nrow = d$N, ncol = 1)
		first <- rep(1, d$N)
		last <- rep(1, d$N)
	}
	if (method == "marginal")
	{
		# Every grid year, equally likely inside the window and impossible outside
		grid <- seq(d$time_ref_min, d$time_ref_min + d$time_ref_range, by = grid.step)
		year <- matrix(grid, nrow = d$N, ncol = length(grid), byrow = TRUE)
		log_p <- matrix(-Inf, nrow = d$N, ncol = length(grid))
		first <- rep(NA, d$N) # First and last grid column inside each window
		last <- rep(NA, d$N)
		for (i in 1:d$N)
		{
			inside <- grid >= sim$Start_date[i] & grid <= sim$End_date[i]
			# A window narrower than the grid step gets the grid year nearest its middle
			if (!any(inside)) inside[which.min(abs(grid - (sim$Start_date[i] + sim$End_date[i]) / 2))] <- TRUE
			log_p[i, inside] <- log(1 / sum(inside))
			first[i] <- min(which(inside))
			last[i] <- max(which(inside))
		}
	}

	dat <- list(n = d$N, n_years = ncol(year), y = sim$Value, year = year,
	            log_p = log_p, first = first, last = last, centre = (d$period_start + d$period_end) / 2)
	fit <- model$sample(data = dat, chains = 4, parallel_chains = 1,
	                    iter_warmup = 500, iter_sampling = 500, adapt_delta = 0.95,
	                    seed = d$seed, refresh = 0, show_messages = FALSE,
	                    show_exceptions = FALSE)
	posterior <- as.data.frame(fit$draws(c("intercept", "slope", "sigma"), format = "draws_df"))

	# Posterior median, its error, and whether the 50% and 90% intervals contain the truth
	# (04_recovery_plots.R adds the design settings back by dataset_id)
	out <- data.frame(dataset_id = d$dataset_id, model = method)
	for (p in c("intercept", "slope", "sigma"))
	{
		q <- quantile(posterior[[p]], c(0.05, 0.25, 0.5, 0.75, 0.95))
		out[[paste0(p, "_med")]] <- q[3]
		out[[paste0(p, "_err")]] <- q[3] - d[[p]]
		out[[paste0(p, "_cov50")]] <- d[[p]] >= q[2] & d[[p]] <= q[4]
		out[[paste0(p, "_cov90")]] <- d[[p]] >= q[1] & d[[p]] <= q[5]
		out[[paste0(p, "_width90")]] <- q[5] - q[1]
	}
	out$n_divergent <- sum(fit$diagnostic_summary(quiet = TRUE)$num_divergent)
	out$max_rhat <- max(fit$summary(c("intercept", "slope", "sigma"))$rhat)
	out
}, mc.cores = n.workers)

# A failed fit comes back as an error message instead of a data.frame
failed <- !sapply(results, is.data.frame)
cat(sum(failed), "of", length(results), "fits failed\n")
results <- do.call(rbind, results[!failed])
write.csv(results, output.file, row.names = FALSE)
