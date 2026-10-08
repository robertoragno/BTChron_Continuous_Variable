# Case 3 recovery study, time as response: every dataset in data/design.csv is
# simulated and fitted twice with shared/models/linear_dates_response.stan:
#   midpoint  centre of each phase window
#   marginal  the whole window, every year in it equally likely
# For a flat range the median equals the midpoint, so it is not fitted.
# Writes one row per fit.
#
# Run with "test" on the command line to fit 20 datasets spread across the
# design and write output/recovery_results_test.csv instead:
#   Rscript .../03_recovery_study.R test
library(here)
library(cmdstanr)
library(parallel)

source(here("Simulations", "shared", "dating", "case3_overlapping_phases.R"))

n.workers <- 47 # Fits run at the same time
output.file <- here("Simulations", "time_as_response", "Case3_OverlappingPhases", "output", "recovery_results.csv")

design <- read.csv(here("Simulations", "time_as_response", "Case3_OverlappingPhases", "data", "design.csv"))
model <- cmdstan_model(here("Simulations", "shared", "models", "linear_dates_response.stan"))

if (length(commandArgs(trailingOnly = TRUE)) > 0 && commandArgs(trailingOnly = TRUE)[1] == "test")
{
	design <- design[round(seq(1, nrow(design), length.out = 20)), ]
	output.file <- sub("\\.csv$", "_test.csv", output.file)
}

# One job per dataset and method ("marginal" is the full distribution)
jobs <- expand.grid(row = 1:nrow(design), method = c("midpoint", "marginal"),
                    stringsAsFactors = FALSE)

results <- mclapply(1:nrow(jobs), function(j)
{
	d <- design[jobs$row[j], ]
	method <- jobs$method[j]

	# Values first, then the true dates from the trend, then the phases (which
	# cover the range of these dates) and the window of each find. The seed is
	# set once here and not again inside simulate_overlap().
	set.seed(d$seed)
	x <- runif(d$N, d$x_min, d$x_max)
	true.date <- rnorm(d$N, d$intercept + d$slope * x, d$sigma)
	sim <- simulate_overlap(d$N, 0, 0, 1, d$K, d$alpha_conc,
	                        overlap = d$overlap, assign_p = d$assign_p,
	                        true_date = true.date)

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
		# Grid from the earliest to the latest window of this dataset. Every grid
		# year is equally likely inside a window and impossible outside it.
		grid <- seq(floor(min(sim$Start_date) / d$grid_step) * d$grid_step,
		            ceiling(max(sim$End_date) / d$grid_step) * d$grid_step, by = d$grid_step)
		year <- matrix(grid, nrow = d$N, ncol = length(grid), byrow = TRUE)
		log_p <- matrix(-Inf, nrow = d$N, ncol = length(grid))
		first <- rep(NA, d$N) # First and last grid column inside each window
		last <- rep(NA, d$N)
		for (i in 1:d$N)
		{
			inside <- grid >= sim$Start_date[i] & grid <= sim$End_date[i]
			# A phase narrower than the grid step gets the grid year nearest its middle
			if (!any(inside)) inside[which.min(abs(grid - (sim$Start_date[i] + sim$End_date[i]) / 2))] <- TRUE
			log_p[i, inside] <- log(1 / sum(inside))
			first[i] <- min(which(inside))
			last[i] <- max(which(inside))
		}
	}

	dat <- list(n = d$N, n_years = ncol(year), x = x, year = year,
	            log_p = log_p, first = first, last = last,
	            centre = (d$period_start + d$period_end) / 2, x_centre = (d$x_min + d$x_max) / 2,
	            half_cell = ifelse(method == "marginal", d$grid_step / 2, 0))
	# Every chain starts near a flat trend through the centre, with sigma about
	# 100 yr. Stan's default random starts can put a chain so far from the data
	# that it never comes back, or trap it in a minor peak with a tiny sigma.
	inits <- lapply(1:4, function(ch) list(a = runif(1, -0.5, 0.5), b = runif(1, -0.001, 0.001), s = runif(1, 0.5, 1.5)))
	fit <- model$sample(data = dat, chains = 4, parallel_chains = 1, init = inits,
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
if (any(failed)) print(results[failed][1:min(3, sum(failed))])
results <- do.call(rbind, results[!failed])
dir.create(dirname(output.file), showWarnings = FALSE, recursive = TRUE)
write.csv(results, output.file, row.names = FALSE)
