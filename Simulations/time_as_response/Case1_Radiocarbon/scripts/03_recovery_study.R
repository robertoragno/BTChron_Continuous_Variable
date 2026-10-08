# Case 1 recovery study, time as response: every dataset in data/design.csv is
# simulated, calibrated, and fitted three times with
# shared/models/linear_dates_response.stan:
#   midpoint  centre of the 95% HPD range of each calibrated date
#   median    median of each calibrated date
#   marginal  the whole calibrated distribution
# Writes one row per fit.
library(here)
library(cmdstanr)
library(parallel)
library(rcarbon)

source(here("Simulations", "shared", "dating", "case1_radiocarbon.R"))

n.workers <- 47 # Fits run at the same time
output.file <- here("Simulations", "time_as_response", "Case1_Radiocarbon", "output", "recovery_results.csv")

design <- read.csv(here("Simulations", "time_as_response", "Case1_Radiocarbon", "data", "design.csv"))
model <- cmdstan_model(here("Simulations", "shared", "models", "linear_dates_response.stan"))

# One job per dataset and method ("marginal" is the full distribution)
jobs <- expand.grid(row = 1:nrow(design), method = c("midpoint", "median", "marginal"),
                    stringsAsFactors = FALSE)

results <- mclapply(1:nrow(jobs), function(j)
{
	d <- design[jobs$row[j], ]
	method <- jobs$method[j]

	# Values first, then the true dates from the trend, then radiocarbon ages
	set.seed(d$seed)
	x <- runif(d$N, d$x_min, d$x_max)
	true.date <- round(rnorm(d$N, d$intercept + d$slope * x, d$sigma))
	sim <- simulate_c14(d$N, 0, 0, 1, c(d$window_start, d$window_end), d$lab_error,
	                    true_date = true.date)

	# Calibrate over the whole grid, then add up the yearly probabilities into
	# 5-year cells: one row per find, one column per grid year
	grid <- seq(d$grid_min, d$grid_max, by = d$grid_step)
	cal <- calibrate(sim$CRA, sim$Error, calMatrix = TRUE, verbose = FALSE,
	                 timeRange = 1950 - c(min(grid) - d$grid_step, max(grid) + d$grid_step))
	years <- 1950 - as.numeric(rownames(cal$calmatrix))
	cell <- round((years - grid[1]) / d$grid_step) + 1
	on.grid <- cell >= 1 & cell <= length(grid)
	prob <- t(rowsum(cal$calmatrix[on.grid, ], cell[on.grid]))
	prob <- prob / rowSums(prob)

	if (method == "marginal")
	{
		# Each find uses the shortest run of grid years holding 99.99% of its probability
		year <- matrix(grid, nrow = d$N, ncol = length(grid), byrow = TRUE)
		log_p <- log(prob)
		first <- rep(NA, d$N)
		last <- rep(NA, d$N)
		for (i in 1:d$N)
		{
			cum <- cumsum(prob[i, ])
			first[i] <- which(cum >= 0.00005)[1]
			last[i] <- which(cum >= 0.99995)[1]
		}
	} else {
		# One candidate year per find, with probability 1
		point.date <- rep(NA, d$N)
		for (i in 1:d$N)
		{
			p <- prob[i, ]
			if (method == "median") point.date[i] <- grid[which(cumsum(p) >= 0.5)[1]]
			if (method == "midpoint")
			{
				# 95% HPD: the most probable years until they hold 95% of the probability
				cut <- sort(p, decreasing = TRUE)[which(cumsum(sort(p, decreasing = TRUE)) >= 0.95)[1]]
				point.date[i] <- (min(grid[p >= cut]) + max(grid[p >= cut])) / 2
			}
		}
		year <- matrix(point.date, ncol = 1)
		log_p <- matrix(0, nrow = d$N, ncol = 1)
		first <- rep(1, d$N)
		last <- rep(1, d$N)
	}

	dat <- list(n = d$N, n_years = ncol(year), x = x, year = year,
	            log_p = log_p, first = first, last = last,
	            centre = (d$window_start + d$window_end) / 2, x_centre = (d$x_min + d$x_max) / 2,
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
results <- do.call(rbind, results[!failed])
dir.create(dirname(output.file), showWarnings = FALSE, recursive = TRUE)
write.csv(results, output.file, row.names = FALSE)
