# Fit the pipe stem contexts with time as the response, date = a + b * mean bore.
# Two dating treatments, on two subsets:
#   midpoint  the centre of each date window, one candidate year
#   marginal  the whole window, every year in it equally likely
#   all       the 60 contexts from 01_prepare.R
#   pre1780   the 58 whose window ends by 1780, Binford's own limit on his formula
# 60 contexts, so the fits take seconds and need no tmux session.
#
# Produces: fits.rds, fit_summary.csv, diagnostics.txt
library(here)
library(cmdstanr)
library(posterior)

data.file <- here("Real_Data", "dataset_2", "data", "pipe_stems_subset.csv")
output.dir <- here("Real_Data", "dataset_2", "output")
model <- cmdstan_model(here("Simulations", "shared", "models", "linear_dates_response.stan"))

GRID.STEP <- 1 # One candidate year per year of the window

contexts <- read.csv(data.file)

# Build the Stan data for one subset and one dating treatment ------------------

stan.data <- function(d, method)
{
	if (method == "midpoint")
	{
		year <- matrix((d$date_start + d$date_end) / 2, ncol = 1)
		log_p <- matrix(0, nrow = nrow(d), ncol = 1)
		first <- rep(1, nrow(d))
		last <- rep(1, nrow(d))
	}
	if (method == "marginal")
	{
		grid <- seq(min(d$date_start), max(d$date_end), by = GRID.STEP)
		year <- matrix(grid, nrow = nrow(d), ncol = length(grid), byrow = TRUE)
		log_p <- matrix(-Inf, nrow = nrow(d), ncol = length(grid))
		first <- rep(NA, nrow(d))
		last <- rep(NA, nrow(d))
		for (i in 1:nrow(d))
		{
			inside <- grid >= d$date_start[i] & grid <= d$date_end[i]
			log_p[i, inside] <- log(1 / sum(inside))
			first[i] <- min(which(inside))
			last[i] <- max(which(inside))
		}
		# The narrowest window is 2 years, so every context keeps several grid years
		stopifnot(all(is.finite(first)), all(last - first >= 1))
	}
	list(n = nrow(d), n_years = ncol(year), x = d$mean_bore, year = year,
	     log_p = log_p, first = first, last = last,
	     centre = (min(d$date_start) + max(d$date_end)) / 2,
	     x_centre = (min(d$mean_bore) + max(d$mean_bore)) / 2,
	     half_cell = ifelse(method == "marginal", GRID.STEP / 2, 0))
}

# Fit -------------------------------------------------------------------------

subsets <- list(all = contexts, pre1780 = contexts[!contexts$past_1780, ])
jobs <- expand.grid(subset = names(subsets), method = c("midpoint", "marginal"),
                    stringsAsFactors = FALSE)

fits <- list()
summaries <- list()
for (j in 1:nrow(jobs))
{
	d <- subsets[[jobs$subset[j]]]
	label <- paste(jobs$subset[j], jobs$method[j], sep = "_")
	cat("fitting", label, "(", nrow(d), "contexts )\n")

	# Chains start near a flat trend through the centre with sigma about 100 yr,
	# as in the simulations: random starts can strand a chain far from the data
	inits <- lapply(1:4, function(ch)
		list(a = runif(1, -0.5, 0.5), b = runif(1, -0.001, 0.001), s = runif(1, 0.5, 1.5)))
	fit <- model$sample(data = stan.data(d, jobs$method[j]), chains = 4, parallel_chains = 4,
	                    init = inits, iter_warmup = 1000, iter_sampling = 1000,
	                    adapt_delta = 0.95, seed = 42, refresh = 0, show_messages = FALSE)
	fits[[label]] <- fit

	draws <- as.data.frame(fit$draws(c("intercept", "slope", "sigma"), format = "draws_df"))
	s <- data.frame(subset = jobs$subset[j], method = jobs$method[j], n = nrow(d))
	for (p in c("intercept", "slope", "sigma"))
	{
		q <- quantile(draws[[p]], c(0.05, 0.5, 0.95))
		s[[paste0(p, "_med")]] <- q[2]
		s[[paste0(p, "_lo")]] <- q[1]
		s[[paste0(p, "_hi")]] <- q[3]
	}
	summaries[[label]] <- s
}

fit.summary <- do.call(rbind, summaries)
write.csv(fit.summary, file.path(output.dir, "fit_summary.csv"), row.names = FALSE)
saveRDS(fits, file.path(output.dir, "fits.rds"))

# Diagnostics: check these before reading any of the figures -------------------

sink(file.path(output.dir, "diagnostics.txt"))
for (label in names(fits))
{
	cat("---", label, "---\n")
	print(fits[[label]]$diagnostic_summary())
	print(fits[[label]]$summary(c("intercept", "slope", "sigma"))[, c("variable", "rhat", "ess_bulk")])
}
sink()

print(fit.summary, row.names = FALSE)
