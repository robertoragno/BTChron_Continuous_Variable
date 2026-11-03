# Fit the GINI houses with time as the predictor, log roofed area = a + b * date.
# Two dating treatments, with the model shared with the simulations:
#   midpoint  the centre of each date window, one candidate year
#   marginal  the whole window, every grid year in it equally likely
# Windows reach 915 years, so the grid is in 5-year steps as in the simulations.
#
# Produces: fits.rds, fit_summary.csv, diagnostics.txt
library(here)
library(cmdstanr)
library(posterior)

data.file <- here("Real_Data", "dataset_1", "data", "gini_gb_filtered.csv")
output.dir <- here("Real_Data", "dataset_1", "output")
model <- cmdstan_model(here("Simulations", "shared", "models", "linear_dates.stan"))

GRID.STEP <- 5

houses <- read.csv(data.file)

stan.data <- function(d, method)
{
	if (method == "midpoint")
	{
		year <- matrix((d$Start_date + d$End_date) / 2, ncol = 1)
		log_p <- matrix(0, nrow = nrow(d), ncol = 1)
		first <- rep(1, nrow(d))
		last <- rep(1, nrow(d))
	}
	if (method == "marginal")
	{
		grid <- seq(floor(min(d$Start_date) / GRID.STEP) * GRID.STEP,
		            ceiling(max(d$End_date) / GRID.STEP) * GRID.STEP, by = GRID.STEP)
		year <- matrix(grid, nrow = nrow(d), ncol = length(grid), byrow = TRUE)
		log_p <- matrix(-Inf, nrow = nrow(d), ncol = length(grid))
		first <- rep(NA, nrow(d))
		last <- rep(NA, nrow(d))
		for (i in 1:nrow(d))
		{
			inside <- which(grid >= d$Start_date[i] & grid <= d$End_date[i])
			# One house has a 3-year window that falls between two grid years;
			# it keeps the nearest one rather than being dropped
			if (length(inside) == 0)
				inside <- which.min(abs(grid - (d$Start_date[i] + d$End_date[i]) / 2))
			log_p[i, inside] <- log(1 / length(inside))
			first[i] <- min(inside)
			last[i] <- max(inside)
		}
		stopifnot(all(is.finite(first)))
	}
	list(n = nrow(d), n_years = ncol(year), y = d$Value, year = year,
	     log_p = log_p, first = first, last = last,
	     centre = (min(d$Start_date) + max(d$End_date)) / 2)
}

fits <- list()
summaries <- list()
for (method in c("midpoint", "marginal"))
{
	cat("fitting", method, "(", nrow(houses), "houses )\n")
	fit <- model$sample(data = stan.data(houses, method), chains = 4, parallel_chains = 4,
	                    iter_warmup = 1000, iter_sampling = 1000, adapt_delta = 0.95,
	                    seed = 42, refresh = 200)
	fits[[method]] <- fit

	draws <- as.data.frame(fit$draws(c("intercept", "slope", "sigma"), format = "draws_df"))
	s <- data.frame(method = method, n = nrow(houses))
	for (p in c("intercept", "slope", "sigma"))
	{
		q <- quantile(draws[[p]], c(0.05, 0.5, 0.95))
		s[[paste0(p, "_med")]] <- q[2]
		s[[paste0(p, "_lo")]] <- q[1]
		s[[paste0(p, "_hi")]] <- q[3]
	}
	summaries[[method]] <- s
}

fit.summary <- do.call(rbind, summaries)
write.csv(fit.summary, file.path(output.dir, "fit_summary.csv"), row.names = FALSE)
saveRDS(fits, file.path(output.dir, "fits.rds"))

sink(file.path(output.dir, "diagnostics.txt"))
for (method in names(fits))
{
	cat("---", method, "---\n")
	print(fits[[method]]$diagnostic_summary())
	print(fits[[method]]$summary(c("intercept", "slope", "sigma"))[, c("variable", "rhat", "ess_bulk")])
}
sink()

print(fit.summary, row.names = FALSE)
