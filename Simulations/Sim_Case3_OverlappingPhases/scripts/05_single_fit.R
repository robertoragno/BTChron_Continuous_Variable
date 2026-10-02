# Case 3: one simulated dataset, fitted twice with shared/models/linear_dates.stan
# (the model used in the recovery study):
#   A. the midpoint of each dating window
#   B. the whole window, every year in it equally likely
library(cmdstanr)
library(ggplot2)
library(here)

source(here("Simulations", "Sim_Case3_OverlappingPhases", "scripts", "simulate.R"))

true.intercept <- 8 # Value at year 0
true.slope <- 0.02 # Change in value per year
true.sigma <- 1 # Noise around the trend
n <- 200 # Number of dated finds
period <- c(100, 900) # Study period
pad <- 249 # Stretched windows run past the period (as in 01_design.R)
grid.step <- 5 # Years between candidate years
title <- "Case 3: one dataset, phases overlapping by half"

sim <- simulate_overlap(n, true.intercept, true.slope, true.sigma, K = 6, alpha_conc = 1,
                        overlap = 0.5, assign_p = 0.5, period_start = period[1],
                        period_end = period[2], seed = 3)

# The same model fitted twice. The midpoint is one candidate year with
# probability 1; the full window is every grid year inside it, equally likely.
midpoint <- (sim$Start_date + sim$End_date) / 2
grid <- seq(period[1] - pad, period[2] + pad, by = grid.step)
log_p <- matrix(-Inf, nrow = n, ncol = length(grid))
first <- rep(NA, n)
last <- rep(NA, n)
for (i in 1:n)
{
	inside <- grid >= sim$Start_date[i] & grid <= sim$End_date[i]
	if (!any(inside)) inside[which.min(abs(grid - midpoint[i]))] <- TRUE
	log_p[i, inside] <- log(1 / sum(inside))
	first[i] <- min(which(inside))
	last[i] <- max(which(inside))
}
inputs <- list(
	midpoint = list(year = matrix(midpoint, ncol = 1), log_p = matrix(0, n, 1),
	                first = rep(1, n), last = rep(1, n)),
	full = list(year = matrix(grid, nrow = n, ncol = length(grid), byrow = TRUE),
	            log_p = log_p, first = first, last = last))

# Fitted trend and 95% interval for each, from the posterior draws
model <- cmdstan_model(here("Simulations", "shared", "models", "linear_dates.stan"))
model.names <- c(midpoint = "A  Midpoint of each window", full = "B  Full window")
pred.years <- seq(period[1], period[2], by = 5)
trend <- data.frame()
for (m in names(inputs))
{
	dat <- c(list(n = n, n_years = ncol(inputs[[m]]$year), y = sim$Value, centre = mean(period)),
	         inputs[[m]])
	fit <- model$sample(data = dat, chains = 4, parallel_chains = 4, seed = 1, refresh = 0)
	posterior <- as.data.frame(fit$draws(c("intercept", "slope", "sigma"), format = "draws_df"))
	cat(m, ": slope", round(median(posterior$slope), 4), "(true", true.slope, ")",
	    "sigma", round(median(posterior$sigma), 2), "(true", true.sigma, ")\n")

	predmatrix <- matrix(NA, nrow = nrow(posterior), ncol = length(pred.years))
	for (i in 1:nrow(posterior))
	{
		predmatrix[i, ] <- posterior$intercept[i] + posterior$slope[i] * pred.years
	}
	trend <- rbind(trend, data.frame(model = model.names[m], x = pred.years,
	                                 fit = apply(predmatrix, 2, median),
	                                 lo = apply(predmatrix, 2, quantile, 0.025),
	                                 hi = apply(predmatrix, 2, quantile, 0.975)))
}

# What each fit sees: the midpoints in A, the windows in B
points <- data.frame(model = model.names["midpoint"], x = midpoint, y = sim$Value)
windows <- data.frame(model = model.names["full"], xmin = sim$Start_date,
                      xmax = sim$End_date, y = sim$Value)

p <- ggplot() +
	geom_segment(data = windows, aes(x = xmin, xend = xmax, y = y, yend = y), colour = "grey70") +
	geom_point(data = points, aes(x = x, y = y), colour = "grey40") +
	geom_ribbon(data = trend, aes(x = x, ymin = lo, ymax = hi), fill = "grey50", alpha = 0.4) +
	geom_line(data = trend, aes(x = x, y = fit)) +
	geom_abline(intercept = true.intercept, slope = true.slope, linetype = 2) +
	facet_wrap(~model) +
	labs(x = "Calendar year", y = "Value", title = title,
	     subtitle = "Solid line and band: fitted trend with 95% interval. Dashed line: true trend.") +
	theme_classic()
ggsave(here("Simulations", "Sim_Case3_OverlappingPhases", "figures", "single_fit_comparison.png"),
       p, width = 9, height = 4.5, dpi = 300)
