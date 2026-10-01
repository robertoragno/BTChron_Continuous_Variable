# Case 1: one simulated radiocarbon dataset, fitted three ways.
# A value changes linearly with calendar date, but each date is only known as a
# calibrated radiocarbon date. The trend is fitted with shared/models/linear_dates.stan
# (the model used in the recovery study) on:
#   1. the midpoint of each 95% calibrated range
#   2. the median calibrated date
#   3. the whole calibrated distribution
# Run for the Hallstatt plateau and for a steep section of the curve.
library(rcarbon)
library(cmdstanr)
library(ggplot2)
library(here)

true.intercept <- 8 # Value at year 0
true.slope <- 0.02 # Change in value per year
true.sigma <- 2 # Noise around the trend
n <- 50 # Number of dated finds
lab.error <- 25 # Lab error of each 14C measurement

windows <- list(c(-800, -400), c(-1600, -1200)) # Calendar windows (BCE is negative)
window.names <- c("the Hallstatt plateau", "the steep section")
file.names <- c("single_fit_comparison.png", "single_fit_comparison_steep.png")

model <- cmdstan_model(here("Simulations", "shared", "models", "linear_dates.stan"))

for (k in 1:2)
{
	window <- windows[[k]]
	centre <- mean(window)

	# Simulate the data
	set.seed(1)
	true.date <- round(runif(n, window[1], window[2]))
	# Radiocarbon age at the true date, plus lab and calibration-curve error
	curve <- uncalibrate(1950 - true.date, verbose = FALSE)
	cra <- round(rnorm(n, curve$ccCRA, sqrt(lab.error^2 + curve$ccError^2)))
	value <- round(true.intercept + true.slope * true.date + rnorm(n, 0, true.sigma), 1)

	# Calibrate, with 600 years either side of the window so no date is cut off
	cal <- calibrate(cra, rep(lab.error, n), calMatrix = TRUE,
	                 timeRange = 1950 - c(window[1] - 600, window[2] + 600), verbose = FALSE)
	years <- 1950 - as.numeric(rownames(cal$calmatrix)) # cal BP to BCE/CE
	prob <- t(cal$calmatrix) # One row per find, one column per year
	prob <- prob / rowSums(prob)

	# Median calibrated date, and midpoint of the 95% HPD range (the most
	# probable years until they hold 95% of the probability)
	median.date <- rep(NA, n)
	midpoint.date <- rep(NA, n)
	for (i in 1:n)
	{
		p <- prob[i, ]
		median.date[i] <- years[which(cumsum(p) >= 0.5)[1]]
		cut <- sort(p, decreasing = TRUE)[which(cumsum(sort(p, decreasing = TRUE)) >= 0.95)[1]]
		midpoint.date[i] <- (min(years[p >= cut]) + max(years[p >= cut])) / 2
	}

	# The same model fitted three times. A point date is one candidate year with
	# probability 1; the full distribution is every year with its calibrated
	# probability, trimmed to the run of years holding 99.99% of it.
	first <- rep(NA, n)
	last <- rep(NA, n)
	for (i in 1:n)
	{
		cum <- cumsum(prob[i, ])
		first[i] <- which(cum >= 0.00005)[1]
		last[i] <- which(cum >= 0.99995)[1]
	}
	inputs <- list(
		midpoint = list(year = matrix(midpoint.date, ncol = 1), log_p = matrix(0, n, 1),
		                first = rep(1, n), last = rep(1, n)),
		median = list(year = matrix(median.date, ncol = 1), log_p = matrix(0, n, 1),
		              first = rep(1, n), last = rep(1, n)),
		full = list(year = matrix(years, nrow = n, ncol = length(years), byrow = TRUE),
		            log_p = log(prob), first = first, last = last))

	# Fitted trend and 95% interval for each, from the posterior draws
	model.names <- c(midpoint = "A  Midpoint of 95% range", median = "B  Median calibrated date",
	                 full = "C  Full distribution")
	pred.years <- seq(window[1], window[2], by = 5)
	trend <- data.frame()
	for (m in names(inputs))
	{
		dat <- c(list(n = n, n_years = ncol(inputs[[m]]$year), y = value, centre = centre), inputs[[m]])
		fit <- model$sample(data = dat, chains = 4, parallel_chains = 4, seed = 1, refresh = 0)
		posterior <- as.data.frame(fit$draws(c("intercept", "slope", "sigma"), format = "draws_df"))
		cat(m, ": slope", round(median(posterior$slope), 4), "(true", true.slope, ")\n")

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
	points <- data.frame(model = rep(model.names[1:2], each = n),
	                     x = c(midpoint.date, median.date), y = c(value, value))

	# Calibrated distributions drawn at the height of each find's value
	shapes <- data.frame()
	for (i in 1:n)
	{
		inside <- prob[i, ] > 0.0001
		shapes <- rbind(shapes, data.frame(id = i, x = years[inside], y = value[i],
		                                   top = value[i] + prob[i, inside] * 60))
	}
	shapes$model <- model.names[3]

	p <- ggplot() +
		geom_ribbon(data = shapes, aes(x = x, ymin = y, ymax = top, group = id),
		            fill = "lightblue", alpha = 0.5) +
		geom_point(data = points, aes(x = x, y = y), colour = "grey40") +
		geom_ribbon(data = trend, aes(x = x, ymin = lo, ymax = hi), fill = "grey50", alpha = 0.4) +
		geom_line(data = trend, aes(x = x, y = fit)) +
		geom_abline(intercept = true.intercept, slope = true.slope, linetype = 2) +
		facet_wrap(~model) +
		coord_cartesian(xlim = window) +
		labs(x = "Calendar year", y = "Value",
		     title = paste("Case 1: one dataset on", window.names[k]),
		     subtitle = "Solid line and band: fitted trend with 95% interval. Dashed line: true trend.") +
		theme_classic()
	ggsave(here("Simulations", "Sim_Case1_Radiocarbon", "figures", file.names[k]),
	       p, width = 12, height = 4.5, dpi = 300)
}
