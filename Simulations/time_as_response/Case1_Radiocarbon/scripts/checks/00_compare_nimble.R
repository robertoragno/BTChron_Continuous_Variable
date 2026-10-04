# Does linear_dates_response.stan give the same answer as Crema's Nimble model?
# Follows measurement_error_example.R (ercrema/statistical_modelling_review):
# same simulation, same Nimble model and settings. Each dataset is also fitted
# with our Stan model, which sums over each date's candidate years (5-year grid)
# instead of sampling the dates. 10 datasets on the Hallstatt plateau (Crema's
# values) and 10 on the steep section (same trend, moved 785 years earlier).
#
# One Nimble fit takes about 40 minutes, so the datasets are run side by side,
# one R process each, by naming them on the command line:
#   Rscript .../00_compare_nimble.R plateau 3
# Each run writes output/compare_nimble/<window>_<dataset>.csv. Run with no
# arguments to only collect the finished rows into the table and the figure.
library(rcarbon)
library(nimbleCarbon)
library(coda)
library(cmdstanr)
library(ggplot2)
library(here)

true.beta <- 1/3.7 # True slope
true.sigma <- 70 # True error
n <- 100 # Sample size
windows <- c(plateau = 2700, steep = 3485) # True intercept (BP) for each window
n.datasets <- 10 # Datasets per window

# Stan model, compiled once
stan.model <- cmdstan_model(here("Simulations", "shared", "models", "linear_dates_response.stan"))
grid <- seq(-2600, 400, by = 5) # Candidate years (BCE/CE) for the Stan model

# Regression model in nimble
model <- nimbleCode({
	for (i in 1:n)
	{
		mu[i]  <- alpha +  beta * d[i]
		theta[i] ~ dnorm(mean=mu[i],sd=sigma)
		c14age[i] <- interpLin(z=theta[i],x=calBP[],y=C14BP[]);
		sigmaCurve[i] <- interpLin(z=theta[i],x=calBP[],y=C14err[]);
		sigmaDate[i] <- (cra.error[i]^2+sigmaCurve[i]^2)^(1/2);
		cra[i] ~ dnorm(mean=c14age[i],sd=sigmaDate[i])
	}
	alpha ~ dunif(1000,5000)
	beta ~ dnorm(0,1)
	sigma ~ dexp(0.1)
})
data(intcal20)

args <- commandArgs(trailingOnly = TRUE)
row.dir <- here("Simulations", "time_as_response", "Case1_Radiocarbon", "output", "compare_nimble")
dir.create(row.dir, showWarnings = FALSE, recursive = TRUE)
for (w in names(windows))
{
	for (k in 1:n.datasets)
	{
		if (length(args) == 0 || args[1] != w || as.integer(args[2]) != k) next

		# Simulate Regression Data (as in Crema's script)
		set.seed(k)
		true.alpha <- windows[w]
		d <- runif(n=n,min=0,max=1000) |> round()
		calendar.dates <- rnorm(n=n,mean=true.alpha - true.beta*d,sd=true.sigma) |> round()
		c14ages <- uncalibrate(calendar.dates)[,4] #back-calibrate calendar date in 14C age
		c14ages.error <- rep(30,n) #assign errors
		median.calibrated <- calibrate(c14ages,c14ages.error,verbose=FALSE) |> medCal()

		# Regression on median calibrated dates
		fitmed <- lm(median.calibrated~d)

		# Bayesian EIV via Nimble (Crema's settings)
		dat <- list(cra = c14ages, cra.error = c14ages.error)
		constants <- list(d = d, calBP = intcal20$CalBP, C14BP = intcal20$C14Age,
		                  C14err = intcal20$C14Age.sigma, n = n)
		inits <- list(theta = median.calibrated, alpha = 3000, beta = 1/2, sigma = 50)
		time.nimble <- system.time(
			out <- nimbleMCMC(code=model,constants=constants,inits=inits,data=dat,nchains=4,niter=100000,nburnin=50000,samplesAsCodaMCMC=TRUE,monitors=c('theta','alpha','beta','sigma'))
		)[3]
		posterior <- do.call(rbind.data.frame,out)

		# Our Stan model: calibrated probabilities summed into 5-year cells.
		# Stan works in calendar years (BCE negative), so its slope has the
		# opposite sign to Crema's beta, which is in years BP.
		cal <- calibrate(c14ages, c14ages.error, calMatrix = TRUE, verbose = FALSE,
		                 timeRange = 1950 - c(min(grid) - 5, max(grid) + 5))
		years <- 1950 - as.numeric(rownames(cal$calmatrix))
		cell <- round((years - grid[1]) / 5) + 1
		on.grid <- cell >= 1 & cell <= length(grid)
		prob <- t(rowsum(cal$calmatrix[on.grid, ], cell[on.grid]))
		prob <- prob / rowSums(prob)
		first <- rep(NA, n)
		last <- rep(NA, n)
		for (i in 1:n)
		{
			cum <- cumsum(prob[i, ])
			first[i] <- which(cum >= 0.00005)[1]
			last[i] <- which(cum >= 0.99995)[1]
		}
		stan.data <- list(n = n, n_years = length(grid), x = d,
		                  year = matrix(grid, nrow = n, ncol = length(grid), byrow = TRUE),
		                  log_p = log(prob), first = first, last = last,
		                  centre = 1950 - mean(median.calibrated), x_centre = 500, half_cell = 2.5)
		time.stan <- system.time(
			fit <- stan.model$sample(data = stan.data, chains = 4, parallel_chains = 1,
			                         seed = k, refresh = 0, show_messages = FALSE)
		)[3]
		stan.posterior <- as.data.frame(fit$draws(c("intercept", "slope", "sigma"), format = "draws_df"))

		# The dates are not parameters in the Stan model, but their posterior can be
		# drawn after the fit: for every posterior draw, each sample gets one year,
		# picked with weight (calibrated probability x density of the trend there).
		# This is the distribution Nimble samples theta from.
		theta.stan <- matrix(NA, nrow = nrow(stan.posterior), ncol = n)
		for (i in 1:n)
		{
			cols <- first[i]:last[i]
			for (s in 1:nrow(stan.posterior))
			{
				weight <- prob[i, cols] * dnorm(grid[cols], stan.posterior$intercept[s] + stan.posterior$slope[s] * d[i],
				                                stan.posterior$sigma[s])
				theta.stan[s, i] <- grid[cols][sample.int(length(cols), 1, prob = weight)]
			}
		}
		theta.stan <- 1950 - theta.stan # to years BP, as in Nimble
		theta.nimble <- as.matrix(posterior[, paste0("theta[", 1:n, "]")])

		# Each sample's date: true value, median and 95% interval from each model
		write.csv(data.frame(
			window = w, dataset = k, sample = 1:n, true.date = calendar.dates,
			nimble.median = apply(theta.nimble, 2, median),
			nimble.lo = apply(theta.nimble, 2, quantile, 0.025), nimble.hi = apply(theta.nimble, 2, quantile, 0.975),
			stan.median = apply(theta.stan, 2, median),
			stan.lo = apply(theta.stan, 2, quantile, 0.025), stan.hi = apply(theta.stan, 2, quantile, 0.975)),
			file.path(row.dir, paste0("theta_", w, "_", k, ".csv")), row.names = FALSE)

		# Every draw of the first six samples of dataset 1, to compare the shapes
		if (k == 1)
		{
			write.csv(rbind(
				data.frame(window = w, model = "Nimble", sample = rep(1:6, each = nrow(theta.nimble)), date = c(theta.nimble[, 1:6])),
				data.frame(window = w, model = "Stan", sample = rep(1:6, each = nrow(theta.stan)), date = c(theta.stan[, 1:6]))),
				file.path(row.dir, paste0("draws_", w, ".csv")), row.names = FALSE)
		}

		write.csv(data.frame(
			window = w, dataset = k, true.beta = true.beta,
			median.beta = -coef(fitmed)[2],
			nimble.beta = -median(posterior$beta),
			nimble.lo = -quantile(posterior$beta, 0.975), nimble.hi = -quantile(posterior$beta, 0.025),
			stan.beta = median(stan.posterior$slope),
			stan.lo = quantile(stan.posterior$slope, 0.025), stan.hi = quantile(stan.posterior$slope, 0.975),
			nimble.sigma = median(posterior$sigma), stan.sigma = median(stan.posterior$sigma),
			nimble.seconds = time.nimble, stan.seconds = time.stan),
			file.path(row.dir, paste0(w, "_", k, ".csv")), row.names = FALSE)
		cat(w, k, "done\n")
	}
}

# Collect every finished dataset
results <- do.call(rbind, lapply(list.files(row.dir, pattern = "^(plateau|steep)_.*csv$", full.names = TRUE), read.csv))
if (is.null(results)) quit()
print(results, digits = 3)

# Slope from each model, with 95% intervals, one point per dataset
p <- ggplot(results, aes(x = nimble.beta, y = stan.beta)) +
	geom_abline(intercept = 0, slope = 1, colour = "grey70") +
	geom_hline(yintercept = true.beta, linetype = 2) +
	geom_vline(xintercept = true.beta, linetype = 2) +
	geom_errorbar(aes(ymin = stan.lo, ymax = stan.hi), width = 0, colour = "grey50") +
	geom_errorbar(aes(xmin = nimble.lo, xmax = nimble.hi), width = 0, orientation = "y", colour = "grey50") +
	geom_point() +
	facet_wrap(~window) +
	labs(x = "Slope, Nimble (Crema's model)", y = "Slope, Stan (dates summed over)",
	     subtitle = "Bars: 95% intervals. Dashed: true slope. Grey line: the two models agree.") +
	theme_classic()
dir.create(here("Simulations", "time_as_response", "Case1_Radiocarbon", "figures", "checks"),
           showWarnings = FALSE, recursive = TRUE)
ggsave(here("Simulations", "time_as_response", "Case1_Radiocarbon", "figures", "checks", "compare_nimble.png"),
       p, width = 8, height = 4, dpi = 300)

# Dates: Stan against Nimble, every sample of every dataset (median and 95% interval)
theta <- do.call(rbind, lapply(list.files(row.dir, pattern = "^theta_", full.names = TRUE), read.csv))
p <- ggplot(theta, aes(x = nimble.median, y = stan.median)) +
	geom_abline(intercept = 0, slope = 1, colour = "grey70") +
	geom_errorbar(aes(ymin = stan.lo, ymax = stan.hi), width = 0, colour = "grey80") +
	geom_errorbar(aes(xmin = nimble.lo, xmax = nimble.hi), width = 0, orientation = "y", colour = "grey80") +
	geom_point(size = 0.6) +
	facet_wrap(~window, scales = "free") +
	labs(x = "Date of each sample, Nimble (years BP)", y = "Date of each sample, Stan (years BP)",
	     subtitle = "Posterior median and 95% interval of every sample's date. Grey line: the two models agree.") +
	theme_classic()
ggsave(here("Simulations", "time_as_response", "Case1_Radiocarbon", "figures", "checks", "compare_nimble_dates.png"),
       p, width = 8, height = 4, dpi = 300)

# Shape of the date posterior for six samples of dataset 1
draws <- do.call(rbind, lapply(list.files(row.dir, pattern = "^draws_", full.names = TRUE), read.csv))
# Stan dates sit on the 5-year grid; spread them across their cell so the density is smooth
is.stan <- draws$model == "Stan"
draws$date[is.stan] <- draws$date[is.stan] + runif(sum(is.stan), -2.5, 2.5)
p <- ggplot(draws, aes(x = date, colour = model)) +
	geom_density(adjust = 0.5) +
	scale_x_reverse() +
	facet_wrap(window ~ sample, scales = "free", ncol = 6, labeller = label_both) +
	labs(x = "Date (years BP)", y = "Posterior density", colour = NULL) +
	theme_classic() +
	theme(legend.position = "top")
ggsave(here("Simulations", "time_as_response", "Case1_Radiocarbon", "figures", "checks", "compare_nimble_date_shapes.png"),
       p, width = 12, height = 5, dpi = 300)
