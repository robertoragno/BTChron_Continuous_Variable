# Do the datasets that did not converge in Stan converge in Crema's Nimble model?
# 20 plateau datasets whose full-distribution fit had R-hat above 1.05 in the
# recovery study (12 with a real slope, 8 with a zero slope), refitted with the
# Nimble model and settings of 00_compare_nimble.R (Crema's
# measurement_error_example.R). The data are rebuilt exactly as in
# 02_recovery_study.R.
#
# One Nimble fit takes about an hour (more for N = 500), so each dataset runs as
# its own R process, named on the command line:
#   Rscript .../01_nimble_nonconverged.R 1093
# Each run writes output/nimble_nonconverged/<dataset_id>.csv. Run with no
# arguments to choose the datasets (written to picks.csv) and, once the fits are
# done, to collect them into one table.
library(rcarbon)
library(nimbleCarbon)
library(coda)
library(here)

source(here("Simulations", "shared", "dating", "case1_radiocarbon.R"))

case.dir <- here("Simulations", "time_as_response", "Case1_Radiocarbon")
out.dir <- file.path(case.dir, "output", "nimble_nonconverged")
dir.create(out.dir, showWarnings = FALSE, recursive = TRUE)
design <- read.csv(file.path(case.dir, "data", "design.csv"))
stan <- read.csv(file.path(case.dir, "output", "recovery_results.csv"))
stan <- stan[stan$model == "marginal", ]
args <- commandArgs(trailingOnly = TRUE)

if (length(args) == 0)
{
	if (!file.exists(file.path(out.dir, "picks.csv")))
	{
		# Choose the datasets: 4 real-slope datasets at each of N = 50, 100, 200,
		# 2 zero-slope datasets at each of N = 50, 100, 200, 500
		bad <- merge(stan[stan$max_rhat >= 1.05, c("dataset_id", "max_rhat")], design)
		set.seed(2026)
		picks <- c()
		for (n in c(50, 100, 200)) picks <- c(picks, sample(bad$dataset_id[bad$N == n & bad$slope_condition == "random"], 4))
		for (n in c(50, 100, 200, 500)) picks <- c(picks, sample(bad$dataset_id[bad$N == n & bad$slope_condition == "zero"], 2))
		write.csv(data.frame(dataset_id = picks), file.path(out.dir, "picks.csv"), row.names = FALSE)
		cat("picked", length(picks), "datasets\n")
	}
	# Collect the finished fits
	files <- list.files(out.dir, pattern = "^[0-9]+\\.csv$", full.names = TRUE)
	if (length(files) > 0)
	{
		results <- do.call(rbind, lapply(files, read.csv))
		results <- merge(results, stan[, c("dataset_id", "slope_med", "max_rhat")])
		names(results)[names(results) == "slope_med"] <- "stan_slope"
		names(results)[names(results) == "max_rhat"] <- "stan_rhat"
		print(results[order(results$slope_condition, results$N), ], digits = 3)
		write.csv(results, file.path(case.dir, "output", "nimble_nonconverged.csv"), row.names = FALSE)
	}
	quit()
}

d <- design[design$dataset_id == as.integer(args[1]), ]

# The same data as in 02_recovery_study.R
set.seed(d$seed)
x <- runif(d$N, d$x_min, d$x_max)
true.date <- round(rnorm(d$N, d$intercept + d$slope * x, d$sigma))
sim <- simulate_c14(d$N, 0, 0, 1, c(d$window_start, d$window_end), d$lab_error, true_date = true.date)

# Crema's Nimble model: dates in years BP, so its beta has the opposite sign to
# our slope (calendar years per unit of x)
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
median.calibrated <- calibrate(sim$CRA, sim$Error, verbose = FALSE) |> medCal()
dat <- list(cra = sim$CRA, cra.error = sim$Error)
constants <- list(d = x, calBP = intcal20$CalBP, C14BP = intcal20$C14Age,
                  C14err = intcal20$C14Age.sigma, n = d$N)
inits <- list(theta = median.calibrated, alpha = 3000, beta = 1/2, sigma = 50)
seconds <- system.time(
	out <- nimbleMCMC(code = model, constants = constants, inits = inits, data = dat, nchains = 4,
	                  niter = 100000, nburnin = 50000, samplesAsCodaMCMC = TRUE,
	                  monitors = c("alpha", "beta", "sigma"), setSeed = d$seed * 10 + 1:4)
)[3]

# Slope and sigma, R-hat, and each chain's median slope
rhat <- gelman.diag(out[, c("alpha", "beta", "sigma")], autoburnin = FALSE, multivariate = FALSE)$psrf[, 1]
posterior <- do.call(rbind.data.frame, out)
slope <- -posterior$beta
q <- quantile(slope, c(0.05, 0.5, 0.95))
chain.slopes <- sapply(out, function(ch) median(-ch[, "beta"]))
write.csv(data.frame(dataset_id = d$dataset_id, N = d$N, slope_condition = d$slope_condition,
                     true_slope = d$slope, true_sigma = d$sigma,
                     nimble_slope = q[2], nimble_cov90 = d$slope >= q[1] & d$slope <= q[3],
                     nimble_sigma = median(posterior$sigma), nimble_rhat = max(rhat),
                     chain_slopes = paste(round(chain.slopes, 3), collapse = " "), seconds = seconds),
          file.path(out.dir, paste0(d$dataset_id, ".csv")), row.names = FALSE)
cat("dataset", d$dataset_id, "done in", round(seconds / 60), "min\n")
