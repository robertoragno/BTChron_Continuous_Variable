# Reversing Case 1: the calendar date as the response, not the predictor
#
# Case 1 puts the date on the x axis: Value = intercept + slope * date + noise,
# so the calibrated date is a predictor measured with error. That is classical
# errors-in-variables, and it attenuates the slope (0.80-0.83 on the steep
# section of the curve).
#
# This script reverses the experiment, following
# https://github.com/ercrema/statistical_modelling_review/blob/main/measurement_error_example.R
# The calendar date is now the response of a known covariate d:
#
#   calendar date ~ Normal(alpha - beta * d, sigma)
#
# and only the date is dated by radiocarbon. Error in the response is Berkson
# rather than classical, so the prediction is that the slope is recovered with
# little or no attenuation, and that the cost of the dating shows up in sigma
# instead. This is the radiocarbon analogue of what the other cases already
# show for midpoints of a uniform range.
#
# It also asks three further questions:
#   - do the radiocarbon dates change shape once the trend is modelled, the way
#     the typochronology dates did
#   - does the agreement index (nimbleCarbon) say the new shape still agrees
#     with the plain calibrated date
#   - what happens to the likelihood when one date is fixed to a constant year
#     instead of being estimated
#
# Written as one linear script, in the style of Crema's example. Nothing here
# feeds the recovery study; it is a separate test on a single dataset.
#
#   Rscript Simulations/Sim_Case1_Radiocarbon/scripts/diagnostics/08_reversed_berkson.R

library(rcarbon)
library(nimbleCarbon)
library(coda)
library(cmdstanr)

out_dir <- "Simulations/Sim_Case1_Radiocarbon/output/reversed"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# The Stan port does not fit this model correctly yet: it returns a slope near
# zero and a sigma three times the truth, while Nimble recovers both. It is kept
# here but switched off, so the figures come from Nimble alone.
#   REVERSED_STAN=1 Rscript .../08_reversed_berkson.R   to run it anyway
run_stan <- Sys.getenv("REVERSED_STAN", "0") == "1"

# Simulate regression data ----
# Crema's values, so the numbers here can be read against his figure 2.
set.seed(12332)
true.alpha <- 2700          # intercept, in cal BP
true.beta  <- 1 / 3.7       # slope, cal BP per unit of d
true.sigma <- 70            # spread of the dates around the trend
n <- 100
d <- round(runif(n = n, min = 0, max = 1000))
calendar.dates <- round(rnorm(n = n, mean = true.alpha - true.beta * d, sd = true.sigma))

# Back-calibrate to a radiocarbon age. Crema takes uncalibrate()[,4], which is a
# draw carrying the curve error only. Case 1 adds the lab error on top, because
# calibrate() assumes both are there, so that is what is used here.
c14ages.error <- rep(30, n)
curve   <- uncalibrate(calendar.dates, verbose = FALSE)
c14ages <- round(rnorm(n, curve$ccCRA, sqrt(c14ages.error^2 + curve$ccError^2)))

median.calibrated <- medCal(calibrate(c14ages, c14ages.error, verbose = FALSE))

# Regression on the true calendar dates, and on the median calibrated dates ----
fitcal <- lm(calendar.dates ~ d)
fitmed <- lm(median.calibrated ~ d)

# Bayesian EIV in Nimble ----
# Crema's model: the date is a latent parameter drawn from the trend,
# and the radiocarbon age is the calibration curve read at that date.
dat <- list(cra = c14ages, cra.error = c14ages.error)

data(intcal20)
constants <- list(d = d, n = n, calBP = intcal20$CalBP,
                  C14BP = intcal20$C14Age, C14err = intcal20$C14Age.sigma)

inits <- list(theta = median.calibrated, alpha = 3000, beta = 1 / 2, sigma = 50)

model <- nimbleCode({
  for (i in 1:n)
  {
    mu[i] <- alpha + beta * d[i]
    theta[i] ~ dnorm(mean = mu[i], sd = sigma)
    c14age[i] <- interpLin(z = theta[i], x = calBP[], y = C14BP[]);
    sigmaCurve[i] <- interpLin(z = theta[i], x = calBP[], y = C14err[]);
    sigmaDate[i] <- (cra.error[i]^2 + sigmaCurve[i]^2)^(1/2);
    cra[i] ~ dnorm(mean = c14age[i], sd = sigmaDate[i])
  }
  alpha ~ dunif(1000, 5000)
  beta ~ dnorm(0, 1)
  sigma ~ dexp(0.1)
})

out <- nimbleMCMC(code = model, constants = constants, inits = inits, data = dat,
                  nchains = 4, niter = 100000, nburnin = 50000,
                  samplesAsCodaMCMC = TRUE,
                  monitors = c('theta', 'alpha', 'beta', 'sigma'))

# Only theta is expected to fail the Rhat check, because a calibrated date is
# not normal and can have several peaks. EC's script says to ignore it.
rhats <- coda::gelman.diag(out, multivariate = FALSE)
cat("Parameters above Rhat 1.01:",
    paste(rownames(rhats[[1]])[rhats[[1]][, 1] > 1.01], collapse = " "), "\n")

posterior <- do.call(rbind.data.frame, out)
theta.nimble <- as.matrix(posterior[, grep('theta', colnames(posterior))])

# The same model in Stan ----
# Stan cannot use a parameter as an index, so the calibration curve is read with
# a hat function: sum(y * max(0, 1 - |position - j|)) over the grid is the same
# linear interpolation interpLin does. The grid is the 5-yr one the other cases
# use, over a range wide enough to hold every date, because the cost of the hat
# function is one pass over the grid for every date at every step.
grid_step  <- 5
grid_years <- seq(1800, 3200, by = grid_step)
grid_c14   <- approx(intcal20$CalBP, intcal20$C14Age, grid_years)$y
grid_err   <- approx(intcal20$CalBP, intcal20$C14Age.sigma, grid_years)$y

# This failed and needs to be fixed: the hat function is not differentiable at the edges, so Stan's autodiff fails. The model is still correct, but the sampler cannot explore the posterior. 
# A fix is to use a smooth approximation to the hat function, such as a Gaussian kernel or a triangular kernel
stan_code <- "
data {
  int<lower=1> n;
  vector[n] cra;
  vector[n] cra_error;
  vector[n] d;
  int<lower=1> G;
  real grid_start;
  real grid_step;
  vector[G] grid_c14;
  vector[G] grid_err;
  int<lower=0,upper=1> fix_first;   // 1 fixes theta[1] to theta_fixed
  real theta_fixed;
}
transformed data {
  vector[G] grid_index;
  for (j in 1:G) grid_index[j] = j;
}
parameters {
  real<lower=1000,upper=5000> alpha;
  real beta;
  real<lower=0> sigma;
  vector<lower=1800,upper=3200>[n - fix_first] theta_free;
}
transformed parameters {
  vector[n] theta;
  vector[n] c14age;
  vector[n] sigma_date;
  if (fix_first == 1) {
    theta[1] = theta_fixed;
    theta[2:n] = theta_free;
  } else {
    theta = theta_free;
  }
  for (i in 1:n) {
    vector[G] w = fmax(0, 1 - abs((theta[i] - grid_start) / grid_step + 1 - grid_index));
    c14age[i] = dot_product(w, grid_c14);
    sigma_date[i] = sqrt(square(cra_error[i]) + square(dot_product(w, grid_err)));
  }
}
model {
  beta ~ normal(0, 1);
  sigma ~ exponential(0.1);
  theta ~ normal(alpha + beta * d, sigma);
  cra ~ normal(c14age, sigma_date);
}
generated quantities {
  // The two terms shared by both versions of the model, so the fixed-date run
  // and the free-date run can be compared on the same scale. lp__ cannot be
  // compared directly: fixing a date removes a parameter.
  real ll_total = normal_lpdf(theta | alpha + beta * d, sigma)
                + normal_lpdf(cra | c14age, sigma_date);
  real ll_first = normal_lpdf(theta[1] | alpha + beta * d[1], sigma)
                + normal_lpdf(cra[1] | c14age[1], sigma_date[1]);
}
"

if (run_stan) {
stan_file <- write_stan_file(stan_code, dir = out_dir, basename = "reversed_eiv")
stan_model <- cmdstan_model(stan_file)

stan_data <- list(n = n, cra = c14ages, cra_error = c14ages.error, d = d,
                  G = length(grid_years), grid_start = min(grid_years),
                  grid_step = grid_step, grid_c14 = grid_c14, grid_err = grid_err,
                  fix_first = 0, theta_fixed = 0)

stan_fit <- stan_model$sample(data = stan_data, chains = 4, parallel_chains = 4,
                              iter_warmup = 1000, iter_sampling = 2000,
                              adapt_delta = 0.95, refresh = 500, seed = 1)

stan_draws   <- stan_fit$draws(format = "df")
theta.stan   <- as.matrix(stan_draws[, grep("^theta\\[", colnames(stan_draws))])
}

# Do the two samplers agree ----
slopes <- data.frame(
  source = c("true calendar dates (lm)", "median calibrated dates (lm)",
             "EIV, Nimble"),
  beta = c(-coef(fitcal)[2], -coef(fitmed)[2], -mean(posterior$beta)),
  lower = c(-confint(fitcal)[2, 2], -confint(fitmed)[2, 2],
            -HPDinterval(as.mcmc(posterior$beta))[2]),
  upper = c(-confint(fitcal)[2, 1], -confint(fitmed)[2, 1],
            -HPDinterval(as.mcmc(posterior$beta))[1])
)
if (run_stan)
  slopes <- rbind(slopes, data.frame(
    source = "EIV, Stan", beta = -mean(stan_draws$beta),
    lower = -HPDinterval(as.mcmc(stan_draws$beta))[2],
    upper = -HPDinterval(as.mcmc(stan_draws$beta))[1]))
# The model is written as alpha + beta * d with dates in cal BP, which run
# backwards, so the fitted beta is negative where the true slope is positive.
# Signs are flipped above so every row is comparable with true.beta.
slopes$ratio <- slopes$beta / true.beta

sigmas <- data.frame(
  source = c("true calendar dates (lm)", "median calibrated dates (lm)",
             "EIV, Nimble"),
  sigma = c(summary(fitcal)$sigma, summary(fitmed)$sigma, mean(posterior$sigma))
)
if (run_stan)
  sigmas <- rbind(sigmas, data.frame(source = "EIV, Stan",
                                     sigma = mean(stan_draws$sigma)))

write.csv(slopes, file.path(out_dir, "slopes.csv"), row.names = FALSE)
write.csv(sigmas, file.path(out_dir, "sigmas.csv"), row.names = FALSE)
print(slopes)
print(sigmas)

# Does the shape of each date change once the trend is modelled ----
# agreementIndex compares the plain calibrated date with the posterior of theta.
# It is the same index OxCal reports, where 60 is the usual threshold.
agr.nimble <- agreementIndex(c14ages, c14ages.error, theta = theta.nimble, verbose = FALSE)

agreement <- data.frame(date = 1:n, nimble = agr.nimble$agreement,
                        cal_median = median.calibrated,
                        post_median = apply(theta.nimble, 2, median),
                        true_date = calendar.dates)
if (run_stan) {
  agr.stan <- agreementIndex(c14ages, c14ages.error, theta = theta.stan, verbose = FALSE)
  agreement$stan <- agr.stan$agreement
  agreement$stan_post_median <- apply(theta.stan, 2, median)
}
agreement$shift <- agreement$post_median - agreement$cal_median
write.csv(agreement, file.path(out_dir, "agreement.csv"), row.names = FALSE)

cat("Overall agreement, Nimble:", round(agr.nimble$overall.agreement, 1), "\n")
if (run_stan)
  cat("Overall agreement, Stan:  ", round(agr.stan$overall.agreement, 1), "\n")
cat("Dates below 60:", sum(agr.nimble$agreement < 60), "\n")

# Fixing one date to a constant ----
# The date closest to the trend line is left free in one run and fixed to a
# round year in the other. If the fixed year is wrong the shared likelihood
# should drop.
# NOTE: this is NOT the check E. asked for at the meeting. That one is about the
# marginal model, where a date is data, and it lives in 09_fixed_date_likelihood.R.
# Here the date is a parameter, so this asks whether a single latent date is
# identified in the reversed model. It needs Stan, so it runs only with
# REVERSED_STAN=1.
if (run_stan) {
target <- which.min(abs(median.calibrated - (true.alpha - true.beta * d)))
fixed_years <- round(c(median(theta.nimble[, target]),
                       median(theta.nimble[, target]) - 200,
                       median(theta.nimble[, target]) + 200))

fixed_runs <- lapply(fixed_years, function(y) {
  dd <- stan_data
  dd$fix_first <- 1
  dd$theta_fixed <- y
  # theta[1] is the fixed one, so the target date is moved to the front
  ord <- c(target, setdiff(1:n, target))
  dd$cra <- dd$cra[ord]; dd$cra_error <- dd$cra_error[ord]; dd$d <- dd$d[ord]
  f <- stan_model$sample(data = dd, chains = 4, parallel_chains = 4,
                         iter_warmup = 1000, iter_sampling = 2000,
                         adapt_delta = 0.95, refresh = 0, seed = 1)
  dr <- f$draws(format = "df")
  data.frame(fixed_year = y, ll_total = mean(dr$ll_total),
             ll_first = mean(dr$ll_first), beta = -mean(dr$beta),
             sigma = mean(dr$sigma))
})

# The same run with the date left free, for the comparison
ord <- c(target, setdiff(1:n, target))
free_data <- stan_data
free_data$cra <- free_data$cra[ord]; free_data$cra_error <- free_data$cra_error[ord]
free_data$d <- free_data$d[ord]
free_fit <- stan_model$sample(data = free_data, chains = 4, parallel_chains = 4,
                              iter_warmup = 1000, iter_sampling = 2000,
                              adapt_delta = 0.95, refresh = 0, seed = 1)
free_draws <- free_fit$draws(format = "df")

fixed <- rbind(
  data.frame(fixed_year = NA, ll_total = mean(free_draws$ll_total),
             ll_first = mean(free_draws$ll_first), beta = -mean(free_draws$beta),
             sigma = mean(free_draws$sigma)),
  do.call(rbind, fixed_runs)
)
fixed$true_date <- calendar.dates[target]
write.csv(fixed, file.path(out_dir, "reversed_fixed_theta.csv"), row.names = FALSE)
print(fixed)
}

# Figures ----
# Crema's figure 2, with the Stan fit added as a fourth panel.
pred <- data.frame(d = -100:1200)
pred$true  <- predict(fitcal, newdata = pred)
pred$m     <- predict(fitmed, newdata = pred)
pred$m.lo  <- predict(fitmed, newdata = pred, interval = 'confidence')[, 2]
pred$m.hi  <- predict(fitmed, newdata = pred, interval = 'confidence')[, 3]

band <- function(a, b) {
  m <- outer(a, rep(1, length(pred$d))) + outer(b, pred$d)
  list(mean = colMeans(m),
       lo = apply(m, 2, function(x) HPDinterval(mcmc(x))[1]),
       hi = apply(m, 2, function(x) HPDinterval(mcmc(x))[2]))
}
b.nimble <- band(posterior$alpha, posterior$beta)
if (run_stan) b.stan <- band(stan_draws$alpha, stan_draws$beta)

cl <- calibrate(c14ages, c14ages.error, calMatrix = TRUE,
                timeRange = c(3000, 2000), verbose = FALSE)

n_panels <- if (run_stan) 4 else 3
png(file.path(out_dir, "reversed_regression.png"), width = 3 * n_panels,
    height = 3.5, units = "in", res = 300)
par(mfrow = c(1, n_panels), mar = c(4, 4, 2, 1), lend = 2)

plot(d, calendar.dates, xlim = c(0, 1000), ylim = c(2850, 2100), pch = 19,
     ylab = 'cal BP', xlab = 'd', col = adjustcolor('black', 0.6),
     main = 'True calendar dates')
abline(a = true.alpha, b = -true.beta, lty = 2, lwd = 2)
lines(pred$d, pred$true, lwd = 2, col = '#780000')

plot(NA, xlim = c(0, 1000), ylim = c(2850, 2100), ylab = 'cal BP', xlab = 'd',
     main = 'What radiocarbon sees')
for (i in 1:n) {
  polygon(c(d[i] + cl$calmatrix[, i] * 1000, d[i] - rev(cl$calmatrix[, i]) * 1000),
          c(3000:2000, 2000:3000), border = NA, col = 'lightblue')
}
segments(d, calendar.dates, d, median.calibrated, lty = 3)
points(d, calendar.dates, pch = 20)
points(d, median.calibrated, pch = 4)

plot(NA, xlim = c(0, 1000), ylim = c(2850, 2100), ylab = 'cal BP', xlab = 'd',
     main = 'Median date against EIV')
abline(a = true.alpha, b = -true.beta, lty = 2, lwd = 2)
polygon(c(pred$d, rev(pred$d)), c(pred$m.lo, rev(pred$m.hi)), border = NA,
        col = adjustcolor('darkorange', 0.3))
lines(pred$d, pred$m, lwd = 2, col = 'darkorange')
polygon(c(pred$d, rev(pred$d)), c(b.nimble$lo, rev(b.nimble$hi)), border = NA,
        col = adjustcolor('darkgreen', 0.3))
lines(pred$d, b.nimble$mean, lwd = 2, col = 'darkgreen')
legend('bottomleft', legend = c('True relationship', 'Median calibrated date', 'EIV, Nimble'),
       lwd = c(1, 8, 8), col = c('black', 'darkorange', 'darkgreen'), bty = 'n',
       cex = 0.9, lty = c(2, 1, 1))

if (run_stan) {
plot(NA, xlim = c(0, 1000), ylim = c(2850, 2100), ylab = 'cal BP', xlab = 'd',
     main = 'Nimble against Stan')
abline(a = true.alpha, b = -true.beta, lty = 2, lwd = 2)
polygon(c(pred$d, rev(pred$d)), c(b.nimble$lo, rev(b.nimble$hi)), border = NA,
        col = adjustcolor('darkgreen', 0.3))
lines(pred$d, b.nimble$mean, lwd = 2, col = 'darkgreen')
lines(pred$d, b.stan$mean, lwd = 2, col = '#780000', lty = 2)
lines(pred$d, b.stan$lo, lwd = 1, col = '#780000', lty = 3)
lines(pred$d, b.stan$hi, lwd = 1, col = '#780000', lty = 3)
legend('bottomleft', legend = c('True relationship', 'EIV, Nimble', 'EIV, Stan'),
       lwd = c(1, 8, 2), col = c('black', 'darkgreen', '#780000'), bty = 'n',
       cex = 0.9, lty = c(2, 1, 2))
}
dev.off()

# The six dates whose shape changed most, plain calibration against posterior
png(file.path(out_dir, "date_shapes.png"), width = 9, height = 6,
    units = "in", res = 300)
par(mfrow = c(2, 3), mar = c(4, 4, 3, 1))
worst <- order(agr.nimble$agreement)[1:6]
for (i in worst) {
  g <- cl$calmatrix[, i]
  years <- 3000:2000
  dens <- density(theta.nimble[, i], from = 2000, to = 3000, n = 1001)
  plot(years, g / max(g), type = 'l', xlim = c(2900, 2100), ylim = c(0, 1.05),
       xlab = 'cal BP', ylab = 'scaled density', col = 'lightblue4',
       main = sprintf('Date %d, agreement %.0f', i, agr.nimble$agreement[i]))
  polygon(c(years, rev(years)), c(g / max(g), rep(0, length(years))),
          border = NA, col = adjustcolor('lightblue', 0.7))
  lines(dens$x, dens$y / max(dens$y), lwd = 2, col = '#780000')
  abline(v = calendar.dates[i], lty = 2)
}
dev.off()

# Agreement against how far the date moved
png(file.path(out_dir, "agreement.png"), width = 6, height = 5, units = "in", res = 300)
par(mar = c(4, 4, 2, 1))
plot(agreement$shift, agreement$nimble, pch = 19, col = adjustcolor('#780000', 0.6),
     xlab = 'posterior median minus calibrated median (yr)', ylab = 'agreement index',
     main = 'Do the modelled dates still agree')
abline(h = 60, lty = 2)
dev.off()

cat("Done. Output in", out_dir, "\n")
