# Do radiocarbon dates still agree with their calibration once the trend is fitted?
# One plateau dataset per direction, fitted with the Stan models using the full
# calibrated distributions. For each date, its posterior date is compared with
# its calibrated date by the agreement index (the one OxCal reports, 60 being
# the usual threshold), computed with nimbleCarbon::agreementIndex.
#   time as response   Crema's measurement_error_example.R dataset (100 dates)
#   time as predictor  100 dates on the plateau, value = 8 + 0.02 * date + noise
# Writes output/agreement_index.csv, figures/checks/agreement_index.png and
# figures/checks/agreement_dates.png (the three dates with the lowest agreement
# in each direction).
library(rcarbon)
library(nimbleCarbon)
library(cmdstanr)
library(ggplot2)
library(here)

source(here("Simulations", "shared", "dating", "case1_radiocarbon.R"))
case.dir <- here("Simulations", "time_as_response", "Case1_Radiocarbon")
step <- 5
n <- 100

# Calibrated probabilities of each date in 5-year cells, as in 02_recovery_study.R
calibrated_cells <- function(cra, error, grid)
{
	cal <- calibrate(cra, error, calMatrix = TRUE, verbose = FALSE,
	                 timeRange = 1950 - c(min(grid) - step, max(grid) + step))
	years <- 1950 - as.numeric(rownames(cal$calmatrix))
	cell <- round((years - grid[1]) / step) + 1
	on.grid <- cell >= 1 & cell <= length(grid)
	prob <- t(rowsum(cal$calmatrix[on.grid, ], cell[on.grid]))
	prob / rowSums(prob)
}

# First and last cell holding 99.99% of each date
date_range <- function(prob)
{
	first <- apply(prob, 1, function(p) which(cumsum(p) >= 0.00005)[1])
	last <- apply(prob, 1, function(p) which(cumsum(p) >= 0.99995)[1])
	list(first = first, last = last)
}

# Time as response: Crema's dataset, as in diagnostics of September 2026
set.seed(12332)
d <- round(runif(n, 0, 1000))
true.bp <- round(rnorm(n, 2700 - d / 3.7, 70))
curve <- uncalibrate(true.bp, verbose = FALSE)
cra.r <- round(rnorm(n, curve$ccCRA, sqrt(30^2 + curve$ccError^2)))
grid.r <- seq(-2600, 400, by = step)
prob.r <- calibrated_cells(cra.r, rep(30, n), grid.r)
range.r <- date_range(prob.r)
model.r <- cmdstan_model(here("Simulations", "shared", "models", "linear_dates_response.stan"))
inits <- lapply(1:4, function(ch) list(a = runif(1, -0.5, 0.5), b = runif(1, -0.001, 0.001), s = runif(1, 0.5, 1.5)))
fit.r <- model.r$sample(data = list(n = n, n_years = length(grid.r), x = d,
                                    year = matrix(grid.r, nrow = n, ncol = length(grid.r), byrow = TRUE),
                                    log_p = log(prob.r), first = range.r$first, last = range.r$last,
                                    centre = -600, x_centre = 500, half_cell = step / 2),
                        chains = 4, parallel_chains = 4, iter_warmup = 500, iter_sampling = 500,
                        adapt_delta = 0.95, init = inits, seed = 1, refresh = 0, show_messages = FALSE)
draws.r <- as.data.frame(fit.r$draws(c("intercept", "slope", "sigma"), format = "draws_df"))

# Posterior date of each find: calibrated probability times the trend's share of
# each cell, averaged over 400 draws
post.r <- matrix(0, n, length(grid.r))
for (k in round(seq(1, nrow(draws.r), length.out = 400)))
{
	gap <- matrix(grid.r, n, length(grid.r), byrow = TRUE) - (draws.r$intercept[k] + draws.r$slope[k] * d)
	w <- prob.r * (pnorm((-abs(gap) + step / 2) / draws.r$sigma[k]) - pnorm((-abs(gap) - step / 2) / draws.r$sigma[k]))
	post.r <- post.r + w / rowSums(w)
}
post.r <- post.r / rowSums(post.r)

# Time as predictor: 100 dates on the plateau
grid.p <- seq(-800 - c14_grid_pad(c(-800, -400), 30), -400 + c14_grid_pad(c(-800, -400), 30), by = step)
sim <- simulate_c14(n, 8, 0.02, 2, c(-800, -400), 30, seed = 1)
prob.p <- calibrated_cells(sim$CRA, sim$Error, grid.p)
range.p <- date_range(prob.p)
model.p <- cmdstan_model(here("Simulations", "shared", "models", "linear_dates.stan"))
fit.p <- model.p$sample(data = list(n = n, n_years = length(grid.p), y = sim$Value,
                                    year = matrix(grid.p, nrow = n, ncol = length(grid.p), byrow = TRUE),
                                    log_p = log(prob.p), first = range.p$first, last = range.p$last,
                                    centre = -600),
                        chains = 4, parallel_chains = 4, iter_warmup = 500, iter_sampling = 500,
                        adapt_delta = 0.95, seed = 1, refresh = 0, show_messages = FALSE)
draws.p <- as.data.frame(fit.p$draws(c("intercept", "slope", "sigma"), format = "draws_df"))

# Posterior date of each find: calibrated probability times how well the trend
# at that year fits the find's value, averaged over 400 draws
post.p <- matrix(0, n, length(grid.p))
for (k in round(seq(1, nrow(draws.p), length.out = 400)))
{
	mu <- draws.p$intercept[k] + draws.p$slope[k] * grid.p
	w <- prob.p * t(sapply(sim$Value, function(y) dnorm(y, mu, draws.p$sigma[k])))
	post.p <- post.p + w / rowSums(w)
}
post.p <- post.p / rowSums(post.p)

# Agreement index: sample 2000 posterior dates per find (spread within each
# 5-year cell) and compare with the calibrated date, in years BP
agreement <- data.frame()
for (direction in c("time as response", "time as predictor"))
{
	if (direction == "time as response") { post <- post.r; grid <- grid.r; cra <- cra.r; true <- 1950 - true.bp }
	if (direction == "time as predictor") { post <- post.p; grid <- grid.p; cra <- sim$CRA; true <- sim$True_date }
	theta <- sapply(1:n, function(i) 1950 - (sample(grid, 2000, replace = TRUE, prob = post[i, ]) + runif(2000, -2.5, 2.5)))
	agr <- agreementIndex(cra, rep(30, n), theta = theta, verbose = FALSE)
	cal <- calibrated_cells(cra, rep(30, n), grid)
	agreement <- rbind(agreement, data.frame(
		direction = direction, date = 1:n, agreement = agr$agreement,
		cal_median = apply(cal, 1, function(p) grid[which(cumsum(p) >= 0.5)[1]]),
		post_median = apply(post, 1, function(p) grid[which(cumsum(p) >= 0.5)[1]]),
		true_date = true))
	cat(direction, ": overall agreement", round(agr$overall.agreement, 1), "| dates below 60:",
	    sum(agr$agreement < 60), "of", n, "\n")
}
agreement$shift <- agreement$post_median - agreement$cal_median
write.csv(agreement, file.path(case.dir, "output", "agreement_index.csv"), row.names = FALSE)

# Agreement against how far the trend moved each date
p <- ggplot(agreement, aes(x = shift, y = agreement)) +
	geom_hline(aes(yintercept = 60, linetype = "OxCal threshold (60)"), colour = "grey40") +
	geom_point(colour = "#780000", alpha = 0.7) +
	scale_linetype_manual(values = c("OxCal threshold (60)" = "dashed")) +
	facet_wrap(~direction) +
	labs(title = "Agreement of each date with its calibration, Hallstatt plateau",
	     x = "Posterior median minus calibrated median (years)", y = "Agreement index", linetype = NULL) +
	theme_classic() +
	theme(legend.position = "top")
ggsave(file.path(case.dir, "figures", "checks", "agreement_index.png"), p, width = 8, height = 4, dpi = 300)

# The three dates with the lowest agreement in each direction: calibrated and posterior date
shapes <- data.frame()
truth <- data.frame()
for (direction in c("time as response", "time as predictor"))
{
	a <- agreement[agreement$direction == direction, ]
	if (direction == "time as response") { post <- post.r; grid <- grid.r; cra <- cra.r }
	if (direction == "time as predictor") { post <- post.p; grid <- grid.p; cra <- sim$CRA }
	cal <- calibrated_cells(cra, rep(30, n), grid)
	for (i in order(a$agreement)[1:3])
	{
		label <- sprintf("%s, date %d (agreement %.0f)", direction, i, a$agreement[i])
		keep <- grid >= -1000 & grid <= -200
		shapes <- rbind(shapes,
			data.frame(panel = label, year = grid[keep], density = cal[i, keep] / max(cal[i, ]), what = "Calibrated date"),
			data.frame(panel = label, year = grid[keep], density = post[i, keep] / max(post[i, ]), what = "Posterior date"))
		truth <- rbind(truth, data.frame(panel = label, year = a$true_date[i]))
	}
}
p <- ggplot(shapes, aes(x = year, y = density)) +
	geom_area(data = shapes[shapes$what == "Calibrated date", ], aes(fill = what), alpha = 0.6) +
	geom_line(data = shapes[shapes$what == "Posterior date", ], aes(colour = what), linewidth = 0.7) +
	geom_vline(data = truth, aes(xintercept = year, linetype = "True date")) +
	scale_fill_manual(values = c("Calibrated date" = "lightblue")) +
	scale_colour_manual(values = c("Posterior date" = "#780000")) +
	scale_linetype_manual(values = c("True date" = "dashed")) +
	facet_wrap(~panel, ncol = 3) +
	labs(title = "Lowest-agreement dates, Hallstatt plateau", x = "Year (BCE negative)",
	     y = "Scaled density", fill = NULL, colour = NULL, linetype = NULL) +
	theme_classic() +
	theme(legend.position = "top", strip.text = element_text(size = 8))
ggsave(file.path(case.dir, "figures", "checks", "agreement_dates.png"), p, width = 10, height = 5.5, dpi = 300)
