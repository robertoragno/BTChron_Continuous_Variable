# Case 3 tables and figures from output/recovery_results.csv.
# For every group of fits (model x setting):
#   slope bias          mean of posterior median minus true slope
#   calibration slope   estimated slopes regressed on true slopes (1 = one-for-one)
#   90% coverage        how often the 90% interval contains the true slope
#   90% interval width  mean width of that interval
#   bias in sigma       mean of posterior median minus true sigma
# Bars are 95% intervals: +/- 2 standard errors for the means, the lm confidence
# interval for the calibration slope, a Jeffreys interval for coverage.
library(here)
library(ggplot2)

source(here("Simulations", "Sim_Case3_OverlappingPhases", "scripts", "simulate.R"))

results <- read.csv(here("Simulations", "Sim_Case3_OverlappingPhases", "output", "recovery_results.csv"))
design <- read.csv(here("Simulations", "Sim_Case3_OverlappingPhases", "data", "design.csv"))
r <- merge(results, design)
r$model <- ifelse(r$model == "midpoint", "Midpoint / Median", "Full distribution")
fig.dir <- here("Simulations", "Sim_Case3_OverlappingPhases", "figures")

cat("fits:", nrow(r), "| R-hat above 1.01:", sum(r$max_rhat > 1.01),
    "| with a divergent transition:", sum(r$n_divergent > 0), "\n")

# Zero-slope datasets: how often the 90% interval wrongly excludes zero
zero <- r[r$slope_condition == "zero", ]
cat("false-positive rate on zero-slope datasets:\n")
print(tapply(!zero$slope_cov90, zero$model, mean))

# The groups of fits each figure compares, one row per level of the setting
# (sorted, so the figure rows run from the smallest setting to the largest)
core <- r[r$sweep == "core" & r$slope_condition == "random", ]
core <- core[order(core$N), ]
factor.sweep <- r[r$sweep == "factor", ]
overlap.sweep <- factor.sweep[factor.sweep$assign_p == 0.5, ]
overlap.sweep <- overlap.sweep[order(overlap.sweep$overlap), ]
assign.sweep <- factor.sweep[factor.sweep$overlap == 0.5, ]
assign.sweep <- assign.sweep[order(assign.sweep$assign_p), ]
sets <- rbind(data.frame(core, figure = "recovery_summary", level = "all"),
              data.frame(core, figure = "checks/recovery_by_n", level = paste("N =", core$N)),
              data.frame(overlap.sweep, figure = "recovery_by_overlap",
                         level = paste("overlap", overlap.sweep$overlap)),
              data.frame(assign.sweep, figure = "checks/recovery_by_assign_p",
                         level = paste("assign_p", assign.sweep$assign_p)))
metrics <- data.frame()
for (f in unique(sets$figure))
{
	for (l in unique(sets$level[sets$figure == f]))
	{
		for (m in unique(sets$model))
		{
			s <- sets[sets$figure == f & sets$level == l & sets$model == m, ]
			n <- nrow(s)
			k <- sum(s$slope_cov90)
			calib <- lm(slope_med ~ slope, data = s)
			se.bias <- sd(s$slope_err) / sqrt(n)
			se.width <- sd(s$slope_width90) / sqrt(n)
			se.sigma <- sd(s$sigma_err) / sqrt(n)
			metrics <- rbind(metrics, data.frame(
				figure = f, level = l, model = m, n = n,
				metric = c("Slope bias", "Calibration slope", "90% coverage",
				           "90% interval width", "Bias in sigma"),
				value = c(mean(s$slope_err), coef(calib)[2], k / n,
				          mean(s$slope_width90), mean(s$sigma_err)),
				lo = c(mean(s$slope_err) - 2 * se.bias, confint(calib)[2, 1],
				       qbeta(0.025, k + 0.5, n - k + 0.5),
				       mean(s$slope_width90) - 2 * se.width, mean(s$sigma_err) - 2 * se.sigma),
				hi = c(mean(s$slope_err) + 2 * se.bias, confint(calib)[2, 2],
				       qbeta(0.975, k + 0.5, n - k + 0.5),
				       mean(s$slope_width90) + 2 * se.width, mean(s$sigma_err) + 2 * se.sigma)))
		}
	}
}
write.csv(metrics, here("Simulations", "Sim_Case3_OverlappingPhases", "output", "recovery_metrics.csv"),
          row.names = FALSE)
print(metrics[metrics$figure == "recovery_summary", c("model", "metric", "value", "lo", "hi")])

# Figures: the three headline metrics, dotted lines at the value each should hit
headline <- c("Calibration slope", "90% coverage", "Bias in sigma")
targets <- data.frame(metric = factor(headline, levels = headline), target = c(1, 0.9, 0))
for (f in unique(metrics$figure))
{
	x <- metrics[metrics$figure == f & metrics$metric %in% headline, ]
	x$metric <- factor(x$metric, levels = headline)
	x$level <- factor(x$level, levels = unique(x$level))
	p <- ggplot(x, aes(x = value, y = model, colour = model)) +
		geom_vline(data = targets, aes(xintercept = target), linetype = "dotted", colour = "grey60") +
		geom_errorbar(aes(xmin = lo, xmax = hi), width = 0, orientation = "y") +
		geom_point(size = 2.4) +
		scale_colour_manual(values = c("Midpoint / Median" = "grey55", "Full distribution" = "#780000")) +
		facet_grid(level ~ metric, scales = "free_x") +
		scale_x_continuous(n.breaks = 4) +
		labs(title = "Case 3: overlapping phases", x = NULL, y = NULL, colour = NULL) +
		theme_classic() +
		theme(legend.position = "top", axis.text.y = element_blank(), axis.ticks.y = element_blank(),
		      panel.spacing.x = unit(1.5, "lines"))
	ggsave(file.path(fig.dir, paste0(f, ".png")), p, width = 8,
	       height = 2 + 1.4 * length(unique(x$level)), dpi = 300)
}

# Dataset anatomy: one simulated dataset per overlap level. Each find is drawn at
# its value, with its dating window as a grey line and its true date as a red dot.
anatomy <- data.frame()
for (ov in c(0, 0.25, 0.5))
{
	sim <- simulate_overlap(70, 8, 0.02, 1, K = 6, alpha_conc = 1.5, overlap = ov, assign_p = 0.2, seed = 3)
	sim$panel <- paste("overlap", ov)
	anatomy <- rbind(anatomy, sim)
}
p <- ggplot(anatomy) +
	geom_abline(intercept = 8, slope = 0.02, linetype = "dashed") +
	geom_segment(aes(x = Start_date, xend = End_date, y = Value, yend = Value), colour = "grey65") +
	geom_point(aes(x = (Start_date + End_date) / 2, y = Value), shape = 15, size = 1, colour = "grey45") +
	geom_point(aes(x = True_date, y = Value), size = 0.9, colour = "#780000") +
	facet_wrap(~panel) +
	labs(x = "Calendar year", y = "Value",
	     subtitle = "Grey line: dating window. Grey square: midpoint. Red dot: true date. Dashed: true trend.") +
	theme_classic()
ggsave(file.path(fig.dir, "dataset_anatomy.png"), p, width = 9, height = 3.6, dpi = 300)
