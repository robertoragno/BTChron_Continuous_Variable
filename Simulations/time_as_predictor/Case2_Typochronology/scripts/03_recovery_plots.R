# Case 2 tables and figures from output/recovery_results.csv.
# For every group of fits (model x setting):
#   slope bias          mean of posterior median minus true slope
#   slope ratio   estimated slopes regressed on true slopes (1 = one-for-one)
#   90% coverage        how often the 90% interval contains the true slope
#   90% interval width  mean width of that interval
#   relative interval width  that width over the size of the true slope
#   bias in sigma       mean of posterior median minus true sigma
# Bars are 95% intervals: +/- 2 standard errors for the means, the lm confidence
# interval for the slope ratio, a Jeffreys interval for coverage.
library(here)
library(ggplot2)

source(here("Simulations", "shared", "dating", "case2_typochronology.R"))

results <- read.csv(here("Simulations", "time_as_predictor", "Case2_Typochronology", "output", "recovery_results.csv"))
design <- read.csv(here("Simulations", "time_as_predictor", "Case2_Typochronology", "data", "design.csv"))
r <- merge(results, design)
r$model <- ifelse(r$model == "midpoint", "Midpoint / Median", "Full distribution")
fig.dir <- here("Simulations", "time_as_predictor", "Case2_Typochronology", "figures")

cat("fits:", nrow(r), "| R-hat above 1.01:", sum(r$max_rhat > 1.01),
    "| with a divergent transition:", sum(r$n_divergent > 0), "\n")

# Zero-slope datasets: how often the 90% interval wrongly excludes zero
zero <- r[r$slope_condition == "zero", ]
cat("false-positive rate on zero-slope datasets:\n")
print(tapply(!zero$slope_cov90, zero$model, mean))

# The groups of fits each figure compares, one row per level of the setting
# Each factor sweep is read against the reference datasets (N = 200, random slope),
# which sit in the core sweep.
# (sorted, so the figure rows run from the smallest setting to the largest)
core <- r[r$sweep == "core" & r$slope_condition == "random", ]
core <- core[order(core$N), ]
reference <- core[core$N == 200, ]
prop <- rbind(reference, r[r$sweep == "prop", ])
prop <- prop[order(prop$prop_coarse_samples), ]
width <- rbind(reference, r[r$sweep == "width", ])
width <- width[order(width$coarse_frac), ]
deposition <- rbind(reference, r[r$sweep == "deposition", ])
deposition <- deposition[order(deposition$growth_ratio), ]
precision <- rbind(reference, r[r$sweep == "precision", ])
precision <- precision[order(precision$precision_trend), ]
# Noise compared with the total change the trend makes across the period
core$noise_ratio <- core$sigma / (abs(core$slope) * (core$period_end - core$period_start))
core$noise_band <- cut(core$noise_ratio, c(0, 0.1, 0.3, 1, Inf), labels = c("< 0.1", "0.1-0.3", "0.3-1", "> 1"))
noise <- core[order(core$noise_band), ]
sets <- rbind(data.frame(core[, names(r)], figure = "recovery_summary", level = "all"),
              data.frame(core[, names(r)], figure = "checks/recovery_by_n", level = paste("N =", core$N)),
              data.frame(prop, figure = "recovery_by_prop_coarse",
                         level = paste(100 * prop$prop_coarse_samples, "% coarse", sep = "")),
              data.frame(width, figure = "recovery_by_coarse_width",
                         level = paste("width", width$coarse_frac)),
              data.frame(deposition, figure = "checks/recovery_by_deposition",
                         level = paste("growth", deposition$growth_ratio)),
              data.frame(precision, figure = "checks/recovery_by_precision",
                         level = paste("trend", precision$precision_trend)),
              data.frame(noise[, names(r)], figure = "checks/recovery_by_noise",
                         level = paste("noise", noise$noise_band)))
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
			# Dividing the width by the size of the true slope removes its units
			# and the length of the period, so the width can be read next to the
			# other cases and the other direction.
			slope.scale <- mean(abs(s$slope))
			metrics <- rbind(metrics, data.frame(
				figure = f, level = l, model = m, n = n,
				metric = c("Slope bias", "Slope ratio", "90% coverage",
				           "90% interval width", "Relative interval width", "Bias in sigma"),
				value = c(mean(s$slope_err), coef(calib)[2], k / n,
				          mean(s$slope_width90), mean(s$slope_width90) / slope.scale,
				          mean(s$sigma_err)),
				lo = c(mean(s$slope_err) - 2 * se.bias, confint(calib)[2, 1],
				       qbeta(0.025, k + 0.5, n - k + 0.5),
				       mean(s$slope_width90) - 2 * se.width,
				       (mean(s$slope_width90) - 2 * se.width) / slope.scale,
				       mean(s$sigma_err) - 2 * se.sigma),
				hi = c(mean(s$slope_err) + 2 * se.bias, confint(calib)[2, 2],
				       qbeta(0.975, k + 0.5, n - k + 0.5),
				       mean(s$slope_width90) + 2 * se.width,
				       (mean(s$slope_width90) + 2 * se.width) / slope.scale,
				       mean(s$sigma_err) + 2 * se.sigma)))
		}
	}
}
write.csv(metrics, here("Simulations", "time_as_predictor", "Case2_Typochronology", "output", "recovery_metrics.csv"),
          row.names = FALSE)
print(metrics[metrics$figure == "recovery_summary", c("model", "metric", "value", "lo", "hi")])

# Figures: the three headline metrics, dotted lines at the value each should hit
headline <- c("Slope ratio", "90% coverage", "Bias in sigma")
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
		labs(title = "Case 2: finds with mixed dating precision", x = NULL, y = NULL, colour = NULL) +
		theme_classic() +
		theme(legend.position = "top", axis.text.y = element_blank(), axis.ticks.y = element_blank(),
		      panel.spacing.x = unit(1.5, "lines"))
	ggsave(file.path(fig.dir, paste0(f, ".png")), p, width = 8,
	       height = 2 + 1.4 * length(unique(x$level)), dpi = 300)
}

# Dataset anatomy: one simulated dataset per proportion of coarsely dated finds,
# before any fitting. Each find is drawn at its value, with its window, its
# midpoint and its true date.
anatomy <- data.frame()
for (share in c(0, 0.5, 1))
{
	sim <- simulate_typo(60, 8, 0.02, 1, prop_coarse_samples = share, seed = 7)
	sim$panel <- paste0(100 * share, "% coarsely dated")
	anatomy <- rbind(anatomy, sim)
}
anatomy$panel <- factor(anatomy$panel, levels = unique(anatomy$panel))
p <- ggplot(anatomy) +
	geom_abline(aes(intercept = 8, slope = 0.02, linetype = "True trend")) +
	geom_segment(aes(x = Start_date, xend = End_date, y = Value, yend = Value, colour = "Dating window")) +
	geom_point(aes(x = (Start_date + End_date) / 2, y = Value, colour = "Midpoint"), shape = 15, size = 1) +
	geom_point(aes(x = True_date, y = Value, colour = "True date"), size = 0.9) +
	scale_colour_manual(values = c("Dating window" = "grey70", "Midpoint" = "grey40", "True date" = "#780000")) +
	scale_linetype_manual(values = c("True trend" = "dashed")) +
	facet_wrap(~panel) +
	labs(title = "Case 2: one dataset per proportion of coarsely dated finds",
	     x = "Calendar year", y = "Value", colour = NULL, linetype = NULL) +
	theme_classic() +
	theme(legend.position = "top")
ggsave(file.path(fig.dir, "dataset_anatomy.png"), p, width = 9, height = 3.8, dpi = 300)
