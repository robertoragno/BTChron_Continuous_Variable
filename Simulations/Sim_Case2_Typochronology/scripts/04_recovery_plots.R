# Case 2 tables and figures from output/recovery_results.csv.
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

results <- read.csv(here("Simulations", "Sim_Case2_Typochronology", "output", "recovery_results.csv"))
design <- read.csv(here("Simulations", "Sim_Case2_Typochronology", "data", "design.csv"))
r <- merge(results, design)
r$model <- ifelse(r$model == "midpoint", "Midpoint / Median", "Full distribution")
fig.dir <- here("Simulations", "Sim_Case2_Typochronology", "figures")

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
                         level = paste("share coarse", prop$prop_coarse_samples)),
              data.frame(width, figure = "recovery_by_coarse_width",
                         level = paste("coarse window", width$coarse_frac, "of period")),
              data.frame(deposition, figure = "checks/recovery_by_deposition",
                         level = paste("growth ratio", deposition$growth_ratio)),
              data.frame(precision, figure = "checks/recovery_by_precision",
                         level = paste("precision trend", precision$precision_trend)),
              data.frame(noise[, names(r)], figure = "checks/recovery_by_noise",
                         level = paste("noise ratio", noise$noise_band)))
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
write.csv(metrics, here("Simulations", "Sim_Case2_Typochronology", "output", "recovery_metrics.csv"),
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
		labs(title = "Case 2: finds with mixed dating precision", x = NULL, y = NULL, colour = NULL) +
		theme_classic() +
		theme(legend.position = "top", axis.text.y = element_blank(), axis.ticks.y = element_blank())
	ggsave(file.path(fig.dir, paste0(f, ".png")), p, width = 8,
	       height = 2 + 1.4 * length(unique(x$level)), dpi = 300)
}

