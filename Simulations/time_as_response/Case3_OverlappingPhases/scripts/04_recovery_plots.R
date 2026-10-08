# Case 3 (time as response) tables and figures from output/recovery_results.csv.
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

results <- read.csv(here("Simulations", "time_as_response", "Case3_OverlappingPhases", "output", "recovery_results.csv"))
design <- read.csv(here("Simulations", "time_as_response", "Case3_OverlappingPhases", "data", "design.csv"))
r <- merge(results, design)
r$model <- ifelse(r$model == "midpoint", "Midpoint / Median", "Full distribution")
fig.dir <- here("Simulations", "time_as_response", "Case3_OverlappingPhases", "figures")
dir.create(file.path(fig.dir, "checks"), showWarnings = FALSE, recursive = TRUE)

cat("fits:", nrow(r), "| R-hat above 1.01:", sum(r$max_rhat > 1.01),
    "| with a divergent transition:", sum(r$n_divergent > 0), "\n")

# A dataset counts as converged when its full-distribution fit has R-hat below
# 1.05. Its point-date fits are grouped with it, so every group compares the
# same datasets.
full <- r[r$model == "Full distribution", ]
r$converged <- full$max_rhat[match(r$dataset_id, full$dataset_id)] < 1.05
full$converged <- full$max_rhat < 1.05
cat("datasets whose full-distribution fit did not converge:", sum(!full$converged), "of", nrow(full), "\n")
# What kind of datasets they are: size, and the trend's span and noise as shares of the period
full$span_share <- abs(full$slope) * (full$x_max - full$x_min) / (full$period_end - full$period_start)
full$sigma_share <- full$sigma / (full$period_end - full$period_start)
print(aggregate(cbind(N, span_share, sigma_share) ~ converged, data = full, FUN = mean))
print(table(full$slope_condition, full$converged))

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
sets <- rbind(data.frame(core, figure = "recovery_summary", level = ifelse(core$converged, "converged", "not converged")),
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
write.csv(metrics, here("Simulations", "time_as_response", "Case3_OverlappingPhases", "output", "recovery_metrics.csv"),
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
		labs(title = "Case 3: overlapping phases, time as response", x = NULL, y = NULL, colour = NULL) +
		theme_classic() +
		theme(legend.position = "top", axis.text.y = element_blank(), axis.ticks.y = element_blank(),
		      panel.spacing.x = unit(1.5, "lines"))
	ggsave(file.path(fig.dir, paste0(f, ".png")), p, width = 8,
	       height = 2 + 1.4 * length(unique(x$level)), dpi = 300)
}
