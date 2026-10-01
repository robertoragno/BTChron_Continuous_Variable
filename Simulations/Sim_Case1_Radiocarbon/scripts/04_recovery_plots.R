# Case 1 tables and figures from output/recovery_results.csv.
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

results <- read.csv(here("Simulations", "Sim_Case1_Radiocarbon", "output", "recovery_results.csv"))
design <- read.csv(here("Simulations", "Sim_Case1_Radiocarbon", "data", "design.csv"))
r <- merge(results, design)
r$model <- c(midpoint = "Midpoint", median = "Calibrated median", marginal = "Full distribution")[r$model]
r$window <- ifelse(r$window == "plateau", "Hallstatt plateau", "Steep section")
fig.dir <- here("Simulations", "Sim_Case1_Radiocarbon", "figures")

cat("fits:", nrow(r), "| R-hat above 1.01:", sum(r$max_rhat > 1.01),
    "| with a divergent transition:", sum(r$n_divergent > 0), "\n")

# Zero-slope datasets: how often the 90% interval wrongly excludes zero
zero <- r[r$slope_condition == "zero", ]
cat("false-positive rate on zero-slope datasets:\n")
print(tapply(!zero$slope_cov90, zero$model, mean))

# The groups of fits each figure compares, one row per level of the setting
# Headline results use even deposition only, all sweeps, one row per window.
# (sorted, so the figure rows run from the smallest setting to the largest)
reference <- r[r$growth_ratio == 1, ]
core <- reference[reference$sweep == "core" & reference$slope_condition == "random", ]
core <- core[order(core$N), ]
lab <- reference[reference$sweep == "factor", ]
lab <- lab[order(lab$lab_error), ]
# Deposition: even vs 4 times denser at the end, lab error 30, one figure per window
deposition <- r[r$sweep %in% c("factor", "deposition") & r$lab_error == 30 & r$slope_condition == "random", ]
deposition <- deposition[order(deposition$growth_ratio), ]
plateau <- deposition[deposition$window == "Hallstatt plateau", ]
steep <- deposition[deposition$window == "Steep section", ]
sets <- rbind(data.frame(reference, figure = "recovery_summary", level = reference$window),
              data.frame(core, figure = "checks/recovery_by_n", level = paste("N =", core$N)),
              data.frame(lab, figure = "recovery_by_lab_error", level = paste("lab error", lab$lab_error)),
              data.frame(plateau, figure = "checks/recovery_by_deposition_plateau",
                         level = paste("growth ratio", plateau$growth_ratio)),
              data.frame(steep, figure = "checks/recovery_by_deposition_steep",
                         level = paste("growth ratio", steep$growth_ratio)))
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
write.csv(metrics, here("Simulations", "Sim_Case1_Radiocarbon", "output", "recovery_metrics.csv"),
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
		scale_colour_manual(values = c("Midpoint" = "grey55", "Calibrated median" = "grey25", "Full distribution" = "#780000")) +
		facet_grid(level ~ metric, scales = "free_x") +
		labs(title = "Case 1: radiocarbon", x = NULL, y = NULL, colour = NULL) +
		theme_classic() +
		theme(legend.position = "top", axis.text.y = element_blank(), axis.ticks.y = element_blank())
	ggsave(file.path(fig.dir, paste0(f, ".png")), p, width = 8,
	       height = 2 + 1.4 * length(unique(x$level)), dpi = 300)
}

