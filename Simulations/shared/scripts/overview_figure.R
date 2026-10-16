# Overview of the recovery studies: slope ratio, slope coverage and the 90%
# interval width for every case and both directions, read from each case's
# output/recovery_metrics.csv. Writes figures/overview.png.
library(here)
library(ggplot2)

cases <- c("Case1_Radiocarbon", "Case2_Typochronology", "Case3_OverlappingPhases", "Case4_MergedPhases")
x <- data.frame()
for (direction in c("time_as_predictor", "time_as_response"))
{
	for (case in cases)
	{
		m <- read.csv(here("Simulations", direction, case, "output", "recovery_metrics.csv"))
		# Headline rows only; the not-converged rows of time as response are left out
		m <- m[m$figure == "recovery_summary" &
		       m$metric %in% c("Slope ratio", "90% coverage", "Relative interval width") &
		       m$level %in% c("Hallstatt plateau", "Steep section", "all", "converged"), ]
		m$case <- paste("Case", substr(case, 5, 5))
		m$case[m$level == "Hallstatt plateau"] <- "Case 1, plateau"
		m$case[m$level == "Steep section"] <- "Case 1, steep"
		m$direction <- gsub("_", " ", direction)
		x <- rbind(x, m)
	}
}
x$model[x$model == "Midpoint / Median"] <- "Midpoint"
x$model <- factor(x$model, levels = c("Midpoint", "Calibrated median", "Full distribution"))
x$metric <- ifelse(x$metric == "90% coverage", "Slope coverage (90%)", x$metric)
x$metric <- factor(x$metric, levels = c("Slope ratio", "Slope coverage (90%)", "Relative interval width"))
x$case <- factor(x$case, levels = rev(c("Case 1, plateau", "Case 1, steep", "Case 2", "Case 3", "Case 4")))
targets <- data.frame(metric = factor(c("Slope ratio", "Slope coverage (90%)"), levels = levels(x$metric)),
                      target = c(1, 0.9))

p <- ggplot(x, aes(x = value, y = case, colour = model)) +
	geom_vline(data = targets, aes(xintercept = target, linetype = "Target"), colour = "grey60") +
	geom_errorbar(aes(xmin = lo, xmax = hi), width = 0, orientation = "y",
	              position = position_dodge(width = 0.6)) +
	geom_point(size = 2, position = position_dodge(width = 0.6)) +
	scale_colour_manual(values = c("Midpoint" = "grey60", "Calibrated median" = "grey25",
	                               "Full distribution" = "#780000")) +
	scale_linetype_manual(values = c("Target" = "dotted")) +
	facet_grid(direction ~ metric, scales = "free_x") +
	labs(title = "Recovery of the slope, all cases", x = NULL, y = NULL, colour = NULL, linetype = NULL,
	     caption = paste("Error bars: regression CI (slope ratio), Jeffreys (coverage),",
	                     "+/- 2 SE (relative width).")) +
	theme_classic() +
	theme(legend.position = "top", panel.spacing.x = unit(1.5, "lines"))
ggsave(here("Simulations", "shared", "figures", "overview.png"), p, width = 11, height = 6, dpi = 300)
