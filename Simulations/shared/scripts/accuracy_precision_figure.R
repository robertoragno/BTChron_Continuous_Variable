# Accuracy against precision, one panel per direction and method, read from
# each case's output/recovery_metrics.csv. Writes figures/accuracy_precision.png.
#
# x: how far the slope ratio is from 1 (0 = the slope is recovered)
# y: the 90% interval over the size of the true slope (smaller = tighter)
# The bottom left corner is the best a method can do: accurate and tight.
#
# Only the pooled rows are used, one per case (two for Case 1, which behaves
# differently on the plateau and on the steep section). The sweeps within a
# case are left out: their interval width moves with the dataset size, and the
# Case 1 lab-error sweep mixes the two windows.
library(here)
library(ggplot2)

cases <- c("Case1_Radiocarbon", "Case2_Typochronology", "Case3_OverlappingPhases", "Case4_MergedPhases")
x <- data.frame()
for (direction in c("time_as_predictor", "time_as_response"))
{
	for (case in cases)
	{
		m <- read.csv(here("Simulations", direction, case, "output", "recovery_metrics.csv"))
		m <- m[m$figure == "recovery_summary" &
		       m$metric %in% c("Slope ratio", "Relative interval width") &
		       m$level %in% c("Hallstatt plateau", "Steep section", "all", "converged"), ]
		m$direction <- gsub("_", " ", direction)
		m$case <- paste("Case", substr(case, 5, 5))
		m$case[m$level == "Hallstatt plateau"] <- "Case 1, plateau"
		m$case[m$level == "Steep section"] <- "Case 1, steep"
		x <- rbind(x, m)
	}
}

# Case 1 also has the calibrated median, a third method, which the overview
# figure shows; here only the midpoint is kept. In Cases 2-4 the windows are
# flat, so their median is the midpoint and nothing is left out.
x <- x[x$model != "Calibrated median", ]
x$method <- ifelse(x$model == "Full distribution", "Full distribution", "Midpoint")

# One row per case, direction and method, with the two measures side by side
x$group <- paste(x$direction, x$case, x$method, sep = " | ")
d <- data.frame()
for (g in unique(x$group))
{
	s <- x[x$group == g, ]
	d <- rbind(d, data.frame(direction = s$direction[1], case = s$case[1], method = s$method[1],
	                         accuracy = abs(s$value[s$metric == "Slope ratio"] - 1),
	                         precision = s$value[s$metric == "Relative interval width"]))
}
d$method <- factor(d$method, levels = c("Midpoint", "Full distribution"))
print(d[order(d$direction, d$case), ], row.names = FALSE)

p <- ggplot(d, aes(x = accuracy, y = precision, colour = case)) +
	geom_point(size = 3) +
	scale_colour_manual(values = c("Case 1, plateau" = "#780000", "Case 1, steep" = "#bc4749",
	                               "Case 2" = "grey25", "Case 3" = "#c1844d", "Case 4" = "grey60")) +
	# Labels and long ticks every 0.2, a short tick on the grid line between them
	scale_x_continuous(breaks = seq(0, 0.8, 0.2), minor_breaks = seq(0, 0.8, 0.1),
	                   guide = guide_axis(minor.ticks = TRUE)) +
	scale_y_continuous(breaks = seq(0, 0.8, 0.2), minor_breaks = seq(0, 0.8, 0.1),
	                   guide = guide_axis(minor.ticks = TRUE)) +
	facet_grid(direction ~ method) +
	labs(title = "Accuracy against precision",
	     x = "Distance of the slope ratio from 1\n← more accurate",
	     # Rotating the y title flips the order of the lines, so the arrow is
	     # written first to end up outside the title, as it is on the x axis
	     y = "← more precise\nInterval width over the true slope",
	     colour = NULL) +
	theme_classic() +
	theme(legend.position = "top", panel.spacing = unit(1.2, "lines"),
	      panel.grid.major = element_line(colour = "grey90"),
	      panel.grid.minor = element_line(colour = "grey96"),
	      axis.ticks.length = unit(6, "pt"),
	      axis.minor.ticks.length = rel(0.5))
ggsave(here("Simulations", "shared", "figures", "accuracy_precision.png"), p,
       width = 8, height = 6.5, dpi = 300)
