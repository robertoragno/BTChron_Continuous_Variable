# Build the analysis subset of the McMillan (2010) pipe stem contexts and the
# exploratory panel. Time as response: the date of a context is the response,
# its mean stem bore the predictor known exactly.
#
# Dropped from the 64 transcribed contexts:
#   3 with no date window   two mean ceramic dates and one documentary year.
#                           A mean ceramic date is a point estimate with error
#                           around it, not a bracket holding the date.
#   1 nested duplicate      St John's Circular Pit appears as "Phase I" and as
#                           "total" over a longer window; the total is kept.
# Flagged but not dropped (02_fit.R fits with and without):
#   2 contexts running past 1780, the limit Binford set for his own formula.
#
# Produces: pipe_stems_subset.csv, dataset_2_exploratory_panel.png
library(here)
library(ggplot2)
library(patchwork)

data.dir <- here("Real_Data", "dataset_2", "data")
figure.file <- here("Real_Data", "dataset_2", "figures", "dataset_2_exploratory_panel.png")

ACCENT <- "#780000"

panel.theme <- theme_classic(base_size = 11) +
	theme(plot.title = element_text(size = 11, face = "bold", hjust = 0),
	      legend.position = "bottom", legend.title = element_blank(),
	      plot.margin = margin(5, 10, 2, 10))

contexts <- read.csv(file.path(data.dir, "pipe_stems_mcmillan2010.csv"))

# 1. Subset ---------------------------------------------------------------

contexts <- contexts[contexts$date_end > contexts$date_start, ]
contexts <- contexts[contexts$context != "Circular Pit Phase I", ]

contexts$width <- contexts$date_end - contexts$date_start
contexts$midpoint <- (contexts$date_start + contexts$date_end) / 2
contexts$past_1780 <- contexts$date_end > 1780

subset.file <- file.path(data.dir, "pipe_stems_subset.csv")
write.csv(contexts, subset.file, row.names = FALSE)
cat("subset:", nrow(contexts), "contexts,", sum(contexts$n), "stems,",
    sum(contexts$past_1780), "running past 1780\n")
cat("windows:", min(contexts$width), "-", max(contexts$width),
    "years, median", median(contexts$width), "\n")

# 2. Exploratory panel ----------------------------------------------------
# The first panel is the data as the model sees it: one segment per context,
# spanning the years its date could take, at the bore measured on it.

contexts$limit <- ifelse(contexts$past_1780, "ends after 1780", "within Binford's limit")

p.data <- ggplot(contexts, aes(x = mean_bore, y = midpoint, colour = limit)) +
	geom_linerange(aes(ymin = date_start, ymax = date_end), linewidth = 0.6) +
	scale_colour_manual(values = c("ends after 1780" = ACCENT,
	                               "within Binford's limit" = "grey30")) +
	labs(title = "Dates and bore diameters", x = "mean bore diameter (1/64 inch)",
	     y = "date window (CE)") +
	panel.theme

p.width <- ggplot(contexts, aes(x = width)) +
	geom_histogram(binwidth = 5, fill = "grey30", colour = "white") +
	labs(title = "Width of the date window", x = "years", y = "contexts") +
	panel.theme

p.stems <- ggplot(contexts, aes(x = n)) +
	geom_histogram(bins = 20, fill = "grey30", colour = "white") +
	scale_x_log10() +
	labs(title = "Measured stems per context", x = "stems (log scale)", y = "contexts") +
	panel.theme

panel <- p.data / (p.width | p.stems) + plot_layout(heights = c(1.6, 1))
ggsave(figure.file, panel, width = 9, height = 8, dpi = 300, bg = "white")
cat("wrote", basename(figure.file), "\n")
