# Build the analysis subset of the GINI house database and the exploratory
# panel. Time as predictor: the roofed area of a house is the response, its
# date the predictor, known only as a window.
#
# Subset: Great Britain, nuclear-family houses, positive roofed area, a dated
# window (EndHouse > BeginHouse) and a midpoint in [-100, 1000].
#
# Produces: gini_gb_filtered.csv, dataset_1_exploratory_panel.png
library(here)
library(ggplot2)
library(patchwork)

data.dir <- here("Real_Data", "dataset_1", "data")
figure.file <- here("Real_Data", "dataset_1", "figures", "dataset_1_exploratory_panel.png")

panel.theme <- theme_classic(base_size = 11) +
	theme(plot.title = element_text(size = 11, face = "bold", hjust = 0),
	      plot.margin = margin(5, 10, 2, 10))

raw <- read.csv(file.path(data.dir, "gini_database_all_records_20240721.csv"))

# 1. Subset ---------------------------------------------------------------

# which() rather than a plain logical index: a record missing a date or an area
# gives NA in the condition, which would otherwise be kept as a row of NAs
houses <- raw[which(raw$Region == "Great Britain" & raw$RoofedArea > 0 &
                    raw$Housetypedesc == "Nuclear family" & raw$EndHouse > raw$BeginHouse), ]
houses$Midpoint <- (houses$BeginHouse + houses$EndHouse) / 2
houses <- houses[which(houses$Midpoint >= -100 & houses$Midpoint <= 1000), ]

houses <- data.frame(ID = seq_len(nrow(houses)), Site = houses$Site,
                     Start_date = houses$BeginHouse, End_date = houses$EndHouse,
                     Midpoint = houses$Midpoint, RoofedArea = houses$RoofedArea,
                     Value = log(houses$RoofedArea))
houses$width <- houses$End_date - houses$Start_date

write.csv(houses, file.path(data.dir, "gini_gb_filtered.csv"), row.names = FALSE)
cat("subset:", nrow(houses), "houses from", length(unique(houses$Site)), "sites\n")
cat("span:", min(houses$Start_date), "-", max(houses$End_date), "CE\n")
cat("windows:", min(houses$width), "-", max(houses$width),
    "years, median", median(houses$width), "\n")

# 2. Exploratory panel ----------------------------------------------------
# The first panel is the data as the model sees it: one segment per house,
# spanning the years its date could take, at the roofed area measured for it.
# Too many houses to draw them all, so a sample is shown.

set.seed(1)
shown <- houses[sample(nrow(houses), 300), ]

p.data <- ggplot(shown, aes(y = Value)) +
	geom_errorbar(aes(xmin = Start_date, xmax = End_date), orientation = "y",
	              width = 0, linewidth = 0.3, colour = "grey30", alpha = 0.6) +
	labs(title = "Dates and roofed areas (300 houses sampled)",
	     x = "date window (CE)", y = "log roofed area") +
	panel.theme

p.width <- ggplot(houses, aes(x = width)) +
	geom_histogram(binwidth = 25, fill = "grey30", colour = "white") +
	labs(title = "Width of the date window", x = "years", y = "houses") +
	panel.theme

p.value <- ggplot(houses, aes(x = Value)) +
	geom_histogram(bins = 30, fill = "grey30", colour = "white") +
	labs(title = "Roofed area", x = "log roofed area", y = "houses") +
	panel.theme

panel <- p.data / (p.width | p.value) + plot_layout(heights = c(1.6, 1))
ggsave(figure.file, panel, width = 9, height = 8, dpi = 300, bg = "white")
cat("wrote", basename(figure.file), "\n")
