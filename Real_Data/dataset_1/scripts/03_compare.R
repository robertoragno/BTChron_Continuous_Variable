# Compare the two dating treatments fitted in 02_fit.R: the trend they imply,
# and their posteriors for the slope and for the noise around it.
#
# Produces: dataset_1_model_comparison.png
library(here)
library(ggplot2)
library(patchwork)
library(posterior)

data.file <- here("Real_Data", "dataset_1", "data", "gini_gb_filtered.csv")
output.dir <- here("Real_Data", "dataset_1", "output")
figure.file <- here("Real_Data", "dataset_1", "figures", "dataset_1_model_comparison.png")

METHOD.LABS <- c(midpoint = "midpoint date", marginal = "whole window")
METHOD.COLS <- c("midpoint date" = "#0b5394", "whole window" = "#780000")

# Negative years are BCE, as in the simulation figures
bce.axis <- scale_x_continuous(
	labels = function(x) ifelse(x < 0, paste0(abs(x), " BCE"), paste(x, "CE")))

panel.theme <- theme_classic(base_size = 11) +
	theme(plot.title = element_text(size = 11, face = "bold", hjust = 0),
	      legend.position = "bottom", legend.title = element_blank(),
	      plot.margin = margin(5, 10, 2, 10))

houses <- read.csv(data.file)
fits <- readRDS(file.path(output.dir, "fits.rds"))

draws <- do.call(rbind, lapply(names(fits), function(m)
{
	d <- as.data.frame(fits[[m]]$draws(c("intercept", "slope", "sigma"), format = "draws_df"))
	data.frame(method = METHOD.LABS[[m]], intercept = d$intercept,
	           slope = d$slope, sigma = d$sigma)
}))

# 1. The two trends over the data ----------------------------------------

grid <- seq(min(houses$Start_date), max(houses$End_date), length.out = 100)
band <- do.call(rbind, lapply(unique(draws$method), function(m)
{
	s <- draws[draws$method == m, ]
	k <- sample(nrow(s), 500)
	mu <- outer(s$intercept[k], rep(1, length(grid))) + outer(s$slope[k], grid)
	data.frame(method = m, date = grid,
	           med = apply(mu, 2, median),
	           lo = apply(mu, 2, quantile, 0.05),
	           hi = apply(mu, 2, quantile, 0.95))
}))

set.seed(1)
shown <- houses[sample(nrow(houses), 300), ]

p.trend <- ggplot() +
	geom_errorbar(data = shown, aes(y = Value, xmin = Start_date, xmax = End_date),
	              orientation = "y", width = 0, linewidth = 0.3,
	              colour = "grey70", alpha = 0.6) +
	geom_ribbon(data = band, aes(x = date, ymin = lo, ymax = hi, fill = method), alpha = 0.25) +
	geom_line(data = band, aes(x = date, y = med, colour = method), linewidth = 0.8) +
	# the key is drawn once, by the density panels below
	scale_colour_manual(values = METHOD.COLS, guide = "none") +
	scale_fill_manual(values = METHOD.COLS, guide = "none") +
	bce.axis +
	labs(title = "Trend in roofed area", x = "date", y = "log roofed area") +
	panel.theme

# 2. Slope and noise posteriors ------------------------------------------

draws$slope100 <- draws$slope * 100

p.slope <- ggplot(draws, aes(x = slope100, fill = method)) +
	geom_density(alpha = 0.5, colour = NA) +
	scale_fill_manual(values = METHOD.COLS) +
	labs(title = "Slope", x = "change in log roofed area per century", y = "density") +
	panel.theme

p.sigma <- ggplot(draws, aes(x = sigma, fill = method)) +
	geom_density(alpha = 0.5, colour = NA) +
	scale_fill_manual(values = METHOD.COLS) +
	labs(title = "Noise around the trend", x = "sigma (log roofed area)", y = "density") +
	panel.theme

panel <- p.trend / (p.slope | p.sigma) +
	plot_layout(heights = c(1.5, 1), guides = "collect") &
	theme(legend.position = "bottom")
ggsave(figure.file, panel, width = 9, height = 8, dpi = 300, bg = "white")
cat("wrote", basename(figure.file), "\n")

print(aggregate(cbind(slope100, sigma) ~ method, draws, function(v) round(quantile(v, c(0.05,0.5,0.95)), 4)))
