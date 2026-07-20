# Why the slope ties between EIV and midpoint (see README): phase dating is
# Berkson error, not classical measurement error, and only classical error
# biases a linear fit. Classical: a true value blurred by independent noise
# (a bad ruler). Berkson: a window is fixed first, true value uniform inside
# it (a phase) -- what sim_linear.stan actually assumes. Same amount of
# x-uncertainty in each panel below, opposite effect on the fitted slope.
#
# Each fit line is labelled directly via a leader (not a legend), since in
# panel B the two lines nearly coincide. Teal instead of the repo's usual
# EIV-red/Midpoint-grey, since this is about the two error types, not a
# model comparison.
#
# Produces: figures/berkson_vs_classical.png

library(ggplot2)
library(patchwork)

set.seed(42)
N <- 90

true_colour <- "black"
used_colour <- "#1B7A8C"   # teal -- deliberately not the repo's EIV red / Midpoint grey

panel_theme <- theme_classic(base_size = 13) +
  theme(plot.title = element_text(size = 14, face = "bold", hjust = 0),
        plot.margin = margin(8, 16, 8, 8),
        axis.title = element_text(size = 12),
        panel.grid = element_blank(),
        legend.position = "none")

figure_path <- function(name) here::here("Simulations", "Sim_Linear", "figures", name)

# One panel: point cloud (true vs used, joined by thin segments) + two fit
# lines, each with its own leader label.
make_panel <- function(x_true, x_used, y, used_label, xlab, title, subtitle,
                       slope_label, box_x, box_y, x_limits, vlines = NULL,
                       callout_x, label_nudge) {
  fit_true <- lm(y ~ x_true)
  fit_used <- lm(y ~ x_used)
  slope_true <- coef(fit_true)[2]
  slope_used <- coef(fit_used)[2]

  grid_x <- seq(x_limits[1], x_limits[2], length.out = 100)
  line_true <- data.frame(x = grid_x, y = predict(fit_true, data.frame(x_true = grid_x)))
  line_used <- data.frame(x = grid_x, y = predict(fit_used, data.frame(x_used = grid_x)))

  points <- rbind(data.frame(x = x_true, y = y, Type = "true"),
                  data.frame(x = x_used, y = y, Type = "used"))

  # Label position: where each fit sits at callout_x, nudged apart vertically
  # so the two labels never overlap.
  y_true_at <- predict(fit_true, data.frame(x_true = callout_x))
  y_used_at <- predict(fit_used, data.frame(x_used = callout_x))
  label_true_y <- y_true_at + label_nudge
  label_used_y <- y_used_at - label_nudge

  p <- ggplot() +
    { if (!is.null(vlines)) geom_vline(xintercept = vlines, colour = "grey88", linewidth = 0.3) } +
    geom_segment(data = data.frame(x1 = x_true, x2 = x_used, y = y),
                 aes(x = x1, xend = x2, y = y, yend = y),
                 colour = "grey82", linewidth = 0.3) +
    geom_point(data = points[points$Type == "true", ], aes(x, y),
               shape = 16, size = 1.8, colour = true_colour) +
    geom_point(data = points[points$Type == "used", ], aes(x, y),
               shape = 15, size = 1.8, colour = used_colour) +
    geom_line(data = line_true, aes(x, y), linewidth = 0.9, linetype = "dashed",
              colour = true_colour) +
    geom_line(data = line_used, aes(x, y), linewidth = 1, colour = used_colour) +
    # leader segments, label back to line
    annotate("segment", x = callout_x, xend = callout_x,
             y = y_true_at, yend = label_true_y - sign(label_nudge) * 0.6,
             colour = true_colour, linewidth = 0.4) +
    annotate("segment", x = callout_x, xend = callout_x,
             y = y_used_at, yend = label_used_y + sign(label_nudge) * 0.6,
             colour = used_colour, linewidth = 0.4) +
    annotate("text", x = callout_x, y = label_true_y, label = "true values",
             colour = true_colour, fontface = "bold", size = 3.6, hjust = 0.5) +
    annotate("text", x = callout_x, y = label_used_y, label = used_label,
             colour = used_colour, fontface = "bold", size = 3.6, hjust = 0.5) +
    annotate("label", x = box_x, y = box_y, hjust = 0,
             fill = "white", alpha = 0.95, size = 4.1,
             label = sprintf(slope_label, slope_true, slope_used,
                             100 * slope_used / slope_true)) +
    coord_cartesian(xlim = x_limits, clip = "off") +
    labs(title = title, subtitle = subtitle, x = xlab, y = "y") +
    panel_theme
  list(plot = p, slope_true = slope_true, slope_used = slope_used)
}

# Classical error: true x blurred by independent noise.
x_true_c <- runif(N, 0, 12)
y_c      <- 2 + 3 * x_true_c + rnorm(N, 0, 2)
x_obs_c  <- x_true_c + rnorm(N, 0, 4)

res_a <- make_panel(
  x_true_c, x_obs_c, y_c, used_label = "observed values", xlab = "x",
  title = "A. Classical error", subtitle = "true value is known, but blurred by noise",
  slope_label = "True-date slope: %.2f\nObserved slope: %.2f  (%.0f%% of true)",
  box_x = -1.5, box_y = 44, x_limits = c(-2, 16),
  callout_x = 10.5, label_nudge = 7)

# Berkson error: a phase window is fixed first, true x uniform inside it.
bounds <- seq(0, 30, by = 5)
mid    <- (bounds[-length(bounds)] + bounds[-1]) / 2
win    <- sample(seq_len(6), N, replace = TRUE)
x_true_b <- runif(N, bounds[win], bounds[win + 1])
y_b      <- 2 + 3 * x_true_b + rnorm(N, 0, 2)
x_mid_b  <- mid[win]

res_b <- make_panel(
  x_true_b, x_mid_b, y_b, used_label = "phase midpoints", xlab = "x",
  title = "B. Berkson error with phase dating",
  subtitle = "known window (grey lines), true value random inside it",
  slope_label = "True-date slope: %.2f\nMidpoint slope: %.2f  (%.0f%% of true)",
  box_x = -1.5, box_y = 88, x_limits = c(-2, 32), vlines = bounds,
  callout_x = 15, label_nudge = 10)

final <- res_a$plot | res_b$plot

ggsave(figure_path("berkson_vs_classical.png"), final,
       width = 12.5, height = 5.8, dpi = 300, bg = "white")

cat(sprintf("Classical: true-date slope %.2f, observed slope %.2f (%.0f%% of true)\n",
            res_a$slope_true, res_a$slope_used, 100 * res_a$slope_used / res_a$slope_true))
cat(sprintf("Berkson:   true-date slope %.2f, midpoint slope %.2f (%.0f%% of true)\n",
            res_b$slope_true, res_b$slope_used, 100 * res_b$slope_used / res_b$slope_true))
