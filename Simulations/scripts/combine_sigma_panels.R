# Cross-simulation "where midpoint fails" figure for the paper.
#
# Each simulation's recovery study writes output/recovery_results.csv with the
# same sigma columns (model, mean_width, sigma_err). This script rebuilds the
# sigma-error-vs-window panel from each of those tables and stacks them into one
# figure with panels A, B, C — one per simulation. Panels are rebuilt from the
# data (not stitched from the per-simulation PNGs) so the styling stays uniform.
#
# Linear and Changepoint are present now; GP is included automatically once its
# recovery_results.csv exists. Missing simulations are skipped.

library(here)
library(readr)
library(dplyr)
library(ggplot2)
library(patchwork)

# Each simulation in display order. label is the panel subtitle.
simulations <- list(
  list(key = "Sim_Linear",      label = "Linear trend"),
  list(key = "Sim_Changepoint", label = "Changepoint trend"),
  list(key = "Sim_GP",          label = "Gaussian-process trend")
)

label_model   <- function(code) factor(code, c("latent", "midpoint"),
                                        c("Latent-date", "Midpoint"))
model_colours <- c(`Latent-date` = "#780000", Midpoint = "grey25")
model_shapes  <- c(`Latent-date` = 16, Midpoint = 17)

window_breaks   <- seq(0, 500, by = 100)
precision_bands <- tibble(
  band_start = window_breaks[-6], band_end = window_breaks[-1],
  band_centre = (window_breaks[-6] + window_breaks[-1]) / 2,
  label = c("Precise", "Fairly precise", "Moderate", "Vague", "Very vague"),
  shade = c("grey88", "grey91", "grey94", "grey97", "white"))

panel_theme <- theme_classic(base_size = 11) +
  theme(plot.title = element_text(size = 12, face = "bold"),
        plot.subtitle = element_text(size = 9),
        strip.background = element_blank(),
        legend.position = "top")

# Read one simulation's converged recovery results, or NULL if not run yet.
read_results <- function(sim) {
  csv <- here("Simulations", sim$key, "output", "recovery_results.csv")
  if (!file.exists(csv)) return(NULL)
  read_csv(csv, show_col_types = FALSE) %>%
    filter(is.na(error), max_rhat <= 1.05, n_divergent == 0) %>%
    mutate(model = label_model(model))
}

# Build one sigma-vs-window panel. y_limits is shared across panels so the
# magnitude of the bias is directly comparable between simulations.
build_panel <- function(sim, results, y_limits) {
  binned <- results %>%
    mutate(bin_centre = (floor(pmin(mean_width, 499) / 100) + 0.5) * 100) %>%
    group_by(model, bin_centre) %>%
    summarise(median_error = median(sigma_err),
              q25 = quantile(sigma_err, 0.25),
              q75 = quantile(sigma_err, 0.75), .groups = "drop") %>%
    mutate(x = bin_centre + ifelse(model == "Latent-date", -12, 12))

  label_y <- y_limits[1] + 0.94 * diff(y_limits)

  ggplot() +
    geom_rect(data = precision_bands,
              aes(xmin = band_start, xmax = band_end, ymin = -Inf, ymax = Inf,
                  fill = shade)) +
    scale_fill_identity() +
    geom_hline(yintercept = 0, linetype = "dotted", colour = "grey45") +
    geom_point(data = results,
               aes(mean_width, sigma_err, colour = model, shape = model),
               size = 0.35, alpha = 0.25) +
    geom_errorbar(data = binned,
                  aes(x = x, ymin = q25, ymax = q75, colour = model),
                  width = 16, linewidth = 0.7) +
    geom_point(data = binned, aes(x = x, y = median_error, colour = model,
                                  shape = model), size = 2) +
    geom_text(data = precision_bands, aes(x = band_centre, y = label_y,
                                          label = label),
              colour = "grey30", size = 2.6) +
    scale_colour_manual(values = model_colours, name = NULL) +
    scale_shape_manual(values = model_shapes, name = NULL) +
    guides(colour = guide_legend(override.aes = list(size = 2.5, alpha = 1))) +
    coord_cartesian(xlim = c(20, 500), ylim = y_limits, expand = FALSE) +
    scale_x_continuous(breaks = seq(100, 400, 100)) +
    labs(x = "Mean timespan (years)",
         y = expression(sigma[err] == sigma[est] - sigma[true]),
         subtitle = sim$label) +
    panel_theme
}

# Read every available simulation, then set one shared y-axis spanning them all.
available <- Filter(function(s) !is.null(s$results),
                    lapply(simulations, function(s) c(s, list(results = read_results(s)))))

if (length(available) == 0) stop("No recovery_results.csv found in any simulation.")

shared_y <- range(sapply(available,
                         function(s) quantile(s$results$sigma_err, c(0.004, 0.996))))

panels <- lapply(available, function(s) build_panel(s, s$results, shared_y))

combined <- wrap_plots(panels, nrow = 1, guides = "collect") +
  plot_annotation(
    tag_levels = "A",
    title = "Recovery of the noise scale under latent and midpoint dating",
    theme = theme(plot.title = element_text(size = 13, face = "bold"),
                  legend.position = "top")
  ) &
  theme(plot.tag = element_text(size = 12, face = "bold"))

out <- here("Simulations", "figures", "sigma_vs_window_combined.png")
ggsave(out, combined, width = 6 * length(panels), height = 5, dpi = 300,
       bg = "white")

cat(sprintf("Wrote %s (%d panel%s: %s)\n", out, length(panels),
            ifelse(length(panels) == 1, "", "s"),
            paste(sapply(available, `[[`, "label"), collapse = ", ")))
