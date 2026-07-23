#' Why do some change-point fits fail to converge? This reads the fits already
#' recorded by 02_recovery_study.R and shows the convergence rate
#' against two properties of each dataset's generating parameters: the resolution
#' of the periodisation (Shannon entropy H), and how detectable the change-point
#' is relative to the noise. The goal is a descriptive explanation of which
#' datasets tend to fail: if convergence mainly falls with coarser periodisations
#' (low H, where one long phase dominates), but not with lower detectability, that
#' suggests non-convergence is driven more by difficult posterior geometry than by
#' intrinsically weak change-point signal. This helps show that the fits dropped by
#' the convergence filter in 03 are not a biased subset defined by undetectable
#' changepoints.

library(ggplot2)
library(patchwork)
library(here)

results <- read.csv(here("Simulations", "Sim_Changepoint", "output",
                         "recovery_results.csv"))
results <- results[is.na(results$error), ]

# A fit counts as converged if Rhat is small and no divergent transitions occurred.
results$converged <- results$max_rhat <= 1.05 & results$n_divergent == 0

# Change-point detectability: slope change x shorter_arm / noise. We multiply by
# shorter_arm because the two slopes separate more as you move away from the
# change-point. We use the shorter side because if the change-point is near one
# edge, there is only a short distance on that side to reveal the bend. We divide
# by noise because the same bend is harder to see when the data are noisier.
shorter_arm <- pmin(results$true_cp - 100, 900 - results$true_cp)
results$detectability <- abs(results$true_slope_2 - results$true_slope_1) *
                         shorter_arm / results$true_sigma

results$model <- factor(results$model, c("latent", "midpoint"),
                        c("EIV", "Midpoint"))

# Convergence rate within quintiles of a chosen property, separately per model.
# Returns the rate and the median property value in each quintile (for the x axis).
rate_by_quintile <- function(values, property_label) {
  binned <- data.frame(
    model     = results$model,
    converged = results$converged,
    value     = values,
    quintile  = cut(values, quantile(values, 0:5 / 5), include.lowest = TRUE)
  )
  # aggregate() summarises one row per quintile/model group: mean(converged) is
  # the proportion that converged (because converged is a logical vector), 
  # median(value) gives the group's x-position,
  # and merge() puts those two summaries together for plotting.
  rate   <- aggregate(converged ~ quintile + model, binned, mean)
  centre <- aggregate(value     ~ quintile + model, binned, median)
  out    <- merge(rate, centre, by = c("quintile", "model"))
  names(out)[names(out) == "quintile"] <- "group"
  out$property <- property_label
  out
}

plot_data <- rbind(
  rate_by_quintile(results$H,             "Shannon entropy of the periodisation (H)"),
  rate_by_quintile(results$detectability, "Change-point detectability (kink / noise)")
)
# Keep entropy as the left panel (it is the property that matters).
plot_data$property <- factor(plot_data$property,
                             c("Shannon entropy of the periodisation (H)",
                               "Change-point detectability (kink / noise)"))

model_colours <- c(`EIV` = "#780000", Midpoint = "grey25")
model_shapes  <- c(`EIV` = 16, Midpoint = 17)

panel_theme <- theme_classic(base_size = 11) +
  theme(plot.title = element_text(size = 12, face = "bold"),
        plot.subtitle = element_text(size = 9),
        legend.position = "top")

build_panel <- function(data, panel_title, x_label, show_y = TRUE, bands = NULL) {
  p <- ggplot(data, aes(value, 100 * converged, colour = model, shape = model)) +
    {if (!is.null(bands)) geom_rect(
      data = bands,
      inherit.aes = FALSE,
      aes(xmin = band_start, xmax = band_end, ymin = -Inf, ymax = Inf, fill = shade)
    )} +
    {if (!is.null(bands)) geom_text(
      data = bands,
      inherit.aes = FALSE,
      aes(x = band_centre, y = 3.5, label = label),
      colour = "grey30",
      size = 3
    )} +
    geom_line(linewidth = 0.7) +
    geom_point(size = 2.2) +
    scale_fill_identity() +
    scale_colour_manual(values = model_colours, name = NULL) +
    scale_shape_manual(values = model_shapes, name = NULL) +
    scale_y_continuous(limits = c(0, 100)) +
    labs(x = x_label,
         y = if (show_y) "Fits that converged (%)" else NULL,
         subtitle = panel_title) +
    panel_theme

  if (!show_y) {
    p <- p + theme(axis.text.y = element_blank(),
                   axis.ticks.y = element_blank())
  }

  p
}

entropy_panel <- build_panel(
  subset(plot_data, property == "Shannon entropy of the periodisation (H)"),
  "Periodisation resolution",
  "Shannon entropy (H)"
)

detect_panel <- build_panel(
  subset(plot_data, property == "Change-point detectability (kink / noise)"),
  "Change-point detectability",
  "Kink / noise",
  show_y = FALSE
)

diagnostic_plot <- wrap_plots(entropy_panel, detect_panel, nrow = 1, guides = "collect") +
  plot_annotation(
    tag_levels = "A",
    title = "Non-convergence tracks coarse periodisation (low H), not weak signal",
    subtitle = "Convergence rate within quintiles of each generating property",
    theme = theme(plot.title = element_text(size = 13, face = "bold"),
                  plot.subtitle = element_text(size = 9),
                  legend.position = "top")
  ) &
  theme(plot.tag = element_text(size = 12, face = "bold"))

ggsave(here("Simulations", "Sim_Changepoint", "figures",
           "convergence_diagnostic.png"),
       diagnostic_plot, width = 9, height = 4.5, dpi = 300, bg = "white")

# Console summary for the manuscript footnote.
# Converged
print(table(results$model, ifelse(results$converged, "kept", "dropped")))

# Percent of fits dropped overall:
print(round(100 * mean(!results$converged), 1))

# Percent of fits dropped by model:
print(round(100 * tapply(!results$converged, results$model, mean), 1))

# Convergence rate by entropy quartile:
entropy_quartile <- cut(results$H, quantile(results$H, 0:4 / 4),
                        include.lowest = TRUE)
print(round(tapply(results$converged, entropy_quartile, mean), 2))