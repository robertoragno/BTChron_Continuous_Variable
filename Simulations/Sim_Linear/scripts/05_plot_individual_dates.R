#' Purpose: Per-sample diagnostic panels (greyscale).
#'   Panel 1 — Linear trend with CI, sample value + date range highlighted
#'   Panel 2 — Posterior date density with TPQ-TAQ and true date lines
#' Adapted from MSCA Model_1_GP 05_diagnostics_on_single_sample.R

library(here)
library(tidyverse)
library(cmdstanr)
library(posterior)
library(patchwork)

# Load data and fit

sim_data <- read_csv(here("Simulations", "Sim_Linear", "data", "simulated_data.csv"),
                     show_col_types = FALSE)

fit <- readRDS(here("Simulations", "Sim_Linear", "output", "fit_linear.rds"))

# Shared theme

theme_panel <- theme_classic(base_size = 11) +
  theme(
    plot.title    = element_text(size = 11, face = "bold", hjust = 0),
    plot.subtitle = element_text(size = 9, colour = "grey40", hjust = 0),
    plot.margin   = margin(5, 10, 2, 10)
  )

# Extract draws

draws <- fit$draws(format = "df")

pred_grid <- seq(min(sim_data$Start_date), max(sim_data$End_date), by = 1)
N_pred    <- length(pred_grid)

trend_cols <- paste0("trend_pred[", 1:N_pred, "]")
trend_mat  <- as.matrix(draws[, trend_cols])

trend_summary <- tibble(
  Year     = pred_grid,
  Median   = apply(trend_mat, 2, median),
  Lower_90 = apply(trend_mat, 2, quantile, 0.05),
  Upper_90 = apply(trend_mat, 2, quantile, 0.95)
)

N_obs <- nrow(sim_data)

# Diagnostic function

diagnose_sample <- function(sample_id) {

  info <- sim_data[sample_id, ]
  label <- paste0(info$ID, " (value = ", round(info$Value, 2), ")")

  # Panel 1: trend + sample highlight
  p_trend <- ggplot(trend_summary) +
    geom_ribbon(aes(x = Year, ymin = Lower_90, ymax = Upper_90),
                fill = "grey80") +
    geom_line(aes(x = Year, y = Median), linewidth = 0.6, colour = "black") +
    annotate("rect",
             xmin = info$Start_date, xmax = info$End_date,
             ymin = -Inf, ymax = Inf,
             fill = "grey50", alpha = 0.15) +
    geom_hline(yintercept = info$Value,
               linetype = "dashed", colour = "black", linewidth = 0.5) +
    scale_x_continuous(expand = c(0.01, 0)) +
    labs(title = paste0(label, " — trend vs observed value"),
         x = "Date (CE)", y = "Value") +
    theme_panel

  # Panel 2: posterior date density
  date_col   <- paste0("true_date_actual[", sample_id, "]")
  date_draws <- draws[[date_col]]

  date_df <- tibble(estimated_date = date_draws)

  p_date <- ggplot(date_df, aes(x = estimated_date)) +
    geom_density(fill = "grey60", colour = "black", linewidth = 0.4, alpha = 0.5) +
    geom_vline(xintercept = info$Start_date,
               linetype = "dashed", colour = "grey40", linewidth = 0.5) +
    geom_vline(xintercept = info$End_date,
               linetype = "dashed", colour = "grey40", linewidth = 0.5) +
    geom_vline(xintercept = info$True_date,
               linetype = "solid", colour = "black", linewidth = 0.6) +
    labs(title = paste0(label, " — posterior date estimate"),
         x = "Estimated date (CE)", y = "Density") +
    theme_panel

  p_trend | p_date
}

# Generate for 6 random samples

set.seed(123)
sample_ids <- sort(sample(1:N_obs, 6))

plots <- map(sample_ids, function(sid) diagnose_sample(sid))

combined <- wrap_plots(plots, ncol = 1) +
  plot_annotation(
    caption = paste0(
      "Left panels: posterior median trend (black line) and 90% CI (grey ribbon). ",
      "Shaded rectangle marks the sample's TPQ–TAQ date range; ",
      "dashed horizontal line marks the observed value.\n",
      "Right panels: posterior density for the estimated date. ",
      "Dashed vertical lines: TPQ–TAQ boundaries. Solid vertical line: true generating date."
    ),
    theme = theme(
      plot.caption = element_text(hjust = 0, size = 8, colour = "grey40")
    )
  )

ggsave(here("Simulations", "Sim_Linear", "figures", "individual_date_posteriors.png"),
       combined, width = 12, height = 18, dpi = 300, bg = "white")

cat("Saved to", here("Simulations", "Sim_Linear", "figures",
                      "individual_date_posteriors.png"), "\n")
