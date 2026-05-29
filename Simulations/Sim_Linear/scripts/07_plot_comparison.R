#' Purpose: Side-by-side comparison of the latent-date model vs the midpoint model.
#'   Row 1 — Trend recovery (both models)
#'   Row 2 — Posterior parameter distributions (both models)
#'   Row 3 — Date recovery (latent model only, since midpoint has no date inference)

library(here)
library(tidyverse)
library(cmdstanr)
library(posterior)
library(patchwork)

# ── Load data ────────────────────────────────────────────────────────────────

sim_data     <- read_csv(here("Simulations", "Sim_Linear", "data", "simulated_data.csv"),
                         show_col_types = FALSE)
ground_truth <- read_csv(here("Simulations", "Sim_Linear", "data", "ground_truth.csv"),
                         show_col_types = FALSE)
true_params  <- read_csv(here("Simulations", "Sim_Linear", "data", "generating_parameters.csv"),
                         show_col_types = FALSE)

fit_latent  <- readRDS(here("Simulations", "Sim_Linear", "output", "fit_linear.rds"))
fit_midpoint <- readRDS(here("Simulations", "Sim_Linear", "output", "fit_midpoint.rds"))

true_baseline <- true_params$value[true_params$parameter == "baseline"]
true_slope    <- true_params$value[true_params$parameter == "slope"]
true_sigma    <- true_params$value[true_params$parameter == "sigma_noise"]

# ── Shared theme ─────────────────────────────────────────────────────────────

theme_panel <- theme_classic(base_size = 11) +
  theme(
    plot.title   = element_text(size = 11, face = "bold", hjust = 0),
    plot.margin  = margin(5, 10, 2, 10)
  )

# ── Helper: extract trend summary from a fit ─────────────────────────────────

pred_grid <- seq(min(sim_data$Start_date), max(sim_data$End_date), by = 1)
N_pred    <- length(pred_grid)

extract_trend <- function(fit, model_label) {
  draws     <- fit$draws(format = "df")
  trend_mat <- as.matrix(draws[, paste0("trend_pred[", 1:N_pred, "]")])
  tibble(
    Year     = pred_grid,
    Median   = apply(trend_mat, 2, median),
    Lower_90 = apply(trend_mat, 2, quantile, 0.05),
    Upper_90 = apply(trend_mat, 2, quantile, 0.95),
    Lower_50 = apply(trend_mat, 2, quantile, 0.25),
    Upper_50 = apply(trend_mat, 2, quantile, 0.75),
    Model    = model_label
  )
}

trend_both <- bind_rows(
  extract_trend(fit_latent,   "Latent dates"),
  extract_trend(fit_midpoint, "Midpoint dates")
)

# ── Helper: extract parameter posteriors ─────────────────────────────────────

extract_params <- function(fit, model_label) {
  draws <- fit$draws(format = "df")
  tibble(
    Baseline = draws$baseline_original,
    Slope    = draws$slope_original,
    Sigma    = draws$sigma
  ) %>%
    pivot_longer(everything(), names_to = "Parameter", values_to = "Value") %>%
    mutate(
      Parameter = factor(Parameter, levels = c("Baseline", "Slope", "Sigma")),
      Model     = model_label
    )
}

params_both <- bind_rows(
  extract_params(fit_latent,   "Latent dates"),
  extract_params(fit_midpoint, "Midpoint dates")
)

# ── Row 1: Trend recovery ───────────────────────────────────────────────────

make_trend_panel <- function(df, title_label) {
  ggplot(df) +
    geom_ribbon(aes(x = Year, ymin = Lower_90, ymax = Upper_90), fill = "grey80") +
    geom_ribbon(aes(x = Year, ymin = Lower_50, ymax = Upper_50), fill = "grey60") +
    geom_line(aes(x = Year, y = Median), linewidth = 0.6, colour = "black") +
    geom_line(data = ground_truth, aes(x = Year, y = True_value),
              linewidth = 0.7, colour = "black", linetype = "dashed") +
    scale_x_continuous(breaks = seq(100, 900, 200), expand = c(0.01, 0)) +
    labs(title = title_label, x = "Year (CE)", y = "Value") +
    theme_panel
}

p_trend_latent <- make_trend_panel(
  trend_both %>% filter(Model == "Latent dates"),
  expression(bold("A.") ~ "Trend recovery — latent dates")
)

p_trend_mid <- make_trend_panel(
  trend_both %>% filter(Model == "Midpoint dates"),
  expression(bold("B.") ~ "Trend recovery — midpoint dates")
)

# ── Row 2: Parameter posteriors (overlaid) ───────────────────────────────────

true_lines <- tibble(
  Parameter = factor(c("Baseline", "Slope", "Sigma"),
                     levels = c("Baseline", "Slope", "Sigma")),
  True      = c(true_baseline, true_slope, true_sigma)
)

p_params <- ggplot(params_both, aes(x = Value, fill = Model)) +
  geom_density(alpha = 0.45, colour = "black", linewidth = 0.3) +
  geom_vline(data = true_lines, aes(xintercept = True),
             linetype = "dashed", linewidth = 0.6, colour = "black") +
  scale_fill_manual(
    name   = NULL,
    values = c("Latent dates" = "grey40", "Midpoint dates" = "grey80")
  ) +
  facet_wrap(~ Parameter, scales = "free", nrow = 1) +
  labs(
    title = expression(bold("C.") ~ "Posterior distributions — both models"),
    x = NULL, y = "Density"
  ) +
  theme_panel +
  theme(
    strip.background = element_blank(),
    strip.text       = element_text(face = "bold", size = 10),
    legend.position  = "bottom",
    legend.text      = element_text(size = 9),
    legend.key.size  = unit(12, "pt")
  )

# ── Row 3: Date recovery (latent model only) ────────────────────────────────

N_obs     <- nrow(sim_data)
draws_lat <- fit_latent$draws(format = "df")
date_mat  <- as.matrix(draws_lat[, paste0("true_date_actual[", 1:N_obs, "]")])

date_recovery <- tibble(
  True_date    = sim_data$True_date,
  Inferred_med = apply(date_mat, 2, median),
  Inferred_lo  = apply(date_mat, 2, quantile, 0.05),
  Inferred_hi  = apply(date_mat, 2, quantile, 0.95)
)

p_dates <- ggplot(date_recovery, aes(x = True_date, y = Inferred_med)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed",
              colour = "black", linewidth = 0.5) +
  geom_linerange(aes(ymin = Inferred_lo, ymax = Inferred_hi),
                 linewidth = 0.15, colour = "grey50", alpha = 0.5) +
  geom_point(shape = 16, size = 0.5, colour = "black", alpha = 0.6) +
  scale_x_continuous(breaks = seq(100, 900, 200), expand = c(0.01, 0)) +
  scale_y_continuous(breaks = seq(100, 900, 200), expand = c(0.01, 0)) +
  coord_equal() +
  labs(
    title = expression(bold("D.") ~ "Date recovery — latent model"),
    subtitle = "Each point = one sample. Dashed line = perfect 1:1 recovery.",
    x = "True generating date (CE)",
    y = "Posterior median date (CE)"
  ) +
  theme_panel +
  theme(plot.subtitle = element_text(size = 9, colour = "grey40", hjust = 0))

# ── Combine ──────────────────────────────────────────────────────────────────

p <- (p_trend_latent | p_trend_mid) /
  p_params /
  p_dates +
  plot_layout(heights = c(1, 0.7, 1)) +
  plot_annotation(
    caption = paste0(
      "Panels A–B: posterior median (solid) and 50%/90% credible intervals ",
      "vs true generating trend (dashed). ",
      "Panel C: overlaid posterior densities; dashed lines = true values. ",
      "Panel D: each dot is one sample; vertical bars = 90% CI on the inferred date."
    ),
    theme = theme(
      plot.caption = element_text(hjust = 0, size = 8, colour = "grey40",
                                  lineheight = 1.2)
    )
  )

ggsave(here("Simulations", "Sim_Linear", "figures", "model_comparison.png"),
       p, width = 12, height = 12, dpi = 300, bg = "white")

cat("Saved to", here("Simulations", "Sim_Linear", "figures",
                      "model_comparison.png"), "\n")
