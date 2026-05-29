#' Purpose: Plot model results — trend recovery, latent-date recovery,
#'          and posterior distributions of the generating parameters.

library(here)
library(tidyverse)
library(cmdstanr)
library(posterior)
library(patchwork)

# ── Load data and fit ────────────────────────────────────────────────────────

sim_data     <- read_csv(here("Simulations", "Sim_Linear", "data", "simulated_data.csv"),
                         show_col_types = FALSE)
ground_truth <- read_csv(here("Simulations", "Sim_Linear", "data", "ground_truth.csv"),
                         show_col_types = FALSE)
true_params  <- read_csv(here("Simulations", "Sim_Linear", "data", "generating_parameters.csv"),
                         show_col_types = FALSE)

fit <- readRDS(here("Simulations", "Sim_Linear", "output", "fit_linear.rds"))

# ── Shared theme (matches exploratory panel) ─────────────────────────────────

theme_panel <- theme_classic(base_size = 11) +
  theme(
    plot.title   = element_text(size = 11, face = "bold", hjust = 0),
    plot.margin  = margin(5, 10, 2, 10),
    axis.title.x = element_blank()
  )

# ── Extract draws ────────────────────────────────────────────────────────────

draws <- fit$draws(format = "df")

pred_grid <- seq(min(sim_data$Start_date), max(sim_data$End_date), by = 1)
N_pred    <- length(pred_grid)

trend_cols <- paste0("trend_pred[", 1:N_pred, "]")
trend_mat  <- as.matrix(draws[, trend_cols])

trend_summary <- tibble(
  Year     = pred_grid,
  Median   = apply(trend_mat, 2, median),
  Lower_90 = apply(trend_mat, 2, quantile, 0.05),
  Upper_90 = apply(trend_mat, 2, quantile, 0.95),
  Lower_50 = apply(trend_mat, 2, quantile, 0.25),
  Upper_50 = apply(trend_mat, 2, quantile, 0.75)
)

# ── Panel A: Recovered trend vs true trend ───────────────────────────────────

p_trend <- ggplot(trend_summary) +
  geom_ribbon(aes(x = Year, ymin = Lower_90, ymax = Upper_90, fill = "90% CI")) +
  geom_ribbon(aes(x = Year, ymin = Lower_50, ymax = Upper_50, fill = "50% CI")) +
  geom_line(aes(x = Year, y = Median), linewidth = 0.6, colour = "black") +
  geom_line(data = ground_truth, aes(x = Year, y = True_value),
            linewidth = 0.7, colour = "black", linetype = "dashed") +
  geom_point(data = sim_data,
             aes(x = (Start_date + End_date) / 2, y = Value),
             shape = 16, size = 0.4, colour = "grey40", alpha = 0.4) +
  scale_fill_manual(
    name   = NULL,
    values = c("90% CI" = "grey80", "50% CI" = "grey60"),
    guide  = guide_legend(override.aes = list(alpha = 1))
  ) +
  scale_x_continuous(breaks = seq(100, 900, 100), expand = c(0.01, 0)) +
  labs(
    title = expression(bold("A.") ~ "Recovered trend vs true trend"),
    x = "Year CE",
    y = "Value"
  ) +
  theme_panel +
  theme(
    axis.title.x  = element_text(),
    legend.position = c(0.05, 0.95),
    legend.justification = c(0, 1),
    legend.background = element_rect(fill = alpha("white", 0.8), colour = NA),
    legend.key.size = unit(10, "pt"),
    legend.text = element_text(size = 8)
  )

# ── Panel B: Latent date recovery ────────────────────────────────────────────

N_obs     <- nrow(sim_data)
date_cols <- paste0("true_date_actual[", 1:N_obs, "]")
date_mat  <- as.matrix(draws[, date_cols])

date_recovery <- tibble(
  True_date      = sim_data$True_date,
  Inferred_med   = apply(date_mat, 2, median),
  Inferred_lo_90 = apply(date_mat, 2, quantile, 0.05),
  Inferred_hi_90 = apply(date_mat, 2, quantile, 0.95)
)

p_dates <- ggplot(date_recovery, aes(x = True_date, y = Inferred_med)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed",
              colour = "black", linewidth = 0.5) +
  geom_linerange(aes(ymin = Inferred_lo_90, ymax = Inferred_hi_90),
                 linewidth = 0.15, colour = "grey50", alpha = 0.5) +
  geom_point(shape = 16, size = 0.5, colour = "black", alpha = 0.6) +
  scale_x_continuous(breaks = seq(100, 900, 100), expand = c(0.01, 0)) +
  scale_y_continuous(breaks = seq(100, 900, 100), expand = c(0.01, 0)) +
  coord_equal() +
  labs(
    title = expression(bold("B.") ~ "Latent date recovery"),
    x = "True date (CE)",
    y = "Inferred date (CE)"
  ) +
  theme_panel +
  theme(axis.title.x = element_text())

# ── Panel C: Posterior distributions of generating parameters ────────────────

true_baseline <- true_params$value[true_params$parameter == "baseline"]
true_slope    <- true_params$value[true_params$parameter == "slope"]
true_sigma    <- true_params$value[true_params$parameter == "sigma_noise"]

post_df <- tibble(
  Baseline = draws$baseline_original,
  Slope    = draws$slope_original,
  Sigma    = draws$sigma
) %>%
  pivot_longer(everything(), names_to = "Parameter", values_to = "Value") %>%
  mutate(Parameter = factor(Parameter, levels = c("Baseline", "Slope", "Sigma")))

true_lines <- tibble(
  Parameter = factor(c("Baseline", "Slope", "Sigma"),
                     levels = c("Baseline", "Slope", "Sigma")),
  True      = c(true_baseline, true_slope, true_sigma)
)

p_post <- ggplot(post_df, aes(x = Value)) +
  geom_histogram(fill = "grey70", colour = "black", linewidth = 0.3, bins = 40) +
  geom_vline(data = true_lines, aes(xintercept = True),
             linetype = "dashed", linewidth = 0.6, colour = "black") +
  facet_wrap(~ Parameter, scales = "free", nrow = 1) +
  labs(
    title = expression(bold("C.") ~ "Posterior distributions vs true values"),
    y = "Count"
  ) +
  theme_panel +
  theme(
    axis.title.x = element_text(),
    strip.background = element_blank(),
    strip.text = element_text(face = "bold", size = 10)
  )

# ── Combine ──────────────────────────────────────────────────────────────────

p <- (p_trend | p_dates) / p_post +
  plot_layout(heights = c(1, 0.7)) +
  plot_annotation(
    caption = "Dashed lines: true generating values. Shaded bands (A): 50% and 90% credible intervals.",
    theme = theme(plot.caption = element_text(hjust = 0, size = 8, colour = "grey40"))
  )

ggsave(here("Simulations", "Sim_Linear", "figures", "model_results_panel.png"), p,
       width = 12, height = 8, dpi = 300, bg = "white")

cat("Saved to", here("Simulations", "Sim_Linear", "figures", "model_results_panel.png"), "\n")
