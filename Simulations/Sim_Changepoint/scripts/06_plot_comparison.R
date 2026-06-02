#' Purpose: Side-by-side comparison of latent-date vs midpoint changepoint model.
#'   Row 1 — Trend recovery (A: latent dates, B: midpoint dates)
#'   Row 2 — Posterior parameter distributions for both models
#'
# caption = "A–B: posterior median (solid) and 50%/90% credible intervals vs
#   true generating trend (dashed). C: posterior histograms with graduated
#   shading by credible interval; dashed lines mark the true generating values."

library(here)
library(tidyverse)
library(cmdstanr)
library(posterior)
library(patchwork)

# Load data

sim_data     <- read_csv(here("Simulations", "Sim_Changepoint", "data",
                              "simulated_data.csv"),
                         show_col_types = FALSE)
ground_truth <- read_csv(here("Simulations", "Sim_Changepoint", "data",
                              "ground_truth.csv"),
                         show_col_types = FALSE)
true_params  <- read_csv(here("Simulations", "Sim_Changepoint", "data",
                              "generating_parameters.csv"),
                         show_col_types = FALSE)

fit_latent   <- readRDS(here("Simulations", "Sim_Changepoint", "output",
                             "fit_changepoint.rds"))
fit_midpoint <- readRDS(here("Simulations", "Sim_Changepoint", "output",
                             "fit_midpoint.rds"))

true_vals <- setNames(true_params$value, true_params$parameter)

# Shared theme

theme_panel <- theme_classic(base_size = 11) +
  theme(
    plot.title   = element_text(size = 11, face = "bold", hjust = 0),
    plot.margin  = margin(5, 10, 2, 10)
  )

pred_grid <- seq(min(sim_data$Start_date), max(sim_data$End_date), by = 1)
N_pred    <- length(pred_grid)

# Helper: extract trend summary

extract_trend <- function(fit) {
  draws   <- fit$draws(format = "df")
  mu_mat  <- as.matrix(draws[, paste0("mu_pred[", 1:N_pred, "]")])
  tibble(
    Year     = pred_grid,
    Median   = apply(mu_mat, 2, median),
    Lower_90 = apply(mu_mat, 2, quantile, 0.05),
    Upper_90 = apply(mu_mat, 2, quantile, 0.95),
    Lower_50 = apply(mu_mat, 2, quantile, 0.25),
    Upper_50 = apply(mu_mat, 2, quantile, 0.75)
  )
}

trend_latent   <- extract_trend(fit_latent)
trend_midpoint <- extract_trend(fit_midpoint)

# Helper: extract parameter posteriors with CI regions

extract_params <- function(fit) {
  draws <- fit$draws(format = "df")
  post_df <- tibble(
    Baseline    = draws$baseline_original,
    `Slope 1`   = draws$slope1_original,
    `Slope 2`   = draws$slope2_original,
    Changepoint = draws$cp_actual,
    Sigma       = draws$sigma
  ) %>%
    pivot_longer(everything(), names_to = "Parameter", values_to = "Value") %>%
    mutate(Parameter = factor(
      Parameter,
      levels = c("Baseline", "Slope 1", "Slope 2", "Changepoint", "Sigma")
    ))

  ci_df <- post_df %>%
    group_by(Parameter) %>%
    summarise(
      q05 = quantile(Value, 0.05),
      q25 = quantile(Value, 0.25),
      q75 = quantile(Value, 0.75),
      q95 = quantile(Value, 0.95),
      .groups = "drop"
    )

  post_df %>%
    left_join(ci_df, by = "Parameter") %>%
    mutate(
      Region = case_when(
        Value >= q25 & Value <= q75 ~ "50% CI",
        Value >= q05 & Value <= q95 ~ "90% CI",
        TRUE                        ~ "Tail"
      ),
      Region = factor(Region, levels = c("50% CI", "90% CI", "Tail"))
    )
}

params_latent   <- extract_params(fit_latent)   %>%
  mutate(Model = "Latent dates")
params_midpoint <- extract_params(fit_midpoint) %>%
  mutate(Model = "Midpoint dates")

params_both <- bind_rows(params_latent, params_midpoint) %>%
  mutate(Model = factor(Model, levels = c("Latent dates", "Midpoint dates")))

# True value reference lines

true_lines <- tibble(
  Parameter = factor(
    c("Baseline", "Slope 1", "Slope 2", "Changepoint", "Sigma"),
    levels = c("Baseline", "Slope 1", "Slope 2", "Changepoint", "Sigma")
  ),
  True = c(
    true_vals["baseline"],
    true_vals["slope_1"],
    true_vals["slope_2"],
    true_vals["changepoint"],
    true_vals["sigma_noise"]
  )
)

# Row 1: Trend recovery

make_trend_panel <- function(df, title_label) {
  ggplot(df) +
    geom_ribbon(aes(x = Year, ymin = Lower_90, ymax = Upper_90),
                fill = "grey80") +
    geom_ribbon(aes(x = Year, ymin = Lower_50, ymax = Upper_50),
                fill = "grey60") +
    geom_line(aes(x = Year, y = Median), linewidth = 0.6, colour = "black") +
    geom_line(data = ground_truth, aes(x = Year, y = True_value),
              linewidth = 0.7, colour = "black", linetype = "dashed") +
    scale_x_continuous(breaks = seq(100, 900, 200), expand = c(0.01, 0)) +
    labs(title = title_label, x = "Year CE", y = "Value") +
    theme_panel
}

p_trend_latent <- make_trend_panel(
  trend_latent,
  expression(bold("A.") ~ "Trend recovery — latent dates")
)

p_trend_mid <- make_trend_panel(
  trend_midpoint,
  expression(bold("B.") ~ "Trend recovery — midpoint dates")
)

# Row 2: Parameter posteriors

p_post <- ggplot(params_both, aes(x = Value, fill = Region)) +
  geom_histogram(colour = "black", linewidth = 0.2, bins = 40) +
  geom_vline(data = true_lines, aes(xintercept = True),
             linetype = "dashed", linewidth = 0.6, colour = "black") +
  scale_fill_manual(
    name   = NULL,
    values = c("50% CI" = "grey40", "90% CI" = "grey65", "Tail" = "grey88")
  ) +
  facet_grid(Model ~ Parameter, scales = "free_x") +
  labs(
    title = expression(bold("C.") ~ "Posterior distributions"),
    x = NULL, y = "Count"
  ) +
  theme_panel +
  theme(
    strip.background = element_blank(),
    strip.text       = element_text(face = "bold", size = 10),
    legend.position  = "bottom",
    legend.text      = element_text(size = 8),
    legend.key.size  = unit(10, "pt")
  )

# Combine

p <- (p_trend_latent | p_trend_mid) /
  p_post +
  plot_layout(heights = c(1, 1.2)) +
  plot_annotation(theme = theme(plot.margin = margin(5, 5, 5, 5)))

ggsave(
  here("Simulations", "Sim_Changepoint", "figures", "model_comparison.png"),
  p, width = 14, height = 10, dpi = 300, bg = "white"
)

cat("Saved to",
    here("Simulations", "Sim_Changepoint", "figures",
         "model_comparison.png"), "\n")
