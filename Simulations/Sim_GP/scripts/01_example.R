# Worked example for the GP simulation.
#
# This mirrors the 01_example.R scripts in the linear and changepoint folders:
# write one dataset, fit both models, and produce the single-dataset figures.

library(here)
library(tidyverse)
library(cmdstanr)
library(posterior)
library(patchwork)

source(here("Simulations", "Sim_GP", "scripts", "simulate.R"))

sim_data <- simulate_gp_dataset()
ground_truth <- ground_truth_grid()
true_params <- generating_parameters()
write_example_files(sim_data)

glimpse(sim_data)
summary(sim_data$Value)
cat("Wrote simulated_data.csv, ground_truth.csv, generating_parameters.csv, and evaluation_targets.csv\n")

# Exploratory panel

sim_data_expl <- sim_data %>%
  mutate(Midpoint = (Start_date + End_date) / 2) %>%
  arrange(Midpoint) %>%
  mutate(rank = row_number(), rank_spaced = rank * 1.8)

theme_panel <- theme_classic(base_size = 11) +
  theme(
    plot.title = element_text(size = 11, face = "bold", hjust = 0),
    plot.margin = margin(5, 10, 2, 10),
    axis.title.x = element_blank()
  )

p_hist <- ggplot(sim_data_expl, aes(x = Value)) +
  geom_histogram(binwidth = 2, fill = "grey70", colour = "black", linewidth = 0.3) +
  labs(
    title = expression(bold("A.") ~ "Distribution of values"),
    y = "Count"
  ) +
  theme_panel +
  theme(axis.title.x = element_text())

tick_h <- 0.5

p_dur <- ggplot(sim_data_expl, aes(y = rank_spaced)) +
  geom_segment(aes(x = Start_date, xend = End_date, yend = rank_spaced),
               linewidth = 0.2, colour = "grey40") +
  geom_segment(aes(x = Start_date, xend = Start_date,
                   y = rank_spaced - tick_h, yend = rank_spaced + tick_h),
               linewidth = 0.3, colour = "grey30") +
  geom_segment(aes(x = End_date, xend = End_date,
                   y = rank_spaced - tick_h, yend = rank_spaced + tick_h),
               linewidth = 0.3, colour = "grey30") +
  geom_point(aes(x = Midpoint), shape = 15, size = 0.35, colour = "black") +
  scale_x_continuous(breaks = seq(100, 900, 100), expand = c(0.01, 0)) +
  scale_y_continuous(expand = expansion(mult = 0.02)) +
  labs(
    title = expression(bold("B.") ~ "Sampled date ranges"),
    x = "Year CE",
    y = "Sample (ordered)"
  ) +
  theme_panel +
  theme(
    axis.title.x = element_text(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank()
  )

p_val <- ggplot(sim_data_expl, aes(x = Midpoint, y = Value)) +
  geom_linerange(aes(xmin = Start_date, xmax = End_date),
                 linewidth = 0.15, colour = "grey60", alpha = 0.5) +
  geom_point(shape = 15, size = 0.6, colour = "black", alpha = 0.7) +
  geom_line(data = ground_truth, aes(x = Year, y = True_value),
            linewidth = 0.7, colour = "black", linetype = "dashed") +
  scale_x_continuous(breaks = seq(100, 900, 100), expand = c(0.01, 0)) +
  labs(
    title = expression(bold("C.") ~ "Distribution of values across time"),
    x = "Year CE",
    y = "Value"
  ) +
  theme_panel +
  theme(axis.title.x = element_text())

p_exploratory <- (p_hist | p_dur | p_val) +
  plot_layout(widths = c(0.8, 1.8, 0.8)) +
  plot_annotation(
    caption = "Dashed line (C): true Gaussian bell-curve trend used to generate the simulated values.",
    theme = theme(plot.caption = element_text(hjust = 0, size = 8, colour = "grey40"))
  )

ggsave(here("Simulations", "Sim_GP", "figures", "exploratory_panel.png"),
       p_exploratory, width = 14, height = 4.5, dpi = 300, bg = "white")
cat("Saved to", here("Simulations", "Sim_GP", "figures", "exploratory_panel.png"), "\n")

# Fit the latent-date model

pred_grid <- seq(min(sim_data$Start_date), max(sim_data$End_date), by = 1)
stan_data <- list(
  N = nrow(sim_data),
  y = sim_data$Value,
  start_date = sim_data$Start_date,
  end_date = sim_data$End_date,
  N_pred = length(pred_grid),
  x_pred = pred_grid,
  M = 20,
  c = 1.5
)

model_latent <- cmdstan_model(here("Simulations", "Sim_GP", "models", "sim_hsgp.stan"))
latent_fit_path <- here("Simulations", "Sim_GP", "output", "fit_hsgp.rds")

if (file.exists(latent_fit_path)) {
  cat("Loading existing latent-date fit.\n")
  fit_latent <- readRDS(latent_fit_path)
  ran_latent_fit <- FALSE
} else {
  fit_latent <- model_latent$sample(
    data = stan_data,
    chains = 4,
    parallel_chains = 4,
    iter_warmup = 1000,
    iter_sampling = 1000,
    seed = 42,
    adapt_delta = 0.99,
    max_treedepth = 12
  )
  fit_latent$save_object(latent_fit_path)
  ran_latent_fit <- TRUE
}

if (ran_latent_fit) fit_latent$cmdstan_diagnose()
fit_latent$summary(variables = c("mu", "alpha", "rho", "sigma", "rho_actual"))

cat("\n-- Generating parameters --\n")
print(true_params)

# Latent-date model results

draws_latent <- fit_latent$draws(format = "df")
N_pred <- length(pred_grid)
trend_cols <- paste0("mu_pred[", seq_len(N_pred), "]")
trend_mat <- as.matrix(draws_latent[, trend_cols])

trend_summary <- tibble(
  Year = pred_grid,
  Median = apply(trend_mat, 2, median),
  Lower_90 = apply(trend_mat, 2, quantile, 0.05),
  Upper_90 = apply(trend_mat, 2, quantile, 0.95),
  Lower_50 = apply(trend_mat, 2, quantile, 0.25),
  Upper_50 = apply(trend_mat, 2, quantile, 0.75)
)

p_trend <- ggplot(trend_summary) +
  geom_ribbon(aes(x = Year, ymin = Lower_90, ymax = Upper_90, fill = "90% CI")) +
  geom_ribbon(aes(x = Year, ymin = Lower_50, ymax = Upper_50, fill = "50% CI")) +
  geom_line(aes(x = Year, y = Median), linewidth = 0.6, colour = "black") +
  geom_line(data = ground_truth, aes(x = Year, y = True_value),
            linewidth = 0.7, colour = "black", linetype = "dashed") +
  geom_rug(data = sim_data, aes(y = Value),
           sides = "r", colour = "grey50", alpha = 0.3, length = grid::unit(3, "pt")) +
  scale_fill_manual(
    name = NULL,
    values = c("90% CI" = "grey80", "50% CI" = "grey60"),
    guide = guide_legend(override.aes = list(alpha = 1))
  ) +
  scale_x_continuous(breaks = seq(100, 900, 100), expand = c(0.01, 0)) +
  labs(
    title = expression(bold("A.") ~ "Recovered trend vs true trend"),
    x = "Year CE",
    y = "Value"
  ) +
  theme_panel +
  theme(
    axis.title.x = element_text(),
    legend.position = c(0.05, 0.95),
    legend.justification = c(0, 1),
    legend.background = element_rect(fill = alpha("white", 0.8), colour = NA),
    legend.key.size = grid::unit(10, "pt"),
    legend.text = element_text(size = 8)
  )

date_cols <- paste0("true_date_actual[", seq_len(nrow(sim_data)), "]")
date_mat <- as.matrix(draws_latent[, date_cols])

date_recovery <- tibble(
  True_date = sim_data$True_date,
  Inferred_med = apply(date_mat, 2, median),
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
    subtitle = "Each point = one sample. Dashed line = perfect 1:1 recovery.",
    x = "True generating date (CE)",
    y = "Posterior median date (CE)"
  ) +
  theme_panel +
  theme(
    axis.title.x = element_text(),
    plot.subtitle = element_text(size = 9, colour = "grey40", hjust = 0)
  )

post_df <- tibble(
  Sigma = draws_latent$sigma,
  Mu = draws_latent$mu,
  Rho_actual = draws_latent$rho_actual
) %>%
  pivot_longer(everything(), names_to = "Parameter", values_to = "Value") %>%
  mutate(Parameter = factor(Parameter, levels = c("Sigma", "Mu", "Rho_actual"),
                            labels = c("sigma", "mu", "rho (CE)")))

ci_df <- post_df %>%
  group_by(Parameter) %>%
  summarise(
    q05 = quantile(Value, 0.05),
    q25 = quantile(Value, 0.25),
    q75 = quantile(Value, 0.75),
    q95 = quantile(Value, 0.95),
    .groups = "drop"
  )

post_df <- post_df %>%
  left_join(ci_df, by = "Parameter") %>%
  mutate(
    Region = case_when(
      Value >= q25 & Value <= q75 ~ "50% CI",
      Value >= q05 & Value <= q95 ~ "90% CI",
      TRUE ~ "Tail"
    ),
    Region = factor(Region, levels = c("50% CI", "90% CI", "Tail"))
  )

true_sigma_line <- tibble(
  Parameter = factor("sigma", levels = c("sigma", "mu", "rho (CE)")),
  True = TRUE_SIGMA
)

p_post <- ggplot(post_df, aes(x = Value, fill = Region)) +
  geom_histogram(colour = "black", linewidth = 0.2, bins = 40) +
  geom_vline(data = true_sigma_line, aes(xintercept = True),
             linetype = "dashed", linewidth = 0.6, colour = "black") +
  scale_fill_manual(
    name = NULL,
    values = c("50% CI" = "grey40", "90% CI" = "grey65", "Tail" = "grey88")
  ) +
  facet_wrap(~ Parameter, scales = "free", nrow = 1) +
  labs(
    title = expression(bold("C.") ~ "Posterior distributions"),
    caption = "Dashed line (sigma only): true generating value.",
    y = "Count"
  ) +
  theme_panel +
  theme(
    axis.title.x = element_text(),
    strip.background = element_blank(),
    strip.text = element_text(face = "bold", size = 10),
    legend.position = "bottom",
    legend.text = element_text(size = 8),
    legend.key.size = grid::unit(10, "pt"),
    plot.caption = element_text(hjust = 0, size = 8, colour = "grey40")
  )

p_model_results <- (p_trend | p_dates) / p_post +
  plot_layout(heights = c(1, 0.7)) +
  plot_annotation(theme = theme(plot.margin = margin(5, 5, 5, 5)))

ggsave(here("Simulations", "Sim_GP", "figures", "model_results_panel.png"),
       p_model_results, width = 12, height = 8, dpi = 300, bg = "white")
cat("Saved to", here("Simulations", "Sim_GP", "figures", "model_results_panel.png"), "\n")

# Fit the midpoint model

model_midpoint <- cmdstan_model(here("Simulations", "Sim_GP", "models", "sim_hsgp_midpoint.stan"))
midpoint_fit_path <- here("Simulations", "Sim_GP", "output", "fit_hsgp_midpoint.rds")

if (file.exists(midpoint_fit_path)) {
  cat("Loading existing midpoint fit.\n")
  fit_midpoint <- readRDS(midpoint_fit_path)
  ran_midpoint_fit <- FALSE
} else {
  fit_midpoint <- model_midpoint$sample(
    data = stan_data,
    chains = 4,
    parallel_chains = 4,
    iter_warmup = 1000,
    iter_sampling = 1000,
    seed = 42,
    adapt_delta = 0.99,
    max_treedepth = 12
  )
  fit_midpoint$save_object(midpoint_fit_path)
  ran_midpoint_fit <- TRUE
}

if (ran_midpoint_fit) fit_midpoint$cmdstan_diagnose()
fit_midpoint$summary(variables = c("mu", "alpha", "rho", "sigma", "rho_actual"))

# Model comparison

extract_trend <- function(fit) {
  draws <- fit$draws(format = "df")
  trend_mat <- as.matrix(draws[, paste0("mu_pred[", seq_len(N_pred), "]")])
  tibble(
    Year = pred_grid,
    Median = apply(trend_mat, 2, median),
    Lower_90 = apply(trend_mat, 2, quantile, 0.05),
    Upper_90 = apply(trend_mat, 2, quantile, 0.95),
    Lower_50 = apply(trend_mat, 2, quantile, 0.25),
    Upper_50 = apply(trend_mat, 2, quantile, 0.75)
  )
}

trend_latent <- extract_trend(fit_latent)
trend_midpoint <- extract_trend(fit_midpoint)

make_trend_panel <- function(df, title_label) {
  ggplot(df) +
    geom_ribbon(aes(x = Year, ymin = Lower_90, ymax = Upper_90), fill = "grey80") +
    geom_ribbon(aes(x = Year, ymin = Lower_50, ymax = Upper_50), fill = "grey60") +
    geom_line(aes(x = Year, y = Median), linewidth = 0.6, colour = "black") +
    geom_line(data = ground_truth, aes(x = Year, y = True_value),
              linewidth = 0.7, colour = "black", linetype = "dashed") +
    scale_x_continuous(breaks = seq(100, 900, 200), expand = c(0.01, 0)) +
    labs(title = title_label, x = "Year CE", y = "Value") +
    theme_panel
}

p_trend_latent <- make_trend_panel(
  trend_latent,
  expression(bold("A.") ~ "Trend recovery - latent dates")
)
p_trend_mid <- make_trend_panel(
  trend_midpoint,
  expression(bold("B.") ~ "Trend recovery - midpoint dates")
)

extract_sigma <- function(fit, model_label) {
  draws <- fit$draws(format = "df")
  sigma <- draws$sigma
  q05 <- quantile(sigma, 0.05)
  q25 <- quantile(sigma, 0.25)
  q75 <- quantile(sigma, 0.75)
  q95 <- quantile(sigma, 0.95)

  tibble(sigma = sigma, Model = model_label) %>%
    mutate(
      Region = case_when(
        sigma >= q25 & sigma <= q75 ~ "50% CI",
        sigma >= q05 & sigma <= q95 ~ "90% CI",
        TRUE ~ "Tail"
      ),
      Region = factor(Region, levels = c("50% CI", "90% CI", "Tail"))
    )
}

sigma_df <- bind_rows(
  extract_sigma(fit_latent, "Latent dates"),
  extract_sigma(fit_midpoint, "Midpoint dates")
) %>%
  mutate(Model = factor(Model, levels = c("Latent dates", "Midpoint dates")))

p_sigma <- ggplot(sigma_df, aes(x = sigma, fill = Region)) +
  geom_histogram(colour = "black", linewidth = 0.2, bins = 40) +
  geom_vline(xintercept = TRUE_SIGMA,
             linetype = "dashed", linewidth = 0.6, colour = "black") +
  scale_fill_manual(
    name = NULL,
    values = c("50% CI" = "grey40", "90% CI" = "grey65", "Tail" = "grey88")
  ) +
  facet_wrap(~ Model, scales = "free_y", nrow = 1) +
  labs(
    title = expression(bold("C.") ~ "Posterior sigma - latent vs midpoint"),
    caption = "Dashed line: true generating sigma (2.5).",
    x = "sigma",
    y = "Count"
  ) +
  theme_panel +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(face = "bold", size = 10),
    legend.position = "bottom",
    legend.text = element_text(size = 8),
    legend.key.size = grid::unit(10, "pt"),
    plot.caption = element_text(hjust = 0, size = 8, colour = "grey40")
  )

p_comparison <- (p_trend_latent | p_trend_mid) / p_sigma +
  plot_layout(heights = c(1, 0.8)) +
  plot_annotation(theme = theme(plot.margin = margin(5, 5, 5, 5)))

ggsave(here("Simulations", "Sim_GP", "figures", "model_comparison.png"),
       p_comparison, width = 12, height = 9, dpi = 300, bg = "white")
cat("Saved to", here("Simulations", "Sim_GP", "figures", "model_comparison.png"), "\n")
