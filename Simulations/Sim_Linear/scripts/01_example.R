# Worked example: ONE dataset drawn from the shared generating process
# (simulate.R), shown and fitted, purely to illustrate what a single dataset and
# fit look like. The systematic analysis is in 02_recovery_study.R, which draws
# many datasets from the same generator.
#
# Produces: exploratory_panel.png, model_comparison_single_fit.png,
#           individual_date_posteriors.png
# Fits are kept in memory (fast for one dataset), not cached to disk.

library(here)
library(tidyverse)
library(cmdstanr)
library(posterior)
library(patchwork)

source(here("Simulations", "Sim_Linear", "scripts", "simulate.R"))

figure_path <- function(name) here("Simulations", "Sim_Linear", "figures", name)

panel_theme <- theme_classic(base_size = 11) +
  theme(plot.title = element_text(size = 11, face = "bold", hjust = 0),
        plot.margin = margin(5, 10, 2, 10))

# 1. One illustrative dataset --------------------------------------------------

intercept <- 5
slope <- 0.015
sigma <- 1.5
n_phases <- 6        # K: number of periodisation phases
alpha_conc <- 1      # Dirichlet concentration (moderately even phases)
n_observations <- 300

sim_raw <- simulate_linear(N = n_observations, intercept = intercept,
                           slope = slope, sigma = sigma,
                           K = n_phases, alpha_conc = alpha_conc, seed = 42)
entropy_H <- attr(sim_raw, "H")

sim_data <- sim_raw %>%
  mutate(ID = row_number(),
         Site_name = sample(c("Site A", "Site B", "Site C", "Site D"),
                            n_observations, replace = TRUE)) %>%
  select(ID, Site_name, Start_date, End_date, True_date, Value)

ground_truth <- tibble(Year = seq(TMIN, TMAX),
                       True_value = intercept + slope * seq(TMIN, TMAX))

write_csv(sim_data, here("Simulations", "Sim_Linear", "data", "simulated_data.csv"))
write_csv(ground_truth, here("Simulations", "Sim_Linear", "data", "ground_truth.csv"))
write_csv(tibble(parameter = c("baseline", "slope", "sigma_noise",
                               "n_phases", "alpha_conc", "shannon_H"),
                 value = c(intercept, slope, sigma,
                           n_phases, alpha_conc, entropy_H)),
          here("Simulations", "Sim_Linear", "data", "generating_parameters.csv"))

# 2. Exploratory panel ---------------------------------------------------------

exploratory <- sim_data %>%
  mutate(Midpoint = (Start_date + End_date) / 2) %>%
  arrange(Midpoint) %>%
  mutate(row_position = row_number() * 1.8)

values_histogram <- ggplot(exploratory, aes(Value)) +
  geom_histogram(binwidth = 2, fill = "grey70", colour = "black", linewidth = 0.3) +
  labs(title = expression(bold("A.") ~ "Distribution of values"), x = "Value",
       y = "Count") + panel_theme

date_ranges_plot <- ggplot(exploratory, aes(y = row_position)) +
  geom_segment(aes(x = Start_date, xend = End_date, yend = row_position),
               linewidth = 0.2, colour = "grey40") +
  geom_point(aes(x = Midpoint), shape = 15, size = 0.35, colour = "black") +
  scale_x_continuous(breaks = seq(100, 900, 100), expand = c(0.01, 0)) +
  labs(title = expression(bold("B.") ~ "Sampled date ranges"), x = "Year CE",
       y = "Sample (ordered)") +
  panel_theme + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())

values_over_time <- ggplot(exploratory, aes(Midpoint, Value)) +
  geom_linerange(aes(xmin = Start_date, xmax = End_date), linewidth = 0.15,
                 colour = "grey60", alpha = 0.5) +
  geom_point(shape = 15, size = 0.6, colour = "black", alpha = 0.7) +
  geom_line(data = ground_truth, aes(Year, True_value), linewidth = 0.7,
            colour = "black", linetype = "dashed") +
  scale_x_continuous(breaks = seq(100, 900, 100), expand = c(0.01, 0)) +
  labs(title = expression(bold("C.") ~ "Values across time"), x = "Year CE",
       y = "Value") + panel_theme

ggsave(figure_path("exploratory_panel.png"),
       (values_histogram | date_ranges_plot | values_over_time) +
         plot_layout(widths = c(0.8, 1.8, 0.8)),
       width = 14, height = 4.5, dpi = 300, bg = "white")

# 3. Fit both models -----------------------------------------------------------

prediction_grid <- seq(min(sim_data$Start_date), max(sim_data$End_date), by = 1)
n_prediction_points <- length(prediction_grid)
stan_data <- list(N = n_observations, y = sim_data$Value,
                  start_date = sim_data$Start_date, end_date = sim_data$End_date,
                  N_pred = n_prediction_points, x_pred = prediction_grid)

fit_latent <- cmdstan_model(here("Simulations", "Sim_Linear", "models",
                                 "sim_linear.stan"))$sample(
  data = stan_data, chains = 4, parallel_chains = 4, iter_warmup = 1000,
  iter_sampling = 1000, seed = 42, adapt_delta = 0.95, max_treedepth = 12,
  refresh = 0, show_messages = FALSE)

fit_midpoint <- cmdstan_model(here("Simulations", "Sim_Linear", "models",
                                   "sim_linear_midpoint.stan"))$sample(
  data = stan_data, chains = 4, parallel_chains = 4, iter_warmup = 1000,
  iter_sampling = 1000, seed = 42, adapt_delta = 0.95, max_treedepth = 12,
  refresh = 0, show_messages = FALSE)

# 4. Model comparison (trend recovery + parameter posteriors) ------------------

summarise_trend <- function(fit) {
  trend_matrix <- as.matrix(
    fit$draws(format = "df")[, paste0("mu_pred[", 1:n_prediction_points, "]")])
  tibble(Year = prediction_grid,
         median   = apply(trend_matrix, 2, median),
         lower_90 = apply(trend_matrix, 2, quantile, 0.05),
         upper_90 = apply(trend_matrix, 2, quantile, 0.95),
         lower_50 = apply(trend_matrix, 2, quantile, 0.25),
         upper_50 = apply(trend_matrix, 2, quantile, 0.75))
}

# Periodisation phases, shared by both trend panels for the top strip.
# Alternating grey just separates neighbouring phases visually.
phase_bounds <- sim_data %>%
  distinct(Start_date, End_date) %>%
  arrange(Start_date) %>%
  mutate(shade = factor(row_number() %% 2))

phase_strip <- function(title_label) {
  ggplot(phase_bounds) +
    geom_rect(aes(xmin = Start_date, xmax = End_date, ymin = 0, ymax = 1,
                  fill = shade), colour = "white", linewidth = 0.3) +
    scale_fill_manual(values = c("0" = "grey85", "1" = "grey65"), guide = "none") +
    scale_x_continuous(limits = range(prediction_grid), expand = c(0.01, 0)) +
    scale_y_continuous(expand = c(0, 0)) +
    labs(title = title_label) +
    theme_void(base_size = 11) +
    theme(plot.title = element_text(size = 11, face = "bold", hjust = 0),
          plot.margin = margin(2, 10, 0, 10))
}

make_trend_panel <- function(fit, title_label) {
  trend <- ggplot(summarise_trend(fit)) +
    geom_ribbon(aes(Year, ymin = lower_90, ymax = upper_90), fill = "grey80") +
    geom_ribbon(aes(Year, ymin = lower_50, ymax = upper_50), fill = "grey60") +
    geom_line(aes(Year, median), linewidth = 0.6, colour = "black") +
    geom_line(data = ground_truth, aes(Year, True_value), linewidth = 0.7,
              colour = "black", linetype = "dashed") +
    scale_x_continuous(breaks = seq(100, 900, 200), expand = c(0.01, 0),
                       limits = range(prediction_grid)) +
    labs(x = "Year (CE)", y = "Value") + panel_theme
  phase_strip(title_label) / trend + plot_layout(heights = c(0.06, 1))
}

parameter_posteriors <- function(fit, model_label) {
  draws <- fit$draws(format = "df")
  tibble(Baseline = draws$baseline_original, Slope = draws$slope_original,
         Sigma = draws$sigma) %>%
    pivot_longer(everything(), names_to = "Parameter", values_to = "Value") %>%
    group_by(Parameter) %>%
    mutate(Region = case_when(
      Value >= quantile(Value, .25) & Value <= quantile(Value, .75) ~ "50% CI",
      Value >= quantile(Value, .05) & Value <= quantile(Value, .95) ~ "90% CI",
      TRUE ~ "Tail")) %>% ungroup() %>%
    mutate(Parameter = factor(Parameter, c("Baseline", "Slope", "Sigma")),
           Model = model_label)
}

both_posteriors <- bind_rows(
  parameter_posteriors(fit_latent, "EIV"),
  parameter_posteriors(fit_midpoint, "Midpoint dates")) %>%
  mutate(Model = factor(Model, c("EIV", "Midpoint dates")),
         Region = factor(Region, c("50% CI", "90% CI", "Tail")))

true_value_lines <- tibble(
  Parameter = factor(c("Baseline", "Slope", "Sigma"),
                     c("Baseline", "Slope", "Sigma")),
  True = c(intercept, slope, sigma))

posteriors_panel <- ggplot(both_posteriors, aes(Value, fill = Region)) +
  geom_histogram(colour = "black", linewidth = 0.2, bins = 40) +
  geom_vline(data = true_value_lines, aes(xintercept = True),
             linetype = "dashed", linewidth = 0.6) +
  scale_fill_manual(NULL, values = c("50% CI" = "grey40", "90% CI" = "grey65",
                                     "Tail" = "grey88")) +
  facet_grid(Model ~ Parameter, scales = "free_x") +
  labs(title = expression(bold("C.") ~ "Posterior distributions"), x = NULL,
       y = "Count") +
  panel_theme + theme(strip.background = element_blank(),
                      strip.text = element_text(face = "bold", size = 10),
                      legend.position = "bottom")

comparison_figure <-
  (make_trend_panel(fit_latent, expression(bold("A.") ~ "Trend — EIV")) |
   make_trend_panel(fit_midpoint, expression(bold("B.") ~ "Trend — midpoint dates"))) /
  posteriors_panel + plot_layout(heights = c(1, 1.2))

ggsave(figure_path("model_comparison_single_fit.png"), comparison_figure,
       width = 12, height = 10, dpi = 300, bg = "white")

# 5. Per-sample date posteriors (latent model) ---------------------------------

latent_draws <- fit_latent$draws(format = "df")
latent_trend_matrix <- as.matrix(
  latent_draws[, paste0("mu_pred[", 1:n_prediction_points, "]")])
latent_trend <- tibble(Year = prediction_grid,
                       median   = apply(latent_trend_matrix, 2, median),
                       lower_90 = apply(latent_trend_matrix, 2, quantile, 0.05),
                       upper_90 = apply(latent_trend_matrix, 2, quantile, 0.95))

make_date_panel <- function(sample_index) {
  sample_row   <- sim_data[sample_index, ]
  sample_label <- paste0(sample_row$ID, " (value = ", round(sample_row$Value, 2), ")")

  trend_with_sample <- ggplot(latent_trend) +
    geom_ribbon(aes(Year, ymin = lower_90, ymax = upper_90), fill = "grey80") +
    geom_line(aes(Year, median), linewidth = 0.6) +
    annotate("rect", xmin = sample_row$Start_date, xmax = sample_row$End_date,
             ymin = -Inf, ymax = Inf, fill = "grey50", alpha = 0.15) +
    geom_hline(yintercept = sample_row$Value, linetype = "dashed",
               linewidth = 0.5) +
    labs(title = paste0(sample_label, " — trend vs observed value"),
         x = "Date (CE)", y = "Value") + panel_theme

  estimated_dates <- latent_draws[[paste0("true_date_actual[", sample_index, "]")]]
  date_density <- ggplot(tibble(estimated_date = estimated_dates),
                         aes(estimated_date)) +
    geom_density(fill = "grey60", colour = "black", linewidth = 0.4, alpha = 0.5) +
    geom_vline(xintercept = c(sample_row$Start_date, sample_row$End_date),
               linetype = "dashed", colour = "grey40", linewidth = 0.5) +
    geom_vline(xintercept = sample_row$True_date, linewidth = 0.6) +
    labs(title = paste0(sample_label, " — posterior date estimate"),
         x = "Estimated date (CE)", y = "Density") + panel_theme

  trend_with_sample | date_density
}

set.seed(123)
date_panels <- map(sort(sample(seq_len(n_observations), 6)), make_date_panel)
ggsave(figure_path("individual_date_posteriors.png"),
       wrap_plots(date_panels, ncol = 1),
       width = 12, height = 18, dpi = 300, bg = "white")

