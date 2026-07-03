# Worked example: ONE dataset drawn from the shared generating process
# (simulate.R), shown and fitted, to illustrate what a single GP dataset and fit
# look like. The systematic analysis is in 02_recovery_study.R, which draws many
# datasets from the same generator.
#
# Produces: exploratory_panel.png, model_comparison_single_fit.png,
#           individual_date_posteriors.png
# Fits are kept in memory (fast enough for one dataset), not cached to disk.

library(here)
library(tidyverse)
library(cmdstanr)
library(posterior)
library(patchwork)

source(here("Simulations", "Sim_GP", "scripts", "simulate.R"))

figure_path <- function(name) here("Simulations", "Sim_GP", "figures", name)

panel_theme <- theme_classic(base_size = 11) +
  theme(plot.title = element_text(size = 11, face = "bold", hjust = 0),
        plot.margin = margin(5, 10, 2, 10))

# HSGP settings. M basis functions and boundary factor c; kept as data so the
# same model file serves the example and the study. M around 25 is enough to
# represent curves on the wiggle scale here; more modes only make the latent-date
# geometry harder to sample, not the fit better.
M_BASIS <- 20
C_BOUNDARY <- 1.5

# 1. One illustrative dataset --------------------------------------------------

# A fixed gently undulating curve so the figure is stable between runs. Two waves
# with periods longer than a typical dating window.
curve <- list(
  baseline = 8,
  amp      = c(6, 4),
  period   = c(520, 380),
  phase    = c(0, 1.5)
)
sigma <- 2.5
n_phases <- 6         # K: number of periodisation phases
alpha_conc <- 1       # Dirichlet concentration (moderately even phases)
n_observations <- 250

sim_raw <- simulate_gp(N = n_observations, sigma = sigma, K = n_phases,
                       alpha_conc = alpha_conc, curve = curve, seed = 42)
entropy_H <- attr(sim_raw, "H")

sim_data <- sim_raw %>%
  mutate(ID = row_number(),
         Site_name = sample(c("Site A", "Site B", "Site C", "Site D"),
                            n_observations, replace = TRUE)) %>%
  select(ID, Site_name, Start_date, End_date, True_date, Value)

ground_truth <- truth_grid(curve)

write_csv(sim_data, here("Simulations", "Sim_GP", "data", "simulated_data.csv"))
write_csv(ground_truth, here("Simulations", "Sim_GP", "data", "ground_truth.csv"))
write_csv(tibble(parameter = c("baseline", "sigma_noise", "n_components",
                               "n_phases", "alpha_conc", "shannon_H"),
                 value = c(curve$baseline, sigma, N_COMPONENTS,
                           n_phases, alpha_conc, entropy_H)),
          here("Simulations", "Sim_GP", "data", "generating_parameters.csv"))

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
                  N_pred = n_prediction_points, x_pred = prediction_grid,
                  M = M_BASIS, c = C_BOUNDARY)

fit_latent <- cmdstan_model(here("Simulations", "Sim_GP", "models",
                                 "sim_hsgp.stan"))$sample(
  data = stan_data, chains = 4, parallel_chains = 4, iter_warmup = 1000,
  iter_sampling = 1000, seed = 42, adapt_delta = 0.99, max_treedepth = 12,
  refresh = 0, show_messages = FALSE)

fit_midpoint <- cmdstan_model(here("Simulations", "Sim_GP", "models",
                                   "sim_hsgp_midpoint.stan"))$sample(
  data = stan_data, chains = 4, parallel_chains = 4, iter_warmup = 1000,
  iter_sampling = 1000, seed = 42, adapt_delta = 0.99, max_treedepth = 12,
  refresh = 0, show_messages = FALSE)

# 4. Model comparison (curve recovery + parameter posteriors) ------------------

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

# Proportion of the true curve caught by the 90% band, over the prediction grid.
# This is the accuracy measure the recovery study uses for the whole curve.
curve_coverage <- function(fit) {
  trend <- summarise_trend(fit)
  truth <- f_true(trend$Year, curve)
  mean(truth >= trend$lower_90 & truth <= trend$upper_90)
}

make_trend_panel <- function(fit, title_label) {
  covered <- curve_coverage(fit)
  ggplot(summarise_trend(fit)) +
    geom_ribbon(aes(Year, ymin = lower_90, ymax = upper_90), fill = "grey80") +
    geom_ribbon(aes(Year, ymin = lower_50, ymax = upper_50), fill = "grey60") +
    geom_line(aes(Year, median), linewidth = 0.6, colour = "black") +
    geom_line(data = ground_truth, aes(Year, True_value), linewidth = 0.7,
              colour = "black", linetype = "dashed") +
    scale_x_continuous(breaks = seq(100, 900, 200), expand = c(0.01, 0)) +
    labs(title = title_label, x = "Year (CE)", y = "Value",
         subtitle = sprintf("Proportion of the true curve inside the 90%% CI: %.0f%%",
                            100 * covered)) +
    panel_theme + theme(plot.subtitle = element_text(size = 9))
}

# sigma has a known true value; mu and rho are properties of the fitted curve
# (no simple generator truth), shown for context.
parameter_posteriors <- function(fit, model_label) {
  draws <- fit$draws(format = "df")
  tibble(Sigma = draws$sigma, Mu = draws$mu, Rho = draws$rho_actual) %>%
    pivot_longer(everything(), names_to = "Parameter", values_to = "Value") %>%
    group_by(Parameter) %>%
    mutate(Region = case_when(
      Value >= quantile(Value, .25) & Value <= quantile(Value, .75) ~ "50% CI",
      Value >= quantile(Value, .05) & Value <= quantile(Value, .95) ~ "90% CI",
      TRUE ~ "Tail")) %>% ungroup() %>%
    mutate(Parameter = factor(Parameter, c("Sigma", "Mu", "Rho")),
           Model = model_label)
}

both_posteriors <- bind_rows(
  parameter_posteriors(fit_latent, "EIV"),
  parameter_posteriors(fit_midpoint, "Midpoint dates")) %>%
  mutate(Model = factor(Model, c("EIV", "Midpoint dates")),
         Region = factor(Region, c("50% CI", "90% CI", "Tail")))

# Only sigma has a generator truth to mark.
true_value_lines <- tibble(Parameter = factor("Sigma", c("Sigma", "Mu", "Rho")),
                           True = sigma)

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
  (make_trend_panel(fit_latent, expression(bold("A.") ~ "Curve — EIV")) |
   make_trend_panel(fit_midpoint, expression(bold("B.") ~ "Curve — midpoint dates"))) /
  posteriors_panel + plot_layout(heights = c(1, 1.2))

ggsave(figure_path("model_comparison_single_fit.png"), comparison_figure,
       width = 12, height = 10, dpi = 300, bg = "white")

# 5. Per-sample date posteriors (EIV model) ------------------------------------

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
    labs(title = paste0(sample_label, " — curve vs observed value"),
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
