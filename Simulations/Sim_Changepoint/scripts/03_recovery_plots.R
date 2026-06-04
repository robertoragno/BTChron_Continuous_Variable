# Figures for the recovery study (greyscale with dark-red/grey model coding):
#   recovery.png        - estimated vs true value for each parameter, both models
#   coverage.png        - how often each 50% / 90% interval contained the truth
#   sigma_vs_window.png - midpoint over-estimates the noise as dating gets vaguer
# Also writes recovery_table.csv (recovery rates) for the paper's summary table.

library(here)
library(readr)
library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)

# The changepoint posterior is harder to sample than the linear one: about a
# fifth of fits do not fully converge (max Rhat > 1.05 or divergences), spread
# across both models and tracking wider dating windows / higher noise rather than
# any one regime. Restrict to converged fits so the coverage rates are not
# distorted by unconverged chains; the rates are essentially unchanged whether or
# not this filter is applied, so it does not drive the conclusions.
all_results <- read_csv(
  here("Simulations", "Sim_Changepoint", "output", "recovery_results.csv"),
  show_col_types = FALSE
) %>% filter(is.na(error))

recovery_results <- all_results %>%
  filter(max_rhat <= 1.05, n_divergent == 0)

dropped <- all_results %>%
  group_by(model) %>%
  summarise(total = n(), kept = sum(max_rhat <= 1.05 & n_divergent == 0),
            .groups = "drop")
cat("Converged fits kept (max Rhat <= 1.05, no divergences):\n")
print(as.data.frame(dropped), row.names = FALSE)
cat("\n")

figure_path <- function(name) here("Simulations", "Sim_Changepoint", "figures", name)

label_model   <- function(code) factor(code, c("latent", "midpoint"),
                                        c("Latent-date", "Midpoint"))
model_colours <- c(`Latent-date` = "#780000", Midpoint = "grey25")
model_shapes  <- c(`Latent-date` = 16, Midpoint = 17)

param_levels <- c("Baseline", "Slope1", "Slope2", "CP", "Sigma")
param_labels <- c("Baseline", "Slope 1", "Slope 2", "Changepoint", "Sigma")

panel_theme <- theme_classic(base_size = 11) +
  theme(plot.title = element_text(size = 11, face = "bold"),
        plot.subtitle = element_text(size = 9),
        strip.background = element_blank(),
        strip.text = element_text(face = "bold"),
        legend.position = "top")

# recovery.png : estimated (posterior median) vs true value

estimates_vs_truth <- recovery_results %>%
  transmute(model = label_model(model),
            Baseline_estimate = baseline_med, Baseline_true = true_baseline,
            Slope1_estimate   = slope1_med,   Slope1_true   = true_slope_1,
            Slope2_estimate   = slope2_med,   Slope2_true   = true_slope_2,
            CP_estimate       = cp_med,       CP_true       = true_cp,
            Sigma_estimate    = sigma_med,    Sigma_true    = true_sigma) %>%
  pivot_longer(-model, names_to = c("parameter", ".value"), names_sep = "_") %>%
  mutate(parameter = factor(parameter, param_levels,
                            labels = param_labels))

recovery_plot <- ggplot(estimates_vs_truth,
                        aes(true, estimate, colour = model, shape = model)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey40") +
  geom_point(size = 0.7, alpha = 0.45) +
  facet_wrap(~ parameter, scales = "free", nrow = 1) +
  scale_colour_manual(values = model_colours, name = NULL) +
  scale_shape_manual(values = model_shapes, name = NULL) +
  guides(colour = guide_legend(override.aes = list(size = 2.5, alpha = 1))) +
  labs(x = "True value", y = "Estimated (posterior median)",
       title = "Recovery: estimated vs true value across random datasets",
       subtitle = "Points on the dashed 1:1 line = recovered") +
  panel_theme

ggsave(figure_path("recovery.png"), recovery_plot,
       width = 14, height = 4, dpi = 300, bg = "white")

# coverage.png : how often the interval contained the truth

coverage_rates <- recovery_results %>%
  group_by(model) %>%
  summarise(across(c(baseline_cov50, slope1_cov50, slope2_cov50, cp_cov50,
                     sigma_cov50, baseline_cov90, slope1_cov90, slope2_cov90,
                     cp_cov90, sigma_cov90), mean),
            .groups = "drop") %>%
  pivot_longer(-model, names_to = c("parameter", "interval"),
               names_sep = "_cov") %>%
  mutate(model = label_model(model),
         parameter = factor(parameter, tolower(param_levels),
                            labels = param_labels),
         interval = paste0(interval, "% interval"),
         target   = ifelse(grepl("50", interval), 0.5, 0.9))

coverage_plot <- ggplot(coverage_rates, aes(parameter, value, fill = model)) +
  geom_col(position = position_dodge(0.7), width = 0.65, colour = "black",
           linewidth = 0.2) +
  geom_hline(aes(yintercept = target), linetype = "dashed", colour = "grey40") +
  facet_wrap(~ interval) +
  scale_fill_manual(values = c("Latent-date" = "grey45", Midpoint = "grey82"),
                    name = NULL) +
  scale_y_continuous(labels = scales::percent, limits = c(0, 1)) +
  labs(x = NULL, y = "Datasets where the interval contained the truth",
       title = "Recovery rate (dashed line = target)") +
  panel_theme

ggsave(figure_path("coverage.png"), coverage_plot,
       width = 10, height = 4.5, dpi = 300, bg = "white")

# sigma_vs_window.png : midpoint mis-states the noise more as dating gets vaguer.
# The "where midpoint fails" panel for the changepoint simulation (matches the
# linear case): midpoint mistakes within-window date spread for measurement noise
# and over-estimates sigma, increasingly so as dating windows widen.

sigma_by_window <- recovery_results %>% mutate(model = label_model(model))

window_breaks  <- seq(0, 500, by = 100)
precision_bands <- tibble(
  band_start = window_breaks[-6], band_end = window_breaks[-1],
  band_centre = (window_breaks[-6] + window_breaks[-1]) / 2,
  label = c("Precise", "Fairly precise", "Moderate", "Vague", "Very vague"),
  shade = c("grey88", "grey91", "grey94", "grey97", "white"))

binned_sigma <- sigma_by_window %>%
  mutate(bin_centre = (floor(pmin(mean_width, 499) / 100) + 0.5) * 100) %>%
  group_by(model, bin_centre) %>%
  summarise(median_error = median(sigma_err),
            q25 = quantile(sigma_err, 0.25),
            q75 = quantile(sigma_err, 0.75), .groups = "drop") %>%
  mutate(x = bin_centre + ifelse(model == "Latent-date", -12, 12))

y_limits_sg <- quantile(sigma_by_window$sigma_err, c(0.004, 0.996))
label_y_sg  <- y_limits_sg[1] + 0.04 * diff(y_limits_sg)

sigma_window_plot <- ggplot() +
  geom_rect(data = precision_bands,
            aes(xmin = band_start, xmax = band_end, ymin = -Inf, ymax = Inf,
                fill = shade)) +
  scale_fill_identity() +
  geom_hline(yintercept = 0, linetype = "dotted", colour = "grey45") +
  geom_point(data = sigma_by_window,
             aes(mean_width, sigma_err, colour = model, shape = model),
             size = 0.55, alpha = 0.65) +
  geom_errorbar(data = binned_sigma,
                aes(x = x, ymin = q25, ymax = q75, colour = model),
                width = 16, linewidth = 0.7) +
  geom_point(data = binned_sigma, aes(x = x, y = median_error, colour = model,
                                      shape = model), size = 2) +
  geom_text(data = precision_bands, aes(x = band_centre, y = label_y_sg,
                                        label = label),
            colour = "grey30", size = 3) +
  scale_colour_manual(values = model_colours, name = NULL) +
  scale_shape_manual(values = model_shapes, name = NULL) +
  guides(colour = guide_legend(override.aes = list(size = 2.5, alpha = 1))) +
  coord_cartesian(xlim = c(20, 500), ylim = y_limits_sg, expand = FALSE) +
  scale_x_continuous(breaks = seq(100, 400, 100)) +
  labs(x = "Mean timespan (years)",
       y = expression(sigma[err] == sigma[est] - sigma[true]),
       title = "Midpoint over-estimates the noise as dating gets vaguer",
       caption = "Markers: bin median with interquartile range (25th-75th percentile).") +
  panel_theme

ggsave(figure_path("sigma_vs_window.png"), sigma_window_plot,
       width = 8.5, height = 5.5, dpi = 300, bg = "white")

# recovery_table.csv : recovery rates for the cross-simulation paper table

recovery_rate_table <- recovery_results %>%
  pivot_longer(matches("_cov(50|90)$"), names_to = c("parameter", "interval"),
               names_sep = "_cov", values_to = "contained_truth") %>%
  group_by(parameter, interval, model) %>%
  summarise(percent = round(100 * mean(contained_truth), 1), .groups = "drop") %>%
  pivot_wider(names_from = model, values_from = percent) %>%
  transmute(Simulation = "Changepoint",
            Parameter = factor(parameter,
                               c("baseline", "slope1", "slope2", "cp", "sigma"),
                               c("Baseline", "Slope 1", "Slope 2", "Changepoint", "Sigma")),
            Interval = paste0(interval, "%"),
            `Latent-date (%)` = latent, `Midpoint (%)` = midpoint) %>%
  arrange(Parameter, Interval)

write_csv(recovery_rate_table,
          here("Simulations", "Sim_Changepoint", "output", "recovery_table.csv"))

cat("Wrote recovery.png, coverage.png, sigma_vs_window.png, recovery_table.csv\n\n")
print(as.data.frame(recovery_rate_table), row.names = FALSE)
