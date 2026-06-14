# Recovery / calibration plots for the GP simulation.
#
# Reads the repeated-study output from 02_recovery_study.R and plots posterior
# rank histograms for sigma, together with the observed 50% and 90% coverage.

library(here)
library(tidyverse)
library(patchwork)

source(here("Simulations", "Sim_GP", "scripts", "simulate.R"))

input_dir <- Sys.getenv(
  "GP_RECOVERY_OUTPUT_DIR",
  here("Simulations", "Sim_GP", "output")
)
figure_dir <- Sys.getenv(
  "GP_RECOVERY_FIGURE_DIR",
  here("Simulations", "Sim_GP", "figures")
)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

results_latent <- read_csv(file.path(input_dir, "calibration_results.csv"),
                           show_col_types = FALSE)
results_midpoint <- read_csv(file.path(input_dir, "calibration_results_midpoint.csv"),
                             show_col_types = FALSE)
feature_results <- read_csv(file.path(input_dir, "feature_recovery_results.csv"),
                            show_col_types = FALSE)

label_model <- function(code) factor(code, c("latent", "midpoint"),
                                     c("Latent-date", "Midpoint"))
model_colours <- c(`Latent-date` = "#780000", Midpoint = "grey25")
model_shapes <- c(`Latent-date` = 16, Midpoint = 17)

target_labels <- c(
  peak_year = "Peak year",
  peak_value = "Peak value",
  fwhm_years = "Width (FWHM)",
  sigma_noise = "Sigma"
)

compute_coverage <- function(res) {
  tibble(
    CI = c("50%", "90%"),
    Coverage = c(
      mean(TRUE_SIGMA >= res$sigma_q25 & TRUE_SIGMA <= res$sigma_q75),
      mean(TRUE_SIGMA >= res$sigma_q05 & TRUE_SIGMA <= res$sigma_q95)
    )
  )
}

cov_latent <- compute_coverage(results_latent) %>% mutate(Model = "Latent-date")
cov_midpoint <- compute_coverage(results_midpoint) %>% mutate(Model = "Midpoint")
coverage <- bind_rows(cov_latent, cov_midpoint)

theme_panel <- theme_classic(base_size = 11) +
  theme(
    plot.title = element_text(size = 11, face = "bold", hjust = 0),
    plot.margin = margin(5, 10, 5, 10)
  )

n_draws <- 4000
n_bins <- 20
expected_count <- nrow(results_latent) / n_bins

rank_df <- bind_rows(
  tibble(Rank = results_latent$sigma_rank, Model = "Latent-date"),
  tibble(Rank = results_midpoint$sigma_rank, Model = "Midpoint")
) %>%
  mutate(Model = factor(Model, levels = c("Latent-date", "Midpoint")))

cov_labels <- coverage %>%
  mutate(label = sprintf("%s CI: %.0f%%", CI, Coverage * 100)) %>%
  group_by(Model) %>%
  summarise(text = paste(label, collapse = "\n"), .groups = "drop") %>%
  mutate(Model = factor(Model, levels = c("Latent-date", "Midpoint")))

p <- ggplot(rank_df, aes(x = Rank)) +
  geom_histogram(bins = n_bins, fill = "grey60", colour = "black", linewidth = 0.2) +
  geom_hline(yintercept = expected_count, linetype = "dashed", colour = "black") +
  geom_text(data = cov_labels, aes(label = text),
            x = n_draws * 0.95, y = Inf, hjust = 1, vjust = 1.5,
            size = 3, lineheight = 1.1) +
  facet_wrap(~ Model, nrow = 2) +
  labs(
    title = "Posterior rank histogram - sigma",
    subtitle = "Dashed line = expected bin count under uniform ranks",
    x = "Rank of true value among posterior draws",
    y = "Count"
  ) +
  theme_panel +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(face = "bold", size = 10)
  )

ggsave(file.path(figure_dir, "calibration_coverage.png"),
       p, width = 8, height = 7, dpi = 300, bg = "white")

feature_results <- feature_results %>%
  filter(is.na(error)) %>%
  mutate(
    model = label_model(model),
    target = factor(target, levels = names(target_labels), labels = unname(target_labels))
  )

feature_summary <- feature_results %>%
  group_by(model, target) %>%
  summarise(
    median_estimate = median(estimate_med),
    q25 = quantile(estimate_med, 0.25),
    q75 = quantile(estimate_med, 0.75),
    .groups = "drop"
  ) %>%
  # We nudge the summary markers sideways so the latent-date and midpoint model
  # summaries do not sit directly on top of one another inside each facet.
  mutate(x = ifelse(model == "Latent-date", 0.88, 2.12))

true_value_lines <- feature_results %>%
  distinct(target, true_value)

feature_recovery_plot <- ggplot(feature_results, aes(x = model, y = estimate_med, colour = model, shape = model)) +
  geom_hline(data = true_value_lines, aes(yintercept = true_value),
             inherit.aes = FALSE, linetype = "dashed", colour = "grey45") +
  geom_point(position = position_jitter(width = 0.09, height = 0), size = 0.7, alpha = 0.35) +
  geom_errorbar(
    data = feature_summary,
    aes(x = x, ymin = q25, ymax = q75, colour = model),
    inherit.aes = FALSE,
    width = 0.12,
    linewidth = 0.7
  ) +
  geom_point(
    data = feature_summary,
    aes(x = x, y = median_estimate, colour = model, shape = model),
    inherit.aes = FALSE,
    size = 2
  ) +
  facet_wrap(~ target, scales = "free_y", nrow = 1) +
  scale_colour_manual(values = model_colours, name = NULL) +
  scale_shape_manual(values = model_shapes, name = NULL) +
  guides(colour = guide_legend(override.aes = list(size = 2.5, alpha = 1))) +
  labs(
    x = NULL,
    y = "Estimated (posterior median)",
    title = "GP feature recovery under fixed date windows",
    subtitle = "Dashed line = true value; markers = median with interquartile range"
  ) +
  theme_panel +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(face = "bold", size = 10),
    legend.position = "top"
  )

ggsave(file.path(figure_dir, "feature_recovery.png"),
       feature_recovery_plot, width = 12, height = 4.6, dpi = 300, bg = "white")

z_95 <- qnorm(0.975)

# Wilson intervals behave better than a plain Wald interval when coverage is
# close to 0 or 1, which matters here because midpoint sigma coverage is
# effectively zero.
coverage_rates <- feature_results %>%
  pivot_longer(c(cov50, cov90), names_to = "interval_key", values_to = "contained_truth") %>%
  mutate(interval = ifelse(interval_key == "cov50", "50% interval", "90% interval")) %>%
  group_by(model, target, interval) %>%
  summarise(successes = sum(contained_truth), n = n(), value = mean(contained_truth), .groups = "drop") %>%
  mutate(
    denom = 1 + z_95^2 / n,
    centre = (value + z_95^2 / (2 * n)) / denom,
    half = z_95 * sqrt((value * (1 - value) + z_95^2 / (4 * n)) / n) / denom,
    lower = pmax(0, centre - half),
    upper = pmin(1, centre + half),
    target_rate = ifelse(interval == "50% interval", 0.5, 0.9)
  )

feature_coverage_plot <- ggplot(coverage_rates, aes(target, value, colour = model, shape = model)) +
  geom_hline(aes(yintercept = target_rate), linetype = "dashed", colour = "grey45") +
  geom_errorbar(aes(ymin = lower, ymax = upper),
                position = position_dodge(width = 0.4),
                width = 0, linewidth = 0.7) +
  geom_point(position = position_dodge(width = 0.4), size = 2.3) +
  facet_wrap(~ interval) +
  scale_colour_manual(values = model_colours, name = NULL) +
  scale_shape_manual(values = model_shapes, name = NULL) +
  guides(colour = guide_legend(override.aes = list(size = 2.5, alpha = 1))) +
  scale_y_continuous(labels = scales::percent, limits = c(0, 1)) +
  labs(
    x = NULL,
    y = "Fits where the interval contained the truth",
    title = "Empirical coverage of GP recovery targets",
    subtitle = "Error bars show 95% Wilson intervals"
  ) +
  theme_panel +
  theme(
    axis.text.x = element_text(angle = 15, hjust = 1),
    strip.background = element_blank(),
    strip.text = element_text(face = "bold", size = 10),
    legend.position = "top"
  )

ggsave(file.path(figure_dir, "feature_coverage.png"),
       feature_coverage_plot, width = 10.5, height = 4.8, dpi = 300, bg = "white")

feature_table <- coverage_rates %>%
  transmute(
    Target = as.character(target),
    Interval = ifelse(interval == "50% interval", "50%", "90%"),
    Model = as.character(model),
    Percent = round(100 * value, 1)
  ) %>%
  pivot_wider(names_from = Model, values_from = Percent)

write_csv(feature_table, file.path(input_dir, "feature_recovery_table.csv"))

cat("Wrote calibration_coverage.png, feature_recovery.png, feature_coverage.png, and feature_recovery_table.csv\n")
