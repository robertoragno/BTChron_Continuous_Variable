# Compare the latent-date and midpoint fits on the GB subset, in the same style
# as the simulation. No ground-truth overlays: real data has no known trend, so
# the contrast between the two models is the point. Run after 02_fit.R.
#
# Produces: dataset_1_model_comparison.png, dataset_1_individual_date_posteriors.png

library(here)
library(tidyverse)
library(posterior)
library(patchwork)

clean_path <- here("Real_Data", "dataset_1", "data", "gini_gb_filtered.csv")
output_dir <- here("Real_Data", "dataset_1", "output")
figure_path <- function(name) here("Real_Data", "dataset_1", "figures", name)

panel_theme <- theme_classic(base_size = 11) +
  theme(plot.title = element_text(size = 11, face = "bold", hjust = 0),
        plot.margin = margin(5, 10, 2, 10))

gb_data <- read_csv(clean_path, show_col_types = FALSE)
fit_latent   <- readRDS(file.path(output_dir, "fit_latent.rds"))
fit_midpoint <- readRDS(file.path(output_dir, "fit_midpoint.rds"))

prediction_grid <- seq(min(gb_data$Start_date), max(gb_data$End_date), by = 1)
n_prediction_points <- length(prediction_grid)
year_breaks <- seq(-200, 1000, 200)

# 1. Trend recovery + parameter posteriors -------------------------------------

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

trend_latent   <- summarise_trend(fit_latent)
trend_midpoint <- summarise_trend(fit_midpoint)

# Shared y-limits so the two trend panels are directly comparable.
shared_ylim <- range(trend_latent$lower_90, trend_latent$upper_90,
                     trend_midpoint$lower_90, trend_midpoint$upper_90)

make_trend_panel <- function(trend, title_label) {
  ggplot(trend) +
    geom_ribbon(aes(Year, ymin = lower_90, ymax = upper_90), fill = "grey80") +
    geom_ribbon(aes(Year, ymin = lower_50, ymax = upper_50), fill = "grey60") +
    geom_line(aes(Year, median), linewidth = 0.6, colour = "black") +
    scale_x_continuous(breaks = year_breaks, expand = c(0.01, 0)) +
    coord_cartesian(ylim = shared_ylim) +
    labs(title = title_label, x = "Year", y = "log(Roofed area)") + panel_theme
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
  parameter_posteriors(fit_latent, "Latent dates"),
  parameter_posteriors(fit_midpoint, "Midpoint dates")) %>%
  mutate(Model = factor(Model, c("Latent dates", "Midpoint dates")),
         Region = factor(Region, c("50% CI", "90% CI", "Tail")))

posteriors_panel <- ggplot(both_posteriors, aes(Value, fill = Region)) +
  geom_histogram(colour = "black", linewidth = 0.2, bins = 40) +
  scale_fill_manual(NULL, values = c("50% CI" = "grey40", "90% CI" = "grey65",
                                     "Tail" = "grey88")) +
  facet_grid(Model ~ Parameter, scales = "free_x") +
  labs(title = expression(bold("C.") ~ "Posterior distributions"), x = NULL,
       y = "Count") +
  panel_theme + theme(strip.background = element_blank(),
                      strip.text = element_text(face = "bold", size = 10),
                      legend.position = "bottom")

comparison_figure <-
  (make_trend_panel(trend_latent, expression(bold("A.") ~ "Trend — latent dates")) |
   make_trend_panel(trend_midpoint, expression(bold("B.") ~ "Trend — midpoint dates"))) /
  posteriors_panel + plot_layout(heights = c(1, 1.2))

ggsave(figure_path("dataset_1_model_comparison.png"), comparison_figure,
       width = 12, height = 10, dpi = 300, bg = "white")

# 2. Per-house date posteriors (latent model) ----------------------------------
# No true date to mark; the window and the posterior pulled within it are shown.

latent_draws <- fit_latent$draws(format = "df")

make_date_panel <- function(sample_index) {
  sample_row   <- gb_data[sample_index, ]
  sample_label <- paste0("House ", sample_row$ID,
                         " (value = ", round(sample_row$Value, 2), ")")

  trend_with_sample <- ggplot(trend_latent) +
    geom_ribbon(aes(Year, ymin = lower_90, ymax = upper_90), fill = "grey80") +
    geom_line(aes(Year, median), linewidth = 0.6) +
    annotate("rect", xmin = sample_row$Start_date, xmax = sample_row$End_date,
             ymin = -Inf, ymax = Inf, fill = "grey50", alpha = 0.15) +
    geom_hline(yintercept = sample_row$Value, linetype = "dashed",
               linewidth = 0.5) +
    scale_x_continuous(breaks = year_breaks, expand = c(0.01, 0)) +
    labs(title = paste0(sample_label, " — trend vs observed value"),
         x = "Year", y = "log(Roofed area)") + panel_theme

  estimated_dates <- latent_draws[[paste0("date_actual[", sample_index, "]")]]
  date_density <- ggplot(tibble(estimated_date = estimated_dates),
                         aes(estimated_date)) +
    geom_density(fill = "grey60", colour = "black", linewidth = 0.4, alpha = 0.5) +
    geom_vline(xintercept = c(sample_row$Start_date, sample_row$End_date),
               linetype = "dashed", colour = "grey40", linewidth = 0.5) +
    labs(title = paste0(sample_label, " — posterior date estimate"),
         x = "Estimated date", y = "Density") + panel_theme

  trend_with_sample | date_density
}

set.seed(123)
date_panels <- map(sort(sample(seq_len(nrow(gb_data)), 6)), make_date_panel)
ggsave(figure_path("dataset_1_individual_date_posteriors.png"),
       wrap_plots(date_panels, ncol = 1),
       width = 12, height = 18, dpi = 300, bg = "white")

cat("Wrote dataset_1_model_comparison.png, dataset_1_individual_date_posteriors.png\n")
