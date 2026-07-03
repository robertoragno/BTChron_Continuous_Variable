# Figures for the GP recovery study:
#   accuracy_vs_entropy.png   - accuracy vs Shannon entropy, curve and sigma
#   precision_vs_entropy.png  - band / interval width vs Shannon entropy
#   precision_vs_n.png        - width vs sample size N
#   sigma_vs_entropy.png      - midpoint sigma error grows as dating coarsens
#   curve_recovery.png        - distribution of curve accuracy and curve RMSE
# Also writes recovery_table.csv (accuracy) for the paper's cross-sim table.
#
# The GP targets are the curve (proportion of the true curve inside the band,
# averaged over a time grid) and sigma (a scalar). Curve accuracy is a proportion
# per dataset; sigma accuracy is a hit / miss as in the other simulations. Both
# are summarised as the mean across datasets.

library(here)
library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)

figure_path <- function(name) here("Simulations", "Sim_GP", "figures", name)

label_model   <- function(code) factor(code, c("latent", "midpoint"), c("EIV", "Midpoint"))
model_colours <- c(`EIV` = "#780000", Midpoint = "grey25")
model_shapes  <- c(`EIV` = 16, Midpoint = 17)
model_fills   <- c(`EIV` = "grey45", Midpoint = "grey82")

panel_theme <- theme_classic(base_size = 11) +
  theme(plot.title = element_text(size = 11, face = "bold"),
        plot.subtitle = element_text(size = 9),
        strip.background = element_blank(),
        strip.text = element_text(face = "bold"),
        legend.position = "top")

recovery_results_all <- read_csv(
  here("Simulations", "Sim_GP", "output", "recovery_results.csv"),
  show_col_types = FALSE
) %>% filter(is.na(error))

# The latent-date GP is the hardest of the three simulations to sample, and the
# curve and the noise scale can fail to mix independently, so convergence is
# checked per target: curve figures keep fits whose grid curve (mu_pred) mixed;
# the sigma figure keeps fits whose sigma mixed. Report the kept counts so the
# drop is visible.
curve_results <- recovery_results_all %>%
  filter(curve_rhat <= 1.05, n_divergent == 0) %>% mutate(model = label_model(model))
sigma_results <- recovery_results_all %>%
  filter(sigma_rhat <= 1.05, n_divergent == 0) %>% mutate(model = label_model(model))

recovery_results_all %>%
  group_by(model) %>%
  summarise(total = n(),
            curve_kept = sum(curve_rhat <= 1.05 & n_divergent == 0),
            sigma_kept = sum(sigma_rhat <= 1.05 & n_divergent == 0),
            .groups = "drop") %>%
  { cat("Converged fits kept per target (Rhat <= 1.05, no divergences):\n"); print(as.data.frame(.)) }

# Helper: bin a metric by entropy quintile and average within model x bin.
bin_by_entropy <- function(df, value_col) {
  df %>%
    mutate(H_bin = cut(H, breaks = quantile(H, seq(0, 1, 0.2)), include.lowest = TRUE)) %>%
    group_by(target, model, H_bin) %>%
    summarise(H_mid = median(H), value = mean(.data[[value_col]]), .groups = "drop")
}

# accuracy_vs_entropy.png : mean accuracy against the entropy of each dataset's
# periodisation. Curve accuracy is the mean proportion of the true curve inside
# the band; sigma accuracy is the fraction of datasets whose interval held the
# truth. Higher H (finer phases) should raise accuracy; midpoint should sit below
# EIV. Curve and sigma use their own converged subsets.

coverage_long <- bind_rows(
  curve_results %>% transmute(model, H, `Curve (90%)` = curve_cov90,
                              `Curve (50%)` = curve_cov50) %>%
    pivot_longer(starts_with("Curve"), names_to = "target", values_to = "coverage"),
  sigma_results %>% transmute(model, H, `Sigma (90%)` = as.numeric(sigma_cov90),
                              `Sigma (50%)` = as.numeric(sigma_cov50)) %>%
    pivot_longer(starts_with("Sigma"), names_to = "target", values_to = "coverage")
) %>%
  mutate(target = factor(target, c("Curve (90%)", "Curve (50%)",
                                   "Sigma (90%)", "Sigma (50%)")))

coverage_by_entropy <- coverage_long %>%
  rename(value_in = coverage) %>% bin_by_entropy("value_in") %>% rename(coverage = value)

nominal_lines <- tibble(
  target = factor(c("Curve (90%)", "Curve (50%)", "Sigma (90%)", "Sigma (50%)"),
                  c("Curve (90%)", "Curve (50%)", "Sigma (90%)", "Sigma (50%)")),
  nominal = c(0.9, 0.5, 0.9, 0.5))

accuracy_entropy_plot <- ggplot(coverage_by_entropy,
                                aes(H_mid, coverage, colour = model, shape = model)) +
  geom_hline(data = nominal_lines, aes(yintercept = nominal),
             linetype = "dashed", colour = "grey40") +
  geom_line(linewidth = 0.7) +
  geom_point(size = 1.8) +
  facet_wrap(~ target, nrow = 1) +
  scale_colour_manual(values = model_colours, name = NULL) +
  scale_shape_manual(values = model_shapes, name = NULL) +
  scale_y_continuous(labels = scales::percent, limits = c(0, 1)) +
  labs(x = "Shannon entropy of the periodisation (H)", y = "Accuracy",
       title = "Accuracy against dating resolution",
       subtitle = "Curve panels: proportion of the true curve inside the CI. Sigma panels: proportion of datasets whose CI held the true value.") +
  panel_theme

ggsave(figure_path("accuracy_vs_entropy.png"), accuracy_entropy_plot,
       width = 12, height = 4, dpi = 300, bg = "white")

# precision_vs_entropy.png : band / interval width against entropy. Curve width is
# the mean 90% band width over the grid; sigma width is the 90% interval width.
# Own units, so free y-axes; own converged subsets.

width_long <- bind_rows(
  curve_results %>% transmute(model, H, target = "Curve", width90 = curve_width90),
  sigma_results %>% transmute(model, H, target = "Sigma", width90 = sigma_width90)
) %>% mutate(target = factor(target, c("Curve", "Sigma")))

width_by_entropy <- width_long %>%
  mutate(H_bin = cut(H, breaks = quantile(H, seq(0, 1, 0.2)), include.lowest = TRUE)) %>%
  group_by(target, model, H_bin) %>%
  summarise(H_mid = median(H), median_width = median(width90),
            q25 = quantile(width90, 0.25), q75 = quantile(width90, 0.75),
            .groups = "drop")

precision_entropy_plot <- ggplot(width_long, aes(H, width90)) +
  geom_point(aes(colour = model, shape = model), size = 0.5, alpha = 0.4) +
  geom_line(data = width_by_entropy, aes(H_mid, median_width, colour = model),
            linewidth = 0.7) +
  geom_point(data = width_by_entropy, aes(H_mid, median_width, colour = model,
                                          shape = model), size = 1.8) +
  facet_wrap(~ target, scales = "free_y") +
  scale_colour_manual(values = model_colours, name = NULL) +
  scale_shape_manual(values = model_shapes, name = NULL) +
  guides(colour = guide_legend(override.aes = list(size = 2.5, alpha = 1))) +
  labs(x = "Shannon entropy of the periodisation (H)",
       y = "Width of the 90% band / interval",
       title = "Precision against dating resolution",
       caption = "Lines and markers: per-bin median width.") +
  panel_theme

ggsave(figure_path("precision_vs_entropy.png"), precision_entropy_plot,
       width = 9, height = 4, dpi = 300, bg = "white")

# precision_vs_n.png : width against sample size N, the swept factor. Precision
# should improve (width shrink) with more data for both models.

width_by_n <- bind_rows(
  curve_results %>% transmute(model, N = factor(N), target = "Curve", width90 = curve_width90),
  sigma_results %>% transmute(model, N = factor(N), target = "Sigma", width90 = sigma_width90)
) %>% mutate(target = factor(target, c("Curve", "Sigma")))

n_summary <- width_by_n %>%
  group_by(target, model, N) %>%
  summarise(median_width = median(width90),
            q25 = quantile(width90, 0.25), q75 = quantile(width90, 0.75),
            .groups = "drop")

dodge <- position_dodge(width = 0.5)
precision_n_plot <- ggplot(n_summary, aes(N, median_width, colour = model, shape = model)) +
  geom_line(aes(group = model), position = dodge, linewidth = 0.7) +
  geom_errorbar(aes(ymin = q25, ymax = q75), position = dodge, width = 0.2,
                linewidth = 0.6) +
  geom_point(position = dodge, size = 2) +
  facet_wrap(~ target, scales = "free_y") +
  scale_colour_manual(values = model_colours, name = NULL) +
  scale_shape_manual(values = model_shapes, name = NULL) +
  labs(x = "Sample size (N)", y = "Width of the 90% band / interval",
       title = "Precision against sample size",
       caption = "Markers: per-N median; bars: interquartile range.") +
  panel_theme

ggsave(figure_path("precision_vs_n.png"), precision_n_plot,
       width = 9, height = 4, dpi = 300, bg = "white")

# sigma_vs_entropy.png : midpoint mistakes within-phase date spread for
# measurement noise and over-estimates sigma; the error grows as the periodisation
# coarsens (low H). Same figure as the linear and changepoint studies, for the
# cross-simulation story. Uses the sigma-converged subset.

sigma_by_H <- sigma_results
H_breaks <- seq(0, 2.5, by = 0.5)
resolution_bands <- tibble(
  band_start = H_breaks[-length(H_breaks)], band_end = H_breaks[-1],
  band_centre = (H_breaks[-length(H_breaks)] + H_breaks[-1]) / 2,
  label = c("Very coarse", "Coarse", "Moderate", "Fine", "Very fine"),
  shade = c("white", "grey97", "grey94", "grey91", "grey88"))

binned_sigma <- sigma_by_H %>%
  mutate(bin_centre = (floor(pmin(H, 2.499) / 0.5) + 0.5) * 0.5) %>%
  group_by(model, bin_centre) %>%
  summarise(median_error = median(sigma_err),
            q25 = quantile(sigma_err, 0.25),
            q75 = quantile(sigma_err, 0.75), .groups = "drop") %>%
  mutate(x = bin_centre + ifelse(model == "EIV", -0.04, 0.04))

y_limits_sg <- quantile(sigma_by_H$sigma_err, c(0.004, 0.996))
label_y_sg  <- y_limits_sg[2] - 0.04 * diff(y_limits_sg)

sigma_entropy_plot <- ggplot() +
  geom_rect(data = resolution_bands,
            aes(xmin = band_start, xmax = band_end, ymin = -Inf, ymax = Inf,
                fill = shade)) +
  scale_fill_identity() +
  geom_hline(yintercept = 0, linetype = "dotted", colour = "grey45") +
  geom_point(data = sigma_by_H, aes(H, sigma_err, colour = model, shape = model),
             size = 0.35, alpha = 0.35) +
  geom_errorbar(data = binned_sigma,
                aes(x = x, ymin = q25, ymax = q75, colour = model),
                width = 0.06, linewidth = 0.7) +
  geom_point(data = binned_sigma, aes(x, median_error, colour = model,
                                      shape = model), size = 2) +
  geom_text(data = resolution_bands, aes(band_centre, label_y_sg, label = label),
            colour = "grey20", size = 3.3, fontface = "bold") +
  scale_colour_manual(values = model_colours, name = NULL) +
  scale_shape_manual(values = model_shapes, name = NULL) +
  guides(colour = guide_legend(override.aes = list(size = 2.5, alpha = 1))) +
  coord_cartesian(xlim = c(0, max(H_breaks)), ylim = y_limits_sg, expand = FALSE) +
  scale_x_continuous(breaks = seq(0, 2.5, 0.5)) +
  labs(x = "Shannon entropy (H)",
       y = expression(sigma[err] == sigma[est] - sigma[true]),
       title = "Midpoint over-estimates the noise as dating gets coarser",
       caption = "Markers: bin median; bars: interquartile range.") +
  panel_theme

ggsave(figure_path("sigma_vs_entropy.png"), sigma_entropy_plot,
       width = 8.5, height = 5.5, dpi = 300, bg = "white")

# curve_recovery.png : how well the curve is recovered overall. Panel A is the
# distribution of curve accuracy (proportion of the true curve inside the 90% CI);
# B is the curve RMSE (posterior median curve vs truth). Curve-converged subset.

coverage_hist <- ggplot(curve_results, aes(curve_cov90, fill = model)) +
  geom_histogram(position = "identity", alpha = 0.55, bins = 25,
                 colour = "black", linewidth = 0.2) +
  geom_vline(xintercept = 0.9, linetype = "dashed", colour = "grey40") +
  scale_fill_manual(values = model_fills, name = NULL) +
  scale_x_continuous(labels = scales::percent) +
  labs(title = expression(bold("A.") ~ "Curve accuracy"),
       subtitle = "Proportion of the true curve inside the 90% CI",
       x = "Proportion of the true curve inside the 90% CI", y = "Count") +
  panel_theme

rmse_box <- ggplot(curve_results, aes(model, curve_rmse, fill = model)) +
  geom_boxplot(outlier.size = 0.4, linewidth = 0.3, width = 0.6) +
  scale_fill_manual(values = model_fills, name = NULL, guide = "none") +
  labs(title = expression(bold("B.") ~ "Curve error"),
       x = NULL, y = "RMSE (median curve vs truth)") +
  panel_theme

ggsave(figure_path("curve_recovery.png"),
       (coverage_hist | rmse_box) + plot_layout(widths = c(1.4, 1)),
       width = 10, height = 4, dpi = 300, bg = "white")

# recovery_table.csv : accuracy per target and model. Curve = mean % of the true
# curve inside the interval; Sigma = % of datasets whose interval held the truth.
# Each target uses its own converged subset. Appends into the cross-sim table.

curve_rows <- curve_results %>%
  transmute(model, `Curve_50` = curve_cov50, `Curve_90` = curve_cov90) %>%
  pivot_longer(-model, names_to = c("parameter", "interval"), names_sep = "_",
               values_to = "coverage")
sigma_rows <- sigma_results %>%
  transmute(model, `Sigma_50` = as.numeric(sigma_cov50), `Sigma_90` = as.numeric(sigma_cov90)) %>%
  pivot_longer(-model, names_to = c("parameter", "interval"), names_sep = "_",
               values_to = "coverage")

accuracy_table <- bind_rows(curve_rows, sigma_rows) %>%
  group_by(parameter, interval, model) %>%
  summarise(percent = round(100 * mean(coverage), 1), .groups = "drop") %>%
  pivot_wider(names_from = model, values_from = percent) %>%
  transmute(Simulation = "GP",
            Parameter = factor(parameter, c("Curve", "Sigma")),
            Interval = paste0(interval, "%"),
            `EIV (%)` = EIV, `Midpoint (%)` = Midpoint) %>%
  arrange(Parameter, Interval)

write_csv(accuracy_table,
          here("Simulations", "Sim_GP", "output", "recovery_table.csv"))

print(as.data.frame(accuracy_table), row.names = FALSE)
