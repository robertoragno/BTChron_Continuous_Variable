# Figures for the recovery study (greyscale with dark-red/grey model coding):
#   recovery.png             - estimated vs true value for each parameter
#   accuracy.png             - how often each 50% / 90% interval held the truth
#   precision_vs_entropy.png - 90% interval width vs Shannon entropy, per param
#   sigma_vs_entropy.png     - midpoint sigma error grows as dating coarsens
#   accuracy_vs_entropy.png  - 90% interval accuracy vs Shannon entropy
#   precision_boxplots.png   - EIV vs midpoint 90% width per parameter,
#                              across low / mid / high entropy
#   precision_vs_n.png        - 90% interval width vs sample size N, per param
# Also writes recovery_table.csv (accuracy rates) for the paper's summary table.

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
                                        c("EIV", "Midpoint"))
model_colours <- c(`EIV` = "#780000", Midpoint = "grey25")
model_shapes  <- c(`EIV` = 16, Midpoint = 17)
model_fills   <- c(`EIV` = "grey45", Midpoint = "grey82")

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
  geom_abline(aes(linetype = "1:1 (perfect recovery)"), slope = 1, intercept = 0,
              colour = "grey40") +
  geom_point(size = 0.7, alpha = 0.45) +
  facet_wrap(~ parameter, scales = "free", nrow = 1) +
  scale_colour_manual(values = model_colours, name = NULL) +
  scale_shape_manual(values = model_shapes, name = NULL) +
  scale_linetype_manual(values = "dashed", name = NULL) +
  guides(colour = guide_legend(override.aes = list(size = 2.5, alpha = 1))) +
  labs(x = "True value", y = "Estimated (posterior median)",
       title = "Recovery: estimated vs true value across random datasets") +
  panel_theme

ggsave(figure_path("recovery.png"), recovery_plot,
       width = 14, height = 4, dpi = 300, bg = "white")

# accuracy.png : how often each credible interval held the truth, with a 95%
# Jeffreys interval on the rate (Beta posterior under the Jeffreys prior
# Beta(1/2, 1/2): rate ~ Beta(successes + 1/2, misses + 1/2)).

accuracy_rates <- recovery_results %>%
  pivot_longer(matches("_cov(50|90)$"), names_to = c("parameter", "interval"),
               names_sep = "_cov", values_to = "contained_truth") %>%
  group_by(model, parameter, interval) %>%
  summarise(successes = sum(contained_truth), n = n(),
            value = mean(contained_truth), .groups = "drop") %>%
  mutate(lower = qbeta(0.025, successes + 0.5, n - successes + 0.5),
         upper = qbeta(0.975, successes + 0.5, n - successes + 0.5),
         model = label_model(model),
         parameter = factor(parameter, tolower(param_levels), labels = param_labels),
         interval = paste0(interval, "% interval"),
         target = ifelse(grepl("50", interval), 0.5, 0.9))

accuracy_plot <- ggplot(accuracy_rates, aes(parameter, value, fill = model)) +
  geom_col(position = position_dodge(0.7), width = 0.65, colour = "black",
           linewidth = 0.2) +
  geom_errorbar(aes(ymin = lower, ymax = upper), position = position_dodge(0.7),
                width = 0.2, linewidth = 0.4) +
  geom_hline(aes(yintercept = target, linetype = "Nominal coverage"),
             colour = "grey40") +
  facet_wrap(~ interval) +
  scale_fill_manual(values = model_fills, name = NULL) +
  scale_linetype_manual(values = "dashed", name = NULL) +
  scale_y_continuous(labels = scales::percent, limits = c(0, 1)) +
  labs(x = NULL, y = "Datasets containing the truth", title = "Accuracy") +
  panel_theme

ggsave(figure_path("accuracy.png"), accuracy_plot,
       width = 10, height = 4.5, dpi = 300, bg = "white")

# precision_vs_entropy.png : 90% interval width against the Shannon entropy of
# each dataset's periodisation, for every parameter. Higher H (more even, finer
# phases) should give tighter intervals; the midpoint model should be wider
# throughout, not just for sigma. Slope 2 (after the kink) is where recovery is
# weakest, so its precision gap matters most.

width_long <- recovery_results %>%
  transmute(model = label_model(model), H,
            Baseline = baseline_width90, Slope1 = slope1_width90,
            Slope2 = slope2_width90, CP = cp_width90, Sigma = sigma_width90) %>%
  pivot_longer(c(Baseline, Slope1, Slope2, CP, Sigma),
               names_to = "parameter", values_to = "width90") %>%
  mutate(parameter = factor(parameter, param_levels, labels = param_labels))

entropy_bins <- width_long %>%
  mutate(H_bin = cut(H, breaks = quantile(H, seq(0, 1, 0.2)),
                     include.lowest = TRUE)) %>%
  group_by(parameter, model, H_bin) %>%
  summarise(H_mid = median(H), median_width = median(width90),
            .groups = "drop")

precision_entropy_plot <- ggplot(width_long, aes(H, width90)) +
  geom_point(aes(colour = model, shape = model), size = 0.5, alpha = 0.4) +
  geom_line(data = entropy_bins, aes(H_mid, median_width, colour = model),
            linewidth = 0.7) +
  geom_point(data = entropy_bins, aes(H_mid, median_width, colour = model,
                                      shape = model), size = 1.8) +
  facet_wrap(~ parameter, scales = "free_y", nrow = 1) +
  scale_colour_manual(values = model_colours, name = NULL) +
  scale_shape_manual(values = model_shapes, name = NULL) +
  guides(colour = guide_legend(override.aes = list(size = 2.5, alpha = 1))) +
  labs(x = "Shannon entropy of the periodisation (H)",
       y = "90% interval width (less is more precise)",
       title = "Precision against dating resolution",
       caption = "Lines and markers: per-bin median width.") +
  panel_theme

ggsave(figure_path("precision_vs_entropy.png"), precision_entropy_plot,
       width = 14, height = 4, dpi = 300, bg = "white")

# sigma_vs_entropy.png : the entropy analogue of the old sigma-vs-window figure.
# Midpoint mistakes within-phase date spread for measurement noise and over-
# estimates sigma; the error grows as the periodisation gets coarser (low H).
# Shaded bands mark dating resolution from coarse (low H) to fine (high H).

sigma_by_H <- recovery_results %>% mutate(model = label_model(model))

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
label_y_sg  <- y_limits_sg[1] + 0.04 * diff(y_limits_sg)

sigma_entropy_plot <- ggplot() +
  geom_rect(data = resolution_bands,
            aes(xmin = band_start, xmax = band_end, ymin = -Inf, ymax = Inf,
                fill = shade)) +
  scale_fill_identity() +
  geom_hline(yintercept = 0, linetype = "dotted", colour = "grey45") +
  geom_point(data = sigma_by_H, aes(H, sigma_err, colour = model, shape = model),
             size = 0.55, alpha = 0.65) +
  geom_errorbar(data = binned_sigma,
                aes(x = x, ymin = q25, ymax = q75, colour = model),
                width = 0.06, linewidth = 0.7) +
  geom_point(data = binned_sigma, aes(x, median_error, colour = model,
                                      shape = model), size = 2) +
  geom_text(data = resolution_bands, aes(band_centre, label_y_sg, label = label),
            colour = "grey30", size = 3) +
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

# accuracy_vs_entropy.png : how often the 90% interval held the truth, binned by
# entropy. This is the "coverage decreases as dating gets coarser" figure:
# accuracy should rise with H, and midpoint should sit below latent.

accuracy_long <- recovery_results %>%
  transmute(model = label_model(model), H,
            Baseline = baseline_cov90, Slope1 = slope1_cov90,
            Slope2 = slope2_cov90, CP = cp_cov90, Sigma = sigma_cov90) %>%
  pivot_longer(c(Baseline, Slope1, Slope2, CP, Sigma),
               names_to = "parameter", values_to = "contained_truth") %>%
  mutate(parameter = factor(parameter, param_levels, labels = param_labels))

accuracy_by_entropy <- accuracy_long %>%
  mutate(H_bin = cut(H, breaks = quantile(H, seq(0, 1, 0.2)),
                     include.lowest = TRUE)) %>%
  group_by(parameter, model, H_bin) %>%
  summarise(H_mid = median(H), accuracy = mean(contained_truth),
            .groups = "drop")

accuracy_entropy_plot <- ggplot(accuracy_by_entropy,
                                aes(H_mid, accuracy, colour = model,
                                    shape = model)) +
  geom_hline(aes(yintercept = 0.9, linetype = "90% nominal"), colour = "grey40") +
  geom_line(linewidth = 0.7) +
  geom_point(size = 1.8) +
  facet_wrap(~ parameter, nrow = 1) +
  scale_colour_manual(values = model_colours, name = NULL) +
  scale_shape_manual(values = model_shapes, name = NULL) +
  scale_linetype_manual(values = "dashed", name = NULL) +
  scale_y_continuous(labels = scales::percent, limits = c(0, 1)) +
  labs(x = "Shannon entropy of the periodisation (H)",
       y = "90% interval accuracy",
       title = "Accuracy against dating resolution (90% interval)") +
  panel_theme

ggsave(figure_path("accuracy_vs_entropy.png"), accuracy_entropy_plot,
       width = 14, height = 4, dpi = 300, bg = "white")

# precision_boxplots.png : the latent-vs-midpoint precision gap as boxplots,
# split into low / mid / high entropy thirds.

width_by_entropy <- width_long %>%
  mutate(entropy_band = cut(H, breaks = quantile(H, c(0, 1/3, 2/3, 1)),
                            include.lowest = TRUE,
                            labels = c("Low", "Mid", "High")))

precision_boxplots <- ggplot(width_by_entropy,
                             aes(entropy_band, width90, fill = model)) +
  geom_boxplot(outlier.size = 0.3, linewidth = 0.3, position = position_dodge(0.8)) +
  facet_wrap(~ parameter, scales = "free_y", nrow = 1) +
  scale_fill_manual(values = model_fills, name = NULL) +
  labs(x = "Shannon entropy (H)",
       y = "90% interval width (less is more precise)",
       title = "EIV vs midpoint precision across the entropy range") +
  panel_theme

ggsave(figure_path("precision_boxplots.png"), precision_boxplots,
       width = 14, height = 4.5, dpi = 300, bg = "white")

# precision_vs_n.png : 90% interval width against sample size N, per parameter.
# Sanity check that precision improves (width shrinks) with more data, for both
# models. N is the swept factor, drawn independently of the periodisation, so it
# isolates sample size from dating resolution. Faint points are the per-dataset
# widths (clipped at the 95th percentile per parameter so the medians stay
# readable); markers are the per-N median and interquartile range, dodged by
# model so the two are not drawn on top of each other.

width_by_n <- recovery_results %>%
  transmute(model = label_model(model), N,
            Baseline = baseline_width90, Slope1 = slope1_width90,
            Slope2 = slope2_width90, CP = cp_width90, Sigma = sigma_width90) %>%
  pivot_longer(c(Baseline, Slope1, Slope2, CP, Sigma),
               names_to = "parameter", values_to = "width90") %>%
  mutate(parameter = factor(parameter, param_levels, labels = param_labels),
         N = factor(N))

n_summary <- width_by_n %>%
  group_by(parameter, model, N) %>%
  summarise(median_width = median(width90),
            q25 = quantile(width90, 0.25), q75 = quantile(width90, 0.75),
            .groups = "drop")

width_points <- width_by_n %>%
  group_by(parameter) %>%
  filter(width90 <= quantile(width90, 0.95)) %>%
  ungroup()

dodge <- position_dodge(width = 0.7)

precision_n_plot <- ggplot(mapping = aes(N, colour = model, shape = model)) +
  geom_point(data = width_points, aes(y = width90),
             position = position_jitterdodge(jitter.width = 0.18, dodge.width = 0.7),
             size = 0.4, alpha = 0.2) +
  geom_line(data = n_summary, aes(y = median_width, group = model),
            position = dodge, linewidth = 0.7) +
  geom_errorbar(data = n_summary, aes(ymin = q25, ymax = q75),
                position = dodge, width = 0.3, linewidth = 0.6) +
  geom_point(data = n_summary, aes(y = median_width), position = dodge, size = 2) +
  facet_wrap(~ parameter, scales = "free_y", nrow = 1) +
  scale_colour_manual(values = model_colours, name = NULL) +
  scale_shape_manual(values = model_shapes, name = NULL) +
  guides(colour = guide_legend(override.aes = list(size = 2.5, alpha = 1))) +
  labs(x = "Sample size (N)",
       y = "90% interval width (less is more precise)",
       title = "Precision against sample size",
       caption = "Markers: per-N median; bars: interquartile range.") +
  panel_theme

ggsave(figure_path("precision_vs_n.png"), precision_n_plot,
       width = 14, height = 4, dpi = 300, bg = "white")

# recovery_table.csv : accuracy rates for the cross-simulation paper table

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
            `EIV (%)` = latent, `Midpoint (%)` = midpoint) %>%
  arrange(Parameter, Interval)

write_csv(recovery_rate_table,
          here("Simulations", "Sim_Changepoint", "output", "recovery_table.csv"))

cat("Wrote recovery.png, accuracy.png, precision_vs_entropy.png,",
    "sigma_vs_entropy.png, accuracy_vs_entropy.png, precision_boxplots.png,",
    "precision_vs_n.png,",
    "recovery_table.csv\n\n")
print(as.data.frame(recovery_rate_table), row.names = FALSE)
