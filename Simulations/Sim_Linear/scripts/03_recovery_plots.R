# Figures for the recovery study:
#   accuracy_precision_composite.png - accuracy and precision against dating
#                                       resolution (Shannon entropy H), stacked
#   precision_vs_n.png               - 90% interval width against sample size N
# Also writes recovery_table.csv (accuracy rates) for the paper's summary table.
#
# Six other cuts of this same H-resolved accuracy/precision relationship
# (recovery.png, accuracy.png, accuracy_vs_entropy.png, precision_vs_entropy.png,
# sigma_vs_entropy.png, precision_boxplots.png) were dropped from the main
# script on 2026-07-15 as redundant with the composite above; the code that
# produces them still lives in
# archive/recovery_figures_trim_20260715/recovery_supplementary_figures.R
# if any of them is wanted again.

library(here)
library(readr)
library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)
library(patchwork)

recovery_results <- read_csv(
  here("Simulations", "Sim_Linear", "output", "recovery_results.csv"),
  show_col_types = FALSE
) %>% filter(is.na(error), max_rhat <= 1.05, n_divergent == 0)

figure_path <- function(name) here("Simulations", "Sim_Linear", "figures", name)

label_model  <- function(code) factor(code, c("latent", "midpoint"),
                                       c("EIV", "Midpoint"))
model_colours <- c(`EIV` = "#780000", Midpoint = "grey25")
model_shapes  <- c(`EIV` = 16, Midpoint = 17)
model_fills   <- c(`EIV` = "grey45", Midpoint = "grey82")

panel_theme <- theme_classic(base_size = 11) +
  theme(plot.title = element_text(size = 11, face = "bold"),
        plot.subtitle = element_text(size = 9),
        strip.background = element_blank(),
        strip.text = element_text(face = "bold"),
        legend.position = "top")

# Anchor markers on the entropy axis: where the well-specified scenario cases
# of 05 sit on this continuum (mean H per case). The skewed cases are left out
# on purpose -- they break the uniform-within-phase assumption, so they live
# off this map. Only drawn if the scenario study has been run.
scenario_table_file <- here("Simulations", "Sim_Linear", "output",
                            "scenario_table.csv")
anchor_layers <- if (file.exists(scenario_table_file)) {
  scenario_anchors <- read_csv(scenario_table_file, show_col_types = FALSE) %>%
    filter(case %in% c("fine", "coarse", "diagnostic")) %>%
    group_by(case) %>%
    summarise(H = mean(H), .groups = "drop") %>%
    mutate(label = str_to_sentence(case))
  list(
    geom_vline(data = scenario_anchors, aes(xintercept = H),
               linetype = "dotted", colour = "grey45", linewidth = 0.4),
    geom_text(data = scenario_anchors, aes(H, Inf, label = label),
              inherit.aes = FALSE, angle = 90, hjust = 1.15, vjust = -0.4,
              size = 2.7, colour = "grey35"))
} else NULL

# precision_vs_entropy (panel B of the composite) : 90% interval width against
# the Shannon entropy of each dataset's periodisation, for every parameter.

width_long <- recovery_results %>%
  transmute(model = label_model(model), H,
            Intercept = intercept_width90,
            Slope     = slope_width90,
            Sigma     = sigma_width90) %>%
  pivot_longer(c(Intercept, Slope, Sigma),
               names_to = "parameter", values_to = "width90") %>%
  mutate(parameter = factor(parameter, c("Intercept", "Slope", "Sigma")))

entropy_bins <- width_long %>%
  mutate(H_bin = cut(H, breaks = quantile(H, seq(0, 1, 0.2)),
                     include.lowest = TRUE)) %>%
  group_by(parameter, model, H_bin) %>%
  summarise(H_mid = median(H), median_width = median(width90),
            q25 = quantile(width90, 0.25), q75 = quantile(width90, 0.75),
            .groups = "drop")

precision_entropy_plot <- ggplot(width_long, aes(H, width90)) +
  geom_point(aes(colour = model, shape = model), size = 0.5, alpha = 0.4) +
  geom_line(data = entropy_bins, aes(H_mid, median_width, colour = model),
            linewidth = 0.7) +
  geom_point(data = entropy_bins, aes(H_mid, median_width, colour = model,
                                      shape = model), size = 1.8) +
  facet_wrap(~ parameter, scales = "free_y") +
  scale_colour_manual(values = model_colours, name = NULL) +
  scale_shape_manual(values = model_shapes, name = NULL) +
  guides(colour = guide_legend(override.aes = list(size = 2.5, alpha = 1))) +
  labs(x = "Shannon entropy of the periodisation (H)",
       y = "Width of the 90% interval",
       title = "Precision against dating resolution",
       caption = "Lines and markers: per-bin median width.") +
  panel_theme +
  anchor_layers

# accuracy_vs_entropy (panel A of the composite) : how often the 90% interval
# held the truth, binned by entropy.

accuracy_long <- recovery_results %>%
  transmute(model = label_model(model), H,
            Intercept = intercept_cov90, Slope = slope_cov90,
            Sigma = sigma_cov90) %>%
  pivot_longer(c(Intercept, Slope, Sigma),
               names_to = "parameter", values_to = "contained_truth") %>%
  mutate(parameter = factor(parameter, c("Intercept", "Slope", "Sigma")))

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
  facet_wrap(~ parameter) +
  scale_colour_manual(values = model_colours, name = NULL) +
  scale_shape_manual(values = model_shapes, name = NULL) +
  scale_linetype_manual(values = "dashed", name = NULL) +
  scale_y_continuous(labels = scales::percent, limits = c(0, 1)) +
  labs(x = "Shannon entropy of the periodisation (H)",
       y = "% of 90% intervals containing the true value",
       title = "Accuracy against dating resolution") +
  panel_theme +
  anchor_layers

# accuracy_precision_composite.png : accuracy (A) over precision (B), per
# parameter, so the two must be read together. A narrow interval that misses
# the truth cannot pass for a good one, and a narrow interval that still
# covers reads as the win it is.

acc_panel  <- accuracy_entropy_plot + labs(title = NULL, x = NULL) +
  facet_wrap(~ parameter, scales = "free_y")
prec_panel <- precision_entropy_plot +
  labs(title = NULL, caption = NULL) +
  guides(colour = "none", shape = "none") +
  theme(strip.text = element_blank())

accuracy_precision_composite <- (acc_panel / prec_panel) +
  plot_layout(guides = "collect", axes = "collect_x") +
  plot_annotation(tag_levels = "A",
                  title = "Accuracy and precision against dating resolution") &
  theme(legend.position = "top",
        plot.title = element_text(size = 11, face = "bold"),
        plot.tag = element_text(size = 15, face = "bold"),
        plot.tag.location = "margin")

ggsave(figure_path("accuracy_precision_composite.png"), accuracy_precision_composite,
       width = 10, height = 7.5, dpi = 300, bg = "white")

# precision_vs_n.png : 90% interval width against sample size N, per parameter.
# Sanity check that precision improves (width shrinks) with more data, for both
# models. N is the swept factor, drawn independently of the periodisation, so it
# isolates sample size from dating resolution. Faint points are the per-dataset
# widths (clipped at the 95th percentile per parameter so the medians stay
# readable); markers are the per-N median and interquartile range, dodged by
# model so the two are not drawn on top of each other.

width_by_n <- recovery_results %>%
  transmute(model = label_model(model), N,
            Intercept = intercept_width90, Slope = slope_width90,
            Sigma = sigma_width90) %>%
  pivot_longer(c(Intercept, Slope, Sigma),
               names_to = "parameter", values_to = "width90") %>%
  mutate(parameter = factor(parameter, c("Intercept", "Slope", "Sigma")),
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
  facet_wrap(~ parameter, scales = "free_y") +
  scale_colour_manual(values = model_colours, name = NULL) +
  scale_shape_manual(values = model_shapes, name = NULL) +
  guides(colour = guide_legend(override.aes = list(size = 2.5, alpha = 1))) +
  labs(x = "Sample size (N)",
       y = "Width of the 90% interval",
       title = "Precision against sample size",
       caption = "Markers: per-N median; bars: interquartile range.") +
  panel_theme

ggsave(figure_path("precision_vs_n.png"), precision_n_plot,
       width = 10, height = 4, dpi = 300, bg = "white")

# recovery_table.csv : accuracy (% of datasets where the interval caught the
# truth) per parameter and model. Appends into the cross-simulation paper table.

accuracy_table <- recovery_results %>%
  pivot_longer(matches("_cov(50|90)$"), names_to = c("parameter", "interval"),
               names_sep = "_cov", values_to = "contained_truth") %>%
  group_by(parameter, interval, model) %>%
  summarise(percent = round(100 * mean(contained_truth), 1), .groups = "drop") %>%
  pivot_wider(names_from = model, values_from = percent) %>%
  transmute(Simulation = "Linear",
            Parameter = factor(str_to_title(parameter),
                               c("Intercept", "Slope", "Sigma")),
            Interval = paste0(interval, "%"),
            `EIV (%)` = latent, `Midpoint (%)` = midpoint) %>%
  arrange(Parameter, Interval)

write_csv(accuracy_table,
          here("Simulations", "Sim_Linear", "output", "recovery_table.csv"))


print(as.data.frame(accuracy_table), row.names = FALSE)

#' recovery_width_summary.csv : how much narrower EIV's interval is, per
#' parameter. Paired by dataset_id (both models fit the same simulated data)
#' and summarised by median, not mean, because width90 is heavy-tailed 
#' (some datasets have a very wide interval for one model), so a mean
#' ratio is pulled around by those few datasets. 
#' We can use this number in the README.

width_pairs <- recovery_results %>%
  select(dataset_id, model, intercept_width90, slope_width90, sigma_width90) %>%
  pivot_longer(c(intercept_width90, slope_width90, sigma_width90),
               names_to = "parameter", names_pattern = "(.*)_width90",
               values_to = "width90") %>%
  pivot_wider(names_from = model, values_from = width90) %>%
  filter(!is.na(latent), !is.na(midpoint))

width_summary <- width_pairs %>%
  group_by(parameter) %>%
  summarise(median_eiv = median(latent), median_midpoint = median(midpoint),
            eiv_narrower_pct = round(100 * (1 - median(latent) / median(midpoint)), 1),
            .groups = "drop") %>%
  mutate(parameter = factor(str_to_title(parameter),
                            c("Intercept", "Slope", "Sigma"))) %>%
  arrange(parameter)

write_csv(width_summary,
          here("Simulations", "Sim_Linear", "output", "recovery_width_summary.csv"))
print(as.data.frame(width_summary), row.names = FALSE)
