# Figures for the scenario study (05_scenarios.R).
#
# Nine fixed cases, many datasets each, both models. Seven are general-purpose
# (periodisation width, between-phase sampling, within-phase skew with no
# assumption about which periods are better-dated -- see the header of
# 05_scenarios.R); two (coarse_strong_early/late) are a worked example, not a
# general claim. Five figures and a table:
#   scenario_accuracy_precision.png - periodisation width and sampling only
#                                     (fine, coarse, diagnostic): accuracy and
#                                     precision, never bias
#   scenario_skew_bias.png          - slope bias vs skew strength, faceted by
#                                     phase width, with no width~position
#                                     correlation assumed (the general result)
#   position_bias_example.png       - worked example: what happens IF phase
#                                     width correlates with position, in
#                                     either direction (illustrative only)
#   scenario_sigma.png              - sigma error across the seven
#                                     general-purpose cases
#   scenario_table.csv              - the quotable per-case numbers, all nine

library(here)
library(readr)
library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)
library(patchwork)

results <- read_csv(
  here("Simulations", "Sim_Linear", "output", "scenarios.csv"),
  show_col_types = FALSE) %>%
  filter(is.na(error), max_rhat <= 1.05, n_divergent == 0)

figure_path <- function(name) here("Simulations", "Sim_Linear", "figures", name)

# House style, shared with 03.
label_model   <- function(code) factor(code, c("latent", "midpoint"),
                                       c("EIV", "Midpoint"))
model_colours <- c(EIV = "#780000", Midpoint = "grey25")
model_shapes  <- c(EIV = 16, Midpoint = 17)
model_fills   <- c(EIV = "grey45", Midpoint = "grey82")

# Cases in a fixed order, with readable labels. coarse_strong_early/late are
# not part of the general-purpose case grid (see 05_scenarios.R) -- they only
# feed position_bias_example.png -- but are included here so scenario_table.csv
# reports them too.
case_levels <- c("fine", "coarse", "diagnostic",
                 "fine_mild", "fine_strong", "coarse_mild", "coarse_strong",
                 "coarse_strong_early", "coarse_strong_late")
case_labs   <- c("Fine\ndating", "Coarse\ndating", "Diagnostic\nsampling",
                 "Fine,\nmild skew", "Fine,\nstrong skew",
                 "Coarse,\nmild skew", "Coarse,\nstrong skew",
                 "Coarse, strong skew\n(width~position, early; illustrative)",
                 "Coarse, strong skew\n(width~position, late; illustrative)")
label_case  <- function(code) factor(code, case_levels, case_labs)

# The seven general-purpose cases, excluding the two illustrative ones --
# used to keep scenario_sigma.png (and other general figures) to the cases
# that make no assumption about width-position correlation.
general_case_levels <- setdiff(case_levels, c("coarse_strong_early", "coarse_strong_late"))

# The skew cases decompose into a width factor (fine/coarse) and a
# skew-strength factor (none/mild/strong); fine and coarse themselves are the
# skew = "none" reference points (skew_shape = 1, i.e. uniform deposition).
skew_width_levels <- c("fine", "coarse")
skew_width_labs   <- c("Fine dating", "Coarse dating")
skew_level_levels <- c("none", "mild", "strong")
label_skew_width <- function(case) {
  factor(dplyr::case_when(
    case %in% c("fine", "fine_mild", "fine_strong")     ~ "fine",
    case %in% c("coarse", "coarse_mild", "coarse_strong") ~ "coarse"),
    skew_width_levels, skew_width_labs)
}
label_skew_level <- function(case) {
  factor(dplyr::case_when(
    case %in% c("fine", "coarse")               ~ "none",
    grepl("_mild$", case)                       ~ "mild",
    grepl("_strong$", case)                     ~ "strong"),
    skew_level_levels)
}

param_levels <- c("Intercept", "Slope", "Sigma")

panel_theme <- theme_classic(base_size = 11) +
  theme(plot.title = element_text(size = 11, face = "bold"),
        strip.background = element_blank(),
        strip.text = element_text(face = "bold"),
        legend.position = "top")

# scenario_accuracy_precision.png : periodisation width and between-phase
# sampling only (fine, coarse, diagnostic). Neither of those can bias a
# linear slope, so this figure is deliberately restricted to the three cases
# where within-phase deposition is uniform -- the skew cases belong to
# scenario_skew_bias.png instead, where bias is the point.
baseline_cases <- c("fine", "coarse", "diagnostic")

# % of 90% intervals containing the truth, per case and parameter, with a 95%
# Jeffreys interval on the rate.
accuracy <- results %>%
  filter(case %in% baseline_cases) %>%
  pivot_longer(c(intercept_cov90, slope_cov90, sigma_cov90),
               names_to = "parameter", names_pattern = "(.*)_cov90",
               values_to = "hit") %>%
  group_by(case, model, parameter) %>%
  summarise(successes = sum(hit), n = n(), value = mean(hit), .groups = "drop") %>%
  mutate(lower = qbeta(0.025, successes + 0.5, n - successes + 0.5),
         upper = qbeta(0.975, successes + 0.5, n - successes + 0.5),
         model = label_model(model), case = label_case(case),
         parameter = factor(str_to_title(parameter), param_levels))

# Panel A: accuracy as points with a 95% Jeffreys interval, matching the house
# style of accuracy_precision_composite.png (EIV dark red circle, midpoint grey
# triangle, dashed 90% nominal).
accuracy_plot <- ggplot(accuracy, aes(case, value, colour = model,
                                      shape = model)) +
  geom_hline(aes(yintercept = 0.9, linetype = "90% nominal"), colour = "grey40") +
  geom_errorbar(aes(ymin = lower, ymax = upper), position = position_dodge(0.5),
                width = 0.25, linewidth = 0.5) +
  geom_point(position = position_dodge(0.5), size = 2.4) +
  facet_wrap(~ parameter) +
  scale_colour_manual(values = model_colours, name = NULL) +
  scale_shape_manual(values = model_shapes, name = NULL) +
  scale_linetype_manual(values = "dashed", name = NULL) +
  scale_y_continuous(labels = scales::percent, limits = c(0, 1)) +
  labs(x = NULL, y = "Accuracy\n(% covering truth)") +
  panel_theme

# Panel B: precision as the per-dataset widths (faint) with a bold median marker,
# same faint-scatter-plus-marker look as the recovery composite.
precision <- results %>%
  filter(case %in% baseline_cases) %>%
  pivot_longer(c(intercept_width90, slope_width90, sigma_width90),
               names_to = "parameter", names_pattern = "(.*)_width90",
               values_to = "width90") %>%
  mutate(model = label_model(model), case = label_case(case),
         parameter = factor(str_to_title(parameter), param_levels))

precision_medians <- precision %>%
  group_by(case, model, parameter) %>%
  summarise(median_width = median(width90), .groups = "drop")

precision_plot <- ggplot(precision, aes(case, width90, colour = model,
                                        shape = model)) +
  geom_point(position = position_jitterdodge(jitter.width = 0.15,
                                             dodge.width = 0.5),
             size = 0.5, alpha = 0.25) +
  geom_point(data = precision_medians, aes(case, median_width),
             position = position_dodge(0.5), size = 2.4) +
  facet_wrap(~ parameter, scales = "free_y") +
  scale_colour_manual(values = model_colours, name = NULL) +
  scale_shape_manual(values = model_shapes, name = NULL) +
  labs(x = NULL, y = "Precision\n(interval width)") + panel_theme

# scenario_accuracy_precision.png : accuracy (A) over precision (B), per
# parameter, so the two must be read together. A narrow interval that misses the
# truth (midpoint sigma) cannot pass for a good one, and a narrow interval that
# still covers (EIV slope) reads as the win it is. Same layout and palette as
# accuracy_precision_composite.png so the two sit together in the paper.
acc_panel  <- accuracy_plot + labs(x = NULL) +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
prec_panel <- precision_plot +
  guides(colour = "none", shape = "none") +
  theme(strip.text = element_blank())

composite <- (acc_panel / prec_panel) +
  plot_layout(guides = "collect", axes = "collect_x") +
  plot_annotation(tag_levels = "A",
                  title = "Accuracy and precision by scenario") &
  theme(legend.position = "top",
        plot.title = element_text(size = 11, face = "bold"),
        plot.tag = element_text(size = 15, face = "bold"),
        plot.tag.location = "margin")

ggsave(figure_path("scenario_accuracy_precision.png"), composite,
       width = 10, height = 7.5, dpi = 300, bg = "white")

#' scenario_skew_bias.png : Within-phase deposition skew
#' is the only thing that can bias a linear slope; skew = "none" is the
#' fine/coarse cases themselves (skew_shape = 1, uniform), so the same runs
#' used for the accuracy/precision figure double as the zero-skew anchor here.
#' Slope bias is summarised as the ratio of estimated to true slope
#' (median over reps, with the 25/75% quantile band) rather than raw error,
#' because true slopes are drawn near zero and a signed error is not
#' comparable across datasets with different true slope magnitudes; a ratio of
#' 1 is unbiased, below 1 is a slope flattened toward zero.
skew_bias <- results %>%
  filter(case %in% c("fine", "coarse", "fine_mild", "fine_strong",
                     "coarse_mild", "coarse_strong")) %>%
  mutate(model = label_model(model),
         width = label_skew_width(case),
         skew  = label_skew_level(case),
         slope_ratio = slope_med / true_slope) %>%
  group_by(width, skew, model) %>%
  summarise(median_ratio = median(slope_ratio),
            lower = quantile(slope_ratio, 0.25),
            upper = quantile(slope_ratio, 0.75), .groups = "drop")

skew_bias_plot <- ggplot(skew_bias, aes(skew, median_ratio, colour = model,
                                        shape = model, group = model)) +
  geom_hline(aes(yintercept = 1, linetype = "Unbiased"), colour = "grey40") +
  geom_errorbar(aes(ymin = lower, ymax = upper), position = position_dodge(0.3),
                width = 0.15, linewidth = 0.5) +
  geom_line(position = position_dodge(0.3), linewidth = 0.4, alpha = 0.6) +
  geom_point(position = position_dodge(0.3), size = 2.4) +
  facet_wrap(~ width) +
  scale_colour_manual(values = model_colours, name = NULL) +
  scale_shape_manual(values = model_shapes, name = NULL) +
  scale_linetype_manual(values = "dashed", name = NULL) +
  labs(x = "Within-phase deposition skew", y = "Slope ratio\n(estimate / true)",
       title = "Slope bias under skewed within-phase deposition") +
  panel_theme

ggsave(figure_path("scenario_skew_bias.png"), skew_bias_plot,
       width = 8, height = 4.5, dpi = 300, bg = "white")

# position_bias_example.png : a worked example, not a general-purpose result.
# coarse_strong has phase width uncorrelated with position on the timeline
# (the general-purpose default); coarse_strong_early / coarse_strong_late
# turn that correlation on in each direction. This is NOT a claim that real
# periodisations work either way -- which direction applies, if either,
# depends on a specific dataset's evidentiary record (e.g. which periods
# happen to have well-studied typologies or coinage). The point of showing
# both directions side by side is that the bias flips sign with the
# correlation, which is the actual mechanism -- not "old" or "young" -- and
# an analyst should check their own periodisation for this correlation
# before trusting a trend estimate from either model.
position_levels <- c("coarse_strong", "coarse_strong_early", "coarse_strong_late")
position_labs   <- c("No width~position\ncorrelation", "Wide phases\nearly",
                     "Wide phases\nlate")

position_bias <- results %>%
  filter(case %in% position_levels) %>%
  mutate(model = label_model(model),
         position = factor(case, position_levels, position_labs),
         slope_ratio = slope_med / true_slope) %>%
  group_by(position, model) %>%
  summarise(median_ratio = median(slope_ratio),
            lower = quantile(slope_ratio, 0.25),
            upper = quantile(slope_ratio, 0.75), .groups = "drop")

position_bias_plot <- ggplot(position_bias, aes(position, median_ratio,
                                                colour = model, shape = model)) +
  geom_hline(aes(yintercept = 1, linetype = "Unbiased"), colour = "grey40") +
  geom_errorbar(aes(ymin = lower, ymax = upper), position = position_dodge(0.4),
                width = 0.15, linewidth = 0.5) +
  geom_point(position = position_dodge(0.4), size = 2.4) +
  scale_colour_manual(values = model_colours, name = NULL) +
  scale_shape_manual(values = model_shapes, name = NULL) +
  scale_linetype_manual(values = "dashed", name = NULL) +
  labs(x = NULL, y = "Slope ratio\n(estimate / true)",
       title = "Worked example: slope bias when phase width correlates with position",
       subtitle = "Illustrative only -- direction depends on the dataset, not assumed here") +
  panel_theme

ggsave(figure_path("position_bias_example.png"), position_bias_plot,
       width = 7, height = 4.5, dpi = 300, bg = "white")

# scenario_sigma.png : sigma error, general-purpose cases only (excludes the
# two illustrative position cases). The midpoint reads within-phase date
# spread as measurement noise and overestimates sigma; the effect is largest
# where the sample is coarsely dated and smallest under diagnostic sampling.
sigma_error <- results %>%
  filter(case %in% general_case_levels) %>%
  mutate(model = label_model(model), case = label_case(case))

sigma_plot <- ggplot(sigma_error, aes(case, sigma_err, fill = model)) +
  geom_hline(yintercept = 0, linetype = "dotted", colour = "grey45") +
  geom_boxplot(outlier.size = 0.3, linewidth = 0.3,
               position = position_dodge(0.8)) +
  scale_fill_manual(values = model_fills, name = NULL) +
  labs(x = NULL, y = expression(sigma[est] - sigma[true]),
       title = "Sigma error by scenario") +
  panel_theme

ggsave(figure_path("scenario_sigma.png"), sigma_plot,
       width = 11, height = 4.5, dpi = 300, bg = "white")

# scenario_table.csv : per-case numbers.
scenario_table <- results %>%
  group_by(case, model) %>%
  summarise(H = round(mean(H), 2), sampled_width = round(mean(mean_width)),
            slope_cov90 = round(100 * mean(slope_cov90), 1),
            slope_width90 = round(mean(slope_width90), 4),
            # 1 = unbiased; below 1 the estimated slope is flattened toward 0.
            # A median of ratios, robust to the near-zero true slopes.
            slope_ratio = round(median(slope_med / true_slope), 2),
            sigma_cov90 = round(100 * mean(sigma_cov90), 1),
            sigma_err = round(mean(sigma_err), 2),
            n_fit = n(), .groups = "drop") %>%
  mutate(case = factor(case, case_levels), model = label_model(model)) %>%
  arrange(case, model)

write_csv(scenario_table,
          here("Simulations", "Sim_Linear", "output", "scenario_table.csv"))
print(as.data.frame(scenario_table), row.names = FALSE)
