# Figures for the scenario study (05_scenarios.R).
#
# Four fixed cases, many datasets each, both models. The question here is not the
# general relationship (that is 02/03) but a clean side-by-side: in each case,
# does EIV help, and how? Three figures and a table:
#   scenario_accuracy.png   - % of 90% intervals holding the truth, per case
#   scenario_precision.png  - 90% interval width, per case
#   scenario_sigma.png      - sigma error (midpoint reads date spread as noise)
#   scenario_table.csv      - the quotable per-case numbers

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

# Cases in a fixed order, with readable labels.
case_levels <- c("fine", "coarse", "skewed", "diagnostic")
case_labs   <- c("Fine\ndating", "Coarse\ndating", "Skewed\ndeposition",
                 "Diagnostic\nsampling")
label_case  <- function(code) factor(code, case_levels, case_labs)

param_levels <- c("Intercept", "Slope", "Sigma")

panel_theme <- theme_classic(base_size = 11) +
  theme(plot.title = element_text(size = 11, face = "bold"),
        strip.background = element_blank(),
        strip.text = element_text(face = "bold"),
        legend.position = "top")

# scenario_accuracy.png : % of 90% intervals containing the truth, per case and
# parameter, with a 95% Jeffreys interval on the rate.
accuracy <- results %>%
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

# scenario_sigma.png : sigma error. The midpoint reads within-phase date spread
# as measurement noise and overestimates sigma; the effect is largest where the
# sample is coarsely dated and smallest under diagnostic sampling.
sigma_error <- results %>%
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
       width = 8, height = 4.5, dpi = 300, bg = "white")

# scenario_table.csv : the quotable per-case numbers.
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
