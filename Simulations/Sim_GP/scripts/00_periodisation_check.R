# Prior predictive check. Draws generating parameters from the same priors as the
# recovery study, simulates datasets, and shows what the priors imply: the spread
# of dating resolution (H), example phases within a timeline, and example
# datasets together with the wavy curves they were drawn from.

library(here)
library(ggplot2)
library(patchwork)

source(here("Simulations", "Sim_GP", "scripts", "simulate.R"))

figure_path <- function(name) here("Simulations", "Sim_GP", "figures", name)

panel_theme <- theme_classic(base_size = 11) +
  theme(plot.title = element_text(size = 11, face = "bold"),
        legend.position = "top")

set.seed(2026)

# Periodisation priors, matching 02_recovery_study.R and shared with the linear
# and changepoint simulations. K phases with Dirichlet concentration alpha_conc;
# alpha_conc leans low (Beta over [0.1, 10]) so coarser phases are common. N is
# the two-level low/high sweep. Defined as functions so they are not accidentally
# reused as stored values.
draw_K          <- function(n = 1) sample(3:10, n, replace = TRUE)
draw_alpha_conc <- function(n = 1) 10^(-1 + 2 * rbeta(n, 2, 3.5))

draw_prior <- function() list(
  sigma      = runif(1, 0.5, 4),
  K          = draw_K(),
  alpha_conc = draw_alpha_conc(),
  N          = sample(c(50, 200), 1),
  curve      = draw_curve()
)

# Shannon entropy H (~ dating resolution) implied by the K and alpha_conc priors.
n_draws <- 5000
K          <- draw_K(n_draws)
alpha_conc <- draw_alpha_conc(n_draws)
H <- vapply(seq_len(n_draws),
            function(i) attr(partition_timeline(1, K[i], alpha_conc[i]), "H"),
            numeric(1))
res <- data.frame(K = K, H = H)

h_hist <- ggplot(res, aes(H)) +
  geom_histogram(bins = 40, fill = "grey70", colour = "black", linewidth = 0.2) +
  labs(title = "Dating resolution H", x = "H", y = "Count") + panel_theme

h_by_k <- ggplot(res, aes(factor(K), H)) +
  geom_boxplot(fill = "grey80", outlier.size = 0.3, linewidth = 0.3) +
  labs(title = "H by number of phases", x = "Phases (K)", y = "H") + panel_theme

ggsave(figure_path("periodisation_H_distribution.png"),
       h_hist | h_by_k, width = 11, height = 4.5, dpi = 300, bg = "white")

# Example phase layouts. Boundaries are rebuilt from the weights, so each panel
# shows the whole K-phase partition, not just where one observation fell.
example_seeds <- c(1, 7, 13, 21, 42, 99)
phases <- do.call(rbind, lapply(seq_along(example_seeds), function(i) {
  set.seed(example_seeds[i])
  k          <- draw_K()
  alpha_conc <- draw_alpha_conc()
  w    <- attr(partition_timeline(1, k, alpha_conc), "weights")
  bound <- round(TMIN + c(0, cumsum(w)) * (TMAX - TMIN))
  data.frame(example = i, phase = seq_len(k),
             label = sprintf("K=%d, alpha=%.2f, H=%.2f",
                             k, alpha_conc, -sum(w * log(w))),
             start = bound[-length(bound)], end = bound[-1])
}))
phases$panel <- paste0("Example ", phases$example, "  (", phases$label, ")")

layout_plot <- ggplot(phases) +
  geom_rect(aes(xmin = start, xmax = end, ymin = 0, ymax = 1,
                fill = factor(phase %% 2)), colour = "black", linewidth = 0.2) +
  scale_fill_manual(values = c("grey85", "grey60"), guide = "none") +
  facet_wrap(~ panel, ncol = 1) +
  scale_x_continuous(breaks = seq(TMIN, TMAX, 100)) +
  labs(title = "Example periodisations", x = "Year CE", y = NULL) +
  panel_theme +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())

ggsave(figure_path("periodisation_examples.png"), layout_plot,
       width = 9, height = 8, dpi = 300, bg = "white")

# Prior predictive datasets: draw full parameter sets, simulate, and look at the
# data generated. The pooled value distribution flags whether the priors imply
# plausible measurements; the panels show six datasets with the curve behind them.
n_sets <- 500
priors <- replicate(n_sets, draw_prior(), simplify = FALSE)
sims   <- lapply(priors, function(p)
  simulate_gp(p$N, p$sigma, p$K, p$alpha_conc, p$curve))

value_hist <- ggplot(data.frame(Value = unlist(lapply(sims, `[[`, "Value"))),
                     aes(Value)) +
  geom_histogram(bins = 50, fill = "grey70", colour = "black", linewidth = 0.2) +
  labs(title = "Predictive distribution of values",
       x = "Simulated value", y = "Count") + panel_theme

panel_of <- function(i) sprintf("Set %d  (N=%d, H=%.2f)",
                                i, priors[[i]]$N, attr(sims[[i]], "H"))

# Stack the six datasets into one frame, each tagged with its panel label.
ex <- do.call(rbind, lapply(1:6, function(i) {
  d <- sims[[i]]; d$panel <- panel_of(i); d
}))

# The wavy curve each set was generated from, on a yearly grid.
curves <- do.call(rbind, lapply(1:6, function(i) {
  g <- truth_grid(priors[[i]]$curve)
  data.frame(panel = panel_of(i), Year = g$Year, Value = g$True_value)
}))

dataset_plot <- ggplot(ex, aes((Start_date + End_date) / 2, Value)) +
  geom_linerange(aes(xmin = Start_date, xmax = End_date),
                 colour = "grey70", linewidth = 0.2, alpha = 0.5) +
  geom_point(size = 0.5, alpha = 0.6) +
  geom_line(data = curves, aes(Year, Value), colour = "#780000", linewidth = 0.6) +
  facet_wrap(~ panel, scales = "free_y") +
  scale_x_continuous(breaks = seq(TMIN, TMAX, 200)) +
  labs(title = "Example datasets drawn from the priors", x = "Year CE", y = "Value") +
  panel_theme

ggsave(figure_path("prior_predictive.png"),
       value_hist / dataset_plot + plot_layout(heights = c(1, 1.6)),
       width = 11, height = 9, dpi = 300, bg = "white")
