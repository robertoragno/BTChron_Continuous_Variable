# Recovery / calibration plots for the GP simulation.
#
# Reads the repeated-study output from 02_recovery_study.R and plots posterior
# rank histograms for sigma, together with the observed 50% and 90% coverage.

library(here)
library(tidyverse)
library(patchwork)

source(here("Simulations", "Sim_GP", "scripts", "simulate.R"))

results_latent <- read_csv(here("Simulations", "Sim_GP", "output", "calibration_results.csv"),
                           show_col_types = FALSE)
results_midpoint <- read_csv(
  here("Simulations", "Sim_GP", "output", "calibration_results_midpoint.csv"),
  show_col_types = FALSE
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

ggsave(here("Simulations", "Sim_GP", "figures", "calibration_coverage.png"),
       p, width = 8, height = 7, dpi = 300, bg = "white")

cat("Wrote calibration_coverage.png\n")
