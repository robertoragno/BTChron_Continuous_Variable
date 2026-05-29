#' Purpose: Posterior date densities for a random subset of individual samples,
#'          greyscale version of MSCA Model_1_Linear figure 04.

library(here)
library(tidyverse)
library(cmdstanr)
library(posterior)
library(patchwork)

# ── Load data and fit ────────────────────────────────────────────────────────

sim_data <- read_csv(here("Simulations", "Sim_Linear", "data", "simulated_data.csv"),
                     show_col_types = FALSE)

fit <- readRDS(here("Simulations", "Sim_Linear", "output", "fit_linear.rds"))

# ── Shared theme (matches exploratory panel) ─────────────────────────────────

theme_panel <- theme_classic(base_size = 11) +
  theme(
    plot.title    = element_text(size = 11, face = "bold", hjust = 0),
    plot.subtitle = element_text(size = 9, colour = "grey40", hjust = 0),
    plot.margin   = margin(5, 10, 2, 10)
  )

# ── Extract latent date posteriors ───────────────────────────────────────────

draws   <- fit$draws(format = "df")
N_obs   <- nrow(sim_data)

set.seed(123)
sample_ids <- sort(sample(1:N_obs, 6))

date_long <- map_dfr(sample_ids, function(i) {
  col <- paste0("true_date_actual[", i, "]")
  tibble(
    sample_idx     = i,
    estimated_date = draws[[col]],
    Start_date     = sim_data$Start_date[i],
    End_date       = sim_data$End_date[i],
    True_date      = sim_data$True_date[i],
    ID             = sim_data$ID[i]
  )
})

# ── Plot ─────────────────────────────────────────────────────────────────────

p <- ggplot(date_long, aes(x = estimated_date)) +
  geom_density(fill = "grey60", colour = "black", linewidth = 0.4, alpha = 0.5) +
  geom_vline(aes(xintercept = Start_date),
             linetype = "dashed", colour = "grey40", linewidth = 0.5) +
  geom_vline(aes(xintercept = End_date),
             linetype = "dashed", colour = "grey40", linewidth = 0.5) +
  geom_vline(aes(xintercept = True_date),
             linetype = "solid", colour = "black", linewidth = 0.5) +
  facet_wrap(~ ID, scales = "free_y", ncol = 2) +
  labs(
    title    = "Posterior date estimates for individual samples",
    subtitle = "Dashed lines: TPQ–TAQ range. Solid line: true generating date.",
    x = "Estimated date (CE)",
    y = "Density"
  ) +
  theme_panel +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(face = "bold", size = 10)
  )

ggsave(here("Simulations", "Sim_Linear", "figures", "individual_date_posteriors.png"),
       p, width = 10, height = 8, dpi = 300, bg = "white")

cat("Saved to", here("Simulations", "Sim_Linear", "figures", "individual_date_posteriors.png"), "\n")
