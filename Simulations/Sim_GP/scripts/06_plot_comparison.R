#' Purpose: Side-by-side comparison of the latent-date HSGP model vs the midpoint model.
#'   Row 1 — Trend recovery (A: latent dates, B: midpoint dates)
#'   Row 2 — Posterior sigma distributions for both models, sigma is the cleanest
#'            diagnostic: the midpoint model absorbs date uncertainty into residual
#'            variance, so its sigma posterior sits above the true value of 2.5.

library(here)
library(tidyverse)
library(cmdstanr)
library(posterior)
library(patchwork)

sim_data     <- read_csv(here("Simulations", "Sim_GP", "data", "simulated_data.csv"),
                         show_col_types = FALSE)
ground_truth <- read_csv(here("Simulations", "Sim_GP", "data", "ground_truth.csv"),
                         show_col_types = FALSE)
true_params  <- read_csv(here("Simulations", "Sim_GP", "data", "generating_parameters.csv"),
                         show_col_types = FALSE)

fit_latent   <- readRDS(here("Simulations", "Sim_GP", "output", "fit_hsgp.rds"))
fit_midpoint <- readRDS(here("Simulations", "Sim_GP", "output", "fit_hsgp_midpoint.rds"))

true_sigma <- true_params$value[true_params$parameter == "sigma_noise"]

theme_panel <- theme_classic(base_size = 11) +
  theme(
    plot.title  = element_text(size = 11, face = "bold", hjust = 0),
    plot.margin = margin(5, 10, 2, 10)
  )

pred_grid <- seq(min(sim_data$Start_date), max(sim_data$End_date), by = 1)
N_pred    <- length(pred_grid)

extract_trend <- function(fit) {
  draws     <- fit$draws(format = "df")
  trend_mat <- as.matrix(draws[, paste0("mu_pred[", 1:N_pred, "]")])
  tibble(
    Year     = pred_grid,
    Median   = apply(trend_mat, 2, median),
    Lower_90 = apply(trend_mat, 2, quantile, 0.05),
    Upper_90 = apply(trend_mat, 2, quantile, 0.95),
    Lower_50 = apply(trend_mat, 2, quantile, 0.25),
    Upper_50 = apply(trend_mat, 2, quantile, 0.75)
  )
}

trend_latent   <- extract_trend(fit_latent)
trend_midpoint <- extract_trend(fit_midpoint)

make_trend_panel <- function(df, title_label) {
  ggplot(df) +
    geom_ribbon(aes(x = Year, ymin = Lower_90, ymax = Upper_90), fill = "grey80") +
    geom_ribbon(aes(x = Year, ymin = Lower_50, ymax = Upper_50), fill = "grey60") +
    geom_line(aes(x = Year, y = Median), linewidth = 0.6, colour = "black") +
    geom_line(data = ground_truth, aes(x = Year, y = True_value),
              linewidth = 0.7, colour = "black", linetype = "dashed") +
    scale_x_continuous(breaks = seq(100, 900, 200), expand = c(0.01, 0)) +
    labs(title = title_label, x = "Year CE", y = "Value") +
    theme_panel
}

p_trend_latent <- make_trend_panel(
  trend_latent,
  expression(bold("A.") ~ "Trend recovery — latent dates")
)
p_trend_mid <- make_trend_panel(
  trend_midpoint,
  expression(bold("B.") ~ "Trend recovery — midpoint dates")
)

# Sigma posteriors for both models with CI shading and true value

extract_sigma <- function(fit, model_label) {
  draws   <- fit$draws(format = "df")
  sigma   <- draws$sigma
  q05     <- quantile(sigma, 0.05)
  q25     <- quantile(sigma, 0.25)
  q75     <- quantile(sigma, 0.75)
  q95     <- quantile(sigma, 0.95)

  tibble(sigma = sigma, Model = model_label) %>%
    mutate(
      Region = case_when(
        sigma >= q25 & sigma <= q75 ~ "50% CI",
        sigma >= q05 & sigma <= q95 ~ "90% CI",
        TRUE                        ~ "Tail"
      ),
      Region = factor(Region, levels = c("50% CI", "90% CI", "Tail"))
    )
}

sigma_df <- bind_rows(
  extract_sigma(fit_latent,   "Latent dates"),
  extract_sigma(fit_midpoint, "Midpoint dates")
) %>%
  mutate(Model = factor(Model, levels = c("Latent dates", "Midpoint dates")))

p_sigma <- ggplot(sigma_df, aes(x = sigma, fill = Region)) +
  geom_histogram(colour = "black", linewidth = 0.2, bins = 40) +
  geom_vline(xintercept = true_sigma,
             linetype = "dashed", linewidth = 0.6, colour = "black") +
  scale_fill_manual(
    name   = NULL,
    values = c("50% CI" = "grey40", "90% CI" = "grey65", "Tail" = "grey88")
  ) +
  facet_wrap(~ Model, scales = "free_y", nrow = 1) +
  labs(
    title   = expression(bold("C.") ~ "Posterior sigma — latent vs midpoint"),
    caption = "Dashed line: true generating sigma (2.5).",
    x = "sigma",
    y = "Count"
  ) +
  theme_panel +
  theme(
    strip.background = element_blank(),
    strip.text       = element_text(face = "bold", size = 10),
    legend.position  = "bottom",
    legend.text      = element_text(size = 8),
    legend.key.size  = unit(10, "pt"),
    plot.caption     = element_text(hjust = 0, size = 8, colour = "grey40")
  )

p <- (p_trend_latent | p_trend_mid) / p_sigma +
  plot_layout(heights = c(1, 0.8)) +
  plot_annotation(theme = theme(plot.margin = margin(5, 5, 5, 5)))

ggsave(here("Simulations", "Sim_GP", "figures", "model_comparison.png"),
       p, width = 12, height = 9, dpi = 300, bg = "white")

cat("Saved to", here("Simulations", "Sim_GP", "figures", "model_comparison.png"), "\n")
