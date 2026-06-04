#' Purpose: Prior predictive check for the HSGP model. Draw hyperparameters
#'          from the priors, compute GP realisations through the HSGP
#'          parameterisation, and compare the resulting trends and simulated
#'          observations against the observed data range. If the priors are
#'          reasonable the observed values should sit comfortably within the
#'          bulk of the prior predictive distribution.

library(here)
library(tidyverse)
library(patchwork)

set.seed(42)

sim_data <- read_csv(here("Simulations", "Sim_GP", "data", "simulated_data.csv"),
                     show_col_types = FALSE)

# ── HSGP setup (must match sim_hsgp.stan) ────────────────────────────────────

M          <- 20
c_boundary <- 1.5
L          <- c_boundary

time_min   <- min(sim_data$Start_date)
time_max   <- max(sim_data$End_date)
time_range <- time_max - time_min

pred_grid      <- seq(time_min, time_max, by = 5)
x_pred_norm    <- (pred_grid - time_min) / time_range

PHI_fn <- function(x_norm, M, L) {
  m_seq <- seq_len(M)
  sin(outer(pi / (2 * L) * (x_norm + L), m_seq)) / sqrt(L)
}

diagSPD_EQ_fn <- function(alpha, rho, L, M) {
  m_seq <- seq_len(M)
  alpha * sqrt(sqrt(2 * pi) * rho) *
    exp(-0.25 * (rho * pi / (2 * L))^2 * m_seq^2)
}

PHI_mat <- PHI_fn(x_pred_norm, M, L)

# ── Draw from priors ──────────────────────────────────────────────────────────
# Priors from sim_hsgp.stan:
#   mu    ~ normal(0, 15)
#   alpha ~ normal+(0, 10)   [half-normal, truncated at 0]
#   rho   ~ inv_gamma(5, 2.5)
#   sigma ~ exponential(1)
#   z[m]  ~ std_normal()

N_sim <- 200

rinvgamma <- function(n, shape, scale) scale / rgamma(n, shape = shape, rate = 1)

prior_draws <- tibble(
  mu    = rnorm(N_sim, 0, 15),
  alpha = abs(rnorm(N_sim, 0, 10)),
  rho   = rinvgamma(N_sim, shape = 5, scale = 2.5),
  sigma = rexp(N_sim, rate = 1)
)

# ── Compute prior predictive trends and observations ─────────────────────────

trend_list <- vector("list", N_sim)
y_sim_all  <- numeric(0)

for (i in seq_len(N_sim)) {
  z       <- rnorm(M)
  spd     <- diagSPD_EQ_fn(prior_draws$alpha[i], prior_draws$rho[i], L, M)
  beta_gp <- spd * z
  f       <- prior_draws$mu[i] + PHI_mat %*% beta_gp

  trend_list[[i]] <- tibble(
    sim  = i,
    Year = pred_grid,
    f    = as.numeric(f)
  )

  y_sim_all <- c(y_sim_all, rnorm(length(f), f, prior_draws$sigma[i]))
}

trend_df <- bind_rows(trend_list)

# ── Trim extreme draws for legible plot (keep middle 90% by max |f|) ─────────

max_abs_f <- trend_df %>%
  group_by(sim) %>%
  summarise(max_f = max(abs(f)), .groups = "drop") %>%
  filter(max_f <= quantile(max_f, 0.90))

trend_plot_df <- trend_df %>% filter(sim %in% max_abs_f$sim)

# ── Theme ─────────────────────────────────────────────────────────────────────

theme_panel <- theme_classic(base_size = 11) +
  theme(
    plot.title  = element_text(size = 11, face = "bold", hjust = 0),
    plot.margin = margin(5, 10, 5, 10)
  )

# ── Panel A: spaghetti of prior predictive trends ────────────────────────────

obs_range <- range(sim_data$Value)

p_trends <- ggplot() +
  geom_line(data = trend_plot_df,
            aes(x = Year, y = f, group = sim),
            colour = "grey70", linewidth = 0.2, alpha = 0.5) +
  annotate("rect",
           xmin = -Inf, xmax = Inf,
           ymin = obs_range[1], ymax = obs_range[2],
           fill = "steelblue", alpha = 0.12) +
  annotate("text",
           x = time_min + 0.02 * time_range,
           y = mean(obs_range),
           label = "Observed\ndata range",
           hjust = 0, vjust = 0.5, size = 3, colour = "steelblue4") +
  scale_x_continuous(breaks = seq(100, 900, 200), expand = c(0.01, 0)) +
  coord_cartesian(ylim = c(-50, 80)) +
  labs(
    title = expression(bold("A.") ~ "Prior predictive trends (n = 200, trimmed to 90th percentile)"),
    x     = "Year CE",
    y     = "f(t)"
  ) +
  theme_panel

# ── Panel B: density of prior predictive y vs observed y ─────────────────────

density_df <- bind_rows(
  tibble(Source = "Prior predictive", Value = y_sim_all),
  tibble(Source = "Observed",         Value = sim_data$Value)
)

p_density <- ggplot(density_df, aes(x = Value, colour = Source, fill = Source)) +
  geom_density(alpha = 0.25, linewidth = 0.6) +
  scale_colour_manual(values = c("Prior predictive" = "grey40", "Observed" = "steelblue4"),
                      name = NULL) +
  scale_fill_manual(values  = c("Prior predictive" = "grey70", "Observed" = "steelblue"),
                    name = NULL) +
  coord_cartesian(xlim = c(-20, 50)) +
  labs(
    title = expression(bold("B.") ~ "Prior predictive y vs observed values"),
    x     = "Value",
    y     = "Density"
  ) +
  theme_panel +
  theme(
    legend.position = c(0.97, 0.97),
    legend.justification = c(1, 1),
    legend.background = element_rect(fill = alpha("white", 0.8), colour = NA),
    legend.text = element_text(size = 9)
  )

p <- p_trends / p_density

ggsave(here("Simulations", "Sim_GP", "figures", "prior_predictive_check.png"),
       p, width = 10, height = 8, dpi = 300, bg = "white")

cat("Saved to", here("Simulations", "Sim_GP", "figures",
                     "prior_predictive_check.png"), "\n")
