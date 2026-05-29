#' Purpose: Generate simulated data with a known changepoint regression
#'          f(t) = baseline + slope_1*(t - t_min)                    for t <= cp
#'          f(t) = baseline + slope_1*(cp - t_min) + slope_2*(t - cp) for t >  cp
#'          The function is continuous at the changepoint.

library(here)
library(tidyverse)

set.seed(42)

# ── Known parameters (to be recovered) ──────────────────────────────────────

gt_baseline <- 8
gt_slope_1  <- 0.02
gt_slope_2  <- -0.01
gt_cp       <- 500          # changepoint year

f_true <- function(t) {
  t_min <- 100
  ifelse(
    t <= gt_cp,
    gt_baseline + gt_slope_1 * (t - t_min),
    gt_baseline + gt_slope_1 * (gt_cp - t_min) + gt_slope_2 * (t - gt_cp)
  )
}

# ── Generate date ranges ────────────────────────────────────────────────────

n         <- 300
date_grid <- seq(100, 900, by = 25)

sim_data <- tibble(
  ID        = 1:n,
  Site_name = sample(c("Site A", "Site B", "Site C", "Site D"), n, replace = TRUE)
) %>%
  rowwise() %>%
  mutate(
    Start_date = sample(date_grid[date_grid <= 700], 1),
    End_date = {
      possible_ends <- date_grid[date_grid > Start_date & date_grid <= 900]
      if (length(possible_ends) < 2) {
        sample(possible_ends, 1)
      } else {
        weights <- dnorm(
          seq_along(possible_ends),
          mean = length(possible_ends) * 0.35,
          sd   = length(possible_ends) * 0.25
        )
        sample(possible_ends, 1, prob = weights)
      }
    }
  ) %>%
  ungroup()

set.seed(123)
flip <- sample(c(TRUE, FALSE), n, replace = TRUE, prob = c(0.35, 0.65))
sim_data <- sim_data %>%
  mutate(End_date = if_else(flip, End_date - 1L, End_date))

# ── Generate observed values from the true function ─────────────────────────

sim_data <- sim_data %>%
  rowwise() %>%
  mutate(
    True_date  = runif(1, min = Start_date, max = End_date),
    True_value = f_true(True_date),
    Value      = round(True_value + rnorm(1, 0, 1.5), 1),
    Value      = pmax(Value, 0.5)
  ) %>%
  ungroup()

write_csv(
  sim_data %>% select(ID, Site_name, Start_date, End_date, True_date, Value),
  here("Simulations", "Sim_Changepoint", "data", "simulated_data.csv")
)

# ── Ground truth on a fine grid ─────────────────────────────────────────────

ground_truth <- tibble(
  Year       = seq(100, 900, by = 1),
  True_value = f_true(Year)
)
write_csv(ground_truth,
          here("Simulations", "Sim_Changepoint", "data", "ground_truth.csv"))

# ── Save the generating parameters ──────────────────────────────────────────

params <- tibble(
  parameter = c("baseline", "slope_1", "slope_2", "changepoint", "sigma_noise"),
  value     = c(gt_baseline, gt_slope_1, gt_slope_2, gt_cp, 1.5)
)
write_csv(params,
          here("Simulations", "Sim_Changepoint", "data", "generating_parameters.csv"))

glimpse(sim_data)
summary(sim_data$Value)
