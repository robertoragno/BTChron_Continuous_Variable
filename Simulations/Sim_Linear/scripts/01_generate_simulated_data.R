#' Purpose: Generate simulated data with a known linear trend
#'          f(t) = baseline + slope * t
#'          Parameters are stored so they can be recovered by the model.

library(here)
library(tidyverse)

set.seed(42)

# Known parameters (to be recovered from the model later)

true_baseline <- 5
true_slope    <- 0.015

f_true <- function(t) {
  true_baseline + true_slope * t
}

# Generate date ranges

n         <- 300 # number of observations
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

#' Introduce some "censoring" by randomly shortening some end dates by 1 year
#' as the flip of a coin. This simulates the real-world scenario where some observations
#' might end a century (for instance the second century) as 199 CE instead of 200 CE.
censor_flag <- sample(c(TRUE, FALSE), n, replace = TRUE, prob = c(0.35, 0.65))
sim_data <- sim_data %>%
  mutate(End_date = if_else(censor_flag, End_date - 1L, End_date))

# Generate observed values from the true function 

sim_data <- sim_data %>%
  rowwise() %>%
  mutate(
    True_date  = runif(1, min = Start_date, max = End_date), # a random date within the range
    True_value = f_true(True_date), # the true value from the linear function at that date
    Value      = round(True_value + rnorm(1, 0, 1.5), 1), # add some noise to the true value
    Value      = pmax(Value, 0.5) # ensure values are not negative
  ) %>%
  ungroup()

write_csv(
  sim_data %>% select(ID, Site_name, Start_date, End_date, True_date, Value),
  here("Simulations", "Sim_Linear", "data", "simulated_data.csv")
)

# ── Ground truth on a fine grid ─────────────────────────────────────────────

ground_truth <- tibble(
  Year       = seq(100, 900, by = 1),
  True_value = f_true(Year)
)
write_csv(ground_truth,
          here("Simulations", "Sim_Linear", "data", "ground_truth.csv"))

# ── Save the generating parameters ──────────────────────────────────────────

params <- tibble(
  parameter = c("baseline", "slope", "sigma_noise"),
  value     = c(true_baseline, true_slope, 1.5)
)
write_csv(params,
          here("Simulations", "Sim_Linear", "data", "generating_parameters.csv"))

glimpse(sim_data)
summary(sim_data$Value)
