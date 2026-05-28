#' Purpose: Generate simulated data to be tested in three regressions

# Load libraries
# Need to add 'here' for reproducibility
library(tidyverse)

# Set seed for reproducibility
set.seed(42)


n <- 300

date_grid <- seq(100, 900, by = 25)

sim_data <- tibble(
  ID = 1:n,
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
          sd = length(possible_ends) * 0.25
        )
        sample(possible_ends, 1, prob = weights)
      }
    }
  ) %>%
  ungroup()

# Randomly convert some end dates to N-1 style (e.g. 199, 249, 299)
# to mimic how archaeologists sometimes write date ranges
set.seed(123)
flip <- sample(c(TRUE, FALSE), n, replace = TRUE, prob = c(0.35, 0.65))
sim_data <- sim_data %>%
  mutate(End_date = if_else(flip, End_date - 1L, End_date))

# Ground truth: smooth rise-and-fall peaking around 450 CE
# f(t) = baseline + amplitude * exp(-((t - peak) / width)^2)
gt_baseline  <- 8
gt_amplitude <- 12
gt_peak      <- 450
gt_width     <- 200

f_true <- function(t) {
  gt_baseline + gt_amplitude * exp(-((t - gt_peak) / gt_width)^2)
}

sim_data <- sim_data %>%
  rowwise() %>%
  mutate(
    True_date  = runif(1, min = Start_date, max = End_date),
    True_value = f_true(True_date),
    Value      = round(True_value + rnorm(1, 0, 2.5), 1),
    Value      = pmax(Value, 0.5)
  ) %>%
  ungroup() %>%
  select(-True_date, -True_value)

write_csv(sim_data, here::here("Simulations", "data", "simulated_data.csv"))

ground_truth <- tibble(
  Year       = seq(100, 900, by = 1),
  True_value = f_true(Year)
)
write_csv(ground_truth, here::here("Simulations", "data", "ground_truth.csv"))

glimpse(sim_data)
summary(sim_data$Value)
