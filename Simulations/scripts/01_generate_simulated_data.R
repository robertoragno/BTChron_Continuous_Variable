library(tidyverse)

set.seed(42)

n <- 500

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

# Simulated continuous archaeological variable (e.g. vessel capacity in litres)
sim_data <- sim_data %>%
  mutate(
    Value = round(rnorm(n, mean = 15, sd = 5), 1),
    Value = pmax(Value, 0.5)
  )

write_csv(sim_data, here::here("Simulations", "data", "simulated_data.csv"))

glimpse(sim_data)
