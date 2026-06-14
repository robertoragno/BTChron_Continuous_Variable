# Shared data-generating process for the GP simulation.
# Sourced by the worked example and the repeated calibration study so both use
# the same bell-shaped truth and date-window generator.

library(tidyverse)

TMIN <- 100
TMAX <- 900

# These are the settings of the bell-shaped curve we use as the hidden truth.
# The GP is not trying to recover these names directly; they just define the
# example signal we simulate from.
TRUE_BASELINE <- 8
TRUE_AMPLITUDE <- 12
TRUE_PEAK <- 450
TRUE_WIDTH <- 200
TRUE_SIGMA <- 2.5

SITE_NAMES <- c("Site A", "Site B", "Site C", "Site D")
DATE_GRID <- seq(TMIN, TMAX, by = 25)

# The underlying smooth trend. Later we will only observe noisy values at
# unknown dates inside each sample's dating window.
f_true <- function(t,
                   baseline = TRUE_BASELINE,
                   amplitude = TRUE_AMPLITUDE,
                   peak = TRUE_PEAK,
                   width = TRUE_WIDTH) {
  baseline + amplitude * exp(-((t - peak) / width)^2)
}

simulate_date_windows <- function(n, seed = 42) {
  set.seed(seed)

  # First create one dating window per observation.
  sim_data <- tibble(
    ID = seq_len(n),
    Site_name = sample(SITE_NAMES, n, replace = TRUE)
  ) %>%
    rowwise() %>%
    mutate(
      # Start dates can begin anywhere up to 700 so there is still room for the
      # window to extend to the right.
      Start_date = sample(DATE_GRID[DATE_GRID <= 700], 1),
      End_date = {
        possible_ends <- DATE_GRID[DATE_GRID > Start_date & DATE_GRID <= TMAX]
        if (length(possible_ends) < 2) {
          sample(possible_ends, 1)
        } else {
          # This weighting makes shorter and medium windows more common than
          # extreme ones, while still leaving some wide windows in the mix.
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

  set.seed(seed + 81)
  # A small adjustment for some samples so the end date is not always sitting
  # exactly on the 25-year grid.
  censor_flag <- sample(c(TRUE, FALSE), n, replace = TRUE, prob = c(0.35, 0.65))
  sim_data %>%
    mutate(End_date = if_else(censor_flag, End_date - 1L, End_date))
}

simulate_gp_dataset <- function(n = 300, seed_windows = 42, seed_values = 123) {
  sim_data <- simulate_date_windows(n, seed = seed_windows)

  set.seed(seed_values)
  sim_data %>%
    rowwise() %>%
    mutate(
      # The "real" date is hidden from the model and sampled somewhere inside
      # the dating window.
      True_date = runif(1, min = Start_date, max = End_date),
      # This is the noiseless value on the true curve at that hidden date.
      True_value = f_true(True_date),
      # Add Gaussian observation noise around the true curve value, then round to
      # one decimal place to mimic a recorded continuous measurement.
      Value = round(True_value + rnorm(1, 0, TRUE_SIGMA), 1)
    ) %>%
    ungroup()
}

ground_truth_grid <- function(by = 1) {
  # A dense version of the true curve, used only for plotting and comparison.
  tibble(
    Year = seq(TMIN, TMAX, by = by),
    True_value = f_true(Year)
  )
}

generating_parameters <- function() {
  tibble(
    parameter = c("baseline", "amplitude", "peak", "width", "sigma_noise"),
    value = c(TRUE_BASELINE, TRUE_AMPLITUDE, TRUE_PEAK, TRUE_WIDTH, TRUE_SIGMA)
  )
}

evaluation_targets <- function() {
  # For the GP, these are more useful checks than the generator's named
  # parameters: a known scalar target (sigma) and simple features of the curve.
  truth <- ground_truth_grid()
  half_height <- min(truth$True_value) + 0.5 * (max(truth$True_value) - min(truth$True_value))
  above_half <- which(truth$True_value >= half_height)
  fwhm_years <- truth$Year[max(above_half)] - truth$Year[min(above_half)]
# fhwm_years means "full width at half maximum" in years. 
# It is a measure of the width of the peak of the curve, specifically the distance between 
# the two points on the curve where the value is half of the maximum value. 
# In this context, it is calculated by finding the years corresponding to the maximum and 
#  minimum indices of the values that are above half of the maximum value, and then taking the difference between those two years. 
# This gives an indication of how wide or narrow the peak of the curve is.

  tibble(
    target = c("peak_year", "peak_value", "fwhm_years", "sigma_noise"),
    value = c(TRUE_PEAK, f_true(TRUE_PEAK), fwhm_years, TRUE_SIGMA),
    description = c(
      "Year where the true curve reaches its maximum",
      "Maximum value of the true curve",
      "Full width at half maximum of the true curve",
      "Residual noise SD used to generate observations"
    )
  )
}

write_example_files <- function(sim_data) {
  # Keep all the single-example inputs in CSV form so the plotting and fitting
  # scripts can be rerun without rebuilding everything by hand.
  write_csv(
    sim_data %>% select(ID, Site_name, Start_date, End_date, True_date, Value),
    here("Simulations", "Sim_GP", "data", "simulated_data.csv")
  )
  write_csv(
    ground_truth_grid(),
    here("Simulations", "Sim_GP", "data", "ground_truth.csv")
  )
  write_csv(
    generating_parameters(),
    here("Simulations", "Sim_GP", "data", "generating_parameters.csv")
  )
  write_csv(
    evaluation_targets(),
    here("Simulations", "Sim_GP", "data", "evaluation_targets.csv")
  )
}
