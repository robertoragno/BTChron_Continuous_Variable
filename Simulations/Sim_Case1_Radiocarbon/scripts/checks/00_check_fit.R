# End-to-end check on one dataset: simulate radiocarbon dates, calibrate, and fit
# the full-distribution model. OLS on calibrated medians is printed alongside for
# reference.

library(here)
suppressMessages(library(rcarbon))
library(cmdstanr)

source(here("Simulations", "shared", "scripts", "weight_rows.R"))
source(here("Simulations", "Sim_Case1_Radiocarbon", "scripts", "simulate.R"))

WINDOW    <- c(-800, -400)   # Hallstatt plateau, the hard case
LAB_ERROR <- 25
N         <- 50
STEP      <- as.integer(Sys.getenv("MARGINAL_GRID_STEP", "5"))
INTERCEPT <- 8
SLOPE     <- 0.02
SIGMA     <- 2

pad  <- c14_grid_pad(WINDOW, LAB_ERROR)
grid <- seq(WINDOW[1] - pad, WINDOW[2] + pad, by = STEP)
cat(sprintf("window %d..%d, pad %d yr, grid %d..%d (%d years at %d yr steps)\n",
            WINDOW[1], WINDOW[2], pad, min(grid), max(grid), length(grid), STEP))

sim <- simulate_c14(N, INTERCEPT, SLOPE, SIGMA, WINDOW, LAB_ERROR, seed = 1)

x   <- calibrate(sim$CRA, errors = sim$Error, calCurves = "intcal20",
                 calMatrix = TRUE, verbose = FALSE)
cal <- calmatrix_to_calendar(x)
w   <- calibrated_rows(cal, grid)

cat(sprintf("support per date: median %d, max %d grid cells (%d..%d yr)\n",
            median(w$row_n_years), max(w$row_n_years), median(w$row_n_years) * STEP,
            max(w$row_n_years) * STEP))

time_ref_min   <- min(grid)
time_ref_range <- max(grid) - min(grid)
x_pred <- seq(WINDOW[1], WINDOW[2], length.out = 50)

fit <- cmdstan_model(here("Simulations", "shared", "models",
                          "marginal_date.stan"))$sample(
  data = list(N = N, y = sim$Value,
              n_years = length(grid), grid_year = grid,
              n_weights = w$n_weights, log_year_prob_packed = w$log_year_prob_packed, row_first_year = w$row_first_year, row_n_years = w$row_n_years,
              time_ref_min = time_ref_min, time_ref_range = time_ref_range,
              N_pred = length(x_pred), x_pred = x_pred),
  chains = 4, parallel_chains = 4,
  iter_warmup = 1000, iter_sampling = 1000, seed = 1, refresh = 0)

print(fit$summary(c("slope_original", "baseline_original", "sigma")))
cat(sprintf("\ntrue: slope %.3f, intercept %.1f, sigma %.1f\n",
            SLOPE, INTERCEPT, SIGMA))

# Point-date reference: OLS on the calibrated median, the usual practice.
med_bcad <- 1950 - summary(x)$MedianBP
cat(sprintf("OLS on calibrated medians: slope %.5f\n",
            coef(lm(sim$Value ~ med_bcad))[2]))

# Recovered dates against the truth, as a sanity check on the weight rows.
dates <- colMeans(fit$draws("date_actual", format = "draws_matrix"))
cat(sprintf("date posterior means vs true: mean error %.0f yr, correlation %.3f\n",
            mean(dates - sim$True_date), cor(dates, sim$True_date)))
