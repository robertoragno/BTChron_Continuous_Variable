# Does the two-group design (well dated or coarsely dated) drive the Case 2
# results? The reference datasets are simulated twice, with the same true dates
# and values:
#   two groups   as in the study: half at 6.25% of the period, half at 20-30%
#   continuous   every find draws its width uniformly between 6.25% and w_max
# w_max is chosen so both versions have the same mean squared width, which sets
# how much the midpoint model flattens the slope.
#
# Writes output/check_continuous_widths.csv (fits) and
# output/check_continuous_widths_table.csv (summary).
# About 800 fits; run it in the BTChron tmux session.

library(here)
library(cmdstanr)

source(here("Simulations", "shared", "scripts", "weight_rows.R"))
source(here("Simulations", "shared", "scripts", "fit_models.R"))
source(here("Simulations", "shared", "scripts", "run_recovery.R"))
source(here("Simulations", "shared", "scripts", "recovery_summary.R"))
source(here("Simulations", "Sim_Case2_Typochronology", "scripts", "simulate.R"))

GRID_STEP <- 5
out_dir   <- here("Simulations", "Sim_Case2_Typochronology", "output")
fits_csv  <- file.path(out_dir, "check_continuous_widths.csv")
table_csv <- file.path(out_dir, "check_continuous_widths_table.csv")

design <- read.csv(here("Simulations", "Sim_Case2_Typochronology", "data",
                        "design.csv"))

# The reference setting: N = 200, random slope, everything else at its default
reference <- design[design$sweep == "core" & design$N == 200 &
                    design$slope_condition == "random", ]

# Mean squared width of the two-group design, as a share of the period.
# For a uniform draw between a and b, the mean of w^2 is (a^2 + a*b + b^2) / 3.
fine   <- reference$fine_frac[1]
coarse <- reference$coarse_frac[1]
prop   <- reference$prop_coarse_samples[1]
coarse_lo <- 0.8 * coarse
coarse_hi <- 1.2 * coarse
mean_w2 <- (1 - prop) * fine^2 +
  prop * (coarse_lo^2 + coarse_lo * coarse_hi + coarse_hi^2) / 3

# w_max: solve (fine^2 + fine*w_max + w_max^2) / 3 = mean_w2 for w_max
w_max <- (-fine + sqrt(12 * mean_w2 - 3 * fine^2)) / 2
cat(sprintf("continuous widths: uniform between %.4f and %.4f of the period\n",
            fine, w_max))

# Check on a large sample that both versions have the same mean squared width
two_groups <- simulate_typo(100000, 8, 0.02, 1, prop_coarse_samples = prop,
                            coarse_frac = coarse, fine_frac = fine, seed = 1)
span   <- reference$period_end[1] - reference$period_start[1]
w_two  <- (two_groups$End_date - two_groups$Start_date) / span
w_cont <- runif(100000, fine, w_max)
cat(sprintf("mean squared width: two groups %.5f, continuous %.5f\n",
            mean(w_two^2), mean(w_cont^2)))
stopifnot(abs(mean(w_two^2) / mean(w_cont^2) - 1) < 0.01)

# One copy of the reference datasets for each version
design_two <- reference
design_two$widths <- "two groups"
design_cont <- reference
design_cont$widths <- "continuous"
check_design <- rbind(design_two, design_cont)

# run_recovery() calls this once per dataset
prepare_dataset <- function(d) {
  sim <- simulate_typo(d$N, d$intercept, d$slope, d$sigma,
                       prop_coarse_samples = d$prop_coarse_samples,
                       coarse_frac = d$coarse_frac, fine_frac = d$fine_frac,
                       period_start = d$period_start, period_end = d$period_end,
                       seed = d$seed)

  # Continuous version: same true dates and values, new widths and positions
  if (d$widths == "continuous") {
    set.seed(d$seed + 100000)
    width <- runif(d$N, d$fine_frac, w_max) * (d$period_end - d$period_start)
    sim$Start_date <- sim$True_date - runif(d$N) * width
    sim$End_date   <- sim$Start_date + width
    sim$Coarse     <- NA
  }

  grid <- seq(d$time_ref_min, d$time_ref_min + d$time_ref_range, by = GRID_STEP)
  list(sim = sim,
       dates = list(start = sim$Start_date, end = sim$End_date,
                    median = (sim$Start_date + sim$End_date) / 2,
                    grid = grid,
                    weights = uniform_rows(sim$Start_date, sim$End_date, grid)))
}

run_recovery(check_design, prepare_dataset,
             keep = c("dataset_id", "widths", "N", "intercept", "slope", "sigma",
                      "period_start", "period_end"),
             models = c("midpoint", "marginal"),
             output_csv = fits_csv,
             x_pred_cols = c("period_start", "period_end"))

r <- read_recovery(fits_csv, merge_median = TRUE)
by_widths <- summarise_models(r, c("widths", "model"))
write.csv(by_widths, table_csv, row.names = FALSE)

cat("\n")
print(by_widths[, c("widths", "model", "n_fits", "slope_calib", "slope_calib_lo",
                    "slope_calib_hi", "slope_coverage", "slope_rmse",
                    "slope_width90", "sigma_bias", "sigma_coverage")],
      digits = 3, row.names = FALSE)
