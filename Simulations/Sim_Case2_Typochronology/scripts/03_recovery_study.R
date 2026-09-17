# Case 2 recovery study: every model fitted to every dataset in the design.
#
# Reads data/design.csv and re-declares nothing. The study itself is
# run_recovery() in shared/scripts, the same one Cases 1 and 4 use - the only
# difference here is that a window produces a flat weight row instead of a
# calibrated one. GRID_STEP is the spacing of candidate years the marginalised
# model sums over.

library(here)
library(cmdstanr)

source(here("Simulations", "shared", "scripts", "weight_rows.R"))
source(here("Simulations", "shared", "scripts", "fit_models.R"))
source(here("Simulations", "shared", "scripts", "run_recovery.R"))
source(here("Simulations", "Sim_Case2_Typochronology", "scripts", "simulate.R"))

GRID_STEP <- as.integer(Sys.getenv("MARGINAL_GRID_STEP", "5"))

design <- read.csv(here("Simulations", "Sim_Case2_Typochronology", "data",
                        "design.csv"))
limit  <- as.integer(Sys.getenv("RECOVERY_LIMIT", "0"))
if (limit > 0) design <- design[seq_len(min(limit, nrow(design))), ]

# Median is dropped here, not because it is uninteresting but because for a flat
# window it is arithmetically the same estimator as midpoint: (start + end) / 2
# IS the median. Fitting it over the whole design would spend a third of the
# compute re-proving algebra and put a duplicate row in every figure.
#
# It is not taken on trust. 00_check_median_identity.R fits both on a small
# sample and asserts the gap is exactly zero; 04_ labels the surviving row
# "Midpoint / Median".
case_models <- setdiff(MODELS, "median")

prepare_dataset <- function(d) {
  sim <- simulate_typo(d$N, d$intercept, d$slope, d$sigma,
                       prop_coarse_samples = d$prop_coarse_samples,
                       coarse_frac = d$coarse_frac, fine_frac = d$fine_frac,
                       precision_trend = d$precision_trend,
                       growth_ratio = d$growth_ratio,
                       period_start = d$period_start, period_end = d$period_end,
                       seed = d$seed)

  grid <- seq(d$time_ref_min, d$time_ref_min + d$time_ref_range, by = GRID_STEP)
  list(sim = sim,
       dates = list(start = sim$Start_date, end = sim$End_date,
                    median = (sim$Start_date + sim$End_date) / 2,
                    grid = grid,
                    weights = uniform_rows(sim$Start_date, sim$End_date, grid)))
}

# Every design column a figure groups by has to travel into the results file,
# or 04_ cannot draw that sweep. Add to this list whenever the design gains a
# factor. intercept, slope and sigma are kept for the noise ratio in 04_.
run_recovery(design, prepare_dataset,
             keep = c("dataset_id", "sweep", "N", "slope_condition",
                      "prop_coarse_samples", "coarse_frac", "growth_ratio",
                      "precision_trend",
                      "intercept", "slope", "sigma", "period_start",
                      "period_end"),
             models = case_models,
             output_csv = Sys.getenv("RECOVERY_OUT",
                                     here("Simulations", "Sim_Case2_Typochronology",
                                          "output", "recovery_results.csv")),
             x_pred_cols = c("period_start", "period_end"))
