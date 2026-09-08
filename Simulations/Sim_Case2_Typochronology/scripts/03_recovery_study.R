# Case 2 recovery study: every model fitted to every dataset in the design.
#
# Reads data/design.csv and re-declares nothing. The study itself is
# run_recovery() in shared/scripts, the same one Cases 1 and 4 use - the only
# difference here is that a window produces a flat weight row instead of a
# calibrated one.
#
# Note the two meanings of "grid". The design's `grid` column is the 25-yr
# snapping of window boundaries; GRID_STEP below is the spacing of candidate
# years the marginalised model sums over. They are unrelated.

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
  widths <- typo_widths(d$period_start, d$period_end, d$min_frac, d$max_frac,
                        d$grid)
  sim <- simulate_typo(d$N, d$intercept, d$slope, d$sigma,
                       d$period_start, d$period_end, widths,
                       width_prob = typo_width_prob(widths, d$width_mix,
                                                    d$strength),
                       skew_shape = d$skew_shape, grid = d$grid, seed = d$seed)

  grid <- seq(d$time_ref_min, d$time_ref_min + d$time_ref_range, by = GRID_STEP)
  list(sim = sim,
       dates = list(start = sim$Start_date, end = sim$End_date,
                    median = (sim$Start_date + sim$End_date) / 2,
                    grid = grid,
                    weights = uniform_rows(sim$Start_date, sim$End_date, grid)))
}

# Every design column a figure groups by has to travel into the results file.
# max_frac was missing here once: the overhang sweep ran, wrote its fits, and
# 04_ skipped the whole figure without an error because the column it grouped by
# did not exist. Silent, so add to this list whenever the design gains a factor.
run_recovery(design, prepare_dataset,
             keep = c("dataset_id", "sweep", "N", "slope_condition",
                      "width_mix", "skew_shape", "max_frac"),
             models = case_models,
             output_csv = Sys.getenv("RECOVERY_OUT",
                                     here("Simulations", "Sim_Case2_Typochronology",
                                          "output", "recovery_results.csv")),
             x_pred_cols = c("period_start", "period_end"))
