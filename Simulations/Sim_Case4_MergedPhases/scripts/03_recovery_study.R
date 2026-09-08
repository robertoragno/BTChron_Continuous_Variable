# Case 4 recovery study: every model fitted to every dataset in the design.
#
# Reads data/design.csv and re-declares nothing. The study itself is
# run_recovery() in shared/scripts, the same one Cases 1 and 2 use; what is
# specific to Case 4 is that a find's recorded window spans a run of merged fine
# phases rather than a single one.
#
# Median is not fitted. For a flat window (start + end) / 2 is the median, so it
# would duplicate the midpoint row; Case 2's 00_check_median_identity.R proves
# that once for every flat-window case. 04_ labels the surviving row
# "Midpoint / Median".

library(here)
library(cmdstanr)

source(here("Simulations", "shared", "scripts", "weight_rows.R"))
source(here("Simulations", "shared", "scripts", "fit_models.R"))
source(here("Simulations", "shared", "scripts", "run_recovery.R"))
source(here("Simulations", "shared", "scripts", "partition.R"))
source(here("Simulations", "Sim_Case4_MergedPhases", "scripts", "simulate.R"))

GRID_STEP <- as.integer(Sys.getenv("MARGINAL_GRID_STEP", "5"))

design <- read.csv(here("Simulations", "Sim_Case4_MergedPhases", "data",
                        "design.csv"))
limit  <- as.integer(Sys.getenv("RECOVERY_LIMIT", "0"))
if (limit > 0) design <- design[seq_len(min(limit, nrow(design))), ]

prepare_dataset <- function(d) {
  sim  <- simulate_merged(d$N, d$intercept, d$slope, d$sigma, d$K, d$alpha_conc,
                          merge_max = d$merge_max, t_min = d$t_min,
                          t_max = d$t_max, seed = d$seed)
  grid <- seq(d$time_ref_min, d$time_ref_min + d$time_ref_range, by = GRID_STEP)
  list(sim = sim,
       dates = list(start = sim$Start_date, end = sim$End_date,
                    median = (sim$Start_date + sim$End_date) / 2,
                    grid = grid,
                    weights = uniform_rows(sim$Start_date, sim$End_date, grid)))
}

run_recovery(design, prepare_dataset,
             keep = c("dataset_id", "sweep", "N", "slope_condition", "merge_max"),
             models = setdiff(MODELS, "median"),
             output_csv = Sys.getenv("RECOVERY_OUT",
                                     here("Simulations", "Sim_Case4_MergedPhases",
                                          "output", "recovery_results.csv")),
             x_pred_cols = c("t_min", "t_max"))
