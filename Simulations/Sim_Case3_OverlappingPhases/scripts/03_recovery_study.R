# Case 3 recovery study: every model fitted to every dataset in the design.
#
# Reads data/design.csv and re-declares nothing. The study itself is
# run_recovery() in shared/scripts, the same one every case uses. GRID_STEP is
# the spacing of candidate years the full-distribution model sums over.

library(here)
library(cmdstanr)

source(here("Simulations", "shared", "scripts", "weight_rows.R"))
source(here("Simulations", "shared", "scripts", "fit_models.R"))
source(here("Simulations", "shared", "scripts", "run_recovery.R"))
source(here("Simulations", "Sim_Case3_OverlappingPhases", "scripts", "simulate.R"))

GRID_STEP <- as.integer(Sys.getenv("MARGINAL_GRID_STEP", "5"))

design <- read.csv(here("Simulations", "Sim_Case3_OverlappingPhases", "data", "design.csv"))
limit  <- as.integer(Sys.getenv("RECOVERY_LIMIT", "0"))
if (limit > 0) design <- design[seq_len(min(limit, nrow(design))), ]

prepare_dataset <- function(d) {
  sim <- simulate_overlap(d$N, d$intercept, d$slope, d$sigma, d$K, d$alpha_conc,
                          overlap = d$overlap, assign_p = d$assign_p,
                          period_start = d$period_start, period_end = d$period_end,
                          seed = d$seed)

  grid <- seq(d$time_ref_min, d$time_ref_min + d$time_ref_range, by = GRID_STEP)
  list(sim = sim,
       dates = list(start = sim$Start_date, end = sim$End_date,
                    median = (sim$Start_date + sim$End_date) / 2,
                    grid = grid,
                    weights = uniform_rows(sim$Start_date, sim$End_date, grid)))
}

# Median is not fitted: for a flat window it equals the midpoint.
# Figures can only group by the design columns kept here.
run_recovery(design, prepare_dataset,
             keep = c("dataset_id", "sweep", "N", "slope_condition", "overlap", "assign_p"),
             models = setdiff(MODELS, "median"),
             output_csv = Sys.getenv("RECOVERY_OUT",
                                     here("Simulations", "Sim_Case3_OverlappingPhases",
                                          "output", "recovery_results.csv")),
             x_pred_cols = c("period_start", "period_end"))
