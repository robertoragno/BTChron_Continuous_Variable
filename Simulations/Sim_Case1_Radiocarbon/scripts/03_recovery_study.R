# Case 1 recovery study: every model fitted to every dataset in the design.
# Reads data/design.csv; writes one row per fit to output/recovery_results.csv.
# The fitting loop is run_recovery() in shared/scripts.

library(here)
library(cmdstanr)

source(here("Simulations", "shared", "scripts", "weight_rows.R"))
source(here("Simulations", "shared", "scripts", "fit_models.R"))
source(here("Simulations", "shared", "scripts", "run_recovery.R"))
source(here("Simulations", "Sim_Case1_Radiocarbon", "scripts", "simulate.R"))
suppressMessages(library(rcarbon))

design <- read.csv(here("Simulations", "Sim_Case1_Radiocarbon", "data",
                        "design.csv"))
limit  <- as.integer(Sys.getenv("RECOVERY_LIMIT", "0"))  # 0 = whole design
if (limit > 0) design <- design[seq_len(min(limit, nrow(design))), ]

# One dataset: simulate, calibrate, and prepare what each model reads.
#   midpoint           centre of the 95% HPD range (start, end)
#   median             calibrated median
#   full distribution  the whole calibrated distribution on the grid (weights)
# mass_kept is the share of each calibrated distribution inside its 95% range.
prepare_dataset <- function(d) {
  sim <- simulate_c14(d$N, d$intercept, d$slope, d$sigma,
                      c(d$window_start, d$window_end), d$lab_error,
                      growth_ratio = d$growth_ratio, seed = d$seed)

  grid <- seq(d$time_ref_min, d$time_ref_min + d$time_ref_range, by = d$grid_step)
  x    <- calibrate(sim$CRA, errors = sim$Error, calCurves = "intcal20",
                    calMatrix = TRUE, verbose = FALSE)
  cal  <- calmatrix_to_calendar(x)

  s <- calibrated_summaries(cal, grid)
  list(sim = sim,
       dates = list(start = s$start, end = s$end, median = s$median,
                    grid = grid, weights = calibrated_rows(cal, grid)),
       extra = list(mass_kept = mean(s$mass_kept)))
}

run_recovery(design, prepare_dataset,
             keep = c("dataset_id", "sweep", "N", "slope_condition",
                      "lab_error", "window", "growth_ratio"),
             models = MODELS,
             output_csv = Sys.getenv("RECOVERY_OUT",
                                     here("Simulations", "Sim_Case1_Radiocarbon",
                                          "output", "recovery_results.csv")),
             x_pred_cols = c("window_start", "window_end"))
