# Case 1 recovery study: every model fitted to every dataset in the design.
#
# Reads data/design.csv and re-declares nothing. The study itself is
# run_recovery() in shared/scripts; all that is case-specific is below, namely
# how a design row becomes a set of calibrated dates. One row per fit goes to
# output/recovery_results.csv; 04_recovery_plots.R turns that into the tables
# and figures.

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

#' Everything the models need from one dataset, built once and shared by all of
#' them so they see identical dates. Calibration is the expensive part, so the
#' same CalDates object serves the weight rows and the envelope summaries.
#'
#' mass_kept rides along in `extra`: how much of the calibrated posterior the
#' HPD envelope keeps, which is the loss the point-date models take and the
#' marginalised one does not.
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
