# The period model against the whole Case 1 design, alongside the three the
# paper reports.
#
# 06_period_prior_check.R showed that an estimated study period recovers the
# control window, but only there and only against a design with no zero-slope
# cells. Two things were left open and this settles both: how the period model
# behaves on the plateau, where the models actually separate, and whether an
# over-tight period inflates slopes that are really zero.
#
# EXPERIMENT, not part of the paper's pipeline. It writes to
# output/period_full/ and nowhere else, so the reported tables and figures are
# untouched and no number here can be mistaken for one of those. It draws no
# figures at all, for the same reason.
#
# The four models are refitted from scratch rather than the period fits being
# joined onto the existing three. Same design and same seeds would make a join
# valid, but a join is a place to get it wrong silently, and the cost of not
# doing one is compute rather than correctness.
#
# End state: either the period model is promoted into MODELS in fit_models.R and
# this file and 06_ are deleted, or the approach is dropped and they are deleted
# anyway. It is not meant to survive the decision.
#
#   RECOVERY_LIMIT=20 Rscript .../07_period_full_check.R

library(here)
library(cmdstanr)
suppressMessages(library(rcarbon))

source(here("Simulations", "shared", "scripts", "weight_rows.R"))
source(here("Simulations", "shared", "scripts", "fit_models.R"))
source(here("Simulations", "shared", "scripts", "run_recovery.R"))
source(here("Simulations", "shared", "scripts", "recovery_summary.R"))
source(here("Simulations", "Sim_Case1_Radiocarbon", "scripts", "simulate.R"))

out_dir <- here("Simulations", "Sim_Case1_Radiocarbon", "output", "period_full")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

design <- read.csv(here("Simulations", "Sim_Case1_Radiocarbon", "data",
                        "design.csv"))
limit  <- as.integer(Sys.getenv("RECOVERY_LIMIT", "0"))
if (limit > 0) design <- design[seq_len(min(limit, nrow(design))), ]

# Identical to 03_recovery_study.R's, so the only difference between the two
# runs is the fourth model.
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
             models = c(MODELS, "period"),
             output_csv = file.path(out_dir, "recovery_results.csv"),
             x_pred_cols = c("window_start", "window_end"))

# Tables only. The headline figures are built by 04_ from the reported run and
# are deliberately not rebuilt here.
r <- read_recovery(file.path(out_dir, "recovery_results.csv"))
r$window <- factor(unname(c(plateau = "plateau", steep = "control")[r$window]),
                   levels = c("plateau", "control"))
reference <- r[r$growth_ratio == 1, ]

by_model <- summarise_models(reference, c("model", "window"))
write.csv(by_model, file.path(out_dir, "recovery_table.csv"), row.names = FALSE)
print(by_model, row.names = FALSE)

cat("\nfalse-positive rate (zero-slope cells, 90% interval excludes 0; nominal 0.10):\n")
print(false_positive_rate(reference))

cat("\nconvergence:", sum(r$max_rhat > 1.01, na.rm = TRUE), "fits with rhat > 1.01,",
    sum(r$max_rhat > 1.05, na.rm = TRUE), "above 1.05\n")
cat(sprintf("%.1f%% of fits had a divergent transition\n",
            100 * mean(r$n_divergent > 0)))
