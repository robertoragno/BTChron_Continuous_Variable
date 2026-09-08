# Fits every model to every dataset in the ratio sweep design.
#
# Same shape as each case's 03_recovery_study.R, and the same shared fitting
# code. The third model, "oracle", is local to this study: it is midpoint.stan
# handed the TRUE dates, which makes it a control on the prior rather than a
# method anyone could use. See README.md.
#
#   RECOVERY_LIMIT=20 Rscript .../02_sweep.R

library(here)
library(cmdstanr)
library(parallel)

source(here("Simulations", "shared", "scripts", "weight_rows.R"))
source(here("Simulations", "shared", "scripts", "fit_models.R"))
source(here("Simulations", "Sim_Case2_Typochronology", "scripts", "simulate.R"))

design_csv <- here("Simulations", "Sim_Sweep_Ratio", "data", "design.csv")
output_csv <- Sys.getenv("SWEEP_OUT",
                         here("Simulations", "Sim_Sweep_Ratio", "output",
                              "sweep_results.csv"))
n_workers <- as.integer(Sys.getenv("RECOVERY_WORKERS", "24"))
limit     <- as.integer(Sys.getenv("RECOVERY_LIMIT", "0"))
warmup    <- as.integer(Sys.getenv("RECOVERY_WARMUP", "500"))
sampling  <- as.integer(Sys.getenv("RECOVERY_SAMPLING", "500"))
GRID_STEP <- as.integer(Sys.getenv("MARGINAL_GRID_STEP", "5"))
WIDTHS    <- seq(50, 200, by = 25)   # must match 01_design.R
GRID      <- 25

design <- read.csv(design_csv)
if (limit > 0) design <- design[seq_len(min(limit, nrow(design))), ]

# Median is not fitted: for a flat window it is the same estimator as midpoint,
# which Case 2's 00_check_median_identity.R establishes.
SWEEP_MODELS <- c("oracle", "midpoint", "marginal")

jobs <- expand.grid(row = seq_len(nrow(design)), model = SWEEP_MODELS,
                    stringsAsFactors = FALSE)
set.seed(1); jobs <- jobs[sample(nrow(jobs)), ]

compiled <- compile_models(here("Simulations", "shared", "models"),
                           c("midpoint", "marginal"))

run_one_job <- function(i) {
  job  <- jobs[i, ]
  d    <- design[job$row, ]
  keep <- c("dataset_id", "period_span", "coarseness", "rise_level", "N")

  out <- tryCatch({
    sim <- simulate_typo(d$N, d$intercept, d$slope, d$sigma,
                         d$period_start, d$period_end, WIDTHS,
                         grid = GRID, seed = d$seed)
    truth  <- list(intercept = d$intercept, slope = d$slope, sigma = d$sigma)
    x_pred <- c(d$period_start, d$period_end)

    dates <- if (job$model == "oracle") {
      # start = end = the true date, so midpoint.stan regresses on the truth
      list(start = sim$True_date, end = sim$True_date)
    } else {
      grid <- seq(d$time_ref_min, d$time_ref_min + d$time_ref_range,
                  by = GRID_STEP)
      list(start = sim$Start_date, end = sim$End_date, grid = grid,
           weights = uniform_rows(sim$Start_date, sim$End_date, grid))
    }

    stan_model <- if (job$model == "oracle") "midpoint" else job$model
    fit_model(stan_model, compiled,
              model_stan_data(stan_model, dates, sim$Value,
                              d$time_ref_min, d$time_ref_range, x_pred),
              truth = truth, iter_warmup = warmup, iter_sampling = sampling)
  }, error = function(e) data.frame(error = conditionMessage(e)))

  cbind(d[, keep], data.frame(model = job$model), out)
}

if (file.exists(output_csv)) {
  cat("sweep_results.csv already exists; delete it to re-run.\n")
} else {
  dir.create(dirname(output_csv), showWarnings = FALSE, recursive = TRUE)
  cat(sprintf("%d fits (%d datasets x %d models) on %d workers\n",
              nrow(jobs), nrow(design), length(SWEEP_MODELS), n_workers))
  t0 <- Sys.time()
  # mc.preschedule = TRUE for the reason recorded in the cases' 03_ scripts:
  # the FALSE variant forks a process per job and is killed on this machine.
  results <- bind_rows_padded(mclapply(seq_len(nrow(jobs)), run_one_job,
                                       mc.cores = n_workers,
                                       mc.preschedule = TRUE))
  write.csv(results, output_csv, row.names = FALSE)
  cat(sprintf("done in %.1f min, wrote %s\n",
              as.numeric(difftime(Sys.time(), t0, units = "mins")), output_csv))

  failed <- sum(!is.na(results$error))
  if (failed > 0) cat(sprintf("WARNING: %d/%d fits errored\n", failed, nrow(results)))
}
