# Does an estimated study period account for the control-window attenuation?
#
# The control window is where the calibration curve is steep, posteriors are
# narrow and single-humped, and the models cannot differ on the shape of a date
# - yet all three flatten the slope by about a sixth. 05_attenuation_diagnostic.R
# located the cause: each weight row carries mass outside the period the finds
# came from, nothing ties the dates to each other, the fitted spread of dates
# exceeds the real one, and an inflated var(x) flattens a slope. Clipping rows
# to the true window recovers c = 1, but that number comes from the generator.
#
# This asks whether the period can be estimated instead of supplied.
#
#   baseline  marginal_date.stan, as the study runs it
#   period    marginal_date_period.stan, the same model with a shared
#             Normal(mu, tau) over the dates, both estimated from the data
#
# Paired: the two conditions see identical datasets from identical seeds, so a
# difference between them is the prior and nothing else.
#
# What each outcome means:
#   c moves to 1 and period_sd lands near the truth  -> mechanism confirmed, and
#       a trapezoid is the shape worth fitting next
#   c moves past 1                                   -> over-correction; the
#       period is being estimated too tight, check the zero-slope false-positive
#       rate before going further
#   c does not move                                  -> the diagnosis is wrong
#       and the attenuation is something else
#
#   PERIOD_REPS=20 Rscript .../06_period_prior_check.R

library(here)
library(cmdstanr)
library(parallel)
suppressMessages(library(rcarbon))

source(here("Simulations", "shared", "scripts", "weight_rows.R"))
source(here("Simulations", "shared", "scripts", "fit_models.R"))
source(here("Simulations", "shared", "scripts", "recovery_summary.R"))
source(here("Simulations", "Sim_Case1_Radiocarbon", "scripts", "simulate.R"))

output_csv <- Sys.getenv("PERIOD_OUT",
                         here("Simulations", "Sim_Case1_Radiocarbon", "output",
                              "period_prior_check.csv"))
n_workers <- as.integer(Sys.getenv("RECOVERY_WORKERS", "24"))
n_rep     <- as.integer(Sys.getenv("PERIOD_REPS", "120"))
warmup    <- as.integer(Sys.getenv("RECOVERY_WARMUP", "500"))
sampling  <- as.integer(Sys.getenv("RECOVERY_SAMPLING", "500"))

# The control window at the main study's settings, so baseline here and the
# control cells of the recovery study are the same quantity.
WINDOW     <- c(-1600, -1200)
LAB_ERRORS <- c(15, 50)
N          <- 200
GRID_STEP  <- 5
CONDITIONS <- c("baseline", "period")
# Condition names read for the report; these are the models behind them.
CONDITION_MODEL <- c(baseline = "marginal", period = "period")

# True dates are uniform over the window, so the period a correctly behaving
# model should recover has this standard deviation.
TRUE_PERIOD_SD <- diff(WINDOW) / sqrt(12)

d1  <- read.csv(here("Simulations", "Sim_Case1_Radiocarbon", "data", "design.csv"))
d1  <- d1[d1$window == "steep", ][1, ]
PAD <- (d1$time_ref_range - diff(WINDOW)) / 2

set.seed(2026)
design <- expand.grid(rep = seq_len(n_rep), lab_error = LAB_ERRORS,
                      condition = CONDITIONS, stringsAsFactors = FALSE)

# One truth per (rep, lab_error), shared by both conditions, so the comparison
# is paired rather than merely averaged.
truth_key <- unique(design[, c("rep", "lab_error")])
truth_key$seed      <- seq_len(nrow(truth_key))
truth_key$intercept <- runif(nrow(truth_key), 2, 15)
truth_key$sigma     <- runif(nrow(truth_key), 0.5, 4)
truth_key$slope     <- runif(nrow(truth_key), -0.03, 0.03)
design <- merge(design, truth_key, by = c("rep", "lab_error"))

compiled <- compile_models(here("Simulations", "shared", "models"),
                           c("marginal", "period"))

run_one <- function(i) {
  d    <- design[i, ]
  grid <- seq(WINDOW[1] - PAD, WINDOW[2] + PAD, by = GRID_STEP)
  keep <- c("rep", "lab_error", "condition")

  out <- tryCatch({
    sim <- simulate_c14(N, d$intercept, d$slope, d$sigma, WINDOW, d$lab_error,
                        seed = d$seed)
    cal <- calmatrix_to_calendar(
      calibrate(sim$CRA, errors = sim$Error, calCurves = "intcal20",
                calMatrix = TRUE, verbose = FALSE))
    s     <- calibrated_summaries(cal, grid)
    dates <- list(start = s$start, end = s$end, median = s$median,
                  grid = grid, weights = calibrated_rows(cal, grid))

    model <- CONDITION_MODEL[[d$condition]]
    fit_model(model, compiled,
              model_stan_data(model, dates, sim$Value,
                              WINDOW[1] - PAD, diff(WINDOW) + 2 * PAD, WINDOW),
              truth = list(intercept = d$intercept, slope = d$slope,
                           sigma = d$sigma),
              iter_warmup = warmup, iter_sampling = sampling,
              extra_pars = if (d$condition == "period")
                             c("period_centre_original", "period_sd_original"))
  }, error = function(e) data.frame(error = conditionMessage(e)))

  cbind(d[, keep], out)
}

cat(sprintf("%d fits (%d reps x %d lab errors x %d conditions) on %d workers\n",
            nrow(design), n_rep, length(LAB_ERRORS), length(CONDITIONS), n_workers))
t0  <- Sys.time()
res <- bind_rows_padded(mclapply(seq_len(nrow(design)), run_one,
                                 mc.cores = n_workers, mc.preschedule = TRUE))
dir.create(dirname(output_csv), showWarnings = FALSE, recursive = TRUE)
write.csv(res, output_csv, row.names = FALSE)
cat(sprintf("done in %.1f min, wrote %s\n",
            as.numeric(difftime(Sys.time(), t0, units = "mins")), output_csv))

failed <- sum(!is.na(res$error))
if (failed > 0) cat(sprintf("WARNING: %d/%d fits errored\n", failed, nrow(res)))

ok <- res[is.na(res$error), ]
cat(sprintf("\ntrue period sd over a %d-yr uniform window: %.1f yr\n",
            diff(WINDOW), TRUE_PERIOD_SD))
cat("\ncondition   lab_err     n      c  [95% CI]        accuracy  sigma_err  period_sd\n")
for (cond in CONDITIONS) for (le in LAB_ERRORS) {
  s <- ok[ok$condition == cond & ok$lab_error == le, ]
  if (nrow(s) < 3) next
  a  <- attenuation(s)
  sd <- if ("period_sd_original" %in% names(s)) mean(s$period_sd_original, na.rm = TRUE) else NA
  cat(sprintf("%-11s %4d %6d  %.3f  [%.3f, %.3f]     %.3f     %+.3f    %s\n",
              cond, le, nrow(s), a["mean"], a["lo"], a["hi"],
              mean(s$slope_cov90), mean(s$sigma_err),
              if (is.na(sd)) "-" else sprintf("%.1f", sd)))
}
