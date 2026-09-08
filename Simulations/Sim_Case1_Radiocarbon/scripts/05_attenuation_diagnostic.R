# Why does the date-marginalised model still attenuate in the steep window?
#
# The main study leaves one result unexplained. Over the plateau the models
# separate as expected (attenuation 1.049 / 0.824 / 0.946 for midpoint, median
# and marginalised), but over the steep window all three land on c ~ 0.82-0.84
# and are statistically indistinguishable, with accuracy falling to 0.23 at the
# largest lab error. A correctly specified errors-in-variables model should not
# flatten the slope by a sixth, so either the model is not doing what it claims
# or the simulator is not producing what the model assumes.
#
# Four conditions, same seeds throughout, so the comparisons are paired:
#
#   oracle       fitted on the true dates. No dating error at all, so c must
#                come out at 1. This is a control on everything else in the
#                pipeline - if the oracle attenuates, the cause is in the
#                fitting or scoring code and none of the other three mean
#                anything.
#   baseline     the study as it stands. Expect c ~ 0.84.
#   matched_dgp  the determination is drawn with the curve's own error folded
#                in, so the generator matches what calibrate() assumes.
#                simulate.R notes that the study ignores this. If c moves to 1
#                here, the attenuation is a simulator artefact and is fixed by
#                changing the generator.
#   window_prior the weight rows are cut to the true study window and
#                renormalised, telling the model what it is otherwise never
#                told. If c moves to 1 here, the attenuation comes from the
#                model spreading dates over the whole padded grid, and is a
#                real property of the method rather than a bug.
#
# The two hypotheses are not exclusive and the conditions are not nested, so
# both can move. What matters is which one moves c further.
#
#   ATTEN_REPS=20 Rscript .../05_attenuation_diagnostic.R

library(here)
library(cmdstanr)
library(parallel)
suppressMessages(library(rcarbon))

source(here("Simulations", "shared", "scripts", "weight_rows.R"))
source(here("Simulations", "shared", "scripts", "fit_models.R"))
source(here("Simulations", "shared", "scripts", "recovery_summary.R"))
source(here("Simulations", "Sim_Case1_Radiocarbon", "scripts", "simulate.R"))

output_csv <- Sys.getenv("ATTEN_OUT",
                         here("Simulations", "Sim_Case1_Radiocarbon", "output",
                              "attenuation_diagnostic.csv"))
n_workers <- as.integer(Sys.getenv("RECOVERY_WORKERS", "24"))
n_rep     <- as.integer(Sys.getenv("ATTEN_REPS", "120"))
warmup    <- as.integer(Sys.getenv("RECOVERY_WARMUP", "500"))
sampling  <- as.integer(Sys.getenv("RECOVERY_SAMPLING", "500"))

# Held at the main study's steep-window settings so the baseline condition is
# comparable with what 03_ already produced, rather than a fresh scenario.
WINDOW     <- c(-1600, -1200)
LAB_ERRORS <- c(15, 50)
N          <- 200
GRID_STEP  <- 5
CONDITIONS <- c("oracle", "baseline", "matched_dgp", "window_prior")

# Same padding as the design, so the baseline condition here and the steep cells
# in 03_ see an identical grid.
d1  <- read.csv(here("Simulations", "Sim_Case1_Radiocarbon", "data", "design.csv"))
d1  <- d1[d1$window == "steep", ][1, ]
PAD <- (d1$time_ref_range - (WINDOW[2] - WINDOW[1])) / 2

set.seed(2026)
design <- expand.grid(rep = seq_len(n_rep), lab_error = LAB_ERRORS,
                      condition = CONDITIONS, stringsAsFactors = FALSE)
# Truths depend on (rep, lab_error) only, so the four conditions see the same
# datasets and differ solely in what the model is given.
truth_key <- unique(design[, c("rep", "lab_error")])
truth_key$seed      <- seq_len(nrow(truth_key))
truth_key$intercept <- runif(nrow(truth_key), 2, 15)
truth_key$sigma     <- runif(nrow(truth_key), 0.5, 4)
truth_key$slope     <- runif(nrow(truth_key), -0.03, 0.03)
design <- merge(design, truth_key, by = c("rep", "lab_error"))

#' As simulate_c14(), but the determination carries the calibration curve's own
#' error as well as the lab's. rcarbon calibrates against
#' sd = sqrt(lab^2 + curve^2); drawing with only the lab error makes the
#' calibrated posterior wider than the truth warrants.
simulate_c14_matched <- function(N, intercept, slope, sigma, window, lab_error,
                                 seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  true_date <- round(runif(N, window[1], window[2]))
  cc  <- rcarbon::uncalibrate(1950 - true_date, verbose = FALSE)
  cra <- round(rnorm(N, cc$ccCRA, sqrt(lab_error^2 + cc$ccError^2)))
  data.frame(CRA = cra, Error = lab_error, True_date = true_date,
             Value = round(intercept + slope * true_date + rnorm(N, 0, sigma), 1))
}

#' Cut each weight row to the study window and renormalise: the model is told
#' the dates came from a known 400-yr period, which is the one thing the study
#' deliberately withholds from it.
#'
#' Rows are stored flat: log_year_prob_packed concatenated, row_first_year the grid
#' position each row starts at, row_n_years its length - so the rows are walked
#' with a running offset rather than indexed as a matrix. Years outside the
#' window go to -Inf, which log_sum_exp handles and the model's generated
#' quantities already floor.
apply_window_prior <- function(w, grid, window) {
  off <- 0L
  for (i in seq_along(w$row_n_years)) {
    k   <- w$row_n_years[i]
    at  <- off + seq_len(k)
    yrs <- grid[w$row_first_year[i] + seq_len(k) - 1L]
    lp  <- w$log_year_prob_packed[at]
    lp[yrs < window[1] | yrs > window[2]] <- -Inf
    m <- max(lp)
    if (!is.finite(m))
      stop("window prior emptied a weight row; the date lies wholly outside")
    w$log_year_prob_packed[at] <- lp - (m + log(sum(exp(lp - m))))
    off <- off + k
  }
  w
}

compiled <- compile_models(here("Simulations", "shared", "models"),
                           c("midpoint", "marginal"))

run_one <- function(i) {
  d    <- design[i, ]
  grid <- seq(WINDOW[1] - PAD, WINDOW[2] + PAD, by = GRID_STEP)
  keep <- c("rep", "lab_error", "condition")

  out <- tryCatch({
    sim <- if (d$condition == "matched_dgp") {
      simulate_c14_matched(N, d$intercept, d$slope, d$sigma, WINDOW,
                           d$lab_error, seed = d$seed)
    } else {
      simulate_c14(N, d$intercept, d$slope, d$sigma, WINDOW, d$lab_error,
                   seed = d$seed)
    }

    truth <- list(intercept = d$intercept, slope = d$slope, sigma = d$sigma)
    x_pred <- WINDOW

    if (d$condition == "oracle") {
      # start = end = the true date, so midpoint.stan regresses on the truth
      dates <- list(start = sim$True_date, end = sim$True_date)
      fit_model("midpoint", compiled,
                model_stan_data("midpoint", dates, sim$Value,
                                WINDOW[1] - PAD, diff(WINDOW) + 2 * PAD, x_pred),
                truth = truth, iter_warmup = warmup, iter_sampling = sampling)
    } else {
      x   <- calibrate(sim$CRA, errors = sim$Error, calCurves = "intcal20",
                       calMatrix = TRUE, verbose = FALSE)
      cal <- calmatrix_to_calendar(x)
      s   <- calibrated_summaries(cal, grid)
      w   <- calibrated_rows(cal, grid)
      if (d$condition == "window_prior") w <- apply_window_prior(w, grid, WINDOW)

      dates <- list(start = s$start, end = s$end, median = s$median,
                    grid = grid, weights = w)
      fit_model("marginal", compiled,
                model_stan_data("marginal", dates, sim$Value,
                                WINDOW[1] - PAD, diff(WINDOW) + 2 * PAD, x_pred),
                truth = truth, iter_warmup = warmup, iter_sampling = sampling)
    }
  }, error = function(e) data.frame(error = conditionMessage(e)))

  cbind(d[, keep], out)
}

cat(sprintf("%d fits (%d reps x %d lab errors x %d conditions) on %d workers\n",
            nrow(design), n_rep, length(LAB_ERRORS), length(CONDITIONS),
            n_workers))
t0 <- Sys.time()
res <- bind_rows_padded(mclapply(seq_len(nrow(design)), run_one,
                                 mc.cores = n_workers, mc.preschedule = TRUE))
dir.create(dirname(output_csv), showWarnings = FALSE, recursive = TRUE)
write.csv(res, output_csv, row.names = FALSE)
cat(sprintf("done in %.1f min, wrote %s\n",
            as.numeric(difftime(Sys.time(), t0, units = "mins")), output_csv))

failed <- sum(!is.na(res$error))
if (failed > 0) cat(sprintf("WARNING: %d/%d fits errored\n", failed, nrow(res)))

ok <- res[is.na(res$error), ]
cat("\ncondition     lab_err    n     c  [95% CI]        accuracy  sigma_err\n")
for (cond in CONDITIONS) for (le in LAB_ERRORS) {
  s <- ok[ok$condition == cond & ok$lab_error == le, ]
  if (nrow(s) < 3) next
  a <- attenuation(s)
  cat(sprintf("%-13s %4d %6d  %.3f  [%.3f, %.3f]     %.3f     %+.3f\n",
              cond, le, nrow(s), a["mean"], a["lo"], a["hi"],
              mean(s$slope_cov90), mean(s$sigma_err)))
}
