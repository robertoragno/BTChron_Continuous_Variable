# One recovery study: every model fitted to every dataset in a design.
#
# Cases 1, 2 and 4 differ in exactly one thing - how a design row becomes a
# dataset the models can see. That part is the case's own prepare_dataset();
# the job list, the workers, the error handling and the results file are the
# same for all three, so they live here and are written once.
#
# Reduce the work with env vars while developing:
#   RECOVERY_LIMIT=20 RECOVERY_WORKERS=8 Rscript .../03_recovery_study.R

#' Fit every model to every dataset and write one row per fit.
#'
#' @param design      the design grid, already read from data/design.csv
#' @param prepare_dataset  function(d) for one design row, returning
#'                    list(sim = , dates = ) and optionally extra = , a named
#'                    list of per-dataset values to record alongside the fit
#'                    (Case 1 uses it for mass_kept; the others have none)
#' @param keep        design columns carried into the results file. Every column
#'                    a figure groups by has to be in here - see the note in
#'                    Case 2's 03_ script for what happens when one is missing.
#' @param models      which of MODELS to fit
#' @param output_csv  written only if absent, so a re-run after a crash resumes
#'                    rather than discarding hours of fits. Delete it to refit.
run_recovery <- function(design, prepare_dataset, keep, models, output_csv,
                         x_pred_cols,
                         model_dir = here::here("Simulations", "shared", "models"),
                         n_workers = as.integer(Sys.getenv("RECOVERY_WORKERS", "24")),
                         warmup    = as.integer(Sys.getenv("RECOVERY_WARMUP", "500")),
                         sampling  = as.integer(Sys.getenv("RECOVERY_SAMPLING", "500"))) {

  if (file.exists(output_csv)) {
    cat(basename(output_csv), "already exists; delete it to re-run.\n")
    return(invisible(output_csv))
  }
  dir.create(dirname(output_csv), showWarnings = FALSE, recursive = TRUE)

  jobs <- expand.grid(row = seq_len(nrow(design)), model = models,
                      stringsAsFactors = FALSE)
  # Shuffled so prescheduled workers get a similar mix of cheap and expensive
  # jobs; the design is ordered by N, so contiguous blocks would be very uneven.
  set.seed(1); jobs <- jobs[sample(nrow(jobs)), ]

  compiled <- compile_models(model_dir, models)

  run_one_job <- function(i) {
    job <- jobs[i, ]  # jobs is shuffled; i indexes rows, not design order
    d   <- design[job$row, ]

    # Failures are captured, not thrown: one pathological dataset should not
    # take down a run of several thousand fits.
    prepared <- tryCatch(prepare_dataset(d), error = function(e) e)
    if (inherits(prepared, "error"))
      return(cbind(d[, keep], data.frame(model = job$model,
                                         error = conditionMessage(prepared))))

    row <- fit_model(job$model, compiled,
                     model_stan_data(job$model, prepared$dates, prepared$sim$Value,
                                     d$time_ref_min, d$time_ref_range,
                                     as.numeric(unlist(d[, x_pred_cols]))),
                     truth = list(intercept = d$intercept, slope = d$slope,
                                  sigma = d$sigma),
                     iter_warmup = warmup, iter_sampling = sampling)

    out <- cbind(d[, keep], data.frame(model = job$model))
    if (!is.null(prepared$extra))
      out <- cbind(out, as.data.frame(prepared$extra))
    cbind(out, row)
  }

  cat(sprintf("%d fits (%d datasets x %d models) on %d workers\n",
              nrow(jobs), nrow(design), length(models), n_workers))
  t0 <- Sys.time()
  # mc.preschedule = TRUE forks n_workers processes once and splits the jobs
  # between them. The FALSE variant forks a fresh process per job - thousands of
  # them - and on this machine that was killed outright the moment mclapply
  # started, reproducibly, at job counts above roughly fifty while small runs
  # completed fine. Prescheduling loses some load balancing across uneven job
  # costs, which matters little here because the jobs are shuffled above.
  results <- bind_rows_padded(parallel::mclapply(seq_len(nrow(jobs)), run_one_job,
                                                 mc.cores = n_workers,
                                                 mc.preschedule = TRUE))
  write.csv(results, output_csv, row.names = FALSE)
  cat(sprintf("done in %.1f min, wrote %s\n",
              as.numeric(difftime(Sys.time(), t0, units = "mins")), output_csv))

  failed <- sum(!is.na(results$error))
  if (failed > 0) cat(sprintf("WARNING: %d/%d fits errored\n", failed, nrow(results)))
  invisible(results)
}
