# The recovery study loop shared by every case: fit every model to every dataset
# in parallel and write one row per fit. Each case supplies prepare_dataset().
#
# For a quick test run:
#   RECOVERY_LIMIT=20 RECOVERY_WORKERS=8 Rscript .../03_recovery_study.R

#' Fit every model to every dataset and write one row per fit.
#'
#' @param design      the design grid, already read from data/design.csv
#' @param prepare_dataset  function(d) turning one design row into
#'                    list(sim, dates) and optionally extra (values to record)
#' @param keep        design columns copied into the results; figures can only
#'                    group by columns listed here
#' @param models      which of MODELS to fit
#' @param output_csv  if it already exists nothing is fitted; delete it to refit
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
  # Shuffle so each worker gets a mix of small and large datasets
  set.seed(1); jobs <- jobs[sample(nrow(jobs)), ]

  compiled <- compile_models(model_dir, models)

  run_one_job <- function(i) {
    job <- jobs[i, ]  # jobs is shuffled; i indexes rows, not design order
    d   <- design[job$row, ]

    # Record a failure instead of stopping the run
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
  # mc.preschedule = TRUE: one process per worker. FALSE (one process per job)
  # was killed on this server for runs of more than about fifty jobs.
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
