# Fitting one model to one dataset. Every case hands over the same `dates` list:
#
#   start, end   the range the midpoint model uses
#   median       the date the median model uses
#   weights      weight rows from weight_rows.R, for the full-distribution model

# The models the paper reports. latent_date.stan is not among them: for a flat
# window it gives the same answer as "marginal", and checks/00_check_marginal.R
# uses it to prove exactly that.
MODELS <- c("midpoint", "median", "marginal")

#' Posterior summary for one parameter: median, error, and for the 50/80/90/95%
#' intervals whether each contains the true value (accuracy) and its width
#' (precision). Averaged across datasets in recovery_summary.R.
summarise_parameter <- function(draws, true_value, name) {
  med <- median(draws)
  out <- list(med, med - true_value)
  names(out) <- paste0(name, c("_med", "_err"))
  for (level in c(50, 80, 90, 95)) {
    q <- quantile(draws, c(1 - level / 100, 1 + level / 100) / 2, names = FALSE)
    out[[paste0(name, "_cov", level)]]   <- true_value >= q[1] && true_value <= q[2]
    out[[paste0(name, "_width", level)]] <- q[2] - q[1]
  }
  out
}

#' Stan data for one model. The median model reuses midpoint.stan with
#' start = end = median, so (start + end) / 2 is the median.
model_stan_data <- function(model, dates, y, time_ref_min, time_ref_range,
                            x_pred) {
  common <- list(N = length(y), y = y,
                 time_ref_min = time_ref_min, time_ref_range = time_ref_range,
                 N_pred = length(x_pred), x_pred = x_pred)
  switch(model,
    midpoint = c(common, list(start_date = dates$start, end_date = dates$end)),
    median   = c(common, list(start_date = dates$median, end_date = dates$median)),
    # period takes exactly the same data as marginal; the prior is internal.
    period   = ,
    marginal = c(common, list(n_years = length(dates$grid),
                              grid_year = dates$grid,
                              n_weights = dates$weights$n_weights,
                              log_year_prob_packed = dates$weights$log_year_prob_packed,
                              row_first_year = dates$weights$row_first_year,
                              row_n_years = dates$weights$row_n_years)),
    stop("unknown model: ", model))
}

# "period" (marginal_date_period.stan) estimates the study period too. It is a
# prototype used only by Case 1's diagnostics/06_period_prior_check.R.
MODEL_FILE <- c(midpoint = "midpoint.stan", median = "midpoint.stan",
                marginal = "marginal_date.stan",
                period = "marginal_date_period.stan")

#' Compile the models once, before the parallel workers start.
compile_models <- function(model_dir, models = MODELS) {
  files  <- unique(MODEL_FILE[models])
  models <- lapply(files, function(f)
    cmdstanr::cmdstan_model(file.path(model_dir, f)))
  names(models) <- files
  models
}

#' Fit one model to one dataset and return one row of results.
#'
#' An error is recorded in the `error` column instead of stopping the whole run.
#' extra_pars: other parameters to record the posterior mean of (diagnostics only).
fit_model <- function(model, compiled, stan_data, truth,
                      iter_warmup = 500, iter_sampling = 500,
                      extra_pars = NULL) {
  tryCatch({
    fit <- compiled[[MODEL_FILE[model]]]$sample(
      data = stan_data, chains = 4, parallel_chains = 1,
      iter_warmup = iter_warmup, iter_sampling = iter_sampling,
      adapt_delta = 0.95, max_treedepth = 10,
      refresh = 0, show_messages = FALSE, show_exceptions = FALSE)

    pars  <- c("slope_original", "baseline_original", "sigma")
    draws <- fit$draws(variables = pars, format = "draws_matrix")
    diag  <- fit$diagnostic_summary(quiet = TRUE)

    extra <- list()
    if (!is.null(extra_pars)) {
      e <- fit$draws(variables = extra_pars, format = "draws_matrix")
      extra <- as.list(colMeans(e))
      names(extra) <- extra_pars
    }

    data.frame(c(extra,
      summarise_parameter(draws[, "baseline_original"], truth$intercept, "intercept"),
      summarise_parameter(draws[, "slope_original"],    truth$slope,     "slope"),
      summarise_parameter(draws[, "sigma"],             truth$sigma,     "sigma"),
      list(n_divergent = sum(diag$num_divergent),
           max_rhat    = max(fit$summary(pars)$rhat, na.rm = TRUE),
           error       = NA_character_)))
  }, error = function(e) data.frame(error = conditionMessage(e)))
}

#' rbind() for results with different columns (a failed fit has only `error`)
bind_rows_padded <- function(rows) {
  cols <- unique(unlist(lapply(rows, names)))
  do.call(rbind, lapply(rows, function(r) {
    r[setdiff(cols, names(r))] <- NA
    r[cols]
  }))
}
