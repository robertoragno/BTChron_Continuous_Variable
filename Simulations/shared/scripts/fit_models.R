# The comparison models, and the code that fits one of them to one dataset.
#
# Dispatch is on MODEL, never on case. A case reaches this file as a `dates`
# frame with the same columns whatever produced it, so nothing here knows
# whether it is looking at radiocarbon or pottery:
#
#   start, end   the interval the midpoint model uses
#   median       the date column the median model uses
#   weights      ragged rows from weight_rows.R, for the marginalised model
#
# Anything case-specific that appears below is a design smell.

# The models the paper reports.
#
# "latent" is deliberately not among them. It samples a date per observation
# inside one [start, end], which for a calibrated posterior means first
# flattening a multimodal curve into a single interval - a reduction nobody
# performs in practice, so reporting it would be arguing with a straw man. And
# for a flat window it is provably the same estimator as "marginal", which
# 00_check_marginal.R demonstrates to five decimal places, so it would add a
# redundant column in Cases 2-4.
#
# latent_date.stan stays in models/ regardless: that equivalence is what proves
# the marginal model correct, and the check needs something to check against.
MODELS <- c("midpoint", "median", "marginal")

#' Posterior summary for one parameter, in the vocabulary the earlier studies
#' used: ACCURACY is whether the interval contained the truth, PRECISION is how
#' wide it was. Recorded per fit; averaged in recovery_summary.R.
summarise_parameter <- function(draws, true_value, name) {
  q <- quantile(draws, c(0.05, 0.25, 0.5, 0.75, 0.95), names = FALSE)
  out <- list(q[3],
              true_value >= q[2] && true_value <= q[4],
              true_value >= q[1] && true_value <= q[5],
              q[5] - q[1],
              q[3] - true_value)
  names(out) <- paste0(name, c("_med", "_cov50", "_cov90", "_width90", "_err"))
  out
}

#' Stan data for one model.
#'
#' Median reuses midpoint.stan with start = end = median: the file takes
#' (start + end) / 2, which is then the median itself. One Stan file, two
#' point-date models, no duplicated code.
model_stan_data <- function(model, dates, y, time_ref_min, time_ref_range,
                            x_pred) {
  common <- list(N = length(y), y = y,
                 time_ref_min = time_ref_min, time_ref_range = time_ref_range,
                 N_pred = length(x_pred), x_pred = x_pred)
  switch(model,
    midpoint = c(common, list(start_date = dates$start, end_date = dates$end)),
    median   = c(common, list(start_date = dates$median, end_date = dates$median)),
    marginal = c(common, list(n_years = length(dates$grid),
                              grid_year = dates$grid,
                              n_weights = dates$weights$n_weights,
                              log_year_prob_packed = dates$weights$log_year_prob_packed,
                              row_first_year = dates$weights$row_first_year,
                              row_n_years = dates$weights$row_n_years)),
    stop("unknown model: ", model))
}

MODEL_FILE <- c(midpoint = "midpoint.stan", median = "midpoint.stan",
                marginal = "marginal_date.stan")

#' Compile the models once, before forking. cmdstanr compiles to a binary beside
#' the .stan file, so letting workers do it races them onto the same path.
compile_models <- function(model_dir, models = MODELS) {
  files  <- unique(MODEL_FILE[models])
  models <- lapply(files, function(f)
    cmdstanr::cmdstan_model(file.path(model_dir, f)))
  names(models) <- files
  models
}

#' Fit one model to one dataset, and return a one-row data.frame.
#'
#' Failures are captured rather than thrown: one pathological dataset should not
#' take down a run of several thousand fits. The error text is kept so a run can
#' be audited afterwards instead of quietly returning fewer rows.
fit_model <- function(model, compiled, stan_data, truth,
                      iter_warmup = 500, iter_sampling = 500) {
  tryCatch({
    fit <- compiled[[MODEL_FILE[model]]]$sample(
      data = stan_data, chains = 4, parallel_chains = 1,
      iter_warmup = iter_warmup, iter_sampling = iter_sampling,
      adapt_delta = 0.95, max_treedepth = 10,
      refresh = 0, show_messages = FALSE, show_exceptions = FALSE)

    pars  <- c("slope_original", "baseline_original", "sigma")
    draws <- fit$draws(variables = pars, format = "draws_matrix")
    diag  <- fit$diagnostic_summary(quiet = TRUE)

    data.frame(c(
      summarise_parameter(draws[, "baseline_original"], truth$intercept, "intercept"),
      summarise_parameter(draws[, "slope_original"],    truth$slope,     "slope"),
      summarise_parameter(draws[, "sigma"],             truth$sigma,     "sigma"),
      list(n_divergent = sum(diag$num_divergent),
           max_rhat    = max(fit$summary(pars)$rhat, na.rm = TRUE),
           error       = NA_character_)))
  }, error = function(e) data.frame(error = conditionMessage(e)))
}

#' Stack one-row results whose columns may differ, because failed fits return
#' only an error column. rbind() would refuse; this pads instead.
bind_rows_padded <- function(rows) {
  cols <- unique(unlist(lapply(rows, names)))
  do.call(rbind, lapply(rows, function(r) {
    r[setdiff(cols, names(r))] <- NA
    r[cols]
  }))
}
