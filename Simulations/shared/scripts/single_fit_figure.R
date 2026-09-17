# One dataset fitted by two or three models: one trend panel per model, then the
# posterior distributions of baseline, slope and sigma. Used by each case's
# 05_single_fit.R.

library(ggplot2)
library(patchwork)

MODEL_TITLES <- c(midpoint = "Midpoint", median = "Calibrated median",
                  marginal = "Full distribution")

#' Fit the models (midpoint and full distribution by default) to one dataset
fit_single_example <- function(dates, y, time_ref_min, time_ref_range, x_pred,
                               models = c("midpoint", "marginal"),
                               model_dir = here::here("Simulations", "shared", "models"),
                               seed = 1, iter_warmup = 1000, iter_sampling = 1000) {
  compiled <- compile_models(model_dir, models)
  fits <- lapply(models, function(m) {
    sd <- model_stan_data(m, dates, y, time_ref_min, time_ref_range, x_pred)
    compiled[[MODEL_FILE[m]]]$sample(
      data = sd, chains = 4, parallel_chains = 4,
      iter_warmup = iter_warmup, iter_sampling = iter_sampling,
      seed = seed, adapt_delta = 0.95, max_treedepth = 12,
      refresh = 0, show_messages = FALSE, show_exceptions = FALSE)
  })
  names(fits) <- models
  fits
}

#' Trend panels: fitted trend with 50% and 90% bands, true trend dashed.
#' obs: data.frame(x, y) drawn as points (e.g. window midpoints).
#' ranges: data.frame(xmin, xmax, y) drawn as horizontal bars (e.g. dating ranges).
trend_panel <- function(fit, x_pred, truth, title, obs = NULL, ranges = NULL,
                        obs_label = "midpoint of dating range",
                        range_label = "dating range", obs_shape = 16) {
  mu    <- fit$draws("mu_pred", format = "draws_matrix")
  trend <- data.frame(x = x_pred,
                      median = apply(mu, 2, median),
                      lo90 = apply(mu, 2, quantile, 0.05), hi90 = apply(mu, 2, quantile, 0.95),
                      lo50 = apply(mu, 2, quantile, 0.25), hi50 = apply(mu, 2, quantile, 0.75))
  ref   <- data.frame(x = x_pred, y = truth$intercept + truth$slope * x_pred)
  p <- ggplot(trend)
  if (!is.null(ranges))
    p <- p + geom_segment(data = ranges,
                          aes(x = xmin, xend = xmax, y = y, yend = y, colour = range_label),
                          linewidth = 0.3, alpha = 0.5, inherit.aes = FALSE)
  p <- p +
    geom_ribbon(aes(x, ymin = lo90, ymax = hi90, fill = "90% interval"), alpha = 0.7) +
    geom_ribbon(aes(x, ymin = lo50, ymax = hi50, fill = "50% interval"), alpha = 0.7) +
    geom_line(aes(x, median, linetype = "posterior median"), linewidth = 0.6) +
    geom_line(data = ref, aes(x, y, linetype = "true trend"), linewidth = 0.7)
  if (!is.null(obs))
    p <- p + geom_point(data = obs, aes(x, y, shape = obs_label), colour = "grey30",
                        size = 1.3, alpha = 0.7, inherit.aes = FALSE)
  p + scale_fill_manual(NULL, values = c("50% interval" = "grey55", "90% interval" = "grey80"),
                        breaks = c("50% interval", "90% interval")) +
    scale_linetype_manual(NULL, values = c("posterior median" = "solid", "true trend" = "dashed")) +
    scale_shape_manual(NULL, values = setNames(obs_shape, obs_label)) +
    scale_colour_manual(NULL, values = setNames("grey35", range_label)) +
    labs(title = title, x = "Calendar year", y = "Value") +
    theme_classic(base_size = 11) +
    theme(plot.title = element_text(size = 11, face = "bold"),
          legend.position = "bottom")
}

#' Posterior histograms, 50% / 90% / tail shaded, true value dashed
posterior_panel <- function(fits, truth, title) {
  one <- function(fit, label) {
    draws <- fit$draws(format = "df")
    long <- data.frame(Baseline = draws$baseline_original,
                       Slope = draws$slope_original, Sigma = draws$sigma)
    long <- utils::stack(long)
    names(long) <- c("Value", "Parameter")
    long$Region <- with(long, ave(Value, Parameter, FUN = function(v) {
      lo50 <- quantile(v, 0.25); hi50 <- quantile(v, 0.75)
      lo90 <- quantile(v, 0.05); hi90 <- quantile(v, 0.95)
      ifelse(v >= lo50 & v <= hi50, "50% CI",
            ifelse(v >= lo90 & v <= hi90, "90% CI", "Tail"))
    }))
    long$Model <- label
    long
  }
  both <- do.call(rbind, lapply(names(fits), function(m) one(fits[[m]], MODEL_TITLES[m])))
  both$Parameter <- factor(both$Parameter, c("Baseline", "Slope", "Sigma"))
  both$Model     <- factor(both$Model, MODEL_TITLES[names(fits)])
  both$Region    <- factor(both$Region, c("50% CI", "90% CI", "Tail"))

  true_lines <- data.frame(
    Parameter = factor(c("Baseline", "Slope", "Sigma"), c("Baseline", "Slope", "Sigma")),
    True = c(truth$intercept, truth$slope, truth$sigma))

  ggplot(both, aes(Value, fill = Region)) +
    geom_histogram(colour = "black", linewidth = 0.2, bins = 40) +
    geom_vline(data = true_lines, aes(xintercept = True),
               linetype = "dashed", linewidth = 0.6) +
    scale_fill_manual(NULL, values = c("50% CI" = "grey40", "90% CI" = "grey65",
                                       "Tail" = "grey88")) +
    facet_grid(Model ~ Parameter, scales = "free_x") +
    labs(title = title, x = NULL, y = "Count") +
    theme_classic(base_size = 11) +
    theme(plot.title = element_text(size = 11, face = "bold"),
          strip.background = element_blank(),
          strip.text = element_text(face = "bold", size = 10),
          legend.position = "bottom")
}

#' The whole figure, one trend panel per fitted model.
#' obs: data.frame(x, y), points on the midpoint panel (and on the
#'   full-distribution panel if ranges is not given).
#' median_obs: data.frame(x, y), points on the calibrated-median panel.
#' ranges: data.frame(xmin, xmax, y), bars on the full-distribution panel.
single_fit_comparison_figure <- function(fits, x_pred, truth, out_path, title,
                                         obs = NULL, ranges = NULL, median_obs = NULL,
                                         obs_label = "midpoint of dating range",
                                         median_label = "calibrated median",
                                         range_label = "dating range") {
  panels <- list()
  for (i in seq_along(fits)) {
    m <- names(fits)[i]
    panel_title <- paste0(LETTERS[i], "  Trend - ", tolower(MODEL_TITLES[m]))
    if (m == "median")
      panels[[i]] <- trend_panel(fits[[m]], x_pred, truth, panel_title, median_obs,
                                 obs_label = median_label, obs_shape = 15)
    else if (m == "marginal" && !is.null(ranges))
      panels[[i]] <- trend_panel(fits[[m]], x_pred, truth, panel_title,
                                 ranges = ranges, range_label = range_label)
    else
      panels[[i]] <- trend_panel(fits[[m]], x_pred, truth, panel_title, obs,
                                 obs_label = obs_label)
  }
  # Same calendar axis on every trend panel
  xlim <- range(x_pred, obs$x, median_obs$x, ranges$xmin, ranges$xmax)
  for (i in seq_along(panels))
    panels[[i]] <- panels[[i]] + coord_cartesian(xlim = xlim)

  pPost <- posterior_panel(fits, truth,
                           paste0(LETTERS[length(fits) + 1], "  Posterior distributions"))

  top <- patchwork::wrap_plots(panels, nrow = 1) +
    patchwork::plot_layout(guides = "collect") &
    theme(legend.position = "bottom")
  n <- length(fits)
  fig <- patchwork::wrap_elements(top) / pPost +
    patchwork::plot_layout(heights = c(1.1, 0.6 * n)) +
    patchwork::plot_annotation(title = title,
                               theme = theme(plot.title = element_text(size = 13, face = "bold")))
  ggsave(out_path, fig, width = 6 * n, height = 5.5 + 2.5 * n, dpi = 300, bg = "white")
  invisible(fig)
}
