# One dataset fitted by the midpoint and full-distribution models: the fitted
# trends (panels A, B) and the posterior distributions of baseline, slope and
# sigma (panel C). Used by each case's 05_single_fit.R.

library(ggplot2)
library(patchwork)

#' Fit the midpoint and full-distribution models to one dataset
fit_single_example <- function(dates, y, time_ref_min, time_ref_range, x_pred,
                               model_dir = here::here("Simulations", "shared", "models"),
                               seed = 1, iter_warmup = 1000, iter_sampling = 1000) {
  compiled <- compile_models(model_dir, c("midpoint", "marginal"))
  fits <- lapply(c("midpoint", "marginal"), function(m) {
    sd <- model_stan_data(m, dates, y, time_ref_min, time_ref_range, x_pred)
    compiled[[MODEL_FILE[m]]]$sample(
      data = sd, chains = 4, parallel_chains = 4,
      iter_warmup = iter_warmup, iter_sampling = iter_sampling,
      seed = seed, adapt_delta = 0.95, max_treedepth = 12,
      refresh = 0, show_messages = FALSE, show_exceptions = FALSE)
  })
  names(fits) <- c("midpoint", "marginal")
  fits
}

#' Panels A and B: fitted trend with 50% and 90% bands, true trend dashed
trend_panel <- function(fit, x_pred, truth, title, obs = NULL) {
  mu    <- fit$draws("mu_pred", format = "draws_matrix")
  trend <- data.frame(x = x_pred,
                      median = apply(mu, 2, median),
                      lo90 = apply(mu, 2, quantile, 0.05), hi90 = apply(mu, 2, quantile, 0.95),
                      lo50 = apply(mu, 2, quantile, 0.25), hi50 = apply(mu, 2, quantile, 0.75))
  ref   <- data.frame(x = x_pred, y = truth$intercept + truth$slope * x_pred)
  p <- ggplot(trend) +
    geom_ribbon(aes(x, ymin = lo90, ymax = hi90), fill = "grey80") +
    geom_ribbon(aes(x, ymin = lo50, ymax = hi50), fill = "grey60") +
    geom_line(aes(x, median), linewidth = 0.6, colour = "black") +
    geom_line(data = ref, aes(x, y), linewidth = 0.7, colour = "black",
              linetype = "dashed")
  if (!is.null(obs))
    p <- p + geom_point(data = obs, aes(x, y), colour = "grey40", size = 1.1,
                        alpha = 0.6, inherit.aes = FALSE)
  p + labs(title = title, x = "Calendar year", y = "Value") +
    theme_classic(base_size = 11) +
    theme(plot.title = element_text(size = 11, face = "bold"))
}

#' Panel C: posterior histograms, 50% / 90% / tail shaded, true value dashed
posterior_panel <- function(fits, truth) {
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
  both <- rbind(one(fits$midpoint, "Midpoint"), one(fits$marginal, "Full distribution"))
  both$Parameter <- factor(both$Parameter, c("Baseline", "Slope", "Sigma"))
  both$Model     <- factor(both$Model, c("Midpoint", "Full distribution"))
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
    labs(title = "C  Posterior distributions", x = NULL, y = "Count") +
    theme_classic(base_size = 11) +
    theme(plot.title = element_text(size = 11, face = "bold"),
          strip.background = element_blank(),
          strip.text = element_text(face = "bold", size = 10),
          legend.position = "bottom")
}

#' The whole figure. obs, if given, is data.frame(x, y) of points to draw on the
#' trend panels (e.g. window midpoints).
single_fit_comparison_figure <- function(fits, x_pred, truth, out_path, title,
                                         obs = NULL) {
  pA <- trend_panel(fits$midpoint, x_pred, truth, "A  Trend - midpoint", obs)
  pB <- trend_panel(fits$marginal, x_pred, truth, "B  Trend - full distribution", obs)
  pC <- posterior_panel(fits, truth)
  fig <- (pA | pB) / pC +
    patchwork::plot_layout(heights = c(1, 1.2)) +
    patchwork::plot_annotation(title = title,
                               theme = theme(plot.title = element_text(size = 13, face = "bold")))
  ggsave(out_path, fig, width = 12, height = 10, dpi = 300, bg = "white")
  invisible(fig)
}
