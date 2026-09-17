# Recovery metrics, tables and figures shared by every case's 04_ script.
#
#   bias         posterior median minus the true value, averaged over datasets
#   attenuation  how much of the true slope survives (1 = all of it)
#   accuracy     how often the 90% interval contains the true value (target 0.9)
#   precision    how wide that interval is
#
# Greyscale, with dark red for the full-distribution model.

library(ggplot2)

# Metrics carried into recovery_table.csv
METRICS <- c("Slope bias", "Slope attenuation", "Slope accuracy (90%)",
             "Slope precision (90% width)", "Sigma error")

# Model names, colours and shapes used in every figure
MODEL_LEVELS  <- c("midpoint", "median", "marginal")
MODEL_LABELS  <- c(midpoint = "Midpoint", median = "Calibrated median",
                   marginal = "Full distribution")
MODEL_COLOURS <- c(Midpoint = "grey55", `Midpoint / Median` = "grey55",
                   `Calibrated median` = "grey25", `Full distribution` = "#780000")
MODEL_SHAPES  <- c(Midpoint = 17, `Midpoint / Median` = 17, `Calibrated median` = 15,
                   `Full distribution` = 16)

ACCENT <- "#780000"

panel_theme <- theme_classic(base_size = 11) +
  theme(plot.title = element_text(size = 11, face = "bold"),
        strip.background = element_blank(),
        strip.text = element_text(face = "bold"),
        legend.position = "top",
        legend.title = element_blank())

#' Read a results file, drop failed fits and label the models.
#'
#' merge_median: for flat windows the median model is not fitted (it equals the
#' midpoint), so the midpoint row is labelled "Midpoint / Median".
read_recovery <- function(path, merge_median = FALSE) {
  r <- read.csv(path)
  n_failed <- sum(!is.na(r$error))
  if (n_failed > 0)
    cat(sprintf("%d of %d fits errored and are excluded\n", n_failed, nrow(r)))
  r <- r[is.na(r$error), ]
  labels <- MODEL_LABELS
  if (merge_median) {
    if ("median" %in% r$model)
      stop("merge_median = TRUE but median fits are present; ",
           "either drop them or set merge_median = FALSE")
    labels["midpoint"] <- "Midpoint / Median"
  }
  r$model <- factor(labels[r$model], levels = unname(labels[MODEL_LEVELS]))
  droplevels(r)
}

#' Mean and its Monte Carlo standard error (how precisely the mean is known from
#' this many datasets)
mc_mean <- function(x) {
  x <- x[is.finite(x)]
  c(mean = mean(x), se = sd(x) / sqrt(length(x)))
}

#' A proportion (e.g. accuracy) with a 95% Jeffreys interval, which stays
#' inside 0-1 even near 1
jeffreys_prop <- function(x, level = 0.95) {
  x <- x[!is.na(x)]
  n <- length(x)
  k <- sum(x)
  a <- (1 - level) / 2
  c(mean = k / n,
    lo = qbeta(a,     k + 0.5, n - k + 0.5),
    hi = qbeta(1 - a, k + 0.5, n - k + 0.5))
}

#' Attenuation: the slope of estimated slopes regressed on true slopes, with a 95%
#' interval. 1 means the slope is recovered one-for-one, below 1 flattened.
#'
#' Bias alone cannot show this: true slopes are drawn around zero, so flattening
#' positive and negative slopes cancels out in the average error. NA when the
#' true slopes do not vary (zero-slope datasets).
attenuation <- function(d, level = 0.95) {
  truth <- d$slope_med - d$slope_err
  ok    <- is.finite(truth) & is.finite(d$slope_med)
  if (sum(ok) < 3 || sd(truth[ok]) < .Machine$double.eps^0.5)
    return(c(mean = NA_real_, lo = NA_real_, hi = NA_real_))
  fit <- stats::lm(d$slope_med[ok] ~ truth[ok])
  ci  <- stats::confint(fit, 2, level = level)
  c(mean = unname(coef(fit)[2]), lo = ci[1], hi = ci[2])
}

#' One row of metrics per group (e.g. per model, or per model and N)
summarise_models <- function(d, by) {
  stats <- do.call(rbind, lapply(split(seq_len(nrow(d)), d[by], drop = TRUE),
    function(i) {
      s     <- d[i, ]
      bias  <- mc_mean(s$slope_err)
      atten <- attenuation(s)
      acc   <- jeffreys_prop(s$slope_cov90)
      prec  <- mc_mean(s$slope_width90)
      sig   <- mc_mean(s$sigma_err)
      # Intercept and sigma get the same bias / accuracy / precision as the
      # slope. They go to the tables only, not the headline figures.
      int_bias <- mc_mean(s$intercept_err)
      int_acc  <- jeffreys_prop(s$intercept_cov90)
      int_prec <- mc_mean(s$intercept_width90)
      sig_acc  <- jeffreys_prop(s$sigma_cov90)
      sig_prec <- mc_mean(s$sigma_width90)
      data.frame(s[1, by, drop = FALSE], n_fits = nrow(s),
                 intercept_bias         = int_bias["mean"],
                 intercept_bias_se      = int_bias["se"],
                 intercept_accuracy     = int_acc["mean"],
                 intercept_accuracy_lo  = int_acc["lo"],
                 intercept_accuracy_hi  = int_acc["hi"],
                 intercept_precision    = int_prec["mean"],
                 intercept_precision_se = int_prec["se"],
                 sigma_accuracy         = sig_acc["mean"],
                 sigma_accuracy_lo      = sig_acc["lo"],
                 sigma_accuracy_hi      = sig_acc["hi"],
                 sigma_precision        = sig_prec["mean"],
                 sigma_precision_se     = sig_prec["se"],
                 slope_bias         = bias["mean"],
                 slope_bias_se      = bias["se"],
                 slope_atten        = atten["mean"],
                 slope_atten_lo     = atten["lo"],
                 slope_atten_hi     = atten["hi"],
                 slope_accuracy     = acc["mean"],
                 slope_accuracy_lo  = acc["lo"],
                 slope_accuracy_hi  = acc["hi"],
                 slope_precision    = prec["mean"],
                 slope_precision_se = prec["se"],
                 sigma_err          = sig["mean"],
                 sigma_err_se       = sig["se"],
                 row.names = NULL)
    }))
  stats[order(stats$model), ]
}

#' Summary table to long format: one row per group and metric, with bounds
metric_long <- function(stats, keep) {
  long <- do.call(rbind, lapply(seq_len(nrow(stats)), function(i) {
    s <- stats[i, ]
    # Accuracy and attenuation have their own intervals; the rest are +/- 2 SE
    value <- c(s$slope_bias, s$slope_atten, s$slope_accuracy, s$slope_precision,
               s$sigma_err)
    se    <- c(s$slope_bias_se, NA, NA, s$slope_precision_se, s$sigma_err_se)
    lo    <- c(NA, s$slope_atten_lo, s$slope_accuracy_lo, NA, NA)
    hi    <- c(NA, s$slope_atten_hi, s$slope_accuracy_hi, NA, NA)
    data.frame(s[, keep, drop = FALSE], metric = METRICS, value = value,
               lo = ifelse(is.na(se), lo, value - 2 * se),
               hi = ifelse(is.na(se), hi, value + 2 * se),
               row.names = NULL)
  }))
  long$metric <- factor(long$metric, levels = METRICS)
  long
}

# The three metrics shown in the headline figures
HEADLINE_METRICS <- c("Slope attenuation", "Slope accuracy (90%)", "Sigma error")

#' Headline metrics per model, one panel per metric, optionally one row per level
#' of a factor. The dotted line in each panel is the target value.
plot_metrics <- function(long, facet_row = NULL, title) {
  long <- long[long$metric %in% HEADLINE_METRICS, ]
  long$metric <- factor(long$metric, levels = HEADLINE_METRICS)
  # The value each metric should land on, in the same order as HEADLINE_METRICS.
  refs <- data.frame(metric = factor(HEADLINE_METRICS, levels = HEADLINE_METRICS),
                     ref = c(1, 0.9, 0))

  p <- ggplot(long, aes(value, model, colour = model, shape = model)) +
    geom_vline(data = refs, aes(xintercept = ref), colour = "grey60",
               linetype = "dotted", linewidth = 0.4, inherit.aes = FALSE) +
    geom_errorbar(aes(xmin = lo, xmax = hi), orientation = "y", width = 0,
                  linewidth = 0.6, na.rm = TRUE) +
    geom_point(size = 2.4) +
    scale_colour_manual(values = MODEL_COLOURS) +
    scale_shape_manual(values = MODEL_SHAPES) +
    scale_x_continuous(n.breaks = 4) +
    labs(title = title, x = NULL, y = NULL) +
    panel_theme +
    theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())
  if (is.null(facet_row))
    p + facet_wrap(~ metric, scales = "free_x", nrow = 1)
  else
    p + facet_grid(stats::reformulate("metric", facet_row), scales = "free_x")
}

#' Bias, attenuation and accuracy against dataset size
plot_by_n <- function(d, title) {
  by_n <- summarise_models(d, c("model", "N"))
  long <- droplevels(metric_long(by_n, c("model", "N")))
  long <- long[long$metric %in% METRICS[1:3], ]
  refs <- data.frame(metric = factor(METRICS[1:3], levels = levels(long$metric)),
                     ref = c(0, 1, 0.9))

  ggplot(long, aes(N, value, colour = model, shape = model)) +
    geom_hline(data = refs, aes(yintercept = ref), colour = "grey60",
               linetype = "dotted", linewidth = 0.4, inherit.aes = FALSE) +
    geom_errorbar(aes(ymin = lo, ymax = hi), width = 0, linewidth = 0.4,
                  alpha = 0.5, na.rm = TRUE) +
    geom_line(aes(group = model), linewidth = 0.7) +
    geom_point(size = 2) +
    scale_colour_manual(values = MODEL_COLOURS) +
    scale_shape_manual(values = MODEL_SHAPES) +
    facet_wrap(~ metric, scales = "free_y", nrow = 1) +
    scale_x_continuous(trans = "log2", breaks = sort(unique(long$N))) +
    labs(title = title, x = "Sample size", y = NULL) +
    panel_theme
}

#' Table and figure for one sweep: writes recovery_by_<name>.csv and .png.
#' Does nothing if the results do not vary in `column`, and stops if the column
#' is missing (the 03_ script did not keep it).
sweep_figure <- function(d, column, name, title, out_dir, fig_dir) {
  if (!column %in% names(d))
    stop("results have no '", column, "' column; add it to the keep vector in the 03_ script")
  if (length(unique(d[[column]])) < 2) return(invisible(NULL))
  stats <- summarise_models(d, c("model", column))
  write.csv(stats, file.path(out_dir, paste0("recovery_by_", name, ".csv")),
            row.names = FALSE)
  ggsave(file.path(fig_dir, paste0("recovery_by_", name, ".png")),
         plot_metrics(metric_long(stats, c("model", column)), facet_row = column,
                      title = title),
         # One facet row per level, plus a fixed allowance for the strip and axis.
         width = 8, height = 2.2 + 1.6 * length(unique(stats[[column]])),
         dpi = 300, bg = "white")
  invisible(stats)
}

#' Zero-slope datasets: how often the 90% interval excludes zero
false_positive_rate <- function(d) {
  zero <- d[d$slope_condition == "zero", ]
  if (nrow(zero) == 0) return(NULL)
  round(tapply(zero$slope_cov90, zero$model, function(x) 1 - mean(x)), 3)
}

#' One simulated dataset per panel: each find at its measured value, its dating
#' range as a segment, its point date as a square, its true date as a red dot,
#' and the true trend dashed.
#'
#' @param d        one row per find: panel, Start_date, End_date, True_date,
#'                 Value, Point_date
#' @param ribbon   optional phase bands: panel, xmin, xmax
#' @param trend    optional true trends: panel, intercept, slope
#' @param free_x   TRUE when panels cover different calendar ranges
#' @param point_label, span_label  legend names for the square and the segment
dataset_anatomy <- function(d, ribbon = NULL, trend = NULL, title = NULL,
                            free_x = FALSE,
                            point_label = "window midpoint",
                            span_label  = "dating window") {
  need <- c("panel", "Start_date", "End_date", "True_date", "Value", "Point_date")
  if (!all(need %in% names(d)))
    stop("dataset_anatomy() needs columns: ", paste(need, collapse = ", "))
  d <- d[order(d$panel, d$Value), ]

  p <- ggplot(d, aes(y = Value))

  if (!is.null(ribbon)) {
    # Alternate shading of phase bands within each panel
    ribbon$panel <- factor(as.character(ribbon$panel), levels = levels(d$panel))
    ribbon <- do.call(rbind, by(ribbon, ribbon$panel, function(g) {
      g <- g[order(g$xmin), ]
      g$shade <- rep_len(c("a", "b"), nrow(g))
      g
    }))
    p <- p +
      geom_rect(data = ribbon, inherit.aes = FALSE, alpha = 0.5,
                aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf,
                    fill = shade)) +
      scale_fill_manual(values = c(a = "grey92", b = "grey82"), guide = "none")
  }

  if (!is.null(trend))
    p <- p +
      geom_abline(data = trend, aes(intercept = intercept, slope = slope),
                  linetype = "dashed", colour = "grey20", linewidth = 0.4)

  # One legend for the segment and both kinds of point
  keys <- c(point_label, "simulated date", span_label)

  p +
    geom_segment(aes(x = Start_date, xend = End_date, yend = Value,
                     colour = span_label),
                 linewidth = 0.3, alpha = 0.7) +
    geom_point(aes(x = Point_date, shape = point_label, colour = point_label),
               size = 1) +
    geom_point(aes(x = True_date, shape = "simulated date",
                   colour = "simulated date"),
               size = 0.9) +
    scale_shape_manual(breaks = keys, limits = keys,
                       values = stats::setNames(c(15, 16, NA), keys)) +
    scale_colour_manual(breaks = keys, limits = keys,
                        values = stats::setNames(
                          c("grey45", ACCENT, "grey65"), keys)) +
    facet_wrap(~ panel, scales = if (free_x) "free_x" else "fixed") +
    labs(title = title, shape = NULL, colour = NULL,
         x = "calendar year", y = "measured value") +
    # Own theme, so it looks the same whichever script calls it
    theme_classic(base_size = 11) +
    theme(plot.title = element_text(size = 11, face = "bold"),
          strip.background = element_blank(),
          strip.text = element_text(face = "bold"),
          legend.position = "top", legend.title = element_blank())
}
