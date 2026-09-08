# Turning one-row-per-fit results into the paper's result axes.
#
# Shared by every case's 04_ script: the metrics are the same whatever the
# dating is, and only the factor a case is faceted by differs.
#
# Vocabulary follows the earlier studies rather than inventing new terms:
#
#   bias         how far the posterior median sits from the truth, on average
#   attenuation  how much of the true slope survives estimation (target 1)
#   accuracy     how often the 90% interval contained the truth (target 0.9).
#                 This is what the statistical literature calls coverage.
#   precision    how wide that interval was. Narrower is better, but only once
#                 accuracy is on target - a narrow interval that misses the truth
#                 is worse than a wide one that does not.
#
# Bias and accuracy are reported side by side and never combined into a single
# score. A model can be unbiased and still badly calibrated, which is the point
# of Claim A.
#
# Bias alone cannot see attenuation, and in these designs it never will. True
# slopes are drawn symmetrically about zero, so the mean truth is ~0; an error
# that scales with the truth then averages to (c - 1) * 0 = 0 whatever c is,
# because the positive-slope and negative-slope datasets cancel. Case 1's median
# model reports a mean bias of +6e-06 while flattening every slope by 18%.
# Attenuation is therefore reported alongside bias, not instead of it: bias
# catches an error of fixed size, attenuation an error proportional to the
# effect, and the two are blind to each other.
#
# Greyscale with a single dark-red accent for the model being argued for, after
# the palette in the archived Sim_Linear figures. Prints legibly in black and
# white, and shape carries the model as well as colour, so neither is
# load-bearing on its own.

library(ggplot2)

# All five metrics. metric_long() carries every one into recovery_table.csv; the
# headline figures show only the three that carry the argument (HEADLINE_METRICS,
# defined next to plot_metrics below).
METRICS <- c("Slope bias", "Slope attenuation", "Slope accuracy (90%)",
             "Slope precision (90% width)", "Sigma error")

# "marginal" plays the role the earlier studies labelled EIV. The point-date
# baselines recede; the model under test carries the accent.
MODEL_LEVELS  <- c("midpoint", "median", "latent", "marginal")
MODEL_LABELS  <- c(midpoint = "Midpoint", median = "Median",
                   latent = "Latent date", marginal = "Date-marginalised")
MODEL_COLOURS <- c(Midpoint = "grey55", `Midpoint / Median` = "grey55",
                   Median = "grey25", `Latent date` = "grey40",
                   `Date-marginalised` = "#780000")
MODEL_SHAPES  <- c(Midpoint = 17, `Midpoint / Median` = 17, Median = 15,
                   `Latent date` = 18, `Date-marginalised` = 16)

# The one non-grey in every figure: the model being argued for, or the true
# calendar date in the anatomy figures. Kept here so it is defined once.
ACCENT <- "#780000"

panel_theme <- theme_classic(base_size = 11) +
  theme(plot.title = element_text(size = 11, face = "bold"),
        plot.subtitle = element_text(size = 9),
        strip.background = element_blank(),
        strip.text = element_text(face = "bold"),
        legend.position = "top",
        legend.title = element_blank(),
        plot.caption = element_text(size = 8, colour = "grey35", hjust = 0))

#' Drop failed fits and label the models for plotting.
#'
#' merge_median is for the flat-window cases, where median was not fitted
#' because it is the same estimator as midpoint (00_check_median_identity.R).
#' The surviving row is relabelled so a reader is not left wondering which of
#' the two they are looking at. Passed explicitly rather than inferred, so
#' shared code never guesses what a case meant.
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

#' Mean and its Monte Carlo standard error: how precisely this summary is known
#' given only n simulated datasets, so a gap between models can be read against
#' simulation noise rather than eyeballed. Shrinks as 1/sqrt(n_datasets).
mc_mean <- function(x) {
  x <- x[is.finite(x)]
  c(mean = mean(x), se = sd(x) / sqrt(length(x)))
}

#' Accuracy is a proportion, so it gets a Jeffreys interval rather than a
#' normal-approximation standard error.
#'
#' The normal approximation, p +/- 2 sqrt(p(1-p)/n), misbehaves at the
#' boundaries: at accuracy 0.98 with n = 100 it draws a bar running past 1.0,
#' which is meaningless. Accuracy near 1 is the expected result for the
#' date-marginalised model, so the boundary is not a corner case here.
#'
#' Jeffreys puts a Beta(0.5, 0.5) prior on the proportion, giving a posterior
#' Beta(k + 0.5, n - k + 0.5); the interval is its equal-tailed quantiles. It
#' stays inside [0, 1] and goes asymmetric near the edges, which is the honest
#' shape. Asymmetric, so it is returned as bounds rather than as one SE.
jeffreys_prop <- function(x, level = 0.95) {
  x <- x[!is.na(x)]
  n <- length(x)
  k <- sum(x)
  a <- (1 - level) / 2
  c(mean = k / n,
    lo = qbeta(a,     k + 0.5, n - k + 0.5),
    hi = qbeta(1 - a, k + 0.5, n - k + 0.5))
}

#' Attenuation: the share of the true slope that survives estimation.
#'
#' Every dataset in a recovery study has its own true slope, so the estimates
#' can be regressed on the truths that produced them. The fitted line's slope is
#' the attenuation factor c: an estimator that recovers the truth one-for-one
#' gives c = 1, one that systematically flattens gives c < 1, one that
#' exaggerates gives c > 1. Distance from 1 is what is being read, in either
#' direction.
#'
#' The truth is reconstructed as slope_med - slope_err rather than joined back
#' from the design, so this works on a results file alone and cannot go out of
#' step with a design that was rebuilt.
#'
#' Returns NA where the truths do not vary - the zero-slope cells have no
#' regression to fit, and asking for one there is a category error, not a
#' failure.
attenuation <- function(d, level = 0.95) {
  truth <- d$slope_med - d$slope_err
  ok    <- is.finite(truth) & is.finite(d$slope_med)
  if (sum(ok) < 3 || sd(truth[ok]) < .Machine$double.eps^0.5)
    return(c(mean = NA_real_, lo = NA_real_, hi = NA_real_))
  fit <- stats::lm(d$slope_med[ok] ~ truth[ok])
  ci  <- stats::confint(fit, 2, level = level)
  c(mean = unname(coef(fit)[2]), lo = ci[1], hi = ci[2])
}

summarise_models <- function(d, by) {
  stats <- do.call(rbind, lapply(split(seq_len(nrow(d)), d[by], drop = TRUE),
    function(i) {
      s     <- d[i, ]
      bias  <- mc_mean(s$slope_err)
      atten <- attenuation(s)
      acc   <- jeffreys_prop(s$slope_cov90)
      prec  <- mc_mean(s$slope_width90)
      sig   <- mc_mean(s$sigma_err)
      data.frame(s[1, by, drop = FALSE], n_fits = nrow(s),
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

#' Wide summary to one row per (group, metric). rbind() drops factor level
#' order, so the levels are reapplied rather than trusted.
metric_long <- function(stats, keep) {
  long <- do.call(rbind, lapply(seq_len(nrow(stats)), function(i) {
    s <- stats[i, ]
    # Bounds rather than a single se: accuracy and attenuation carry their own
    # asymmetric or regression intervals, the rest are +/- 2 Monte Carlo SE.
    # Both kinds are named in the caption.
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

# The three axes that carry the argument. Bias is ~0 by construction here (true
# slopes are symmetric about zero) and precision is only readable once accuracy
# is on target, so both live in recovery_table.csv rather than the headline
# figure. Cutting from five panels to three is also what makes the strip labels
# and axis numbers fit without colliding.
HEADLINE_METRICS <- c("Slope attenuation", "Slope accuracy (90%)", "Sigma error")

#' The headline metrics per model, optionally faceted by a case's own factor.
#'
#' One panel per metric, each with its own x scale - attenuation sits near 1,
#' accuracy near 0.9, sigma error near 0, and forcing them onto a shared axis
#' would flatten the differences that matter. The dotted line in each panel is
#' that metric's target.
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
    labs(title = title, x = NULL, y = NULL,
         caption = paste("Accuracy: Jeffreys 95% interval.",
                         "Attenuation: 95% interval on the regression of",
                         "estimate on truth.\nSigma error: +/- 2 Monte Carlo",
                         "standard errors. Bias and precision are in",
                         "recovery_table.csv.")) +
    panel_theme +
    theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())
  if (is.null(facet_row))
    p + facet_wrap(~ metric, scales = "free_x", nrow = 1)
  else
    p + facet_grid(stats::reformulate("metric", facet_row), scales = "free_x")
}

#' Bias, attenuation and accuracy against sample size. All three should stay
#' flat as N grows and only the spread should shrink; that is what separates a
#' systematic error from noise. Accuracy falling as N rises is the signature of
#' a systematic error the intervals do not know about - more data narrows them
#' around the wrong value.
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
    labs(title = title, x = "Sample size", y = NULL,
         caption = paste("Accuracy and attenuation: 95% intervals.",
                         "Bias: +/- 2 Monte Carlo standard errors.")) +
    panel_theme
}

#' Was a sweep actually run, and is the column there to tell?
#'
#' Guards the `if` around each optional figure. The distinction matters: a
#' column holding one value means that sweep was not part of this design, which
#' is a fine reason to skip the figure. A column that is absent entirely means
#' the study's 03_ script did not carry it into the results, which is a bug -
#' and testing `length(unique(d$missing)) > 1` silently reports FALSE for both,
#' so the figure just never appears. That happened to Case 2's overhang sweep:
#' 400 fits ran and were scored, and nothing ever plotted them.
has_sweep <- function(d, column) {
  if (!column %in% names(d))
    stop("results have no '", column, "' column, so the figure grouped by it ",
         "cannot be drawn. Add it to the keep vector in the 03_ script and ",
         "re-run, or join it from the design.")
  length(unique(d[[column]])) > 1
}

#' One sweep, as a table and a figure: the models compared across the levels of
#' whatever factor this case varies.
#'
#' Every case asks the same question of its own factor - lab error, dating
#' resolution, deposition skew, window overhang, broad-period width - and the
#' answer is always the same table and the same faceted figure. Only the column,
#' the file name and the title differ, so those are the arguments.
#'
#' @param d       the rows this sweep is read from, already filtered by the case
#' @param column  the design column whose levels the facets are
#' @param name    file stem: writes recovery_by_<name>.csv and .png
sweep_figure <- function(d, column, name, title, out_dir, fig_dir) {
  if (!has_sweep(d, column)) return(invisible(NULL))
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

#' Zero-slope cells: how often a 90% interval excludes the true zero.
false_positive_rate <- function(d) {
  zero <- d[d$slope_condition == "zero", ]
  if (nrow(zero) == 0) return(NULL)
  round(tapply(zero$slope_cov90, zero$model, function(x) 1 - mean(x)), 3)
}

# --- what a dataset looks like before any model sees it -----------------------

#' The structure of a single dataset, in the idiom of the archived Sim_Linear
#' scenario figures: one row per find at its measured value, the dating window a
#' faint horizontal segment, the true calendar date a dark-red dot, the
#' point-date estimate a grey square, the true trend a dashed line. Where the
#' dating is by phases, the phase boundaries are alternating vertical bands
#' behind the cloud, so it is visible which phase each find sits in.
#'
#' One representative dataset per small multiple. What the panels vary - lab
#' error, overlap, merge width, plateau vs control - is the case's own factor,
#' passed in the `panel` column already labelled.
#'
#' The helper never inspects the case. Each case builds the frame it needs:
#' `Point_date` is the calibrated median for radiocarbon and the window midpoint
#' for the flat-window cases, computed by the caller, not here.
#'
#' @param d       one row per find: panel, Start_date, End_date, True_date,
#'                 Value, Point_date.
#' @param ribbon  optional phase bands: panel, xmin, xmax. Alternate shading is
#'                 assigned from left to right within each panel.
#' @param trend   optional true lines: panel, intercept, slope.
#' @param free_x  TRUE when panels sit on different calendar ranges (Case 1's two
#'                 windows do), FALSE when they share one axis.
#' @param point_label,span_label  what the grey square and the grey segment are
#'                 in this case, named in the legend. They are not the same
#'                 object from case to case: for a flat window the square is the
#'                 midpoint and the segment is the recorded window; for
#'                 radiocarbon the square is the calibrated median and the
#'                 segment is the 95% HPD envelope. Naming them in the figure
#'                 rather than only in the caption is the point - the gap
#'                 between the square and the red dot is the dating error, and
#'                 a reader should not have to infer which is which.
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
    # Match the ribbon's panels to the data's, so shading and facet order do not
    # depend on how the caller ordered its rows.
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

  # One legend for all three marks, so colour and shape carry the same breaks.
  # The segment has no shape (NA) and the points have no line; ggplot draws each
  # key from the layer that uses it.
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
    # Its own theme rather than the shared panel_theme: this helper is called
    # from the 00_check scripts, which each define a panel_theme of their own,
    # and the figure should look the same whichever one is in scope.
    theme_classic(base_size = 11) +
    theme(plot.title = element_text(size = 11, face = "bold"),
          strip.background = element_blank(),
          strip.text = element_text(face = "bold"),
          legend.position = "top", legend.title = element_blank())
}
