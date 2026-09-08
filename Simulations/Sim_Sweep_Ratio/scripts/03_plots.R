# Table and figures for the sweep.
#
# Metrics against dating coarseness, one row of panels per signal level.
# Plotted as curves rather than as the cases' dot-and-interval panels, because
# here the x axis is a continuous quantity a reader looks their own study up on.

library(here)
library(ggplot2)

source(here("Simulations", "shared", "scripts", "recovery_summary.R"))

results_csv <- Sys.getenv("SWEEP_OUT",
                          here("Simulations", "Sim_Sweep_Ratio", "output",
                               "sweep_results.csv"))
fig_dir <- here("Simulations", "Sim_Sweep_Ratio", "figures")
out_dir <- here("Simulations", "Sim_Sweep_Ratio", "output")

LABELS  <- c(oracle = "True dates (control)", midpoint = "Midpoint",
             marginal = "Date-marginalised")
COLOURS <- c(`True dates (control)` = "grey75", Midpoint = "grey45",
             `Date-marginalised` = "#780000")
SHAPES  <- c(`True dates (control)` = 1, Midpoint = 17,
             `Date-marginalised` = 16)

r <- read.csv(results_csv)
n_failed <- sum(!is.na(r$error))
if (n_failed > 0) cat(sprintf("%d of %d fits errored and are excluded\n",
                              n_failed, nrow(r)))
r <- r[is.na(r$error), ]
r$model <- factor(LABELS[r$model], levels = unname(LABELS))

stats <- do.call(rbind, lapply(
  split(seq_len(nrow(r)), r[c("model", "coarseness", "rise_level")], drop = TRUE),
  function(i) {
    s <- r[i, ]
    a <- attenuation(s)
    p <- jeffreys_prop(s$slope_cov90)
    sg <- mc_mean(s$sigma_err)
    data.frame(model = s$model[1], coarseness = s$coarseness[1],
               rise_level = s$rise_level[1], period_span = s$period_span[1],
               n_fits = nrow(s),
               atten = a["mean"], atten_lo = a["lo"], atten_hi = a["hi"],
               accuracy = p["mean"], accuracy_lo = p["lo"],
               accuracy_hi = p["hi"],
               sigma_err = sg["mean"], sigma_err_se = sg["se"],
               row.names = NULL)
  }))
stats <- stats[order(stats$model, stats$rise_level, stats$coarseness), ]
stats$signal <- factor(sprintf("total rise +/- %g", stats$rise_level),
                       levels = sprintf("total rise +/- %g",
                                        sort(unique(stats$rise_level))))
write.csv(stats, file.path(out_dir, "sweep_table.csv"), row.names = FALSE)
print(stats, row.names = FALSE, digits = 3)

keep <- c("model", "coarseness", "period_span", "signal")
long <- rbind(
  data.frame(stats[keep], metric = "Slope attenuation (target 1)",
             value = stats$atten, lo = stats$atten_lo, hi = stats$atten_hi),
  data.frame(stats[keep], metric = "Slope accuracy (target 0.9)",
             value = stats$accuracy, lo = stats$accuracy_lo,
             hi = stats$accuracy_hi),
  data.frame(stats[keep], metric = "Sigma error (target 0)",
             value = stats$sigma_err,
             lo = stats$sigma_err - 2 * stats$sigma_err_se,
             hi = stats$sigma_err + 2 * stats$sigma_err_se))
long$metric <- factor(long$metric, levels = unique(long$metric))
refs <- data.frame(metric = factor(levels(long$metric),
                                   levels = levels(long$metric)),
                   ref = c(1, 0.9, 0))

p <- ggplot(long, aes(coarseness, value, colour = model, shape = model)) +
  geom_hline(data = refs, aes(yintercept = ref), colour = "grey60",
             linetype = "dotted", linewidth = 0.4, inherit.aes = FALSE) +
  geom_errorbar(aes(ymin = lo, ymax = hi), width = 0, linewidth = 0.4,
                alpha = 0.6, na.rm = TRUE) +
  geom_line(aes(group = model), linewidth = 0.7) +
  geom_point(size = 2.2) +
  scale_colour_manual(values = COLOURS) +
  scale_shape_manual(values = SHAPES) +
  scale_x_continuous(trans = "log2", breaks = sort(unique(stats$coarseness)),
                     labels = function(x) sprintf("%.2f", x)) +
  facet_grid(signal ~ metric, scales = "free_y") +
  labs(title = "What coarse dating costs, at three signal strengths",
       subtitle = paste("Window widths fixed at 50-200 yr; the studied period",
                        "varies from 200 to 3200 yr"),
       x = "dating coarseness  (mean window width / studied period)", y = NULL,
       caption = paste("Accuracy: Jeffreys 95% interval. Attenuation: 95%",
                       "interval on the regression of estimate on truth.\n",
                       "The control is fitted on the true dates and must sit",
                       "at the target in every cell.")) +
  panel_theme
ggsave(file.path(fig_dir, "sweep_ratio.png"), p,
       width = 9, height = 8, dpi = 300, bg = "white")

cat("\nwrote figures to", fig_dir, "\n")
