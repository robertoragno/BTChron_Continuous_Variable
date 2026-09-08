# Case 2 tables and figures from output/recovery_results.csv.
# Metrics and plots live in shared/scripts/recovery_summary.R; this file only
# says what Case 2 is faceted by, which is dating resolution and deposition skew.

library(here)

source(here("Simulations", "shared", "scripts", "recovery_summary.R"))

results_csv <- Sys.getenv("RECOVERY_OUT",
                          here("Simulations", "Sim_Case2_Typochronology",
                               "output", "recovery_results.csv"))
fig_dir <- here("Simulations", "Sim_Case2_Typochronology", "figures")
out_dir <- here("Simulations", "Sim_Case2_Typochronology", "output")

# median is not fitted for flat windows; see 00_check_median_identity.R
r <- read_recovery(results_csv, merge_median = !("median" %in% read.csv(results_csv)$model))

by_model <- summarise_models(r, "model")
write.csv(by_model, file.path(out_dir, "recovery_table.csv"), row.names = FALSE)
print(by_model, row.names = FALSE)

# The headline for Case 2: the point-date models should be close to unbiased here
# (Berkson), and still miss on accuracy and sigma. Bias and accuracy side by
# side is the whole argument, so the two are never collapsed into one number.
p1 <- plot_metrics(metric_long(by_model, "model"),
                   title = "Case 2: unbiased is not the same as calibrated")
ggsave(file.path(fig_dir, "recovery_summary.png"), p1,
       width = 8, height = 3.4, dpi = 300, bg = "white")

core <- r[r$sweep == "core" & r$slope_condition == "random", ]
if (nrow(core) > 0)
  ggsave(file.path(fig_dir, "recovery_by_n.png"),
         plot_by_n(core, "Case 2: precision improves with N, accuracy does not"),
         width = 8, height = 4.2, dpi = 300, bg = "white")

# Dating resolution: coarser windows should widen intervals and inflate sigma
# without moving the slope much.
sweep_figure(r[r$sweep %in% c("core", "resolution") & r$slope_condition == "random", ],
             "width_mix", "resolution", "Case 2: effect of dating resolution",
             out_dir, fig_dir)

# Skew is where Berkson breaks: with deposition leaning late, the point-date
# models should pick up real slope bias.
sweep_figure(r[r$sweep %in% c("core", "skew") & r$slope_condition == "random", ],
             "skew_shape", "skew",
             "Case 2b: what breaks when deposition is not uniform",
             out_dir, fig_dir)

# Overhang: how much the boundary flattens the slope, as a function of how wide
# windows are relative to the study period. Reported as a result rather than
# engineered away - it is a real property of typological dating.
sweep_figure(r[r$sweep %in% c("core", "overhang") & r$slope_condition == "random" &
               r$N == 200, ],
             "max_frac", "overhang",
             "Case 2: window width relative to the study period",
             out_dir, fig_dir)

fp <- false_positive_rate(r)
if (!is.null(fp)) {
  cat("\nfalse-positive rate (zero-slope cells, 90% interval excludes 0):\n")
  print(fp)
}

cat("wrote figures to", fig_dir, "\n")
