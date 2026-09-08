# Case 1 tables and figures from output/recovery_results.csv.
# The metrics and their plots live in shared/scripts/recovery_summary.R; this
# file only says what Case 1 is faceted by, which is the calendar window.

library(here)

source(here("Simulations", "shared", "scripts", "recovery_summary.R"))

results_csv <- Sys.getenv("RECOVERY_OUT",
                          here("Simulations", "Sim_Case1_Radiocarbon", "output",
                               "recovery_results.csv"))
fig_dir <- here("Simulations", "Sim_Case1_Radiocarbon", "figures")
out_dir <- here("Simulations", "Sim_Case1_Radiocarbon", "output")

r <- read_recovery(results_csv)

# The design names the second window "steep"; the paper and this case's README
# call it the control. Use the reader-facing name in every table and figure.
r$window <- factor(unname(c(plateau = "plateau", steep = "control")[r$window]),
                   levels = c("plateau", "control"))

# The headline is the reference condition, so the deposition sweep is held out
# of it: growth appears in its own figures below and nowhere else. Without this
# the headline silently averages 100 growth datasets per window into the totals
# and stops being comparable with earlier runs.
reference <- r[r$growth_ratio == 1, ]

by_model <- summarise_models(reference, c("model", "window"))
write.csv(by_model, file.path(out_dir, "recovery_table.csv"), row.names = FALSE)
print(by_model, row.names = FALSE)

p1 <- plot_metrics(metric_long(by_model, c("model", "window")), facet_row = "window",
                   title = "Case 1: what each model gets right and wrong")
ggsave(file.path(fig_dir, "recovery_summary.png"), p1, width = 8,
       height = 2.2 + 1.6 * length(unique(by_model$window)), dpi = 300, bg = "white")

core <- reference[reference$sweep == "core" &
                  reference$slope_condition == "random", ]
if (nrow(core) > 0)
  ggsave(file.path(fig_dir, "recovery_by_n.png"),
         plot_by_n(core, "Case 1: more data does not fix a biased estimator"),
         width = 8, height = 4.2, dpi = 300, bg = "white")

# Lab error is Case 1's other factor: wider calibrated posteriors should deepen
# attenuation in the point-date models and leave the marginal model alone.
sweep_figure(reference[reference$sweep == "factor", ], "lab_error", "lab_error",
             "Case 1: effect of lab error", out_dir, fig_dir)

# Deposition: uniform against fourfold growth in find density across the window.
# The reference cells are the `factor` ones at the reference lab error, so the
# only thing that differs between the two rows of each panel is the deposition
# shape. One figure per window rather than one pooled: the whole point is
# whether growth bites on the plateau and not off it, and a pooled figure would
# average that away.
deposition <- r[r$sweep %in% c("factor", "deposition") &
                r$lab_error == 30 & r$slope_condition == "random", ]
for (w in levels(droplevels(deposition$window)))
  sweep_figure(deposition[deposition$window == w, ], "growth_ratio",
               paste0("deposition_", w),
               sprintf("Case 1: uniform vs growing deposition (%s)", w),
               out_dir, fig_dir)

# Convergence before conclusions. read_recovery() has already dropped and counted
# the fits that errored outright; these are the ones that returned an answer but
# may not have earned it. Broken down by sweep so a problem in one block - the
# deposition cells, say - shows up isolated rather than averaged into the total.
cat("\nconvergence:", sum(r$max_rhat > 1.01, na.rm = TRUE), "fits with rhat > 1.01,",
    sum(r$max_rhat > 1.05, na.rm = TRUE), "above 1.05, worst",
    sprintf("%.4f", max(r$max_rhat, na.rm = TRUE)), "\n")
cat(sprintf("%.1f%% of fits had a divergent transition, worst %d\n",
            100 * mean(r$n_divergent > 0), max(r$n_divergent)))
cat("\nworst rhat by sweep and model:\n")
print(round(tapply(r$max_rhat, list(r$sweep, r$model), max, na.rm = TRUE), 3))
cat("\n% of fits with any divergence, by sweep and model:\n")
print(round(tapply(r$n_divergent > 0, list(r$sweep, r$model), mean) * 100, 1))

fp <- false_positive_rate(reference)
if (!is.null(fp)) {
  cat("\nfalse-positive rate (zero-slope cells, 90% interval excludes 0):\n")
  print(fp)
}

cat(sprintf("\nHPD envelope retains %.1f%% of the calibrated mass on average\n",
            100 * mean(r$mass_kept, na.rm = TRUE)))
cat("wrote figures to", fig_dir, "\n")
