# Case 1 tables and figures from output/recovery_results.csv.
# Metrics and plots live in shared/scripts/recovery_summary.R.

library(here)

source(here("Simulations", "shared", "scripts", "recovery_summary.R"))

results_csv <- Sys.getenv("RECOVERY_OUT",
                          here("Simulations", "Sim_Case1_Radiocarbon", "output",
                               "recovery_results.csv"))
fig_dir <- here("Simulations", "Sim_Case1_Radiocarbon", "figures")
chk_dir <- file.path(fig_dir, "checks")
out_dir <- here("Simulations", "Sim_Case1_Radiocarbon", "output")
dir.create(chk_dir, showWarnings = FALSE)

r <- read_recovery(results_csv)

# Readable window names, on two lines so the facet strips fit
WINDOW_LABELS <- c(plateau = "Hallstatt plateau\n(800–400 BCE)",
                   steep   = "Steep section\n(1600–1200 BCE)")
r$window <- factor(unname(WINDOW_LABELS[r$window]), levels = WINDOW_LABELS)

# Headline results use even deposition only; growth has its own figures below
reference <- r[r$growth_ratio == 1, ]

by_model <- summarise_models(reference, c("model", "window"))
write.csv(by_model, file.path(out_dir, "recovery_table.csv"), row.names = FALSE)
print(by_model, row.names = FALSE)

p1 <- plot_metrics(metric_long(by_model, c("model", "window")), facet_row = "window",
                   title = "Case 1 (radiocarbon): parameter recovery by model",
                   caption = paste0(
                     "Each point summarises all simulated datasets; dotted lines mark the target. Bars are 95% intervals:\n",
                     "OLS confidence interval for the calibration slope b in estimated\u1d62 = a + b \u00b7 true\u1d62 + \u03b5\u1d62 ",
                     "(slope in dataset i; estimated = posterior median);\n",
                     "Jeffreys interval for coverage; \u00b12 Monte Carlo SE for bias in \u03c3 (posterior median minus true \u03c3)."))
ggsave(file.path(fig_dir, "recovery_summary.png"), p1, width = 8,
       height = 2.2 + 1.6 * length(unique(by_model$window)), dpi = 300, bg = "white")

core <- reference[reference$sweep == "core" &
                  reference$slope_condition == "random", ]
if (nrow(core) > 0)
  ggsave(file.path(chk_dir, "recovery_by_n.png"),
         plot_by_n(core, "Case 1: more data does not fix a biased estimator"),
         width = 8, height = 4.2, dpi = 300, bg = "white")

# Lab error sweep
sweep_figure(reference[reference$sweep == "factor", ], "lab_error", "lab_error",
             "Case 1: effect of lab error", out_dir, fig_dir)

# Deposition: even vs 4 times denser at the end, lab error 30, one figure per
# window so the plateau and the steep section are not averaged together
deposition <- r[r$sweep %in% c("factor", "deposition") &
                r$lab_error == 30 & r$slope_condition == "random", ]
for (code in names(WINDOW_LABELS))
  sweep_figure(deposition[deposition$window == WINDOW_LABELS[code], ], "growth_ratio",
               paste0("deposition_", code),
               sprintf("Case 1: uniform vs growing deposition, %s",
                       sub("\n", " ", WINDOW_LABELS[code])),
               out_dir, chk_dir)

# Convergence: R-hat and divergent transitions, by sweep and model
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
