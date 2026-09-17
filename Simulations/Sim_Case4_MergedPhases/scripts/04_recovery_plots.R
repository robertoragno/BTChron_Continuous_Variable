# Case 4 tables and figures from output/recovery_results.csv.
# Metrics and plots live in shared/scripts/recovery_summary.R; this file only
# says which sweeps Case 4 has.

library(here)
library(ggplot2)

source(here("Simulations", "shared", "scripts", "recovery_summary.R"))
source(here("Simulations", "Sim_Case4_MergedPhases", "scripts", "simulate.R"))

results_csv <- Sys.getenv("RECOVERY_OUT",
                          here("Simulations", "Sim_Case4_MergedPhases",
                               "output", "recovery_results.csv"))
fig_dir <- here("Simulations", "Sim_Case4_MergedPhases", "figures")
chk_dir <- file.path(fig_dir, "checks")
out_dir <- here("Simulations", "Sim_Case4_MergedPhases", "output")
dir.create(chk_dir, showWarnings = FALSE)

# median is not fitted for flat windows
r <- read_recovery(results_csv, merge_median = TRUE)

core <- r[r$sweep == "core" & r$slope_condition == "random", ]
by_model <- summarise_models(core, "model")
write.csv(by_model, file.path(out_dir, "recovery_table.csv"), row.names = FALSE)
print(by_model, row.names = FALSE)

ggsave(file.path(fig_dir, "recovery_summary.png"),
       plot_metrics(metric_long(by_model, "model"),
                    title = "Case 4: merged phases"),
       width = 8, height = 3.4, dpi = 300, bg = "white")

ggsave(file.path(chk_dir, "recovery_by_n.png"),
       plot_by_n(core, "Case 4: effect of dataset size"),
       width = 8, height = 4.2, dpi = 300, bg = "white")

sweep_figure(r[r$sweep == "merge", ], "merge_max", "merge",
             "Case 4: largest broad period, in phases", out_dir, fig_dir)

# One simulated dataset per merge level, fine phases as alternating bands
merge_levels <- c(2, 3, 4)
panels <- lapply(merge_levels, function(mm) {
  sim <- simulate_merged(70, 8, 0.02, 2, K = 8, alpha_conc = 1,
                         merge_max = mm, seed = 4)
  b <- attr(sim, "bounds")
  sim$panel <- paste("up to", mm, "phases")
  sim$Point_date <- (sim$Start_date + sim$End_date) / 2
  list(sim = sim,
       ribbon = data.frame(panel = sim$panel[1], xmin = b[-length(b)], xmax = b[-1]))
})
labels <- paste("up to", merge_levels, "phases")
anat   <- do.call(rbind, lapply(panels, function(p) p$sim))
ribbon <- do.call(rbind, lapply(panels, function(p) p$ribbon))
anat$panel <- factor(anat$panel, levels = labels)
trend <- data.frame(panel = labels, intercept = 8, slope = 0.02)
ggsave(file.path(fig_dir, "dataset_anatomy.png"),
       dataset_anatomy(anat, ribbon = ribbon, trend = trend),
       width = 9, height = 3.6, dpi = 300, bg = "white")

fp <- false_positive_rate(r)
if (!is.null(fp)) {
  cat("\nfalse-positive rate (zero-slope cells, 90% interval excludes 0):\n")
  print(fp)
}

cat("wrote figures to", fig_dir, "\n")
