# Case 3 tables and figures from output/recovery_results.csv.
# Metrics and plots live in shared/scripts/recovery_summary.R; this file only
# says which sweeps Case 3 has.

library(here)
library(ggplot2)

source(here("Simulations", "shared", "scripts", "recovery_summary.R"))
source(here("Simulations", "Sim_Case3_OverlappingPhases", "scripts", "simulate.R"))

results_csv <- Sys.getenv("RECOVERY_OUT",
                          here("Simulations", "Sim_Case3_OverlappingPhases",
                               "output", "recovery_results.csv"))
fig_dir <- here("Simulations", "Sim_Case3_OverlappingPhases", "figures")
chk_dir <- file.path(fig_dir, "checks")
out_dir <- here("Simulations", "Sim_Case3_OverlappingPhases", "output")
dir.create(chk_dir, showWarnings = FALSE)

# median is not fitted for flat windows
r <- read_recovery(results_csv, merge_median = TRUE)

core <- r[r$sweep == "core" & r$slope_condition == "random", ]
by_model <- summarise_models(core, "model")
write.csv(by_model, file.path(out_dir, "recovery_table.csv"), row.names = FALSE)
print(by_model, row.names = FALSE)

ggsave(file.path(fig_dir, "recovery_summary.png"),
       plot_metrics(metric_long(by_model, "model"),
                    title = "Case 3: overlapping phases"),
       width = 8, height = 3.4, dpi = 300, bg = "white")

ggsave(file.path(chk_dir, "recovery_by_n.png"),
       plot_by_n(core, "Case 3: effect of dataset size"),
       width = 8, height = 4.2, dpi = 300, bg = "white")

factor_rows <- r[r$sweep == "factor", ]
sweep_figure(factor_rows[factor_rows$assign_p == 0.5, ], "overlap", "overlap",
             "Case 3: overlap (assign_p = 0.5)", out_dir, fig_dir)
sweep_figure(factor_rows[factor_rows$overlap == 0.5, ], "assign_p", "assign_p",
             "Case 3: assign_p (overlap = 0.5)", out_dir, chk_dir)

# One simulated dataset per overlap level
overlaps <- c(0, 0.25, 0.5)
anat <- do.call(rbind, lapply(overlaps, function(ov) {
  sim <- simulate_overlap(70, 8, 0.02, 1, K = 6, alpha_conc = 1.5,
                          overlap = ov, assign_p = 0.2, seed = 3)
  sim$panel <- paste("overlap", ov)
  sim$Point_date <- (sim$Start_date + sim$End_date) / 2
  sim
}))
anat$panel <- factor(anat$panel, levels = paste("overlap", overlaps))
trend <- data.frame(panel = levels(anat$panel), intercept = 8, slope = 0.02)
ggsave(file.path(fig_dir, "dataset_anatomy.png"),
       dataset_anatomy(anat, trend = trend),
       width = 9, height = 3.6, dpi = 300, bg = "white")

cat("\nfits with rhat > 1.01:", sum(r$max_rhat > 1.01, na.rm = TRUE),
    "\nfits with a divergent transition:", sum(r$n_divergent > 0), "\n")

fp <- false_positive_rate(r)
if (!is.null(fp)) {
  cat("\nfalse-positive rate (zero-slope cells, 90% interval excludes 0):\n")
  print(fp)
}

cat("wrote figures to", fig_dir, "\n")
