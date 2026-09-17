# Case 2 tables and figures from output/recovery_results.csv.
# Metrics and plots live in shared/scripts/recovery_summary.R; this file only
# says which sweeps Case 2 has.

library(here)

source(here("Simulations", "shared", "scripts", "recovery_summary.R"))

results_csv <- Sys.getenv("RECOVERY_OUT",
                          here("Simulations", "Sim_Case2_Typochronology",
                               "output", "recovery_results.csv"))
fig_dir <- here("Simulations", "Sim_Case2_Typochronology", "figures")
chk_dir <- file.path(fig_dir, "checks")
out_dir <- here("Simulations", "Sim_Case2_Typochronology", "output")
dir.create(chk_dir, showWarnings = FALSE)

# median is not fitted for flat windows; see 00_check_median_identity.R
r <- read_recovery(results_csv, merge_median = !("median" %in% read.csv(results_csv)$model))

by_model <- summarise_models(r[r$sweep == "core" & r$slope_condition == "random", ],
                             "model")
write.csv(by_model, file.path(out_dir, "recovery_table.csv"), row.names = FALSE)
print(by_model, row.names = FALSE)

p1 <- plot_metrics(metric_long(by_model, "model"),
                   title = "Case 2: finds with mixed dating precision")
ggsave(file.path(fig_dir, "recovery_summary.png"), p1,
       width = 8, height = 3.4, dpi = 300, bg = "white")

core <- r[r$sweep == "core" & r$slope_condition == "random", ]
ggsave(file.path(chk_dir, "recovery_by_n.png"),
       plot_by_n(core, "Case 2: effect of dataset size"),
       width = 8, height = 4.2, dpi = 300, bg = "white")

# Every factor sweep is read against the reference cell, which sits in core
reference <- r[r$sweep == "core" & r$slope_condition == "random" & r$N == 200, ]
sweep_rows <- function(name) rbind(reference, r[r$sweep == name, ])

sweep_figure(sweep_rows("prop"), "prop_coarse_samples", "prop_coarse",
             "Case 2: proportion of coarsely dated finds", out_dir, fig_dir)
sweep_figure(sweep_rows("width"), "coarse_frac", "coarse_width",
             "Case 2: coarse window length, share of the period", out_dir, fig_dir)
sweep_figure(sweep_rows("deposition"), "growth_ratio", "deposition",
             "Case 2: finds denser at the end of the period", out_dir, chk_dir)
sweep_figure(sweep_rows("precision"), "precision_trend", "precision",
             "Case 2: coarsely dated finds more common early", out_dir, chk_dir)

# Noise ratio: noise compared with the total change the trend produces across
# the period. Not a design factor, worked out per dataset after the fact.
core$noise_ratio <- core$sigma / (abs(core$slope) * (core$period_end - core$period_start))
core$noise_band  <- cut(core$noise_ratio, c(0, 0.1, 0.3, 1, Inf),
                        labels = c("< 0.1", "0.1-0.3", "0.3-1", "> 1"))
by_noise <- summarise_models(core, c("model", "noise_band"))
write.csv(by_noise, file.path(out_dir, "recovery_by_noise.csv"), row.names = FALSE)

fp <- false_positive_rate(r)
if (!is.null(fp)) {
  cat("\nfalse-positive rate (zero-slope cells, 90% interval excludes 0):\n")
  print(fp)
}

cat("wrote figures to", fig_dir, "\n")
