# Case 4 tables and figures from output/recovery_results.csv.
# The metrics and their plots live in shared/scripts/recovery_summary.R; this
# file only says what Case 4 is faceted by, which is broad-period width.

library(here)
library(ggplot2)

source(here("Simulations", "shared", "scripts", "recovery_summary.R"))
source(here("Simulations", "shared", "scripts", "partition.R"))
source(here("Simulations", "Sim_Case4_MergedPhases", "scripts", "simulate.R"))

results_csv <- Sys.getenv("RECOVERY_OUT",
                          here("Simulations", "Sim_Case4_MergedPhases", "output",
                               "recovery_results.csv"))
fig_dir <- here("Simulations", "Sim_Case4_MergedPhases", "figures")
out_dir <- here("Simulations", "Sim_Case4_MergedPhases", "output")

# median is not fitted for flat windows; see 00_check_median_identity.R
r <- read_recovery(results_csv,
                   merge_median = !("median" %in% read.csv(results_csv)$model))

by_model <- summarise_models(r, "model")
write.csv(by_model, file.path(out_dir, "recovery_table.csv"), row.names = FALSE)
print(by_model, row.names = FALSE)

ggsave(file.path(fig_dir, "recovery_summary.png"),
       plot_metrics(metric_long(by_model, "model"),
                    title = "Case 4: one model across mixed dating resolution"),
       width = 8, height = 3.4, dpi = 300, bg = "white")

core <- r[r$sweep == "core" & r$slope_condition == "random", ]
if (nrow(core) > 0)
  ggsave(file.path(fig_dir, "recovery_by_n.png"),
         plot_by_n(core, "Case 4: precision improves with N, accuracy does not"),
         width = 8, height = 4.2, dpi = 300, bg = "white")

# How coarse the broad periods get. Wider merges should widen the intervals and
# inflate sigma for the point-date without moving the slope.
sweep_figure(r[r$sweep %in% c("core", "merge") & r$slope_condition == "random", ],
             "merge_max", "merge", "Case 4: effect of broad-period width",
             out_dir, fig_dir)

# One example dataset per merge width, in the archived Sim_Linear idiom: the fine
# phases as alternating bands, each find's broad period the segment spanning the
# run of them it was dated to, the true date against the trend. Everything but
# merge_max is held fixed across the panels, so the only thing that changes from
# one to the next is how coarse the broad periods are.
design <- read.csv(here("Simulations", "Sim_Case4_MergedPhases", "data",
                        "design.csv"))
merge_levels <- sort(unique(design$merge_max))
t_min <- design$t_min[1]; t_max <- design$t_max[1]
anat_rows <- lapply(merge_levels, function(mm) {
  sim <- simulate_merged(70, intercept = 8, slope = 0.02, sigma = 2, K = 8,
                         alpha_conc = 1, merge_max = mm,
                         t_min = t_min, t_max = t_max, seed = 4)
  lab <- sprintf("merge up to %d phases", mm)
  b   <- attr(sim, "bounds")
  list(anat = data.frame(panel = lab,
                         Start_date = sim$Start_date, End_date = sim$End_date,
                         True_date = sim$True_date, Value = sim$Value,
                         Point_date = (sim$Start_date + sim$End_date) / 2),
       ribbon = data.frame(panel = lab, xmin = b[-length(b)], xmax = b[-1]),
       trend  = data.frame(panel = lab, intercept = 8, slope = 0.02))
})
lvl  <- sprintf("merge up to %d phases", merge_levels)
bind <- function(k) do.call(rbind, lapply(anat_rows, `[[`, k))
anat <- bind("anat");   anat$panel <- factor(anat$panel, levels = lvl)
rib  <- bind("ribbon"); rib$panel  <- factor(rib$panel,  levels = lvl)
trd  <- bind("trend");  trd$panel  <- factor(trd$panel,  levels = lvl)
ggsave(file.path(fig_dir, "dataset_anatomy.png"),
       dataset_anatomy(anat, ribbon = rib, trend = trd,
                       title = "Case 4: one dataset per broad-period width"),
       width = 9, height = 3.6, dpi = 300, bg = "white")

fp <- false_positive_rate(r)
if (!is.null(fp)) {
  cat("\nfalse-positive rate (zero-slope cells, 90% interval excludes 0):\n")
  print(fp)
}
cat("wrote figures to", fig_dir, "\n")
