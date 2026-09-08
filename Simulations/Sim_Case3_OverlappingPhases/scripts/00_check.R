# What the Case 3 generator produces, before anything is fitted to it.
#
# Four things, in order of how much rides on them:
#
#   1. overlap = 0 reproduces the archived partition_timeline() year for year.
#      This is the unit check the plan (section 12) gates the case on.
#   2. Every recorded window contains its own true date, boundaries stay on
#      whole years, and the phases still cover [t_min, t_max] with no gap.
#   3. The midpoint slope ratio by (overlap, assign_p), OLS only, no Stan. This
#      was meant to decide the design. It decided against it: the ratio is 1.00
#      at every level of both factors (run 2026-09-07). assign_p below 0.5 does
#      thin each phase's early and late shared strips unequally, but by a
#      near-constant number of years for every phase, which moves the regression
#      intercept and not the slope. So overlapping phases are a CALIBRATION
#      problem here, not a slope-bias one. See the README and section 5.1 of the
#      plan; the case is paused pending a mechanism whose displacement varies
#      with calendar position.
#   4. The retained-date density inside one wide phase: the picture of the
#      per-phase Berkson violation that (3) shows does not reach the slope. The
#      red retained-mean line pulls left of the midpoint as assign_p drops.
#
# Stepped through or run whole with Rscript, same result.

library(here)
library(ggplot2)

source(here("Simulations", "shared", "scripts", "partition.R"))
source(here("Simulations", "Sim_Case3_OverlappingPhases", "scripts", "simulate.R"))
source(here("Simulations", "shared", "scripts", "recovery_summary.R"))

# ---- settings to play with ----
T_MIN     <- 100      # matches the archived Sim_Linear span
T_MAX     <- 900
OVERLAPS  <- c(0, 0.25, 0.5)
ASSIGN_PS <- c(0.5, 0.35, 0.2)   # 0.5 = coin flip; lower leans dates later
N_BIG     <- 40000
N_REP     <- 200      # datasets per (overlap, assign_p) cell in the ratio check
TRUE_SLOPE <- 0.02
# -------------------------------

figure_path <- function(name)
  here("Simulations", "Sim_Case3_OverlappingPhases", "figures", name)
panel_theme <- theme_classic(base_size = 11) +
  theme(plot.title = element_text(size = 11, face = "bold"),
        strip.background = element_blank(),
        strip.text = element_text(face = "bold"),
        legend.position = "top", legend.title = element_blank())

midpoint <- function(d) (d$Start_date + d$End_date) / 2

# 1. overlap = 0 must reproduce the archived function exactly. The archived file
# reads TMIN / TMAX as globals off its own top level, so it is sourced into its
# own environment and called from there - a plain source() would clobber the new
# partition_timeline() with the old one.
archived <- new.env()
sys.source(here("Simulations", "archive", "Sim_Linear", "scripts", "simulate.R"),
           envir = archived)

for (s in 1:20) {
  set.seed(s); a <- archived$partition_timeline(300, 6, 1)
  set.seed(s); b <- partition_timeline(300, 6, 1, overlap = 0,
                                       t_min = T_MIN, t_max = T_MAX)
  stopifnot(identical(a$Start_date, b$Start_date),
            identical(a$End_date,   b$End_date),
            identical(a$True_date,  b$True_date))
}
cat("overlap = 0 reproduces the archived partition_timeline() on 20 seeds\n")

# 2. Containment, whole-year boundaries, gap-free coverage.
for (ov in OVERLAPS) for (ap in ASSIGN_PS) {
  d <- simulate_overlap(5000, intercept = 8, slope = TRUE_SLOPE, sigma = 1,
                        K = 6, alpha_conc = 1, overlap = ov, assign_p = ap,
                        t_min = T_MIN, t_max = T_MAX, seed = 1)
  b <- attr(d, "bounds")
  stopifnot(all(d$True_date >= d$Start_date & d$True_date <= d$End_date),
            all(d$Start_date == round(d$Start_date)),
            all(d$End_date == round(d$End_date)),
            b[1] == T_MIN, b[length(b)] == T_MAX,
            all(diff(b) > 0))
}
cat("containment, whole-year boundaries and gap-free coverage hold for every",
    "overlap x assign_p\n")

# 3. Midpoint slope ratio, and the true-date-on-midpoint slope that predicts it.
# Both from OLS, so this is the dating structure and nothing the Stan models do.
ratio_grid <- expand.grid(overlap = OVERLAPS, assign_p = ASSIGN_PS)
ratio_grid <- ratio_grid[!(ratio_grid$overlap == 0 & ratio_grid$assign_p != 0.5), ]

ratio_rows <- lapply(seq_len(nrow(ratio_grid)), function(g) {
  ov <- ratio_grid$overlap[g]; ap <- ratio_grid$assign_p[g]
  est <- vapply(seq_len(N_REP), function(i) {
    s <- simulate_overlap(200, intercept = 8, slope = TRUE_SLOPE,
                          sigma = 1, K = sample(3:10, 1),
                          alpha_conc = 10^(-1 + 2 * rbeta(1, 2, 3.5)),
                          overlap = ov, assign_p = ap,
                          t_min = T_MIN, t_max = T_MAX, seed = 1000 + i)
    m <- midpoint(s)
    c(coef(lm(s$Value ~ m))[2] / TRUE_SLOPE,
      coef(lm(s$True_date ~ m))[2])
  }, numeric(2))
  # A very lumpy partition can drop every date into one phase, leaving the
  # midpoint constant and its OLS slope undefined; those reps drop out.
  data.frame(overlap = ov, assign_p = ap,
             midpoint_slope_ratio = median(est[1, ], na.rm = TRUE),
             truedate_on_midpoint = median(est[2, ], na.rm = TRUE),
             n_used = sum(!is.na(est[1, ])))
})
ratio_tab <- do.call(rbind, ratio_rows)
cat("\nmidpoint slope ratio (< 1 means the trend comes out flattened):\n")
print(ratio_tab, row.names = FALSE, digits = 3)
cat("\ntruedate_on_midpoint is the slope of lm(True_date ~ midpoint) - its",
    "deficit\nfrom 1 is the attenuation the midpoint model should show.\n")

# 4. Retained-date density inside one wide phase, at the widest overlap.
# Fix the periodisation so the same phase is compared across assign_p.
set.seed(7)
probe_K <- 5
probe_bounds <- attr(simulate_overlap(10, 8, TRUE_SLOPE, 1, K = probe_K,
                                      alpha_conc = 2, overlap = 0,
                                      t_min = T_MIN, t_max = T_MAX,
                                      seed = 7), "bounds")
probe_phase <- which.max(diff(probe_bounds)[-c(1, probe_K)]) + 1  # widest interior

dens <- do.call(rbind, lapply(ASSIGN_PS, function(ap) {
  d <- simulate_overlap(N_BIG, 8, TRUE_SLOPE, 1, K = probe_K, alpha_conc = 2,
                        overlap = 0.5, assign_p = ap,
                        t_min = T_MIN, t_max = T_MAX, seed = 7)
  ph <- findInterval(midpoint(d), probe_bounds, rightmost.closed = TRUE,
                     all.inside = TRUE)
  keep <- d[ph == probe_phase, ]
  data.frame(assign_p = sprintf("assign_p = %.2f", ap),
             True_date = keep$True_date,
             mid = mean(c(keep$Start_date[1], keep$End_date[1])),
             retained_mean = mean(keep$True_date))
}))

marks <- unique(dens[, c("assign_p", "mid", "retained_mean")])
p <- ggplot(dens, aes(True_date)) +
  geom_histogram(bins = 40, fill = "grey75", colour = "grey30",
                 linewidth = 0.2) +
  geom_vline(data = marks, aes(xintercept = mid),
             linetype = "dashed", colour = "grey30") +
  geom_vline(data = marks, aes(xintercept = retained_mean),
             colour = "#780000", linewidth = 0.7) +
  facet_wrap(~ assign_p, ncol = 1) +
  labs(title = "Case 3: retained true dates in one wide phase (overlap 0.5)",
       subtitle = paste("dashed: window midpoint;",
                        "red: mean of the dates the phase kept"),
       x = "True date (yr)", y = "Count") +
  panel_theme

ggsave(figure_path("check_depletion.png"), p,
       width = 7.5, height = 6.5, dpi = 300, bg = "white")
cat("\nwrote", figure_path("check_depletion.png"), "\n")

# 5. One example dataset per overlap level, in the archived Sim_Linear idiom.
# No phase ribbon here: the phase windows overlap by construction, so alternating
# bands would paint over each other. The assigned-window segments carry it - at
# overlap 0.5 the segments from adjacent phases share calendar range.
anat <- do.call(rbind, lapply(OVERLAPS, function(ov) {
  d <- simulate_overlap(70, intercept = 8, slope = TRUE_SLOPE, sigma = 1,
                        K = 6, alpha_conc = 1.5, overlap = ov, assign_p = 0.5,
                        t_min = T_MIN, t_max = T_MAX, seed = 3)
  data.frame(panel = factor(sprintf("overlap = %.2f", ov),
                            levels = sprintf("overlap = %.2f", OVERLAPS)),
             Start_date = d$Start_date, End_date = d$End_date,
             True_date = d$True_date, Value = d$Value,
             Point_date = midpoint(d))
}))
trend <- data.frame(panel = factor(sprintf("overlap = %.2f", OVERLAPS),
                                   levels = sprintf("overlap = %.2f", OVERLAPS)),
                    intercept = 8, slope = TRUE_SLOPE)
ggsave(figure_path("dataset_anatomy.png"),
       dataset_anatomy(anat, trend = trend,
                       title = "Case 3: one dataset per overlap level"),
       width = 9, height = 3.6, dpi = 300, bg = "white")
cat("wrote", figure_path("dataset_anatomy.png"), "\n")
