# Builds the Case 3 design grid and writes it to data/design.csv.
# 03_recovery_study.R reads that file and re-derives nothing.
#
#   core    dataset size crossed with random / zero slope, at overlap 0.5 and
#           assign_p 0.2
#   factor  overlap crossed with assign_p, at N = 200. assign_p does nothing
#           when phases do not overlap, so overlap 0 appears once.

library(here)

# Simulating script
source(here("Simulations", "Sim_Case3_OverlappingPhases", "scripts", "simulate.R"))

# add_true_trend(): an id, a seed and a true intercept, slope and sigma per dataset
source(here("Simulations", "shared", "scripts", "design.R"))

# ---- settings ----
PERIOD_START <- 100
PERIOD_END   <- 900
OVERLAPS     <- c(0, 0.25, 0.5)
ASSIGN_PS    <- c(0.5, 0.35, 0.2)   # chance a find in a shared strip goes to the earlier phase
REF_OVERLAP  <- 0.5
REF_ASSIGN_P <- 0.2
K_RANGE      <- 5:10                # phases per dataset
N_LEVELS     <- c(50, 100, 200, 500)
N_FACTOR     <- 200
N_REP        <- 100                 # datasets per setting
OUTPUT_CSV   <- Sys.getenv("DESIGN_CASE3_OUT",
                           here("Simulations", "Sim_Case3_OverlappingPhases",
                                "data", "design.csv"))
# ------------------

set.seed(2026)

# Phase windows run past the period by up to overlap / 2 of the first or last
# phase. Grid padding: the widest stretch seen over 20000 random sets of
# phases, plus 25%. (sigma = 0 so no noise is drawn; this keeps the random
# numbers, and so design.csv, the same as in the original run.)
reach <- replicate(20000, {
  K          <- sample(K_RANGE, 1)
  alpha_conc <- 10^(-1 + 2 * rbeta(1, 2, 3.5))
  sim        <- simulate_overlap(1, 0, 0, 0, K, alpha_conc, overlap = max(OVERLAPS),
                                 period_start = PERIOD_START, period_end = PERIOD_END)
  lengths    <- diff(attr(sim, "bounds"))
  max(OVERLAPS) * max(lengths[1], lengths[K]) / 2
})
pad <- ceiling(1.25 * max(reach))
cat("grid padding", pad, "yr\n")

core <- expand.grid(rep = seq_len(N_REP), N = N_LEVELS,
                    slope_condition = c("random", "zero"),
                    stringsAsFactors = FALSE)
core$overlap  <- REF_OVERLAP
core$assign_p <- REF_ASSIGN_P
core$sweep    <- "core"

settings <- expand.grid(overlap = OVERLAPS, assign_p = ASSIGN_PS)
settings <- settings[settings$overlap > 0 | settings$assign_p == 0.5, ]
factor_sweep <- settings[rep(seq_len(nrow(settings)), each = N_REP), ]
factor_sweep$rep             <- rep(seq_len(N_REP), nrow(settings))
factor_sweep$N               <- N_FACTOR
factor_sweep$slope_condition <- "random"
factor_sweep$sweep           <- "factor"

keep   <- c("sweep", "rep", "N", "slope_condition", "overlap", "assign_p")
design <- rbind(core[, keep], factor_sweep[, keep])
rownames(design) <- NULL

design <- add_true_trend(design)
design$K          <- sample(K_RANGE, nrow(design), replace = TRUE)
design$alpha_conc <- 10^(-1 + 2 * rbeta(nrow(design), 2, 3.5))

design$period_start <- PERIOD_START
design$period_end   <- PERIOD_END

# Calendar axis passed to Stan, the same for every dataset and model
design$time_ref_min   <- PERIOD_START - pad
design$time_ref_range <- (PERIOD_END + pad) - design$time_ref_min

# Generator checks: every window contains its true date and fits in the grid
for (ov in OVERLAPS) {
  d <- simulate_overlap(4000, 8, 0.02, 1, K = sample(K_RANGE, 1), alpha_conc = 0.1,
                        overlap = ov, assign_p = REF_ASSIGN_P,
                        period_start = PERIOD_START, period_end = PERIOD_END, seed = 1)
  stopifnot(all(d$True_date >= d$Start_date & d$True_date <= d$End_date),
            all(d$Start_date >= PERIOD_START - pad & d$End_date <= PERIOD_END + pad))
}
cat("generator checks ok\n")

dir.create(dirname(OUTPUT_CSV), showWarnings = FALSE, recursive = TRUE)
write.csv(design, OUTPUT_CSV, row.names = FALSE)

cat(sprintf("wrote %d datasets to %s\n", nrow(design), OUTPUT_CSV))
print(table(design$overlap, design$assign_p, design$sweep))
