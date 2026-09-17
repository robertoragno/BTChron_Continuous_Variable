# Builds the Case 4 design grid and writes it to data/design.csv.
# 03_recovery_study.R reads that file and re-derives nothing.
#
#   core   dataset size crossed with random / zero slope, broad periods up to
#          3 phases
#   merge  broad periods up to 2, 3 or 4 phases, at N = 200

library(here)

# Simulating script
source(here("Simulations", "Sim_Case4_MergedPhases", "scripts", "simulate.R"))

# add_true_trend(): an id, a seed and a true intercept, slope and sigma per dataset
source(here("Simulations", "shared", "scripts", "design.R"))

# ---- settings ----
PERIOD_START <- 100
PERIOD_END   <- 900
MERGE_LEVELS <- c(2, 3, 4)   # largest number of fine phases in a broad period
MERGE_REF    <- 3
K_RANGE      <- 5:10         # fine phases per dataset
N_LEVELS     <- c(50, 100, 200, 500)
N_FACTOR     <- 200
N_REP        <- 100          # datasets per setting
OUTPUT_CSV   <- Sys.getenv("DESIGN_CASE4_OUT",
                           here("Simulations", "Sim_Case4_MergedPhases",
                                "data", "design.csv"))
# ------------------

set.seed(2026)

core <- expand.grid(rep = seq_len(N_REP), N = N_LEVELS,
                    slope_condition = c("random", "zero"),
                    stringsAsFactors = FALSE)
core$merge_max <- MERGE_REF
core$sweep     <- "core"

merge_sweep <- expand.grid(rep = seq_len(N_REP), merge_max = MERGE_LEVELS)
merge_sweep$N               <- N_FACTOR
merge_sweep$slope_condition <- "random"
merge_sweep$sweep           <- "merge"

keep   <- c("sweep", "rep", "N", "slope_condition", "merge_max")
design <- rbind(core[, keep], merge_sweep[, keep])

design <- add_true_trend(design)
design$K          <- sample(K_RANGE, nrow(design), replace = TRUE)
design$alpha_conc <- 10^(-1 + 2 * rbeta(nrow(design), 2, 3.5))

design$period_start <- PERIOD_START
design$period_end   <- PERIOD_END

# Calendar axis passed to Stan. Merged windows never leave the period, so no
# padding is needed.
design$time_ref_min   <- PERIOD_START
design$time_ref_range <- PERIOD_END - PERIOD_START

# Generator checks: every window contains its true date, stays in the period,
# and each dataset mixes single phases with merged ones
for (mm in MERGE_LEVELS) {
  d <- simulate_merged(4000, 8, 0.02, 1, K = 8, alpha_conc = 1, merge_max = mm,
                       period_start = PERIOD_START, period_end = PERIOD_END, seed = 1)
  stopifnot(all(d$True_date >= d$Start_date & d$True_date <= d$End_date),
            all(d$Start_date >= PERIOD_START & d$End_date <= PERIOD_END),
            min(d$Phases) == 1, max(d$Phases) == mm)
}
cat("generator checks ok\n")

dir.create(dirname(OUTPUT_CSV), showWarnings = FALSE, recursive = TRUE)
write.csv(design, OUTPUT_CSV, row.names = FALSE)

cat(sprintf("wrote %d datasets to %s\n", nrow(design), OUTPUT_CSV))
print(table(design$merge_max, design$sweep))
