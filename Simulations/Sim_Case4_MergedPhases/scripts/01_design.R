# Builds the Case 4 design grid and writes it to data/design.csv.
# 03_recovery_study.R reads that file and re-derives nothing, so a partial
# re-run cannot draw a different design by accident. Edit the settings below and
# re-run to change the study.
#
# Two blocks:
#   core   sample size crossed with the slope condition, at the reference merge
#          width. The N sweep reads precision against sample size; the zero-slope
#          cells give the false-positive rate.
#   merge  how coarse the broad periods get, at fixed N. The reference width is
#          not repeated here; the matching core cell stands in for it.
#
# A few generator checks run at the end, on freshly simulated datasets, before
# anything is written.

library(here)

source(here("Simulations", "shared", "scripts", "partition.R"))
source(here("Simulations", "Sim_Case4_MergedPhases", "scripts", "simulate.R"))
source(here("Simulations", "shared", "scripts", "design.R"))

# ---- settings ----
T_MIN        <- 100     # same span as Case 3, so the phase structure is comparable
T_MAX        <- 900
MERGE_LEVELS <- c(2, 3, 4)
MERGE_REF    <- 3       # reference broad-period width for the core sweep
K_RANGE      <- 5:10    # fine phases per dataset, drawn per dataset
N_LEVELS     <- c(50, 100, 200, 500)
N_FACTOR     <- 200     # sample size the merge sweep runs at
N_REP        <- 100     # datasets per cell
OUTPUT_CSV   <- Sys.getenv("DESIGN_CASE4_OUT",
                           here("Simulations", "Sim_Case4_MergedPhases", "data",
                                "design.csv"))
# ------------------

set.seed(2026)

core <- expand.grid(rep = seq_len(N_REP), N = N_LEVELS,
                    slope_condition = c("random", "zero"),
                    stringsAsFactors = FALSE)
core$merge_max <- MERGE_REF
core$sweep     <- "core"

merge_sweep <- expand.grid(rep = seq_len(N_REP), merge_max = MERGE_LEVELS,
                           stringsAsFactors = FALSE)
merge_sweep$N               <- N_FACTOR
merge_sweep$slope_condition <- "random"
merge_sweep$sweep           <- "merge"

keep   <- c("sweep", "rep", "N", "slope_condition", "merge_max")
design <- rbind(core[, keep], merge_sweep[, keep])

design <- add_nuisance(design)
design$K          <- sample(K_RANGE, nrow(design), replace = TRUE)
design$alpha_conc <- 10^(-1 + 2 * rbeta(nrow(design), 2, 3.5))

# Recorded so the study and any re-analysis use the same generator settings.
design$t_min <- T_MIN
design$t_max <- T_MAX

# The merged windows never leave [T_MIN, T_MAX], so the calendar axis needs no
# padding. Passed to Stan as data, identical for every dataset and model, so the
# prior on the normalised slope means the same thing everywhere.
design$time_ref_min   <- T_MIN
design$time_ref_range <- T_MAX - T_MIN

# ---- generator checks ----
for (mm in MERGE_LEVELS) {
  d <- simulate_merged(4000, intercept = 8, slope = 0.02, sigma = 1,
                       K = 8, alpha_conc = 1, merge_max = mm,
                       t_min = T_MIN, t_max = T_MAX, seed = 1)
  span <- attr(d, "span")
  stopifnot(
    all(d$True_date >= d$Start_date & d$True_date <= d$End_date),
    all(d$Start_date == round(d$Start_date)),
    all(d$End_date == round(d$End_date)),
    all(d$Start_date >= T_MIN & d$End_date <= T_MAX),
    all(span >= 1 & span <= mm),
    min(span) == 1,          # some finds keep a fine phase
    (mm == 1) || max(span) > 1  # and some are merged
  )
}
cat("generator checks passed for merge_max", paste(MERGE_LEVELS, collapse = ", "), "\n")

write.csv(design, OUTPUT_CSV, row.names = FALSE)
cat(sprintf("wrote %d datasets to %s\n", nrow(design), OUTPUT_CSV))
print(table(design$sweep, design$merge_max))
