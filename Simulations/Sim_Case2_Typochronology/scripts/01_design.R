# Builds the Case 2 design grid and writes it to data/design.csv.
# 03_recovery_study.R reads that file and re-derives nothing.
#
# One core sweep plus one sweep per factor. Each factor sweep changes one thing
# and keeps the rest at the reference (N = 200, half the finds coarsely dated,
# coarse windows 25% of the period, even deposition, no precision trend, random
# slope).
#
#   core        dataset size crossed with random / zero slope
#   prop        proportion of coarsely dated finds
#   width       length of a coarse window compared with the period
#   deposition  finds denser at the end of the period
#   precision   coarsely dated finds more common early than late

library(here)

# Simulating script
source(here("Simulations", "Sim_Case2_Typochronology", "scripts", "simulate.R"))

# add_true_trend(): an id, a seed and a true intercept, slope and sigma per dataset
source(here("Simulations", "shared", "scripts", "design.R"))

# ---- settings ----
PERIOD_START    <- 100   # Define the period of interest
PERIOD_END      <- 900
FINE_FRAC       <- 0.0625 # 50 yr of an 800-yr period
COARSE_LEVELS   <- c(0.10, 0.25, 0.45) # typical coarse window, share of the period
COARSE_REF      <- 0.25
N_LEVELS        <- c(50, 100, 200, 500)
N_REF           <- 200
# List of values tried for the proportion of coarsely dated samples in a dataset
PROP_LEVELS     <- c(0, 0.25, 0.5, 0.75, 1)
PROP_REF        <- 0.5
GROWTH_RATIO    <- 4    # finds 4 times denser at the end of the period
# Are coarsely dated finds more common early than late? 0 = no, 1 = yes
# The proportion is the same, what changes is the position of coarse samples in the timeline
PRECISION_TREND <- 0.6 # (80% of coarse samples in the first 20% of the period)
N_REP           <- 200  # datasets per setting
OUTPUT_CSV      <- Sys.getenv("DESIGN_CASE2_OUT",
                              here("Simulations", "Sim_Case2_Typochronology",
                                   "data", "design.csv"))
# ------------------

set.seed(2026)

reference <- data.frame(N = N_REF, slope_condition = "random",
                        prop_coarse_samples = PROP_REF, coarse_frac = COARSE_REF,
                        growth_ratio = 1, precision_trend = 0)

# Repeat the reference, then overwrite the column(s) a sweep varies
sweep_block <- function(name, varied) {
  grid <- expand.grid(c(list(rep = seq_len(N_REP)), varied),
                      stringsAsFactors = FALSE)
  for (col in setdiff(names(reference), names(varied))) grid[[col]] <- reference[[col]]
  grid$sweep <- name
  grid
}

design <- rbind(
  sweep_block("core",       list(N = N_LEVELS, slope_condition = c("random", "zero"))),
  sweep_block("prop",       list(prop_coarse_samples = setdiff(PROP_LEVELS, PROP_REF))),
  sweep_block("width",      list(coarse_frac = setdiff(COARSE_LEVELS, COARSE_REF))),
  sweep_block("deposition", list(growth_ratio = GROWTH_RATIO)),
  sweep_block("precision",  list(precision_trend = PRECISION_TREND))
)

design <- add_true_trend(design)

design$period_start <- PERIOD_START
design$period_end   <- PERIOD_END
design$fine_frac    <- FINE_FRAC

# Calendar axis passed to Stan, the same for every dataset and model. Windows
# can run past the period by up to the longest window in the whole design.
pad <- ceiling(1.2 * max(COARSE_LEVELS) * (PERIOD_END - PERIOD_START))
design$time_ref_min   <- PERIOD_START - pad
design$time_ref_range <- (PERIOD_END + pad) - design$time_ref_min

# Generator checks on large simulated datasets
check <- function(...) simulate_typo(40000, 8, 0.02, 1, ..., fine_frac = FINE_FRAC,
                                     period_start = PERIOD_START,
                                     period_end = PERIOD_END, seed = 1)
d <- check(prop_coarse_samples = 0.5, coarse_frac = max(COARSE_LEVELS))
stopifnot(all(d$True_date >= d$Start_date & d$True_date <= d$End_date),
          abs(mean(d$Coarse) - 0.5) < 0.01,
          abs(mean(d$True_date) - (PERIOD_START + PERIOD_END) / 2) < 5,
          min(d$Start_date) >= design$time_ref_min[1],
          max(d$End_date) <= design$time_ref_min[1] + design$time_ref_range[1])

d <- check(growth_ratio = GROWTH_RATIO)
stopifnot(mean(d$True_date) > (PERIOD_START + PERIOD_END) / 2 + 50)

d <- check(precision_trend = PRECISION_TREND)
early <- d$True_date < PERIOD_START + 200
late  <- d$True_date > PERIOD_END - 200
stopifnot(mean(d$Coarse[early]) > mean(d$Coarse[late]) + 0.3)
cat("generator checks ok\n")

dir.create(dirname(OUTPUT_CSV), showWarnings = FALSE, recursive = TRUE)
write.csv(design, OUTPUT_CSV, row.names = FALSE)

cat(sprintf("wrote %d datasets to %s\n", nrow(design), OUTPUT_CSV))
print(table(design$sweep))
