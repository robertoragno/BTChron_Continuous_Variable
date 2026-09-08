# Builds the Case 1 design grid and writes it to data/design_case1.csv.
#
# Same shape as 01_design_case2.R: written before any fitting, and every
# downstream script reads this file rather than rebuilding the grid.
#
# Two blocks:
#   core    sample size crossed with the slope condition, at the reference
#           lab error and on the plateau. The N sweep reads precision against
#           sample size; the zero-slope cells give the false-positive rate.
#   factor  lab error crossed with calendar window, at fixed N.
#   deposition  growth in find density, at both windows and the reference lab
#           error. Uniform is not repeated here - the matching `factor` cells at
#           the reference lab error are its reference, the way Case 2's core
#           cells stand in for its "even" resolution.

library(here)

source(here("Simulations", "shared", "scripts", "weight_rows.R"))
source(here("Simulations", "shared", "scripts", "design.R"))
source(here("Simulations", "Sim_Case1_Radiocarbon", "scripts", "simulate.R"))

# ---- settings to play with ----
# Two windows of the same length, named for the shape of IntCal20 across them.
# Over the plateau the curve is nearly flat, so the calibrated posterior is
# widest and most multimodal and the median and the HPD envelope midpoint
# disagree most (45 yr on average, 98 at worst, against 9 and 25 off it). Over
# the steep window the curve rises monotonically, posteriors are narrow and
# single-humped, and every summary of a date says much the same thing.
WINDOWS      <- list(plateau = c(-800, -400),
                     steep   = c(-1600, -1200))
LAB_ERRORS   <- c(15, 30, 50)
REF_ERROR    <- 30
REF_WINDOW   <- "plateau"
GRID_STEP    <- 5       # candidate-year spacing; see the equivalence check
N_LEVELS     <- c(50, 100, 200, 500)
N_FACTOR     <- 200     # sample size the factor sweep runs at
# Deposition: how many times denser finds are at the end of the window than at
# the start. 1 is uniform. A ratio, not a rate in years, so it means the same
# thing whatever the window length.
GROWTH_RATIO <- 4
N_REP        <- 100     # datasets per cell
OUTPUT_CSV   <- Sys.getenv("DESIGN_CASE1_OUT",
                           here("Simulations", "Sim_Case1_Radiocarbon", "data",
                                "design.csv"))
# -------------------------------

set.seed(2026)

# One padding constant for the whole case, not one per cell.
#
# This is where the "no window prior" decision becomes a number. The model is
# not told that true dates come from a 400-yr window (see the DECISION note in
# weight_rows.R), so every date's calibrated posterior has to fit on the grid
# whole. Clipping it would impose that window assumption by the back door and
# drag edge dates inward. The padding is measured, not guessed: the widest
# calibrated support that actually occurs.
#
# Taking the maximum across every window and lab error, rather than per cell,
# keeps time_ref_range identical everywhere: the windows are the same length, so
# normal(0, 10) on beta means the same thing on the calendar scale in every
# cell, and cells stay comparable.
pad <- max(sapply(WINDOWS, function(w)
  sapply(LAB_ERRORS, function(e) c14_grid_pad(w, e))))

cat(sprintf("grid padding %d yr (widest calibrated support over %d windows x %d lab errors)\n",
            pad, length(WINDOWS), length(LAB_ERRORS)))

core <- expand.grid(rep = seq_len(N_REP), N = N_LEVELS,
                    slope_condition = c("random", "zero"),
                    stringsAsFactors = FALSE)
core$lab_error <- REF_ERROR
core$window    <- REF_WINDOW
core$sweep     <- "core"

factor_sweep <- expand.grid(rep = seq_len(N_REP), lab_error = LAB_ERRORS,
                            window = names(WINDOWS), stringsAsFactors = FALSE)
factor_sweep$N               <- N_FACTOR
factor_sweep$slope_condition <- "random"
factor_sweep$sweep           <- "factor"

deposition <- expand.grid(rep = seq_len(N_REP), window = names(WINDOWS),
                          stringsAsFactors = FALSE)
deposition$N               <- N_FACTOR
deposition$lab_error       <- REF_ERROR
deposition$slope_condition <- "random"
deposition$growth_ratio    <- GROWTH_RATIO
deposition$sweep           <- "deposition"

core$growth_ratio         <- 1
factor_sweep$growth_ratio <- 1

keep   <- c("sweep", "rep", "N", "slope_condition", "lab_error", "window",
            "growth_ratio")
design <- rbind(core[, keep], factor_sweep[, keep], deposition[, keep])

design <- add_nuisance(design)

design$window_start <- sapply(WINDOWS[design$window], `[`, 1)
design$window_end   <- sapply(WINDOWS[design$window], `[`, 2)
design$pad          <- pad
design$grid_step    <- GRID_STEP

# Fixed reference constants for the calendar axis, passed to Stan as data so the
# normalisation is identical for every dataset and every model.
design$time_ref_min   <- design$window_start - pad
design$time_ref_range <- (design$window_end + pad) - design$time_ref_min

write.csv(design, OUTPUT_CSV, row.names = FALSE)

cat(sprintf("wrote %d datasets to %s\n", nrow(design), OUTPUT_CSV))
print(table(design$sweep, design$window))
print(table(design$sweep, design$growth_ratio))
