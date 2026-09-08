# Builds the Case 2 design grid and writes it to data/design_case2.csv.
#
# Written before any fitting; every downstream script reads this file rather than
# rebuilding the grid from a seed, so a partial re-run cannot quietly draw a
# different design. Edit the settings below and re-run to change the study.
#
# Three blocks:
#   core       sample size crossed with the slope condition, at even dating
#              resolution and no skew. The N sweep reads precision against
#              sample size; the zero-slope cells give the false-positive rate.
#   resolution fine vs coarse dating, at fixed N. "even" is not repeated here,
#              the matching core cell is its reference.
#   skew       deposition leaning to the late edge of the window (case 2b).
#   overhang   how wide windows are relative to the study period. Windows may
#              run past the period, and near the boundary the midpoint stops
#              being the conditional mean of the true date, which flattens the
#              slope. That is a real property of typological dating, not a
#              simulator defect, so it is swept and measured rather than tuned
#              away - shrinking windows until the midpoint looks unbiased would
#              be rigging the comparison.

library(here)

source(here("Simulations", "Sim_Case2_Typochronology", "scripts", "simulate.R"))
source(here("Simulations", "shared", "scripts", "design.R"))

# ---- settings to play with ----
PERIOD_START <- 100     # matches Sim_Linear's TMIN/TMAX so Case 3 stays comparable
PERIOD_END   <- 900
MIN_FRAC     <- 0.05    # narrowest window as a share of the period
MAX_FRAC     <- 0.25    # widest; past ~1/3 the midpoint slope starts flattening
# Levels for the overhang sweep. MAX_FRAC stays the reference: it is a realistic
# upper end for typological dating, and measured attenuation there is ~0.974.
OVERHANG_FRACS <- c(0.10, 0.25, 0.40)
GRID         <- 25      # window boundaries are multiples of this
STRENGTH     <- 4       # how lopsided fine/coarse are (1 = flat)
SKEW_LEVELS  <- c(2, 4) # Beta shapes for 2b, matching Sim_Linear's mild/strong
N_LEVELS     <- c(50, 100, 200, 500)
N_FACTOR     <- 200     # sample size the factor sweeps run at
N_REP        <- 100     # datasets per cell
# The skew cells are the underpowered ones: at 100 datasets the Monte Carlo
# standard error on the slope bias (~6e-05) is the same size as the effect being
# looked for, so nothing there can be resolved. Given more replicates.
N_REP_SKEW   <- 300
OUTPUT_CSV   <- Sys.getenv("DESIGN_CASE2_OUT",
                           here("Simulations", "Sim_Case2_Typochronology", "data",
                                "design.csv"))
# -------------------------------

set.seed(2026)

widths <- typo_widths(PERIOD_START, PERIOD_END, MIN_FRAC, MAX_FRAC, GRID)
cat("width set:", paste(widths, collapse = ", "), "yr  (",
    paste0(round(100 * range(widths) / (PERIOD_END - PERIOD_START)), "%",
           collapse = " to "), "of the period )\n")

core <- expand.grid(rep = seq_len(N_REP), N = N_LEVELS,
                    slope_condition = c("random", "zero"),
                    stringsAsFactors = FALSE)
core$width_mix  <- "even"
core$skew_shape <- 1
core$sweep      <- "core"

resolution <- expand.grid(rep = seq_len(N_REP), width_mix = c("fine", "coarse"),
                          stringsAsFactors = FALSE)
resolution$N               <- N_FACTOR
resolution$slope_condition <- "random"
resolution$skew_shape      <- 1
resolution$sweep           <- "resolution"

skew <- expand.grid(rep = seq_len(N_REP_SKEW), skew_shape = SKEW_LEVELS,
                    stringsAsFactors = FALSE)
skew$N               <- N_FACTOR
skew$slope_condition <- "random"
skew$width_mix       <- "even"
skew$sweep           <- "skew"

overhang <- expand.grid(rep = seq_len(N_REP), max_frac = OVERHANG_FRACS,
                        stringsAsFactors = FALSE)
overhang$N               <- N_FACTOR
overhang$slope_condition <- "random"
overhang$width_mix       <- "even"
overhang$skew_shape      <- 1
overhang$sweep           <- "overhang"

core$max_frac       <- MAX_FRAC
resolution$max_frac <- MAX_FRAC
skew$max_frac       <- MAX_FRAC

keep   <- c("sweep", "rep", "N", "slope_condition", "width_mix", "skew_shape",
            "max_frac")
design <- rbind(core[, keep], resolution[, keep], skew[, keep],
                overhang[, keep])

design <- add_nuisance(design)

# Recorded so every downstream script uses the same generator settings, and so
# the study can be reconstructed from this file alone.
design$period_start <- PERIOD_START
design$period_end   <- PERIOD_END
design$min_frac     <- MIN_FRAC
design$grid         <- GRID
design$strength     <- STRENGTH

# Fixed reference constants for the calendar axis, passed to Stan as data so the
# normalisation is identical for every dataset and model. Derived from the design
# rather than from each dataset's realised dates: with a per-dataset range,
# normal(0, s) on beta is a different prior on the calendar-scale slope in every
# dataset. Windows can overhang the period by up to one width either side.
#
# Padded by the widest window in the WHOLE design, not per cell - the overhang
# sweep varies max_frac, and a per-cell pad would give each level its own
# time_ref_range and so its own prior on the calendar-scale slope, confounding
# the sweep with the centring.
widest <- max(sapply(c(MAX_FRAC, OVERHANG_FRACS), function(f)
  max(typo_widths(PERIOD_START, PERIOD_END, MIN_FRAC, f, GRID))))
design$time_ref_min   <- PERIOD_START - widest
design$time_ref_range <- (PERIOD_END + widest) - design$time_ref_min

write.csv(design, OUTPUT_CSV, row.names = FALSE)

cat(sprintf("wrote %d datasets to %s\n", nrow(design), OUTPUT_CSV))
print(table(design$sweep, design$width_mix))
cat("\nwidest window across the whole design:", widest, "yr\n")
print(table(design$sweep, design$max_frac))
