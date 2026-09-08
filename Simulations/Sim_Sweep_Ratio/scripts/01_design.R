# Builds the sweep design and writes it to data/design.csv.
#
# Two factors, crossed:
#
#   dating coarseness  a typical window width as a share of the studied period.
#                      Widths are held fixed in years and the period varies, so
#                      this moves without changing dating practice.
#   signal strength    the total change in y across the whole period, against a
#                      residual scatter of 0.5-4.
#
# The second factor is here because the first run held the signal at one level
# and produced a much smaller midpoint sigma error than Case 2's overhang sweep
# at comparable coarseness. Sigma inflation from unmodelled date error should
# scale with the slope, so the suspicion is that the two studies differ in
# signal rather than in coarseness. Crossing them settles it.

library(here)

source(here("Simulations", "Sim_Case2_Typochronology", "scripts", "simulate.R"))

# ---- settings to play with ----
# Window widths in YEARS, not as a share of the period. A realistic
# typochronological range, held constant while the period around it changes.
WIDTHS       <- seq(50, 200, by = 25)
GRID         <- 25
PERIOD_START <- 100
# Studied periods, from a single phase to a continental compilation. Chosen to
# roughly double each time, so the ratio axis is evenly spaced on a log scale.
PERIOD_SPANS <- c(200, 400, 800, 1600, 3200)
N            <- 200
# Total change in y across the whole period, drawn uniform on +/- the level.
# Fixing the rise rather than the slope keeps the signal the same size as the
# period grows; see README.md. The levels span faint to obvious against a
# scatter of 0.5-4, and 24 is what Case 2 draws (slope 0.03 over 800 yr).
RISE_LEVELS  <- c(4, 12, 24)
SIGMA_RANGE  <- c(0.5, 4)
N_REP        <- as.integer(Sys.getenv("SWEEP_REPS", "100"))
OUTPUT_CSV   <- Sys.getenv("SWEEP_DESIGN_OUT",
                           here("Simulations", "Sim_Sweep_Ratio", "data",
                                "design.csv"))
# -------------------------------

set.seed(2026)

design <- expand.grid(rep = seq_len(N_REP), period_span = PERIOD_SPANS,
                      rise_level = RISE_LEVELS, stringsAsFactors = FALSE)
design$dataset_id <- seq_len(nrow(design))
design$seed       <- design$dataset_id
design$N          <- N

design$period_start <- PERIOD_START
design$period_end   <- PERIOD_START + design$period_span

design$intercept <- runif(nrow(design), 2, 15)
design$sigma     <- runif(nrow(design), SIGMA_RANGE[1], SIGMA_RANGE[2])
design$rise      <- runif(nrow(design), -design$rise_level, design$rise_level)
design$slope     <- design$rise / design$period_span

# The axis grows with the period, padded by the widest window either side
# because windows may overhang the period. Because the slope shrinks as the
# period grows while the axis grows with it, the NORMALISED slope stays roughly
# constant, so normal(0, 10) on beta means about the same thing in every cell.
# The oracle model is what actually demonstrates this; see README.md.
widest <- max(WIDTHS)
design$time_ref_min   <- design$period_start - widest
design$time_ref_range <- (design$period_end + widest) - design$time_ref_min

# The quantity the whole study is about, recorded so no downstream script has to
# recompute it and risk defining it differently.
design$mean_width <- mean(WIDTHS)
design$coarseness <- design$mean_width / design$period_span

write.csv(design, OUTPUT_CSV, row.names = FALSE)

cat(sprintf("wrote %d datasets to %s\n", nrow(design), OUTPUT_CSV))
cat("window widths (fixed):", paste(WIDTHS, collapse = ", "), "yr, mean",
    mean(WIDTHS), "\n\n")
print(unique(design[, c("period_span", "coarseness")]), row.names = FALSE)
cat("\nsignal: total rise across the period, against scatter",
    paste(SIGMA_RANGE, collapse = "-"), "\n")
print(table(design$rise_level, design$period_span))
