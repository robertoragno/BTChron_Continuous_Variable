# Builds the Case 1 design grid and writes it to data/design.csv.
# 03_recovery_study.R reads that file and re-derives nothing.
#
# One sweep at a time, the rest at the reference (N = 200, lab error 30,
# Hallstatt plateau, even deposition, random slope):
#
#   core        dataset size crossed with random / zero slope, on the plateau
#   factor      lab error crossed with calendar window
#   deposition  finds 4 times denser at the end of the window, both windows

library(here)

source(here("Simulations", "shared", "scripts", "design.R"))
source(here("Simulations", "shared", "scripts", "weight_rows.R"))
source(here("Simulations", "Sim_Case1_Radiocarbon", "scripts", "simulate.R"))

# ---- settings ----
# Two windows of the same length. On the plateau IntCal20 is nearly flat, so
# calibrated dates are wide and multi-peaked. On the steep section they are
# narrow and single-peaked.
WINDOWS      <- list(plateau = c(-800, -400),
                     steep   = c(-1600, -1200))
LAB_ERRORS   <- c(15, 30, 50)
REF_ERROR    <- 30
REF_WINDOW   <- "plateau"
N_LEVELS     <- c(50, 100, 200, 500)
N_FACTOR     <- 200
GROWTH_RATIO <- 4
GRID_STEP    <- 5    # spacing of candidate years for the full-distribution model
N_REP        <- 100  # datasets per setting
OUTPUT_CSV   <- Sys.getenv("DESIGN_CASE1_OUT",
                           here("Simulations", "Sim_Case1_Radiocarbon", "data",
                                "design.csv"))
# ------------------

set.seed(2026)

# How far the calendar grid extends past each window: the widest calibrated date
# over both windows and all lab errors (588 yr). This is the first rcarbon call,
# and loading rcarbon draws random numbers, so moving it before set.seed()
# changes every dataset's intercept, slope and sigma.
PAD <- max(sapply(WINDOWS, function(w)
  sapply(LAB_ERRORS, function(e) c14_grid_pad(w, e))))

reference <- data.frame(N = N_FACTOR, slope_condition = "random",
                        lab_error = REF_ERROR, window = REF_WINDOW,
                        growth_ratio = 1)

# Repeat the reference, then overwrite the column(s) a sweep varies
sweep_block <- function(name, varied) {
  grid <- expand.grid(c(list(rep = seq_len(N_REP)), varied),
                      stringsAsFactors = FALSE)
  for (col in setdiff(names(reference), names(varied))) grid[[col]] <- reference[[col]]
  grid$sweep <- name
  grid[, c("sweep", "rep", names(reference))]
}

design <- rbind(
  sweep_block("core",       list(N = N_LEVELS, slope_condition = c("random", "zero"))),
  sweep_block("factor",     list(lab_error = LAB_ERRORS, window = names(WINDOWS))),
  sweep_block("deposition", list(window = names(WINDOWS), growth_ratio = GROWTH_RATIO))
)

design <- add_true_trend(design)

design$window_start <- sapply(WINDOWS[design$window], `[`, 1)
design$window_end   <- sapply(WINDOWS[design$window], `[`, 2)
design$pad          <- PAD
design$grid_step    <- GRID_STEP

# Calendar axis passed to Stan, the same for every dataset and model
design$time_ref_min   <- design$window_start - PAD
design$time_ref_range <- (design$window_end + PAD) - design$time_ref_min

dir.create(dirname(OUTPUT_CSV), showWarnings = FALSE, recursive = TRUE)
write.csv(design, OUTPUT_CSV, row.names = FALSE)

cat(sprintf("wrote %d datasets to %s\n", nrow(design), OUTPUT_CSV))
print(table(design$sweep, design$window))
