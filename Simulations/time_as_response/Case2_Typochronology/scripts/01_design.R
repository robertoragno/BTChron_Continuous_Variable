# Case 2, time as response: builds data/design.csv, one row per simulated dataset.
# Each dataset follows Case 1 (time as response): values x drawn between 0 and
# 1000, true dates = intercept + slope * x + noise, then each find gets a date
# range that contains its true date (03_recovery_study.R does the simulating).
#
# The trend is set relative to the study period, as in Case 1: across the values
# it spans 25-75% of the period, and the noise sd is 5-20% of the period. Dates
# are centred on the period.
#
# Sweeps, 100 datasets per setting, the rest at the reference (N = 200, half the
# finds coarsely dated, coarse windows 25% of the period, no precision trend,
# random slope):
#   core       dataset size 50, 100, 200, 500 crossed with random / zero slope
#   prop       proportion of coarsely dated finds 0, 0.25, 0.75, 1
#   width      coarse window 10% or 45% of the period
#   precision  coarsely dated finds more common early than late
# The deposition sweep of time as predictor is left out: here the dates come
# from the trend, not from a deposition rate.
library(here)

source(here("Simulations", "shared", "dating", "case2_typochronology.R"))

period <- c(100, 900)
x.range <- c(0, 1000) # Range of the measured values
fine.frac <- 0.0625 # Window of a well-dated find: 50 yr of an 800-yr period
grid.step <- 5
n.rep <- 100

set.seed(2026)
core <- expand.grid(rep = 1:n.rep, N = c(50, 100, 200, 500), slope_condition = c("random", "zero"),
                    prop_coarse_samples = 0.5, coarse_frac = 0.25, precision_trend = 0,
                    stringsAsFactors = FALSE)
core$sweep <- "core"
prop <- expand.grid(rep = 1:n.rep, N = 200, slope_condition = "random",
                    prop_coarse_samples = c(0, 0.25, 0.75, 1), coarse_frac = 0.25, precision_trend = 0,
                    stringsAsFactors = FALSE)
prop$sweep <- "prop"
width <- expand.grid(rep = 1:n.rep, N = 200, slope_condition = "random",
                     prop_coarse_samples = 0.5, coarse_frac = c(0.10, 0.45), precision_trend = 0,
                     stringsAsFactors = FALSE)
width$sweep <- "width"
precision <- expand.grid(rep = 1:n.rep, N = 200, slope_condition = "random",
                         prop_coarse_samples = 0.5, coarse_frac = 0.25, precision_trend = 0.6,
                         stringsAsFactors = FALSE)
precision$sweep <- "precision"
design <- rbind(core, prop, width, precision)
design$dataset_id <- 1:nrow(design)
design$seed <- design$dataset_id

design$period_start <- period[1]
design$period_end <- period[2]
period.length <- diff(period)
centre <- mean(period)

# True trend per dataset: years per unit of x, a random sign, noise in years
span <- runif(nrow(design), 0.25, 0.75) * period.length
direction <- sample(c(-1, 1), nrow(design), replace = TRUE)
design$slope <- ifelse(design$slope_condition == "zero", 0, direction * span / diff(x.range))
design$sigma <- runif(nrow(design), 0.05, 0.2) * period.length
design$intercept <- centre - design$slope * mean(x.range) # date at x = 0

design$x_min <- x.range[1]
design$x_max <- x.range[2]
design$fine_frac <- fine.frac
design$grid_step <- grid.step

# Generator checks on a large dataset at the widest setting
x <- runif(40000, x.range[1], x.range[2])
true.date <- rnorm(40000, centre - 0.75 * period.length / 2 + 0.75 * period.length * x / 1000, 0.2 * period.length)
d <- simulate_typo(40000, 0, 0, 1, prop_coarse_samples = 0.5, coarse_frac = 0.45, fine_frac = fine.frac,
                   period_start = period[1], period_end = period[2], true_date = true.date)
stopifnot(all(d$True_date >= d$Start_date & d$True_date <= d$End_date),
          abs(mean(d$Coarse) - 0.5) < 0.01)
d <- simulate_typo(40000, 0, 0, 1, precision_trend = 0.6, fine_frac = fine.frac,
                   period_start = period[1], period_end = period[2], true_date = true.date)
early <- d$True_date < period[1] + 200
late <- d$True_date > period[2] - 200
stopifnot(mean(d$Coarse[early]) > mean(d$Coarse[late]) + 0.3)
cat("generator checks ok\n")

dir.create(here("Simulations", "time_as_response", "Case2_Typochronology", "data"),
           showWarnings = FALSE, recursive = TRUE)
write.csv(design, here("Simulations", "time_as_response", "Case2_Typochronology", "data", "design.csv"),
          row.names = FALSE)
cat("wrote", nrow(design), "datasets\n")
print(table(design$sweep))
