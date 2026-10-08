# Case 3, time as response: builds data/design.csv, one row per simulated dataset.
# Each dataset follows Cases 1 and 2 (time as response): values x drawn between
# 0 and 1000, true dates = intercept + slope * x + noise, then the dates are
# split into overlapping phases and each find gets the window of its phase
# (03_recovery_study.R does the simulating). The phases cover the range of the
# dataset's true dates, so every window contains its true date.
#
# The trend is set relative to the study period, as in Cases 1 and 2: across the
# values it spans 25-75% of the period, and the noise sd is 5-20% of the period.
# Dates are centred on the period.
#
# Sweeps, 100 datasets per setting, as in Case 3 time as predictor:
#   core    dataset size 50, 100, 200, 500 crossed with random / zero slope,
#           at overlap 0.5 and assign_p 0.2
#   factor  overlap 0, 0.25, 0.5 crossed with assign_p 0.5, 0.35, 0.2, at N = 200
#           (assign_p does nothing without overlap, so overlap 0 appears once)
library(here)

source(here("Simulations", "shared", "dating", "case3_overlapping_phases.R"))

period <- c(100, 900)
x.range <- c(0, 1000) # Range of the measured values
grid.step <- 5
n.rep <- 100

set.seed(2026)
core <- expand.grid(rep = 1:n.rep, N = c(50, 100, 200, 500), slope_condition = c("random", "zero"),
                    overlap = 0.5, assign_p = 0.2, stringsAsFactors = FALSE)
core$sweep <- "core"
factor.sweep <- expand.grid(rep = 1:n.rep, N = 200, slope_condition = "random",
                            overlap = c(0, 0.25, 0.5), assign_p = c(0.5, 0.35, 0.2), stringsAsFactors = FALSE)
factor.sweep <- factor.sweep[factor.sweep$overlap > 0 | factor.sweep$assign_p == 0.5, ]
factor.sweep$sweep <- "factor"
design <- rbind(core, factor.sweep)
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

# Phases per dataset: how many, and how even their lengths are
design$K <- sample(5:10, nrow(design), replace = TRUE)
design$alpha_conc <- 10^(-1 + 2 * rbeta(nrow(design), 2, 3.5))

design$x_min <- x.range[1]
design$x_max <- x.range[2]
design$grid_step <- grid.step

# Generator checks on large datasets at the widest trend: every window contains
# its true date, at every overlap
x <- runif(40000, x.range[1], x.range[2])
true.date <- rnorm(40000, centre - 0.75 * period.length / 2 + 0.75 * period.length * x / 1000, 0.2 * period.length)
for (ov in c(0, 0.25, 0.5))
{
	d <- simulate_overlap(40000, 0, 0, 1, K = 8, alpha_conc = 0.1, overlap = ov, assign_p = 0.2,
	                      true_date = true.date)
	stopifnot(all(d$True_date >= d$Start_date & d$True_date <= d$End_date))
}
cat("generator checks ok\n")

dir.create(here("Simulations", "time_as_response", "Case3_OverlappingPhases", "data"),
           showWarnings = FALSE, recursive = TRUE)
write.csv(design, here("Simulations", "time_as_response", "Case3_OverlappingPhases", "data", "design.csv"),
          row.names = FALSE)
cat("wrote", nrow(design), "datasets\n")
print(table(design$sweep))
