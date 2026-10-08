# Case 4, time as response: builds data/design.csv, one row per simulated dataset.
# Each dataset follows Cases 1 and 2 (time as response): values x drawn between
# 0 and 1000, true dates = intercept + slope * x + noise, then the dates are
# split into fine phases and each find gets a window of 1 to merge_max adjacent
# phases (03_recovery_study.R does the simulating). The phases cover the range
# of the dataset's true dates, so every window contains its true date.
#
# The trend is set relative to the study period, as in Cases 1 and 2: across the
# values it spans 25-75% of the period, and the noise sd is 5-20% of the period.
# Dates are centred on the period.
#
# Sweeps, 100 datasets per setting, as in Case 4 time as predictor:
#   core   dataset size 50, 100, 200, 500 crossed with random / zero slope,
#          windows of up to 3 phases
#   merge  windows of up to 2, 3 or 4 phases, at N = 200
library(here)

source(here("Simulations", "shared", "dating", "case4_merged_phases.R"))

period <- c(100, 900)
x.range <- c(0, 1000) # Range of the measured values
grid.step <- 5
n.rep <- 100

set.seed(2026)
core <- expand.grid(rep = 1:n.rep, N = c(50, 100, 200, 500), slope_condition = c("random", "zero"),
                    merge_max = 3, stringsAsFactors = FALSE)
core$sweep <- "core"
merge.sweep <- expand.grid(rep = 1:n.rep, N = 200, slope_condition = "random",
                           merge_max = c(2, 3, 4), stringsAsFactors = FALSE)
merge.sweep$sweep <- "merge"
design <- rbind(core, merge.sweep)
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

# Fine phases per dataset: how many, and how even their lengths are
design$K <- sample(5:10, nrow(design), replace = TRUE)
design$alpha_conc <- 10^(-1 + 2 * rbeta(nrow(design), 2, 3.5))

design$x_min <- x.range[1]
design$x_max <- x.range[2]
design$grid_step <- grid.step

# Generator checks on large datasets at the widest trend: every window contains
# its true date, and each dataset mixes single phases with merged ones
x <- runif(40000, x.range[1], x.range[2])
true.date <- rnorm(40000, centre - 0.75 * period.length / 2 + 0.75 * period.length * x / 1000, 0.2 * period.length)
for (mm in c(2, 3, 4))
{
	d <- simulate_merged(40000, 0, 0, 1, K = 8, alpha_conc = 1, merge_max = mm, true_date = true.date)
	stopifnot(all(d$True_date >= d$Start_date & d$True_date <= d$End_date),
	          min(d$Phases) == 1, max(d$Phases) == mm)
}
cat("generator checks ok\n")

dir.create(here("Simulations", "time_as_response", "Case4_MergedPhases", "data"),
           showWarnings = FALSE, recursive = TRUE)
write.csv(design, here("Simulations", "time_as_response", "Case4_MergedPhases", "data", "design.csv"),
          row.names = FALSE)
cat("wrote", nrow(design), "datasets\n")
print(table(design$sweep))
