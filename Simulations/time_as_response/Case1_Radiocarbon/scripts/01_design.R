# Case 1, time as response: builds data/design.csv, one row per simulated dataset.
# Each dataset follows measurement_error_example.R: values x drawn between 0 and
# 1000, true dates = intercept + slope * x + noise, then each date is turned into
# a radiocarbon age and calibrated (03_recovery_study.R does the simulating).
#
# The trend is set relative to the window, not in years: across the values it
# spans 25-75% of the window, and the noise sd is 5-20% of the window. Dates are
# centred on the window. Crema's example (270 years of trend, sd 70, 400-year
# window) sits inside these ranges.
#
# Sweeps, 100 datasets per setting, the rest at the reference (N = 200, lab
# error 30, Hallstatt plateau, random slope):
#   core    dataset size 50, 100, 200, 500 crossed with random / zero slope, plateau
#   factor  lab error 15, 30, 50 crossed with window
library(here)

windows <- list(plateau = c(-800, -400), steep = c(-1600, -1200))
x.range <- c(0, 1000) # Range of the measured values
pad <- 950 # Grid padding: widest calibrated date (588 yr) plus 4 noise sd at most
grid.step <- 5
n.rep <- 100

set.seed(2026)
core <- expand.grid(rep = 1:n.rep, N = c(50, 100, 200, 500), slope_condition = c("random", "zero"),
                    lab_error = 30, window = "plateau", stringsAsFactors = FALSE)
core$sweep <- "core"
factor.sweep <- expand.grid(rep = 1:n.rep, N = 200, slope_condition = "random",
                            lab_error = c(15, 30, 50), window = names(windows), stringsAsFactors = FALSE)
factor.sweep$sweep <- "factor"
design <- rbind(core, factor.sweep)
design$dataset_id <- 1:nrow(design)
design$seed <- design$dataset_id

design$window_start <- sapply(windows[design$window], function(w) w[1])
design$window_end <- sapply(windows[design$window], function(w) w[2])
window.length <- design$window_end - design$window_start
centre <- (design$window_start + design$window_end) / 2

# True trend per dataset: years per unit of x, a random sign, noise in years
span <- runif(nrow(design), 0.25, 0.75) * window.length
direction <- sample(c(-1, 1), nrow(design), replace = TRUE)
design$slope <- ifelse(design$slope_condition == "zero", 0, direction * span / diff(x.range))
design$sigma <- runif(nrow(design), 0.05, 0.2) * window.length
design$intercept <- centre - design$slope * mean(x.range) # date at x = 0

design$x_min <- x.range[1]
design$x_max <- x.range[2]
design$grid_min <- design$window_start - pad
design$grid_max <- design$window_end + pad
design$grid_step <- grid.step

dir.create(here("Simulations", "time_as_response", "Case1_Radiocarbon", "data"),
           showWarnings = FALSE, recursive = TRUE)
write.csv(design, here("Simulations", "time_as_response", "Case1_Radiocarbon", "data", "design.csv"),
          row.names = FALSE)
cat("wrote", nrow(design), "datasets\n")
print(table(design$sweep, design$window))
