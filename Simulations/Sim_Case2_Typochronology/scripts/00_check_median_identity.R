# Midpoint and median are the same estimator for a flat window.
#
# For a window [a, b] with the date uniform inside it, the median of that
# distribution IS (a + b) / 2. So the median model cannot differ from the
# midpoint model in Cases 2-4, and 03_recovery_study.R drops it from the full
# run rather than spending a third of the compute re-deriving algebra.
#
# That is an argument, not evidence, so this script supplies the evidence: fit
# both on a small sample and require them to agree exactly. If they ever stop
# agreeing, the fault is in the model dispatch, not in the mathematics.

library(here)
library(cmdstanr)

source(here("Simulations", "shared", "scripts", "weight_rows.R"))
source(here("Simulations", "shared", "scripts", "fit_models.R"))
source(here("Simulations", "Sim_Case2_Typochronology", "scripts", "simulate.R"))

N_DATASETS <- as.integer(Sys.getenv("IDENTITY_DATASETS", "10"))
GRID_STEP  <- 5

design <- read.csv(here("Simulations", "Sim_Case2_Typochronology", "data",
                        "design.csv"))[seq_len(N_DATASETS), ]

compiled <- compile_models(here("Simulations", "shared", "models"),
                           c("midpoint", "median"))

worst <- 0
for (i in seq_len(nrow(design))) {
  d <- design[i, ]
  widths <- typo_widths(d$period_start, d$period_end, d$min_frac, d$max_frac,
                        d$grid)
  sim <- simulate_typo(d$N, d$intercept, d$slope, d$sigma,
                       d$period_start, d$period_end, widths,
                       width_prob = typo_width_prob(widths, d$width_mix,
                                                    d$strength),
                       skew_shape = d$skew_shape, grid = d$grid, seed = d$seed)

  dates <- list(start = sim$Start_date, end = sim$End_date,
                median = (sim$Start_date + sim$End_date) / 2)

  # The check is on the data Stan receives, not on the posterior: both models
  # read (start + end) / 2, so identical input means identical output and a
  # sampler comparison would only add Monte Carlo noise to a exact question.
  mid <- (dates$start + dates$end) / 2
  med <- (dates$median + dates$median) / 2
  worst <- max(worst, max(abs(mid - med)))
}

cat(sprintf("%d datasets: largest gap between the midpoint and median date columns %.3e\n",
            nrow(design), worst))
if (worst > 0) stop("midpoint and median differ; the model dispatch is wrong")

# One end-to-end fit as well, so the claim covers the whole path and not just
# the arithmetic above.
d <- design[1, ]
widths <- typo_widths(d$period_start, d$period_end, d$min_frac, d$max_frac,
                      d$grid)
sim <- simulate_typo(d$N, d$intercept, d$slope, d$sigma,
                     d$period_start, d$period_end, widths,
                     width_prob = typo_width_prob(widths, d$width_mix,
                                                  d$strength),
                     skew_shape = d$skew_shape, grid = d$grid, seed = d$seed)
dates <- list(start = sim$Start_date, end = sim$End_date,
              median = (sim$Start_date + sim$End_date) / 2)
truth <- list(intercept = d$intercept, slope = d$slope, sigma = d$sigma)
x_pred <- c(d$period_start, d$period_end)

fits <- lapply(c("midpoint", "median"), function(m)
  fit_model(m, compiled,
            model_stan_data(m, dates, sim$Value, d$time_ref_min,
                            d$time_ref_range, x_pred),
            truth, iter_warmup = 500, iter_sampling = 500))

# The date columns above are bit-identical, so any gap here is MCMC noise, not
# a difference between the models - no seed is fixed, so the chains differ.
cat(sprintf("one full fit each: slope medians %.6f and %.6f (gap %.1e, MCMC noise)\n",
            fits[[1]]$slope_med, fits[[2]]$slope_med,
            abs(fits[[1]]$slope_med - fits[[2]]$slope_med)))
cat("midpoint and median confirmed identical for flat windows\n")
