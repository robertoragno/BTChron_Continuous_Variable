#' CASE 2: Archaeological finds with dating uncertainty coming from
#' the context of deposition (example of typochronology)
#'
#' We test the effects of finds from well-dated contexts (a short window)
#' and finds from mixed or coarsely dated contexts (a long window).
#' The true trend is linear, but the observed trend is biased by the dating
#' uncertainty. 
#' Memo for the readme/paper:
#' - Each find is assumed to be independent.
#' - Window lengths are shares of the study period, because the result
#' depends on how long a window is compared with the period, not on its
#' length in years (so the results are scale-invariant).

#' Function simulate_typo to simulate a single (Case 2) dataset.
#'
#' @param N                    number of samples within a dataset
#' @param intercept, slope, sigma; the simulated known trend: 
#'        Value = intercept + slope * date + noise
#' @param prop_coarse_samples  proportion of samples that are coarsely dated
#' @param coarse_frac          typical window of a coarsely dated sample, as a share
#'                             of the period. Each sample draws its own, between 0.8
#'                             and 1.2 times this.
#' @param fine_frac            window of a well-dated sample, as a share of the period
#' @param precision_trend      how much higher the chance of being coarsely dated
#'                             is at the start of the period than at the end. 0
#'                             means no change over time (uniform). This is an optional
#'                             parameter to test the effect of a non-uniform distribution of
#'                             coarse samples over time (for instance prehistory vs historical periods). 
#' @param growth_ratio         how many times denser finds are at the end of the
#'                             period than at the start. 1 is uniform (even spread).
#' @param period_start, period_end  the study period
#' @param seed                 random seed
simulate_typo <- function(N, intercept, slope, sigma,
                          prop_coarse_samples = 0.5, coarse_frac = 0.25,
                          fine_frac = 0.0625, precision_trend = 0,
                          growth_ratio = 1, period_start = 100,
                          period_end = 900, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  span <- period_end - period_start

  # When each find was deposited. Same draw as rdeposition() in Case 1.
  if (growth_ratio == 1) {
    true_date <- runif(N, period_start, period_end)
  } else {
    k <- log(growth_ratio) / span
    true_date <- period_start + log(1 + runif(N) * (growth_ratio - 1)) / k
  }

  # Which finds are coarsely dated
  position <- (true_date - period_start) / span   # 0 at the start, 1 at the end
  p_coarse <- prop_coarse_samples + precision_trend * (0.5 - position)
  p_coarse <- pmin(pmax(p_coarse, 0), 1)
  coarse   <- runif(N) < p_coarse

  width <- ifelse(coarse, runif(N, 0.8, 1.2) * coarse_frac, fine_frac) * span

  #' The window contains the true date, which is equally likely anywhere in it.
  #' Windows can run past the period: a date range does not know where the study
  #' stops. 
  #' Note to self: If I force windows within the study period they tend to be placed 
  #' at the center.
  start <- true_date - runif(N) * width

  data.frame(Start_date = start,
             End_date   = start + width,
             True_date  = true_date,
             Coarse     = coarse,
             Value      = intercept + slope * true_date + rnorm(N, 0, sigma))
}
