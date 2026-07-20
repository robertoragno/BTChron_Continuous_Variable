// =============================================================================
// Changepoint Regression with Marginalized Latent Dates
// =============================================================================
// Same broken-stick trend and the same "each bone's true date is unknown,
// uniform within its phase window" assumption as sim_changepoint.stan. The
// difference is HOW that assumption enters the likelihood.
//
// sim_changepoint.stan gives each observation its own free parameter
// (true_date_raw[n]) and lets HMC sample it. For a bone whose window
// straddles the changepoint, that parameter effectively has to pick a side
// (before/after the kink), and different chains can settle on different
// sides -- a multimodality that shows up as poor mixing / high Rhat / dropped
// fits, worse for wide windows and steep slope changes.
//
// Here, instead of sampling a date, each observation's likelihood is the
// date integrated out analytically:
//
//   p(y_n | params) = INTEGRAL over [start_n, end_n] of
//                      (1/width_n) * Normal(y_n | mu(t), sigma) dt
//
// approximated by a fixed grid of M points per observation (midpoint-rule
// quadrature). Because the prior density (1/width) and the quadrature step
// (width/M) cancel, this is just the AVERAGE likelihood across M evenly
// spaced candidate dates -- a log_sum_exp minus log(M), no extra parameters,
// same trick used for interval-censored / rounded data in Stan.
//
// No true_date parameter means no per-observation mode for chains to
// disagree about: cp_norm, beta1, beta2 see the exact (quadrature-accurate)
// marginal likelihood on every iteration.
//
// Per-observation date information is not lost -- generated quantities
// reconstructs one posterior draw of each bone's date per MCMC iteration, by
// sampling from the same grid weights used in the likelihood. Stacked across
// draws these behave like the true_date_actual draws in the latent model, so
// the same downstream plotting code (density of estimated dates against the
// phase window and the true date) applies unchanged.
// =============================================================================

data {
  int<lower=1> N;
  array[N] real y;
  array[N] real start_date;
  array[N] real end_date;
  int<lower=1> N_pred;
  array[N_pred] real x_pred;
}

transformed data {
  real time_min   = min(start_date);
  real time_max   = max(end_date);
  real time_range = time_max - time_min;

  int M = 40;  // quadrature points per observation window
  real neg_log_sqrt_2pi = -0.9189385332046727;  // -0.5 * log(2 * pi())

  vector[N_pred] x_pred_norm;
  for (p in 1:N_pred)
    x_pred_norm[p] = (x_pred[p] - time_min) / time_range;

  // Grid of M candidate dates per observation, normalised to [0, 1], plus
  // the half-width of one grid cell (used later to jitter a sampled grid
  // index back into a continuous-looking date).
  array[N] vector[M] time_grid;
  vector[N] half_cell_norm;
  for (n in 1:N) {
    real s = (start_date[n] - time_min) / time_range;
    real e = (end_date[n]   - time_min) / time_range;
    real cell = (e - s) / M;
    half_cell_norm[n] = cell / 2;
    for (k in 1:M)
      time_grid[n, k] = s + (k - 0.5) * cell;
  }
}

parameters {
  real alpha;
  real beta1;
  real beta2;
  real<lower=0, upper=1> cp_norm;
  real<lower=0> sigma;
}

model {
  alpha   ~ normal(0, 10);
  beta1   ~ normal(0, 10);
  beta2   ~ normal(0, 10);
  cp_norm ~ uniform(0, 1);
  sigma   ~ exponential(1);

  // Each observation's date, marginalized over its window: average the
  // likelihood across M grid points instead of sampling one date. Vectorized
  // over the grid (one fmax + one vector of normal log-densities per
  // observation) rather than a scalar loop over M -- the scalar version costs
  // roughly M Stan-function-call autodiff nodes per observation per leapfrog
  // step; this collapses that to a handful of vector ops.
  for (n in 1:N) {
    vector[M] excess = fmax(time_grid[n] - cp_norm, 0.0);
    vector[M] mu     = alpha + beta1 * time_grid[n] + (beta2 - beta1) * excess;
    vector[M] z      = (y[n] - mu) / sigma;
    vector[M] lp      = neg_log_sqrt_2pi - log(sigma) - 0.5 * square(z);
    target += log_sum_exp(lp) - log(M);
  }
}

generated quantities {
  vector[N] log_lik;
  array[N] real true_date_actual;

  for (n in 1:N) {
    vector[M] excess = fmax(time_grid[n] - cp_norm, 0.0);
    vector[M] mu     = alpha + beta1 * time_grid[n] + (beta2 - beta1) * excess;
    vector[M] z      = (y[n] - mu) / sigma;
    vector[M] lp      = neg_log_sqrt_2pi - log(sigma) - 0.5 * square(z);
    log_lik[n] = log_sum_exp(lp) - log(M);

    // One posterior draw of this bone's date: pick a grid cell weighted by
    // how well it fits (softmax of lp), then jitter uniformly within that
    // cell so stacked draws look like a continuous density, not spikes.
    int k_draw = categorical_logit_rng(lp);
    real t_centre = time_grid[n, k_draw];
    real t_jittered = uniform_rng(t_centre - half_cell_norm[n],
                                  t_centre + half_cell_norm[n]);
    true_date_actual[n] = time_min + t_jittered * time_range;
  }

  array[N_pred] real mu_pred;
  array[N_pred] real y_rep;
  for (p in 1:N_pred) {
    real t     = x_pred_norm[p];
    mu_pred[p] = alpha + beta1 * t + (beta2 - beta1) * fmax(0.0, t - cp_norm);
    y_rep[p]   = normal_rng(mu_pred[p], sigma);
  }

  real cp_actual          = time_min + cp_norm * time_range;
  real baseline_original  = alpha;
  real slope1_original    = beta1 / time_range;
  real slope2_original    = beta2 / time_range;
}
