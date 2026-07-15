// =============================================================================
// Linear Regression using Midpoint Dates (no latent date inference)
// =============================================================================
// Baseline comparison: dates are fixed at (start + end) / 2.
// All dates are mapped to [-1, 1] internally, as in linear_latent.stan.
// =============================================================================

data {
  int<lower=1> N;
  vector[N] y;
  vector[N] start_date;
  vector[N] end_date;

  int<lower=1> N_pred;
  vector[N_pred] x_pred;
}

transformed data {
  real time_min   = min(start_date);
  real time_max   = max(end_date);
  real time_range = time_max - time_min;

  vector[N] mid_norm =
    2 * ((start_date + end_date) / 2 - time_min) / time_range - 1;
  vector[N_pred] x_pred_norm = 2 * (x_pred - time_min) / time_range - 1;
}

parameters {
  real alpha;
  real beta;
  real<lower=0> sigma;
}

model {
  alpha ~ normal(0, 10);
  beta  ~ normal(0, 10);
  sigma ~ exponential(1);

  y ~ normal(alpha + beta * mid_norm, sigma);
}

generated quantities {
  vector[N] log_lik;
  {
    vector[N] mu = alpha + beta * mid_norm;
    for (n in 1:N)
      log_lik[n] = normal_lpdf(y[n] | mu[n], sigma);
  }

  // Expected trend on prediction grid (noise-free)
  vector[N_pred] mu_pred = alpha + beta * x_pred_norm;
  // Posterior predictive draws (includes observation noise)
  array[N_pred] real y_rep = normal_rng(mu_pred, sigma);

  real slope_original    = 2 * beta / time_range;
  real baseline_original = alpha - beta - 2 * beta * time_min / time_range;
}
