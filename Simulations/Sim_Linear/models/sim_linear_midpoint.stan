// =============================================================================
// Linear Regression using Midpoint Dates (no latent date inference)
// =============================================================================
// Baseline comparison: dates are fixed at (start + end) / 2.
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

  vector[N] mid_norm;
  vector[N_pred] x_pred_norm;

  for (n in 1:N)
    mid_norm[n] = ((start_date[n] + end_date[n]) / 2.0 - time_min) / time_range;

  for (p in 1:N_pred)
    x_pred_norm[p] = (x_pred[p] - time_min) / time_range;
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

  for (n in 1:N)
    y[n] ~ normal(alpha + beta * mid_norm[n], sigma);
}

generated quantities {
  vector[N] log_lik;
  for (n in 1:N)
    log_lik[n] = normal_lpdf(y[n] | alpha + beta * mid_norm[n], sigma);

  array[N_pred] real trend_pred;
  array[N_pred] real obs_sim;
  for (p in 1:N_pred) {
    trend_pred[p] = alpha + beta * x_pred_norm[p];
    obs_sim[p]    = normal_rng(trend_pred[p], sigma);
  }

  real slope_original    = beta / time_range;
  real baseline_original = alpha - beta * time_min / time_range;
}
