// =============================================================================
// Simple Linear Regression with Latent Date Inference
// =============================================================================
//
//   trend(x) = alpha + beta * x
//   y_n ~ Normal(trend(date_n), sigma)
//
// Latent dates are inferred uniformly within each [start, end] window.
// All dates are mapped to [-1, 1] internally: centring the covariate keeps
// alpha and beta from being strongly correlated in the posterior.
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

  vector[N] start_norm = 2 * (start_date - time_min) / time_range - 1;
  vector[N] end_norm   = 2 * (end_date   - time_min) / time_range - 1;
  vector[N_pred] x_pred_norm = 2 * (x_pred - time_min) / time_range - 1;
}

parameters {
  vector<lower=0, upper=1>[N] date_raw;

  real alpha;
  real beta;
  real<lower=0> sigma;
}

transformed parameters {
  vector[N] date_norm = start_norm + date_raw .* (end_norm - start_norm);
}

model {
  alpha ~ normal(0, 10);
  beta  ~ normal(0, 10);
  sigma ~ exponential(1);

  y ~ normal(alpha + beta * date_norm, sigma);
}

generated quantities {
  // Pointwise log-likelihood (LOO-CV)
  vector[N] log_lik;
  {
    vector[N] mu = alpha + beta * date_norm;
    for (n in 1:N)
      log_lik[n] = normal_lpdf(y[n] | mu[n], sigma);
  }

  // Expected trend on prediction grid (noise-free)
  vector[N_pred] mu_pred = alpha + beta * x_pred_norm;
  // Posterior predictive draws (includes observation noise)
  array[N_pred] real y_rep = normal_rng(mu_pred, sigma);

  // Back-transform latent dates to original scale
  vector[N] date_actual = time_min + (date_norm + 1) / 2 * time_range;

  // Back-transform trend parameters to original scale
  real slope_original    = 2 * beta / time_range;
  real baseline_original = alpha - beta - 2 * beta * time_min / time_range;
}
