// =============================================================================
// Simple Linear Regression with Latent Date Inference
// =============================================================================
//
//   trend(x) = alpha + beta * x
//   y_n ~ Normal(trend(true_date_n), sigma)
//
// Latent dates are inferred uniformly within each [start, end] window.
// All dates are normalised to [0, 1] internally for numerical stability.
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

  vector[N] start_norm;
  vector[N] end_norm;
  vector[N_pred] x_pred_norm;

  for (n in 1:N) {
    start_norm[n] = (start_date[n] - time_min) / time_range;
    end_norm[n]   = (end_date[n]   - time_min) / time_range;
  }
  for (p in 1:N_pred)
    x_pred_norm[p] = (x_pred[p] - time_min) / time_range;
}

parameters {
  array[N] real<lower=0, upper=1> true_date_raw;

  real alpha;
  real beta;
  real<lower=0> sigma;
}

transformed parameters {
  vector[N] true_date_norm;
  for (n in 1:N)
    true_date_norm[n] = start_norm[n]
                        + true_date_raw[n] * (end_norm[n] - start_norm[n]);
}

model {
  alpha ~ normal(0, 10);
  beta  ~ normal(0, 10);
  sigma ~ exponential(1);

  for (n in 1:N) {
    real trend_n = alpha + beta * true_date_norm[n];
    y[n] ~ normal(trend_n, sigma);
  }
}

generated quantities {
  // Pointwise log-likelihood (LOO-CV)
  vector[N] log_lik;
  for (n in 1:N) {
    real trend_n = alpha + beta * true_date_norm[n];
    log_lik[n] = normal_lpdf(y[n] | trend_n, sigma);
  }

  // Expected trend on prediction grid (noise-free)
  array[N_pred] real mu_pred;
  // Posterior predictive draws (includes observation noise)
  array[N_pred] real y_rep;
  for (p in 1:N_pred) {
    mu_pred[p] = alpha + beta * x_pred_norm[p];
    y_rep[p]   = normal_rng(mu_pred[p], sigma);
  }

  // Back-transform latent dates to original scale
  array[N] real true_date_actual;
  for (n in 1:N)
    true_date_actual[n] = time_min + true_date_norm[n] * time_range;

  // Back-transform trend parameters to original scale
  real slope_original    = beta / time_range;
  real baseline_original = alpha - beta * time_min / time_range;
}
