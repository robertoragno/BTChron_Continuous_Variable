// =============================================================================
// Changepoint Regression using Midpoint Dates (no latent date inference)
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

  for (n in 1:N) {
    real t  = mid_norm[n];
    real mu = alpha + beta1 * t + (beta2 - beta1) * fmax(0.0, t - cp_norm);
    y[n] ~ normal(mu, sigma);
  }
}

generated quantities {
  vector[N] log_lik;
  for (n in 1:N) {
    real t  = mid_norm[n];
    real mu = alpha + beta1 * t + (beta2 - beta1) * fmax(0.0, t - cp_norm);
    log_lik[n] = normal_lpdf(y[n] | mu, sigma);
  }

  // mu_pred: the model's expected value (the trend line) at each prediction point,
  // with no observation noise added. We use this to plot the fitted curve.
  array[N_pred] real mu_pred;
  // y_rep: a simulated observation at each prediction point, drawn from
  // Normal(mu_pred, sigma). Unlike mu_pred, this includes sigma noise, so it
  // represents what a new data point would look like. We use this for predictive intervals.
  array[N_pred] real y_rep;
  for (p in 1:N_pred) {
    real t     = x_pred_norm[p];
    mu_pred[p] = alpha + beta1 * t + (beta2 - beta1) * fmax(0.0, t - cp_norm);
    y_rep[p]   = normal_rng(mu_pred[p], sigma);
  }

  real cp_actual         = time_min + cp_norm * time_range;
  real baseline_original = alpha;
  real slope1_original   = beta1 / time_range;
  real slope2_original   = beta2 / time_range;
}
