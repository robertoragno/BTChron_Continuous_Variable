// =============================================================================
// Changepoint Regression with Latent Date Inference
// =============================================================================
// Broken-stick (segmented) parameterization, continuous changepoint:
//
//   mu(t) = alpha + beta1*t + (beta2 - beta1) * max(0, t - cp)
//
// Equivalent to the piecewise form used to generate the data:
//   f(t) = baseline + slope1*(t - t_min)            for t <= cp
//   f(t) = baseline + slope1*(cp - t_min) + slope2*(t - cp)  for t > cp
//
// Latent dates are inferred uniformly within each [start, end] window.
// All dates are normalised to [0, 1] for numerical stability.
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
  real beta1;
  real beta2;
  real<lower=0, upper=1> cp_norm;
  real<lower=0> sigma;
}

transformed parameters {
  vector[N] true_date_norm;
  for (n in 1:N)
    true_date_norm[n] = start_norm[n]
                        + true_date_raw[n] * (end_norm[n] - start_norm[n]);
}

model {
  alpha   ~ normal(0, 10);
  beta1   ~ normal(0, 10);
  beta2   ~ normal(0, 10);
  cp_norm ~ uniform(0, 1);
  sigma   ~ exponential(1);

  for (n in 1:N) {
    real t  = true_date_norm[n];
    real mu = alpha + beta1 * t + (beta2 - beta1) * fmax(0.0, t - cp_norm);
    y[n] ~ normal(mu, sigma);
  }
}

generated quantities {
  // Pointwise log-likelihood (LOO-CV)
  vector[N] log_lik;
  for (n in 1:N) {
    real t  = true_date_norm[n];
    real mu = alpha + beta1 * t + (beta2 - beta1) * fmax(0.0, t - cp_norm);
    log_lik[n] = normal_lpdf(y[n] | mu, sigma);
  }

  // mu_pred: the model's expected value (the trend line) at each prediction point,
  // with no observation noise added. We use this to plot the fitted curve.
  array[N_pred] real mu_pred;
  // y_rep: a simulated observation at each prediction point, drawn from
  // Normal(mu_pred, sigma). Unlike mu_pred, this includes sigma noise, so it
  // represents what a new data point would look like. We would use this for predictive intervals.
  array[N_pred] real y_rep;
  for (p in 1:N_pred) {
    real t   = x_pred_norm[p];
    mu_pred[p] = alpha + beta1 * t + (beta2 - beta1) * fmax(0.0, t - cp_norm);
    y_rep[p]   = normal_rng(mu_pred[p], sigma);
  }

  // Back-transform latent dates to original scale
  array[N] real true_date_actual;
  for (n in 1:N)
    true_date_actual[n] = time_min + true_date_norm[n] * time_range;

  // Back-transform parameters to original scale
  real cp_actual         = time_min + cp_norm * time_range;
  real baseline_original = alpha;                  // f(t_min): alpha at t_norm = 0
  real slope1_original   = beta1 / time_range;
  real slope2_original   = beta2 / time_range;
}
