// =============================================================================
// Changepoint Regression with Latent Date Inference
// =============================================================================
// Broken-stick (segmented) parameterisation, continuous changepoint.
//
// The mean function is a single formula that behaves like two separate lines:
//
//   mu(t) = alpha + beta1*t + (beta2 - beta1) * max(0, t - cp)
//
// When t < cp: max(0, t - cp) = 0, so the last term vanishes entirely:
//   mu(t) = alpha + beta1*t                        [first slope only]
//
// When t > cp: max(0, t - cp) = t - cp, and beta1 cancels algebraically:
//   mu(t) = alpha + beta2*t + cp*(beta1 - beta2)   [second slope only]
//
// The line is continuous at cp — no jump, just a change of direction (a kink).
//
// All dates are mapped to [-1, 1] internally: centring the covariate keeps
// alpha and the slopes from being strongly correlated in the posterior.
// Slopes and intercepts are therefore on the normalised scale; the
// generated quantities block back-transforms them to the original time scale.
//
// Latent dates: each observation has an uncertain date within a known
// [start_date, end_date] window. These are inferred jointly with all
// other parameters.
//
// Changepoint identification: cp_norm is a single global parameter shared
// across all N observations. Its posterior is driven by all N residuals
// simultaneously — it concentrates wherever placing the kink minimises
// total residual variance across the whole dataset.
// =============================================================================

data {
  int<lower=1> N;                      // number of observations
  vector[N] y;                         // response variable (e.g. LSI values)
  vector[N] start_date;                // earliest possible date for each observation
  vector[N] end_date;                  // latest possible date for each observation
  int<lower=1> N_pred;                 // number of prediction points
  vector[N_pred] x_pred;               // time points to predict at
}

transformed data {
  // Compute the normalisation constants from the full date range.
  // Everything is mapped to [-1, 1] so that HMC priors (normal(0,10) etc.)
  // are on a sensible scale and alpha/beta1/beta2 are not strongly
  // correlated in the posterior. Raw years (e.g. -200 to 400 CE) would
  // produce tiny slopes and poor sampler geometry.
  real time_min   = min(start_date);
  real time_max   = max(end_date);
  real time_range = time_max - time_min;

  vector[N] start_norm = 2 * (start_date - time_min) / time_range - 1;
  vector[N] end_norm   = 2 * (end_date   - time_min) / time_range - 1;
  vector[N_pred] x_pred_norm = 2 * (x_pred - time_min) / time_range - 1;
}

parameters {
  // Latent date for each observation, as a raw uniform(0,1) variable.
  // The <lower=0, upper=1> constraint IS the prior — no explicit ~ statement
  // is needed. The transformed parameters block rescales each date_raw[n]
  // to its actual window [start_norm[n], end_norm[n]].
  // This reparameterisation (sample on [0,1], then stretch) is friendlier
  // to HMC than declaring a uniform with data-dependent bounds directly.
  vector<lower=0, upper=1>[N] date_raw;

  real alpha;                          // intercept at t_norm = -1 (i.e. at time_min)
  real beta1;                          // slope before the changepoint (normalised scale)
  real beta2;                          // slope after the changepoint (normalised scale)

  // Global changepoint, shared across all N observations.
  // The posterior for cp_norm is informed by every observation at once:
  // each observation is implicitly assigned to the first or second slope
  // depending on whether its latent date falls before or after cp_norm,
  // and cp_norm settles where the total score across all N is highest.
  // The <lower=-1, upper=1> constraint IS the prior (flat), matching the
  // normalised time range — no explicit ~ statement is needed.
  real<lower=-1, upper=1> cp_norm;

  real<lower=0> sigma;                 // observation noise
}

transformed parameters {
  // Rescale each raw [0,1] latent date to its actual [start_norm, end_norm] window.
  // date_norm is what enters the likelihood.
  vector[N] date_norm = start_norm + date_raw .* (end_norm - start_norm);
}

model {
  // Priors — all on the normalised [-1, 1] time scale.
  // normal(0, 10) is weakly informative given that t is in [-1, 1].
  alpha ~ normal(0, 10);
  beta1 ~ normal(0, 10);
  beta2 ~ normal(0, 10);
  sigma ~ exponential(1);  // weakly informative, keeps sigma positive

  // Likelihood, vectorised. fmax(0, t - cp_norm) is 0 when t < cp (first
  // slope only) and equals t - cp when t > cp (activating the slope change).
  y ~ normal(alpha + beta1 * date_norm
             + (beta2 - beta1) * fmax(0, date_norm - cp_norm), sigma);
}

generated quantities {
  // --- Pointwise log-likelihood for LOO-CV (loo package) ---
  vector[N] log_lik;
  {
    vector[N] mu = alpha + beta1 * date_norm
                  + (beta2 - beta1) * fmax(0, date_norm - cp_norm);
    for (n in 1:N)
      log_lik[n] = normal_lpdf(y[n] | mu[n], sigma);
  }

  // --- Posterior predictive quantities ---

  // mu_pred: the expected value of the trend line at each prediction point.
  // No noise is added — this is the fitted curve itself, used for plotting
  // the mean trend.
  vector[N_pred] mu_pred = alpha + beta1 * x_pred_norm
                           + (beta2 - beta1) * fmax(0, x_pred_norm - cp_norm);

  // y_rep: a simulated new observation at each prediction point, drawn from
  // Normal(mu_pred, sigma). Unlike mu_pred, this includes observation noise,
  // so it represents what a new data point would look like. Use this for
  // posterior predictive intervals (the wider band around the trend line).
  array[N_pred] real y_rep = normal_rng(mu_pred, sigma);

  // --- Back-transform latent dates to original time scale ---
  vector[N] date_actual = time_min + (date_norm + 1) / 2 * time_range;

  // --- Back-transform parameters to original time scale ---

  // Changepoint in original time units (e.g. years CE)
  real cp_actual = time_min + (cp_norm + 1) / 2 * time_range;

  // Intercept at time_min (t_norm = -1): alpha is the value at t_norm = 0,
  // so the value at t_norm = -1 on the first segment is alpha - beta1.
  real baseline_original = alpha - beta1;

  // Slopes: beta1/beta2 are Δy per unit of normalised time, and the
  // normalised time scale spans 2 units (-1 to 1) per time_range.
  // Dividing by time_range/2 (i.e. multiplying by 2) converts to Δy per
  // unit of original time.
  real slope1_original = 2 * beta1 / time_range;
  real slope2_original = 2 * beta2 / time_range;
}
