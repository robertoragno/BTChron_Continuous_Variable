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
// All dates are normalised to [0, 1] for numerical stability of HMC.
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
  array[N] real y;                     // response variable (e.g. LSI values)
  array[N] real start_date;            // earliest possible date for each observation
  array[N] real end_date;              // latest possible date for each observation
  int<lower=1> N_pred;                 // number of prediction points
  array[N_pred] real x_pred;           // time points to predict at
}

transformed data {
  // Compute the normalisation constants from the full date range.
  // Everything is mapped to [0, 1] so that HMC priors (normal(0,10) etc.)
  // are on a sensible scale. Raw years (e.g. -200 to 400 CE) would produce
  // tiny slopes and poor sampler geometry.
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
  // Latent date for each observation, as a raw uniform(0,1) variable.
  // The <lower=0, upper=1> constraint IS the prior — no explicit ~ statement
  // is needed. The transformed parameters block rescales each true_date_raw[n]
  // to its actual window [start_norm[n], end_norm[n]].
  // This reparameterisation (sample on [0,1], then stretch) is friendlier
  // to HMC than declaring a uniform with data-dependent bounds directly.
  array[N] real<lower=0, upper=1> true_date_raw;

  real alpha;                          // intercept at t_norm = 0 (i.e. at time_min)
  real beta1;                          // slope before the changepoint (normalised scale)
  real beta2;                          // slope after the changepoint (normalised scale)

  // Global changepoint, shared across all N observations.
  // The posterior for cp_norm is informed by every observation at once:
  // each observation is implicitly assigned to the first or second slope
  // depending on whether its latent date falls before or after cp_norm,
  // and cp_norm settles where the total score across all N is highest.
  real<lower=0, upper=1> cp_norm;

  real<lower=0> sigma;                 // observation noise
}

transformed parameters {
  // Rescale each raw [0,1] latent date to its actual [start_norm, end_norm] window.
  // true_date_norm[n] is what enters the likelihood.
  vector[N] true_date_norm;
  for (n in 1:N)
    true_date_norm[n] = start_norm[n]
                        + true_date_raw[n] * (end_norm[n] - start_norm[n]);
}

model {
  // Priors — all on the normalised [0,1] time scale.
  // normal(0, 10) is weakly informative given that t is in [0,1].
  alpha   ~ normal(0, 10);
  beta1   ~ normal(0, 10);
  beta2   ~ normal(0, 10);
  cp_norm ~ uniform(0, 1);   // no prior preference for where the kink falls
  sigma   ~ exponential(1);  // weakly informative, keeps sigma positive

  // Likelihood.
  // For each observation we compute mu using the broken-stick formula.
  // max(0, t - cp_norm) is the key: it is 0 when t < cp (first slope only)
  // and equals t - cp when t > cp (activating the slope change).
  for (n in 1:N) {
    real t  = true_date_norm[n];
    real mu = alpha + beta1 * t + (beta2 - beta1) * fmax(0.0, t - cp_norm);
    y[n] ~ normal(mu, sigma);
  }
}

generated quantities {
  // --- Pointwise log-likelihood for LOO-CV (loo package) ---
  vector[N] log_lik;
  for (n in 1:N) {
    real t  = true_date_norm[n];
    real mu = alpha + beta1 * t + (beta2 - beta1) * fmax(0.0, t - cp_norm);
    log_lik[n] = normal_lpdf(y[n] | mu, sigma);
  }

  // --- Posterior predictive quantities ---

  // mu_pred: the expected value of the trend line at each prediction point.
  // No noise is added — this is the fitted curve itself, used for plotting
  // the mean trend.
  array[N_pred] real mu_pred;

  // y_rep: a simulated new observation at each prediction point, drawn from
  // Normal(mu_pred, sigma). Unlike mu_pred, this includes observation noise,
  // so it represents what a new data point would look like. Use this for
  // posterior predictive intervals (the wider band around the trend line).
  array[N_pred] real y_rep;

  for (p in 1:N_pred) {
    real t     = x_pred_norm[p];
    mu_pred[p] = alpha + beta1 * t + (beta2 - beta1) * fmax(0.0, t - cp_norm);
    y_rep[p]   = normal_rng(mu_pred[p], sigma);
  }

  // --- Back-transform latent dates to original time scale ---
  array[N] real true_date_actual;
  for (n in 1:N)
    true_date_actual[n] = time_min + true_date_norm[n] * time_range;

  // --- Back-transform parameters to original time scale ---

  // Changepoint in original time units (e.g. years CE)
  real cp_actual = time_min + cp_norm * time_range;

  // Intercept at time_min — no rescaling needed because at t_norm = 0
  // (i.e. t = time_min) the normalised and original intercepts coincide.
  real baseline_original = alpha;

  // Slopes: beta1/beta2 are Δy per unit of normalised time.
  // Dividing by time_range converts to Δy per unit of original time.
  real slope1_original = beta1 / time_range;
  real slope2_original = beta2 / time_range;
}
