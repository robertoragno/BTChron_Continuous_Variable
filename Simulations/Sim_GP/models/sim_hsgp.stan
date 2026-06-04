// HSGP regression model with latent date inference.
//
// The GP trend f(t) = mu + PHI(t)' * beta approximates a squared-exponential
// GP via the Hilbert-space low-rank basis expansion (Solin & Sarkka 2020).
// The basis functions are fixed harmonic modes on [-L, L]; the spectral
// weights beta are scaled by the EQ spectral density evaluated at each
// eigenfrequency. This avoids building and inverting an N×N covariance matrix.
//
// Each observation has an uncertain date drawn uniformly within its known
// [start_date, end_date] window. Latent dates are inferred jointly with the
// GP hyperparameters, so posterior uncertainty captures both the noise in
// the signal and temporal imprecision in the dating.
//
// Time is normalised to [0, 1] internally for numerical stability. The
// generated quantities block back-transforms dates and the GP length-scale
// to the original time units.

functions {
  vector diagSPD_EQ(real alpha, real rho, real L, int M) {
    return alpha * sqrt(sqrt(2 * pi()) * rho)
           * exp(-0.25 * (rho * pi() / 2 / L)^2
                 * linspaced_vector(M, 1, M)^2);
  }

  matrix PHI(int N, int M, real L, vector x) {
    return sin(diag_post_multiply(
      rep_matrix(pi() / (2 * L) * (x + L), M),
      linspaced_vector(M, 1, M)
    )) / sqrt(L);
  }
}

data {
  int<lower=1> N;
  array[N] real y;
  array[N] real start_date;
  array[N] real end_date;

  int<lower=1> N_pred;
  array[N_pred] real x_pred;

  int<lower=1> M;         // number of HSGP basis functions
  real<lower=1> c;        // boundary factor; the HSGP domain extends to [-c, c]
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

  // c > 1 ensures the normalised data (in [0, 1]) sits strictly inside [-L, L]
  real L = c;

  // Prediction basis is constant across iterations, so build it here once
  matrix[N_pred, M] PHI_pred = PHI(N_pred, M, L, x_pred_norm);
}

parameters {
  // Each raw value is a uniform draw on [0, 1] that gets stretched into the
  // actual [start_norm, end_norm] window in transformed parameters.
  array[N] real<lower=0, upper=1> true_date_raw;

  real mu;                  // overall mean of the GP process
  real<lower=0> alpha;      // GP marginal SD (sets the scale of wiggles in y)
  real<lower=0> rho;        // GP length-scale on the normalised [0, 1] time axis
  vector[M] z;              // standard-normal basis weights (non-centred form)
  real<lower=0> sigma;      // observation noise
}

transformed parameters {
  vector[N] true_date_norm;
  for (n in 1:N)
    true_date_norm[n] = start_norm[n]
                        + true_date_raw[n] * (end_norm[n] - start_norm[n]);

  // Scale the raw basis weights by the spectral density of the EQ kernel
  vector[M] beta_gp = diagSPD_EQ(alpha, rho, L, M) .* z;
}

model {
  mu    ~ normal(0, 15);      // weakly informative over the observed value range
  alpha ~ normal(0, 10);      // half-normal; amplitude up to ~12 is plausible
  rho   ~ inv_gamma(5, 2.5);  // peaks near 0.5 on [0,1]; rules out near-zero and near-1
  z     ~ std_normal();
  sigma ~ exponential(1);

  // true_date_norm changes every iteration, so PHI_obs must be rebuilt here
  matrix[N, M] PHI_obs = PHI(N, M, L, true_date_norm);

  for (n in 1:N)
    y[n] ~ normal(mu + PHI_obs[n] * beta_gp, sigma);
}

generated quantities {
  vector[N] log_lik;
  {
    matrix[N, M] PHI_obs_gq = PHI(N, M, L, true_date_norm);
    for (n in 1:N)
      log_lik[n] = normal_lpdf(y[n] | mu + PHI_obs_gq[n] * beta_gp, sigma);
  }

  // Noise-free trend evaluated at the prediction grid
  array[N_pred] real mu_pred;
  // Posterior predictive draws that include observation noise
  array[N_pred] real y_rep;
  {
    vector[N_pred] f_grid = mu + PHI_pred * beta_gp;
    for (p in 1:N_pred) {
      mu_pred[p] = f_grid[p];
      y_rep[p]   = normal_rng(f_grid[p], sigma);
    }
  }

  // Latent dates back on the original time scale
  array[N] real true_date_actual;
  for (n in 1:N)
    true_date_actual[n] = time_min + true_date_norm[n] * time_range;

  // rho is sampled on [0,1]; multiply by time_range to recover original units
  real rho_actual = rho * time_range;
}
