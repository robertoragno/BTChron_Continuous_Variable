// HSGP regression model with latent date inference.
//
// The GP trend f(t) = mu + PHI(t)' * beta approximates a squared-exponential
// GP via the Hilbert-space low-rank basis expansion (https://link.springer.com/article/10.1007/s11222-019-09886-w).
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

  // Prediction basis, built once here since the grid never changes across
  // iterations. PHI_pred is [N_pred x M]: row p holds the M basis functions
  // evaluated at grid time x_pred_norm[p]. Multiplied by the weights beta_gp
  // later, it turns the basis into one trend value per grid point.
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
  alpha ~ normal(0, 10);      // half-normal; amplitude of the wiggles
  // rho is the GP length-scale on the normalised [0,1] axis. inv_gamma(5, 0.9)
  // puts the 90% mass on roughly 0.10-0.46, i.e. ~80-370 years on an 800-year
  // span: the signal is expected to change over a century or a few, not to be
  // either flat across the whole period or spikier than the dating can resolve.
  rho   ~ inv_gamma(5, 0.9);
  z     ~ std_normal();
  // half-normal, weakly informative over the generating range (sigma ~ 0.5-4).
  // exponential(1) has mean 1, too tight here: it fights the data and biases
  // sigma low. Widening it also eases the latent-date sampling geometry.
  sigma ~ normal(0, 5);

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

  // --- Posterior predictive quantities on the prediction grid ---
  // The grid is the N_pred times in x_pred where the fitted curve is drawn for
  // plotting, independent of the observed dates.

  // mu_pred: noise-free trend at each grid point (the fitted curve itself).
  array[N_pred] real mu_pred;
  // y_rep: a simulated new observation at each grid point (mu_pred plus noise).
  // Use mu_pred for the mean trend line, y_rep for the wider predictive band.
  array[N_pred] real y_rep;
  {
    // The GP trend at every grid point in one matrix-vector product:
    // PHI_pred is [N_pred x M] and beta_gp is length M, so PHI_pred * beta_gp is
    // a length-N_pred vector; adding the mean mu gives the trend at each point.
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
