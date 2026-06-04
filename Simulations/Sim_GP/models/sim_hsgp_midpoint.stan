// HSGP regression using midpoint dates (no latent date inference).
//
// Dates are fixed at (start_date + end_date) / 2 before sampling begins.
// The GP mean function is identical to sim_gp.stan, so LOO scores from the
// two models are directly comparable. A better ELPD in the latent date
// model means within-window date uncertainty is genuinely informing the
// posterior; if the two are indistinguishable, fixing the midpoint is a
// reasonable simplification.
//
// Because midpoints are constants, the PHI basis matrices for both observations
// and prediction points can be precomputed once in transformed data rather than
// rebuilt on every HMC iteration.

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

  int<lower=1> M;
  real<lower=1> c;
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

  real L = c;

  // Both matrices are fixed constants, so they only need to be computed here
  matrix[N, M]      PHI_mid  = PHI(N,      M, L, mid_norm);
  matrix[N_pred, M] PHI_pred = PHI(N_pred, M, L, x_pred_norm);
}

parameters {
  real mu;
  real<lower=0> alpha;
  real<lower=0> rho;
  vector[M] z;
  real<lower=0> sigma;
}

transformed parameters {
  vector[M] beta_gp = diagSPD_EQ(alpha, rho, L, M) .* z;
}

model {
  mu    ~ normal(0, 15);
  alpha ~ normal(0, 10);
  rho   ~ inv_gamma(5, 2.5);
  z     ~ std_normal();
  sigma ~ exponential(1);

  for (n in 1:N)
    y[n] ~ normal(mu + PHI_mid[n] * beta_gp, sigma);
}

generated quantities {
  vector[N] log_lik;
  for (n in 1:N)
    log_lik[n] = normal_lpdf(y[n] | mu + PHI_mid[n] * beta_gp, sigma);

  // Noise-free trend on the prediction grid
  array[N_pred] real mu_pred;
  // Posterior predictive draws including observation noise
  array[N_pred] real y_rep;
  {
    vector[N_pred] f_grid = mu + PHI_pred * beta_gp;
    for (p in 1:N_pred) {
      mu_pred[p] = f_grid[p];
      y_rep[p]   = normal_rng(f_grid[p], sigma);
    }
  }

  real rho_actual = rho * time_range;
}
