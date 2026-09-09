// =============================================================================
// Linear Regression with Marginalised Latent Dates, Study-Period Prior
// =============================================================================
//
// PROTOTYPE. marginal_date.stan with one term added: the dates are given a
// shared distribution whose location and spread are estimated from the dates
// themselves.
//
// Why. rcarbon calibrates against a uniform prior over the whole calibration
// range, so each weight row carries mass outside the period the finds actually
// came from, and nothing ties the dates to each other. The fitted spread of
// dates then exceeds the real one, and since a slope is roughly
// cov(x, y) / var(x), an inflated var(x) flattens it. That is the suspected
// cause of the control-window attenuation, and this file is the test of it.
//
// A weight row is the calibrated likelihood normalised under a uniform prior,
// so multiplying it by a prior on the year and renormalising gives the
// posterior under that prior instead:
//
//   p(y_n, C14_n) = SUM over t of  w_n(t) * Normal(t | mu, tau) * Normal(y_n | trend(t), sigma)
//
// The normalising constant depends on mu and tau and is deliberately NOT
// divided out - it carries the information about how well the shared period
// accounts for the dates, which is what makes the period estimable.
//
// Normal is the wrong shape: the generator draws dates uniformly over the
// window. It is used because attenuation is driven by var(x) rather than by the
// shape of the date distribution, so recovering the variance should be enough
// to say whether the mechanism is the right one. A trapezoid - ramp, plateau,
// ramp - is the shape to reach for if this confirms the diagnosis, being both
// closer to how deposition behaves and smooth at the boundaries in a way a
// uniform is not.
//
// Not for model comparison as it stands: log_lik below is left as the
// no-prior version, so LOO or WAIC against marginal_date.stan would compare
// unlike quantities.
// =============================================================================
//
//   trend(t) = alpha + beta * t
//   y_n ~ Normal(trend(theta_n), sigma),  theta_n unknown
//
// Each observation carries a probability for every year it could belong to,
// instead of a single [start, end] window, and those years are averaged over
// rather than sampled:
//
//   p(y_n | alpha, beta, sigma) = SUM over years t of prob[n, t] * Normal(y_n | trend(t), sigma)
//
// evaluated as log_sum_exp(log prob[n, .] + log-density). No date parameter is
// introduced, so the posterior stays at three parameters whatever N is, and the
// same file serves every dating case - only the probabilities change:
//
//   Case 1  calibrated radiocarbon posterior on the annual grid (rcarbon calMatrix)
//   Case 2  uniform over the sample's typochronological window
//   Case 3  uniform over the assigned phase's window
//   Case 4  uniform over the union of the merged phases' windows
//
// The probabilities sit on a SHARED annual grid, so they are supplied ragged:
// one flat vector, plus each row's first grid index and its length. Rows are
// trimmed to their support (a 100-yr window out of a 1500-yr axis is 93% zeros)
// and renormalised in R before being passed in.
//
// Calendar years are mapped to [-1, 1] using reference constants passed as data,
// not derived from the realised dates, so the prior on beta is the same prior on
// the calendar-scale slope in every dataset and every model.
// =============================================================================

data {
  int<lower=1> N;
  vector[N] y;

  // Shared candidate-year grid. grid_year[i] is the calendar year of grid
  // position i; this is what turns a row index into a date.
  int<lower=1> n_years;
  vector[n_years] grid_year;

  // Ragged probability rows, already normalised so each row sums to 1.
  //
  // Stan has no ragged array type and every observation covers a different
  // number of years, so the rows are glued end to end into one flat vector and
  // two index arrays say where each row lives inside it:
  //
  //   obs 1: 3 years   obs 2: 4 years   obs 3: 2 years
  //   packed = [ . . . | . . . . | . . ]
  //   position   1 2 3   4 5 6 7   8 9
  //   row_first_year and row_offset point at 1, 4, 8; row_n_years is 3, 4, 2
  //
  // Nothing here is ever modified - this is the data block. The model reads one
  // observation's slice with segment() and calls it candidate_log_prob.
  int<lower=1> n_weights;                  // total stored probabilities
  vector[n_weights] log_year_prob_packed;  // log probabilities, rows concatenated
  array[N] int<lower=1, upper=n_years> row_first_year;  // grid position each row starts at
  array[N] int<lower=1> row_n_years;                    // years covered by each row

  // Calendar axis reference, shared across datasets and models.
  real time_ref_min;
  real<lower=0> time_ref_range;

  int<lower=1> N_pred;
  vector[N_pred] x_pred;
}

transformed data {
  vector[n_years] grid_year_norm =
    2 * (grid_year - time_ref_min) / time_ref_range - 1;
  vector[N_pred] x_pred_norm = 2 * (x_pred - time_ref_min) / time_ref_range - 1;

  real neg_log_sqrt_2pi = -0.9189385332046727;  // -0.5 * log(2 * pi())

  // Where each row begins inside the flat log_year_prob_packed vector.
  array[N] int row_offset;
  {
    int pos = 1;
    for (n in 1:N) {
      row_offset[n] = pos;
      pos += row_n_years[n];
    }
    if (pos - 1 != n_weights)
      reject("row_n_years sums to ", pos - 1,
             " but log_year_prob_packed has length ", n_weights);
  }

  // Rows must be proper distributions over the grid. Note what this does and
  // does not buy: scaling a row by a constant adds a constant to the log
  // target, so it does NOT bias alpha, beta or sigma. It does shift that row's
  // log_lik by log(c), which would corrupt any ELPD comparison against another
  // model, and an off row sum is the symptom of a real construction bug - a
  // mis-trimmed or truncated support - which is what this check is really for.
  for (n in 1:N) {
    real row_sum =
      exp(log_sum_exp(segment(log_year_prob_packed, row_offset[n],
                              row_n_years[n])));
    if (abs(row_sum - 1) > 1e-6)
      reject("probability row ", n, " sums to ", row_sum, ", expected 1");
    if (row_first_year[n] + row_n_years[n] - 1 > n_years)
      reject("probability row ", n, " runs past the end of the grid");
  }

  // Half a grid step, in normalised units, for jittering sampled dates.
  real half_step_norm = n_years > 1
    ? (grid_year_norm[2] - grid_year_norm[1]) / 2
    : 1.0 / time_ref_range;
}

parameters {
  real alpha;
  real beta;
  real<lower=0> sigma;

  // The shared date distribution, on the normalised [-1, 1] calendar axis.
  // The true 400-yr window inside a 1576-yr grid spans about 0.51 normalised
  // units, so a uniform over it has sd near 0.15; the priors below are wide
  // enough not to impose that and weak enough to keep tau off zero.
  real mu_norm;
  real<lower=0> tau_norm;
}

model {
  alpha ~ normal(0, 10);
  // beta is on the [-1, 1] axis: beta = slope_original * time_ref_range / 2.
  // Over a ~1600-yr span the design's steepest slopes (+/- 0.03 y/yr) reach
  // |beta| ~ 24, and normal(0, 10) shrank those toward zero, manufacturing slope
  // attenuation that has nothing to do with dating. normal(0, 40) covers the
  // design range for every case.
  beta  ~ normal(0, 40);
  sigma ~ exponential(1);
  mu_norm  ~ normal(0, 1);
  tau_norm ~ normal(0, 0.5);

  // For each observation, ask how well the current line explains it at every
  // year it could belong to, then average those answers weighted by how likely
  // each year is. Vectorised over a row's years rather than looped: the scalar
  // version costs one autodiff node per year per leapfrog step.
  for (n in 1:N) {
    // the years this observation could belong to, and their probabilities
    vector[row_n_years[n]] candidate_year =
      segment(grid_year_norm, row_first_year[n], row_n_years[n]);
    vector[row_n_years[n]] candidate_log_prob =
      segment(log_year_prob_packed, row_offset[n], row_n_years[n]);

    // how far y sits from the line at each of those years, in units of sigma
    vector[row_n_years[n]] resid_scaled =
      (y[n] - alpha - beta * candidate_year) / sigma;

    // log(probability of the year) + log(normal density of that miss),
    // summed over years on the probability scale
    // The shared period, evaluated at each candidate year. Written out rather
    // than normal_lpdf(), which sums a vector instead of returning one value
    // per element.
    vector[row_n_years[n]] log_period_prior = neg_log_sqrt_2pi - log(tau_norm)
      - 0.5 * square((candidate_year - mu_norm) / tau_norm);

    target += log_sum_exp(candidate_log_prob + log_period_prior
                          + neg_log_sqrt_2pi - log(sigma)
                          - 0.5 * square(resid_scaled));
  }
}

generated quantities {
  // Marginal, not conditional on a latent date, unlike latent_date.stan
  vector[N] log_lik;
  // One posterior draw of each date, stacking across iterations into the same
  // density the latent model's date_actual gives
  vector[N] date_actual;

  for (n in 1:N) {
    vector[row_n_years[n]] candidate_year =
      segment(grid_year_norm, row_first_year[n], row_n_years[n]);
    vector[row_n_years[n]] candidate_log_prob =
      segment(log_year_prob_packed, row_offset[n], row_n_years[n]);
    vector[row_n_years[n]] resid_scaled =
      (y[n] - alpha - beta * candidate_year) / sigma;

    // How well each candidate year accounts for this observation, on the log
    // scale: its own probability times the density of the miss it implies.
    vector[row_n_years[n]] log_year_fit = candidate_log_prob + neg_log_sqrt_2pi
                                          - log(sigma) - 0.5 * square(resid_scaled);

    log_lik[n] = log_sum_exp(log_year_fit);

    // Pick a candidate year weighted by how well it fits, then jitter within
    // the cell so stacked draws read as a density rather than as spikes.
    //
    // A calibrated row trimmed to one contiguous span can hold exact zeros in
    // its interior - the near-empty gaps between the humps of a plateau date -
    // and log(0) is -inf, which log_sum_exp accepts but categorical_logit_rng
    // rejects. Floor the weights 700 log-units below the row maximum: that is
    // a probability around 1e-304, so the year is still never drawn.
    vector[row_n_years[n]] log_year_fit_safe =
      fmax(log_year_fit, max(log_year_fit) - 700);

    real year_norm = candidate_year[categorical_logit_rng(log_year_fit_safe)];
    real year_jittered =
      uniform_rng(year_norm - half_step_norm, year_norm + half_step_norm);
    date_actual[n] = time_ref_min + (year_jittered + 1) / 2 * time_ref_range;
  }

  vector[N_pred] mu_pred = alpha + beta * x_pred_norm;
  array[N_pred] real y_rep = normal_rng(mu_pred, sigma);

  real slope_original    = 2 * beta / time_ref_range;
  real baseline_original = alpha - beta - 2 * beta * time_ref_min / time_ref_range;

  // The estimated period, back on the calendar scale, so it can be read against
  // the window the design actually used.
  real period_centre_original = time_ref_min + (mu_norm + 1) / 2 * time_ref_range;
  real period_sd_original     = tau_norm * time_ref_range / 2;
}
