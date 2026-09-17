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
// The period has no fixed shape: it is a fixed basis of evenly spaced bumps
// whose weights are estimated. Six free-moving bumps were tried first and were
// unusable - 2-4 hours a fit - because their positions and widths are only
// weakly identified. Fixing the basis keeps the flexibility and drops the cost.
//
// Two fixed shapes were tried first and both failed by getting the spread
// wrong, which is what drives slope attenuation. A normal collapsed to sd 93
// against a true 115 on the plateau. A trapezoid reached 98, still short, and
// only 76 against 105 once deposition was skewed - a symmetric flat top cannot
// represent exponential deposition. Heaton (arXiv:2109.15024) makes the general
// argument and fits a Dirichlet process mixture; a finite mixture is the
// standard truncation of that, and marginalising over the year grid here
// removes the multimodality that forced slice sampling there.
//
// Normalised over the grid rather than in closed form, which keeps the prior
// proper whatever the components do. This is the prior's own constant and is
// unrelated to the per-row constant above, which stays undivided.
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
// All finds share one grid of candidate years, but each find only covers some
// of them, so rows have different lengths. They are passed as one long vector
// with rows placed end to end, plus where each row starts on the grid and how
// long it is. Rows are trimmed to the years they cover (a 100-yr window out of
// a 1500-yr axis is 93% zeros) and renormalised in R before being passed in.
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

  // Probability rows, already normalised so each row sums to 1.
  //
  // Every observation covers a different number of years and Stan has no array
  // for rows of different lengths, so the rows are placed end to end in one
  // vector and two index arrays say where each row sits inside it:
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

  // Fixed basis of overlapping bumps, evenly spaced over the grid. Only their
  // weights are estimated, so the basis is built once here rather than every
  // iteration, and the period is one matrix-vector product in the model block.
  //
  // Letting the bumps move as well cost 2-4 hours a fit against 10 minutes for
  // a two-parameter shape: free centres and widths are barely identified, and
  // the sampler spends its time exploring which bump is which rather than what
  // the period looks like. Fixing them makes the density linear in the weights,
  // which is the whole saving.
  // K sets the resolution: the bumps are spaced 2/(K-1) apart and are 0.8 of
  // that wide, so K = 12 gives bumps of 115 yr on this axis - as wide as the
  // period being measured, which put the recovered spread at 173 against a true
  // 115. K = 40 gives 32 yr, fine enough to resolve a 400-yr window.
  // ponytail: fixed K and spacing, raise K if a period needs finer structure
  int K = 40;
  matrix[n_years, K] basis;
  {
    real spacing = 2.0 / (K - 1);
    for (k in 1:K) {
      vector[n_years] col =
        exp(-0.5 * square((grid_year_norm - (-1 + (k - 1) * spacing))
                          / (spacing * 0.8)));
      basis[, k] = col / sum(col);
    }
  }

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

  // How much of the period sits under each basis bump, on the log scale. A
  // random walk along the bumps rather than a simplex prior: see the model
  // block. softmax() makes it a proper density whatever the values do.
  vector[K] log_weight;
  real<lower=0> weight_sd;
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
  // Neighbouring bumps carry similar weight, by an amount the data sets. This
  // is what a period is: a contiguous stretch of time, not a scatter of
  // isolated years.
  //
  // Two dirichlet priors were tried first and neither works. Flat, it behaves
  // like K pseudo-observations spread over the whole calendar axis and the
  // recovered spread came out at 212 against a true 115. Sparse (alpha = 1/K)
  // recovers the spread but costs 120 minutes a fit against 3: below 1 the
  // density spikes wherever a weight approaches zero, and the sampler has to
  // crawl through that. Smoothness buys the same concentration with none of the
  // geometry, because adjacent-and-similar is a weaker demand than
  // few-and-isolated.
  log_weight[1] ~ normal(0, 2);
  log_weight[2:K] ~ normal(log_weight[1:(K - 1)], weight_sd);
  weight_sd ~ normal(0, 1);

  // Floored against underflow: a bump's tail can round to zero several bumps
  // away, and log(0) would take the whole target with it.
  vector[n_years] log_period = log(basis * softmax(log_weight) + 1e-12);

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
    vector[row_n_years[n]] log_period_prior =
      segment(log_period, row_first_year[n], row_n_years[n]);

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

  // The estimated period, back on the calendar scale. Reported as the fitted
  // density's own mean and sd rather than as the bounds: the ramps sit outside
  // the flat top, so the bounds overstate the width that var(x) actually sees.
  real period_centre_original;
  real period_sd_original;
  {
    vector[n_years] p  = basis * softmax(log_weight);
    vector[n_years] yr = time_ref_min + (grid_year_norm + 1) / 2 * time_ref_range;
    period_centre_original = dot_product(p, yr);
    period_sd_original = sqrt(dot_product(p, square(yr - period_centre_original)));
  }
}
