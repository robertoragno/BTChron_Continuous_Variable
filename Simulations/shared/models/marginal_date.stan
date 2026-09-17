// Linear regression with the dates integrated out
//
//   trend(t) = alpha + beta * t
//   y_n ~ Normal(trend(theta_n), sigma),  theta_n unknown
//
// Each observation has a probability for every year it could date to. The model
// averages over those years instead of sampling a date:
//
//   p(y_n | alpha, beta, sigma) = SUM over years t of
//   prob[n, t] * Normal(y_n | trend(t), sigma)
//
// computed as log_sum_exp(log prob[n, .] + log density). There is no date
// parameter, so the model has three parameters whatever N is. The same file is
// used by every case, only the probabilities change:
//
//   Case 1  calibrated radiocarbon date on the grid (rcarbon calMatrix)
//   Case 2  flat over the sample's typochronological window
//   Case 3  flat over the assigned phase's window
//   Case 4  flat over the merged phases' window
//
// All finds share one grid of candidate years, but each find only covers some
// of them, so rows have different lengths (ragged).
// They are passed as one long vector with rows placed end to end, plus where
// each row starts on the grid and how long it is.
// Rows are trimmed to the years they cover (a 100-yr window out of
// a 1500-yr axis is 93% zeros) and renormalised in R before being passed in.
//
// Calendar years are rescaled to [-1, 1] with time_ref_min (start of the axis)
// and time_ref_range (its length), set in each case's 01_design.R. They are
// the same for every dataset in a case, not taken from each dataset's dates,
// so the prior on beta means the same slope in years in every dataset and model.

data {
  int<lower=1> N;
  vector[N] y;

  // Grid of candidate years. grid_year[i] is the calendar year at position i.
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
  //   row_offset is 1, 4, 8; row_n_years is 3, 4, 2
  //
  // row_first_year is where each row starts on the grid of years, not in the
  // packed (ragged) vector. The model reads one observation's row with segment().
  int<lower=1> n_weights;                  // total stored probabilities
  vector[n_weights] log_year_prob_packed;  // log probabilities, rows end to end
  array[N] int<lower=1, upper=n_years> row_first_year;  // grid position each row starts at
  array[N] int<lower=1> row_n_years;                    // years covered by each row

  // Calendar axis, the same for every dataset and model
  real time_ref_min;
  real<lower=0> time_ref_range;

  int<lower=1> N_pred;
  vector[N_pred] x_pred;
}

transformed data {
  vector[n_years] grid_year_norm =
    2 * (grid_year - time_ref_min) / time_ref_range - 1;
  vector[N_pred] x_pred_norm = 2 * (x_pred - time_ref_min) / time_ref_range - 1;

  // Constant term of the normal log density, which the model writes out by hand
  // because it needs one density per candidate year. It does not change the
  // fit, but keeps log_lik comparable with midpoint.stan (e.g. for LOO).
  real neg_log_sqrt_2pi = -0.5 * log(2 * pi());

  // Where each row starts inside log_year_prob_packed
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

  // Each row must sum to 1 and stay on the grid. A row that does not sum to 1
  // would not change alpha, beta or sigma, but it shifts that row's log_lik and
  // usually means the row was cut or trimmed wrongly in R.
  for (n in 1:N) {
    real row_sum =
      exp(log_sum_exp(segment(log_year_prob_packed, row_offset[n],
                              row_n_years[n])));
    if (abs(row_sum - 1) > 1e-6)
      reject("probability row ", n, " sums to ", row_sum, ", expected 1");
    if (row_first_year[n] + row_n_years[n] - 1 > n_years)
      reject("probability row ", n, " runs past the end of the grid");
  }

  // Half a grid step on the [-1, 1] scale, to spread sampled dates within a cell
  real half_step_norm = n_years > 1
    ? (grid_year_norm[2] - grid_year_norm[1]) / 2
    : 1.0 / time_ref_range;
}

parameters {
  real alpha;
  real beta;
  real<lower=0> sigma;
}

model {
  alpha ~ normal(0, 10);
  // Time runs from -1 to 1, but beta can take any value: it is the change in y
  // per unit of rescaled time, beta = slope_original * time_ref_range / 2.
  // Over a ~1600-yr axis the steepest slopes in the design (+/- 0.03 per yr)
  // give |beta| ~ 24. normal(0, 10) pulled those towards zero and flattened the
  // slope for reasons unrelated to dating. normal(0, 40) covers every case.
  beta  ~ normal(0, 40);
  sigma ~ exponential(1);

  // For each observation, how well the line fits it at each year it could date
  // to, averaged with each year's probability as weight. Written over the whole
  // row at once rather than year by year, which is much faster in Stan.
  for (n in 1:N) {
    // candidate years for this observation and their log probabilities
    vector[row_n_years[n]] candidate_year =
      segment(grid_year_norm, row_first_year[n], row_n_years[n]);
    vector[row_n_years[n]] candidate_log_prob =
      segment(log_year_prob_packed, row_offset[n], row_n_years[n]);

    // distance of y from the line at each candidate year, in units of sigma
    vector[row_n_years[n]] resid_scaled =
      (y[n] - alpha - beta * candidate_year) / sigma;

    // log(probability of the year) + log(normal density of the distance),
    // summed over years on the probability scale
    target += log_sum_exp(candidate_log_prob + neg_log_sqrt_2pi - log(sigma)
                          - 0.5 * square(resid_scaled));
  }
}

generated quantities {
  // log likelihood with the date averaged out (latent_date.stan keeps the date)
  vector[N] log_lik;
  // One date per observation per draw. Over all draws these give the same
  // distribution as date_actual in latent_date.stan.
  vector[N] date_actual;

  for (n in 1:N) {
    vector[row_n_years[n]] candidate_year =
      segment(grid_year_norm, row_first_year[n], row_n_years[n]);
    vector[row_n_years[n]] candidate_log_prob =
      segment(log_year_prob_packed, row_offset[n], row_n_years[n]);
    vector[row_n_years[n]] resid_scaled =
      (y[n] - alpha - beta * candidate_year) / sigma;

    // How well each candidate year fits this observation, on the log scale:
    // the year's probability times the normal density of the distance.
    vector[row_n_years[n]] log_year_fit = candidate_log_prob + neg_log_sqrt_2pi
                                          - log(sigma) - 0.5 * square(resid_scaled);

    log_lik[n] = log_sum_exp(log_year_fit);

    // Draw a year weighted by how well it fits, then a uniform position within
    // its grid cell, so the dates form a smooth distribution instead of spikes.
    //
    // A calibrated row can contain exact zeros (the gaps between the peaks of a
    // plateau date). log(0) is -inf, which categorical_logit_rng does not
    // accept, so weights are floored 700 log units below the row maximum. That
    // is a probability of about 1e-304, so those years are still never drawn.
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
}
