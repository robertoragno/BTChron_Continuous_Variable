// Linear trend with time as the response: date = intercept + slope * x + noise,
// where x is a value known exactly and the date of each find is uncertain.
// The same idea as linear_dates.stan, the other way round: the likelihood of each
// find is averaged over its candidate years, weighted by their probability.
//   midpoint or median   one candidate year with probability 1
//   full distribution    every year on the grid, with the calibrated
//                        probability (Case 1) or a flat window (Cases 2-4)
data {
  int n;                        // number of finds
  int n_years;                  // candidate years per find
  vector[n] x;                  // measured values, known exactly
  matrix[n, n_years] year;      // candidate years of each find
  matrix[n, n_years] log_p;     // log probability of each candidate year
  array[n] int first;           // each find only uses the candidate years
  array[n] int last;            // from column first to column last
  real centre;                  // middle of the study period (a year)
  real x_centre;                // middle of the measured values
}
parameters {
  // Dates are counted in centuries from the centre of the period, so the
  // parameters are on similar scales and sampling is fast
  real a;                       // date at x_centre
  real b;                       // change in date per unit of x
  real<lower=0> s;              // noise around the trend
}
model {
  a ~ normal(0, 5);
  b ~ normal(0, 5);
  s ~ exponential(1);
  for (i in 1:n) {
    row_vector[last[i] - first[i] + 1] resid =
      ((year[i, first[i]:last[i]] - centre) / 100 - a - b * (x[i] - x_centre)) / s;
    target += log_sum_exp(log_p[i, first[i]:last[i]] - log(s) - 0.5 * square(resid));
  }
}
generated quantities {
  real slope = 100 * b;                                // years per unit of x
  real sigma = 100 * s;                                // noise, in years
  real intercept = centre + 100 * (a - b * x_centre);  // date at x = 0
}
