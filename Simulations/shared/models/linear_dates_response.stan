// Linear trend with time as the response: date = intercept + slope * x + noise,
// where x is a value known exactly and the date of each find is uncertain.
// The same idea as linear_dates.stan, the other way round: the likelihood of each
// find is averaged over its candidate years, weighted by their probability.
//   midpoint or median   one candidate year with probability 1
//   full distribution    every year on the grid, with the calibrated
//                        probability (Case 1) or a flat window (Cases 2-4)
// For the full distribution each grid year stands for a cell of the grid
// (half_cell years either side), and the trend scores it by the share of its
// normal distribution inside the cell. Scoring by the density at the grid year
// lets sigma collapse to zero on a grid year and the likelihood blow up.
// Point dates have no grid (half_cell = 0) and use the density.
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
  real<lower=0> half_cell;      // half the grid step in years; 0 for point dates
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
    // distance of each candidate year from the trend, in centuries
    row_vector[last[i] - first[i] + 1] gap =
      (year[i, first[i]:last[i]] - centre) / 100 - a - b * (x[i] - x_centre);
    if (half_cell > 0) {
      // The share is the same above or below the trend, so it is computed on
      // the lower side (-abs), where both Phi values are small and precise.
      // Above the trend both would be close to 1 and their difference noisy.
      row_vector[last[i] - first[i] + 1] cell =
        Phi((-abs(gap) + half_cell / 100) / s) - Phi((-abs(gap) - half_cell / 100) / s);
      target += log_sum_exp(log_p[i, first[i]:last[i]] + log(fmax(cell, 1e-300)));
    } else {
      target += log_sum_exp(log_p[i, first[i]:last[i]] - log(s) - 0.5 * square(gap / s));
    }
  }
}
generated quantities {
  real slope = 100 * b;                                // years per unit of x
  real sigma = 100 * s;                                // noise, in years
  real intercept = centre + 100 * (a - b * x_centre);  // date at x = 0
}
