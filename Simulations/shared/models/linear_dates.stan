// Linear trend when the date of each find is uncertain.
// Each find has a set of candidate years, each with a probability. The
// likelihood of a value is averaged over the candidate years, weighted by
// their probability. The same model is used for every method:
//   midpoint or median   one candidate year with probability 1
//   full distribution    every year on the grid, with the calibrated
//                        probability (Case 1) or a flat window (Cases 2-4)
data {
  int n;                        // number of finds
  int n_years;                  // candidate years per find
  vector[n] y;                  // measured values
  matrix[n, n_years] year;      // candidate years of each find
  matrix[n, n_years] log_p;     // log probability of each candidate year
  array[n] int first;           // each find only uses the candidate years
  array[n] int last;            // from column first to column last
  real centre;                  // middle of the study period
}
parameters {
  real a;                       // value at the centre of the period
  real slope100;                // change in value per 100 years (per year the
                                // slope is so small that sampling is slow)
  real<lower=0> sigma;          // noise around the trend
}
model {
  a ~ normal(0, 10);
  slope100 ~ normal(0, 5);
  sigma ~ exponential(1);
  for (i in 1:n) {
    row_vector[last[i] - first[i] + 1] resid =
      (y[i] - a - slope100 * (year[i, first[i]:last[i]] - centre) / 100) / sigma;
    target += log_sum_exp(log_p[i, first[i]:last[i]] - log(sigma) - 0.5 * square(resid));
  }
}
generated quantities {
  real slope = slope100 / 100;          // change in value per year
  real intercept = a - slope * centre;  // value at year 0

  // One plausible date per find and per draw: a candidate year picked with
  // weight (its probability x how well the trend fits there). Not needed for the
  // recovery study and slow to save, so off; uncomment for real data or figures.
  // Years with probability 0 (gaps between plateau peaks) are floored, because
  // categorical_logit_rng does not accept -inf.
  // vector[n] date;
  // for (i in 1:n) {
  //   row_vector[last[i] - first[i] + 1] resid =
  //     (y[i] - a - slope100 * (year[i, first[i]:last[i]] - centre) / 100) / sigma;
  //   vector[last[i] - first[i] + 1] w = (log_p[i, first[i]:last[i]] - 0.5 * square(resid))';
  //   w = fmax(w, max(w) - 700);
  //   date[i] = year[i, first[i] - 1 + categorical_logit_rng(w)];
  // }
}
