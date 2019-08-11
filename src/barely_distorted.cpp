#include <cmath>
#include <vector>

#include "range_mi/barely_distorted.hpp"

void range_mi::barely_distorted::line(
    const unsigned int * const line,
    const double * const vacancy,
    const double * const p_not_measured,
    const double * const width,
    unsigned int num_cells,
    double dtheta,
    double noise_l,
    unsigned int dimension,
    double * const output) {

  // Precompute noise quantities
  double neg_log_noise_l = -std::log(noise_l);

  // Precompute factorials
  std::vector<unsigned int> factorial(dimension);
  factorial[0] = 1;
  for (unsigned int k = 1; k < dimension; k++) {
    factorial[k] = factorial[k - 1] * k;
  }

  // Precompute the expected values
  // E[N^kI(N)]
  std::vector<double> distk_info_noise(dimension);
  double noise_l_inv = 1;
  for (unsigned int k = 0; k < dimension; k++) {
    distk_info_noise[k] =
      -noise_l_inv * (factorial[k+1] + factorial[k] * neg_log_noise_l);
    noise_l_inv /= noise_l;
  }

  // Initialize the expected values
  // E[R^kI(R)] to E[N^kI(N)]
  std::vector<double> distk_info(dimension);
  for (unsigned int k = 0; k < dimension; k++) {
    distk_info[k] = distk_info_noise[k];
  }

  // Initialize the expected values
  // E[N^k] to E[N^k]
  std::vector<double> distk(dimension);
  noise_l_inv = 1;
  for (unsigned int k = 0; k < dimension; k++) {
    distk[k] = noise_l_inv * factorial[k];
    noise_l_inv /= noise_l;
  }

  // Also initialize a space for gamma values
  std::vector<double> gamma(dimension + 1);
  // and powers of w
  std::vector<double> w_to_the(dimension);
  // And negative powers of l
  std::vector<double> l_to_the_neg(dimension);

  // Iterate backwards over the cells in the line
  for (int i = num_cells - 1; i >= 0; i--) {

    // Fetch the vacancy, p_not_measured,
    // and width of the current cell
    unsigned int j = line[i];
    double v = vacancy[j];
    double pnm = p_not_measured[j];
    double w = width[i];

    // Compute the negative log of the vacancy
    double l;
    if (v <= 0) {
      l = noise_l;
    } else {
      l = -std::log(v);
    }

    // Clip the vacancy to the maximum
    l = std::max(l, noise_l);

    // Then compute the negative log of l
    double neg_log_l;
    if (l <= 0) {
      neg_log_l = noise_l;
    } else {
      neg_log_l = -std::log(l);
    }

    // Precompute the probability of missing
    // the constant region
    double p_miss = std::exp(-l * w);

    // Precompute the gamma values
    gamma[0] = 1 - p_miss;
    double wl = 1;
    for (unsigned int k = 1; k <= dimension; k++) {
      // Update (wl)^k
      wl *= w * l;

      // Update the gamma function
      gamma[k] = k * gamma[k - 1] - p_miss * wl;
    }

    // Precompute powers of w
    w_to_the[0] = 1;
    for (unsigned int k = 1; k < dimension; k++) {
      w_to_the[k] = w_to_the[k - 1] * w;
    }

    // precompute
    l_to_the_neg[0] = 1;
    for (unsigned int k = 1; k < dimension; k++) {
      l_to_the_neg[k] = l_to_the_neg[k - 1]/l;
    }

    // Update the expected values
    for (int k = dimension - 1; k >= 0; k--) {

      double miss_info = 0;
      for (int i = 0; i <= k; i++) {
        unsigned int binom =
          factorial[k]/(factorial[i] * factorial[k - i]);
        miss_info +=
          binom * w_to_the[k - i] * (distk_info[i] + pnm * w * l * distk[i]);
      }
      distk_info[i] += p_miss * miss_info;
      distk_info[i] -=
        pnm * l_to_the_neg[k] * (gamma[k + 1] + neg_log_l * gamma[k]);
      distk_info[i] += (1 - pnm) * (1 - p_miss) * distk_info_noise[k];

      double miss_dist = 0;
      for (int i = 0; i <= k; i++) {
        unsigned int binom =
          factorial[k]/(factorial[i] * factorial[k - i]);
        miss_dist +=
          binom * w_to_the[k - i] * distk[i];
      }
      distk[i] += p_miss * miss_dist;
      distk[i] -= l_to_the_neg[k] * gamma[k];
    }

    // MI += E[R^(n-1)I(R)]dtheta
    output[j] += distk_info[dimension - 1] * dtheta;
    // MI -= (1 - log Lambda) E[R^(n-1)]dtheta
    output[j] -= (1 + neg_log_noise_l) * distk[dimension - 1] * dtheta;
  }
}
