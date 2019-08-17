#pragma once

#include <cmath>
#include <array>

namespace range_mi {
namespace barely_distorted {

const double noise_l = 9e100;

template <unsigned int dimension>
void line(
    const unsigned int * const line,
    const double * const vacancy,
    const double * const p_not_measured,
    const double * const width,
    unsigned int num_cells,
    double dtheta,
    double * const output) {

  // Precompute noise quantities
  double neg_log_noise_l = -std::log(noise_l);

  // Precompute factorials
  std::array<unsigned int, dimension + 1> factorial;
  factorial[0] = 1;
  for (unsigned int k = 1; k <= dimension; k++) {
    factorial[k] = factorial[k - 1] * k;
  }

  // Initialize the expected values
  // E[R^kI(R)] to E[N^kI(N)]
  std::array<double, dimension> distk_info;
  double noise_l_inv = 1;
  for (unsigned int k = 0; k < dimension; k++) {
    distk_info[k] =
      noise_l_inv * (factorial[k+1] + factorial[k] * neg_log_noise_l);
    noise_l_inv /= noise_l;
  }

  // Initialize the expected values
  // E[R^k] to E[N^k]
  std::array<double, dimension> distk;
  noise_l_inv = 1;
  for (unsigned int k = 0; k < dimension; k++) {
    distk[k] = noise_l_inv * factorial[k];
    noise_l_inv /= noise_l;
  }

  // Also initialize a space for gamma values
  std::array<double, dimension + 1> gamma;
  // and powers of w
  std::array<double, dimension> w_to_the;
  // And negative powers of l
  std::array<double, dimension> l_to_the_neg;

  // Iterate backwards over the cells in the line
  for (int line_index = num_cells - 1; line_index >= 0; line_index--) {

    // Fetch the vacancy, p_not_measured,
    // and width of the current cell
    unsigned int map_index = line[line_index];
    double v = vacancy[map_index];
    double pnm = p_not_measured[map_index];
    double w = width[line_index];

    // Compute the negative log of the vacancy
    double l;
    if (v <= 0) {
      l = noise_l;
    } else {
      l = -std::log(v);
    }

    // Clip the vacancy to the maximum
    l = std::min(l, noise_l);

    // Then compute the negative log of l
    double neg_log_l;
    if (l <= 0) {
      // When l = 0 this terms is zeroed out by gammas
      neg_log_l = 0;
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
      if (l > 0) {
        l_to_the_neg[k] = l_to_the_neg[k - 1]/l;
      } else {
        // Watch out for division by zero
        // When l = 0 this is zeroed out by gammas
        l_to_the_neg[k] = 1;
      }
    }

    // Update the expected values
    for (int k = dimension - 1; k >= 0; k--) {

      double miss_distk_info = 0;
      double miss_distk = 0;
      for (int i = 0; i <= k; i++) {
        unsigned int binom =
          factorial[k]/(factorial[i] * factorial[k - i]);
        miss_distk_info +=
          binom * w_to_the[k - i] * (distk_info[i] + pnm * w * l * distk[i]);
        miss_distk +=
          binom * w_to_the[k - i] * distk[i];
      }
      distk_info[k] = p_miss * miss_distk_info;
      distk_info[k] +=
        pnm * l_to_the_neg[k] * (gamma[k + 1] + neg_log_l * gamma[k]);
      distk_info[k] += (1 - pnm) * l_to_the_neg[k] * gamma[k] * (1 + neg_log_noise_l);

      distk[k] = p_miss * miss_distk;
      distk[k] += l_to_the_neg[k] * gamma[k];
    }

    // MI += E[R^(n-1)I(R)]dtheta
    output[map_index] += distk_info[dimension - 1] * dtheta;
    // MI -= (1 - log Lambda) E[R^(n-1)]dtheta
    output[map_index] -= (1 + neg_log_noise_l) * distk[dimension - 1] * dtheta;
  }
}

}}
