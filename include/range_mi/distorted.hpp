#pragma once

#include <cmath>
#include <array>

namespace range_mi {
namespace distorted {

double normal_pdf(
    double x,
    double mean,
    double std_dev) {
  double x_norm = (x - mean)/std_dev;
  return std::exp(-0.5 * x_norm * x_norm)/(std_dev * std::sqrt(2 * M_PI));
}

void range_pdf(
    const unsigned int * const line,
    const double * const vacancy,
    const double * const width,
    unsigned int num_cells,
    double noise_dev,
    double noise_half_width,
    double integration_step,
    double pdf_width,
    double * const pdf,
    unsigned int & pdf_size) {

  double pdf_decay = 1;
  double width_sum = 0;
  pdf_size = 0;

  // Iterate over the cells
  for (unsigned int line_index = 0; line_index < num_cells; line_index++) {
    // If the width is too large the cells are of no use
    if (width_sum > pdf_width + 2 * noise_half_width) break;

    unsigned int map_index = line[line_index];
    double v = vacancy[map_index];
    double w = width[line_index];

    // Compute the negative log of the vacancy
    double l;
    if (v <= 0) {
      l = 9999999999999999;
    } else {
      l = -std::log(v);
    }

    // Pre-compute the probability of missing
    // the constant region
    double p_miss = std::exp(-l * w);

    // The probability mass of a hit
    double mass = 1 - p_miss;

    // The center of mass
    // (1 - exp(-l*w))/2 = 1 - exp(-l*com)
    double center_of_mass;
    if (l > 0) {
      center_of_mass = -std::log(0.5 * (1 + p_miss))/l;
    } else {
      center_of_mass = 0;
    }

    double z = -noise_half_width;
    unsigned int z_index = width_sum/integration_step;
    while (z < w + noise_half_width and z + width_sum < pdf_width + noise_half_width) {

      // The value of the convolution is the normal
      // centered at the center of mass, multiplied
      double value = pdf_decay * mass * normal_pdf(z, center_of_mass, noise_dev);

      // If this is the first time the cell
      // has been touched, clear it
      if (z_index == pdf_size) {
        pdf[z_index] = 0;
        pdf_size++;
      }

      // Accumulate
      pdf[z_index] += value;

      z += integration_step;
      z_index++;
    }

    // Update width and decay
    width_sum += w;
    pdf_decay *= p_miss;
  }
}

template <unsigned int dimension>
void line(
    const unsigned int * const line,
    const double * const vacancy,
    const double * const width,
    unsigned int num_cells,
    double noise_dev,
    double noise_half_width,
    double integration_step,
    double dtheta,
    double * const pdf,
    double * const output) {

  // Pre-compute factorials
  std::array<unsigned int, dimension + 1> factorial;
  factorial[0] = 1;
  for (unsigned int k = 1; k <= dimension; k++) {
    factorial[k] = factorial[k - 1] * k;
  }

  // Pre-compute E[N^kI(N)]
  std::array<double, dimension> distk_info_noise;
  switch(dimension) {
    case 3:
      distk_info_noise[2] =
        noise_dev * noise_dev *
        std::log(noise_dev * std::sqrt(2 * M_PI * std::pow(M_E,3)));
      [[fallthrough]];
    case 2:
      distk_info_noise[1] = 0;
      [[fallthrough]];
    case 1:
      distk_info_noise[0] = std::log(noise_dev * std::sqrt(2 * M_PI * M_E));
  }

  // Initialize the expected values
  // E[Z^kI(Z^k)] to E[N^kI(N^k)]
  std::array<double, dimension> distk_info;
  for (unsigned int k = 0; k < dimension; k++) {
    distk_info[k] = distk_info_noise[k];
  }

  // Initialize the expected values E[R^k]
  std::array<double, dimension> distk;
  distk[0] = 1;
  for (unsigned int k = 1; k < dimension; k++) {
    distk[k] = 0;
  }

  // Also initialize a space for gamma values
  std::array<double, dimension + 1> gamma;
  // and powers of w
  std::array<double, dimension> w_to_the;
  // And negative powers of l
  std::array<double, dimension> l_to_the_neg;

  // Iterate backwards over the cells in the line
  for (int line_index = num_cells - 1; line_index >= 0; line_index--) {

    // Fetch the vacancy
    // and width of the current cell
    unsigned int map_index = line[line_index];
    double v = vacancy[map_index];
    double w = width[line_index];

    // Compute the negative log of the vacancy
    double l;
    if (v <= 0) {
      l = 9999999999999999;
    } else {
      l = -std::log(v);
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

    // Precompute the density functions
    unsigned int hit_pdf_size;
    range_pdf(
        line + line_index,
        vacancy,
        width + line_index,
        num_cells - line_index,
        noise_dev,
        noise_half_width,
        integration_step,
        w,
        pdf,
        hit_pdf_size);
    unsigned int miss_pdf_size;
    range_pdf(
        line + line_index + 1,
        vacancy,
        width + line_index + 1,
        num_cells - line_index - 1,
        noise_dev,
        noise_half_width,
        integration_step,
        0,
        pdf + hit_pdf_size,
        miss_pdf_size);

    // Update the expected values
    for (int k = dimension - 1; k >= 0; k--) {

      double miss_distk_info = 0;
      double miss_distk = 0;
      for (int i = 0; i <= k; i++) {
        unsigned int binom =
          factorial[k]/(factorial[i] * factorial[k - i]);
        miss_distk_info +=
          binom * w_to_the[k - i] * (distk_info[i] + w * l * distk[i]);
        if (i == 2) {
          miss_distk_info +=
            binom * w_to_the[k - i] * w * l * noise_dev * noise_dev;
        }
        miss_distk +=
          binom * w_to_the[k - i] * distk[i];
      }
      distk_info[k] = p_miss * miss_distk_info;

      // Perform the numerical integration
      for (unsigned int i = 0; i < hit_pdf_size; i++) {
        double r = integration_step * i - noise_half_width;
        if (pdf[i] > 0) {
          distk_info[k] -=
            pdf[i] * std::pow(r, k) *
            std::log(pdf[i]) * integration_step;
        }
      }
      for (unsigned int i = 0; i < miss_pdf_size; i++) {
        double r = integration_step * i - noise_half_width;
        if (pdf[hit_pdf_size + i] > 0) {
          distk_info[k] +=
            p_miss * pdf[hit_pdf_size + i] * std::pow(r + w, k) *
            (std::log(pdf[hit_pdf_size + i]) - l * w) * integration_step;
        }
      }

      distk[k] = p_miss * miss_distk;
      distk[k] += l_to_the_neg[k] * gamma[k];
    }

    // MI += E[R^(n-1)I(R)]dtheta
    double mi = distk_info[dimension - 1] * dtheta;
    // MI -= E[I(N)]E[R^(n-1)]dtheta
    mi -= distk_info_noise[0] * distk[dimension - 1] * dtheta;
    if (dimension == 3) {
      // MI -= E[N^2I(N)]dtheta
      mi -= distk_info_noise[2] * dtheta;
    }
    if (mi > 0) output[map_index] += mi;
  }
}

}}
