#pragma once

#include <cmath>

namespace range_mi {
namespace p_not_measured {

template<unsigned int dimension>
void line(
    const unsigned int * const line,
    const double * const vacancy,
    const double * const width,
    unsigned int num_cells,
    double dtheta,
    double * const output_) {

  double width_sum = 0;
  double miss_p_product = 1;
  double pnm = 0;
  for (unsigned int i = 0; i < num_cells; i++) {
    unsigned int j = line[i];

    // Scale the probability by r^d to account
    // for radial overlap as well as the width
    // to account for aliasing.
    output_[j] +=
      width[i] * dtheta *
      std::pow(width_sum, dimension - 1) *
      pnm;

    // The probability of the next cell not being
    // measured is the probability some previous
    // combination of cells were missed followed
    // by a hit.
    double miss_p = std::pow(vacancy[j], width[i]);
    pnm += miss_p_product * (1 - miss_p);
    width_sum += width[i];
    miss_p_product *= miss_p;
  }
}

}}
