#pragma once

namespace range_mi {
namespace barely_distorted {

void line(
    const unsigned int * const line,
    const double * const vacancy,
    const double * const p_not_measured,
    const double * const width,
    unsigned int num_cells,
    double dtheta,
    double noise_l,
    unsigned int dimension,
    double * const output);

}}
