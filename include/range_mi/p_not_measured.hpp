#pragma once

namespace range_mi {
namespace p_not_measured {

void line(
    const unsigned int * const line,
    const double * const vacancy,
    const double * const width,
    unsigned int num_cells,
    unsigned int dimension,
    double dtheta,
    double * const output_);

}}
