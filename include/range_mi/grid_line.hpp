#pragma once

namespace range_mi {
namespace grid_line {

void draw(
    unsigned int height,
    unsigned int widht,
    double x,
    double y,
    double theta,
    unsigned int * const line,
    double * const widths,
    unsigned int & num_cells);

void sample(
    unsigned int height,
    unsigned int width,
    double theta,
    double & x,
    double & y,
    double & spatial_interpolation);

}}
