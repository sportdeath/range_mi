#pragma once

#include <algorithm>

namespace range_entropy {

class GridLine {

  public:

    GridLine() {}

    GridLine(unsigned int height_, unsigned int width_,
        unsigned int spatial_jitter_, unsigned int num_beams_)
      : height(height_), width(width_),
        spatial_jitter(spatial_jitter_), num_beams(num_beams_)
    {}

    void draw(
        double x,
        double y,
        double theta,
        unsigned int * const line,
        double * const widths,
        unsigned int & num_cells) const;

    void sample(
        double & x,
        double & y,
        double & theta,
        double & spatial_interpolation,
        double & angular_interpolation) const;

    unsigned int size() const {
      return 2 * std::max(height, width);
    }

    unsigned int height, width;
    unsigned int spatial_jitter, num_beams;
};

}
