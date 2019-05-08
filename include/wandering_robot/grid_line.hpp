#pragma once

#include <random>
#include <algorithm>

namespace wandering_robot {

class GridLine {

  public:

    GridLine() {}

    GridLine(unsigned int height_, unsigned int width_)
      : height(height_), width(width_) {

      // Initialize randomness
      std::random_device random_device;
      gen = std::mt19937(random_device());
      dist = std::uniform_real_distribution<double>(0., 1.);
    }

    void draw(
        double col,
        double row,
        double theta,
        unsigned int * const line,
        double * const widths,
        unsigned int & num_cells) const;

    void draw(
        unsigned int cell,
        double theta,
        unsigned int * const line,
        double * const widths,
        unsigned int & num_cells);

    void sample_regularly(
        double & x,
        double & y,
        double & theta,
        double & spatial_interpolation,
        double & angular_interpolation,
        unsigned int spatial_jitter,
        unsigned int num_beams) const;

    unsigned int size() {
      return 2 * std::max(height, width);
    }

    unsigned int height, width;

  private:
    std::mt19937 gen;
    std::uniform_real_distribution<double> dist;
};

}
