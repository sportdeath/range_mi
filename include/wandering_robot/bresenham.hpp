#pragma once

#include <random>

namespace wandering_robot {

class Bresenham {

  public:

    Bresenham() {}

    Bresenham(unsigned int height_, unsigned int width_)
      : height(height_), width(width_) {

      // Initialize randomness
      std::random_device random_device;
      gen = std::mt19937(random_device());
      dist = std::uniform_real_distribution<double>(0., 1.);
    }

    void line(
        double row,
        double col,
        double theta,
        unsigned int * const line,
        unsigned int & num_cells);

    void sample(
        double & x, 
        double & y, 
        double & theta);

  private:
    unsigned int height, width;
    std::mt19937 gen;
    std::uniform_real_distribution<double> dist;
};

}
