#pragma once

#include <vector>

#include "range_entropy/grid_line.hpp"

namespace range_entropy {

class GridExpected {

private:
  // The computed surface
  std::vector<double> surface_;
  std::vector<double> p_not_measured_;
  std::vector<double> p_not_measured_single;

  // A computing device
  GridLine grid_line;

  // Temporary storage
  std::vector<unsigned int> line;
  std::vector<double> widths;
  double x, y, theta;
  unsigned int num_cells;

public:
  GridExpected() {}

  GridExpected(
      unsigned int height,
      unsigned int width,
      unsigned int spatial_jitter,
      unsigned int num_beams);

  void compute_surface(
      const double * const p_free,
      bool information,
      unsigned int dimension,
      double noise_dev=0,
      double noise_width=0,
      double noise_step_size=0);

  void compute_surface_beam(
      double & spatial_interpolation,
      double & angular_interpolation,
      const double * const p_free,
      bool information,
      unsigned int dimension,
      double noise_dev=0,
      double noise_width=0,
      double noise_step_size=0);

  void condition(
      double x,
      double y,
      const double * const vacancy,
      unsigned int condition_steps,
      unsigned int dimension);

  void reset_surface() {
    std::fill(surface_.begin(), surface_.end(), 0);
  }
  void reset_p_not_measured() {
    std::fill(p_not_measured_.begin(), p_not_measured_.end(), 1);
  }

  const std::vector<double> & surface() const {
    return surface_;
  }
  const std::vector<double> & p_not_measured() const {
    return p_not_measured_;
  }
};

}
