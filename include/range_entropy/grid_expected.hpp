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
      void f(const unsigned int * const,
      const double * const,
      const double * const,
      const double * const,
      unsigned int,
      double * const));

  void compute_surface_beam(
      double & spatial_interpolation,
      double & angular_interpolation,
      const double * const p_free,
      void f(const unsigned int * const,
      const double * const,
      const double * const,
      const double * const,
      unsigned int,
      double * const));

  void condition(
      double x,
      double y,
      const double * const vacancy);

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
