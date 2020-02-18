#pragma once

#include <vector>

namespace range_mi {

class GridMI {

private:
  // The computed mutual information
  std::vector<double> mi_;
  std::vector<double> p_not_measured_;
  std::vector<double> p_not_measured_single;

  // Map parameters
  unsigned int height, width;

  // Noise parameters
  double noise_dev, noise_half_width, integration_step;
  std::vector<double> pdf;

  // Temporary storage for lines
  std::vector<unsigned int> line;
  std::vector<double> widths;
  double x, y;
  unsigned int num_cells;

public:
  GridMI() {}

  GridMI(
      unsigned int height,
      unsigned int width,
      double noise_dev=0,
      double noise_half_width=0,
      double integration_step=0);

  void compute_mi_beam(
      const double * const vacancy,
      double theta,
      double dtheta,
      double & spatial_interpolation);

  void compute_mi(
      const double * const vacancy,
      unsigned int num_beams);

  void condition(
      const double * const vacancy,
      double x,
      double y,
      double theta_min,
      double theta_max,
      double dtheta);

  void reset_mi() {
    std::fill(mi_.begin(), mi_.end(), 0);
  }
  void reset_p_not_measured() {
    std::fill(p_not_measured_.begin(), p_not_measured_.end(), 1);
  }

  const std::vector<double> & mi() const {
    return mi_;
  }
  const std::vector<double> & p_not_measured() const {
    return p_not_measured_;
  }

  static const unsigned int dimension = 2;
  static const bool lower_bound = true;
};

}
