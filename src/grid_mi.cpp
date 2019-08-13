#include <cmath>

#include "range_mi/barely_distorted.hpp"
#include "range_mi/p_not_measured.hpp"
#include "range_mi/grid_line.hpp"
#include "range_mi/grid_mi.hpp"

range_mi::GridMI::GridMI(
    unsigned int height_,
    unsigned int width_) :
  height(height_),
  width(width_) {

  // Set the mutual information to zero
  mi_ = std::vector<double>(height * width, 0);

  // Initialize the conditioning map
  p_not_measured_ = std::vector<double>(height * width, 1);

  // Initialize the storage vectors
  line = std::vector<unsigned int>(2 * std::max(height, width));
  widths = std::vector<double>(2 * std::max(height, width));
  p_not_measured_single = std::vector<double>(height * width, 0);
}

void range_mi::GridMI::compute_mi_beam(
    const double * const vacancy,
    double theta,
    double dtheta,
    double & spatial_interpolation) {

  // Convert the interpolation parameters to
  // x, y, theta
  range_mi::grid_line::sample(
      height, width,
      theta, x, y,
      spatial_interpolation);

  // Compute the intersections of
  // the line with the grid
  range_mi::grid_line::draw(
      height, width,
      x, y, theta,
      line.data(),
      widths.data(),
      num_cells);

  // Accumulate the mutual information
  range_mi::barely_distorted::line<dimension>(
      line.data(),
      vacancy,
      p_not_measured_.data(),
      widths.data(),
      num_cells,
      dtheta,
      mi_.data());
}

void range_mi::GridMI::condition(
    const double * const vacancy,
    double x,
    double y,
    double theta_min,
    double theta_max,
    double dtheta) {

  // Clear the old p_measured
  std::fill(p_not_measured_single.begin(), p_not_measured_single.end(), 0);

  // Compute the new one for the given point
  for (double theta = theta_min; theta < theta_max; theta+=dtheta) {
    // Draw a line
    range_mi::grid_line::draw(
        height, width,
        x, y, theta,
        line.data(),
        widths.data(),
        num_cells);

    // Accumulate p_not_measured along the line
    range_mi::p_not_measured::line<dimension>(
        line.data(),
        vacancy,
        widths.data(),
        num_cells,
        dtheta,
        p_not_measured_single.data());
  }

  // Update the probabilities
  for (unsigned int i = 0; i < p_not_measured_.size(); i++) {
    p_not_measured_[i] = std::min(p_not_measured_[i], p_not_measured_single[i]);
  }
}
