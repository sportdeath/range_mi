#include <cmath>

#include "range_entropy/expected.hpp"
#include "range_entropy/grid_line.hpp"
#include "range_entropy/grid_expected.hpp"

using namespace range_entropy;

range_entropy::GridExpected::GridExpected(
    unsigned int height,
    unsigned int width,
    unsigned int spatial_jitter,
		unsigned int num_beams) {

  // Set the expected value to zero
  surface_ = std::vector<double>(height * width, 0);

  // Initialize the conditioning map
  p_not_measured_ = std::vector<double>(height * width, 1);

  // Initialize the beam sampler
  grid_line = GridLine(height, width, spatial_jitter, num_beams);

  // Initialize the storage vectors
  line = std::vector<unsigned int>(grid_line.size());
  widths = std::vector<double>(grid_line.size());
  p_not_measured_single = std::vector<double>(height * width, 0);
}

void range_entropy::GridExpected::compute_surface(
    const double * const p_free,
    void f(const unsigned int * const,
    const double * const,
    const double * const,
    const double * const,
    unsigned int,
    double * const)) {

  reset_surface();

  // Iterate over the steps
  double spatial_interpolation = 0;
  double angular_interpolation = 0;
  while (angular_interpolation < 1)
    compute_surface_beam(
        spatial_interpolation,
        angular_interpolation,
        p_free, f);
}

void range_entropy::GridExpected::compute_surface_beam(
    double & spatial_interpolation,
    double & angular_interpolation,
    const double * const p_free,
    void f(const unsigned int * const,
    const double * const,
    const double * const,
    const double * const,
    unsigned int,
    double * const)) {

  // Convert the interpolation parameters to
  // x, y, theta
  grid_line.sample(
      x, y, theta,
      spatial_interpolation,
      angular_interpolation);

  // Compute the intersections of
  // the line with the grid
  grid_line.draw(
      x, y, theta,
      line.data(),
      widths.data(),
      num_cells);

  // Accumulate value along the line
  f(line.data(),
    p_free,
    p_not_measured_.data(),
    widths.data(),
    num_cells,
    surface_.data());
}

void range_entropy::GridExpected::condition(
    double x,
    double y,
    const double * const vacancy) {

  // Clear the old p_measured
  std::fill(p_not_measured_single.begin(), p_not_measured_single.end(), 0);

  // Compute the new one for the given point
  double theta = 0;
  double theta_step = 2 * M_PI /((double) 100 * grid_line.num_beams);
  for (unsigned int i = 0; i < 100 * grid_line.num_beams; i++) {
    // Draw a line
    grid_line.draw(
        x, y, theta,
        line.data(),
        widths.data(),
        num_cells);

    // Accumulate p_measured along the line
    range_entropy::expected::p_not_measured(
        line.data(),
        vacancy,
        widths.data(),
        num_cells,
        2,
        p_not_measured_single.data());

    theta += theta_step;
  }

  // Update the probabilities
  for (unsigned int i = 0; i < p_not_measured_.size(); i++) {
    p_not_measured_[i] *= p_not_measured_single[i] * theta_step;
  }
}
