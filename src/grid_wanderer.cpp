#include "wandering_robot/mutual_information.hpp"
#include "wandering_robot/grid_line.hpp"
#include "wandering_robot/grid_wanderer.hpp"

void wandering_robot::GridWanderer::set_map(
    const std::vector<wandering_robot::OccupancyState> & states_,
    unsigned int height,
    unsigned int width) {

  // Store the map states
  states = states_;

  // Set the MI to zero
  mi_ = std::vector<double>(height * width, 0);

  // Set the probability of measurement to 1
  p_not_measured_ = std::vector<double>(height * width, 1);

  // Initialize the beam sampler
  grid_line = wandering_robot::GridLine(height, width);

  // Initialize the storage vectors
  line = std::vector<unsigned int>(grid_line.size());
  widths = std::vector<double>(grid_line.size());
  condition_distances = std::vector<double>(height * width);
}

void wandering_robot::GridWanderer::accrue_mi(
    double spatial_interpolation,
    double angular_interpolation) {

  // Randomly sample a point
  grid_line.sample_regularly(x, y, theta,
      spatial_interpolation,
      angular_interpolation);

  // Compute the intersections of
  // the line with the grid
  grid_line.draw(
      x, y, theta,
      line.data(),
      widths.data(),
      num_cells);

  // Make vectors of the states, etc.
  if (beam_independence) {
    mi_computer.d2(
        states.data(),
        widths.data(),
        p_not_measured_.data(),
        line.data(),
        num_cells,
        mi_.data());
  } else {
    mi_computer.d1(
        states.data(),
        widths.data(),
        p_not_measured_.data(),
        line.data(),
        num_cells,
        mi_.data());
  }
}

void wandering_robot::GridWanderer::condition(double x, double y, unsigned int theta_steps) {
  // Empty the condition distances
  // These values will represent the distance 
  // through unknown space that a beam, originating
  // through cell (x, y), travels.
  std::fill(condition_distances.begin(), condition_distances.end(), 0);

  double theta = 0;
  double theta_step = (2 * M_PI)/theta_steps;

  for (unsigned int i = 0; i < theta_steps; i++) {
    // Compute the intersections of
    // the line with the grid
    grid_line.draw(
        x, y, theta,
        line.data(),
        widths.data(),
        num_cells);

    double total_length = 0;
    double unknown_length = 0;
    for (unsigned int j = 0; j < num_cells; j++) {
      condition_distances[line[j]] += widths[j] * total_length * unknown_length * theta_step;
      total_length += widths[j];
      if (states[line[j]] == OccupancyState::unknown)
        unknown_length += widths[j];
    }

    theta += theta_step;
  }

  // Given the distances apply the conditioning formula
  for (unsigned int i = 0; i < p_not_measured_.size(); i++) {
    mi_computer.condition(p_not_measured_[i], condition_distances[i]);
  }
}
