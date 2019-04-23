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
}

void wandering_robot::GridWanderer::iterate_mi() {
  // Randomly sample a point
  grid_line.sample(x, y, theta);

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
