#include "wandering_robot/mutual_information.hpp"
#include "wandering_robot/grid_line.hpp"
#include "wandering_robot/grid_wanderer.hpp"

wandering_robot::GridWanderer::GridWanderer(
    unsigned int height,
    unsigned int width,
    double poisson_rate,
    bool beam_independence_)
  : beam_independence(beam_independence_) {

  // Store the map states
  states_ = std::vector<wandering_robot::OccupancyState>(
      height * width,
      wandering_robot::OccupancyState::unknown);

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

  mi_computer = wandering_robot::MutualInformation(poisson_rate);
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
        states_.data(),
        widths.data(),
        p_not_measured_.data(),
        line.data(),
        num_cells,
        mi_.data());
  } else {
    mi_computer.d1(
        states_.data(),
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
      if (states_[line[j]] == OccupancyState::unknown)
        unknown_length += widths[j];
    }

    theta += theta_step;
  }

  // Given the distances apply the conditioning formula
  for (unsigned int i = 0; i < p_not_measured_.size(); i++) {
    mi_computer.condition(p_not_measured_[i], condition_distances[i]);
  }
}

void wandering_robot::GridWanderer::apply_scan(
    double x, double y, const std::vector<double> & scan) {
  double theta = -M_PI;

  for (size_t i = 0; i < scan.size(); i++) {
    // Draw a line in the direction of each scan
    grid_line.draw(
        x, y, theta,
        line.data(),
        widths.data(),
        num_cells);

    // Iterate over the line
    double length = 0;
    for (unsigned int j = 0; j < num_cells; j++) {
      // If we've reached the end of the beam mark occupied
      if (length >= scan[i]) {
        states_[line[j]] = wandering_robot::OccupancyState::occupied;
        break;
      }

      // Otherwise mark free
      states_[line[j]] = wandering_robot::OccupancyState::free;

      // Accumulate the length
      length += widths[j];
    }

    theta += 2 * M_PI/scan.size();
  }
}

std::vector<double> wandering_robot::GridWanderer::make_scan(
    double x, double y,
    unsigned int num_beams) {

  std::vector<double> scan(num_beams);
  double theta = -M_PI;
  for (unsigned int i = 0; i < num_beams; i++) {
    // Draw a line in the direction of each scan
    grid_line.draw(
        x, y, theta,
        line.data(),
        widths.data(),
        num_cells);

    double length = 0;
    for (unsigned int j = 0; j < num_cells; j++) {
      if (map[line[j]] == wandering_robot::OccupancyState::occupied)
        break;

      length += widths[j];
    }

    scan[i] = length;
    theta += (2 * M_PI)/num_beams;
  }

  return scan;
}
