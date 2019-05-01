#include <limits>
#include <algorithm>

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

void wandering_robot::GridWanderer::condition(unsigned int cell, unsigned int angular_steps) {
  // Empty the condition distances
  // These values will represent the distance 
  // through unknown space that a beam, originating
  // through cell (x, y), travels.
  std::fill(condition_distances.begin(), condition_distances.end(), 0);

  double theta = 0;
  double theta_step = (2 * M_PI)/angular_steps;

  for (unsigned int i = 0; i < angular_steps; i++) {
    // Compute the intersections of
    // the line with the grid
    grid_line.draw(
        cell, theta,
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

void wandering_robot::GridWanderer::make_scan(
    unsigned int cell,
    unsigned int num_beams) {

  std::vector<double> scan(num_beams);
  double theta = -M_PI;
  for (unsigned int i = 0; i < num_beams; i++) {
    // Draw a line in the direction of each scan
    grid_line.draw(
        cell, theta,
        line.data(),
        widths.data(),
        num_cells);

    for (unsigned int j = 0; j < num_cells; j++) {
      // If we hit an occupied cell, stop
      if (map[line[j]] != wandering_robot::OccupancyState::free) {
        states_[line[j]] = wandering_robot::OccupancyState::occupied;
        break;
      }

      // Otherwise mark free
      states_[line[j]] = wandering_robot::OccupancyState::free;
    }

    theta += (2 * M_PI)/num_beams;
  }
}

void wandering_robot::GridWanderer::dijkstra(
    unsigned int start, 
    std::vector<double> & distances,
    std::vector<unsigned int> & parents) {

  // Create arrays of distances and parents
  distances = std::vector<double>(states_.size(), std::numeric_limits<double>::max());
  parents = std::vector<unsigned int>(states_.size());

  // Initialize the nodes and set the distance to zero
  std::vector<bool> closed_set(states_.size(), false);
  std::vector<unsigned int> open_set = {start};
  distances[start] = 0;

  // Initialize a comparator and make the heap
  auto cmp = [&distances](unsigned int left, unsigned int right) mutable {
    return distances[left] > distances[right];
  };

  while (not open_set.empty()) {
    // Remake into a heap
    std::make_heap(open_set.begin(), open_set.end(), cmp);

    // Find the element with the lowest value.
    // Remove the node and don't visit it again
    unsigned int current = open_set.front();
    std::pop_heap(open_set.begin(), open_set.end(), cmp);
    open_set.pop_back();
    if (closed_set[current]) continue;
    closed_set[current] = true;

    // Convert current to x, y
    int y = current/grid_line.width;
    int x = current - y * grid_line.width;

    // Iterate around neighbors
    for (int x_ = x - 1; x_ <= x + 1; x_++) {
      for (int y_ = y - 1; y_ <= y + 1; y_++) {
        // Check if out of bounds
        if (x_ < 0 or y_ < 0 or
            x_ >= (int) grid_line.width or
            y_ >= (int) grid_line.height)
          continue;

        // Convert back to cell
        unsigned int neighbor = y_ * grid_line.width + x_;

        // Ignore, it's already been visited
        if (closed_set[neighbor]) continue;

        // Ignore if the cell is not free
        if (states_[neighbor] != OccupancyState::free) continue;

        // Update the new distance
        double new_distance = distances[current] +
          std::sqrt(std::abs(x - x_) + std::abs(y - y_));

        // If this new distance is better than before
        if (new_distance < distances[neighbor]) {
          // Update it and the parent
          distances[neighbor] = new_distance;
          parents[neighbor] = current;

          // Add it to the open set
          open_set.push_back(neighbor);
        }
      }
    }
  }
}
