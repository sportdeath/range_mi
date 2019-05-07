#include <limits>
#include <algorithm>

#include "wandering_robot/grid_line.hpp"
#include "wandering_robot/grid_mapper.hpp"

wandering_robot::GridMapper::GridMapper(
    unsigned int height,
    unsigned int width) {

  // Store the map states
  map_ = std::vector<wandering_robot::OccupancyState>(
      height * width,
      wandering_robot::OccupancyState::unknown);

  // Initialize the beam sampler
  grid_line = wandering_robot::GridLine(height, width);

  // Initialize the storage vectors
  line = std::vector<unsigned int>(grid_line.size());
  widths = std::vector<double>(grid_line.size());
}

void wandering_robot::GridMapper::make_scan(
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
      if (map_gt[line[j]] != wandering_robot::OccupancyState::free) {
        map_[line[j]] = wandering_robot::OccupancyState::occupied;
        break;
      }

      // Otherwise mark free
      map_[line[j]] = wandering_robot::OccupancyState::free;
    }

    theta += (2 * M_PI)/num_beams;
  }
}

void wandering_robot::GridMapper::dijkstra(
    unsigned int start, 
    std::vector<double> & distances,
    std::vector<int> & parents) const {

  // Create arrays of distances and parents
  distances = std::vector<double>(map_.size(), std::numeric_limits<double>::max());
  parents = std::vector<int>(map_.size());

  // Initialize the nodes and set the distance to zero
  std::vector<bool> closed_set(map_.size(), false);
  std::vector<unsigned int> open_set = {start};
  distances[start] = 0;
  parents[start] = -1;

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
        if (map_[neighbor] != OccupancyState::free) continue;

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

std::vector<unsigned int> wandering_robot::GridMapper::dijkstra_path(
    int end,
    const std::vector<int> & parents) const {

  std::vector<unsigned int> path;
  while (end >= 0) {
    path.push_back(end);
    end = parents[end];
  }

  std::reverse(path.begin(), path.end());

  return path;
}
