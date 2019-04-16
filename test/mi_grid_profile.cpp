#include <vector>

#include <wandering_robot/bresenham.hpp>
#include <wandering_robot/mutual_information.hpp>

#include "random_inputs.hpp"

// Map constants
double poisson_rate = 0.1234;
unsigned int height = 500;
unsigned int width = 500;
unsigned int num_iterations = 10000;

int main() {
  // Initialize states and probabilities
  std::vector<OccupancyState> states(height * width);
  std::vector<double> p_not_measured(height * width);
  std::vector<double> mi(height * width, 0);

  // Make random states
  random_states(states);
  random_probabilities(p_not_measured);

  // Initialize the beam sampler and mutual information
  wandering_robot::MutualInformation mi_(poisson_rate);
  wandering_robot::Bresenham bresenham(height, width);

  // Initialize the line
  unsigned int line_length = std::max(height, width);
  std::vector<unsigned int> line(line_length);

  double x, y, theta;
  for (unsigned int i = 0; i < num_iterations; i++) {
    // Randomly sample a point
    bresenham.sample(x, y, theta);

    // Compute Bresenham's line
    unsigned int num_cells;
    bresenham.line(
        y, x, theta,
        line.data(),
        num_cells);

    // Make vectors of the states, etc.
    mi_.d2_grid(
        states.data(),
        p_not_measured.data(),
        line.data(),
        theta,
        num_cells,
        mi.data());
  }
}
