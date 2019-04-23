#include <vector>

#include <wandering_robot/grid_wanderer.hpp>

#include "random_inputs.hpp"

// Map constants
double poisson_rate = 0.1234;
unsigned int height = 500;
unsigned int width = 500;
unsigned int num_iterations = 10000;

int main() {
  // Initialize the wanderer
  wandering_robot::GridWanderer w(poisson_rate);

  // Initialize random states
  std::vector<OccupancyState> states(height * width);
  random_states(states);
  w.set_map(states, height, width);

  for (unsigned int i = 0; i < num_iterations; i++) {
    w.iterate_mi();
  }
}
