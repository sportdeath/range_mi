#include <vector>

#include <wandering_robot/grid_mutual_information.hpp>

#include "random_inputs.hpp"

// Map constants
double poisson_rate = 0.1234;
unsigned int height = 500;
unsigned int width = 500;
unsigned int num_beams = 100;
unsigned int jitter = 4;
bool beam_independence = true;

int main() {
  // Initialize random states
  std::vector<OccupancyState> states(height * width);
  random_states(states);

  // Initialize the computation
  wandering_robot::GridMutualInformation mi(
      &states,
      height,
      width,
      poisson_rate,
      beam_independence);

  // Compute!
  mi.compute_mi_surface(jitter, num_beams);
}
