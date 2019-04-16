#include <random>
#include <vector>

#include <wandering_robot/occupancy_state.hpp>

using namespace wandering_robot;

double p_occupied = 0.05;
double p_unknown = 0.5;

// Initialize random generator
std::random_device random_device;
std::mt19937 gen(random_device());
std::uniform_real_distribution<double> dist(0.,1.);

void random_states(std::vector<OccupancyState> & states) {

  for (size_t i = 0; i < states.size(); i++) {
    // Construct the occupancy states
    // according to their probability
    double probability = dist(gen);
    if (probability < p_occupied) {
      states[i] = OccupancyState::occupied;
    } else if (probability < p_occupied + p_unknown) {
      states[i] = OccupancyState::unknown;
    } else {
      states[i] = OccupancyState::free;
    }
  }
}

void random_probabilities(std::vector<double> & p) {
  for (size_t i = 0; i < p.size(); i++) {
    p[i] = dist(gen);
  }
}
