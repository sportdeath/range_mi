#include <random>
#include <vector>
#include <iostream>

#include <wandering_robot/occupancy_state.hpp>
#include <wandering_robot/mutual_information.hpp>

using namespace wandering_robot;

// Initialize constants
double poisson_rate = 1.;
int num_iterations = 100;
size_t num_cells = 1000;
double p_occupied = 0.05;
double p_unknown = 0.5;

// Initialize random generator
std::random_device random_device;
std::mt19937 gen(random_device());
std::uniform_real_distribution<double> dist(0.,1.);

void random_cells(
    std::vector<OccupancyState> & states,
    std::vector<double> & widths,
    std::vector<double> & p_not_measured) {

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

    // Make the cell widths random numbers in [0,1]
    widths[i] = dist(gen);

    // Make the probability that a cell has not already
    // been measured uniform in [0,1]
    p_not_measured[i] = dist(gen);
  }
}

int main() {
  // Initialize the inputs and outputs
  std::vector<OccupancyState> states(num_cells);
  std::vector<double> widths(num_cells);
  std::vector<double> p_not_measured(num_cells);
  std::vector<double> mi_1d(num_cells);

  // Initialize the mutual information
  MutualInformation mi(poisson_rate);

  for (int i = 0; i < num_iterations; i++) {
    // Generate random cells
    random_cells(states, widths, p_not_measured);

    // Compute the mutual information all at once
    mi.d1(states.data(), widths.data(), p_not_measured.data(), num_cells, mi_1d.data());

    // Compute the mutual information one at a time
    for (unsigned int j = 0; j < num_cells; j++) {
      double mi_1d_cell = mi.d1(states.data() + j, widths.data() + j, p_not_measured.data() + j, num_cells);

      double error = std::abs(mi_1d_cell - mi_1d[j]);
      if (error > 0.000001) {
        std::cout <<
          "MI from single computation: " << mi_1d_cell <<
          ", MI from joint computation " << mi_1d[j] << std::endl;
        return -1;
      }
    }
  }

  return 0;
}
