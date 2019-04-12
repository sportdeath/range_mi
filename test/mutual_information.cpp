#include <random>
#include <vector>
#include <iostream>
#include <chrono>

#include <wandering_robot/occupancy_state.hpp>
#include <wandering_robot/mutual_information.hpp>

using namespace wandering_robot;

// Initialize constants
double poisson_rate = 0.1234;
int num_iterations = 1000;
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
  std::vector<double> mi_1d(num_cells), mi_2d(num_cells);

  // Initialize timing
  std::chrono::time_point<std::chrono::system_clock> start, end;
  double time_1d_joint = 0, time_1d_individual = 0;
  double time_2d_joint = 0, time_2d_individual = 0;

  // Initialize the mutual information
  MutualInformation mi(poisson_rate);

  for (int i = 0; i < num_iterations; i++) {
    // Generate random cells
    random_cells(states, widths, p_not_measured);

    // Compute the mutual information all at once

    // In 1 dimension ...
    start = std::chrono::system_clock::now();
    mi.d1(states.data(), widths.data(), p_not_measured.data(), num_cells, mi_1d.data());
    end = std::chrono::system_clock::now();
    time_1d_joint += std::chrono::duration<double>(end - start).count();

    // ... and in 2 dimensions.
    start = std::chrono::system_clock::now();
    mi.d2(states.data(), widths.data(), p_not_measured.data(), num_cells, mi_2d.data());
    end = std::chrono::system_clock::now();
    time_2d_joint += std::chrono::duration<double>(end - start).count();

    // Compute the mutual information one at a time
    for (unsigned int j = 0; j < num_cells; j++) {
      // In 1 dimension ...
      start = std::chrono::system_clock::now();
      double mi_1d_cell = mi.d1(states.data() + j, widths.data() + j, p_not_measured.data() + j, num_cells);
      end = std::chrono::system_clock::now();
      time_1d_individual += std::chrono::duration<double>(end - start).count();

      // ... and in 2 dimensions.
      start = std::chrono::system_clock::now();
      double mi_2d_cell = mi.d2(states.data() + j, widths.data() + j, p_not_measured.data() + j, num_cells);
      end = std::chrono::system_clock::now();
      time_2d_individual += std::chrono::duration<double>(end - start).count();

      double error_1d = std::abs(mi_1d_cell - mi_1d[j]);
      if (error_1d > 0.000001) {
        std::cout <<
          "MI from single computation: " << mi_1d_cell <<
          ", MI from joint computation " << mi_1d[j] << std::endl;
        return -1;
      }
      double error_2d = std::abs(mi_2d_cell - mi_2d[j]);
      if (error_2d > 0.000001) {
        std::cout <<
          "MI from single computation 2D: " << mi_2d_cell <<
          ", MI from joint computation 2D: " << mi_2d[j] << std::endl;
        return -1;
      }
    }
  }

  std::cout << "1D individual time: " << time_1d_individual  << " seconds" << std::endl;
  std::cout << "1D joint time: " << time_1d_joint  << " seconds" << std::endl;
  std::cout << "2D individual time: " << time_2d_individual  << " seconds" << std::endl;
  std::cout << "2D joint time: " << time_2d_joint  << " seconds" << std::endl;

  return 0;
}
