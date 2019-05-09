#include <vector>
#include <iostream>
#include <chrono>

#include <wandering_robot/occupancy_state.hpp>
#include <wandering_robot/mutual_information.hpp>

#include "random_inputs.hpp"

using namespace wandering_robot;

// Initialize constants
double poisson_rate = 0.1234;
int num_iterations = 100;
size_t num_cells = 1000;

int main() {
  // Initialize the inputs and outputs
  std::vector<OccupancyState> states(num_cells);
  std::vector<double> widths(num_cells);
  std::vector<double> p_not_measured(num_cells);
  std::vector<double> mi_1d(num_cells), mi_2d(num_cells);
  std::vector<unsigned int> line(num_cells);

  // Make the line 0...num_cells -1
  std::iota(line.begin(), line.end(), 0);

  // Initialize timing
  std::chrono::time_point<std::chrono::system_clock> start, end;
  double time_1d_joint = 0, time_1d_individual = 0;
  double time_2d_joint = 0, time_2d_individual = 0;

  // Initialize the mutual information
  MutualInformation mi(poisson_rate);

  double max_error_1d = 0;
  double max_error_2d = 0;
  for (int i = 0; i < num_iterations; i++) {
    // Generate random cells
    random_states(states);
    random_probabilities(widths);
    random_probabilities(p_not_measured);

    // Clear MI
    std::fill(mi_1d.begin(), mi_1d.end(), 0);
    std::fill(mi_2d.begin(), mi_2d.end(), 0);

    // Compute the mutual information all at once

    // In 1 dimension ...
    start = std::chrono::system_clock::now();
    mi.d1(states.data(), widths.data(), p_not_measured.data(), line.data(), num_cells, mi_1d.data());
    end = std::chrono::system_clock::now();
    time_1d_joint += std::chrono::duration<double>(end - start).count();

    // ... and in 2 dimensions.
    start = std::chrono::system_clock::now();
    mi.d2(states.data(), widths.data(), p_not_measured.data(), line.data(), num_cells, mi_2d.data());
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
      if (error_1d > max_error_1d) max_error_1d = error_1d;
      double error_2d = std::abs(mi_2d_cell - mi_2d[j]);
      if (error_2d > max_error_2d) max_error_2d = error_2d;
    }
  }

  std::cout << "1D individual time: " << time_1d_individual  << " seconds" << std::endl;
  std::cout << "1D joint time: " << time_1d_joint  << " seconds" << std::endl;
  std::cout << "Max error between 1d individual and joint computation: " << max_error_1d << std::endl;
  std::cout << std::endl;
  std::cout << "2D individual time: " << time_2d_individual  << " seconds" << std::endl;
  std::cout << "2D joint time: " << time_2d_joint  << " seconds" << std::endl;
  std::cout << "Max error between 2d individual and joint computation: " << max_error_2d << std::endl;

  return 0;
}
