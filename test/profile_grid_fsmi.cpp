/**
 * Compute the time to compute the mutual information
 * in a grid at a variety of map sizes.
 */

#include <vector>
#include <chrono>
#include <fstream>
#include <iostream>
#include <range_mi/grid_line.hpp>

#include "helpers.hpp"

double func(double delta, double vacancy) {
  vacancy = std::max(vacancy, 0.00000001);
  vacancy = std::min(vacancy, 0.99999999);
  double r = (1 - vacancy)/vacancy;
  return std::log((r + 1)/(r + 1/delta)) - std::log(delta)/(r * delta + 1);
}

int main() {

unsigned int num_beams = 200;
unsigned int max_side_length = 1000;
unsigned int max_area = max_side_length * max_side_length;

// FSMI parameters
double delta_occ = 1.5;
double delta_emp = 1/delta_occ;

// Randomize the occupancy grid
std::vector<double> vacancy(max_area);
random_p(vacancy);

// Initialize a place for mutual information
std::vector<double> mi(max_area);

// Initialize data
std::vector<unsigned int> line(2*max_side_length);
std::vector<double> widths(2*max_side_length);

// Output file
std::ofstream f;
f.open("profile_grid_fsmi_results.txt");

for (unsigned int side_length = 1; side_length < max_side_length; side_length++) {
  std::cout << side_length << std::endl;

  double theta;
  unsigned int num_cells;
  double p_previous_empty;
  double p_i_first_non_empty;
  double info_gain, info_loss;

  auto start = std::chrono::high_resolution_clock::now();
  for (unsigned int x = 0; x < side_length; x++) {
    for (unsigned int y = 0; y < side_length; y++) {
      for (unsigned int b = 0; b < num_beams; b++) {
        theta = b/num_beams;

        // Compute the intersections of
        // the line with the grid
        range_mi::grid_line::draw(
            side_length, side_length,
            x, y, theta,
            line.data(),
            widths.data(),
            num_cells);

        // Compute the mutual information
        p_previous_empty = 1;
        info_loss = 0;
        mi[x + y * side_length] = 0;
        for (unsigned int i = 0; i < num_cells; i++) {
          p_i_first_non_empty = p_previous_empty * (1 - vacancy[line[i]]);
          p_previous_empty *= vacancy[line[i]]; 

          info_gain = info_loss + func(delta_occ, vacancy[line[i]]);
          info_loss += func(delta_emp, vacancy[line[i]]);

          // Assume the noise width is zero to be fair,
          // so no inner loop
          mi[x + y * side_length] += p_i_first_non_empty * info_gain;
        }
      }
    }
  }
  auto end = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double> total = end - start;
  
  f << side_length * side_length << " " << total.count() << std::endl;
}

f.close();
}
