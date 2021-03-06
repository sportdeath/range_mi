/**
 * This test computes the mutual information
 * of a large number of random beams to determine
 * the average time spent computing mutual information
 * per cell
 */

#include <vector>
#include <random>
#include <chrono>
#include <iostream>

#include <range_mi/barely_distorted.hpp>

#include "helpers.hpp"

// Define constants
unsigned int num_iterations = 100000;
unsigned int num_cells = 100;
double dtheta = 0.1;
const unsigned int dimension = 1;

int main() {
  // Initialize the inputs and outputs
  std::vector<double> vacancy(num_cells);
  std::vector<double> p_not_measured(num_cells);
  std::vector<double> width(num_cells);
  std::vector<double> mi(num_cells, 0);

  // Make a line
  std::vector<unsigned int> line(num_cells);
  std::iota(line.begin(), line.end(), 0);

  // Set the clock to zero
  std::chrono::duration<double> elapsed(0);

  for (unsigned int i = 0; i < num_iterations; i++) {

    // Randomize the inputs
    random_p(vacancy);
    random_p(p_not_measured);
    random_p(width);

    auto start = std::chrono::high_resolution_clock::now();

    // Compute the entropy
    range_mi::barely_distorted::line<dimension>(
        line.data(),
        vacancy.data(),
        p_not_measured.data(),
        width.data(),
        num_cells,
        dtheta,
        mi.data());

    auto end = std::chrono::high_resolution_clock::now();
    elapsed += end - start;
  }

  std::cout << "Computed " << num_iterations << " beams " <<
    " with " << num_cells << " cells in " <<
    elapsed.count() << " seconds." << std::endl;

  std::cout << "That makes an average time of " <<
    (elapsed.count() * std::nano::den)/(num_iterations * num_cells) <<
    " nanoseconds per cell" << std::endl;
}
