#include <cmath>
#include <vector>
#include <iostream>
#include <random>

#include <range_entropy/expected.hpp>

/**
 * This script compares the output of FSMI
 * with the output of the expected range entropy.
 *
 * To do this we chop the map into tiny cells
 * with constant width so that FSMI can be applied.
 * Then we run the algorithm with zero noise.
 */

double cell_size = 0.00001;
unsigned int num_cells = 100;

void fsmi_iteration(
    double p_free,
    double & information,
    double & miss_p,
    double & miss_info) {

  // The probability of hitting the cell
  // is the probability that everything
  // prior was missed times the probability
  // of hitting the current cell.
  double hit_p = miss_p * (1 - p_free);

  // The info from hitting is the info
  // from missing everything plus the info
  // from hitting
  // lim delta -> infinity f(r, delta)
  double hit_info = miss_info - std::log(1 - p_free);

  // The information is the probability
  // of hitting the cell times information
  // gained by hitting
  information += hit_p * hit_info;
  
  // Update the miss probability
  miss_p = miss_p * p_free;
  // Update the miss info
  // lim delta -> 0 f(r, delta)
  miss_info += -std::log(p_free);
}

// Initialize random generator
std::random_device random_device;
std::mt19937 gen(random_device());
std::uniform_real_distribution<double> dist(0.,1.);

// ... and a way to make random vectors
void random_p(std::vector<double> & p) {
  for (size_t i = 0; i < p.size(); i++) {
    p[i] = dist(gen);
  }
}

int main() {
  // Initialize the inputs and outputs
  std::vector<double> p_free(num_cells);
  std::vector<double> p_not_measured(num_cells,1);
  std::vector<double> width(num_cells);
  // Randomize then
  random_p(p_free);
  random_p(width);

  // Make a line
  std::vector<unsigned int> line(num_cells);
  std::iota(line.begin(), line.end(), 0);

  // Compute FSMI
  double fsmi = 0;
  double miss_p = 1;
  double miss_info = 0;
  double width_sum = 0;
  double r = 0;
  unsigned int i = 0;
  while (i < num_cells) {
    // The probability that the cell is free
    // is p^cell_size
    double p_free_cell = std::pow(p_free[i], cell_size);

    // Accumulate FSMI
    fsmi_iteration(
        p_free_cell,
        fsmi,
        miss_p,
        miss_info);

    // Make a step
    r += cell_size;
    if (r - width_sum > width[i]) {
      // Cell completed, move on
      width_sum += width[i];
      i++;
    }
  }

  // Because FSMI returns mutual information,
  // not differential entropy, we need to subtract
  // away the noise added by the discretization.
  //
  // The entropy of a uniform distribution of
  // width cell_size is -log(cell_size)
  double fsmi_entropy = fsmi - (-std::log(cell_size));

  std::vector<double> expected_information1(num_cells, 0);
  range_entropy::expected::information1(
      line.data(), p_free.data(), p_not_measured.data(), width.data(), num_cells,
      expected_information1.data());

  std::cout << " FSMI entropy: " << fsmi_entropy << std::endl;
  std::cout << "Range entropy: " << expected_information1[0] << std::endl;
}
