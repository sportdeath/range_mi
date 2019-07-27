#include <vector>
#include <random>
#include <chrono>
#include <iostream>

#include <range_entropy/expected.hpp>

using namespace range_entropy;

// Define constants
unsigned int num_iterations = 100000;
unsigned int num_cells = 1000;

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
  std::vector<double> p_not_measured(num_cells);
  std::vector<double> width(num_cells);
  std::vector<double> expected_information2(num_cells, 0);
  // Make a line
  std::vector<unsigned int> line(num_cells);
  std::iota(line.begin(), line.end(), 0);

  std::chrono::duration<double> elapsed(0);

  for (unsigned int i = 0; i < num_iterations; i++) {

    // Randomize then
    random_p(p_free);
    random_p(p_not_measured);
    random_p(width);

    auto start = std::chrono::high_resolution_clock::now();

    // Compute the entropy
    expected::line(
        line.data(), p_free.data(), p_not_measured.data(), width.data(), num_cells,
        true, 2,
        expected_information2.data());

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
