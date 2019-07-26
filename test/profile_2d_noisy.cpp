#include <vector>
#include <random>
#include <chrono>
#include <iostream>

#include <range_entropy/expected_noisy.hpp>

using namespace range_entropy;

// Define constants
double integration_step = 0.01;
double noise_dev = 2;
double noise_half_width = 3 * noise_dev;
unsigned int num_iterations = 1;
unsigned int num_cells = 100;

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
  std::vector<double> width(num_cells);
  std::vector<double> expected_information2(num_cells, 0);
  std::vector<double> hit_pdf(num_cells/integration_step);
  std::vector<double> miss_pdf(num_cells/integration_step);
  // Make a line
  std::vector<unsigned int> line(num_cells);
  std::iota(line.begin(), line.end(), 0);

  std::chrono::duration<double> elapsed(0);

  for (unsigned int i = 0; i < num_iterations; i++) {
    // Randomize then
    random_p(p_free);
    random_p(width);

    auto start = std::chrono::high_resolution_clock::now();

    // Compute the entropy
    expected_noisy::line(
        line.data(), p_free.data(), width.data(), num_cells,
        noise_dev, noise_half_width, integration_step,
        true, 2,
        hit_pdf.data(), miss_pdf.data(), expected_information2.data());

    auto end = std::chrono::high_resolution_clock::now();
    elapsed += end - start;
  }

  std::cout << "Computed " << num_iterations << " beams " <<
    " with " << num_cells << " cells in " <<
    elapsed.count() << " seconds." << std::endl;
}
