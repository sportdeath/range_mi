/**
 * This is meant to be a human-readable example
 * of how this library can be used to compute
 * mutual information surfaces.
 */

#include <vector>
#include <cmath>
#include <iostream>
#include <iomanip>

#include <range_mi/grid_mi.hpp>

int main() {

// Create a map
unsigned int height = 7;
unsigned int width = 7;
// The map has a vacant region in the bottom
// center. This region is is bordered on the left
// and right by occupied "walls". Elsewhere there
// is relatively unknown space.
std::vector<double> vacancy = {
  0.8,  0.8,  0.8,  0.8,  0.8,  0.8, 0.8,
  0.8,  0.8,  0.8,  0.8,  0.8,  0.8, 0.8,
  0.8,  0.8,  0.8,  0.8,  0.8,  0.8, 0.8,
  0.8, 0.01, 0.99, 0.99, 0.99, 0.01, 0.8,
  0.8, 0.01, 0.99, 0.99, 0.99, 0.01, 0.8,
  0.8, 0.01, 0.99, 0.99, 0.99, 0.01, 0.8,
  0.8, 0.01, 0.99, 0.99, 0.99, 0.01, 0.8,
};
           
// The number of beams to compute mutual information with.
unsigned int num_beams = 200;

// Initialize the mutual information computation device
range_mi::GridMI mi_computer(height, width);

// Compute the mutual information
mi_computer.compute_mi(vacancy.data(), num_beams);

// Print out the mutual information map
std::cout << "Mutual information surface: " << std::endl;
std::cout << std::fixed;
std::cout << std::setprecision(0);
for (unsigned int i = 0; i < height; i++) {
  for (unsigned int j = 0; j < width; j++) {
    std::cout << std::setw(6) << mi_computer.mi()[j + i * width];
  }
  std::cout << std::endl;
}

}
