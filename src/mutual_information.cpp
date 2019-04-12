#include <cmath>

#include "wandering_robot/occupancy_state.hpp"
#include "wandering_robot/mutual_information.hpp"

double wandering_robot::MutualInformation::d1(
    const OccupancyState * const states,
    const double * const widths,
    unsigned int num_cells) {

  double mi = 0;
  double unknown_length = 0;
  for (unsigned int i = 0; i < num_cells; i++) {
    if (states[i] == OccupancyState::unknown) {
      // Update the mutual information with a piece-wise integral over the region
      mi +=
        std::exp(-poisson_rate * unknown_length) *
        (1. - std::exp(-poisson_rate * widths[i]));

      // Compute the total length of the unknown cells that the
      // beam passes through before it hits an occupied cell
      unknown_length += widths[i];
    } else if (states[i] == OccupancyState::occupied) {
      break;
    }
  }

  // Use that to compute mutual information
  return mi;
}

void wandering_robot::MutualInformation::d1(
    const OccupancyState * const states,
    const double * const widths,
    unsigned int num_cells,
    double * const mutual_information) {

  // Loop backwards from the end of the sequence,
  // keeping track of the length of unknown cells that
  // the beams pass through before they hit occupied cells
  double mi = 0;
  for (int i = num_cells - 1; i >= 0; i--) {
    if (states[i] == OccupancyState::unknown) {
      mi = 1 + std::exp(-poisson_rate * widths[i]) * (mi - 1);
    } else if (states[i] == OccupancyState::occupied) {
      mi = 0;
    }

    // Compute the mutual information
    mutual_information[i] = mi;
  }
}
