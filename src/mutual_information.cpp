#include <cmath>

#include "wandering_robot/occupancy_state.hpp"
#include "wandering_robot/mutual_information.hpp"

double wandering_robot::MutualInformation::d1(
    const OccupancyState * const states,
    const double * const widths,
    unsigned int num_cells) {

  // Compute the total length of the unknown cells that the
  // beam passes through before it hits an occupied cell
  double unknown_length = 0;
  for (unsigned int i = 0; i < num_cells; i++) {
    if (states[i] == OccupancyState::unknown) {
      unknown_length += widths[i];
    } else if (states[i] == OccupancyState::occupied) {
      break;
    }
  }

  // Use that to compute mutual information
  return 1 - std::exp(-poisson_rate * unknown_length);
}

void wandering_robot::MutualInformation::d1(
    const OccupancyState * const states,
    const double * const widths,
    unsigned int num_cells,
    double * const mutual_information) {

  // Loop backwards from the end of the sequence,
  // keeping track of the length of unknown cells that
  // the beams pass through before they hit occupied cells
  double unknown_length = 0;
  for (int i = num_cells - 1; i >= 0; i--) {
    if (states[i] == OccupancyState::unknown) {
      unknown_length += widths[i];
    } else if (states[i] == OccupancyState::occupied) {
      unknown_length = 0;
    }

    // Compute the mutual information
    mutual_information[i] = 1 - std::exp(-poisson_rate * unknown_length);
  }
}
