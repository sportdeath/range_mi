#include <cmath>
#include <algorithm>

#include "wandering_robot/occupancy_state.hpp"
#include "wandering_robot/mutual_information.hpp"

double wandering_robot::MutualInformation::d1(
    const OccupancyState * const states,
    const double * const widths,
    const double * const p_not_measured,
    unsigned int num_cells) {

  double mi = 0;
  double unknown_length = 0;
  for (unsigned int i = 0; i < num_cells; i++) {
    if (states[i] == OccupancyState::unknown) {
      // Update the mutual information with a piece-wise integral over the region
      mi +=
        p_not_measured[i] *
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
    const double * const p_not_measured,
    unsigned int num_cells,
    double * const mutual_information) {

  // Loop backwards from the end of the sequence,
  double mutual_information_local = 0;
  for (int i = num_cells - 1; i >= 0; i--) {
    d1_update(
        states[i],
        widths[i],
        p_not_measured[i],
        -1,
        mutual_information_local,
        mutual_information[i]);
  }
}

void wandering_robot::MutualInformation::d1_grid(
    const OccupancyState * const states,
    const double * const p_not_measured,
    const unsigned int * const line,
    double theta,
    unsigned int num_cells,
    double * const mutual_information) {

  double width = theta_to_cell_width(theta);

  // Compute the probability that there is not a hit
  // within the cell, i.
  double p_no_hit = std::exp(-poisson_rate * width);

  // Loop backwards from the end of the sequence,
  double mutual_information_local = 0;
  for (int j = num_cells - 1; j >= 0; j--) {
    unsigned int i = line[j];
    d1_update(
        states[i],
        width,
        p_not_measured[i],
        p_no_hit,
        mutual_information_local,
        mutual_information[i]);
  }
}

void wandering_robot::MutualInformation::d1_update(
    OccupancyState state,
    double width,
    double p_not_measured,
    double p_no_hit,
    double & mutual_information_local,
    double & mutual_information) {

  if (state == OccupancyState::unknown) {
    if (p_no_hit < 0)
      p_no_hit = std::exp(-poisson_rate * width);

    // Multiply the contribution from the current cell,
    // (1 - e^(-lambda w)), with the probability that
    // a cell has not already been hit. Then add this
    // to the previous mutual information, scaled by the
    // probability that the beam reaches that far.
    mutual_information_local =
      p_not_measured * (1 - p_no_hit) + p_no_hit * mutual_information_local;
  } else if (state == OccupancyState::occupied) {
    mutual_information_local = 0;
  }

  mutual_information += mutual_information_local;
}


double wandering_robot::MutualInformation::d2(
    const OccupancyState * const states,
    const double * const widths,
    const double * const p_not_measured,
    unsigned int num_cells) {

  double mi = 0;
  double total_length = 0;
  double unknown_length = 0;
  for (unsigned int i = 0; i < num_cells; i++) {

    if (states[i] == OccupancyState::unknown) {

      double e_width = std::exp(-poisson_rate * widths[i]);
      double e_unknown = std::exp(-poisson_rate * unknown_length);

      // Update the mutual information with a
      // piece-wise integral over the region
      mi +=
        p_not_measured[i] *
        e_unknown * (
          1 - e_width * (
            1 + poisson_rate * widths[i]
          )
          + total_length * poisson_rate * (
            1 - e_width
          )
        );

      // Increase the length of unknown regions the beam
      // has traveled through
      unknown_length += widths[i];
    } else if (states[i] == OccupancyState::occupied) {
      break;
    }

    // Increase the total length that the beam has traveled
    total_length += widths[i];
  }

  // Use that to compute mutual information
  return mi;
}

void wandering_robot::MutualInformation::d2(
    const OccupancyState * const states,
    const double * const widths,
    const double * const p_not_measured,
    unsigned int num_cells,
    double * const mutual_information) {

  // Loop backwards from the end of the sequence,
  double a = 0, b = 0, c = 0;
  for (int i = num_cells - 1; i >= 0; i--) {
    d2_update(a, b, c, states[i], widths[i], p_not_measured[i], -1, mutual_information[i]);
  }
}

void wandering_robot::MutualInformation::d2_grid(
    const OccupancyState * const states,
    const double * const p_not_measured,
    const unsigned int * const line,
    double theta,
    unsigned int num_cells,
    double * const mutual_information) {

  double width = theta_to_cell_width(theta);

  // Compute the hit probability
  double p_no_hit = std::exp(-poisson_rate * width);

  // Loop backwards from the end of the sequence,
  double a = 0, b = 0, c = 0;
  for (int j = num_cells - 1; j >= 0; j--) {
    unsigned int i = line[j];
    d2_update(a, b, c, states[i], width, p_not_measured[i], p_no_hit, mutual_information[i]);
  }
}

void wandering_robot::MutualInformation::d2_update(
    double & a, double & b, double & c,
    OccupancyState state,
    double width,
    double p_not_measured,
    double p_no_hit,
    double & mutual_information) {

  if (state == OccupancyState::unknown) {
    if (p_no_hit < 0)
      p_no_hit = std::exp(-poisson_rate * width);
    a = p_not_measured * (1 - p_no_hit * (1 + poisson_rate * width)) + p_no_hit * a;
    b = p_no_hit * (b + width * c);
    c = p_not_measured * poisson_rate * (1 - p_no_hit) + p_no_hit * c;
  } else if (state == OccupancyState::free) {
    b += width * c;
  } else {
    a = 0; b = 0; c = 0;
  }

  // Compute the mutual information
  mutual_information += a + b;
}

double wandering_robot::MutualInformation::theta_to_cell_width(
    double theta) {
  return 1./std::max(std::abs(std::sin(theta)), std::abs(std::cos(theta)));
}
