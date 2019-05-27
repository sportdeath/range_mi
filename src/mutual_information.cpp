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
    const unsigned int * const line,
    unsigned int num_cells,
    double * const mutual_information) {

  // Loop backwards from the end of the sequence,
  double mutual_information_local = 0;
  double p_no_hit;
  for (int j = num_cells - 1; j >= 0; j--) {
    unsigned int i = line[j];

    if (states[i] == OccupancyState::unknown) {
      p_no_hit = std::exp(-poisson_rate * widths[j]);

      // Multiply the contribution from the current cell,
      // (1 - e^(-lambda w)), with the probability that
      // a cell has not already been hit. Then add this
      // to the previous mutual information, scaled by the
      // probability that the beam reaches that far.
      mutual_information_local =
        p_not_measured[i] * (1 - p_no_hit) + p_no_hit * mutual_information_local;
    } else if (states[i] == OccupancyState::occupied) {
      mutual_information_local = 0;
    }

    mutual_information[i] += mutual_information_local;
  }
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
    const unsigned int * const line,
    unsigned int num_cells,
    double * const mutual_information) {

  // Loop backwards from the end of the sequence,
  double a = 0, b = 0, c = 0;
  double p_no_hit;

  for (int j = num_cells - 1; j >= 0; j--) {
    unsigned int i = line[j];
    if (states[i] == OccupancyState::unknown) {
      p_no_hit = std::exp(-poisson_rate * widths[j]);
      a = p_not_measured[i] * (1 - p_no_hit * (1 + poisson_rate * widths[j])) + p_no_hit * a;
      b = p_no_hit * (b + widths[j] * c);
      c = p_not_measured[i] * poisson_rate * (1 - p_no_hit) + p_no_hit * c;
    } else if (states[i] == OccupancyState::free) {
      b += widths[j] * c;
    } else {
      a = 0; b = 0; c = 0;
    }

    // Compute the mutual information
    mutual_information[i] += a + b;
    }
}

void wandering_robot::MutualInformation::condition(
    double & p_not_measured,
    double unknown_length) {
  p_not_measured *= (1 - std::exp(-poisson_rate * unknown_length));
}

double wandering_robot::MutualInformation::d1_fsmi(
    const OccupancyState * const states,
    unsigned int num_cells,
    double width,
    double standard_deviation,
    unsigned int gaussian_width,
    double delta_emp,
    double detla_occ) {

  // Initialize the total mutual information to zero
  double mi = 0;

  // The probability of hitting a cell
  // decreases over time
  double probability_of_hitting_cell = 1;

  // The MI contribution
  double mi_contribution = 0;
  std::vector<double> mi_contributions(num_cells);
  for (unsigned int i = 0; i < num_cells; i++) {
    mi_contributions[i] = mi_contribution + f(delta_occ, state[i]);
    mi_contribution += f(delta_emp, state[i]) ;

    // Update the mutual information contribtion
    double mi_contribution;
    fsmi_contribution(delta, state) {
      if (state == free) {
        // r = 0
        // f = log(delta) (1 - 1/(delta + 1))
      } else if (state == occupied) {
        // r == infinity
        // f = 0
      } else if (state == unknown) {
        // r = 1? or exp(-poisson_rate)/(1 - exp(-poisson_rate)
      }
    }
  }

  // Loop over the cells
  for (unsigned int i = 0; i < num_cells; i++) {

    // Loop over the Gaussian width
    for (unsigned int j = -gaussian_width; j <= gaussian_width; j++) {
      // Make sure we don't go over the boundary
      if (i + j < 0 or i + j >= num_cells) continue;

      // Compute the integral of the Gaussian centered
      // at (i + 0.5), across cells i + j to i + j + 1
      // 
      // As the std deviation -> 0 this becomes 1
      // for cell i and 0 elsewhere.
      double mean = width * (i + 0.5);
      double left = width * (i + j);
      double right = width * (i + j + 1);
      double gaussian_weight = 0.5 * (
          std::erf((right - mean)/(std::sqrt(2) * standard_deviation)) -
          std::erf((left - mean)/(std::sqrt(2) * standard_deviation)));

      // Accumulate the mutual information
      mi +=
        probability_of_hitting_cell *
        gaussian_weight *
        mi_contributions[i + j];
    }

    // Update the probability of hitting the next cell
    if (states[i] == wandering_robot::OccupancyState::occupied) {
      // Go to zero!
      probability_of_hitting_cell = 0;
    } else if (states[i] == wandering_robot::OccupancyState::unknown) {
      // The probability decreases according to the Poisson distribution
      // If poisson_rate = log(2) then this is equal to
      // 1 - 0.5 ^ width
      probability_of_hitting_cell *=
        (1 - std::exp(-poisson_rate * width));
    }

  }

  return mi;
}
