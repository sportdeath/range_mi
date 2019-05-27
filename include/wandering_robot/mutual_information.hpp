#pragma once

#include "wandering_robot/occupancy_state.hpp"

namespace wandering_robot {

/**
 * A class for computing the mutual
 * information between lidar beams and
 * probabilistic occupancy cells.
 */
class MutualInformation {

  public:

    /**
     * Initialize
     */
    MutualInformation() {}
    MutualInformation(double poisson_rate_)
      : poisson_rate(poisson_rate_) {}

    /**
     * Compute the mutual information between the
     * occupancy states and a range measurement taken
     * from cell 0 pointing towards cells 1, 2, etc.
     *
     * The mutual information takes into account the
     * probability that each cell has already been
     * measured by another beam.
     *
     * The events are all considered to be independent
     */
    double d1(
        const OccupancyState * const states,
        const double * const widths,
        const double * const p_not_measured,
        unsigned int num_cells);

    /**
     * The same as above using the FSMI formula.
     */
    double d1_fsmi(
        const OccupancyState * const states,
        unsigned int num_cells,
        double width, // FSMI does not support multiple widths
        double gaussian_width, // In number of cells
        double delta_emp, // The log probability of occupied
        double detla_occ);

    /**
     * Compute the mutual information between the
     * occupancy states and range measurements taken
     * from each cell i pointing towards cells
     * i+1, i+2, etc.
     */
    void d1(
        const OccupancyState * const states,
        const double * const widths,
        const double * const p_not_measured,
        const unsigned int * const line,
        unsigned int num_cells,
        double * const mutual_information);

    /**
     * The same as above but in 2D,
     * i.e. accounting for the fact area
     * spanned by the beams increases linearly
     * with the distance away from the source.
     */
    double d2(
        const OccupancyState * const states,
        const double * const widths,
        const double * const p_not_measured,
        unsigned int num_cells);
    void d2(
        const OccupancyState * const states,
        const double * const widths,
        const double * const p_not_measured,
        const unsigned int * const line,
        unsigned int num_cells,
        double * const mutual_information);

    /**
     * Given a previous value for p_not_measured
     * at some point in space, compute a new value
     * for p_not_measured at that point, given that
     * a measurement has been made in that direction,
     * passing through an unknown region of length
     * unknown_length.
     */
    void condition(
        double & p_not_measured,
        double unknown_length);

  private:

    /**
     * The rate at which obstacle boundaries are expected
     * to occur in unknown space.
     */
    double poisson_rate;
};

}
