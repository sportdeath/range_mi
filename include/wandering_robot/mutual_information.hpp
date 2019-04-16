#pragma once

#include "wandering_robot/occupancy_state.hpp"

namespace wandering_robot {

class MutualInformation {

  public:

    /**
     * Initialize
     */
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
     * Compute the mutual information between the
     * occupancy states and range measurements taken
     * from each cell i pointing towards cells
     * i+1, i+2, etc.
     */
    void d1(
        const OccupancyState * const states,
        const double * const widths,
        const double * const p_not_measured,
        unsigned int num_cells,
        double * const mutual_information);

    double d2(
        const OccupancyState * const states,
        const double * const widths,
        const double * const p_not_measured,
        unsigned int num_cells);

    void d2(
        const OccupancyState * const states,
        const double * const widths,
        const double * const p_not_measured,
        unsigned int num_cells,
        double * const mutual_information);

    /**
     * Compute mutual information for a grid
     */
    void d2_grid(
        const OccupancyState * const states,
        const double * const p_not_measured,
        const unsigned int * const line,
        double theta,
        unsigned int num_cells,
        double * const mutual_information);

    /**
     * Update step for mutual information computation
     */
    void d2_update(
        double & a, double & b, double & c,
        double state,
        double width,
        double p_not_measured,
        double p_no_hit,
        double & mutual_information);

  private:

    /**
     * The rate at which obstacle boundaries are expected
     * to occur in unknown space.
     */
    double poisson_rate;
};

}
