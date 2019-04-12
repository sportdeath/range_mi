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
     */
    double d1(
        const OccupancyState * const states,
        const double * const widths,
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
        unsigned int num_cells,
        double * const mutual_information);

  private:

    /**
     * The rate at which obstacle boundaries are expected
     * to occur in unknown space.
     */
    double poisson_rate;
};

}
