#pragma once

#include "wandering_robot/mutual_information.hpp"
#include "wandering_robot/grid_line.hpp"

namespace wandering_robot {

class GridMutualInformation {

  public:

    GridMutualInformation() {}

    GridMutualInformation(
        const std::vector<wandering_robot::OccupancyState> * states,
        unsigned int height,
        unsigned int width,
        double poisson_rate,
        bool beam_independence_=true);

    /**
     * Compute the mutual information between the map
     * and a measurement made from the cell, "cell",
     * emanating "angular_steps" beams.
     */
    double compute_mi(unsigned int cell, unsigned int angular_steps);

    /**
     * Compute the mutual information between the map
     * and measurements made from each cell on the map.
     */
    void compute_mi_surface(unsigned int spatial_jitter, unsigned int num_beams);
    const std::vector<double> & mi_surface() const {return mi_;}
    void reset_mi_surface() {std::fill(mi_.begin(), mi_.end(), 0);}
    void compute_mi_surface_beam(
        double & spatial_interpolation, double & angular_interpolation,
        unsigned int spatial_jitter, unsigned int num_beams);

    void reset_p_not_measured() {std::fill(p_not_measured_.begin(), p_not_measured_.end(), 1);}
    void condition(unsigned int cell, unsigned int angular_steps);

    const std::vector<double> & p_not_measured() const {return p_not_measured_;}

  private:

    // Parameters
    bool beam_independence;

    // The mutual information and other maps
    const std::vector<wandering_robot::OccupancyState> * states;
    std::vector<double> mi_;
    std::vector<double> p_not_measured_;

    // Computation devices
    wandering_robot::GridLine grid_line;
    wandering_robot::MutualInformation mi_computer;

    // Reused storage
    double x, y, theta;
    unsigned int num_cells;
    std::vector<unsigned int> line;
    std::vector<double> widths;
    std::vector<double> condition_distances;
};

}
