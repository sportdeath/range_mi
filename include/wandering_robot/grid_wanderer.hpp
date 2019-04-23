#pragma once

#include "wandering_robot/mutual_information.hpp"
#include "wandering_robot/grid_line.hpp"

namespace wandering_robot {

class GridWanderer {

  public:

    GridWanderer() {}

    GridWanderer(
        double poisson_ratio,
        bool beam_independence_=true)
      : beam_independence(beam_independence_) {
      mi_computer = wandering_robot::MutualInformation(poisson_ratio);
    }

    void set_map(
        const std::vector<wandering_robot::OccupancyState> & states_,
        unsigned int height,
        unsigned int width);

    void iterate_mi();

    // Getters
    const std::vector<double> & mi() const {return mi_;};
    const std::vector<double> & p_not_measured() const {return p_not_measured_;};

  private:

    // Parameters
    bool beam_independence;

    // Initialize a map and a line
    std::vector<wandering_robot::OccupancyState> states;
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
};

}
