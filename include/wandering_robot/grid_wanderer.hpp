#pragma once

#include "wandering_robot/mutual_information.hpp"
#include "wandering_robot/grid_line.hpp"

namespace wandering_robot {

class GridWanderer {

  public:

    GridWanderer() {}

    GridWanderer(
        unsigned int height,
        unsigned int width,
        double poisson_rate,
        bool beam_independence_=true);

    void set_map(const std::vector<wandering_robot::OccupancyState> & map_) {
      map = map_;
    }

    void make_scan(unsigned int cell, unsigned int num_beams);
    void reset_states() {std::fill(states_.begin(), states_.end(), wandering_robot::OccupancyState::unknown);}

    void reset_mi() {std::fill(mi_.begin(), mi_.end(), 0);}
    void accrue_mi(double spatial_interpolation, double angular_interpolation);

    void condition(unsigned int cell, unsigned int angular_steps);
    void reset_p_not_measured() {std::fill(p_not_measured_.begin(), p_not_measured_.end(), 1);}

    void dijkstra(unsigned int start, std::vector<double> & distances, std::vector<unsigned int> & parents);

    const std::vector<double> & mi() const {return mi_;}
    const std::vector<double> & p_not_measured() const {return p_not_measured_;}
    const std::vector<wandering_robot::OccupancyState> & states() const {return states_;}

  private:

    // Parameters
    bool beam_independence;

    // Initialize a map and a line
    std::vector<wandering_robot::OccupancyState> states_;
    std::vector<wandering_robot::OccupancyState> map;
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
