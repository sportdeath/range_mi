#pragma once

#include "wandering_robot/mutual_information.hpp"
#include "wandering_robot/grid_line.hpp"
#include "wandering_robot/grid_mapper.hpp"

namespace wandering_robot {

class GridMapper {

  public:

    GridMapper() {}

    GridMapper(
        unsigned int height,
        unsigned int width);

    void set_ground_truth(const std::vector<wandering_robot::OccupancyState> & map_gt_) {
      map_gt = map_gt_;
    }

    void make_scan(unsigned int cell, unsigned int num_beams);
    void reset_map() {std::fill(map_.begin(), map_.end(), wandering_robot::OccupancyState::unknown);}

    void dijkstra(unsigned int start, std::vector<double> & distances, std::vector<int> & parents) const;
    std::vector<unsigned int> dijkstra_path(int end, const std::vector<int> & parents) const;

    const std::vector<wandering_robot::OccupancyState> & map() const {return map_;}

  private:

    // Initialize a map and a line
    std::vector<wandering_robot::OccupancyState> map_;
    std::vector<wandering_robot::OccupancyState> map_gt;

    // Computation devices
    wandering_robot::GridLine grid_line;

    // Reused storage
    double x, y, theta;
    unsigned int num_cells;
    std::vector<unsigned int> line;
    std::vector<double> widths;
};

}
