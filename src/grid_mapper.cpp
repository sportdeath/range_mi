#include <vector>
#include "wandering_robot/occupancy_state.hpp"

namespace wandering_robot {

class GridMapper {
  public:
    GridMapper(unsigned int height, unsigned int width) {
      states_ = std::vector<wandering_robot::OccupancyState>(height * width);
    }

    void apply_scan(double x, double y, const std::vector<double> & scan) {
      // Draw a line in the direction of each scan
      
      // For each point if the distance is less than
      // the measured distance, make it free.
      // Otherwise make it occupied.
    }

    const std::vector<wandering_robot::OccupancyState> states();

  private:

    std::vector<wandering_robot::OccupancyState> states_;
};

}
