#include <ros/ros.h>

#include <nav_msgs/OccupancyGrid.h>

#include "wandering_robot/bresenham.hpp"
#include "wandering_robot/mutual_information.hpp"

std::string map_topic = "/map";
std::string mi_topic = "/mi";
double poisson_rate = 0.05;
unsigned int num_iterations = 10000000;
unsigned int draw_rate = 100000;
double epsilon = 0.1;

class MutualInformationVisualizer {

  public:

    MutualInformationVisualizer() {
      // Construct a publisher for mutual information
      mi_pub = n.advertise<nav_msgs::OccupancyGrid>(mi_topic, 1, true);

      // Subscribe to maps
      map_sub = n.subscribe(map_topic, 1, &MutualInformationVisualizer::map_callback, this);
    }

    void map_callback(const nav_msgs::OccupancyGrid & map_msg) {
      unsigned int height = map_msg.info.height;
      unsigned int width = map_msg.info.width;

      // Convert to states
      std::vector<wandering_robot::OccupancyState> states(height * width);
      for (unsigned int i = 0; i < states.size(); i++) {
        double value = map_msg.data[i]/100.;
        if (value < 0) {
          states[i] = wandering_robot::OccupancyState::unknown;
        } else if (value < 0.5 - epsilon) {
          states[i] = wandering_robot::OccupancyState::free;
        } else if (value > 0.5 + epsilon) {
          states[i] = wandering_robot::OccupancyState::occupied;
        } else {
          states[i] = wandering_robot::OccupancyState::unknown;
        }
      }

      // Initialize the beam sampler and mutual information
      wandering_robot::Bresenham bresenham(height, width);
      wandering_robot::MutualInformation mi_(poisson_rate);

      // Initialize a map and a line
      unsigned int line_length = std::max(height, width);
      std::vector<unsigned int> line(line_length);
      std::vector<double> line_widths(line_length), line_p_not_measured(line_length);
      std::vector<double> line_mi(line_length);
      std::vector<wandering_robot::OccupancyState> line_states(line_length);
      std::vector<double> mi(height * width, 0);

      double x, y, theta;
      for (unsigned int i = 0; i < num_iterations; i++) {
        if (i % draw_rate == 0) {
          draw_map(map_msg, mi);
        }

        // Randomly sample a point
        bresenham.sample(x, y, theta);

        // Compute Bresenham's line
        unsigned int num_cells;
        bresenham.line(
            y, x, theta,
            line.data(),
            num_cells);

        // Compute width
        double sc = std::max(std::abs(std::sin(theta)), std::abs(std::cos(theta)));
        double width = 1./sc;

        for (unsigned int j = 0; j < num_cells; j++) {
          line_states[j] = states[line[j]];
          line_widths[j] = width;
          line_p_not_measured[j] = 1;
        }

        // Make vectors of the states, etc.
        mi_.d2(
            line_states.data(),
            line_widths.data(),
            line_p_not_measured.data(),
            num_cells,
            line_mi.data());

        for (unsigned int j = 0; j < num_cells; j++) {
          mi[line[j]] += line_mi[j];
        }
      }
    }

    void draw_map(const nav_msgs::OccupancyGrid & map_msg, const std::vector<double> & mi) {
      // Convert to int8
      nav_msgs::OccupancyGrid mi_msg;
      mi_msg.data = std::vector<int8_t>(map_msg.info.height * map_msg.info.width);
      double mi_max = *std::max_element(mi.begin(), mi.end());
      for (size_t i = 0; i < mi_msg.data.size(); i++) {
        mi_msg.data[i] = 100 * (1 - mi[i]/mi_max);
      }

      // Add info
      mi_msg.header = map_msg.header;
      mi_msg.header.stamp = ros::Time::now();
      mi_msg.info = map_msg.info;

      // Publish
      mi_pub.publish(mi_msg);
    }

  private:
    ros::NodeHandle n;
    ros::Subscriber map_sub;
    ros::Publisher mi_pub;
};

int main(int argc, char ** argv) {
  ros::init(argc, argv, "mutual_information_visualizer");
  MutualInformationVisualizer mi;
  ros::spin();
  return 0;
}
