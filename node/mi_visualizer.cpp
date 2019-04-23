#include <ros/ros.h>

#include <nav_msgs/OccupancyGrid.h>
#include <geometry_msgs/PointStamped.h>

#include "wandering_robot/grid_wanderer.hpp"

class MutualInformationVisualizer {

  public:

    MutualInformationVisualizer() {
      // Initialize the node handle
      n = ros::NodeHandle("~");

      // Fetch the ROS parameters
      std::string mi_topic, mi_max_topic, map_topic, click_topic;
      double poisson_rate;
      bool beam_independence;
      n.getParam("mi_topic", mi_topic);
      n.getParam("mi_max_topic", mi_max_topic);
      n.getParam("map_topic", map_topic);
      n.getParam("click_topic", click_topic);
      n.getParam("num_beams", num_beams);
      n.getParam("beams_per_draw", beams_per_draw);
      n.getParam("unknown_threshold", unknown_threshold);
      n.getParam("poisson_rate", poisson_rate);
      n.getParam("beam_independence", beam_independence);

      // Initialize the wanderer
      w = wandering_robot::GridWanderer(poisson_rate, beam_independence);

      // Construct a publisher for mutual information
      mi_pub = n.advertise<nav_msgs::OccupancyGrid>(mi_topic, 1, true);
      mi_max_pub = n.advertise<geometry_msgs::PointStamped>(mi_max_topic, 1);

      // Subscribe to maps and clicked points
      map_sub = n.subscribe(map_topic, 1, &MutualInformationVisualizer::map_callback, this);
      click_sub = n.subscribe(click_topic, 1, &MutualInformationVisualizer::click_callback, this);
    }

    void map_callback(const nav_msgs::OccupancyGrid & map_msg) {
      // Store map information
      map_info = map_msg.info;
      map_frame_id = map_msg.header.frame_id;

      // Convert to states
      std::vector<wandering_robot::OccupancyState> states(map_info.height * map_info.width);
      for (unsigned int i = 0; i < states.size(); i++) {
        double value = map_msg.data[i]/100.;
        if (value < 0) {
          states[i] = wandering_robot::OccupancyState::unknown;
        } else if (value < 0.5 - unknown_threshold) {
          states[i] = wandering_robot::OccupancyState::free;
        } else if (value > 0.5 + unknown_threshold) {
          states[i] = wandering_robot::OccupancyState::occupied;
        } else {
          states[i] = wandering_robot::OccupancyState::unknown;
        }
      }

      // Initialize the map for computation
      w.set_map(states, map_info.height, map_info.width);

      // Compute the MI, drawing the map regularly
      for (int i = 0; i < num_beams; i++) {
        if (i % beams_per_draw == 0) {
          draw_map();
          if (not ros::ok()) break;
        }
        w.iterate_mi();
      }
      draw_map();
    }

    void click_callback(const geometry_msgs::PointStamped & click_msg) {
      std::cout << "Click!" << std::endl;
    }

    void draw_map() {
      // Convert to int8
      nav_msgs::OccupancyGrid mi_msg;
      mi_msg.data = std::vector<int8_t>(map_info.height * map_info.width);
      auto mi_max = std::max_element(w.mi().begin(), w.mi().end());
      for (size_t i = 0; i < mi_msg.data.size(); i++) {
        mi_msg.data[i] = 100 * (1 - w.mi()[i]/(*mi_max));
      }

      // Add info
      mi_msg.header.frame_id = map_frame_id;
      mi_msg.header.stamp = ros::Time::now();
      mi_msg.info = map_info;

      // Publish
      mi_pub.publish(mi_msg);

      // Plot the position of the maximum mutual information
      geometry_msgs::PointStamped mi_max_msg;
      mi_max_msg.header = mi_msg.header;
      unsigned int mi_max_cell = std::distance(w.mi().begin(), mi_max);
      unsigned int y = mi_max_cell/map_info.width;
      unsigned int x = mi_max_cell - y * map_info.width;
      mi_max_msg.point.x = x * map_info.resolution;
      mi_max_msg.point.y = y * map_info.resolution;
      mi_max_pub.publish(mi_max_msg);
    }

  private:
    ros::NodeHandle n;
    ros::Subscriber map_sub, click_sub;
    ros::Publisher mi_pub, mi_max_pub;

    // Parameters
    double unknown_threshold;
    int num_beams, beams_per_draw;

    // Map data
    nav_msgs::MapMetaData map_info;
    std::string map_frame_id;

    // Computation devices
    wandering_robot::GridWanderer w;

};

int main(int argc, char ** argv) {
  ros::init(argc, argv, "mutual_information_visualizer");
  MutualInformationVisualizer mi;
  ros::spin();
  return 0;
}
