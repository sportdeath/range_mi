#include <ros/ros.h>

#include <nav_msgs/OccupancyGrid.h>
#include <geometry_msgs/PointStamped.h>

#include "wandering_robot/grid_mutual_information.hpp"

class MutualInformationVisualizer {

  public:
    MutualInformationVisualizer() {
      // Initialize the node handle
      n = ros::NodeHandle("~");

      // Fetch the ROS parameters
      std::string mi_topic, p_not_measured_topic, click_topic, map_topic;
      n.getParam("mi_topic", mi_topic);
      n.getParam("p_not_measured_topic", p_not_measured_topic);
      n.getParam("click_topic", click_topic);
      n.getParam("map_topic", map_topic);
      n.getParam("spatial_jitter", spatial_jitter);
      n.getParam("num_beams", num_beams);
      n.getParam("condition_steps", condition_steps);
      n.getParam("unknown_threshold", unknown_threshold);
      n.getParam("poisson_rate", poisson_rate);
      n.getParam("beam_independence", beam_independence);

      // Construct a publisher for mutual information
      mi_pub = n.advertise<nav_msgs::OccupancyGrid>(mi_topic, 1, true);
      p_not_measured_pub = n.advertise<nav_msgs::OccupancyGrid>(p_not_measured_topic, 1, true);

      // Subscribe to maps and clicked points
      map_sub = n.subscribe(map_topic, 1, &MutualInformationVisualizer::map_callback, this);
      click_sub = n.subscribe(click_topic, 1, &MutualInformationVisualizer::click_callback, this);
    }

    void map_callback(const nav_msgs::OccupancyGrid & map_msg) {
      // Store map information
      map_info = map_msg.info;
      map_frame_id = map_msg.header.frame_id;

      // Convert to occupancy states
      map = std::vector<wandering_robot::OccupancyState>(map_info.height * map_info.width);
      for (unsigned int i = 0; i < map.size(); i++) {
        double value = map_msg.data[i]/100.;
        if (value < 0 or value > 1) {
          map[i] = wandering_robot::OccupancyState::unknown;
        } else if (value < 0.5 - unknown_threshold) {
          map[i] = wandering_robot::OccupancyState::free;
        } else if (value > 0.5 + unknown_threshold) {
          map[i] = wandering_robot::OccupancyState::occupied;
        } else {
          map[i] = wandering_robot::OccupancyState::unknown;
        }
      }
      
      // Initialize mutual information computation on the grid
      mi = wandering_robot::GridMutualInformation(
          &map,
          map_info.height,
          map_info.width,
          poisson_rate,
          beam_independence);

      compute_mi(true);
    }

    void compute_mi(bool first) {
      mi.reset_mi_surface();

      double spatial_interpolation = 0;
      double angular_interpolation = 0;
      while (angular_interpolation < 1) {
        // Add the beam
        mi.compute_mi_surface_beam(
            spatial_interpolation,
            angular_interpolation,
            spatial_jitter,
            num_beams);

        // Draw every time a spatial section is completed
        if (spatial_interpolation == 0) {
          if (first) {
            mi_max = *std::max_element(mi.mi_surface().begin(), mi.mi_surface().end());
          }
          draw_map();
          mi.reset_mi_surface();
          if (not ros::ok()) break;
        }
      }
    }

    void click_callback(const geometry_msgs::PointStamped & click_msg) {
      // Convert to map cell
      unsigned int cell = ((int) click_msg.point.y) * map_info.width + ((int) click_msg.point.x);

      // Condition the map on the clicked point
      mi.condition(cell, condition_steps);

      // Update the mutual information
      compute_mi(false);
    }

    void draw_map() {
      // Construct a message for the mutual information surface
      nav_msgs::OccupancyGrid mi_msg;
      mi_msg.data = std::vector<int8_t>(map_info.height * map_info.width);
      for (size_t i = 0; i < mi_msg.data.size(); i++) {
        mi_msg.data[i] = 100 * (1 - mi.mi_surface()[i]/(mi_max));
      }

      // Add info
      mi_msg.header.frame_id = map_frame_id;
      mi_msg.header.stamp = ros::Time::now();
      mi_msg.info = map_info;

      // Publish
      mi_pub.publish(mi_msg);

      // Do the same for p_not measured
      nav_msgs::OccupancyGrid p_not_measured_msg = mi_msg;
      for (size_t i = 0; i < p_not_measured_msg.data.size(); i++)
        p_not_measured_msg.data[i] = 100 * (1 - mi.p_not_measured()[i]);
      p_not_measured_pub.publish(p_not_measured_msg);
    }

  private:
    ros::NodeHandle n;
    ros::Subscriber map_sub, click_sub;
    ros::Publisher mi_pub, p_not_measured_pub;

    // Parameters
    double poisson_rate, unknown_threshold;
    bool beam_independence;
    int spatial_jitter, num_beams;
    int condition_steps;
    double mi_max;

    // Map data
    nav_msgs::MapMetaData map_info;
    std::string map_frame_id;
    std::vector<wandering_robot::OccupancyState> map;

    // Computation devices
    wandering_robot::GridMutualInformation mi;
};

int main(int argc, char ** argv) {
  ros::init(argc, argv, "mutual_information_visualizer");
  MutualInformationVisualizer mi;
  ros::spin();
  return 0;
}
