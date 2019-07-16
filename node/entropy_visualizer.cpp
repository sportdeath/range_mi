#include <ros/ros.h>

#include <nav_msgs/OccupancyGrid.h>
#include <geometry_msgs/PointStamped.h>

#include "range_entropy/grid_expected.hpp"
#include "range_entropy/expected.hpp"

class EntropyVisualizer {

  public:
    EntropyVisualizer() {
      // Initialize the node handle
      n = ros::NodeHandle("~");

      // Fetch the ROS parameters
      std::string entropy_topic, p_not_measured_topic, click_topic, map_topic;
      n.getParam("entropy_topic", entropy_topic);
      n.getParam("p_not_measured_topic", p_not_measured_topic);
      n.getParam("click_topic", click_topic);
      n.getParam("map_topic", map_topic);
      n.getParam("spatial_jitter", spatial_jitter);
      n.getParam("num_beams", num_beams);
      n.getParam("condition_steps", condition_steps);

      // Construct a publisher for mutual information
      entropy_pub = n.advertise<nav_msgs::OccupancyGrid>(entropy_topic, 1, true);
      p_not_measured_pub = n.advertise<nav_msgs::OccupancyGrid>(p_not_measured_topic, 1, true);

      // Subscribe to maps and clicked points
      map_sub = n.subscribe(map_topic, 1, &EntropyVisualizer::map_callback, this);
      click_sub = n.subscribe(click_topic, 1, &EntropyVisualizer::click_callback, this);
    }

    void map_callback(const nav_msgs::OccupancyGrid & map_msg) {
      // Store map information
      map_info = map_msg.info;
      map_frame_id = map_msg.header.frame_id;

      // Convert to probability
      p_free = std::vector<double>(map_info.height * map_info.width);
      for (unsigned int i = 0; i < p_free.size(); i++) {
        p_free[i] = 1 - map_msg.data[i]/100.;
        if (p_free[i] < 0 or p_free[i] > 1) {
          std::cout << "Map contains invalid data" << std::endl;
          p_free[i] = 0.5;
        }

        p_free[i] = std::pow(p_free[i], 0.1);
      }

      // Initialize mutual information computation on the grid
      grid_caster = range_entropy::GridExpected(
          map_info.height,
          map_info.width,
          spatial_jitter,
          num_beams);

      compute_entropy(true);
    }

    void compute_entropy(bool first) {
      grid_caster.reset_surface();

      double spatial_interpolation = 0;
      double angular_interpolation = 0;
      while (angular_interpolation < 1) {
        // Add the beam
        grid_caster.compute_surface_beam(
            spatial_interpolation,
            angular_interpolation,
            p_free.data(),
            range_entropy::expected::information2);

        // Draw every time a spatial section is completed
        if (spatial_interpolation == 0) {
          if (first) {
            entropy_max = *std::max_element(grid_caster.surface().begin(), grid_caster.surface().end());
            entropy_min = *std::min_element(grid_caster.surface().begin(), grid_caster.surface().end());
          }
          draw_map();
          //grid_caster.reset_surface();
          if (not ros::ok()) break;
        }

      }
    }

    void click_callback(const geometry_msgs::PointStamped & click_msg) {
      // Condition the map on the clicked point
      double x = click_msg.point.x;
      double y = click_msg.point.y;
      grid_caster.condition(x, y, p_free.data());

      // Update the mutual information
      compute_entropy(false);
    }

    void draw_map() {
      // Construct a message for the mutual information surface
      nav_msgs::OccupancyGrid entropy_msg;
      entropy_msg.data = std::vector<int8_t>(grid_caster.surface().size());
      for (size_t i = 0; i < grid_caster.surface().size(); i++) {
        double normalized = (grid_caster.surface()[i] - entropy_min)/(entropy_max - entropy_min);
        entropy_msg.data[i] = 100 * (1 - normalized);
      }

      // Add info
      entropy_msg.header.frame_id = map_frame_id;
      entropy_msg.header.stamp = ros::Time::now();
      entropy_msg.info = map_info;

      // Publish
      entropy_pub.publish(entropy_msg);

      // Do the same for p_not measured
      nav_msgs::OccupancyGrid p_not_measured_msg = entropy_msg;
      for (size_t i = 0; i < p_not_measured_msg.data.size(); i++)
        p_not_measured_msg.data[i] = 100 * (1 - grid_caster.p_not_measured()[i]);
      p_not_measured_pub.publish(p_not_measured_msg);
    }

  private:
    ros::NodeHandle n;
    ros::Subscriber map_sub, click_sub;
    ros::Publisher entropy_pub, p_not_measured_pub;

    // Parameters
    int spatial_jitter, num_beams;
    int condition_steps;

    // Map data
    nav_msgs::MapMetaData map_info;
    std::string map_frame_id;
    std::vector<double> p_free;
    double entropy_min, entropy_max;

    // Computation devices
    range_entropy::GridExpected grid_caster;
};

int main(int argc, char ** argv) {
  ros::init(argc, argv, "entropy_visualizer");
  EntropyVisualizer ev;
  ros::spin();
  return 0;
}
