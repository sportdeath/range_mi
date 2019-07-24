#include <ros/ros.h>

#include <nav_msgs/OccupancyGrid.h>
#include <range_entropy/EntropyMap.h>
#include <geometry_msgs/PointStamped.h>

#include "range_entropy/grid_expected.hpp"
#include "range_entropy/expected.hpp"

class EntropyVisualizer {

  public:
    EntropyVisualizer() {
      // Initialize the node handle
      n = ros::NodeHandle("~");

      // Fetch the ROS parameters
      std::string map_topic, entropy_topic, condition_topic, reset_condition_topic;
      n.getParam("map_topic", map_topic);
      n.getParam("entropy_topic", entropy_topic);
      n.getParam("condition_topic", condition_topic);
      n.getParam("reset_condition_topic", reset_condition_topic);
      // Ray tracing parameters
      n.getParam("spatial_jitter", spatial_jitter);
      n.getParam("num_beams", num_beams);
      n.getParam("condition_steps", condition_steps);
      // Noise parameters
      n.getParam("noise_dev", noise_dev);
      n.getParam("noise_width", noise_width);
      n.getParam("noise_step_size", noise_step_size);
      // Visualization vvv
      std::string click_condition_topic, entropy_viz_topic, p_not_measured_viz_topic;
      n.getParam("visualize", visualize);
      n.getParam("visualize_more", visualize_more);
      n.getParam("click_condition_topic", click_condition_topic);
      n.getParam("entropy_viz_topic", entropy_viz_topic);
      n.getParam("p_not_measured_viz_topic", p_not_measured_viz_topic);

      // Construct a publisher for mutual information
      entropy_pub = n.advertise<range_entropy::EntropyMap>(entropy_topic, 1, true);
      entropy_map_pub = n.advertise<nav_msgs::OccupancyGrid>(entropy_viz_topic, 1, true);
      p_not_measured_pub = n.advertise<nav_msgs::OccupancyGrid>(p_not_measured_viz_topic, 1, true);

      // Subscribe to maps and clicked points
      map_sub = n.subscribe(map_topic, 1, &EntropyVisualizer::map_callback, this);
      click_sub = n.subscribe(click_condition_topic, 1, &EntropyVisualizer::click_callback, this);
    }

    void map_callback(const nav_msgs::OccupancyGrid & map_msg) {
      // Store map information
      map_info = map_msg.info;
      map_header = map_msg.header;

      // Convert to probability
      vacancy = std::vector<double>(map_info.height * map_info.width);
      for (unsigned int i = 0; i < vacancy.size(); i++) {
        vacancy[i] = 1 - map_msg.data[i]/100.;
        if (vacancy[i] < 0  or vacancy[i] > 1) {
          std::cout << "THIS IS BAD" << std::endl;
          vacancy[i] = 0;
        }

        vacancy[i] = std::pow(vacancy[i], 0.1);
      }

      // Initialize mutual information computation on the grid
      grid_caster = range_entropy::GridExpected(
          map_info.height,
          map_info.width,
          spatial_jitter,
          num_beams);

      compute_entropy();
    }

    void compute_entropy() {
      grid_caster.reset_surface();

      double spatial_interpolation = 0;
      double angular_interpolation = 0;
      while (angular_interpolation < 1) {
        // Add the beam
        grid_caster.compute_surface_beam(
            spatial_interpolation,
            angular_interpolation,
            vacancy.data(),
            true,
            2,
            noise_dev,
            noise_width,
            noise_step_size);

        // Draw every time a spatial section is completed
        if (spatial_interpolation == 0) {
          if (visualize and visualize_more) {
            draw_map();
          }
          //grid_caster.reset_surface();
          if (not ros::ok()) break;
        }
      }

      if (visualize) draw_map();
      publish_entropy();
    }

    void click_callback(const geometry_msgs::PointStamped & click_msg) {
      // Condition the map on the clicked point
      double x = click_msg.point.x;
      double y = click_msg.point.y;
      grid_caster.condition(
          x, y, vacancy.data(), 
          condition_steps, 2);

      // Update the mutual information
      compute_entropy();
    }

    void publish_entropy() {
      // Construct a message for the mutual information surface
      range_entropy::EntropyMap entropy_msg;
      entropy_msg.header = map_header;
      entropy_msg.data = grid_caster.surface();
      entropy_msg.height = map_info.height;
      entropy_msg.width = map_info.width;

      // Publish
      entropy_pub.publish(entropy_msg);
    }

    void draw_map() {
      // Convert to an occupancy map for visualization
      nav_msgs::OccupancyGrid entropy_map_msg;
      entropy_map_msg.header = map_header;
      entropy_map_msg.info = map_info;
      entropy_map_msg.data = std::vector<int8_t>(grid_caster.surface().size());
      double entropy_max = *std::max_element(grid_caster.surface().begin(), grid_caster.surface().end());
      double entropy_min = *std::min_element(grid_caster.surface().begin(), grid_caster.surface().end());
      for (size_t i = 0; i < grid_caster.surface().size(); i++) {
        // Normalize between 0 and 1
        double normalized = (grid_caster.surface()[i] - entropy_min)/(entropy_max - entropy_min);
        // Change to 100%
        entropy_map_msg.data[i] = 100 * (1 - normalized);
      }
      entropy_map_pub.publish(entropy_map_msg);

      // Do the same with p not measured
      nav_msgs::OccupancyGrid p_not_measured_msg = entropy_map_msg;
      for (size_t i = 0; i < p_not_measured_msg.data.size(); i++)
        p_not_measured_msg.data[i] = 100 * (1 - grid_caster.p_not_measured()[i]);
      p_not_measured_pub.publish(p_not_measured_msg);
    }

  private:
    ros::NodeHandle n;
    ros::Subscriber map_sub, click_sub;
    ros::Publisher entropy_pub,
      entropy_map_pub,
      p_not_measured_pub;

    // Parameters
    int spatial_jitter, num_beams;
    int condition_steps;
    bool visualize, visualize_more;
    double noise_dev, noise_width, noise_step_size;

    // Map data
    nav_msgs::MapMetaData map_info;
    std_msgs::Header map_header;
    std::vector<double> vacancy;

    // Computation devices
    range_entropy::GridExpected grid_caster;
};

int main(int argc, char ** argv) {
  ros::init(argc, argv, "entropy_visualizer");
  EntropyVisualizer ev;
  ros::spin();
  return 0;
}
