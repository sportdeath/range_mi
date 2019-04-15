#include <ros/ros.h>

#include <nav_msgs/OccupancyGrid.h>
#include <geometry_msgs/PointStamped.h>

#include "wandering_robot/bresenham.hpp"

unsigned int height = 100;
unsigned int width = 200;
double sample_rate = 0.000001;
double refresh_rate = 0.04;

class SampleBeams {

  public:

    SampleBeams() {
      // Initialize the map
      map = std::vector<double>(height * width, 0);
      line = std::vector<unsigned int>(std::max(height, width));

      // Initialize the beam sampler
      bresenham = wandering_robot::Bresenham(height, width);

      // Construct a publisher
      map_pub = n.advertise<nav_msgs::OccupancyGrid>("/map", 1, true);

      // Start a timer
      sample_timer = n.createTimer(ros::Duration(sample_rate), &SampleBeams::sample_callback, this);

      // Refresh the map
      map_timer = n.createTimer(ros::Duration(refresh_rate), &SampleBeams::draw_map, this);
    }

    void sample_callback(const ros::TimerEvent & event) {
      // Randomly sample a point
      double x, y, theta;
      bresenham.sample(x, y, theta);

      // Compute Bresenham's line
      unsigned int num_cells;
      bresenham.line(
          y, x, theta,
          line.data(),
          num_cells);
      // Add it to the map
      for (unsigned int i = 0; i < num_cells; i++) {
        map[line[i]] += 1;
      }
    }

    void draw_map(const ros::TimerEvent & event) {
      // Convert to int8
      nav_msgs::OccupancyGrid map_msg;
      map_msg.data = std::vector<int8_t>(height * width);
      double map_max = *std::max_element(map.begin(), map.end());
      for (size_t i = 0; i < map_msg.data.size(); i++) {
        map_msg.data[i] = 100 * map[i]/map_max;
      }

      // Add info
      map_msg.header.stamp = ros::Time::now();
      map_msg.header.frame_id = "map";
      map_msg.info.resolution = 1;
      map_msg.info.width = width;
      map_msg.info.height = height;
      map_msg.info.origin.orientation.w = 1;

      // Publish
      map_pub.publish(map_msg);
    }

  private:
    double x, y, theta;
    std::vector<double> map;
    std::vector<unsigned int> line;
    wandering_robot::Bresenham bresenham;

    ros::NodeHandle n;
    ros::Timer sample_timer, map_timer;
    ros::Publisher map_pub;
};

int main(int argc, char ** argv) {
  ros::init(argc, argv, "bresenham_test");
  SampleBeams bt;
  ros::spin();
  return 0;
}
