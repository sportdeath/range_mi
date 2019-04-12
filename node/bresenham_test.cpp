#include <ros/ros.h>

#include <nav_msgs/OccupancyGrid.h>
#include <geometry_msgs/PointStamped.h>

#include "wandering_robot/bresenham.hpp"

unsigned int height = 100;
unsigned int width = 200;
double theta_per_second = 1;
double refresh_rate = 0.04;

class BresenhamTest {

  public:

    BresenhamTest() {
      // Initialize the line point
      x = 0; y = 0; theta = 0;

      // Construct a publisher
      map_pub = n.advertise<nav_msgs::OccupancyGrid>("/map", 1, true);

      // Subscribe to clicked points
      click_sub = n.subscribe("/clicked_point", 1, &BresenhamTest::click_callback, this);

      // Start a timer
      theta_timer = n.createTimer(ros::Duration(refresh_rate), &BresenhamTest::theta_timer_callback, this);
    }

    void click_callback(const geometry_msgs::PointStamped::ConstPtr & msg) {
      x = msg -> point.x;
      y = msg -> point.y;
    }

    void theta_timer_callback(const ros::TimerEvent & event) {
      theta = ros::Time::now().toSec() * theta_per_second;
      draw_bresenham();
    }

    void draw_bresenham() {
      // Initialize an empty map
      nav_msgs::OccupancyGrid map_msg;
      map_msg.data = std::vector<int8_t>(height * width);
      for (size_t i = 0; i < map_msg.data.size(); i++) {
        map_msg.data[i] = 0;
      }

      // Add info
      map_msg.header.stamp = ros::Time::now();
      map_msg.header.frame_id = "map";
      map_msg.info.resolution = 1;
      map_msg.info.width = width;
      map_msg.info.height = height;
      map_msg.info.origin.orientation.w = 1;

      // Compute Bresenham
      unsigned int num_cells;
      std::vector<unsigned int> line(std::max(height, width));
      wandering_robot::bresenham(
          y, x,
          theta,
          height,
          width,
          line.data(),
          num_cells);

      // Draw on the map
      for (unsigned int i = 0; i < num_cells; i++) {
        map_msg.data[line[i]] = 100;
      }

      map_pub.publish(map_msg);
    }

  private:
    double x, y, theta;

    ros::NodeHandle n;
    ros::Timer theta_timer;
    ros::Subscriber click_sub;
    ros::Publisher map_pub;
};

int main(int argc, char ** argv) {
  ros::init(argc, argv, "bresenham_test");
  BresenhamTest bt;
  ros::spin();
  return 0;
}
