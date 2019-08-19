#include <ros/ros.h>

#include <random>

#include <nav_msgs/OccupancyGrid.h>
#include <range_mi/MIGrid.h>

#include <range_mi/grid_line.hpp>

unsigned int num_iterations = 999;

class MonteCarloMI {

  public:
    MonteCarloMI() {
      // Initialize the uniform distribution
      gen = std::mt19937(random_device());
      dist = std::uniform_real_distribution<double>(0, 1);

      // Initialize the node handle
      n = ros::NodeHandle("~");

      // Fetch the ROS parameters
      std::string map_topic, mi_topic;
      n.getParam("map_topic", map_topic);
      n.getParam("mi_topic", mi_topic);
      // Ray tracing parameters
      n.getParam("num_beams", num_beams);
      // Visualization vvv
      std::string mi_map_topic,
                  binary_map_topic, conditional_map_topic;
      n.getParam("visualize", visualize);
      n.getParam("visualize_more", visualize_more);
      n.getParam("mi_map_topic", mi_map_topic);
      n.getParam("binary_map_topic", binary_map_topic);
      n.getParam("conditional_map_topic", conditional_map_topic);

      // Construct a publisher for mutual information
      mi_pub = n.advertise<range_mi::MIGrid>(mi_topic, 1, true);
      mi_map_pub = n.advertise<nav_msgs::OccupancyGrid>(mi_map_topic, 1, true);
      binary_map_pub = n.advertise<nav_msgs::OccupancyGrid>(binary_map_topic, 1, true);
      conditional_map_pub = n.advertise<nav_msgs::OccupancyGrid>(conditional_map_topic, 1, true);

      // Subscribe to the map
      map_sub = n.subscribe(map_topic, 1, &MonteCarloMI::map_callback, this);
    }

    void map_callback(const nav_msgs::OccupancyGrid & map_msg) {
      // Store map information
      map_info = map_msg.info;
      map_header = map_msg.header;

      // Convert to probability
      vacancy = std::vector<double>(map_info.height * map_info.width);
      for (unsigned int i = 0; i < vacancy.size(); i++) {
        vacancy[i] = 1 - map_msg.data[i]/99.;
        if (vacancy[i] < 0  or vacancy[i] > 1) {
          vacancy[i] = 0;
        }

        vacancy[i] = std::pow(vacancy[i], map_info.resolution);
      }

      // Also initialize a place for randomized and conditional vacancies
      binary_vacancy = std::vector<double>(map_info.height * map_info.width);
      conditional_vacancy = std::vector<double>(map_info.height * map_info.width);
      cconditional_vacancy = std::vector<double>(map_info.height * map_info.width);
      conditional_vacancy_viz = std::vector<double>(map_info.height * map_info.width);

      // And a place for drawing lines
      line = std::vector<unsigned int>(2 * std::max(map_info.height, map_info.width));
      widths = std::vector<double>(2 * std::max(map_info.height, map_info.width));

      // And the mutual information map
      mi = std::vector<double>(map_info.height * map_info.width, 0);

      compute_mi();
    }

    static double entropy(const std::vector<double> & v) {
      double e = 0;
      for (unsigned int i = 0; i < v.size(); i++) {
        if (v[i] > 0 and v[i] < 1) {
          // Accumulate the binary entropy
          e -=
            v[i] * std::log(v[i]) +
            (1 - v[i]) * std::log(1 - v[i]);
        }
      }

      return e;
    }

    void compute_mi() {
      // Clear the mutual information
      std::fill(mi.begin(), mi.end(), 0);

      // Compute the entropy of the map
      double map_entropy = entropy(vacancy);

      for (unsigned int it = 0; it < num_iterations and ros::ok(); it++) {
        // Randomly assign binary vacancy values based
        // on the vacancy probabilities
        for (unsigned int cell = 0; cell < vacancy.size(); cell++) {
          double random = dist(gen);
          if (random < vacancy[cell]) {
            binary_vacancy[cell] = 1;
          } else {
            binary_vacancy[cell] = 0;
          }
        }

        // Initialize the map to equal the original vacancies
        for (unsigned int i = 0; i < vacancy.size(); i++) {
          cconditional_vacancy[i] = vacancy[i];
        }

        // If we are conditioning
        if (false) {
          // Compute a scan at the conditioned point
          for (int beam = 0; beam < num_beams; beam++) {
            double theta = beam * (2 * M_PI)/(num_beams);
            range_mi::grid_line::draw(
                map_info.height,
                map_info.width,
                73.5, 44.5, theta,
                line.data(),
                widths.data(),
                num_cells);

            // Iterate over the beam
            for (unsigned int i = 0; i < num_cells; i++) {
              unsigned int line_cell = line[i];

              cconditional_vacancy[line_cell] = binary_vacancy[line_cell];

              // Stop if we have reached an occupied cell
              if (binary_vacancy[line_cell] == 0) break;
            }
          }

          // Compute the entropy of the map
          map_entropy = entropy(cconditional_vacancy);
        }

        // Given the randomized "ground truth" values compute
        // the result of a scan in the map and how the entropy changed.
        for (unsigned int cell = 0; cell < vacancy.size(); cell++) {
          // Initialize the map to equal the conditioned vacancies
          for (unsigned int i = 0; i < vacancy.size(); i++) {
            conditional_vacancy[i] = cconditional_vacancy[i];
          }

          // For each beam cast rays and update the map
          for (int beam = 0; beam < num_beams; beam++) {
            double theta = beam * (2 * M_PI)/(num_beams);

            double y = std::floor(cell/map_info.width);
            double x = cell - y * map_info.width;
            y += 0.5;
            x += 0.5;

            // Draw a beam
            range_mi::grid_line::draw(
                map_info.height,
                map_info.width,
                x, y, theta,
                line.data(),
                widths.data(),
                num_cells);

            // Iterate over the beam
            for (unsigned int i = 0; i < num_cells; i++) {
              unsigned int line_cell = line[i];

              conditional_vacancy[line_cell] = binary_vacancy[line_cell];

              // Stop if we have reached an occupied cell
              if (binary_vacancy[line_cell] == 0) break;
            }
          }

          // Choose the center cell for visualization purposes
          if (cell == (map_info.height * map_info.width)/2 + map_info.width/2) {
            for (unsigned int i = 0; i < vacancy.size(); i++) {
              conditional_vacancy_viz[i] = conditional_vacancy[i];
            }
            conditional_vacancy_viz[cell] = 0;
          }

          // Compute the entropy of the new map
          double conditional_map_entropy = entropy(conditional_vacancy);

          // Accumulate
          mi[cell] += (map_entropy - conditional_map_entropy);
        }

        // Visualize
        if (visualize and visualize_more) draw_map();

        // Plot
        publish_mi();
      }
    }

    void publish_mi() {
      // Construct a message for the mutual information
      range_mi::MIGrid mi_msg;
      mi_msg.header = map_header;
      mi_msg.data = mi;
      mi_msg.height = map_info.height;
      mi_msg.width = map_info.width;

      // Publish
      mi_pub.publish(mi_msg);
    }

    void draw_map() {
      // Convert to an occupancy map for visualization
      nav_msgs::OccupancyGrid mi_map_msg;
      mi_map_msg.header = map_header;
      mi_map_msg.info = map_info;
      mi_map_msg.data = std::vector<int8_t>(mi.size());
      double mi_max = *std::max_element(mi.begin(), mi.end());
      for (size_t i = 0; i < mi.size(); i++) {
        // Normalize between 0 and 1
        double normalized = mi[i]/mi_max;
        // Change to 100%
        mi_map_msg.data[i] = 100 * (1 - normalized);
      }
      mi_map_pub.publish(mi_map_msg);

      // Publish the binary map
      nav_msgs::OccupancyGrid binary_map_msg = mi_map_msg;
      for (size_t i = 0; i < mi.size(); i++) {
        binary_map_msg.data[i] = 100 * (1 - binary_vacancy[i]);
      }
      binary_map_pub.publish(binary_map_msg);

      // Publish the conditioned map
      nav_msgs::OccupancyGrid conditional_map_msg = mi_map_msg;
      for (size_t i = 0; i < mi.size(); i++) {
        conditional_map_msg.data[i] = 100 * (1 - conditional_vacancy_viz[i]);
      }
      conditional_map_pub.publish(conditional_map_msg);
    }

  private:
    ros::NodeHandle n;
    ros::Subscriber map_sub;
    ros::Publisher mi_pub;
    ros::Publisher
      mi_map_pub,
      binary_map_pub,
      conditional_map_pub;

    // Random number generator
    std::random_device random_device;
    std::mt19937 gen;
    std::uniform_real_distribution<double> dist;

    // Parameters
    int num_beams;
    bool visualize, visualize_more;

    // Map data
    nav_msgs::MapMetaData map_info;
    std_msgs::Header map_header;
    std::vector<double> mi;
    std::vector<double> vacancy, binary_vacancy, conditional_vacancy, cconditional_vacancy, conditional_vacancy_viz;

    // Line data
    std::vector<unsigned int> line;
    std::vector<double> widths;
    unsigned int num_cells;
};

int main(int argc, char ** argv) {
  ros::init(argc, argv, "monte_carlo_mi");
  MonteCarloMI mc_mi;
  ros::spin();
  return 0;
}
