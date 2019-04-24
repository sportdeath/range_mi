#include <cmath>
#include <utility>

#include "wandering_robot/grid_line.hpp"

void wandering_robot::GridLine::draw(
    double x,
    double y,
    double theta,
    unsigned int * const line,
    double * const widths,
    unsigned int & num_cells) {

  // Pre-compute sin and cos
  double sin = std::sin(theta);
  double cos = std::cos(theta);

  // Reset the line
  num_cells = 0;

  // Initialize the cell indices
  int floor_x = std::floor(x);
  int floor_y = std::floor(y);

  // Initialize the distances along the line to the
  // next horizontal cell boundary (dx) or vertical
  // cell boundary (dy). The maximum distance for
  // either is dx/y_unit.
  double dx_unit = 1/std::abs(cos);
  double dy_unit = 1/std::abs(sin);
  double dx = dx_unit * ((cos > 0) ? floor_x + 1 - x : x - floor_x);
  double dy = dy_unit * ((sin > 0) ? floor_y + 1 - y : y - floor_y);

  // Compute the sign of the steps taken
  int x_step = (cos > 0) ? 1 : -1;
  int y_step = (sin > 0) ? 1 : -1;

  // While we are inside the map
  while (
      floor_x >= 0 and
      floor_y >= 0 and
      floor_x < (int) width and
      floor_y < (int) height) {

    // Add the cell index
    line[num_cells] = floor_y * width + floor_x;

    // The width is the minimum distance to either cell
    widths[num_cells] = std::min(dx, dy);

    // Subtract the width from the line boundaries
    // and the distance to the boundary
    dx -= widths[num_cells];
    dy -= widths[num_cells];

    // Replenish if necessary
    if (dx <= 0) {
      dx = dx_unit;
      floor_x += x_step;
    }
    if (dy <= 0) {
      dy = dy_unit;
      floor_y += y_step;
    }

    // Increment the cell number
    num_cells++;
  }
}

void wandering_robot::GridLine::sample(
    double & x,
    double & y,
    double & theta) {
  sample_regularly(x, y, theta, dist(gen), dist(gen));
}

void wandering_robot::GridLine::sample_regularly(
    double & x,
    double & y,
    double & theta,
    double spatial_interpolation,
    double angular_interpolation) {

  // Calculate theta
  theta = 2 * M_PI * angular_interpolation;

  // Pre-compute trig values
  double s = std::sin(theta);
  double c = std::cos(theta);
  double s_abs = std::abs(s);
  double c_abs = std::abs(c);

  // Compute the normalized spatial dimensions
  double normalized_width = s_abs * width;
  double normalized_height = c_abs * height;

  double perimeter =
    (normalized_width + normalized_height) *
    spatial_interpolation;

  if (perimeter < normalized_width) {
    x = perimeter/s_abs;
    if (s < 0) {
      // Facing down
      y = height - 0.0000001;
    } else {
      y = 0;
    }
  } else {
    // Horizontal
    y = (perimeter - normalized_width)/c_abs;
    if (c < 0) {
      // Facing left
      x = width - 0.0000001;
    } else {
      x = 0;
    }
  }
}
