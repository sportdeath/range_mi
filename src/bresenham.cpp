#include <cmath>
#include <utility>

#include "wandering_robot/bresenham.hpp"

void wandering_robot::Bresenham::line(
    double col,
    double row,
    double theta,
    unsigned int * const line,
    unsigned int & num_cells) {

  // Initialize constants
  num_cells = 0;
  double s = std::sin(theta);
  double c = std::cos(theta);
  int col_step = (c > 0) ? 1 : -1;
  int row_step = (s > 0) ? 1 : -1;

  // Swap for symmetry
  double error, derror;
  double * col_ = &col;
  double * row_ = &row;
  if (std::abs(s) < 1/std::sqrt(2)) {
    // Predominant motion is along the x axis
    error = std::fmod(row, 1.);
    derror = std::abs(s/c);
  } else {
    // Predominant motion is along the y axis
    error = std::fmod(col, 1.);
    derror = std::abs(c/s);
    std::swap(col_, row_);
    std::swap(row_step, col_step);
  }

  while (true) {

    if (row < 0 or
        col < 0 or
        row >= height or
        col >= width)
      // Outside the map, stop casting!
      break;

    // Add the hit cell
    line[num_cells++] = ((unsigned int) row) * width + ((unsigned int) col);

    // Move to the next cell
    error += derror;
    *col_ += col_step;
    if (error > 1) {
      *row_ += row_step;
      error -= 1;
    }
  }
}

void wandering_robot::Bresenham::sample(
    double & x,
    double & y,
    double & theta) {

  // Only do the first quadrant for now
  theta = (2 * M_PI) * dist(gen);
  double s = std::sin(theta);
  double c = std::cos(theta);
  double s_abs = std::abs(s);
  double c_abs = std::abs(c);

  double normalized_width = s_abs * width;
  double normalized_height = c_abs * height;

  double choice = (normalized_width + normalized_height) * dist(gen);
  if (choice < normalized_width) {
    // Vertical
    x = choice/s_abs;
    if (s < 0) {
      // Facing down
      y = height - 1;
    } else {
      y = 0;
    }
  } else {
    // Horizontal
    y = (choice - normalized_width)/c_abs;
    if (c < 0) {
      // Facing left
      x = width - 1;
    } else {
      x = 0;
    }
  }
}
