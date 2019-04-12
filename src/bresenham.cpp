#include <cmath>
#include <utility>

#include "wandering_robot/bresenham.hpp"

void wandering_robot::bresenham(
    double row,
    double col,
    double theta,
    unsigned int height,
    unsigned int width,
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

    if (row < 0 ||
        col < 0 ||
        row >= height ||
        col >= width)
      // Outside the map, stop casting!
      break;

    // Compute the hit entry
    unsigned int hit_cell = std::floor(row) * width + std::floor(col);

    // Add the hit cell
    line[num_cells] = hit_cell;
    num_cells++;

    // Move to the next cell
    error += derror;
    *col_ += col_step;
    if (error > 1) {
      *row_ += row_step;
      error -= 1;
    }
  }
}
