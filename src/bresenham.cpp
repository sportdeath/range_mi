#include <cmath>
#include <utility>

#include "wandering_robot/bresenham.hpp"

void wandering_robot::Bresenham::line(
    double row,
    double col,
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

void wandering_robot::Bresenham::sample(
    double & x,
    double & y,
    double & theta) {

  // Choose a random point on the perimeter
  double perimeter = dist_perimeter(gen);

  if (perimeter < width) {
    // On the bottom
    theta = dist_theta(gen);
    x = perimeter;
    y = 0;
  } else if (perimeter < 2 * width) {
    // On the top
    theta = dist_theta(gen) + M_PI;
    x = perimeter - width;
    y = height - 1;
  } else if (perimeter < 2 * width + height) {
    // On the left
    theta = dist_theta(gen) - M_PI/2.;
    y = perimeter - 2 * width;
    x = 0;
  } else {
    // On the right
    theta = dist_theta(gen) + M_PI/2.;
    y = perimeter - 2 * width - height;
    x = width - 1;
  }
}
