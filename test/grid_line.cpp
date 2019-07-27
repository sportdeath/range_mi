#include <vector>
#include <cassert>
#include <cmath>

#include <range_entropy/grid_line.hpp>

using namespace range_entropy;

unsigned int height = 10;
unsigned int width = 10;
unsigned int num_beams = 100;
double eps = 0.0001;

unsigned int max_cells = height * width * 2;
std::vector<unsigned int> line(max_cells);
std::vector<double> widths(max_cells);
unsigned int num_cells;

void axis() {
  GridLine grid_line(height, width, num_beams);

  grid_line.draw(0, 0, 0, line.data(), widths.data(), num_cells);

  for (unsigned int i = 0; i < num_cells; i++) {
    assert(std::abs(widths[i] - 1) < eps);
    assert(line[i] == i);
  }

  grid_line.draw(10 - (eps * eps), 0, M_PI, line.data(), widths.data(), num_cells);

  for (unsigned int i = 0; i < num_cells; i++) {
    assert(std::abs(widths[i] - 1) < eps);
    assert(line[i] == 9 - i);
  }
}

void diagonal() {
  GridLine grid_line(height, width, num_beams);
  grid_line.draw(0, 0, M_PI/4., line.data(), widths.data(), num_cells);

  unsigned int ii = 0;
  for (unsigned int i = 0; i < num_cells; i++) {
    if (widths[i] > eps) {
      assert(std::abs(widths[i] - std::sqrt(2)) < eps);
      assert(line[i] == 10 * ii + ii);
      ii++;
    }
  }
  assert(ii == 10);

  grid_line.draw(10 - (eps * eps), 10 - (eps * eps), - 3. * M_PI/4., line.data(), widths.data(), num_cells);

  ii = 0;
  for (unsigned int i = 0; i < num_cells; i++) {
    if (widths[i] > eps) {
      assert(std::abs(widths[i] - std::sqrt(2)) < eps);
      assert(line[i] == 10 * (10 - 1 - ii) + (10 - 1 - ii));
      ii++;
    }
  }
  assert(ii == 10);

  grid_line.draw(0, 10 - (eps * eps), - M_PI/4., line.data(), widths.data(), num_cells);

  ii = 0;
  for (unsigned int i = 0; i < num_cells; i++) {
    if (widths[i] > eps) {
      assert(std::abs(widths[i] - std::sqrt(2)) < eps);
      assert(line[i] == 10 * (10 - 1 - ii) + ii);
      ii++;
    }
  }
  assert(ii == 10);

  grid_line.draw(10 - (eps * eps), 0, 3 * M_PI/4., line.data(), widths.data(), num_cells);

  ii = 0;
  for (unsigned int i = 0; i < num_cells; i++) {
    if (widths[i] > eps) {
      assert(std::abs(widths[i] - std::sqrt(2)) < eps);
      assert(line[i] == 10 * ii + (10 - 1 - ii));
      ii++;
    }
  }
  assert(ii == 10);
}

int main() {
  diagonal();
  axis();
}
