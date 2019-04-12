#pragma once

namespace wandering_robot {

void bresenham(
    double row,
    double col,
    double theta,
    unsigned int height,
    unsigned int width,
    unsigned int * const line,
    unsigned int & num_cells);

}
