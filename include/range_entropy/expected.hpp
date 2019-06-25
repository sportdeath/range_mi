#pragma once

#include <vector>

/**
 * Computes a variety of quantities in expectation
 */

namespace range_entropy {
namespace expected {

struct local {
  double width;
  double p_not_measured;
  double miss_p;
  double miss_odds;
  double miss_info;
  double miss_info_inv;
  double hit_p;
  double zero_info;
};

/**
 * Pre-computes quantities for
 * a region of size width, constant
 * occupancy probability p_occupied,
 * and probability of not being measured
 * by another independent sensor measurement
 * p_not_measured.
 */
void update_local(
    double p_free,
    double p_not_measured,
    double width,
    local & l);

/**
 * The first, second and third moment of
 * distance that the range measurement returns.
 */
double distance1(
    const local & l,
    double distance1_);
double distance2(
    const local & l,
    double distance1_,
    double distance2_);

/**
 * The expected information gain of
 * a range measurement for 1d, 2d and
 * 3d beams.
 *
 * The integrations for 2d and 3d beams
 * have multiplicative terms r and r^2
 * respectively to account for radial
 * expansion.
 */
double information1(
    const local & l,
    double information1);

/**
 * Functions for the above that apply to
 * a vector of values.
 */
std::vector<double> distance1(
    const std::vector<double> & p_free,
    const std::vector<double> & p_not_measured,
    const std::vector<double> & width);
std::vector<double> distance2(
    const std::vector<double> & p_free,
    const std::vector<double> & p_not_measured,
    const std::vector<double> & width);
std::vector<double> information1(
    const std::vector<double> & p_free,
    const std::vector<double> & p_not_measured,
    const std::vector<double> & width);

}}
