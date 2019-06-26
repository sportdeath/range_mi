#pragma once

/**
 * Computes a variety of quantities in expectation
 */

namespace range_entropy {
namespace expected {

struct local {
  double width;
  double p_not_measured;
  double miss_p;
  double miss_info;
  double miss_info_inv;
  double hit_p;
  double zero_info;
};

/**
 * Pre-computes quantities for
 * a region of size width, constant
 * volumetric probability of being free p_free
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
    double information1_);
double information2(
    const local & l,
    double distance1_,
    double information1_,
    double information2_);
double information3(
    const local & l,
    double distance1_,
    double distance2_,
    double information1_,
    double information2_,
    double information3_);

/**
 * Functions for the above that apply to
 * a vector of values.
 */
void distance1(
    const unsigned int * const line,
    const double * const p_free,
    const double * const p_not_measured,
    const double * const p_width,
    unsigned int num_cells,
    double * const distance1_);
void distance2(
    const unsigned int * const line,
    const double * const p_free,
    const double * const p_not_measured,
    const double * const p_width,
    unsigned int num_cells,
    double * const distance2_);
void information1(
    const unsigned int * const line,
    const double * const p_free,
    const double * const p_not_measured,
    const double * const p_width,
    unsigned int num_cells,
    double * const information1_);
void information2(
    const unsigned int * const line,
    const double * const p_free,
    const double * const p_not_measured,
    const double * const p_width,
    unsigned int num_cells,
    double * const information2_);
void information3(
    const unsigned int * const line,
    const double * const p_free,
    const double * const p_not_measured,
    const double * const p_width,
    unsigned int num_cells,
    double * const information3_);

}}
