#pragma once

#include <cmath>

/**
 * Computes a variety of quantities in expectation
 */

namespace range_entropy {
namespace expected {

const double vacancy_min = 0.000000001;
const double pdf_max = -std::log(vacancy_min);
const double entropy_min = -std::log(pdf_max);

/**
 * The functions this class aims to
 * evaluate in expectation over range measurements
 * r with probability density functions f.
 */
double p_not_measured(double r, double f);
double distance1(double r, double f);
double distance2(double r, double f);
double neg_log(double f);
double information1(double r, double f);
double information2(double r, double f);
double information3(double r, double f);

struct local {
  double width;
  double p_not_measured;
  double miss_p;
  double miss_info;
  double neg_log_p_free_inv;
  double hit_p;
  double zero_info;
};

struct expected {
  double information3;
  double information2;
  double information1;
  double distance2;
  double distance1;
  double p_not_measured;
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
 * The expected value of the quantity if a
 * hit occurs times p_hit
 */
double p_not_measured_hit(const local & l);
double distance1_hit(const local & l);
double distance2_hit(const local & l);
double information1_hit(const local & l);
double information2_hit(const local & l);
double information3_hit(const local & l);

/**
 * The expected value the quantity if a
 * miss occurs (not times p_miss)
 */
double p_not_measured_miss(const local & l, const expected & exp);
double distance1_miss(const local & l, const expected & exp);
double distance2_miss(const local & l, const expected & exp);
double information1_miss(const local & l, const expected & exp);
double information2_miss(const local & l, const expected & exp);
double information3_miss(const local & l, const expected & exp);

/**
 * The expected value of a quantity if
 * a hit or miss occurs. The hit and miss
 * functions for the particular value are
 * passed as parameters
 */
double hit_or_miss(
    const local & l,
    const expected & exp,
    double (*hit)(const local &),
    double (*miss)(const local &, const expected &));

/**
 * The expected values of a quantity
 * along all cells in a line
 */
void line(
    const unsigned int * const line,
    const double * const p_free,
    const double * const p_not_measured,
    const double * const width,
    unsigned int num_cells,
    bool entropy,
    unsigned int dimension,
    double * const output);

void p_not_measured(
    const unsigned int * const line,
    const double * const vacancy,
    const double * const width,
    unsigned int num_cells,
    unsigned int dimension,
    double * const p_not_measured1_);

}}
