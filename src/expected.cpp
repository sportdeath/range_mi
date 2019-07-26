#include <cmath>

#include "range_entropy/expected.hpp"

using namespace range_entropy;
using namespace expected;

double range_entropy::expected::p_not_measured(double r, double f) {
  (void) f;
  (void) r;
  return 1;
}
double range_entropy::expected::distance1(double r, double f) {
  (void) f;
  return r;
}
double range_entropy::expected::distance2(double r, double f) {
  (void) f;
  return r*r;
}
double range_entropy::expected::neg_log(double f) {
  if (f <= 0) {
    return 999999;
  } else {
    return -std::log(f);
  }
}
double range_entropy::expected::information1(double r, double f) {
  (void) r;
  return neg_log(f);
}
double range_entropy::expected::information2(double r, double f) {
  return r * neg_log(f);
}
double range_entropy::expected::information3(double r, double f) {
  return r * r * neg_log(f);
}

void range_entropy::expected::update_local(
    double p_free,
    double p_not_measured,
    double width,
    local & l) {

  // Clip the inputs to be defined.
  if (p_free >= 1) {
    p_free = 0.9999999;
  } else if (p_free <= 0) {
    p_free = 0.0000001;
  }

  if (width <= 0) {
    width = 0.0000001;
  }

  l.width = width;
  l.p_not_measured = p_not_measured;

  // Compute once
  double log_p_free = std::log(p_free);

  // The probability that the range measurement
  // does not hit any occupied state within the
  // width.
  //
  //      w
  // 1 - int f(r) dr = p^w
  //      0
  //
  // Note that f(r) = -p^r log(p)
  // 
  l.miss_p = std::exp(log_p_free * width);

  // The converse of p_miss
  l.hit_p = 1 - l.miss_p;

  // The amount of information learned if
  // a miss occurs
  //
  // -log(p_miss) = -w log p
  // 
  l.miss_info = -width * log_p_free;

  // The inverse of the information
  // precomputed for optimization
  l.miss_info_inv = 1./l.miss_info;

  // The information gained by hitting
  // an infinitesimal region.
  //
  //                      w
  //  lim (1 - p^w)^(-1) int - f(r) log(f(r)) dr
  // w->0                 0
  //
  //        = -log(-log(p))
  l.zero_info = -std::log(-log_p_free);
}

double range_entropy::expected::p_not_measured_hit(
    const local & l) {
  return l.hit_p;
}

double range_entropy::expected::distance1_hit(
    const local & l) {
  return l.width * (l.hit_p * l.miss_info_inv - l.miss_p);
}

double range_entropy::expected::distance2_hit(
    const local & l) {
  return l.width * l.width * (
      2 * l.miss_info_inv * l.miss_info_inv * l.hit_p - l.miss_p * (l.miss_info + 2) * l.miss_info_inv);
}

double range_entropy::expected::information1_hit(
    const local & l) {
  return l.hit_p * (l.zero_info + 1) - l.miss_p * l.miss_info;
}

double range_entropy::expected::information2_hit(
    const local & l) {
  return l.width * (
      l.hit_p * (2 + l.zero_info) * l.miss_info_inv - l.miss_p * (2 + l.miss_info + l.zero_info));
}

double range_entropy::expected::information3_hit(
    const local & l) {
  return l.width * l.width * (
      l.hit_p * (6 + 2 * l.zero_info) * l.miss_info_inv * l.miss_info_inv -
      l.miss_p * l.miss_info_inv * (6 + l.miss_info + (l.miss_info + l.zero_info) * (l.miss_info + 2)));
}

double range_entropy::expected::p_not_measured_miss(
    const local & l,
    const expected & exp) {
  (void) l;
  return exp.p_not_measured;
}

double range_entropy::expected::distance1_miss(
    const local & l,
    const expected & exp) {
  // If the range misses the local cell,
  // the expected distance is simply the old
  // expected distance plus the width of the cell
  return exp.distance1 + exp.p_not_measured * l.width;
}

double range_entropy::expected::distance2_miss(
    const local & l,
    const expected & exp) {
  return exp.distance2 + 2 * l.width * exp.distance1 + exp.p_not_measured * l.width * l.width;
}

double range_entropy::expected::information1_miss(
    const local & l,
    const expected & exp) {
  // If the range measurement misses
  // we gain all the information from the far
  // measurement plus the information gained
  // by the fact that we missed.
  return exp.information1 + exp.p_not_measured * l.miss_info;
}

double range_entropy::expected::information2_miss(
    const local & l,
    const expected & exp) {
  return
    exp.information2 + l.miss_info * exp.distance1 +
    l.width * (exp.information1 + exp.p_not_measured * l.miss_info);
}

double range_entropy::expected::information3_miss(
    const local & l,
    const expected & exp) {
  return
    exp.information3 + l.miss_info * exp.distance2 +
    2 * l.width * (exp.information2 + l.miss_info * exp.distance1) +
    l.width * l.width * (exp.information1 + exp.p_not_measured * l.miss_info);
}

double range_entropy::expected::hit_or_miss(
    const local & l,
    const expected & exp,
    double (*hit)(const local &),
    double (*miss)(const local &, const expected &)) {

  return l.p_not_measured * hit(l) + l.miss_p * miss(l, exp);
}

void range_entropy::expected::line(
    const unsigned int * const line,
    const double * const p_free,
    const double * const p_not_measured,
    const double * const width,
    unsigned int num_cells,
    bool information,
    unsigned int dimension,
    double * const output) {

  local l;
  expected exp;

  // Set everything to zero
  exp.information3 = 0;
  exp.information2 = 0;
  exp.information1 = 0;
  exp.distance2 = 0;
  exp.distance1 = 0;
  exp.p_not_measured = 1;

  // Iterate backwards over the cells in the line
  for (int i = num_cells - 1; i >= 0; i--) {

    // Pre-compute
    unsigned int j = line[i];
    update_local(p_free[j], p_not_measured[j], width[i], l);

    // Update the expected information
    if (information) {
      switch (dimension) {
        case 3:
          exp.information3 = hit_or_miss(l, exp, information3_hit, information3_miss);
          [[fallthrough]];
        case 2:
          exp.information2 = hit_or_miss(l, exp, information2_hit, information2_miss);
          [[fallthrough]];
        case 1:
          exp.information1 = hit_or_miss(l, exp, information1_hit, information1_miss);
      }
    }

    // Update the expected distance
    switch (dimension) {
      case 3:
        exp.distance2 = hit_or_miss(l, exp, distance2_hit, distance2_miss);
        [[fallthrough]];
      case 2:
        exp.distance1 = hit_or_miss(l, exp, distance1_hit, distance1_miss);
    }

    // Update the average p_not_measured
    exp.p_not_measured = hit_or_miss(l, exp, p_not_measured_hit, p_not_measured_miss);

    // Update the output
    if (information) {
      switch (dimension) {
        case 3:
          output[j] += exp.information3;
          break;
        case 2:
          output[j] += exp.information2;
          break;
        case 1:
          output[j] += exp.information1;
          break;
      }
    } else {
      switch (dimension) {
        case 3:
          output[j] += exp.distance2;
          break;
        case 2:
          output[j] += exp.distance1;
          break;
      }
    }
  }
}

void range_entropy::expected::p_not_measured(
    const unsigned int * const line,
    const double * const vacancy,
    const double * const width,
    unsigned int num_cells,
    unsigned int dimension,
    double * const p_not_measured_) {

  double width_sum = 0;
  double miss_p_product = 1;
  double p_not_measured1_ = 0;
  for (unsigned int i = 0; i < num_cells; i++) {
    unsigned int j = line[i];
    double miss_p = std::pow(vacancy[j], width[i]);

    // Scale the probability by r^d to account
    // for radial overlap as well as the width
    // to account for aliasing.
    p_not_measured_[j] +=
      width[i] *
      std::pow(width_sum, dimension - 1) *
      p_not_measured1_;

    // The probability of the next cell not being
    // measured is the probability some previous
    // combination of cells were missed followed
    // by a hit.
    //
    // Also equal to the integral of the PDF of
    // a range measurement from 0 to the cell.
    p_not_measured1_ += miss_p_product * (1 - miss_p);
    width_sum += width[i];
    miss_p_product *= miss_p;
  }
}
