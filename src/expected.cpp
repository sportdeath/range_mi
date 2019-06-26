#include <cmath>

#include "range_entropy/expected.hpp"

using namespace range_entropy;
using namespace expected;

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
  l.miss_p = std::pow(p_free, width);

  // The converse of p_miss
  l.hit_p = 1 - l.miss_p;

  // The amount of information learned if
  // a miss occurs
  //
  // -log(p_miss) = -w log p
  // 
  l.miss_info = - width * std::log(p_free);

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
  l.zero_info = -std::log(-std::log(p_free));
}

double range_entropy::expected::distance1(
    const local & l,
    double distance1_) {

  // If the range misses the local cell,
  // the expected distance is simply the old
  // expected distance plus the width of the cell
  double miss_distance1 = distance1_ + l.width;

  // If the range hits the cell the
  // expected distance can be evaluated analytically
  double hit_p_hit_distance1 = l.width * (l.hit_p * l.miss_info_inv - l.miss_p);

  // Return the expected value
  return hit_p_hit_distance1 + l.miss_p * miss_distance1;
}

double range_entropy::expected::distance2(
    const local & l,
    double distance1_,
    double distance2_) {

  // If the range misses the local cell,
  // the expected distance squared is as follows
  double miss_distance2 =
    distance2_ + 2 * l.width * distance1_ + l.width * l.width;

  double hit_p_hit_distance2 = l.width * l.width * (
      2 * l.miss_info_inv * l.miss_info_inv * l.hit_p - l.miss_p * (l.miss_info + 2) * l.miss_info_inv);
  
  // Return the expected value
  return hit_p_hit_distance2 + l.miss_p * miss_distance2;
}

double range_entropy::expected::information1(
    const local & l,
    double information1_) {
  // If the range measurement misses
  // we gain all the information from the far
  // measurement plus the information gained
  // by the fact that we missed.
  double miss_info = information1_ + l.miss_info;

  // If the measurement hits then in expectation
  // we gain the following information
  double hit_p_hit_info = l.hit_p * (l.zero_info + 1) - l.miss_p * l.miss_info;

  return hit_p_hit_info + l.miss_p * miss_info;
}

double range_entropy::expected::information2(
    const local & l,
    double distance1_,
    double information1_,
    double information2_) {

  double miss_info =
    information2_ + l.miss_info * distance1_ +
    l.width * (information1_ + l.miss_info);

  double hit_p_hit_info = l.width * (
      l.hit_p * (2 + l.zero_info) * l.miss_info_inv - l.miss_p * (2 + l.miss_info + l.zero_info));

  return hit_p_hit_info + l.miss_p * miss_info;
}

double range_entropy::expected::information3(
    const local & l,
    double distance1_,
    double distance2_,
    double information1_,
    double information2_,
    double information3_) {

  double miss_info =
    information3_ + l.miss_info * distance2_ +
    2 * l.width * (information2_ + l.miss_info * distance1_) +
    l.width * l.width * (information1_ + l.miss_info);

  double hit_p_hit_info = l.width * l.width * (
      l.hit_p * (6 + 2 * l.zero_info) * l.miss_info_inv * l.miss_info_inv -
      l.miss_p * l.miss_info_inv * (6 + l.miss_info + (l.miss_info + l.zero_info) * (l.miss_info + 2)));

  return hit_p_hit_info + l.miss_p * miss_info;
}

void range_entropy::expected::distance1(
    const unsigned int * const line,
    const double * const p_free,
    const double * const p_not_measured,
    const double * const width,
    unsigned int num_cells,
    double * const distance1_) {

  local l;

  double d1 = 0;
  for (int i = num_cells - 1; i >= 0; i--) {
    unsigned int j = line[i];
    update_local(p_free[j], p_not_measured[j], width[i], l);
    d1 = distance1(l, d1);
    distance1_[j] += d1;
  }
}

void range_entropy::expected::distance2(
    const unsigned int * const line,
    const double * const p_free,
    const double * const p_not_measured,
    const double * const width,
    unsigned int num_cells,
    double * const distance2_) {

  local l;

  double d1 = 0, d2 = 0;
  for (int i = num_cells - 1; i >= 0; i--) {
    unsigned int j = line[i];
    update_local(p_free[j], p_not_measured[j], width[i], l);
    d2 = distance2(l, d1, d2);
    d1 = distance1(l, d1);
    distance2_[j] += d2;
  }
}

void range_entropy::expected::information1(
    const unsigned int * const line,
    const double * const p_free,
    const double * const p_not_measured,
    const double * const width,
    unsigned int num_cells,
    double * const information1_) {

  local l;

  double i1 = 0;
  for (int i = num_cells - 1; i >= 0; i--) {
    unsigned int j = line[i];
    update_local(p_free[j], p_not_measured[j], width[i], l);
    i1 = information1(l, i1);
    information1_[j] += i1;
  }
}

void range_entropy::expected::information2(
    const unsigned int * const line,
    const double * const p_free,
    const double * const p_not_measured,
    const double * const width,
    unsigned int num_cells,
    double * const information2_) {

  local l;

  double d1 = 0, i1 = 0, i2 = 0;
  for (int i = num_cells - 1; i >= 0; i--) {
    unsigned int j = line[i];
    update_local(p_free[j], p_not_measured[j], width[i], l);
    i2 = information2(l, d1, i1, i2);
    i1 = information1(l, i1);
    d1 = distance1(l, d1);
    information2_[j] += i2;
  }
}

void range_entropy::expected::information3(
    const unsigned int * const line,
    const double * const p_free,
    const double * const p_not_measured,
    const double * const width,
    unsigned int num_cells,
    double * const information3_) {

  local l;

  double d1 = 0, d2 = 0, i1 = 0, i2 = 0, i3 = 0;
  for (int i = num_cells - 1; i >= 0; i--) {
    unsigned int j = line[i];
    update_local(p_free[j], p_not_measured[j], width[i], l);
    i3 = information3(l, d1, d2, i1, i2, i3);
    i2 = information2(l, d1, i1, i2);
    i1 = information1(l, i1);
    d2 = distance2(l, d1, d2);
    d1 = distance1(l, d1);
    information3_[j] += i3;
  }
}
