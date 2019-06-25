#include <cmath>
#include <vector>

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

  // The odds of missing
  l.miss_odds = l.miss_p/l.hit_p;

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
  double hit_distance1 = l.width * (l.miss_info_inv - l.miss_odds);

  // Return the expected value
  return l.hit_p * hit_distance1 + l.miss_p * miss_distance1;
}

double range_entropy::expected::distance2(
    const local & l,
    double distance1_,
    double distance2_) {

  // If the range misses the local cell,
  // the expected distance squared is as follows
  double miss_distance2 =
    distance2_ + 2 * l.width * distance1_ + l.width * l.width;

  double hit_distance2 = l.width * l.width * (
      2 * l.miss_info_inv * l.miss_info_inv - l.miss_odds * (l.miss_info + 2) * l.miss_info_inv);
  
  // Return the expected value
  return l.hit_p * hit_distance2 + l.miss_p * miss_distance2;
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
  double hit_info = l.zero_info + 1 - l.miss_odds * l.miss_info;

  return l.hit_p * hit_info + l.miss_p * miss_info;
}

std::vector<double> range_entropy::expected::distance1(
    const std::vector<double> & p_free,
    const std::vector<double> & p_not_measured,
    const std::vector<double> & width) {

  std::vector<double> out(p_free.size());
  local l;

  double d1 = 0;
  for (int i = p_free.size() - 1; i >= 0; i--) {
    update_local(p_free[i], p_not_measured[i], width[i], l);
    d1 = distance1(l, d1);
    out[i] = d1;
  }

  return out;
}

std::vector<double> range_entropy::expected::distance2(
    const std::vector<double> & p_free,
    const std::vector<double> & p_not_measured,
    const std::vector<double> & width) {

  std::vector<double> out(p_free.size());
  local l;

  double d1 = 0, d2 = 0;
  for (int i = p_free.size() - 1; i >= 0; i--) {
    update_local(p_free[i], p_not_measured[i], width[i], l);
    d2 = distance2(l, d1, d2);
    d1 = distance1(l, d1);
    out[i] = d2;
  }

  return out;
}

std::vector<double> range_entropy::expected::information1(
    const std::vector<double> & p_free,
    const std::vector<double> & p_not_measured,
    const std::vector<double> & width) {

  std::vector<double> out(p_free.size());
  local l;

  double i1 = 0;
  for (int i = p_free.size() - 1; i >= 0; i--) {
    update_local(p_free[i], p_not_measured[i], width[i], l);
    i1 = information1(l, i1);
    out[i] = i1;
  }

  return out;
}
