#pragma once

#include "range_entropy/expected.hpp"

namespace range_entropy {
namespace expected_noisy {

double normal_pdf(double x, double dev);

/**
 * Calculate the probability density of a range
 * measurement with additive noise. The noise
 * is normally distributed with deviation noise_dev,
 * and mean zero.
 *
 * For computational purposes the noise is truncated
 * to zero beyond the noise_width. The PDF is computed
 * at an intervals of step_size.
 *
 * The pdf is calculated from -noise_width to 
 * pdf_width + noise_width. Zero is the start of the
 * first cell.
 */
void pdf(
    const unsigned int * const line,
    const double * const vacancy,
    const double * const width,
    unsigned int num_cells,
    const double * const noise,
    unsigned int noise_size,
    double integration_step,
    double pdf_width,
    double * const pdf,
    unsigned int & pdf_size);

/**
 * Compute the expected value
 * in the hit and close miss regions.
 */
double hit(
    const double * const pdf,
    unsigned int pdf_size,
    unsigned int noise_size,
    double integration_step,
    double (*value)(double, double));
double miss(
    const range_entropy::expected::local & l,
    const double * const pdf,
    unsigned int pdf_size,
    unsigned int noise_size,
    double integration_step,
    double (*value)(double, double));

/**
 * Compute the expected value over the
 * entire beam.
 */
double hit_or_miss(
    const range_entropy::expected::local & l,
    const range_entropy::expected::expected & exp,
    double hit_integral,
    double miss_integral,
    double (*miss)(
      const range_entropy::expected::local &,
      const range_entropy::expected::expected &));

void line(
    const unsigned int * const line,
    const double * const vacancy,
    const double * const width,
    unsigned int num_cells,
    const double * const noise,
    unsigned int noise_size,
    double integration_step,
    bool information,
    unsigned int dimension,
    double * const hit_pdf,
    double * const miss_pdf,
    double * const output);

}}
