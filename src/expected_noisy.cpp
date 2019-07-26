#include <cmath>
#include <algorithm>
#include <iostream>

#include "range_entropy/expected.hpp"
#include "range_entropy/expected_noisy.hpp"

double range_entropy::expected_noisy::normal_pdf(
    double x,
    double dev) {
  double x_scaled = x/dev;
  return 1./(std::sqrt(2 * M_PI) * dev) * std::exp(-(x_scaled * x_scaled)/2.);
}

void range_entropy::expected_noisy::pdf(
    const unsigned int * const line,
    const double * const vacancy,
    const double * const width,
    unsigned int num_cells,
    const double * const noise,
    unsigned int noise_size,
    double integration_step,
    double pdf_width,
    double * const pdf,
    unsigned int & pdf_size) {

  double noise_width = integration_step * noise_size;
  double noise_half_width = noise_width/2.;

  double pdf_decay = 1;
  double width_sum = 0;
  double r = 0;
  unsigned int r_index = 0;
  pdf_size = 0;

  // Iterate over the cells
  for (unsigned int line_cell = 0; line_cell < num_cells; line_cell++) {
    // If the width is too large the cells are of no use
    if (width_sum > pdf_width + noise_width) break;

    // Fetch vacancy and width
    unsigned int cell = line[line_cell];
    double v = vacancy[cell];
    double w = width[line_cell];
    double neg_log_v = -std::log(v);

    // Until we hit the end of the cell
    while (r < width_sum + w) {

      // Compute the value of the PDF
      double pdf_value = pdf_decay * std::pow(v, r) * neg_log_v;

      // And convolve with a normal distribution
      double n = -noise_half_width;
      for (unsigned int n_index = 0; n_index < noise_size; n_index++) {
        // This is the convolved position
        double z = r + n;

        // Don't compute past necessary
        if (z > pdf_width + noise_half_width) break;

        // If this is the first time we're visiting
        // the cell, clear it.
        if (r_index + n_index >= pdf_size) {
          pdf[r_index + n_index] = 0;
          pdf_size++;
        }

        // Accumulated to convolute
        pdf[r_index + n_index] += integration_step * noise[n_index] * pdf_value;

        // Make a step
        n += integration_step;
      }

      // Make a step
      r += integration_step;
      r_index++;
    }

    // Move to the next cell
    width_sum += w;
    pdf_decay *= std::pow(v, w);
  }
}

double range_entropy::expected_noisy::hit(
    const double * const pdf,
    unsigned int pdf_size,
    unsigned int noise_size,
    double integration_step,
    double (*value)(double, double)) {

  double integral = 0;
  double r = -(integration_step * noise_size)/2.;
  for (unsigned int i = 0; i < pdf_size; i++) {
    integral += integration_step * pdf[i] * value(r, pdf[i]);
    r += integration_step;
  }

  return integral;
}

double range_entropy::expected_noisy::miss(
    const range_entropy::expected::local & l,
    const double * const pdf,
    unsigned int pdf_size,
    unsigned int noise_size,
    double integration_step,
    double (*value)(double, double)) {

  double integral = 0;
  double r = -(integration_step * noise_size)/2. + l.width;
  for (unsigned int i = 0; i < pdf_size; i++) {
    integral += integration_step * pdf[i] * value(r, pdf[i] * l.miss_p);
    r += integration_step;
  }

  return integral;
}

double range_entropy::expected_noisy::hit_or_miss(
    const range_entropy::expected::local & l,
    const range_entropy::expected::expected & exp,
    double hit_integral,
    double miss_integral,
    double (*miss)(
      const range_entropy::expected::local &,
      const range_entropy::expected::expected &)) {

  return hit_integral + l.miss_p * (miss(l, exp) - miss_integral);
}

void range_entropy::expected_noisy::line(
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
    double * const output) {

  range_entropy::expected::local l;
  range_entropy::expected::expected exp;

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
    update_local(vacancy[j], 1, width[i], l);

    // Compute the PDF along
    // [-noise_width, width[i] + noise_width]
    unsigned int hit_pdf_size;
    range_entropy::expected_noisy::pdf(
        line + i,
        vacancy,
        width + i,
        num_cells - i,
        noise,
        noise_size,
        integration_step,
        width[i],
        hit_pdf,
        hit_pdf_size);
    // Compute the PDF along
    // [-noise_width, noise_width]
    // of the previous cell
    unsigned int miss_pdf_size;
    range_entropy::expected_noisy::pdf(
        line + i + 1,
        vacancy,
        width + i + 1,
        num_cells - i - 1,
        noise,
        noise_size,
        integration_step,
        0,
        miss_pdf,
        miss_pdf_size);

    double hit_integral, miss_integral;
    
    // Update the expected information
    if (information) {
      switch (dimension) {
        case 3:
          hit_integral = hit(hit_pdf, hit_pdf_size, noise_size, integration_step,
              range_entropy::expected::information3);
          miss_integral = miss(l, miss_pdf, miss_pdf_size, noise_size, integration_step,
              range_entropy::expected::information3);
          exp.information3 = hit_or_miss(l, exp, hit_integral, miss_integral,
              range_entropy::expected::information3_miss);
          [[fallthrough]];
        case 2:
          hit_integral = hit(hit_pdf, hit_pdf_size, noise_size, integration_step,
              range_entropy::expected::information2);
          miss_integral = miss(l, miss_pdf, miss_pdf_size, noise_size, integration_step,
              range_entropy::expected::information2);
          exp.information2 = hit_or_miss(l, exp, hit_integral, miss_integral,
              range_entropy::expected::information2_miss);
          [[fallthrough]];
        case 1:
          hit_integral = hit(hit_pdf, hit_pdf_size, noise_size, integration_step,
              range_entropy::expected::information1);
          miss_integral = miss(l, miss_pdf, miss_pdf_size, noise_size, integration_step,
              range_entropy::expected::information1);
          exp.information1 = hit_or_miss(l, exp, hit_integral, miss_integral,
              range_entropy::expected::information1_miss);
      }
    }

    // Update the expected distance
    switch (dimension) {
      case 3:
        hit_integral = hit(hit_pdf, hit_pdf_size, noise_size, integration_step,
            range_entropy::expected::distance2);
        miss_integral = miss(l, miss_pdf, miss_pdf_size, noise_size, integration_step,
            range_entropy::expected::distance2);
        exp.distance2 = hit_or_miss(l, exp, hit_integral, miss_integral,
            range_entropy::expected::distance2_miss);
        [[fallthrough]];
      case 2:
        hit_integral = hit(hit_pdf, hit_pdf_size, noise_size, integration_step,
            range_entropy::expected::distance1);
        miss_integral = miss(l, miss_pdf, miss_pdf_size, noise_size, integration_step,
            range_entropy::expected::distance1);
        exp.distance1 = hit_or_miss(l, exp, hit_integral, miss_integral,
            range_entropy::expected::distance1_miss);
    }

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
