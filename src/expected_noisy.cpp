#include <cmath>
#include <algorithm>
#include <iostream>

#include "range_entropy/expected.hpp"
#include "range_entropy/expected_noisy.hpp"

double range_entropy::expected_noisy::normal_cdf(
    double x,
    double mean,
    double std_dev) {
  return 0.5 * (1 + std::erf((x - mean)/(std_dev * std::sqrt(2))));
}

void range_entropy::expected_noisy::pdf(
    const unsigned int * const line,
    const double * const vacancy,
    const double * const width,
    unsigned int num_cells,
    double std_dev,
    double truncation_width,
    double step_size,
    double pdf_width,
    double * const pdf,
    unsigned int & pdf_size) {

  // Zero the probability density function. 
  pdf_size = 0;
  for (
      double r = -truncation_width;
      r < pdf_width + truncation_width;
      r+= step_size) {
    pdf[pdf_size++] = 0;
  }

  double variance = std_dev * std_dev;

  double pdf_decay = 1;
  double width_sum = 0;

  // Iterate over the cells
  for (unsigned int line_cell = 0; line_cell < num_cells; line_cell++) {
    // Stop computing when the density function
    // of a cell no longer overlaps with
    // [-truncation_width, width[0] + truncation_width]
    if (width_sum > pdf_width + 2 * truncation_width) break;

    // Pre-compute
    unsigned int cell = line[line_cell];
    double v = vacancy[cell];
    double w = width[line_cell];
    double neg_log_v = -std::log(v);

    // The center of the normal
    // distribution is ahead of zero
    double normal_center = variance * neg_log_v;

    double normalization_constant =
      std::exp(0.5 * variance * neg_log_v * neg_log_v);

    unsigned int i_min = std::floor(width_sum/step_size);
    unsigned int i_max = std::ceil((width_sum + w + 2 * truncation_width)/step_size);
    for (unsigned int i = i_min; i < std::min(i_max, pdf_size); i++) {

      // Convert to a range measurement
      double r = i * step_size - width_sum - truncation_width;

      // A box with blurred edges is the
      // difference between normal CDFs
      double normal_scaling =
        normal_cdf(r, normal_center, std_dev) -
        normal_cdf(r, normal_center + w, std_dev);

      // Put it all together
      double pdf_noiseless = std::pow(v, r) * neg_log_v;

      double value = pdf_decay * pdf_noiseless * normal_scaling * normalization_constant;
      pdf[i] += value;
    }

    // Update width and decay
    width_sum += w;
    pdf_decay *= std::pow(v, w);
  }
}

double range_entropy::expected_noisy::hit(
    const double * const pdf,
    unsigned int pdf_size,
    double step_size,
    double noise_width,
    double (*value)(double, double)) {

  double integral = 0;
  for (unsigned int i = 0; i < pdf_size; i++) {
    double r = step_size * i - noise_width;
    integral += step_size * pdf[i] * value(r, pdf[i]);
  }

  return integral;
}

double range_entropy::expected_noisy::miss(
    const range_entropy::expected::local & l,
    const double * const pdf,
    unsigned int pdf_size,
    double step_size,
    double noise_width,
    double (*value)(double, double)) {

  double integral = 0;
  for (unsigned int i = 0; i < pdf_size; i++) {
    double r = step_size * i - noise_width + l.width;
    integral += step_size * pdf[i] * value(r, pdf[i] * l.miss_p);
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
    double noise_std_dev,
    double noise_width,
    double step_size,
    bool information,
    unsigned int dimension,
    double * const * const pdfs,
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

  double * const hit_pdf = pdfs[0];
  double * const miss_pdf = pdfs[1];

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
        noise_std_dev,
        noise_width,
        step_size,
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
        noise_std_dev,
        noise_width,
        step_size,
        0,
        miss_pdf,
        miss_pdf_size);

    double hit_integral, miss_integral;
    
    // Update the expected information
    if (information) {
      switch (dimension) {
        case 3:
          hit_integral = hit(hit_pdf, hit_pdf_size, step_size, noise_width,
              range_entropy::expected::information3);
          miss_integral = miss(l, miss_pdf, miss_pdf_size, step_size, noise_width,
              range_entropy::expected::information3);
          exp.information3 = hit_or_miss(l, exp, hit_integral, miss_integral,
              range_entropy::expected::information3_miss);
          [[fallthrough]];
        case 2:
          hit_integral = hit(hit_pdf, hit_pdf_size, step_size, noise_width,
              range_entropy::expected::information2);
          miss_integral = miss(l, miss_pdf, miss_pdf_size, step_size, noise_width,
              range_entropy::expected::information2);
          exp.information2 = hit_or_miss(l, exp, hit_integral, miss_integral,
              range_entropy::expected::information2_miss);
          [[fallthrough]];
        case 1:
          hit_integral = hit(hit_pdf, hit_pdf_size, step_size, noise_width,
              range_entropy::expected::information1);
          miss_integral = miss(l, miss_pdf, miss_pdf_size, step_size, noise_width,
              range_entropy::expected::information1);
          exp.information1 = hit_or_miss(l, exp, hit_integral, miss_integral,
              range_entropy::expected::information1_miss);
      }
    }

    // Update the expected distance
    switch (dimension) {
      case 3:
        hit_integral = hit(hit_pdf, hit_pdf_size, step_size, noise_width,
            range_entropy::expected::distance2);
        miss_integral = miss(l, miss_pdf, miss_pdf_size, step_size, noise_width,
            range_entropy::expected::distance2);
        exp.distance2 = hit_or_miss(l, exp, hit_integral, miss_integral,
            range_entropy::expected::distance2_miss);
        [[fallthrough]];
      case 2:
        hit_integral = hit(hit_pdf, hit_pdf_size, step_size, noise_width,
            range_entropy::expected::distance1);
        miss_integral = miss(l, miss_pdf, miss_pdf_size, step_size, noise_width,
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
