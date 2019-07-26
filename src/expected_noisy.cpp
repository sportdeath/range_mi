#include <cmath>
#include <algorithm>

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
    double noise_dev,
    double noise_half_width,
    double integration_step,
    double pdf_width,
    double * const pdf,
    unsigned int & pdf_size) {

  double variance = noise_dev * noise_dev;

  double pdf_decay = 1;
  double width_sum = 0;
  pdf_size = 0;

  // Iterate over the cells
  for (unsigned int line_cell = 0; line_cell < num_cells; line_cell++) {
    // If the width is too large the cells are of no use
    if (width_sum > pdf_width + 2 * noise_half_width) break;

    unsigned int cell = line[line_cell];
    double v = vacancy[cell];
    double w = width[line_cell];
    double neg_log_v = -std::log(v);

    // The center of the normal
    // distribution is ahead of zero
    double normal_center = variance * neg_log_v;

    double normalization_constant =
      std::exp(0.5 * variance * neg_log_v * neg_log_v);

    // Until we hit the end of the cell (plus noise)
    double z = -noise_half_width;
    unsigned int z_index = width_sum/integration_step;
    while (z < w + noise_half_width and z + width_sum < pdf_width + noise_half_width) {

      // Compute the value of the PDF
      double noiseless_value = pdf_decay * std::pow(v, z) * neg_log_v;

      // The CDF
      double normal_window =
        normal_cdf(z, normal_center, noise_dev) -
        normal_cdf(z, normal_center + w, noise_dev);

      // Compute the full value
      double value = normalization_constant * normal_window * noiseless_value;

      // If this is the first time the cell
      // has been touched, clear it
      if (z_index == pdf_size) {
        pdf[z_index] = 0;
        pdf_size++;
      }

      // Accumulate
      pdf[z_index] += value;

      z += integration_step;
      z_index++;
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
