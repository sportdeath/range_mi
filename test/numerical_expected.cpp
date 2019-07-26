#include <vector>
#include <iostream>
#include <cmath>
#include <random>

#include <range_entropy/expected.hpp>
#include <range_entropy/expected_noisy.hpp>

using namespace range_entropy;

// Define constants
double integration_step = 0.01;
double noise_dev = 2;
double noise_width = 4 * noise_dev;
unsigned int num_cells = 100;

void numerical_pdf(
    const unsigned int * const line,
    const double * const p_free,
    const double * const p_not_measured,
    const double * const width,
    unsigned int num_cells,
    double * const pdf,
    double * const p_not_measured_stepped,
    unsigned int & pdf_size) {

  unsigned int i = 0;
  double r = 0;
  double width_sum = 0;
  double pdf_decay = 1;
  pdf_size = 0;
  while (i < num_cells) {
    unsigned int j = line[i];

    // Compute the pdf
    p_not_measured_stepped[pdf_size] = p_not_measured[j];
    pdf[pdf_size++] = 
      pdf_decay * (
          -std::pow(p_free[j], r - width_sum) *
          std::log(p_free[j]));

    // Make a step
    r += integration_step;
    if (r - width_sum > width[i]) {
      // Cell completed,
      // move to the next cell!
      pdf_decay *= std::pow(p_free[j], width[i]);
      width_sum += width[i];
      i++;
    }
  }
}

double numerical_expected(
    const double * const pdf,
    const double * const p_not_measured_stepped,
    unsigned pdf_size,
    double (*value)(double,double)) {

  double integral = 0;
  for (unsigned int i = 0; i < pdf_size; i++) {
    double r = i * integration_step;

    // Accumulate
    integral +=
      p_not_measured_stepped[i] * 
      pdf[i] * value(r, pdf[i]) * integration_step;
  }

  return integral;
}

double numerical_expected_noisy(
    const double * const pdf,
    unsigned pdf_size,
    double * const convolve,
    double (*value)(double,double)) {

  // Clear the convolution
  unsigned int convolve_size = pdf_size + (2 * noise_width)/integration_step + 1;
  for (unsigned int i = 0; i < convolve_size; i++) {
    convolve[i] = 0;
  }

  // Perform the convolution
  for (unsigned int i = 0; i < pdf_size; i++) {
    double r = i * integration_step;

    for (double n = -noise_width; n < noise_width; n+=integration_step) {
      unsigned int j = (r + n + noise_width)/integration_step;

      double noise = 1/std::sqrt(2 * M_PI * noise_dev * noise_dev) * std::exp(-(n * n)/(2 * noise_dev * noise_dev));

      convolve[j] += pdf[i] * noise * integration_step;
    }
  }
  
  // Perform the integration
  double integral = 0;
  for (unsigned int i = 0; i < convolve_size; i++) {
    double z = i * integration_step - noise_width;
    integral +=
      convolve[i] * value(z, convolve[i]) * integration_step;
  }

  return integral;
}

// Initialize random generator
std::random_device random_device;
std::mt19937 gen(random_device());
std::uniform_real_distribution<double> dist(0.,1.);

// ... and a way to make random vectors
void random_p(std::vector<double> & p) {
  for (size_t i = 0; i < p.size(); i++) {
    p[i] = dist(gen);
  }
}

int main() {
  // Initialize the inputs and outputs
  std::vector<double> p_free(num_cells);
  std::vector<double> p_not_measured(num_cells);
  std::vector<double> width(num_cells);
  // Randomize then
  random_p(p_free);
  random_p(p_not_measured);
  random_p(width);

  // Make a line
  std::vector<unsigned int> line(num_cells);
  std::iota(line.begin(), line.end(), 0);

  // Compute the pdf
  unsigned int pdf_size;
  std::vector<double> pdf(num_cells/integration_step);
  std::vector<double> p_not_measured_stepped(num_cells/integration_step);
  std::vector<double> convolution((num_cells + 2 * noise_width)/integration_step);
  numerical_pdf(
      line.data(), p_free.data(), p_not_measured.data(), width.data(), num_cells,
      pdf.data(), p_not_measured_stepped.data(), pdf_size);

  std::cout << "Computing noiseless information numerically and recursively..." << std::endl;

  // Compute the values numerically
  double numerical_expected_distance1 =
    numerical_expected(
        pdf.data(), p_not_measured_stepped.data(), pdf_size,
        range_entropy::expected::distance1);
  double numerical_expected_distance2 =
    numerical_expected(
        pdf.data(), p_not_measured_stepped.data(), pdf_size,
        range_entropy::expected::distance2);
  double numerical_expected_information1 =
    numerical_expected(
        pdf.data(), p_not_measured_stepped.data(), pdf_size,
        range_entropy::expected::information1);
  double numerical_expected_information2 =
    numerical_expected(
        pdf.data(), p_not_measured_stepped.data(), pdf_size,
        range_entropy::expected::information2);
  double numerical_expected_information3 =
    numerical_expected(
        pdf.data(), p_not_measured_stepped.data(), pdf_size,
        range_entropy::expected::information3);

  // Compute the values exactly
  std::vector<double> expected_distance1(num_cells, 0),
                      expected_distance2(num_cells, 0),
                      expected_information1(num_cells, 0),
                      expected_information2(num_cells, 0),
                      expected_information3(num_cells, 0);
  expected::line(
      line.data(), p_free.data(), p_not_measured.data(), width.data(), num_cells,
      false, 2,
      expected_distance1.data());
  expected::line(
      line.data(), p_free.data(), p_not_measured.data(), width.data(), num_cells,
      false, 3,
      expected_distance2.data());
  expected::line(
      line.data(), p_free.data(), p_not_measured.data(), width.data(), num_cells,
      true, 1,
      expected_information1.data());
  expected::line(
      line.data(), p_free.data(), p_not_measured.data(), width.data(), num_cells,
      true, 2,
      expected_information2.data());
  expected::line(
      line.data(), p_free.data(), p_not_measured.data(), width.data(), num_cells,
      true, 3,
      expected_information3.data());

  std::cout << "d1: " << numerical_expected_distance1 << ", " << expected_distance1[0] << std::endl;
  std::cout << "d2: " << numerical_expected_distance2 << ", " << expected_distance2[0] << std::endl;
  std::cout << "i1: " << numerical_expected_information1 << ", " << expected_information1[0] << std::endl;
  std::cout << "i2: " << numerical_expected_information2 << ", " << expected_information2[0] << std::endl;
  std::cout << "i3: " << numerical_expected_information3 << ", " << expected_information3[0] << std::endl;

  std::cout << std::endl;
  std::cout << "Computing noisy information numerically and recursively..." << std::endl;

  // Compute the values numerically
  double numerical_expected_noisy_distance1 =
    numerical_expected_noisy(
        pdf.data(), pdf_size, convolution.data(),
        range_entropy::expected::distance1);
  double numerical_expected_noisy_distance2 =
    numerical_expected_noisy(
        pdf.data(), pdf_size, convolution.data(),
        range_entropy::expected::distance2);
  double numerical_expected_noisy_information1 =
    numerical_expected_noisy(
        pdf.data(), pdf_size, convolution.data(),
        range_entropy::expected::information1);
  double numerical_expected_noisy_information2 =
    numerical_expected_noisy(
        pdf.data(), pdf_size, convolution.data(),
        range_entropy::expected::information2);
  double numerical_expected_noisy_information3 =
    numerical_expected_noisy(
        pdf.data(), pdf_size, convolution.data(),
        range_entropy::expected::information3);

  // Compute the values via noisy computation
  std::vector<double> expected_noisy_distance1(num_cells, 0),
                      expected_noisy_distance2(num_cells, 0),
                      expected_noisy_information1(num_cells, 0),
                      expected_noisy_information2(num_cells, 0),
                      expected_noisy_information3(num_cells, 0);
  std::vector<double> hit_pdf(num_cells/integration_step),
                      miss_pdf(num_cells/integration_step);
  expected_noisy::line(
      line.data(), p_free.data(), width.data(), num_cells,
      noise_dev, noise_width, integration_step,
      false, 2,
      hit_pdf.data(), miss_pdf.data(),
      expected_noisy_distance1.data());
  expected_noisy::line(
      line.data(), p_free.data(), width.data(), num_cells,
      noise_dev, noise_width, integration_step,
      false, 3,
      hit_pdf.data(), miss_pdf.data(),
      expected_noisy_distance2.data());
  expected_noisy::line(
      line.data(), p_free.data(), width.data(), num_cells,
      noise_dev, noise_width, integration_step,
      true, 1,
      hit_pdf.data(), miss_pdf.data(),
      expected_noisy_information1.data());
  expected_noisy::line(
      line.data(), p_free.data(), width.data(), num_cells,
      noise_dev, noise_width, integration_step,
      true, 2,
      hit_pdf.data(), miss_pdf.data(),
      expected_noisy_information2.data());
  expected_noisy::line(
      line.data(), p_free.data(), width.data(), num_cells,
      noise_dev, noise_width, integration_step,
      true, 3,
      hit_pdf.data(), miss_pdf.data(),
      expected_noisy_information3.data());

  std::cout << "d1: " << numerical_expected_noisy_distance1 << ", " << expected_noisy_distance1[0] << std::endl;
  std::cout << "d2: " << numerical_expected_noisy_distance2 << ", " << expected_noisy_distance2[0] << std::endl;
  std::cout << "i1: " << numerical_expected_noisy_information1 << ", " << expected_noisy_information1[0] << std::endl;
  std::cout << "i2: " << numerical_expected_noisy_information2 << ", " << expected_noisy_information2[0] << std::endl;
  std::cout << "i3: " << numerical_expected_noisy_information3 << ", " << expected_noisy_information3[0] << std::endl;

}
