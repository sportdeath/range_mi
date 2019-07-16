#include <vector>
#include <iostream>
#include <cmath>
#include <random>

#include <range_entropy/expected.hpp>

using namespace range_entropy;

// Define constants
double integration_step = 0.00001;
unsigned int num_cells = 100;

double numerical_expected(
    const unsigned int * const line,
    const double * const p_free,
    const double * const p_not_measured,
    const double * const width,
    unsigned int num_cells,
    double (*value)(double,double)) {

  unsigned int i = 0;
  double r = 0;
  double width_sum = 0;
  double pdf_decay = 1;
  double cdf = 0;
  double integral = 0;
  while (i < num_cells) {
    unsigned int j = line[i];

    // Compute the pdf
    double pdf = 
      pdf_decay * (
          -std::pow(p_free[j], r - width_sum) *
          std::log(p_free[j]));

    // Accumulate
    integral +=
      p_not_measured[j] * 
      pdf * value(r, pdf) * integration_step;

    // Integrate the cdf
    cdf += pdf * integration_step;

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

  // Ignore probability of miss
  // because it is very low

  return integral;
}

double distance1(double r, double pdf) {
  return r;
}
double distance2(double r, double pdf) {
  return r * r;
}
double information1(double r, double pdf) {
  return -std::log(pdf);
}
double information2(double r, double pdf) {
  return -r * std::log(pdf);
}
double information3(double r, double pdf) {
  return -r * r * std::log(pdf);
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

  // Compute the values numerically
  double numerical_expected_distance1 =
    numerical_expected(
        line.data(), p_free.data(), p_not_measured.data(), width.data(), num_cells,
        distance1);
  double numerical_expected_distance2 =
    numerical_expected(
        line.data(), p_free.data(), p_not_measured.data(), width.data(), num_cells,
        distance2);
  double numerical_expected_information1 =
    numerical_expected(
        line.data(), p_free.data(), p_not_measured.data(), width.data(), num_cells,
        information1);
  double numerical_expected_information2 =
    numerical_expected(
        line.data(), p_free.data(), p_not_measured.data(), width.data(), num_cells,
        information2);
  double numerical_expected_information3 =
    numerical_expected(
        line.data(), p_free.data(), p_not_measured.data(), width.data(), num_cells,
        information3);

  // Compute the values exactly
  std::vector<double> expected_distance1(num_cells, 0),
                      expected_distance2(num_cells, 0),
                      expected_information1(num_cells, 0),
                      expected_information2(num_cells, 0),
                      expected_information3(num_cells, 0);
  expected::distance1(
      line.data(), p_free.data(), p_not_measured.data(), width.data(), num_cells,
      expected_distance1.data());
  expected::distance2(
      line.data(), p_free.data(), p_not_measured.data(), width.data(), num_cells,
      expected_distance2.data());
  expected::information1(
      line.data(), p_free.data(), p_not_measured.data(), width.data(), num_cells,
      expected_information1.data());
  expected::information2(
      line.data(), p_free.data(), p_not_measured.data(), width.data(), num_cells,
      expected_information2.data());
  expected::information3(
      line.data(), p_free.data(), p_not_measured.data(), width.data(), num_cells,
      expected_information3.data());

  std::cout << "d1: " << numerical_expected_distance1 << ", " << expected_distance1[0] << std::endl;
  std::cout << "d2: " << numerical_expected_distance2 << ", " << expected_distance2[0] << std::endl;
  std::cout << "i1: " << numerical_expected_information1 << ", " << expected_information1[0] << std::endl;
  std::cout << "i2: " << numerical_expected_information2 << ", " << expected_information2[0] << std::endl;
  std::cout << "i3: " << numerical_expected_information3 << ", " << expected_information3[0] << std::endl;
}
