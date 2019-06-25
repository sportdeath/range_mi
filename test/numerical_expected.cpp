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
    const std::vector<double> & p_free,
    const std::vector<double> & width,
    double (*value)(double,double)) {

  unsigned int i = 0;
  double r = 0;
  double width_sum = 0;
  double pdf_decay = 1;
  double cdf = 0;
  double integral = 0;
  while (i < p_free.size()) {
    // Compute the pdf
    double pdf = 
      pdf_decay * (
          -std::pow(p_free[i], r - width_sum) *
          std::log(p_free[i]));

    // Accumulate
    integral +=
      pdf * value(r, pdf) * integration_step;

    // Integrate the cdf
    cdf += pdf * integration_step;

    // Make a step
    r += integration_step;
    if (r - width_sum > width[i]) {
      // Cell completed,
      // move to the next cell!
      pdf_decay *= std::pow(p_free[i], width[i]);
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

  // Compute the values numerically
  double numerical_expected_distance1 =
    numerical_expected(
        p_free, width, distance1);
  double numerical_expected_distance2 =
    numerical_expected(
        p_free, width, distance2);
  double numerical_expected_information1 =
    numerical_expected(
        p_free, width, information1);
  double numerical_expected_information2 =
    numerical_expected(
        p_free, width, information2);
  double numerical_expected_information3 =
    numerical_expected(
        p_free, width, information3);

  // Compute the values exactly
  double expected_distance1 =
    expected::distance1(
        p_free, p_not_measured, width)[0];
  double expected_distance2 =
    expected::distance2(
        p_free, p_not_measured, width)[0];
  double expected_information1 =
    expected::information1(
        p_free, p_not_measured, width)[0];
  double expected_information2 =
    expected::information2(
        p_free, p_not_measured, width)[0];
  double expected_information3 =
    expected::information3(
        p_free, p_not_measured, width)[0];

  std::cout << "d1: " << numerical_expected_distance1 << ", " << expected_distance1 << std::endl;
  std::cout << "d2: " << numerical_expected_distance2 << ", " << expected_distance2 << std::endl;
  std::cout << "i1: " << numerical_expected_information1 << ", " << expected_information1 << std::endl;
  std::cout << "i2: " << numerical_expected_information2 << ", " << expected_information2 << std::endl;
  std::cout << "i3: " << numerical_expected_information3 << ", " << expected_information3 << std::endl;
}
