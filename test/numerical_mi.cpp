#include <vector>
#include <iostream>
#include <cmath>
#include <random>
#include <fstream>

#include <range_mi/barely_distorted.hpp>

#include "helpers.hpp"

// Define constants
double integration_step = 0.000001;
double vacancy_scaling = 0.1;
double dtheta = 0.1;
unsigned int num_cells = 100;
const unsigned int num_dimensions = 5;

double numerical_mi(
    const double * const pdf,
    unsigned pdf_size,
    unsigned int dimension,
    double dtheta,
    double integration_step) {

  double h_Z = 0;
  double h_Z_given_M = 0;
  for (unsigned int i = 0; i < pdf_size; i++) {
    double r = i * integration_step;

    h_Z +=
      pdf[i] *
      std::pow(r, dimension - 1) *
      -std::log(pdf[i]) *
      integration_step *
      dtheta;

    h_Z_given_M +=
      (1 - std::log(range_mi::barely_distorted::noise_l)) *
      pdf[i] *
      std::pow(r, dimension - 1) *
      integration_step *
      dtheta;
  }

  return h_Z - h_Z_given_M;
}

int main() {
  // Initialize the inputs and outputs
  std::vector<double> vacancy(num_cells);
  std::vector<double> p_not_measured(num_cells, 1);
  std::vector<double> width(num_cells);
  // Randomize then
  random_p(vacancy);
  random_p(width);

  // Scale vacancy to make everything generally lower
  for (unsigned int i = 0; i < num_cells; i++) {
    vacancy[i] = std::pow(vacancy[i], vacancy_scaling);
  }

  // Make a line
  std::vector<unsigned int> line(num_cells);
  std::iota(line.begin(), line.end(), 0);

  std::cout << "Computing the barely distorted pdf" << std::endl;

  // Compute the pdf
  unsigned int pdf_size;
  std::vector<double> pdf(num_cells/integration_step);
  numerical_pdf(
      line.data(),
      vacancy.data(),
      width.data(),
      num_cells,
      integration_step,
      pdf.data(),
      pdf_size);

  std::cout << "Computing the barely distorted mutual information numerically..." << std::endl;

  std::vector<double> numerical_mi_(num_dimensions);
  for (unsigned int i = 0; i < num_dimensions; i++) {
    numerical_mi_[i] = numerical_mi(pdf.data(), pdf_size, i + 1, dtheta, integration_step);
  }

  std::cout << "Computing the barely distorted mutual information exactly..." << std::endl;

  std::vector<double> exact_mi(num_dimensions, 0);
  std::vector<double> mi(num_cells);

  // Clear the output
  std::fill(mi.begin(), mi.end(), 0);
  range_mi::barely_distorted::line<1>(
      line.data(), vacancy.data(), p_not_measured.data(),
      width.data(), num_cells, dtheta, mi.data());
  exact_mi[0] = mi[0];
  std::fill(mi.begin(), mi.end(), 0);
  range_mi::barely_distorted::line<2>(
      line.data(), vacancy.data(), p_not_measured.data(),
      width.data(), num_cells, dtheta, mi.data());
  exact_mi[1] = mi[0];
  std::fill(mi.begin(), mi.end(), 0);
  range_mi::barely_distorted::line<3>(
      line.data(), vacancy.data(), p_not_measured.data(),
      width.data(), num_cells, dtheta, mi.data());
  exact_mi[2] = mi[0];
  std::fill(mi.begin(), mi.end(), 0);
  range_mi::barely_distorted::line<4>(
      line.data(), vacancy.data(), p_not_measured.data(),
      width.data(), num_cells, dtheta, mi.data());
  exact_mi[3] = mi[0];
  std::fill(mi.begin(), mi.end(), 0);
  range_mi::barely_distorted::line<5>(
      line.data(), vacancy.data(), p_not_measured.data(),
      width.data(), num_cells, dtheta, mi.data());
  exact_mi[4] = mi[0];

  std::cout << std::endl;
  for (unsigned int i = 0; i < num_dimensions; i++) {
    std::cout << "d" << i + 1 << ": " << numerical_mi_[i] << ", " << exact_mi[i] <<
      ", difference: " << std::abs(numerical_mi_[i] - exact_mi[i]) <<
      std::endl;
  }
}
