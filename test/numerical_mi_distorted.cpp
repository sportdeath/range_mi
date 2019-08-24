#include <vector>
#include <iostream>
#include <cmath>
#include <random>
#include <fstream>

#include <range_mi/distorted.hpp>

#include "helpers.hpp"

// Define constants
double integration_step = 0.01;
double dtheta = 0.1;
unsigned int num_cells = 100;
double vacancy_scaling = 10;
double noise_dev = 4;
double noise_half_width = noise_dev * 4;
double noise_integration_step = 0.01;
unsigned int num_dimensions = 3;

double numerical_mi_distorted(
    const double * const pdf,
    const double * const pdf_distorted,
    unsigned pdf_size,
    unsigned pdf_distorted_size,
    unsigned int dimension,
    double dtheta,
    double noise_dev,
    double integration_step) {

  double h_Z = 0;
  for (unsigned int i = 0; i < pdf_distorted_size; i++) {
    double r = i * integration_step;

    h_Z +=
      pdf_distorted[i] *
      std::pow(r, dimension - 1) *
      -std::log(pdf_distorted[i]) *
      integration_step *
      dtheta;
  }

  double h_Z_given_M = 0;
  for (unsigned int i = 0; i < pdf_size; i++) {
    double r = i * integration_step;

    h_Z_given_M +=
      std::log(noise_dev * std::sqrt(2 * M_PI * M_E)) *
      pdf[i] *
      std::pow(r, dimension - 1) *
      integration_step *
      dtheta;
  }

  if (dimension == 3) {
    h_Z_given_M += 
      noise_dev * noise_dev *
      std::log(noise_dev * std::sqrt(2 * M_PI * M_E * M_E * M_E))
      * dtheta;
  }

  return h_Z - h_Z_given_M;
}

int main() {
  // Initialize the inputs and outputs
  std::vector<double> vacancy(num_cells);
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

  std::cout << "Computing the pdf" << std::endl;

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

  std::cout << "Convolving pdf with noise" << std::endl;

  unsigned int noise_size = (2 * noise_half_width)/integration_step;
  unsigned int pdf_distorted_size = pdf_size + noise_size - 1;
  std::vector<double> pdf_distorted(pdf_distorted_size, 0);
  for (unsigned int i = 0; i < pdf_size; i++) {
    for (unsigned int j = 0; j < noise_size; j++) {
      double r = j * integration_step - noise_half_width;
      pdf_distorted[i + j] +=
        pdf[i] *
        range_mi::distorted::normal_pdf(r, 0, noise_dev) *
        integration_step;
    }
  }

  std::cout << "Computing the barely distorted mutual information numerically..." << std::endl;

  std::vector<double> numerical_mi_(num_dimensions);
  for (unsigned int i = 0; i < num_dimensions; i++) {
    numerical_mi_[i] = numerical_mi_distorted(
        pdf.data(),
        pdf_distorted.data(),
        pdf_size,
        pdf_distorted_size,
        i + 1,
        dtheta,
        noise_dev,
        integration_step);
  }

  std::cout << "Computing the barely distorted mutual information exactly..." << std::endl;

  std::vector<double> exact_mi(num_dimensions, 0);
  std::vector<double> mi(num_cells);
  std::vector<double> pdf_((num_cells + 2 * noise_half_width)/noise_integration_step);

  // Clear the output
  std::fill(mi.begin(), mi.end(), 0);
  range_mi::distorted::line<1>(
      line.data(), vacancy.data(), width.data(),
      num_cells, noise_dev, noise_half_width, noise_integration_step,
      dtheta, pdf_.data(), mi.data());
  exact_mi[0] = mi[0];
  std::fill(mi.begin(), mi.end(), 0);
  range_mi::distorted::line<2>(
      line.data(), vacancy.data(), width.data(),
      num_cells, noise_dev, noise_half_width, noise_integration_step,
      dtheta, pdf_.data(), mi.data());
  exact_mi[1] = mi[0];
  std::fill(mi.begin(), mi.end(), 0);
  range_mi::distorted::line<3>(
      line.data(), vacancy.data(), width.data(),
      num_cells, noise_dev, noise_half_width, noise_integration_step,
      dtheta, pdf_.data(), mi.data());
  exact_mi[2] = mi[0];

  std::cout << std::endl;
  for (unsigned int i = 0; i < num_dimensions; i++) {
    std::cout << "d" << i + 1 << ": " << numerical_mi_[i] << ", " << exact_mi[i] << std::endl;
  }
}
