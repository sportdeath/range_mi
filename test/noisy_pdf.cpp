#include <vector>
#include <numeric>
#include <fstream>
#include <iostream>

#include <range_entropy/expected_noisy.hpp>

// The map
std::vector<double> vacancy = {0.2, 0.7, 0.4, 0.2};
std::vector<double> width = {0.1875, 0.5, 0.25, 0.375};
std::vector<double> p_not_measured = {1, 0.5, 1, 0.25};

double noise_dev = 0.05;
double noise_half_width = noise_dev * 5;
double integration_step = 0.001;

int main() {

  // Make a line
  std::vector<unsigned int> line(vacancy.size());
  std::iota(line.begin(), line.end(), 0);

  // Determine how long the density function should be
  double pdf_width = 0;
  for (unsigned int i = 0; i < width.size(); i++) {
    pdf_width += width[i];
  }
  std::cout << "noise_half_width: " << noise_half_width << std::endl;
  std::cout << "pdf_width: " << pdf_width << std::endl;

  // Initialize a place to put the convolution
  std::vector<double> pdf((pdf_width + 2 * noise_half_width)/integration_step + 1);
  std::cout << "pdf size: " << pdf.size() << std::endl;

  // Construct the noise
  unsigned int noise_size = (noise_half_width * 2)/integration_step;
  std::vector<double> noise(noise_size);
  double r = -noise_half_width;
  for (unsigned int i = 0; i < noise_size; i++) {
    noise[i] = range_entropy::expected_noisy::normal_pdf(r, noise_dev);
    r += integration_step;
  }

  unsigned int pdf_size;
  range_entropy::expected_noisy::pdf(
      line.data(),
      vacancy.data(),
      width.data(),
      line.size(),
      noise.data(),
      noise_size,
      integration_step,
      pdf_width,
      pdf.data(),
      pdf_size);

  std::cout << "noise size: " << noise_size << std::endl;
  std::cout << "convolve size: " << pdf_size << std::endl;

  std::ofstream f;
  f.open ("noisy_pdf.txt");
  // Plot the density function
  r = -noise_half_width;
  for (unsigned int i = 0; i < pdf_size; i++) {
    f << r << " " << pdf[i] << std::endl;
    r += integration_step;
  }
  f.close();
}
