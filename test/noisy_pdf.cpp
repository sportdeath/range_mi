#include <vector>
#include <numeric>
#include <fstream>

#include <range_entropy/expected_noisy.hpp>

// The map
std::vector<double> vacancy = {0.2, 0.7, 0.4, 0.2};
std::vector<double> width = {0.1875, 0.5, 0.25, 0.375};
std::vector<double> p_not_measured = {1, 0.5, 1, 0.25};

double noise_dev = 0.035;
double noise_width = noise_dev * 5;
double step_size = 0.001;

int main() {
  // Make a line
  std::vector<unsigned int> line(vacancy.size());
  std::iota(line.begin(), line.end(), 0);

  double pdf_width = 0;
  for (unsigned int i = 0; i < width.size(); i++) {
    pdf_width += width[i];
  }

  // Initialize a place to put the pdf
  std::vector<double> pdf(2 * (pdf_width + 2 * noise_width)/step_size);

  unsigned int pdf_size;
  range_entropy::expected_noisy::pdf(
      line.data(),
      vacancy.data(),
      width.data(),
      line.size(),
      noise_dev,
      noise_width,
      step_size,
      pdf_width,
      pdf.data(),
      pdf_size);

  std::ofstream f;
  f.open ("noisy_pdf.txt");
  // Plot the density function
  for (unsigned int i = 0; i < pdf_size; i++) {
    f << i * step_size - noise_width << " " << pdf[i] << std::endl;
  }
  f.close();
}
