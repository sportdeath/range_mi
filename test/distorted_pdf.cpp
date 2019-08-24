#include <vector>
#include <fstream>

#include <range_mi/distorted.hpp>

std::vector<double> vacancy = {0.2, 0.7, 0.4, 0.2};
std::vector<double> width = {0.75, 2, 1, 1.5};
std::vector<unsigned int> line = {0, 1, 2, 3};
//double noise_dev = 0.06;
double noise_dev = 0.2;
double noise_half_width = noise_dev * 4;
double integration_step = 0.001;

double width_scaling = 0.5;

int main() {
  
  double width_sum = 0;
  for (unsigned int i = 0; i < width.size(); i++) {
    width[i] *= width_scaling;
    width_sum += width[i];
  }

  std::vector<double> pdf(2 * (width_sum + 2 * noise_half_width)/integration_step);
  unsigned int pdf_size;
  range_mi::distorted::range_pdf(
      line.data(),
      vacancy.data(),
      width.data(),
      vacancy.size(),
      noise_dev,
      noise_half_width,
      integration_step,
      width_sum,
      pdf.data(),
      pdf_size);

  // Write it out
  std::ofstream f;
  f.open("noise.txt");
  for (unsigned int i = 0; i < pdf_size; i++) {
    f << pdf[i] << std::endl;
  }
  f.close();
}
