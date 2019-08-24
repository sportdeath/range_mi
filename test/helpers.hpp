#include <random>

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

void numerical_pdf(
    const unsigned int * const line,
    const double * const vacancy,
    const double * const width,
    unsigned int num_cells,
    double integration_step,
    double * const pdf,
    unsigned int & pdf_size) {

  unsigned int i = 0;
  double r = 0;
  double width_sum = 0;
  double pdf_decay = 1;
  pdf_size = 0;
  while (i < num_cells) {
    unsigned int j = line[i];

    // Compute the pdf
    pdf[pdf_size++] = 
      pdf_decay * (
          -std::pow(vacancy[j], r - width_sum) *
          std::log(vacancy[j]));

    // Make a step
    r += integration_step;
    if (r - width_sum > width[i]) {
      // Cell completed,
      // move to the next cell!
      pdf_decay *= std::pow(vacancy[j], width[i]);
      width_sum += width[i];
      i++;
    }
  }
}
