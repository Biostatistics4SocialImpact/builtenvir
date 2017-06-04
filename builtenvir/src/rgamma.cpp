
#include <random>

using namespace std;

#include "rgamma.h"


double rgamma(const double shape, const double scale) {
  static default_random_engine rgenGamma;
  gamma_distribution< double > G(shape, scale);
  return (G(rgenGamma));
}
