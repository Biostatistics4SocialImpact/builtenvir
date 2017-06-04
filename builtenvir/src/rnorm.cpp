
#include <random>

using namespace std;

#include "rnorm.h"


double rnorm(const double mu, const double sigma) {
  static default_random_engine rgenNormal;
  static normal_distribution< double > z(0, 1);
  return (z(rgenNormal) * sigma + mu);
}
