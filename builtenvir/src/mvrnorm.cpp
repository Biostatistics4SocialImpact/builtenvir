
#include <RcppEigen.h>

using namespace Eigen;

#include "rnorm.h"


VectorXd mvrnorm(const VectorXd &mu, const MatrixXd &Sigma) {
  const MatrixXd L = Sigma.llt().matrixL();
  const int k = mu.size();
  VectorXd z(k);
  for (int i = 0; i < k; i++)
    z(i) = rnorm();
  return (mu + L * z);
}
