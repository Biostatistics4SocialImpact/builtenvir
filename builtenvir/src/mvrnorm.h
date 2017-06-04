
#ifndef _builtenvr_MVRNORM_
#define _builtenvr_MVRNORM_

#include <RcppEigen.h>


VectorXd mvrnorm(const VectorXd &mu, const MatrixXd &Sigma);

#endif  // _builtenvr_MVRNORM_
