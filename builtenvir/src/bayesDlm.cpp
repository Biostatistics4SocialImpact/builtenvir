
#define _USE_MATH_DEFINES
#include <cmath>
#include <cstdlib>
#include <algorithm>
#include <Rcpp.h>
#include <RcppEigen.h>

using namespace Rcpp;
using namespace Eigen;

#include "rnorm.h"
#include "mvrnorm.h"
#include "rgamma.h"



// =============================================================================
// Defines a routine to do Bayesian MCMC simulation from the full conditionals
// of the parameters of Distributed Lag Models of the type discussed by Baek et
// al. (2016) Epidemiology. Conjugate family priors are assumed. Computations
// use the Eigen linear algebra library. R shared object interface details
// abstracted with the Rcpp library.
//
//                                                   Andrew Whiteman - May, 2017
// =============================================================================






// perform argument checking and fit the bayesian DLM
List fitBayesDlm(
  const Map< VectorXd > &y,     // outcome vector
  const Map< MatrixXd > &D,     // entire design matrix
  const Map< VectorXd > &g,     // random coefficient diagonal
  const double alphaSigmaPrior, // prior on alpha_sigma^2 (residual variance)
  const double betaSigmaPrior,  // prior on beta_sigma^2
  const double alphaTauPrior,   // prior on alpha_tau^2 (random effect variance)
  const double betaTauPrior,    // prior on beta_tau^2
  const int nSims,              // total number of simulations (approximate)
  const int nSave,              // number of simulations to save
  const int burnin,             // number of simulations to discard
  const int seed                // RNG seed
);


// compute -2 * log-likelihood of the data
double computeNegLogLik(const int n, const double sigmaSq, const double rss);



// R interface
extern "C" SEXP bayesDlm(
  const SEXP ys,
  const SEXP Ds,
  const SEXP gs,
  const SEXP alphaSigmaPriorS,
  const SEXP betaSigmaPriorS,
  const SEXP alphaTauPriorS,
  const SEXP betaTauPriorS,
  const SEXP nSimS,
  const SEXP nSaveS,
  const SEXP burninS,
  const SEXP seedS
) {
  try {
    const Map< VectorXd > y(as< Map< VectorXd > >(ys));
    const Map< MatrixXd > D(as< Map< MatrixXd > >(Ds));
    const Map< VectorXd > g(as< Map< VectorXd > >(gs));
    const double alphaSigmaPrior = as< double >(alphaSigmaPriorS);
    const double betaSigmaPrior = as< double >(betaSigmaPriorS);
    const double alphaTauPrior = as< double >(alphaTauPriorS);
    const double betaTauPrior = as< double >(betaTauPriorS);
    const int nSims = as< int >(nSimS);
    const int nSave = as< int >(nSaveS);
    const int burnin = as< int >(burninS);
    const int seed = as< int >(seedS);

    return (wrap(fitBayesDlm(y, D, g,
			     alphaTauPrior, betaTauPrior,
			     alphaSigmaPrior, betaSigmaPrior,
			     nSims, nSave, burnin, seed
			     )));
  }
  catch (std::exception& __ex__) {
    forward_exception_to_r(__ex__);
  }
  catch (...) {
    ::Rf_error("C++ exception (unkown cause)");
  }

  return (R_NilValue);  // not reached
}






// Function Definitions
// -------------------------------------------------------------------

List fitBayesDlm(
  const Map< VectorXd > &y,
  const Map< MatrixXd > &D,
  const Map< VectorXd > &g,
  const double alphaSigmaPrior,
  const double betaSigmaPrior,
  const double alphaTauPrior,
  const double betaTauPrior,
  const int nSims,
  const int nSave,
  const int burnin,
  const int seed
) {
  const int k = D.cols(), n = y.size(), L = g.sum();
  const int thin = std::max((nSims - burnin) / nSave, 1);

  // argument checking
  if (D.rows() != n)
    throw (std::domain_error("Outcome/design dimension mismatch"));

  if (g.size() != k)
    throw (std::domain_error("Design/G dimension mismatch"));

  if (alphaSigmaPrior <= 0 || betaSigmaPrior <= 0 ||
      alphaTauPrior <= 0   || betaTauPrior <= 0)
    throw (std::domain_error("All Inv-Gamma priors must be > 0"));

  if (nSims < nSave + burnin || nSims <= 0 || nSave <= 0 || burnin < 0)
    throw (std::domain_error("Simulation parameters not viable"));

  // sufficient statistics
  const MatrixXd DtD = D.transpose() * D, G = g.asDiagonal();
  const VectorXd Dty = D.transpose() * y;

  // parameter storage
  VectorXd tauSqSims(nSave), sigmaSqSims(nSave), logLik(nSave);
  MatrixXd thetaSims(nSave, k), SigmaTheta(k, k);
  VectorXd fitted(n), residuals(n);
  int saveCount = 0;
  double rss;
  List result;

  // randomly initialize theta, tau^2
  std::srand(seed);
  VectorXd theta(k);
  for (int i = 0; i < k; i++)
    theta(i) = rnorm();
  double tauSq = 1 / rgamma(alphaTauPrior, 1 / betaTauPrior);
  double sigmaSq;

  // gibbs sampler - full conditionals
  for (int i = 0; saveCount < nSave; i++) {
    fitted = D * theta;
    residuals = y - fitted;
    rss = residuals.transpose() * residuals;
    sigmaSq = 1 / rgamma(alphaSigmaPrior + n / 2,
		       2 / (rss + 2 * betaSigmaPrior)
		       );
    SigmaTheta = (DtD / sigmaSq + G / tauSq).inverse();
    theta = mvrnorm(SigmaTheta * Dty / sigmaSq, SigmaTheta);
    tauSq = 1 / rgamma(alphaTauPrior + double(L - 2) / 2,
			 2 / (theta.transpose() * G * theta +
			      2 * betaTauPrior)
			 );

    // save simulations
    if (i > burnin && (i % thin) == 0) {
      thetaSims.row(saveCount) = theta;
      sigmaSqSims(saveCount) = sigmaSq;
      tauSqSims(saveCount) = tauSq;
      logLik(saveCount) = computeNegLogLik(n, sigmaSq, rss);
      saveCount++;
    }
  }

  // compute posterior means
  theta = thetaSims.colwise().mean();
  sigmaSq = sigmaSqSims.mean();
  tauSq = tauSqSims.mean();
  fitted = D * theta;
  residuals = y - fitted;
  rss = residuals.transpose() * residuals;

  result = List::create(
    Named("theta") = theta,
    Named("sigmaSq") = sigmaSq,
    Named("tauSq") = tauSq,
    Named("fitted") = fitted,
    Named("residuals") = residuals,
    Named("thetaSims") = thetaSims,
    Named("tauSqSims") = tauSqSims,
    Named("sigmaSqSims") = sigmaSqSims,
    Named("logLik") = -computeNegLogLik(n, sigmaSq, rss) / 2,
    Named("DIC") = 2 * logLik.mean() - computeNegLogLik(n, sigmaSq, rss)
    );

  return (result);
}


double computeNegLogLik(const int n, const double sigmaSq, const double rss) {
  return (n * std::log(sigmaSq) + rss / sigmaSq + n * std::log(2 * M_PI));
}
