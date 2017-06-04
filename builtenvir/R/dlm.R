
#' @title Distributed Lag Models
#'
#' @description
#' Fit distributed lag models to distance-profiled data.
#' \code{frequentist} method relies on \code{nlme::lme}.
#' \code{bayesian} method uses Gibbs to sample from full conditionals,
#' with sampling parameters named in the \code{control} list.
#'
#' @param y model response variable
#' @param X matrix of lag count data
#' @param Z covariate matrix  --> TODO: make optional?
#' @param lags lag distances for columns of \code{X}
#' @param intercept \code{c(TRUE, FALSE)}; \code{TRUE} if an intercept
#'   should be included in the model
#' @param method \code{c("frequentist", "bayesian")}; method used to fit the
#'   DLM. Partial matching and capitalization allowed
#' @param constrol \code{list(...)}; a list of control parameters.
#'   See 'Details'
#'
#' @details
#' \code{control} list specifies additional optional \code{"bayes"} method
#' arguments, and may include: \code{n.sims}, the total number of simulations
#' to run; \code{n.save}, the total number of simulations to save;
#' \code{burnin}, number of simulations to discard from the start of the chain;
#' \code{alpha.tau.prior}, prior (shape) parameter \eqn{\alpha_{\tau^2}};
#' \code{beta.tau.prior}, prior (rate) parameter \eqn{\beta_{\tau^2}};
#' \code{alpha.sigma.prior}, prior (shape) parameter \eqn{\alpha_{\sigma^2}};
#' \code{beta.sigma.prior}, prior (rate) parameter \eqn{\beta_{\sigma^2}}.
#'
#' The prior distribution heirarchy we assume in the Bayesian
#' formulation of the DLM is as follows:
#'
#' \deqn{
#'   y \sim N(D \theta, \tau^2 I_n) \\
#'   \theta \sim N(\mu_{\theta}, \Sigma_{\theta}) \\
#'   \theta_{(lag)} \sim N(\mu_{\theta_{(lag)}}, \sigma^2 I_L) \\
#'   \tau^2 \sim Inv-Gamma(\alpha_{\tau^2}, \frac{1}{\beta_{\tau^2}}) \\
#'   \sigma^2 \sim Inv-Gamma(\alpha_{\sigma^2}, \frac{1}{\beta_{\sigma^2}})
#' }
#'
#'
#' @return
#' An S4 \code{Dlm} object containing the results of the fitted model
#'
#' @references Baek J, et al. (2016) Epidemiology 27(1):116-24.
#' (\href{https://www.ncbi.nlm.nih.gov/pubmed/26414942}{PubMed})
#'
#'
#' @examples
#' data (simdata)
#' y <- simdata[, 1]
#' Z <- as.matrix(simdata[, 2:3])
#' X <- as.matrix(simdata[, 4:ncol(simdata)])
#' lg <- seq(0.1, 10, length.out = ncol(X))
#' fit <- dlm(y, X, Z, lg)
#' fit
#'


## library (nlme)

#' @export
dlm <- function(y, X, Z, lags,
                intercept = TRUE,
                method = c("frequentist", "bayesian", "Frequentist", "Bayesian"),
                control = list(),
                ...
                ) {
  method <- tolower(match.arg(method))
  con <- list(n.sims = 5000, n.save = 1000, burnin = 2000,
              alpha.tau.prior = 0.1, beta.tau.prior = 1e-6,
              alpha.sigma.prior = 0.1, beta.sigma.prior = 1e-6,
              rng.seed = as.integer(Sys.time())
              )
  con[(namc <- names(control))] <- control
  if (length(noNms <- namc[!(namc %in% names(con))]))
    warning ("unknown names in control: ", paste(noNms, collapse = ", "))
  if (con$n.sims < con$n.save + con$burnin) {
    cat ("N total simulation draws:", con$n.sims, "\n")
    cat ("N draws to save:", con$n.save, "\n")
    cat ("N burn-in iterations:", con$burnin, "\n")
    stop ("simulation parameters not viable")
  }
  if (intercept)
    Z <- cbind("(Intercept)" = 1, Z)
  basis <- lag.basis(lags)
  X_ <- X %*% basis@C0
  Z_ <- X %*% basis@K1
  K <- list(covariates = ncol(Z),
            fixed = ncol(Z) + ncol(X_),
            lag = ncol(X_) + ncol(Z_),
            random = ncol(Z_),
            total = ncol(Z) + ncol(X_) + ncol(Z_)
            )
  D <- cbind(Z, X_, Z_)
  colnames (D) <- c(.remake.names(Z), .remake.names(X_), .remake.names(Z_))
  ## probably a more robust way to handle missing data. Check this and re-code
  ## if necessary
  ok <- !(is.na(rowSums(D)) | is.na(y))
  D <- D[ok, ]
  if (!any(ok))
    stop ("no non-missing cases")
  if (method == "frequentist") {
    group <- rep(1, length(y))
    fit <- lme(y ~ Z + X_ - 1, random = list(group = pdIdent(~ Z_ - 1)))
  }
  else if (method == "bayesian") {
    y <- y[ok]
    g <- as.double(rep(0:1, c(K$fixed, K$random)))
    fit <- .Call("bayesDlm",
                 y, D, g, con$alpha.sigma.prior, con$beta.sigma.prior,
                 con$alpha.tau.prior, con$beta.tau.prior, con$n.sims,
                 con$n.save, con$burnin, con$rng.seed,
                 PACKAGE = "builtenvir"
                 )
  }
  DlmConstruct(method, fit, basis, match.call(), D, K)
}





## Construct Dlm object from lme or "bayesDlm" output
## May make this exportable at some point

DlmConstruct <- function(method, fit, basis, call, D, K) {
  tryCatch ({
    g <- rep(0:1, c(K$fixed, K$random))
    if (method == "frequentist") {
      sig <- sigma(fit)
      tau <- sig * exp(unname(unlist(fit$modelStruct)))
      .FreqDlm(.Dlm(method = method, call = call,
                    coefficients = unlist(coef(fit)),
                    sigma = sig,
                    tau = tau,
                    logLik = as.numeric(logLik(fit)),
                    vcov = qr.solve(crossprod(D) / sig^2 + diag(g) / tau^2),
                    basis = basis,
                    N = length(residuals(fit)),
                    K = K
                    ),
               lme = fit
               )
    }
    else if (method == "bayesian") {
      names (fit$theta) <- colnames (fit$thetaSims) <- colnames (D)
      .BayesDlm(.Dlm(method = method, call = call,
                     coefficients = fit$theta,
                     sigma = sqrt(fit$sigmaSq),
                     tau = sqrt(fit$tauSq),
                     logLik = fit$logLik,
                     vcov = cov(fit$thetaSims),
                     basis = basis,
                     N = nrow(D),
                     K = K
                     ),
                fitted = fit$fitted,
                DIC = fit$DIC,
                residuals = fit$residuals,
                sims = list(coefficients = fit$thetaSims,
                            sigma = sqrt(fit$sigmaSqSims),
                            tau = sqrt(fit$tauSqSims)
                            )
                )

    }
    else stop ("unrecognized method '", method, "'")
  },
  error = function(e) {
    warning (e)
    warning ("returning component list")
    list (method = method, fit = fit, basis = basis, call = call)
  })
}
