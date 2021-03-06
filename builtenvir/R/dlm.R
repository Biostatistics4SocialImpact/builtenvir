
#' @title Distributed Lag Models
#'
#' @description
#' Fit distributed lag models to distance-profiled data.
#' \code{frequentist} method relies on \code{\link[nlme]{lme}}.
#' \code{bayesian} method uses Gibbs to sample from full conditionals,
#' with sampling parameters named in the \code{control} list.
#'
#' @param formula
#'   an object of class \code{"\link[stats]{formula}"}:
#'   a symbolic description of the model to be fitted. See Details.
#' @param data
#'   an optional data frame, list, or environment containing the
#'   variables of the fitted model.
#' @param subset
#'   optional vector specifying a subset of observations to be used in the
#'   fitting process.
#' @param na.action
#'   optional function that indicates what should happen when the data contains
#'   \code{NA}'s.
#' @param method
#'   method used to fit the DLM. Partial matching and capitalization allowed.
#' @param control
#'   a list of simulation control parameters. See Details.
#' @param ...
#'   Ignored for now.
#'
#' @details
#' Models are specified using familiar R \code{formula} syntax with one set of
#' lag terms returned by a given smoothing function (e.g. see \code{\link{cr}}).
#' The smoothing function can be any that returns a \code{\link{SmoothLag}} basis
#' object. Multiple lag terms or interactions with lag terms are not allowed, nor
#' are \code{NA} and missing-value lag terms supported. See Examples
#' for a basic call to \code{dlm} using the formula interface, and a cubic
#' radial lag basis specified via \code{cr}.
#'
#' The \code{control} list specifies additional optional \code{"bayes"} method
#' arguments, and may include: \code{n.sims}, the total number of simulations
#' to run (default = 5000); \code{n.save}, the total number of simulations to
#' save (default = 1000); \code{burnin}, number of simulations to discard from
#' the start of the chain (default = 2000);
#' \code{alpha.tau.prior}, prior (shape) parameter
#'   \eqn{\alpha_{\tau^2}}{\alpha[\tau^2]} (default = 0.1);
#' \code{beta.tau.prior}, prior (rate) parameter
#'   \eqn{\beta_{\tau^2}}{\alpha[\tau^2]} (default = 1e-6);
#' \code{alpha.sigma.prior}, prior (shape) parameter
#'   \eqn{\alpha_{\sigma^2}}{\alpha[\sigma^2]} (default = 0.1);
#' \code{beta.sigma.prior}, prior (rate) parameter
#'   \eqn{\beta_{\sigma^2}}{\beta[\sigma^2]} (default = 1e-6).
#'
#' The prior distribution heirarchy we assume in the Bayesian
#' formulation of the DLM is as follows:
#'
#' \deqn{y \sim N(D \theta, \sigma^2 I_n)}{y ~ N(D * \theta, \sigma^2 * In)}
#' \deqn{\theta \sim N(\mu_{\theta}, \Sigma_{\theta})}{\theta ~ N(\mu, \Sigma)}
#' \deqn{\theta_l \sim N(\mu_l, \tau^2)}{\theta(l) ~ N(\mu(l), \tau^2)}
#' \deqn{\sigma^2 \sim Inv-Gamma(\alpha_{\sigma^2}, \frac{1}{\beta_{\sigma^2}})}{\sigma^2 ~ Inv-Gamma(\alpha[\sigma^2], 1 / \beta[\sigma^2])}
#' \deqn{\tau^2 \sim Inv-Gamma(\alpha_{\tau^2}, \frac{1}{\beta_{\tau^2}})}{\tau^2 ~ Inv-Gamma(\alpha[\tau^2], 1 / \beta[\tau^2])}
#'
#' Where \eqn{l \in L}{l :- L} indexes the set of lag coefficients.
#'
#' @return
#' An S4 object that inherits from \code{"\link[=Dlm-class]{Dlm}"} containing
#' the results of the fitted model. If construction of this object fails,
#' \code{dlm} will issue a warning, and as a last resort attempt to return a
#' list with components of the fitted model.
#'
#'
#' @references Baek J, et al. (2016) Epidemiology 27(1):116-24.
#' (\href{https://www.ncbi.nlm.nih.gov/pubmed/26414942}{PubMed})
#'
#'
#' @examples
#' data (simdata)
#'
#' ## Setup distance count matrix and corresponding lag distances
#' X <- as.matrix(simdata[, -(1:3)])
#' lag <- seq(0.1, 10, length.out = ncol(X))
#'
#' fit <- dlm(Y ~ Age + Gender + cr(lag, X), data = simdata)
#' summary (fit)
#'
#' @seealso \code{\link[nlme]{lme}}, \code{\link{cr}},
#'   \code{\link{Dlm-class}}, \code{\link{FreqDlm}}, \code{\link{BayesDlm}}
#'

dlm <- function(formula, data, subset, na.action,
                method = c("frequentist", "bayesian", "Frequentist", "Bayesian"),
                control = list(),
                ...
                ) {
  mf <- match.call(expand.dots = FALSE)
  m <- match(c("formula", "data", "subset", "na.action"), names(mf), 0)
  mf <- mf[c(1, m)]
  mf$drop.unused.levels <- TRUE
  mf[[1]] <- quote(stats::model.frame)
  method <- tolower(match.arg(method))
  con <- list(n.sims = 5000, n.save = 1000, burnin = 2000,
              alpha.tau.prior = 0.1, beta.tau.prior = 1e-6,
              alpha.sigma.prior = 0.1, beta.sigma.prior = 1e-6,
              rng.seed = as.integer(Sys.time())
              )
  nmsC <- names(con)
  con[(namc <- names(control))] <- control
  if (length(noNms <- namc[!(namc %in% nmsC)]))
    warning ("unknown names in control: ", paste(noNms, collapse = ", "))
  if (con$n.sims < con$n.save + con$burnin) {
    cat ("N total simulation draws:", con$n.sims, "\n")
    cat ("N draws to save:", con$n.save, "\n")
    cat ("N burn-in iterations:", con$burnin, "\n")
    stop ("simulation parameters not viable")
  }
  mf <- eval(mf, parent.frame())
  mt <- attr(mf, "terms")
  if (is.empty.model(mt))
    stop ("empty model")
  which.lag <- which(vapply(mf, inherits, logical(1), "SmoothLag"))
  if (length(which.lag) != 1)
    stop ("there should be exactly 1 lag set")
  ## extract fixed effects and reorder lag terms to last
  tl <- attr(mt, "term.labels")
  sm.ndx <- grep(mf[[which.lag]]@signature, tl, fixed = TRUE)
  if (length(sm.ndx) > 1)
    stop ("interactions with smooth lag terms are not allowed")
  terms.ord <- c((1:length(tl))[-sm.ndx], sm.ndx)
  attr(mt, "factors") <- attr(mt, "factors")[, terms.ord, drop = FALSE]
  y <- model.response(mf, "numeric")
  X <- model.matrix(mt, mf)     # fixed effects
  Z_ <- mf[[which.lag]]@random  # random effects
  ## stats::model.frame applies na.action only to the resultant data.frame,
  ## whereas subsets are handled more intuitively by indexing on the variables
  if (!is.null(attr(mf, "na.action")))
    Z_ <- Z_[-attr(mf, "na.action"), ]
  K <- list(covariates = ncol(X) - ncol(mf[[which.lag]]),
            fixed = ncol(X),
            lag = ncol(mf[[which.lag]]) + ncol(Z_),
            random = ncol(Z_),
            total = ncol(X) + ncol(Z_)
            )
  D <- cbind(X, Z_)
  colnames (D)[(K$covariates + 1):K$total] <- paste(
    mf[[which.lag]]@signature, mf[[which.lag]]@basis@x, sep = ""
    )
  if (method == "frequentist") {
    group <- rep(1, length(y))
    fit <- lme(y ~ X - 1, random = list(group = pdIdent(~ Z_ - 1)))
  }
  else if (method == "bayesian") {
    g <- as.double(rep(0:1, c(K$fixed, K$random)))
    fit <- .Call("bayesDlm",
                 y, D, g, con$alpha.sigma.prior, con$beta.sigma.prior,
                 con$alpha.tau.prior, con$beta.tau.prior, con$n.sims,
                 con$n.save, con$burnin, con$rng.seed,
                 PACKAGE = "builtenvir"
                 )
  }
  DlmConstruct(method, fit, mf[[which.lag]]@basis, match.call(), D, K)
}





## Construct Dlm object from lme or "bayesDlm" output
## May make this exportable at some point

DlmConstruct <- function(method, fit, basis, call, D, K) {
  tryCatch ({
    g <- rep(0:1, c(K$fixed, K$random))
    if (method == "frequentist") {
      sig <- sigma(fit)
      tau <- sig * exp(unname(unlist(fit$modelStruct)))
      coefs <- unlist(coef(fit))
      names (coefs) <- colnames (D)
      .FreqDlm(.Dlm(method = method, call = call,
                    coefficients = coefs,
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




