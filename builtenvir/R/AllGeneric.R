

## New Generics
## -------------------------------------------------------------------



#' @title Extract Lag Basis Matrix
#'
#' @description
#' Extract lag basis matrix, \eqn{\Omega}
#'
#' @return
#' A square numeric matrix
#'
#' @details
#' \code{omega} is S4 generic.
setGeneric("omega", function(object, ...) standardGeneric("omega"))



#' @title Extract Fixed Effects
#'
#' @description
#' Returns a vector of fixed effects coefficients from a fitted
#' mixed effects model.
#'
#' @return
#' A numeric vector
#'
#' @details
#' \code{fixef} is S4 generic.
#'
setGeneric("fixef",
           function(object) standardGeneric("fixef")
           )



#' @title Extract Random Effects
#'
#' @description
#' Returns a vector of random effects coefficients from a fitted
#' mixed effects model.
#'
#' @return
#' A numeric vector
#'
#' @details
#' \code{ranef} is S4 generic.
#'
setGeneric("ranef",
           function(object) standardGeneric("ranef")
           )



#' @title Extract Standard Errors of Fixed Effects
#'
#' @description
#' Returns a vector of fixed effects standard errors from a fitted
#' mixed effects model.
#'
#' @return
#' A numeric vector
#'
#' @details
#' \code{fixef} is S4 generic.
#'
setGeneric("se.fixef",
           function(object) standardGeneric("se.fixef")
           )


#' @title Extract Standard Errors of Random Effects
#'
#' @description
#' Returns a vector of random effects standard errors from a fitted
#' mixed effects model.
#'
#' @return
#' A numeric vector
#'
#' @details
#' \code{fixef} is S4 generic.
#'
setGeneric("se.ranef",
           function(object) standardGeneric("se.ranef")
           )



#' @title Extract Distributed Lag Scale Matrix
#'
#' @description
#' Return lag coefficient scale matrix, S, such that the distributed lag
#' coefficients fit by the model are obtained via the transformation
#' \eqn{\beta = S \theta}{\beta = S * \theta}. S should be invertable.
#'
#' @param object a fitted model object
#'
#' @return A square numeric matrix
#'
#' @details
#' \code{scaleMat} is S4 generic.
#'
setGeneric("scaleMat",
           function(object, ...) standardGeneric("scaleMat")
           )



#' @title Extract Natural Scale Model Parameters 'Theta'
#'
#' @description
#' Extract \eqn{\theta} coefficients from a fitted model
#'
#' @param object a fitted model object
#'
#' @return A numeric vector
#'
#' @details
#' Returns the natural scale coefficients
#' \eqn{\theta = S^{-1} \beta}{\theta = S^-1 * \beta} given some invertable
#' scale matrix, S.
#'
#' \code{theta} is S4 generic.
#'
setGeneric("theta",
           function(object, ...) standardGeneric("theta")
           )


#' @title Extract Natural Scale Coefficient Covariance Matrix
#'
#' @description
#' Extract \eqn{var \theta} matrix from a fitted model
#'
#' @param object a fitted model object
#'
#' @return A square numeric matrix
#'
#' @details
#' Returns the covariance matrix for the natural scale coefficients,
#' \eqn{\theta}.
#'
#' \code{vcovTheta} is S4 generic.
#'
setGeneric("vcovTheta",
           function(object, ...) standardGeneric("vcovTheta")
           )


#' @title Extract Coefficient Simulations
#'
#' @description
#' Extract simulated coefficients from a model fitted with MCMC
#'
#' @param object a fitted model object
#'
#' @return A numeric matrix
#'
#' @details
#' \code{coef.sims} is S4 generic.
setGeneric("coef.sims", function(object, ...) standardGeneric("coef.sims"))





## Methods for LagBasis class
## -------------------------------------------------------------------


#' @describeIn omega Method for "\code{LagBasis}" objects
setMethod("omega", signature = "LagBasis",
          function(object, ...)  cbind(object@C0, object@K1)
          )


#' @title Predict New Values for Fitted Lag Basis
#'
#' @param object A \code{\link{LagBasis}} object
#'
#' @description
#' Not yet implemented
#'
predict.LagBasis <- function(object, ...) {
  stop ("Not yet Implemented")
}








## Methods for Dlm class
## -------------------------------------------------------------------

#' @describeIn fixef S4 Method for \code{"\link[=Dlm-class]{Dlm}"} Objects
setMethod("fixef", signature = "Dlm",
          function(object) coef(object)[1:object@K$fixed]
          )




#' @describeIn ranef S4 Method for \code{"\link[=Dlm-class]{Dlm}"} Objects
setMethod("ranef", signature = "Dlm",
          function(object) {
            ran.index <- object@K$fixed + 1
            coef(object)[ran.index:object@K$total]
          })



#' @describeIn se.fixef S4 Method for \code{"\link[=Dlm-class]{Dlm}"} Objects
setMethod("se.fixef", signature = "Dlm",
          function(object)
            sqrt(diag(vcov(object))[1:object@K$fixed])
          )


#' @describeIn se.ranef S4 Method for \code{"\link[=Dlm-class]{Dlm}"} Objects
setMethod("se.ranef", signature = "Dlm", function(object) object@tau)



#' @describeIn scaleMat S4 Method for \code{"\link[=Dlm-class]{Dlm}"} Objects
setMethod("scaleMat", signature = "Dlm",
          function(object, ...) {
            k <- object@K$total
            lag.index <- object@K$covariates + 1
            S <- diag(k)
            S[lag.index:k, lag.index:k] <- omega(object@basis)
            return (S)
          })


#' @describeIn theta S4 Method for \code{"\link[=Dlm-class]{Dlm}"} Objects
setMethod("theta", signature = "Dlm", function(object, ...) object@coefficients)


#' @describeIn vcovTheta S4 Method for \code{"\link[=Dlm-class]{Dlm}"} Objects
setMethod("vcovTheta", signature = "Dlm", function(object, ...) object@vcov)


##
setMethod("show", signature = "Dlm",
          function(object) {
            tryCatch ({
              cat ("Call: "); print (object@call)
              method <- if (object@method == "frequentist") {
                  object@lme$method
                } else if (object@method == "bayesian") {
                  "MCMC"
                } else {
                  NULL
                }
              if (!is.null(method))
                cat ("(", object@method, ") distributed lag model fit via ",
                     method, "\n", sep = ""
                     )
            }, error = function(e) NULL
            )
            cat ("Number of Observations:", object@N, "\n")
            cat ("Log-Likelihood:", logLik(object), "\n")
            cat ("Fixed Effects:\n")
            tryCatch ({
              semat <- cbind("Random Effects" = se.ranef(object),
                             "Residuals" = sigma(object)
                             )
              rownames (semat) <- "Std Error"
              printCoefmat(t(fixef(object)), has.Pvalue = FALSE); cat ("\n")
              printCoefmat(semat, has.Pvalue = FALSE)
            }, error = function(e) NULL
            )
          })




#' @title Extract DLM Model Fitted Values
#'
#' @description
#' Extract the predicted \code{y} values from a fitted by \code{\link{dlm}}
#'
#' @param object a fitted \code{"\link[=Dlm-class]{Dlm}"} object
#'
#' @return A numeric vector
#'
fitted.Dlm <- function(object) {
  tryCatch ({
    object@fitted
  },
  error = function(e) {
    tryCatch ({
      fitted(object@lme)
    },
    error = function(e) {
      warning ("object not properly initialized")
      NA_real_
    }) })
}



#' @title Extract DLM Model Residuals
#'
#' @description
#' Extract residuals from a model fitted by \code{\link{dlm}}
#'
#' @param object a fitted \code{"\link[=Dlm-class]{Dlm}"} object
#'
#' @return A numeric vector
#'
residuals.Dlm <- function(object) {
  tryCatch ({
    object@residuals
  },
  error = function(e) {
    tryCatch ({
      residuals(object@lme)
    },
    error = function(e) {
      ## warning ("object not properly initialized")
      NA_real_
    }) })
}



#' @title Generate New Predictions from a Fitted DLM
#'
#' @description
#' Not yet implemented
#'
predict.Dlm <- function(object, ...) {
  stop ("not yet implemented")
}


#' @title Extract Distributed Lag Model Coefficients
#'
#' @description
#' Return DLM coefficients
#'
#' @return
#' A numeric vector of fitted coefficients.
#'
#' @details
#'   Note that these are the basis-scaled coefficients. To extract
#'   coefficients on the scale of the original design (\eqn{\theta}),
#'   see \code{\link{theta}}.
#'
coef.Dlm <- function(object, ...) {
  b <- c(scaleMat(object) %*% theta(object))
  names (b) <- names (theta(object))
  return (b)
}


#' @title Extract Residual Standard Deviation 'Sigma'
#'
#' @description
#' Return observed residual standard deviation
#'
#' @return
#' A positive numeric scalar: the residual standard deviation
#'
sigma.Dlm <- function(object, ...) object@sigma




#' @title Extract Log-Likelihood
#'
#' @description
#' Log likelihood of the observed data given the model
#'
#' @return
#' An object of class \code{logLik} containing the log-likelihood of the
#' data given the fitted DLM
#'
logLik.Dlm <- function(object, ...) {
  ll <- object@logLik
  attr(ll, "nall") <- attr(ll, "nobs") <- object@N
  attr(ll, "df") <- object@K$fixed + 2
  ## fixed effects + 2 variance parameters
  class(ll) <- "logLik"
  return (ll)
}



#' @title Akaike Information Criterion
#'
#' @description
#' Compute AIC statistic for fitted DLM
#'
#' @return
#' A numeric scalar
#'
AIC.Dlm <- function(object, ..., k = 2) {
  ll <- logLik(object)
  k <- max(k, attr(ll, "df"))
  c(2 * (k - ll))
}



#' @title "Bayesian" Information Criterion
#'
#' @description
#' Compute BIC statistic for fitted DLM
#'
#' @return
#' A numeric scalar: BIC
#'
BIC.Dlm <- function(object, ...) {
  ll <- logLik(object)
  k <- max(2, attr(ll, "df"))
  c(log(object@N) * k - 2 * ll)
}




#' @title Calculate the Distributed Lag Variance-Covariance Matrix for a
#'   Fitted Model Object
#'
#' @description
#' Returns the variance-covariance matrix of the fitted model coefficients
#' (including the random effects terms).
#'
#' @return
#' A numeric matrix: the coefficient variance-covariance matrix
#'
vcov.Dlm <- function(object, ...) {
  S <- scaleMat(object)
  S %*% vcovTheta(object) %*% t(S)
}






## Methods for Class FreqDlm
## -------------------------------------------------------------------

#' @title Summarize a Fitted Dlm Object
#'
#' @description
#' Computes likelihood based statistics, t-statistic summaries, and correlations
#' between fixed effects terms after the \code{nlme} package. Also computes
#' rough distributional summaries for distributed lag terms (estimates and
#' intervals, with confidence coefficient level determined by the \code{level}
#' parameter).
#'
#' @param object a fitted \code{"\link[=Dlm-class]{Dlm}"} object
#'
#' @param level desired confidence coefficient (or credible interval) level for
#'   summarizing distributed lag terms
#'
#' @param ... ignored
#'
#' @return
#' An S4 object of class \code{"\link{SummaryDlm}"}
#'
#' @name summary.Dlm
NULL

#' @rdname summary.Dlm
summary.FreqDlm <- function(object, level = 0.95, ...) {
  method <- c(object@method, object@lme$method)
  like.stats <- c(AIC = AIC(object), BIC = BIC(object),
                  logLik = c(logLik(object))
                  )
  variances <- c("Random Effects" = object@tau^2,
                 "Residuals" = sigma(object)^2
                 )
  est <- fixef(object)
  se <- se.fixef(object)
  t.val <- est / se
  df <- object@N - object@K$fixed
  tTable <- cbind("Estimate" = est,
                  "Std. Error" = se,
                  "t value" = t.val,
                  "Pr(>|t|)" = 2 * (1 - pt(abs(t.val), df))
                  )
  cor.fixef <- vcov(object)[1:object@K$fixed, 1:object@K$fixed]
  cor.fixef <- diag(sqrt(diag(cor.fixef))^-1) %*% cor.fixef %*%
    diag(sqrt(diag(cor.fixef))^-1)
  rownames (cor.fixef) <- colnames (cor.fixef) <- names(est)
  level <- max(min(level, 1), 0)
  q <- abs(qnorm((1 - level) / 2))
  lag.coefs <- coef(object) +
    (sqrt(diag(vcov(object))) %*% t(c(0, -1, 1) * q))
  lag.coefs <- lag.coefs[-(1:object@K$covariates), ]
  colnames (lag.coefs) <- c("DL Coefficients",
                            sprintf("%.1f%%", pnorm(c(-1, 1) * q) * 100))
  .SummaryDlm(method = method,
              call = object@call,
              coefficients = coef(object),
              likelihood.stats = like.stats,
              variances = variances,
              tTable = tTable,
              lag.coefs = lag.coefs,
              cor.fixef = cor.fixef,
              residuals = residuals(object),
              N = object@N,
              K = object@K
              )
}





## Methods for Class BayesDlm
## -------------------------------------------------------------------

#' @rdname summary.Dlm
summary.BayesDlm <- function(object, level = 0.95, ...) {
  method <- c(object@method, "MCMC")
  like.stats <- c(DIC = object@DIC, AIC = AIC(object),
                  BIC = BIC(object), logLik = c(logLik(object))
                  )
  variances <- c("Random Effects" = object@tau^2,
                 "Residuals" = sigma(object)^2
                 )
  est <- fixef(object)
  se <- se.fixef(object)
  cs <- coef.sims(object)
  p.vals <- apply(cs[, 1:object@K$fixed], 2, sims.pval)
  tTable <- cbind("Estimate" = est,
                  "Std. Error" = se,
                  ## "Pr(|Est| > 0)" = p.vals
                  "Pr(>|Est|)" = p.vals
                  )
  cor.fixef <- vcov(object)[1:object@K$fixed, 1:object@K$fixed]
  cor.fixef <- diag(sqrt(diag(cor.fixef))^-1) %*% cor.fixef %*%
    diag(sqrt(diag(cor.fixef))^-1)
  rownames (cor.fixef) <- colnames (cor.fixef) <- names(est)
  level <- max(min(level, 1), 0)
  probs <- (1 - level) / 2 + c(0, level)
  lag.coefs <- cbind("DL Coefficients" = coef(object),
                     t(apply(cs, 2, quantile, probs = probs))
                     )[-(1:object@K$covariates), ]
  .SummaryDlm(method = method,
              call = object@call,
              coefficients = coef(object),
              likelihood.stats = like.stats,
              variances = variances,
              tTable = tTable,
              lag.coefs = lag.coefs,
              cor.fixef = cor.fixef,
              residuals = residuals(object),
              N = object@N,
              K = object@K
              )
}







#' @title Extract Distributed Lag Coefficient Simulation Draws
#'
#' @description
#' \code{Dlm} objects fit by the Bayesian method store a matrix of
#' simulated coefficients drawn via MCMC from their full conditional
#' distribution. Extract the distributed lag-scaled simulated coefficient
#' matrix.
#'
#' @return
#' A numeric matrix containing the scaled simulated coefficient matrix
#'
setMethod("coef.sims", signature = "BayesDlm",
          function(object, ...) {
            object@sims$coefficients %*% t(scaleMat(object))
          })


## Methods for Class SummaryDlm
## -------------------------------------------------------------------

setMethod("show", signature = "SummaryDlm",
          function(object) {
            vars <- t(object@variances)
            rownames (vars) <- "Var"
            cors <- object@cor.fixef
            cors[upper.tri(cors, TRUE)] <- NA
            if (object@method[1] != "")
              cat ("(", object@method[1], ") ", sep = "")
            cat ("distributed lag model fit via", object@method[2], "\n")
            cat ("Call: "); print (object@call)
            cat ("\nNumber of observations:", object@N, "\n")
            cat ("Standardized Residuals:\n")
            print (summary(object@residuals / sqrt(object@variances[2])))
            cat ("\n")
            printCoefmat(vars, has.Pvalue = FALSE)
            cat ("\nFixed Effects:\n")
            printCoefmat(object@tTable, has.Pvalue = TRUE)
            cat ("\nCorrelation:\n")
            tryCatch ({
              printCoefmat(cors[-1, -ncol(cors)], has.Pvalue = FALSE,
                           na.print = ""
                           )
            }, error = function(e) NULL
            )
            cat ("\n")
          })







## Methods for SmoothLag Class
## -------------------------------------------------------------------

## getGeneric("[")
## ## standardGeneric for "[" defined from package "base"
## ## function (x, i, j, ..., drop = TRUE)
## ## standardGeneric("[", .Primitive("["))

## getGeneric("[<-")
## ## standardGeneric for "[<-" defined from package "base"
## ## function (x, i, j, ..., value)
## ## standardGeneric("[<-", .Primitive("[<-"))

## x[]
setMethod("[", signature = c(x = "SmoothLag", i = "missing", j = "missing",
                             drop = "ANY"),
          function(x, i, j, ..., drop) x
          )


setMethod("[", signature = c(x = "SmoothLag", i = "numeric", j = "missing",
                             drop = "ANY"),
          function(x, i, j, ..., drop) {
            if (!missing(drop))
              if (!is.logical(drop) || drop)
                warning ("argument 'drop' ignored. Only drop = FALSE accepted")
            .SmoothLag(x@.Data[i, , drop = FALSE],
                       basis = x@basis,
                       random = x@random[i, , drop = FALSE],
                       signature = x@signature
                       )
          })




setMethod("[", signature = c(x = "SmoothLag", i = "missing", j = "numeric",
                             drop = "ANY"),
          function(x, i, j, ..., drop) {
            if (!missing(drop))
              if (!is.logical(drop) || drop)
                warning ("argument 'drop' ignored. Only drop = FALSE accepted")
            x@.Data <- x@.Data[, j, drop = FALSE]
            x
          })


setMethod("[", signature = c(x = "SmoothLag", i = "numeric", j = "numeric",
                             drop = "ANY"),
          function(x, i, j, ..., drop) {
            if (!missing(drop))
              if (!is.logical(drop) || drop)
                warning ("argument 'drop' ignored. Only drop = FALSE accepted")
            .SmoothLag(x@.Data[i, j, drop = FALSE],
                       basis = x@basis,
                       random = x@random[i, , drop = FALSE],
                       signature = x@signature
                       )
          })


setMethod("[", signature = c(x = "SmoothLag", i = "logical", j = "missing",
                             drop = "ANY"),
          function(x, i, j, ..., drop) {
            if (!missing(drop))
              if (!is.logical(drop) || drop)
                warning ("argument 'drop' ignored. Only drop = FALSE accepted")
            .SmoothLag(x@.Data[i, , drop = FALSE],
                       basis = x@basis,
                       random = x@random[i, , drop = FALSE],
                       signature = x@signature
                       )
          })

setMethod("[", signature = c(x = "SmoothLag", i = "missing", j = "logical",
                             drop = "ANY"),
          function(x, i, j, ..., drop) {
            if (!missing(drop))
              if (!is.logical(drop) || drop)
                warning ("argument 'drop' ignored. Only drop = FALSE accepted")
            x@.Data <- x@.Data[, j, drop = FALSE]
            x
          })

setMethod("[", signature = c(x = "SmoothLag", i = "logical", j = "logical",
                             drop = "ANY"),
          function(x, i, j, ..., drop) {
            if (!missing(drop))
              if (!is.logical(drop) || drop)
                warning ("argument 'drop' ignored. Only drop = FALSE accepted")
            .SmoothLag(x@.Data[i, j, ..., drop = FALSE],
                       basis = x@basis,
                       random = x@random[i, , drop = FALSE],
                       signature = x@signature
                       )
          })

setMethod("[", signature = c(x = "SmoothLag", i = "numeric", j = "logical",
                             drop = "ANY"),
          function(x, i, j, ..., drop) {
            if (!missing(drop))
              if (!is.logical(drop) || drop)
                warning ("argument 'drop' ignored. Only drop = FALSE accepted")
            .SmoothLag(x@.Data[i, j, ..., drop = FALSE],
                       basis = x@basis,
                       random = x@random[i, , drop = FALSE],
                       signature = x@signature
                       )
          })

setMethod("[", signature = c(x = "SmoothLag", i = "logical", j = "numeric",
                             drop = "ANY"),
          function(x, i, j, ..., drop) {
            if (!missing(drop))
              if (!is.logical(drop) || drop)
                warning ("argument 'drop' ignored. Only drop = FALSE accepted")
            .SmoothLag(x@.Data[i, j, ..., drop = FALSE],
                       basis = x@basis,
                       random = x@random[i, , drop = FALSE],
                       signature = x@signature
                       )
          })


## alternative catch-all
setMethod("[", signature = c(x = "SmoothLag", i = "ANY", j = "ANY",
                             drop = "ANY"),
          function(x, i, j, ..., drop)
            stop ("invalid or not yet implemented dimension set")
          )

