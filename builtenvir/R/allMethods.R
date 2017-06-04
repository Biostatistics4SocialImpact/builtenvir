
library (methods)



## Methods for LagBasis class
## -------------------------------------------------------------------

setGeneric("omega", function(object, ...) standardGeneric("omega"))


#' @title Extract Lag Basis Matrix
#'
#' @description
#' Extract lag basis matrix, \eqn{\Omega}
#'
#' @return
#' A square numeric basis matrix
#'
#' @export
setMethod("omega", signature = "LagBasis",
          function(object, ...)  cbind(object@C0, object@K1)
          )


#' @title Predict New Values for Fitted Lag Basis
#'
#' @description
#' Not yet implemented
#'
setMethod("predict", signature = "LagBasis",
          function(object, ...) stop ("Not yet implemented")
          )





## Methods for Dlm class
## -------------------------------------------------------------------

#' @title Extract Distributed Lag Model Coefficients
#'
#' @description
#' Return DLM coefficients
#'
#' @return
#' A numeric vector of fitted coefficients
#'
#' @export
setMethod("coef", signature = "Dlm",
          function(object) {
            b <- c(scaleMat(object) %*% theta(object))
            names (b) <- names (theta(object))
            return (b)
          })


#' @title Extract Residual Standard Deviation 'Sigma'
#'
#' @description
#' Return observed residual standard deviation
#'
#' @return
#' A positive numeric scalar: the residual standard deviation
#'
#' @export
setMethod("sigma", signature = "Dlm", function(object) object@sigma)


#' @title Extract Log-Likelihood
#'
#' @description
#' Log likelihood of the observed data given the model
#'
#' @return
#' An object of class \code{logLik} containing the log-likelihood of the
#' data given the fitted DLM
#'
#' @export
setMethod("logLik", signature = "Dlm",
          function(object) {
            ll <- object@logLik
            attr(ll, "nall") <- attr(ll, "nobs") <- object@N
            attr(ll, "df") <- object@K$fixed + 2
            ## fixed effects + 2 variance params
            class(ll) <- "logLik"
            return (ll)
          })


#' @title Akaike Information Criterion
#'
#' @description
#' Compute AIC statistic for fitted DLM
#'
#' @return
#' A numeric scalar: AIC
#'
#' @export
setMethod("AIC", signature = "Dlm",
          function(object, ..., k = 2) {
            ll <- logLik(object)
            k <- max(k, attr(ll, "df"))
            c(2 * (k - ll))
          })


#' @title "Bayesian" Information Criterion
#'
#' @description
#' Compute BIC statistic for fitted DLM
#'
#' @return
#' A numeric scalar: BIC
#'
#' @export
setMethod("BIC", signature = "Dlm",
          function(object, ...) {
            ll <- logLik(object)
            k <- max(2, attr(ll, "df"))
            c(log(object@N) * k - 2 * ll)
          })


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
#' @export
setMethod("vcov", signature = "Dlm",
          function(object) {
            S <- scaleMat(object)
            S %*% vcov.theta(object) %*% t(S)
          })


## define Generics
setGeneric("fixef",
           function(object) standardGeneric("fixef")
           )

setGeneric("ranef",
           function(object) standardGeneric("ranef")
           )

setGeneric("se.fixef",
           function(object) standardGeneric("se.fixef")
           )

setGeneric("se.ranef",
           function(object) standardGeneric("se.ranef")
           )

setGeneric("scaleMat",
           function(object, ...) standardGeneric("scaleMat")
           )

setGeneric("theta",
           function(object, ...) standardGeneric("theta")
           )

setGeneric("vcov.theta",
           function(object, ...) standardGeneric("vcov.theta")
           )



#' @title Extract Fixed Effects
#'
#' @description
#' Will return a vector of fixed effects coefficients from the fitted DLM.
#'
#' @return
#' A numeric vector: fixed effect coefficients
#'
#' @export
setMethod("fixef", signature = "Dlm",
          function(object) coef(object)[1:object@K$fixed]
          )


#' @title Extract Random Distributed Lag Effects
#'
#' @description
#' Will return a vector of random effects coefficients
#' from the fitted DLM.
#'
#' @return
#' A numeric vector: random effect coefficients
#'
#' @export
setMethod("ranef", signature = "Dlm",
          function(object) {
            ran.index <- object@K$fixed + 1
            coef(object)[ran.index:object@K$total]
          })


#' @title Extract Fixed Effect Standard Errors
#'
#' @description
#' Will return a vector of standard errors for each fixed effect coefficient.
#'
#' @return
#' A numeric vector: fixed effects coefficients standard errors
#'
#' @export
setMethod("se.fixef", signature = "Dlm",
          function(object)
            sqrt(diag(vcov(object))[1:object@K$fixed])
          )


#' @title Extract Random Effect (Distributed Lag) Standard Errors
#'
#' @description
#' Will return the random effects standard error. Note that this is a
#' single value in the non-heirarchical formulation of this model.
#'
#' @return
#' A numeric scalar: the random effects (distributed lag) coefficient
#' standard error
#'
#' @export
setMethod("se.ranef", signature = "Dlm", function(object) object@tau)



#' @title Extract Lag Coefficient Scale Matrix
#'
#' @description
#' Return lag coefficient scale matrix, S, such that the distributed lag
#' coefficients fit by the model are obtained via the transformation
#' \eqn{\beta = S \theta}.
#'
#' @return
#' A numeric matrix: the block diagonal scale matrix, S
#'
#' @export
setMethod("scaleMat", signature = "Dlm",
          function(object, ...) {
            k <- object@K$total
            lag.index <- object@K$covariates + 1
            S <- diag(k)
            S[lag.index:k, lag.index:k] <- omega(object@basis)
            return (S)
          })



setMethod("theta", signature = "Dlm", function(object, ...) object@coefficients)

setMethod("vcov.theta", signature = "Dlm", function(object, ...) object@vcov)


#' @title Extract Model Fitted Values
#'
#' @description
#' Will return predicted \code{y} values from a fitted DLM.
#'
#' @return
#' A numeric vector: the model fitted values
#'
#' @export
setMethod("fitted", signature = "Dlm",
          function(object) {
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
              }) }) })


#' @title Extract Model Residuals
#'
#' @description
#' Will return the observed residuals from the model.
#'
#' @return
#' A numeric vector: the model residuals
#'
#' @export
setMethod("residuals", signature = "Dlm",
          function(object) {
            tryCatch ({
              object@residuals
            },
            error = function(e) {
              tryCatch ({
                residuals(object@lme)
              },
              error = function(e) {
                warning ("object not properly initialized")
                NA_real_
              }) }) })



#' @title Plot a Fitted DLM
#'
#' @description
#' Coming soon!
#'
#' @export
setMethod("plot", signature = "Dlm",
          function(x, y, ...) {
            summ <- summary(x, ...)
            lag <- x@basis@x
            df <- cbind(Lag = lag, data.frame(summ@lag.coefs))
            names (df)[3:4] <- c("low", "high")
            ggplot(df, aes(x = Lag)) +
              geom_hline(yintercept = 0, linetype = "dashed") +
              geom_line(aes(y = DL.Coefficients)) +
              geom_ribbon(aes(ymin = low, ymax = high), alpha = 0.35)
          })


#' @title Generate New Predictions from a Fitted DLM
#'
#' @description
#' Not yet implemented
#'
setMethod("predict", signature = "Dlm",
          function(object, ...) {
            stop ("Not yet implemented")
          })


#' @title Display Results of a Fitted DLM
#'
#' @description
#' Print basic output from a fitted DLM object. Output structured after
#' \code{nlme::lme.print()}
#'
#' @export
setMethod("show", signature = "Dlm",
          function(object) {
            tryCatch ({
              cat ("Call: "); print (object@call)
              method <- if (object@method == "frequentist") {
                  object@lme$method
                } else if (object@method == "bayesian") {
                  "MCMC"
                } else {
                  "???"
                }
              if (method != "???")
                cat ("(", object@method, ") distributed lag model fit via ",
                     method, "\n", sep = ""
                     )
            }, error = function(e) NULL
            )
            semat <- cbind("Random Effects" = se.ranef(object),
                           "Residuals" = sigma(object)
                           )
            rownames (semat) <- "Std Error"
            cat ("Number of Observations:", object@N, "\n")
            cat ("Log-Likelihood:", logLik(object), "\n")
            cat ("Fixed Effects:\n")
            printCoefmat(t(fixef(object)), has.Pvalue = FALSE); cat ("\n")
            printCoefmat(semat, has.Pvalue = FALSE)
          })




## Methods for Class FreqDlm
## -------------------------------------------------------------------

#' @title Summarize a \code{Dlm} Object
#'
#' @description
#' Summarize a fitted DLM object. Details coming soon!
#'
#' @return
#' A '\code{SummaryDlm}' object with information about the fitted model.
#'
#' @seealso \code{\link{SummaryDlm}}
#'
#' @name Dlm_Summary
#'
#' @export
setMethod("summary", signature = "FreqDlm",
          function(object, level = 0.95, ...) {
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
          })



## Methods for Class BayesDlm
## -------------------------------------------------------------------

#' @describeIn Dlm_Summary
#'
#' @description Compute and store various statistical summaries
#'   in a \code{SummaryDlm} object
#'
#' @export
setMethod("summary", signature = "BayesDlm",
          function(object, level = 0.95, ...) {
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
                            "Pr(|Est| > 0)" = p.vals
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
          })


## define generics
setGeneric("coef.sims", function(object, ...) standardGeneric("coef.sims"))


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
#' @export
setMethod("coef.sims", signature = "BayesDlm",
          function(object, ...) {
            object@sims$coefficients %*% t(scaleMat(object))
          })


## Methods for Class SummaryDlm
## -------------------------------------------------------------------

#' @describeIn Dlm_Summary
#'
#' @description Print contents of a \code{SummaryDlm} object
#'
#' @export
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
            printCoefmat(cors[-1, -ncol(cors)], has.Pvalue = FALSE,
                         na.print = ""
                         )
            cat ("\n")
          })


