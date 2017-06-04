
library (methods)



## LagBasis
## -------------------------------------------------------------------
#' @title Create and Manipulate Lag Basis Functions
#'
#' @description
#' S4 class object to store and query components of
#' lag basis functions. User interface for creating this class can be found in
#' \code{\link{lag.basis}}.
#'
#'
#' @slot x
#'   original lag data
#'
#' @slot x.center
#'   store the value the lag data was centered to. In theory
#'   this is useful for the non-yet-implemented \code{predict} method
#'
#' @slot x.scale
#'   store the value the lag data was scaled by. Again, should be
#'   useful primarily for the \code{predict} method in the future
#'
#' @slot C0
#'   C_0 part of basis matrix
#'
#' @slot K1
#'   K_1 part of basis matrix
#'
#'
#' @references Baek J, et al. (2016) Epidemiology 27(1):116-24.
#' (\href{https://www.ncbi.nlm.nih.gov/pubmed/26414942}{PubMed})
#'
#' @name LagBasis
#'
#' @export
.LagBasis <- setClass("LagBasis",
         slots = c(
           x  = "numeric",
           x.center = "numeric",
           x.scale = "numeric",
           C0 = "matrix",
           K1 = "matrix"
         ),
         prototype = list(
           x = NA_real_,
           x.center = 0,
           x.scale = 1,
           C0 = matrix(),
           K1 = matrix()
         ))



## Dlm
## -------------------------------------------------------------------

#' @title S4 Distributed Lag Models
#'
#' @description
#' Parent class for DLM object representation. Users will interact mainly
#' with daughter classes, \code{\link{FreqDlm}} and \code{\link{BayesDlm}}.
#'
#'
#' @slot method
#'   character identifier for method (frequentist vs. Bayesian)
#'   used to fit the DLM
#'
#' @slot call
#'   language object containing matching function call
#'
#' @slot coefficients
#'   DLM coefficients
#'
#' @slot sigma
#'   residual standard deviation
#'
#' @slot tau
#'   random effects standard deviation
#'
#' @slot logLik
#'   log likelihood
#'
#' @slot vcov
#'   fitted variance-covariance matrix
#'
#' @slot basis
#'   \code{\link{LagBasis}} object containing lag basis representation
#'
#' @slot N
#'   Number of data points the model was fit to
#'
#' @slot K
#'   list of coefficient type lengths. See 'Details.'
#'
#' @details
#' Coming soon!
#'
#'
#' @references
#' Baek J, et al. (2016) Epidemiology 27(1):116-24.
#' (\href{https://www.ncbi.nlm.nih.gov/pubmed/26414942}{PubMed})
#'
#' @seealso \code{\link{dlm}}, \code{\link{LagBasis}}
#'
#' @name Dlm
#'
#' @export
.Dlm <- setClass("Dlm",
         slots = c(
           method = "character",
           call = "call",
           coefficients = "numeric",
           sigma = "numeric",
           tau = "numeric",
           logLik = "numeric",
           vcov = "matrix",
           basis = "LagBasis",
           N = "numeric",
           K = "list"
         ),
         prototype = list(
           method = "",
           coefficients = NA_real_,
           sigma = NA_real_,
           tau = NA_real_,
           logLik = NA_real_,
           vcov = matrix(),
           basis = .LagBasis(),
           N = NA_real_,
           K = list(covariates = NA, fixed = NA,
                    lag = NA, random = NA, total = NA
                    )
         ))



## FreqDlm
## -------------------------------------------------------------------


setOldClass("lme")

#' @title Classical Distributed Lag Model
#'
#' @description
#' Parameters and summaries from a classically fit DLM (freqentist;
#' model fit using \code{\link{nlme::lme}}). Inherits from
#' \code{\link{Dlm}}.
#'
#'
#' @slot lme
#' output of fitted \code{lme} object
#'
#'
#' @references Baek J, et al. (2016) Epidemiology 27(1):116-24.
#' (\href{https://www.ncbi.nlm.nih.gov/pubmed/26414942}{PubMed})
#'
#' @seealso \code{\link{dlm}}, \code{\link{lme}}
#'
#' @name FreqDlm
#'
#' @export
.FreqDlm <- setClass("FreqDlm",
         slots = c(lme = "lme"),
         contains = "Dlm"
         )




## BayesDlm
## -------------------------------------------------------------------
#' @title Bayesian Distributed Lag Model
#'
#' @description
#' Parameters and summaries from a Bayesian DLM. Inherits from
#' \code{\link{Dlm}}.
#'
#'
#' @slot fitted
#'   fitted values; conditional expectation of observed \code{y}
#'
#' @slot DIC
#'   Deviance Information Criterion
#'
#' @slot residuals
#'   observed residuals
#'
#' @slot sims
#'   list of posterior simulation draws of model parameters
#'
#'
#' @references Baek J, et al. (2016) Epidemiology 27(1):116-24.
#' (\href{https://www.ncbi.nlm.nih.gov/pubmed/26414942}{PubMed})
#'
#' @seealso \code{\link{dlm}}
#'
#' @name BayesDlm
#'
#' @export
.BayesDlm <- setClass("BayesDlm",
         slots = c(
           fitted = "numeric",
           DIC = "numeric",
           residuals = "numeric",
           sims = "list"
         ),
         prototype = list(
           fitted = NA_real_,
           DIC = NA_real_,
           residuals = NA_real_,
           sims = list()
         ),
         contains = "Dlm"
         )





## SummaryDlm
## -------------------------------------------------------------------
#' @title Summarize a Distributed Lag Model
#'
#' @description
#' Compute and store various statistical summaries of a fitted distributed
#' lag model object
#'
#'
#' @slot method
#'   Character vector with details of the model fitting procedure
#'
#' @slot call
#'   Matching call from \code{dlm}
#'
#' @slot coefficients
#'   Numeric vector of distributed lag coefficients
#'
#' @slot likelihood.stats
#'   Named numeric vector of likelihood based statistical summaries
#'
#' @slot variances
#'   Named numeric vector with estimaed random effects and residual variance
#'   parameters
#'
#' @slot tTable
#'   t-statistic table summaries of fixed effects with (approximate)
#'   significance
#'
#' @slot cor.fixef
#'   Correlation matrix of fixed effects
#'
#' @slot residuals
#'   Numeric vector with the residuals from the fitted model
#'
#' @slot N
#'   Number of data points the model was fit to
#'
#' @slot K
#'   List of numbers of various covariate types
#'
#'
#' @name SummaryDlm
#'
#' @export
.SummaryDlm <- setClass("SummaryDlm",
                        slots = c(
                          method = "character",
                          call = "call",
                          coefficients = "numeric",
                          likelihood.stats = "numeric",
                          variances = "numeric",
                          tTable = "matrix",
                          lag.coefs = "matrix",
                          cor.fixef = "matrix",
                          residuals = "numeric",
                          N = "numeric",
                          K = "list"
                        ),
                        prototype = list(
                          method = c("", ""),
                          likelihood.stats = NA_real_,
                          variances = NA_real_,
                          tTable = matrix(),
                          lag.coefs = matrix(),
                          cor.fixef = matrix(),
                          residuals = NA_real_,
                          N = NA_real_,
                          K = list(covariates = NA, fixed = NA,
                                   lag = NA, random = NA, total = NA
                                   )
                          ))

