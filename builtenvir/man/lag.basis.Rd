% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/lag.basis.R
\name{lag.basis}
\alias{lag.basis}
\title{Cubic Radial Basis Functions for Distributed Lags}
\usage{
\method{lag}{basis}(lag, center = TRUE, scale = FALSE)
}
\arguments{
\item{lag}{Distributed lag to compute basis of}

\item{center}{Either \code{logical} or \code{numeric} value to indicate
if the lag should be mean centered before computing basis
(\code{center = TRUE}), or else giving the value to center lag at.}

\item{scale}{Either \code{logical} or \code{numeric} value to indicate
if the lag should be standard deviation-scaled before computing basis
(\code{center = TRUE}), or else giving the value to scale lag by.}
}
\value{
\code{\link{LagBasis}} object containing the basis matrix.
}
\description{
Compute cubic radial basis for a given lag set.
}
\examples{
l <- seq(0.1, 10, length.out = 100)
lb <- lag.basis(l)

}
