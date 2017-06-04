
#' @title Cubic Radial Basis Functions for Distributed Lags
#'
#' @description
#' Compute cubic radial basis for a given lag set.
#'
#' @param lag
#' Distributed lag to compute basis of
#'
#' @param center
#' Either \code{logical} or \code{numeric} value to indicate
#'   if the lag should be mean centered before computing basis
#'   (\code{center = TRUE}), or else giving the value to center lag at.
#'
#' @param scale
#' Either \code{logical} or \code{numeric} value to indicate
#'   if the lag should be standard deviation-scaled before computing basis
#'   (\code{center = TRUE}), or else giving the value to scale lag by.
#'
#' @return
#' \code{\link{LagBasis}} object containing the basis matrix.
#'
#' @examples
#' l <- seq(0.1, 10, length.out = 100)
#' lb <- lag.basis(l)
#'

#' @export
lag.basis <- function(lag, center = TRUE, scale = FALSE) {
  ## scale and center lag, if desired
  cntr <- 0
  scl <- 1
  if (center && is.logical(center))
    cntr <- mean(lag)
  else if (center)
    cntr <- center

  if (scale && is.logical(scale))
    scl <- sd(lag)
  else if (scale)
    scl <- scale
  clag <- (lag - cntr) / scl

  ## compute basis functions
  C0 <- cbind("lag-const" = 1, "lag-lin" = clag)
  C1 <- abs(outer(clag, clag, "-"))^3
  M1 <- qr.Q(qr(cbind(C0, C1)))[, -(1:2)]
  S <- svd(t(M1) %*% C1 %*% M1)
  M2.inv.sqrt <- S$v %*% diag(1 / sqrt(S$d)) %*% t(S$u)
  new("LagBasis",
      x = lag,
      x.center = cntr,
      x.scale = scl,
      C0 = C0,
      K1 = C1 %*% M1 %*% M2.inv.sqrt
      )
}
