
#' @title Plot a Fitted Dlm
#'
#' @description
#' Plots point estimates and intervals for distributed lag terms from a
#' fitted \code{\link{dlm}}. Default aesthetic relies on
#' \code{\link[ggplot2]{geom_pointrange}} from the \code{ggplot2} package.
#'
#' @param x a fitted \code{"\link{Dlm-class}"} object
#'
#' @param ... additional arguments passed to \code{\link{summary.Dlm}}
#'
#' @return
#' A \code{ggplot2} graphic object
#'

plot.Dlm <- function(x, y, ...) {
  summ <- summary(x, ...)
  df <- cbind(Lag = x@basis@x, data.frame(summ@lag.coefs))
  names (df)[3:4] <- c("low", "high")
  ggplot(df, aes(x = Lag)) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    geom_pointrange(aes(y = DL.Coefficients, ymin = low, ymax = high),
                    size = 0.2
                    )
  ## geom_line(aes(y = DL.Coefficients)) +
  ## geom_ribbon(aes(ymin = low, ymax = high), alpha = 0.35)
}




