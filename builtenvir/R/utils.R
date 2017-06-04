

## Non-Exported Misc Utility Functions for builtenvir Package
## -------------------------------------------------------------------



## temporary solution to make.names() until I find a better way
## to handle this

.remake.names <- function(x) {
  paste(deparse(substitute(x)),
        make.unique(colnames(x, do.NULL = FALSE, prefix = ""), sep = ""),
        sep = ""
        )
}




## compute empirical Pr(x > reference) or Pr(|x| > |reference|)
## from a vector of data, x. Take half influence of ties.

sims.pval <- function(x, ref = 0, two.sided = TRUE) {
  x <- x - ref
  mu <- mean(x)
  p <- 1 - (sum(x * sign(mu) > 0) + sum(x == 0) / 2) / length(x)
  if (two.sided)
    p <- 2 * p
  max(min(p, 1), 1 / length(x))
}
