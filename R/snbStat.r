
# Simulate the snb process n times and return a vector of the process
# values when it stopped.
snb.last <- function(n, s, t, prob) {
  foreach(i=1:n, .combine=c) %do% {
    x <- cumsum(rbinom(t, 1, prob=prob))
    min(s, x[length(x)])
  }
}

#' Fit the Stopped Negative Binomial Distribution
#' 
#' Fit a vector of zero's and ones according to the stopped negative
#' binomial distribution with specified 's' and 't' parameters and
#' a beta prior.
#'
#' @param x The vector of 1's and 0's to fit.
#' @param s The ceiling for the snb process.
#' @param t The maximum number of steps the process can be run.
#' @param prior The two element vector giving the parameters for the beta prior.
#' @export
fitsnb <- function(x, s, t, prior=c(0.5, 0.5), numsim=1000) {
  stop("Not implemented yet")
}
