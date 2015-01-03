#' Fit the Stopped Negative Binomial Distribution
#' 
#' Fit a vector of zero's and ones according to the stopped negative
#' binomial distribution with specified 's' and 't' parameters and
#' a beta prior.
#'
#' @param x The vector of 1's and 0's to fit.
#' @param s The ceiling for the snb process.
#' @param t The right barrier for the snb process.
#' @param prior The two element vector giving the parameters for the beta prior.
#' @export
fit_flips <- function(x, s, t, prior=c(0.5, 0.5)) {
  alpha = sum(x) + prior[1]
  beta = sum(x==0) + prior[2] 
  ret = c(alpha, beta, s, t)
  names(ret) = c("shape1", "shape2", "s", "t")
  class(ret) = "flips"
  ret
}

#' @export
summary.flips = function(object, ...) {
  header = c("s", "t", "Mean of p", "Var of p\n")
  cat(paste(sprintf("%8s", header)))
  alpha = object['shape1']
  beta = object['shape2']
  m = alpha / (alpha + beta)
  v = (alpha * beta) / ( (alpha + beta)^2 * (alpha + beta + 1) )
  s =  c(object['s'], object['t'], m, signif(v, 4))
  cat(paste(sprintf("%8s", c(object['s'], object['t'], m, signif(v, 4)))), "\n")
  s =  c(object['s'], object['t'], m, signif(v, 4))
  names(s) = header
  invisible(s)
}

#' @export
fit_snb = function(x, s, t, prior=c(0.5, 0.5), num_sim=10000) {
  snb = fit_flips(x, s, t, prior)
  if (is.null(getDoParName())) registerDoSEQ()
  ps = rbeta(num_sim, snb['shape1'], snb['shape2'])
  ret = sapply(ps, function(p) tail(snb_flips(1, p, s, t), 1))
  list(s_prob = sum(ret), t_prob = sum(ret == 0), n=num_sim)
  ret = c(sum(ret), sum(ret==0))
  names(ret) = c("s", "t")
  class(ret) = "snb"
  ret
}

#' @export
summary.snb = function(object, ...) {
  header = c("Prob s", "Prob t", "n")
  cat(paste(sprintf("%8s", header)), "\n")
  s = object['s']
  t = object['t']
  ret = c(s/(s+t), t / (s+t), s+t)
  cat(paste(sprintf("%8s", ret)), "\n")
  invisible(ret)
}


