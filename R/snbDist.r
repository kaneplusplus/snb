require(foreach)

dsnb.private <- function(x, p, s, t) {
  y <- choose(x-1, s-1) * p^s * (1-p)^(x-3)
  w <- which(x == t)
  y[w] <- y[w] + sum(dbinom(0:(s-1), t, p))
  y[x > t | x < s] <- 0
  y

}

#' The Stopped Negative Binomial Distribution
#'
#' Density, distribution function, quantile function and random generation
#' for the stopped negative binomial distribution with parameters, 'prob',
#' 's', and 't'.
#' 
#' @rdname snb
#' @aliases psnb snb qsnb rsnb
#' @param x,q vector of quantiles.
#' @param prob probility of success on each trial. 
#' @param s the ceiling for the snb process
#' @param t the length of the the process can run for.
#' @return 'dsnb' give the density, 'psnb' give the distribution function
#' function, 'qsnb' gives the quantile function, 'rsnb' generates random
#' deviates.
#' @export
dsnb <- function(x, prob, s, t) {
  if (length(s) != 1)
    stop("dsnb s-parameter may only have length 1")
  if (length(t) != 1)
    stop("dsnb t-parameter may only have length 1")
  if (s < 1) 
    stop("dsnb s-parameter must be at least 1")
  if (t < 1) 
    stop("dsnb t-parameter must be at least 1")
  if (any(prob > 1) || any(prob < 0))
    stop("dsnb prob-parameter must be between zero and one inclusive")
  ret <- c()
  if (length(x) > 1 && length(prob) > 1) {
    ret <- foreach(xx=x, pp=prob, .combine=c) %do% {
      dsnb.private(xx, pp, s, t)
    }
  } else if (length(x) > 1 && length(prob) == 1) {
    ret <- foreach(xx=x, .combine=c) %do% {
      dsnb.private(xx, prob, s, t)
    }
  } else if (length(x) == 1 && length(prob) > 1 ) {
    ret <- foreach(pp=prob, .combine=c) %do% {
      dsnb.private(x, pp, s, t)
    }
  } else if (length(x) == 1 && length(prob) == 1 ) {
    ret <- dsnb.private(x, prob, s, t)
  }
  ret
}

#' @export
rsnb <- function(n, prob, s, t) {
  if (length(prob) > 1)
    stop("rsnb prob-parameter must have lenght 1")
  # Get the distribution function.
  support <- min(s,t):t
  ps <- dsnb( support, prob, s, t)
  sample(support, n, replace=TRUE, prob=prob)
}

#' @export
psnb <- function(q, prob, s, t) {
  if (length(prob) > 1)
    stop("psnb prob-parameter may only have length 1")
  support <- min(s, t):t
  cdf <- c(rep(0, support[1]-1), cumsum(dsnb(support, prob, s, t)))
  qs <- floor(q)
  qs[qs < support[1]] <- support[1]-1
  qs[qs > support[length(support)]] <- support[length(support)]
  cdf[qs]
}

#' @export
qsnb <- function(p, prob, s, t) {
  if (length(prob) > 1)
    stop("psnb prob-parameter may only have length 1")
  support <- min(s, t):t
  cdf <- c(rep(0, support[1]-1), cumsum(dsnb(support, prob, s, t)))
  ret <- foreach(pr=p, .combine=c) %do% {
    r <- NA
    if (!is.na(pr)) {
      r <- which(pr < cdf)[1]
      if (is.na(r))
        r <- support[length(support)]
    }
    if (pr > 1 || pr < 0)
      r <- NaN
    r
  }
  ret[ret < support[1]-1] <- support[1] - 1
  ret
}


