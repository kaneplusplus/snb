require(foreach)
require(ggplot2)
require(grid)
require(reshape2)

R <- function(k, p, t) {
  choose(k-1, k-t) * p^(k-t) * (1-p)^t
}

N <- function(k, p, s) {
  choose(k-1, s-1) * p^s * (1-p)^(k-s)
}

dsnb.private <- function(x, p, s, t) {
  a <- foreach(k=1:(s+t-1), .combine=c) %do% N(k, p, s)
  b <- foreach(k=1:(s+t-1), .combine=c) %do% R(k, p, t)
  d <- a + b
  d <- d/sum(d)
  ret <- 0
  if (x %in% 1:length(d))
    ret <- d[x]
  ret
}

dsnb.private.stacked <- function(x, p, s, t) {
  a <- foreach(k=1:(s+t-1), .combine=c) %do% N(k, p, s)
  b <- foreach(k=1:(s+t-1), .combine=c) %do% R(k, p, t)
  d <- a + b
  u <- a/sum(d)
  r <- b/sum(d)
  ret <- foreach (i=x, .combine=rbind) %do% {
    v <- c(i, 0, 0)
    if (i %in% 1:length(d)) {
      v[2] <- u[i] 
      v[3] <- r[i] 
    }
    v
  }
  colnames(ret) <- c("x", "t", "r")
  rownames(ret) <- NULL
  ret
}

dsnbStackPlot <- function(x, p, s, t) {
  d <- as.data.frame(
    dsnb.private.stacked(x, p=p, s=s, t=t))
  d <- melt(data=d, id.vars="x") 
  names(d)[names(d) == "variable"] <- "border"
  qplot(x=factor(x), y=value, data=d, fill=border, geom="bar", 
    position="stack", stat="identity", ylab="f(k)", xlab="k")
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
  support <- min(s,t):(t+1-1)
  ps <- dsnb( support, prob, s, t)
  sample(support, n, replace=TRUE, prob=ps)
}

#' @export
snbFlips <- function(n, prob, s, t, drop=TRUE) {
  if (length(prob) > 1)
    stop("rsnb prob-parameter must have lenght 1")
  flips <- foreach(i=1:n) %do% {
    flip <- rbinom(s+t-1, 1, prob=prob)
    path <- c(cumsum(flip), sum(flip))
    m <- which(path >= s)
    if (length(m) == 0) {
      m <- s+t-1
    } else {
      m <- min(m)
    }
    r <- which(path < 0:(s+t-1)-(t-s+1))
    if (length(r) == 0) {
      r <- s+t-1
    } else {
      r <- min(r)
    }
    flip[1:min(m, r)]
  }
  if (n == 1 && drop) {
    flips <- unlist(flips)
  }
  flips
}

#' @export
zplot <- function(flips, s, t) {
  d <- data.frame(k=0:length(flips))
  d$head <- c(0, cumsum(flips))
  d$tail<- c(0, cumsum(!(flips)))
  d$headEnd <- c(d$head[-1], NA)
  d$tailEnd <- c(d$tail[-1], NA)

  ggplot(data=na.omit(d)) +
    scale_x_continuous(breaks=0:t, limits=c(0, t)) +
    scale_y_continuous(breaks=0:s, limits=c(0, s)) +
    geom_segment(mapping=aes(x=tail, y=head, xend=tailEnd,
      yend=headEnd), arrow=arrow()) +
    geom_segment(x=0, y=s, xend=t-1, yend=s, color="red") +
    geom_segment(x=t, y=0, xend=t, yend=s-1, color="green")
}

stairs <- function(p, xstart, xend) {
  x <- c(xstart, rep((xstart+1):xend, each=2))
  y <- rep(0:(xend-xstart), each=2)
  y <- y[-length(y)]
  for (i in 1:(length(x)-1)) {
    p <- p + geom_segment(x=x[i], y=y[i], xend=x[i+1], yend=y[i+1],
      color="green")
  }
  p
}

#' @export
kplot <- function(flips, s, t) {
  d <- data.frame(k=0:length(flips))
  d$head <- c(0, cumsum(flips))
  d$tail<- c(0, cumsum(1-flips))
  d$headEnd <- c(d$head[-1], NA)
  d$tailEnd <- c(d$tail[-1], NA)
  d$path <- c(0, cumsum(flips))
  d$k <- 0:(nrow(d)-1)

  p <- qplot(k, path, data=d, geom="line") +
    scale_x_continuous(breaks=0:(t+s), limits=c(0, t+s)) +
    scale_y_continuous(breaks=0:s, limits=c(0, s)) +
    geom_segment(x=0, y=s, xend=(t+s-1), yend=s, color="red")
  stairs(p, t, s+t-1)
}

#' @export
psnb <- function(q, prob, s, t) {
  if (length(prob) > 1)
    stop("psnb prob-parameter may only have length 1")
  support <- min(s, t):(t+s-1)
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


