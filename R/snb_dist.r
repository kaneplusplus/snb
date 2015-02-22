require(foreach)
require(ggplot2)
require(reshape2)

R = function(k, p, t) {
  choose(k-1, k-t) * p^(k-t) * (1-p)^t
}

N = function(k, p, s) {
  choose(k-1, s-1) * p^s * (1-p)^(k-s)
}

dsnb_private = function(x, p, s, t, tol=1e-7) {
  a = foreach(k=1:(s+t-1), .combine=c) %do% N(k, p, s)
  b = foreach(k=1:(s+t-1), .combine=c) %do% R(k, p, t)
  d = a + b
  if (abs(sum(a+b) - 1) > tol) stop("Density error")
  inds = which(x %in% 1:length(d))
  ret = rep(0, length(x))
  ret[inds] = d[x[inds]]
  ret
}

dsnb_private_stacked = function(x, p, s, t, tol=1e-7) {
  a = foreach(k=1:(s+t-1), .combine=c) %do% N(k, p, s)
  b = foreach(k=1:(s+t-1), .combine=c) %do% R(k, p, t)
  d = a + b
  if (abs(sum(a+b) - 1) > tol) stop("Density error")
  u = a/sum(d)
  r = b/sum(d)
  ret = foreach (i=x, .combine=rbind) %do% {
    v = c(i, 0, 0)
    if (i %in% 1:length(d)) {
      v[2] = a[i] #u[i] 
      v[3] = b[i] #r[i] 
    }
    v
  }
  colnames(ret) = c("x", "t", "r")
  rownames(ret) = NULL
  ret
}

#' The Stopped Negative Binomial p.m.f. Stack-Plot
#'
#' The stacked plot of the probability mass function for the snb showing
#' the contributions from N (the top barrier) and R (the right barrier).
#' @param x the range of the distribution (defaults to min(s,t):(t+s-1)).
#' @param p the probability of a success on each trial.
#' @param s the top barrier for the snb process.
#' @param t the right barrier for the snb process.
#' @return a plot of the probability mass function.
#' @export
dsnb_stack_plot = function(x, p, s, t) {
  d = as.data.frame(
    dsnb.private.stacked(x, p=p, s=s, t=t))
  d = melt(data=d, id.vars="x") 
  names(d)[names(d) == "variable"] = "barrier"
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
dsnb = function(x, prob, s, t) {
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
  ret = c()
  if (length(x) > 1 && length(prob) > 1) {
    ret = foreach(xx=x, pp=prob, .combine=c) %do% {
      dsnb_private(xx, pp, s, t)
    }
  } else if (length(x) > 1 && length(prob) == 1) {
    ret = foreach(xx=x, .combine=c) %do% {
      dsnb_private(xx, prob, s, t)
    }
  } else if (length(x) == 1 && length(prob) > 1 ) {
    ret = foreach(pp=prob, .combine=c) %do% {
      dsnb_private(x, pp, s, t)
    }
  } else if (length(x) == 1 && length(prob) == 1 ) {
    ret = dsnb_private(x, prob, s, t)
  }
  ret
}

#' @export
rsnb = function(n, prob, s, t) {
  if (length(prob) > 1)
    stop("rsnb prob-parameter must have length 1")
  # Get the distribution function.
  support = min(s,t):(t+1-1)
  ps = dsnb(support, prob, s, t)
  sample(support, n, replace=TRUE, prob=ps)
}

#' @export
snb_flips = function(n, prob, s, t, drop=TRUE) {
  if (length(prob) > 1)
    stop("rsnb prob-parameter must have lenght 1")
  flips = foreach(i=1:n) %do% {
    flip = rbinom(s+t-1, 1, prob=prob)
    path = c(cumsum(flip), sum(flip))
    m = which(path >= s)
    if (length(m) == 0) {
      m = s+t-1
    } else {
      m = min(m)
    }
    r = which(path < 0:(s+t-1)-(t-s+1))
    if (length(r) == 0) {
      r = s+t-1
    } else {
      r = min(r)
    }
    flip[1:min(m, r)]
  }
  if (n == 1 && drop) {
    flips = unlist(flips)
  }
  flips
}

flips_to_zplot_df = function(flips) {
  d = data.frame(k=0:length(flips))
  d$head = c(0, cumsum(flips))
  d$tail= c(0, cumsum(!(flips)))
  d$headEnd = c(d$head[-1], NA)
  d$tailEnd = c(d$tail[-1], NA)
  d
}

#' @export
zplot = function(flips, s, t, show_arrows=TRUE, unif_jitter=0.2, xlab=NULL,
                 ylab=NULL) {
  p = NULL
  if (!is.list(flips)) {
    d =flips_to_zplot_df(flips)
    if (show_arrows) {
      p = ggplot(data=na.omit(d)) + 
        geom_segment(mapping=aes(x=tail, y=head, xend=tailEnd,
                                 yend=headEnd), arrow=arrow()) 
    } else {
      p = ggplot(data=na.omit(d)) + 
        geom_segment(mapping=aes(x=tail, y=head, xend=tailEnd,
                                 yend=headEnd)) 
    }
  } else {
    flip_set = lapply(flips, flips_to_zplot_df)
    for (i in 1:length(flip_set)) {
      flip_set[[i]] = na.omit(flip_set[[i]])
      flip_set[[i]]$num = as.factor(i)
      if (tail(flip_set[[i]]$headEnd, 1) == s) {
        # We hit the top barrier. Jitter on the x values
        flip_set[[i]]$tail= flip_set[[i]]$tail + runif(nrow(flip_set[[i]]),
          -unif_jitter, unif_jitter) 
        flip_set[[i]]$tailEnd = flip_set[[i]]$tailEnd + 
          runif(nrow(flip_set[[i]]), -unif_jitter, unif_jitter)
      } else {
        flip_set[[i]]$head = flip_set[[i]]$head + runif(nrow(flip_set[[i]]), 
          -unif_jitter, unif_jitter)
        flip_set[[i]]$headEnd = flip_set[[i]]$headEnd + 
          runif(nrow(flip_set[[i]]), -unif_jitter, unif_jitter)
      }
      # Make sure that the paths "connect".
      for (j in nrow(flip_set[[i]]):2) {
        flip_set[[i]][j, c("head", "tail")] = 
          flip_set[[i]][j-1, c("headEnd", "tailEnd")]
      }
    }
    d = Reduce(rbind, flip_set)
    if (show_arrows) {
      p = ggplot(data=na.omit(d)) + 
        geom_segment(mapping=aes(x=tail, y=head, xend=tailEnd,
                     yend=headEnd, group=num), arrow=arrow())
    } else {
      p = ggplot(data=na.omit(d)) + 
        geom_segment(mapping=aes(x=tail, y=head, xend=tailEnd,
                     yend=headEnd, group=num))
    }
  }
  p = p+scale_x_continuous(breaks=0:t, limits=c(-unif_jitter, t)) +
    scale_y_continuous(breaks=0:s, limits=c(-unif_jitter, s)) +
    geom_segment(x=0, y=s, xend=t-1, yend=s, color="red") +
    geom_segment(x=t, y=0, xend=t, yend=s-1, color="green")
  if (!is.null(xlab))
    p = p + xlab(xlab)
  if (!is.null(ylab))
    p = p + ylab(ylab)
  p
}

stairs = function(p, xstart, xend) {
  x = c(xstart, rep((xstart+1):xend, each=2))
  y = rep(0:(xend-xstart), each=2)
  y = y[-length(y)]
  for (i in 1:(length(x)-1)) {
    p = p + geom_segment(x=x[i], y=y[i], xend=x[i+1], yend=y[i+1],
      color="green")
  }
  p
}

flips_to_kplot_df = function(flips) {
  d = data.frame(k=0:length(flips))
  d$head = c(0, cumsum(flips))
  d$tail= c(0, cumsum(1-flips))
  d$headEnd = c(d$head[-1], NA)
  d$tailEnd = c(d$tail[-1], NA)
  d$path = c(0, cumsum(flips))
  d$k = 0:(nrow(d)-1)
  d
}

#' @export
kplot = function(flips, s, t) {
  if (!is.list(flips)) {
    d = flips_to_kplot_df(flips)
    p = qplot(k, path, data=d, geom="line") +
      scale_x_continuous(breaks=0:(t+s), limits=c(0, t+s)) +
      scale_y_continuous(breaks=0:s, limits=c(0, s)) +
      geom_segment(x=0, y=s, xend=(t+s-1), yend=s, color="red")
    p = stairs(p, t, s+t-1)
  } else {
    flip_set = lapply(flips, flips_to_kplot_df)
    for (i in 1:length(flip_set)) {
      flip_set[[i]]$num = as.factor(i)
      flip_set[[i]]$k = jitter(flip_set[[i]]$k)
      flip_set[[i]]$k[flip_set[[i]]$k < 0] = 0
    }
    d = Reduce(rbind, flip_set)[,-(4:5)]
    p = qplot(k, path, data=d, geom="path", group=num) +
      scale_x_continuous(breaks=0:(t+s), limits=c(0, t+s)) +
      scale_y_continuous(breaks=0:s, limits=c(0, s)) +
      geom_segment(x=0, y=s, xend=(t+s-1), yend=s, color="red") 
    p = stairs(p, t, s+t-1)
    p
  }
  p
}

#' @export
psnb = function(q, prob, s, t) {
  if (length(prob) > 1)
    stop("psnb prob-parameter may only have length 1")
  support = min(s, t):(t+s-1)
  cdf = c(rep(0, support[1]-1), cumsum(dsnb(support, prob, s, t)))
  qs = floor(q)
  qs[qs < support[1]] = support[1]-1
  qs[qs > support[length(support)]] = support[length(support)]
  cdf[qs]
}

#' @export
qsnb = function(p, prob, s, t) {
  if (length(prob) > 1)
    stop("psnb prob-parameter may only have length 1")
  support = min(s, t):t
  cdf = c(rep(0, support[1]-1), cumsum(dsnb(support, prob, s, t)))
  ret = foreach(pr=p, .combine=c) %do% {
    r = NA
    if (!is.na(pr)) {
      r = which(pr < cdf)[1]
      if (is.na(r))
        r = support[length(support)]
    }
    if (pr > 1 || pr < 0)
      r = NaN
    r
  }
  ret[ret < support[1]-1] = support[1] - 1
  ret
}

