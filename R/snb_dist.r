require(foreach)
require(ggplot2)
require(reshape2)

R = function(k, p, t) {
  choose(k-1, k-t) * p^(k-t) * (1-p)^t
}

N = function(k, p, s) {
  choose(k-1, s-1) * p^s * (1-p)^(k-s)
}

Rc = function(k, s, t, shape1, shape2) {
  suppressWarnings({ret = choose(k-1, k-t) * beta(shape1 + k - t, t + shape2) / 
    beta(shape1, shape2)})
  ret[!is.finite(ret)] = 0
  ret
}

Nc = function(k, s, t, shape1, shape2) {
  suppressWarnings({ret = choose(k-1, s-1) * beta(shape1 + s, k - s + shape2) / 
    beta(shape1, shape2)})
  ret[!is.finite(ret)] = 0
  ret
}

dsnb_private = function(x, p, s, t) {
  k=NULL
  a = foreach(k=1:(s+t-1), .combine=c) %do% N(k, p, s)
  b = foreach(k=1:(s+t-1), .combine=c) %do% R(k, p, t)
  d = a + b
  inds = which(x %in% 1:length(d))
  ret = rep(0, length(x))
  ret[inds] = d[x[inds]]
  ret
}

# Remember, shape1 and shape2 are data plus priors.
dsnbc_private = function(x, s, t, shape1, shape2) {
  k = NULL
  a = foreach(k=1:(s+t-1), .combine=c) %do% Nc(k, s, t, shape1, shape2)
  b = foreach(k=1:(s+t-1), .combine=c) %do% Rc(k, s, t, shape1, shape2)
  d = a + b
  inds = which(x %in% 1:length(d))
  ret = rep(0, length(x))
  ret[inds] = d[x[inds]]
  ret
}

#' Stack the distribution by responders and non-responders.
#'
#' Stacked distribution function for the stopped negative binomial distribution.
#' @param x vector of quantiles
#' @param p vector of probabilities
#' @param s the number of heads to stop at.
#' @param t the number of tails to stop at.
#' @export
dsnb_stacked = function(x, p, s, t) {
  k = i = NULL
  a = foreach(k=1:(s+t-1), .combine=c) %do% N(k, p, s)
  b = foreach(k=1:(s+t-1), .combine=c) %do% R(k, p, t)
  ret = foreach (i=x, .combine=rbind) %do% {
    v = c(i, 0, 0)
    if (i %in% 1:max(x)) {
      v[2] = a[i]
      v[3] = b[i]
    }
    v
  }
  colnames(ret) = c("x", "s", "t")
  rownames(ret) = NULL
  ret[is.na(ret[,2]), 2] = 0
  ret[is.na(ret[,3]), 3] = 0
  ret
}

dsnb_private_stacked = dsnb_stacked

dsnbc_private_stacked = function(x, shape1, shape2, s, t) {
  k=i=NULL
  a = foreach(k=1:(s+t-1), .combine=c) %do% Nc(k, s, t, shape1, shape2)
  b = foreach(k=1:(s+t-1), .combine=c) %do% Rc(k, s, t, shape1, shape2)
  d = a + b
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
  colnames(ret) = c("x", "s", "t")
  rownames(ret) = NULL
  ret
}


#' The Stopped Negative Binomial p.m.f. Stack-Plot
#'
#' The stacked plot of the probability mass function for the snb showing
#' the contributions from N (the top barrier) and R (the right barrier).
#' @param p the probability of a success on each coin flip. 
#' @param s the top barrier for the snb process.
#' @param t the right barrier for the snb process.
#' @param x the range of the distribution (defaults to min(s,t):(t+s-1)).
#' @param offset an offset on the domain of the distribution. This is 
#' used when getting the conditional distribution where the domain does 
#' not start at 1.
#' @return a plot of the probability mass function.
#' @export
dsnb_stack_plot = function(p, s, t, x, offset) {
  value = Outcome = k = NULL
  if (missing(x))
    x = min(s,t):(t+s-1)
  d = as.data.frame(
    dsnb_private_stacked(x, p=p, s=s, t=t))
  if (!missing(offset))
    d$x = d$x+offset
  d = melt(data=d, id.vars="x") 
  names(d)[names(d) == "variable"] = "Outcome"
  stop("here")
  ggplot(data=d, aes(x=factor(x), y=value, fill=Outcome)) +
    geom_bar(position="stack", stat="identity") + xlab("k") +
    ylab=("f(k,p,s,t)")
}

#' The Conditional Stopped Negative Binomial Density
#'
#' The conditional stacked snb density function. This function gets
#' the distribution of the stopping time when the binomial process has not 
#' reached one of its endpoints. The success probability is fitted using 
#' fit_flips with specified prior.
#' @param d a sequence of 1's and 0's corresponding to the binomial process.
#' @param s the top barrier for the snb process.
#' @param t the right barrier for the snb process.
#' @param prior the shape parameters of the prior on the success probability.
#' @export
cdsnb = function(d, s, t, prior=c(0.5, 0.5)) {
  fit = fit_flips(d, s, t, prior)
  phat = (fit['shape1'] - 1) / (fit['shape1'] + fit['shape2'] - 2)
  total_range = s+t-1
  cond_range = (length(d):total_range)
  cond_s = s - sum(d)
  cond_t = t - sum(1-d)
  # Use the mode unless we don't have enough data.
  p_hat = if(fit[1] < 1 || fit[2] < 1) {
            fit[1] / (fit[1] + fit[2])
          } else {
            (fit[1] - 1) / (fit[1] + fit[2] - 2)
          }
  ret = as.data.frame(dsnb_stacked(min(cond_s,cond_t):(cond_s+cond_t-1), 
                      p_hat, cond_s, cond_t))
  ret$x = ret$x + length(d)
  ret
}

#' The Conditional Stopped Negative Binomial Density Plot
#' 
#' A plot of the stacked snb density function. The plot shows the distribution
#' of teh stopping time when the binomial process has not reached one of 
#' its endpoints. The success probabilty is fitted using the fit_flips function
#' with specefied prior.
#' @param d a sequence of 1's and 0's corresponding to the binomial process.
#' @param s the top barrier for the snb process.
#' @param t the right barrier for the snb process.
#' @param prior the shape parameters of the prior on the success probability.
#' @importFrom reshape2 melt
#' @importFrom ggplot2 qplot
#' @export
cdsnb_stack_plot = function(d, s, t, prior=c(0.5, 0.5)) {
  value = Outcome = NULL
  x = cdsnb(d, s, t, prior)
  x = melt(data=x, id.vars="x") 
  names(x)[names(x) == "variable"] = "Outcome"
  qplot(x=factor(x), y=value, data=x, fill=Outcome, geom="bar", 
    position="stack", stat="identity", ylab="f(k,p,s,t)", xlab="k")
}

#' Stacked Plot of the Compound Stopped Negative Binomial Density 
#'
#' The stacked plot of the probability mass function for the snb showing
#' the contributions from N (the top barrier) and R (the right barrier).
#' @param d the data, a vector of 0 and 1 values.
#' @param s the top barrier for the snb process.
#' @param t the right barrier for the snb process.
#' @param shape1 the value of the first shape parameter on the prior
#' @param shape2 the value of the second shape parameter on the prior
#' @param x the range of the distribution (defaults to min(s,t):(t+s-1)).
#' @return a plot of the probability mass function.
#' @export
dsnbc_stack_plot = function(d, s, t, shape1=0.5, shape2=0.5,
                            x=min(s,t):(t+s-1)) {
  value = Outcome = NULL
  d = dsnbc_stack(d, s, t, shape1, shape2, x)
  d = melt(data=d, id.vars="x") 
  names(d)[names(d) == "variable"] = "Outcome"
  qplot(x=factor(x), y=value, data=d, fill=Outcome, geom="bar", 
    position="stack", stat="identity", ylab="f(k|x,s,t,alpha,beta)", xlab="k")
}

#' The "Stacked" Compound Negative Binomial Density Function
#' 
#' This function returns the "stacked" density function showing the 
#' contribution from each of the end points to the total mass.
#' @param d the data, a vector of 0 and 1 values.
#' @param s the top barrier for the snb process.
#' @param t the right barrier for the snb process.
#' @param shape1 the value of the first shape parameter on the prior
#' @param shape2 the value of the second shape parameter on the prior
#' @param x the range of the distribution (defaults to min(s,t):(t+s-1)).
#' @return the "stacked" density.
#' @export
dsnbc_stack = function(d, s, t, shape1=0.5, shape2=0.5,
                            x=min(s,t):(t+s-1)) {
  num_heads = sum(d)
  num_flips = length(d)
  as.data.frame(
    dsnbc_private_stacked(x, shape1+num_heads, shape2+num_flips-num_heads,s,t))
}


#' The Stopped Negative Binomial Distribution
#'
#' Density, distribution function, quantile function and random generation
#' for the stopped negative binomial distribution with parameters, 'prob',
#' 's', and 't'.
#' 
#' @rdname snb
#' @aliases psnb dsnb qsnb rsnb
#' @usage dsnb(x, prob, s, t)
#' psnb(q, prob, s, t)
#' qsnb(p, prob, s, t)
#' rsnb(n, prob, s, t)
#' @param x vector of quantiles.
#' @param q vector of quantiles.
#' @param p vector of probabilities.
#' @param n number of observations.
#' @param prob probility of success on each coin flip.
#' @param s the ceiling for the snb process
#' @param t the length of the the process can run for.
#' @return 'dsnb' give the density, 'psnb' give the distribution function
#' function, 'qsnb' gives the quantile function, 'rsnb' generates random
#' deviates.
#' @importFrom foreach foreach %do%
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
    pp=xx=NULL
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

#' Simulate the binomial process
#'
#' Generate coin flip trajectories that stop after either s heads or t tails.
#' @param n the number of trajectories to simulate.
#' @param prob the probability of a head.
#' @param s the number of heads to stop at.
#' @param t the number of tails to stop at.
#' @param drop if TRUE then return the results as a matrix. Otherwise return as a list.
#' @export
snb_flips = function(n, prob, s, t, drop=TRUE) {
  if (length(prob) > 1)
    stop("rsnb prob-parameter must have length 1")
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

#' The Z-Plot for the Binomial Process
#'
#' Visualize the stopped Bernoulli process with horizontal axis counting 
#' successes and vertical axis counting failure.
#'
#' @param flips the sequence of coing flips (1's and 0's) to visualize.
#' Note that this can be a list in which case multiple processes will be 
#' shown.
#' @param s the top barrier for the Bernoulli process.
#' @param t the right barrier for the Bernoulli process.
#' @param show_arrows should arrows be shown in the Bernoullis process path?
#' @param unif_jitter for multiple flip paths, how much jitter to add 
#' (default is 0.2).
#' @param xlab the name of the x axis.
#' @param ylab the name of the y axis.
#' @importFrom ggplot2 ggplot geom_segment scale_y_continuous aes
#' @importFrom grid arrow
#' @examples
#' flips = c(0, 0, 1)
#' zplot(flips, 2, 3)
#' @export
zplot = function(flips, s, t, show_arrows=TRUE, unif_jitter=0.2, xlab=NULL,
                 ylab=NULL) {
  p = tailEnd = headEnd = num = NULL
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

#' The K-Plot for the Binomial Process
#'
#' Visualize the stopped Bernoulli process with a horizontal step axis and a 
#' vertical axis counting the number of successes.
#'
#' @param flips the sequence of coin flips (1's and 0's) to visualize.
#' @param s the top barrier for the Bernoulli process.
#' @param t the right barrier for the Bernoulli process.
#' @importFrom ggplot2 scale_x_continuous scale_y_continuous geom_segment qplot
#' @examples
#' flips = c(0, 0, 1)
#' kplot(flips, 2, 3)
#' @export
kplot = function(flips, s, t) {
  k = path = num = NULL
  if (!is.list(flips)) {
    d = flips_to_kplot_df(flips)
    p = qplot(k, path, data=d, geom="line") +
      scale_x_continuous(breaks=0:(t+s), limits=c(0, t+s)) +
      scale_y_continuous(breaks=0:s, limits=c(0, s)) +
      geom_segment(x=s, y=s, xend=(t+s-1), yend=s, color="red")
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
      geom_segment(x=s, y=s, xend=(t+s-1), yend=s, color="red") 
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
  pr = NULL
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

