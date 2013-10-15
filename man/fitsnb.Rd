\name{fitsnb}
\alias{fitsnb}
\title{Fit the Stopped Negative Binomial Distribution}
\usage{
  fitsnb(x, s, t, prior = c(0.5, 0.5), numsim = 1000)
}
\arguments{
  \item{x}{The vector of 1's and 0's to fit.}

  \item{s}{The ceiling for the snb process.}

  \item{t}{The maximum number of steps the process can be
  run.}

  \item{prior}{The two element vector giving the parameters
  for the beta prior.}
}
\description{
  Fit a vector of zero's and ones according to the stopped
  negative binomial distribution with specified 's' and 't'
  parameters and a beta prior.
}

