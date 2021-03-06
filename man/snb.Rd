% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/snb_dist.r
\name{dsnb}
\alias{dsnb}
\alias{psnb}
\alias{qsnb}
\alias{rsnb}
\title{The Stopped Negative Binomial Distribution}
\usage{
dsnb(x, prob, s, t)
psnb(q, prob, s, t)
qsnb(p, prob, s, t)
rsnb(n, prob, s, t)
}
\arguments{
\item{x}{vector of quantiles.}

\item{prob}{probility of success on each coin flip.}

\item{s}{the ceiling for the snb process}

\item{t}{the length of the the process can run for.}

\item{q}{vector of quantiles.}

\item{p}{vector of probabilities.}

\item{n}{number of observations.}
}
\value{
'dsnb' give the density, 'psnb' give the distribution function
function, 'qsnb' gives the quantile function, 'rsnb' generates random
deviates.
}
\description{
Density, distribution function, quantile function and random generation
for the stopped negative binomial distribution with parameters, 'prob',
's', and 't'.
}

