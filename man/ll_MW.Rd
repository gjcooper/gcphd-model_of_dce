% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/architectures.R, R/old_loglike.R
\name{ll_MW}
\alias{ll_MW}
\title{Max Winner model}
\usage{
ll_MW(rt, A, b_acc, b_rej, t0, drifts, accept)

ll_MW(rt, A, b_acc, b_rej, t0, drifts, accept)
}
\arguments{
\item{rt}{A vector of response times}

\item{A}{Start point variability}

\item{b_acc}{positive evidence threshold}

\item{b_rej}{negative evidence threshold}

\item{t0}{non decision time parameter}

\item{drifts}{An array of drift rates, 4 x length(rt) with the columns being
vectors of drift rates for accept price, reject price, accept rating,
reject rating respectively}

\item{accept}{Whether we are looking at accept or reject trials}
}
\value{
The log likelihood of the rts for the accept or reject trials given
the provided parameter values

The log likelihood of the rts for the accept or reject trials given
the provided parameter values
}
\description{
A model where both accumulators for Price and Rating must terminate before
either the option is accepted or rejected.

A model where both accumulators for Price and Rating must terminate before
either the option is accepted or rejected.
}
