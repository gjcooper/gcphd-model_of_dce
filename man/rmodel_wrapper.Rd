% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/old_loglike.R, R/wrappers.R
\name{rmodel_wrapper}
\alias{rmodel_wrapper}
\title{Wrapper for individual model log likelihood function}
\usage{
rmodel_wrapper(x, data, model, contaminant_prob = 0.02, min_rt = 0, max_rt = 1)

rmodel_wrapper(x, data, model, contaminant_prob = 0.02, min_rt = 0, max_rt = 1)
}
\arguments{
\item{x}{A named vector containing parameter values to test}

\item{data}{The data for a single subject for which the likelihood should be
calculated}

\item{model}{The model to be wrapped and returned}

\item{contaminant_prob}{The probability used for contaminant process in the
modelling. A contaminant process is just a uniform random response in the
allowable time window.}

\item{min_rt}{The smallest possible response time in the data}

\item{max_rt}{The largest possible response time in the data}
}
\value{
The log of the likelihood for the data under parameter values x

The log of the likelihood for the data under parameter values x
}
\description{
This function performs some common steps such as rearrangeing the data
pulling out parameter items into variable and combining the results of
applying to model to different subsets of the data (accept and reject
responses)

This function performs some common steps such as rearrangeing the data
pulling out parameter items into variable and combining the results of
applying to model to different subsets of the data (accept and reject
responses)
}
