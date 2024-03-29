% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plotting.R
\name{get.manlimits}
\alias{get.manlimits}
\title{Get limts of any parameter considering all spict objects in rep$man}
\usage{
get.manlimits(rep, par, CI = 0.95)
}
\arguments{
\item{rep}{A result report as generated by running fit.spict.}

\item{par}{The name of the parameter to be plotted.}

\item{CI}{Confidence intervals to be calculated, e.g. 0.9 for the 90\%
confidence intervals. By default (CI = 0.95), the 95\% confidence
intervals are estimated.}
}
\value{
plotting limits for all reps in rep$man
}
\description{
Get limts of any parameter considering all spict objects in rep$man
}
