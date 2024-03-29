% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/retro.R
\name{retro}
\alias{retro}
\title{Conduct retrospective analysis}
\usage{
retro(rep, nretroyear = 5, reduce_output_size = TRUE, mc.cores = 1)
}
\arguments{
\item{rep}{A valid result from fit.spict.}

\item{nretroyear}{Number of years of data to remove (this is also the total
number of model runs).}

\item{reduce_output_size}{logical, if \code{TRUE} (default) the retro is run
with \code{getReportCovariance} and \code{getJointPrecision} set as
\code{FALSE}}

\item{mc.cores}{Number of cores for \code{parallel::mclapply} function. By
default 1.}
}
\value{
A rep list with the added key retro containing the results of the
    retrospective analysis. Use plotspict.retro() to plot these results.
}
\description{
Conduct retrospective analysis
}
\details{
A retrospective analysis consists of estimating the model with later
    data points removed sequentially one year at a time.
}
\examples{
data(pol)
inp <- pol$albacore
rep <- fit.spict(inp)
rep <- retro(rep, nretroyear=6)
plotspict.retro(rep)
}
