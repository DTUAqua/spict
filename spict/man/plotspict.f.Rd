% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plotting.R
\name{plotspict.f}
\alias{plotspict.f}
\title{Plot estimated fishing mortality.}
\usage{
plotspict.f(
  rep,
  logax = FALSE,
  main = "Absolute fishing mortality",
  ylim = NULL,
  plot.obs = TRUE,
  qlegend = TRUE,
  xlab = "Time",
  ylab = NULL,
  rel.axes = TRUE,
  rel.ci = TRUE,
  stamp = get.version(),
  verbose = TRUE,
  CI = 0.95
)
}
\arguments{
\item{rep}{A result report as generated by running fit.spict.}

\item{logax}{Take log of y-axis? default: FALSE}

\item{main}{Title of plot.}

\item{ylim}{Limits for y-axis.}

\item{plot.obs}{If TRUE observations are plotted.}

\item{qlegend}{If TRUE legend explaining colours of observation data is plotted.}

\item{xlab}{Label of x-axis.}

\item{ylab}{Label of y-axis.}

\item{rel.axes}{Plot secondary y-axis contatning relative level of F.}

\item{rel.ci}{Plot confidence interval for relative level of F.}

\item{stamp}{Stamp plot with this character string.}

\item{verbose}{Should detailed outputs be provided (default: TRUE).}

\item{CI}{Confidence intervals to be calculated, e.g. 0.9 for the 90\%
confidence intervals. By default (CI = 0.95), the 95\% confidence
intervals are estimated.}
}
\value{
Nothing.
}
\description{
Plot estimated fishing mortality.
}
\details{
Plots estimated fishing mortality with Fmsy and associated
    confidence interval.

If no management scenarios are included in \code{rep$man}, the grey vertical
line corresponds to the time of the last observation. If management scenarios
are included in \code{rep$man}, the prediction and confidence intervals of
the base scenario (\code{rep}) are omitted and instead the projections of the
different management scenarios are drawn in different colours. Dotted lines
of the management scenarios reflect the intermediate period, while solid
lines reflect the management period. Additionally, two vertical lines
correspond to the start and end of the management period.
}
\examples{
data(pol)
rep <- fit.spict(pol$albacore)
plotspict.f(rep)

}
