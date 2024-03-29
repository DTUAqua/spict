% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plotting.R
\name{plot2}
\alias{plot2}
\title{Plot summarising spict results (alternative plot composition)}
\usage{
plot2(rep, stamp = get.version(), verbose = TRUE, CI = 0.95, ...)
}
\arguments{
\item{rep}{A result report as generated by running fit.spict.}

\item{stamp}{Stamp plot with this character string.}

\item{verbose}{Should detailed outputs be provided (default: TRUE).}

\item{CI}{Confidence intervals to be calculated, e.g. 0.9 for the 90\%
confidence intervals. By default (CI = 0.95), the 95\% confidence
intervals are estimated.}

\item{...}{additional arguments affecting the summary produced.}
}
\value{
Nothing.
}
\description{
Plot summarising spict results (alternative plot composition)
}
\details{
Create a plot containing the following:
\itemize{
 \item{1. Estimated biomass relative to Bmsy using plotspict.bbmsy().}
 \item{2. Estimated fishing mortality relative to Fmsy using plotspict.ffmsy().}
 \item{3. Observed versus predicted catches using plotspict.catch().}
 \item{4. Estimated F versus estimated B using plotspict.fb().}
}

If no management scenarios are included in \code{rep$man}, the grey vertical
line corresponds to the time of the last observation. If management scenarios
are included in \code{rep$man}, the prediction and confidence intervals of
the base scenario (\code{rep}) are omitted and instead the projections of the
different management scenarios are drawn in different colours. Dotted lines
of the management scenarios reflect the intermediate period, while solid
lines reflect the management period. Additionally, two vertical lines
correspond to the start and end of the management period.

Be aware that potential catch intervals of more a year, e.g. biennial
assessment so that the intermediate period spans two years, or management
period spans two years, are equally split up into annual intervals.

Be aware of the fact that the catches represent intervals, where the length
of the interval is indicated by \code{dtc}, e.g. with \eqn{dtc = 1, C(1990) =
[1990,1990[}. In the plot the catches (and vertical lines) correspond to the
beginning of the catch interval. It might thus seem as if the time of the
vertical lines and the management interval would not align.
}
\examples{
data(pol)
rep <- fit.spict(pol$albacore)
plot2(rep)

}
