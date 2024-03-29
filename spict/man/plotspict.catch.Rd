% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plotting.R
\name{plotspict.catch}
\alias{plotspict.catch}
\title{Plot observed catch and predictions.}
\usage{
plotspict.catch(
  rep,
  main = "Catch",
  ylim = NULL,
  qlegend = TRUE,
  lcol = "blue",
  xlab = "Time",
  ylab = NULL,
  stamp = get.version(),
  verbose = TRUE,
  CI = 0.95
)
}
\arguments{
\item{rep}{A result report as generated by running fit.spict.}

\item{main}{Title of plot.}

\item{ylim}{Limits for y-axis.}

\item{qlegend}{If TRUE legend explaining colours of observation data is plotted.}

\item{lcol}{Colour of prediction lines.}

\item{xlab}{Label of x-axis.}

\item{ylab}{Label of y-axis.}

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
Plot observed catch and predictions.
}
\details{
Plots observed catch and predictions using the current F and Fmsy. The plot
also contains the equilibrium catch if the current F is maintained. If no
management scenarios are included in \code{rep$man}, the grey vertical line
corresponds to the time of the last observation.

If management scenarios are included in \code{rep$man}, the prediction and
confidence intervals of the base scenario (\code{rep}) are omitted and
instead the projections of the different management scenarios are drawn in
different colours. Generally, dotted lines of the management scenarios
reflect the intermediate period, while solid lines reflect the management
period. The catch of management period which are longer than 1 year are split
up equally into annual intervals. Two vertical lines correspond to the start
and end of the management period, respectively. However, there are special
cases in which there is only one or no vertical line drawn, the catch
trajectories are missing completely, or the line of the catch trajectory is
solid even in the intermediate period. These cases and their implications on
the annual catch plot are described in the following:

\itemize{

\item{If the management period is shorter than a year, no catch trajectories
 are drawn and there is only one vertical line indicating the start of the
 assessment period.}

\item{If the management timeline differs between the scenarios in
\code{rep$man}, no vertical lines are drawn as they would be at different
times for each scenario.}

\item{If the management period cannot be split equally into annual intervals,
e.g. because it is 1.5 years long, the uneven remains are not displayed, in
this example only the catch representative of one year is displayed.
Additionally, the second vertical line indicating the end of the management
period is omitted.}

\item{If the intermediate period is shorter or longer than a year, e.g. 0.5
or 1.25 years, the lines of the management period start at the time of the
last observation, because the catch in the intermediate period cannot be
aggregated and displayed correctly. Additionally, the first vertical line
indicating the start of the management period is omitted.}

}

All catches in SPiCT represent intervals, where the length of the interval is
indicated by \code{dtc}, e.g. with \eqn{dtc = 1, C(1990) = [1990,1990[}. In
the plot the catches (and vertical lines) correspond to the beginning of the
catch interval. It might thus seem as if the time of the vertical lines and
the management interval would not align.
}
\examples{
data(pol)
rep <- fit.spict(pol$albacore)
plotspict.catch(rep)

}
