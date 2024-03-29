% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/manage.R
\name{add.man.scenario}
\alias{add.man.scenario}
\alias{get.TAC,}
\alias{make.man.inp}
\alias{get.TAC}
\title{Define management scenario}
\usage{
add.man.scenario(
  rep,
  scenarioTitle = "",
  maninterval = NULL,
  maneval = NULL,
  ffac = NULL,
  fabs = NULL,
  cfac = NULL,
  cabs = NULL,
  fractiles = list(catch = 0.5, bbmsy = 0.5, ffmsy = 0.5),
  breakpointB = 0,
  safeguardB = list(limitB = 0, prob = 0.95),
  intermediatePeriodCatch = NULL,
  intermediatePeriodCatchSDFac = 1,
  intermediatePeriodCatchList = NULL,
  ctol = 0.001,
  evalBreakpointB = 0,
  verbose = TRUE,
  dbg = 0,
  mancheck = TRUE
)

get.TAC(
  rep,
  scenarioTitle = "",
  maninterval = NULL,
  maneval = NULL,
  ffac = NULL,
  fabs = NULL,
  cfac = NULL,
  cabs = NULL,
  fractiles = list(catch = 0.5, bbmsy = 0.5, ffmsy = 0.5),
  breakpointB = 0,
  safeguardB = list(limitB = 0, prob = 0.95),
  intermediatePeriodCatch = NULL,
  intermediatePeriodCatchSDFac = 1,
  intermediatePeriodCatchList = NULL,
  ctol = 0.001,
  evalBreakpointB = 0,
  verbose = TRUE,
  dbg = 0,
  mancheck = TRUE
)

make.man.inp(
  rep,
  scenarioTitle = "",
  maninterval = NULL,
  maneval = NULL,
  ffac = NULL,
  fabs = NULL,
  cfac = NULL,
  cabs = NULL,
  fractiles = list(catch = 0.5, bbmsy = 0.5, ffmsy = 0.5),
  breakpointB = 0,
  safeguardB = list(limitB = 0, prob = 0.95),
  intermediatePeriodCatch = NULL,
  intermediatePeriodCatchSDFac = 1,
  intermediatePeriodCatchList = NULL,
  ctol = 0.001,
  evalBreakpointB = 0,
  verbose = TRUE,
  dbg = 0,
  mancheck = TRUE
)
}
\arguments{
\item{rep}{A result report as generated by running \code{fit.spict}.}

\item{scenarioTitle}{Title of scenario (default: \code{'customScenario_X'},
where X is an integer equal to the number of scenarios
with the same name in \code{rep$man} plus 1, e.g. \code{'customScenario_3'}).}

\item{maninterval}{Two floats representing the start and end of the
management period. Example: \code{maninterval = c(2020.25,2021.25)}.
Default: NULL.}

\item{maneval}{Time at which to evaluate model states. Example: \code{maneval
= 2021.25}. Default: NULL.}

\item{ffac}{Factor to multiply current fishing mortality by (default: NULL).}

\item{fabs}{Absolute fishing mortality for management period (default: NULL).}

\item{cfac}{Factor to multiply current catch by (default: NULL). Please refer to
the details for more information.}

\item{cabs}{Absolute catch for the management period (default: NULL).}

\item{fractiles}{List defining the fractiles of the 3 distributions of
'catch', 'bbmsy', and 'ffmsy'. By default (0.5) median is used for all 3
quantities. Please refer to the details for more information.}

\item{breakpointB}{Breakpoints in terms of \eqn{B/B_{MSY}} for the
hockey-stick HCR. By default (0) no breakpoint is assumed. If one value
is provided, F is reduced linearly to zero, if \eqn{B/B_{MSY}} is below
the breakpoint. If two values ara provided, F is reduced linearly to the
lower of the two provided values, if \eqn{B/B_{MSY}} is below the higher
and above the lower value, and F is zero if \eqn{B/B_{MSY}} is below the
lower value. The higher value corresponds to ICES's \eqn{B_{trigger}} and
the lower to ICES's \eqn{B_{lim}}. Note that the breakpoints are
evaluated at the start of the management period. Please refer to the
details for more information.}

\item{safeguardB}{List defining an optional precautionary buffer by means of
a biomass reference level relative to \eqn{B/B_{MSY}} (\code{'limitB'};
default: 0, i.e. deactivating the PA buffer) and the risk aversion
probability (\code{'prob'}; default: 0.95). Please refer to the details
for more information.}

\item{intermediatePeriodCatch}{Catch during intermediate period, e.g. last
year's TAC (default: \code{NULL}). Please refer to the details for more
information.}

\item{intermediatePeriodCatchSDFac}{Factor for the multiplication of the
standard deviation of the catch during the intermediate period (default:
1). Please refer to the details for more information.}

\item{intermediatePeriodCatchList}{List defining catch in the intermediate
period obtaining the elements 'obsC', 'timeC', and 'dtc' (optional
element 'stdevfacC' which is 1 if not provided). Please refer to the
details for more information.}

\item{ctol}{Tolerance of \code{nlminb} when finding F that leads to
provided target catch (via arguments \code{cfac} or \code{cabs})}

\item{evalBreakpointB}{Time for the evaluation of the hockey-stick component
of the HCR: 0 indicating start of the management period and 1 indicating
the end of the management period (default: 0).}

\item{verbose}{Should detailed outputs be provided (default: TRUE).}

\item{dbg}{Debug flag, dbg=1 some output, dbg=2 more output.}

\item{mancheck}{Should the time-dependent objects in \code{inp} be checked
against the management time and corrected if necessary? (Default: TRUE)}
}
\value{
\code{add.man.scenario} returns the input object \code{rep} with the
specified HCR added to the \code{man} list.
\code{get.TAC} returns the total allowable catch (TAC) based on the
specified scenario.
\code{make.man.inp} returns the updated \code{inp} list based on specified HCR.
}
\description{
Define management scenario
}
\details{
\subsection{Default management scenario}{
The default management scenario is fish at \eqn{F_{MSY}}. This is when
\code{ffac}, \code{cfac}, \code{fabs}, \code{cabs} are all \code{NULL}, and
\code{breakpointB} and \code{safeguardB$limitB} are 0. In practice \code{ffac}
is set equal to \eqn{F_{MSY}/F_m}.
}

\subsection{Catch scenarios}{
Management scenarios can be defined based on a desired catch during the
management period. Common examples include scenarios like "increase catch by
25\%", "keep current catch", or "zero catch". The catch can be relative to the
predicted "previous catch", using the multiplier \code{cfac}, or in absolute
terms using \code{cabs} catch in the same units as the input data. By default,
the respective previous catch corresponds to that part of the previous
year which corresponds to the management interval. For example,
if the management period is \eqn{[1991, 1992[}, the
whole catch from the year \eqn{[1990, 1991[} is being used. If the
management period is \eqn{[1991.5, 1991.75[}, the same interval from the
previous year \eqn{[1990.5, 1990.75[} is being used. If the management
period spans several years, e.g. \eqn{[1991, 1993[}, the whole catch from
the previous year \eqn{[1990, 1991[} is being used two times.
}
\subsection{Harvest Control Rules (HCRs)}{

The combination of the arguments "fractiles", "breakpointB", and
"safeguardB" allow the specification of a number of different harvest control
rules:

\itemize{
\item{MSY hockey-stick rule: Fishing at F_{MSY} above a certain biomass reference
    level (here defined as a fraction of B_{MSY} with \code{breakpointB}).
    Below the reference level, fishing is reduced linearly to 0 as suggested in
    ICES (2017).}
\item{MSY (hockey-stick) rule with additional precautionary buffer: As long
    as the probability of the predicted biomass relative to a reference
    biomass level (e.g. 0.3 B_{MSY}, defined by \code{safeguardB$limitB}) is smaller or
    equal to a specified risk aversion probability (e.g. 95\%, defined by
    \code{safeguardB$prob}), fishing at F_{MSY} or following the hockey-stick rule
    (if \code{breakpoint != 0}), otherwise reduce fishing mortality to meet
    specified risk aversion probability (\code{safeguardB$prob}) as introduced in
    ICES (2018).}
\item{By ICES (2019) recommended MSY hockey-stick rule with 35th percentiles:
    Fishing at 35th percentile of F_{MSY} above the 35th percentile of 0.5
    \eqn{B/B_{MSY}} (\code{breakpointB = 0.5}) and 35th percentile of
    linearly reduced F_{MSY} below the 35th percentile of 0.5
    \eqn{B/B_{MSY}}. TAC corresponds to 35th percentile of predicted catch.
    Rule is applied with \code{fractiles = list(catch=0.35, bbmsy=0.35,
    ffmsy=0.35)}, \code{breakpointB = 0.5}, and \code{safeguardB =
    list(limitB = 0, prob = 0.95)}.}}
}
\subsection{Fractiles}{
By default, the median (fractile of 0.5) is used for the stock status (\eqn{B/B_{MSY}},
\eqn{F/F_{MSY}}) and predicted catch distribution. A more precautionary
approach is to used fractiles lower than the median (0.5) to account for the
estimated uncertainty. The arguments of the `fractiles` are:
\itemize{
\item{catch - Fractile of the predicted catch distribution}
\item{bbmsy - Fractile of the \eqn{B/B_{MSY}} distribution}
\item{ffmsy - Fractile of the \eqn{F/F_{MSY}} distribution}}

Note that the fractile for the \eqn{F/F_{MSY}} distribution is 1 minus the
fractile specified. As the current fishing mortality is divided by the value
of this distribution \eqn{F{y+1} = \frac{F_y}{F_y/F_{MSY}}}, a lower
percentile of the \eqn{F/F_{MSY}} distribution is more conservative than a
larger one. This allows a consistent setting of fractiles among the different
quantities.
}

\subsection{Biomass safeguard}{
The argument list "safeguardB" includes:
\itemize{
\item{limitB - Reference level for the evaluation of the predicted biomass
  defined as fraction of \eqn{B/B_{MSY}}. By default (\code{safeguardB$limitB
  == 0}) the PA buffer is not used. Theoretically, any value smaller than 1
  is meaningful, but an ICES recommended value would be 30\%
  \code{safeguardB$limitB = 0.3} (ICES, 2018).}
\item{prob - Risk aversion probability of the predicted biomass relative to
  specified reference level (\code{safeguardB$limitB}) for all rules with PA
  buffer (\code{safeguardB$limitB != 0}). Default: 0.95 as recommended by ICES
  (2018).}}
}
\subsection{Intermediate period assumptions}{
Dependent on the start of the management period (e.g. advice year), there
might be a time lag between the last observation and the start of the
management period, often referred to as the intermediate period. If this is
the case, an assumption about the catch during intermediate time period (e.g.
assessment year) has to be made. Two meaningful assumptions are:

\itemize{
\item{1: }{The catch in the intermediate period is based on the fishing
mortality which is extrapolated from the previous year. This is the default
assumption;}
\item{2: }{The catch in the intermediate period is directly specified. This
could for example be the TAC recommended in the previous year. The catch can
be specified by means of the argument \code{intermediatePeriodCatch}. Be
aware that this catch might correspond to several years or a fraction of a
year depending on the time between the last observation and the start of the
management period. The function \code{\link{man.timeline}} can help
visualising the default or specified intermediate period in your data. The
argument \code{intermediatePeriodCatchSDFac} allows to specify the factor
with which to multiply the standard deviation of the catch (\eqn{\sigma_C})
with. It is thus a measure of the certainty around the catch in the
intermediate period. The argument \code{intermediatePeriodCatchList} allows
to define a list with catches and their intervals. It is a list with the
elements 'obsC', 'timeC', 'dtc' and the optional element 'stdevfacC' (which
is equal to 1 if not provided).}}
}

\subsection{\code{make.man.inp}}{
Internal function that creates the required input list for the specific HCR.
}
}
\examples{
data(pol)
rep <- fit.spict(pol$albacore)

## Fishing at Fmsy
rep <- add.man.scenario(rep)

## MSY hockey-stick rule
rep <- add.man.scenario(rep, breakpointB = 0.5)

## ICES (2019) recommended HCR
rep <- add.man.scenario(rep, fractiles = list(catch=0.35, bbmsy=0.35, ffmsy=0.35), breakpointB=0.5)

## Get the TAC for the ICES (2020) recommended HCR (as used in WKMSYSPICT)
rep <- add.man.scenario(rep, fractiles = list(catch=0.35), breakpointB = c(0.3, 0.5))

## Now `rep` includes 3 management scenarios

## Get the TAC when fishing mortality is equal to Fmsy
get.TAC(rep)

## Get TAC for the MSY hockey-stick rule (only using Btrigger)
get.TAC(rep, breakpointB = 0.5)

## Get TAC for the MSY hockey-stick rule (with Btrigger and Blim)
get.TAC(rep, breakpointB = c(0.3, 0.5))

## Get the TAC for the ICES (2019) recommended HCR
get.TAC(rep, fractiles = list(catch=0.35, bbmsy=0.35, ffmsy=0.35), breakpointB=0.5)

## Get the TAC for the ICES (2020) recommended HCR (as used in WKMSYSPICT)
get.TAC(rep, fractiles = list(catch=0.35), breakpointB = c(0.3, 0.5))
}
\references{
ICES. 2017. Report of the Workshop on the Development of the ICES
approach to providing MSY advice for category 3 and 4 stocks
(WKMSYCat34), 6-10 March 2017, Copenhagen, Denmark. ICES CM 2017/
ACOM:47. 53 pp.

ICES. 2018. Report of the Eighth Workshop on the Development of
Quantitative Assessment Methodologies based on LIFE-history traits,
exploitation characteristics, and other relevant parameters for
data-limited stocks (WKLIFE VIII), 8-12 October 2018, Lisbon,
Portugal. ICES CM 2018/ACOM:40. 172 pp.

ICES.2019. Ninth Workshop on the Development of Quantitative Assessment
Methodologies based on LIFE-history traits, exploitation characteristics,
and other relevant parameters for data-limited stocks (WKLIFE IX).
ICES Scientific Reports. 1:77. 131 pp.http://doi.org/10.17895/ices.pub.5550

ICES 2020. Report of the Ninth Workshop on the Development of
Quantitative Assessment Methodologies based on LIFE-history traits,
exploitation characteristics, and other relevant parameters for
data-limited stocks (WKLIFE X), ICES Scientific Reports. 2:98. 72 pp.
http://doi.org/10.17895/ices.pub.5985
}
