---
title: "Guidelines for the stochastic production model in continuous time (SPiCT)"
author: "T.K. Mildenberger, A. Kokkalis, C.W. Berg"
date: "06 September, 2021"
output:
  rmarkdown::html_vignette:
    number_sections: false
    toc: false
  pdf_document:
    number_sections: false
    toc: false
    keep_tex: false
bibliography: bibliography.bib
vignette: >
  %\VignetteIndexEntry{SPiCT Guidelines}
  %\VignetteEngine{knitr::rmarkdown}
  \usepackage[utf8]{inputenc}
---


This is a living document, be sure to check for the latest update
([link](https://github.com/DTUAqua/spict/commits/master)). The SPiCT package
is actively developed, check for the most recent package version and
report problems at ([link](https://github.com/DTUAqua/spict/releases)).

| Document version | SPiCT version | Comments                              |
| :--------------- | :-----------: | ---------------------------:          |
|  v 1.0           | v 1.2.7       | Initial version (2019-09-23)          |
|  v 1.1           | v 1.2.8       | Management section added              |
|  v 1.2           | v 1.3.0       | Adjusted to package update            |
|  v 1.2           | v 1.3.4       | Adjusted after WKMSYSPICT and WKLIFEX |

This document provides specific guidelines for the use of the stochastic
production model in continuous time [SPiCT, @spict] and is divided into three
parts, containing: (i) the main assumptions and data requirements of SPiCT, (ii)
a checklist for the acceptance of a SPiCT assessment, (iii) options for
assessment tuning, and (iv) harvest control rules for SPiCT.

Accompanying sample code that assumes that

- input data is contained in a list called `inp`,
- the results are contained in a list called `fit` after fitting SPiCT to the
input data with `fit <- fit.spict(inp)`, and
- calculating the one-step-ahead residuals with `fit <-
calc.osa.resid(fit)`.


## Main assumptions and input data for SPiCT

- Catch data should be representative of both landings and bycatch. It is also
possible to use landings only, but then the interpretation of the results
changes. If available, seasonal catches should be used as input. Catches are
assumed to be taken over a time interval (e.g. years or quarters), thus the
associated time vector in SPiCT `inp$timeC` should reflect the beginning of each
catch interval (e.g. 2002.25 and 2002.75 for the second and fourth quarter
catches, respectively). Additionally, the vector `inp$dtc` should reflect the
length of each time interval (e.g. 1 for annual and 0.25 for quarterly catches,
respectively).

- Stock size indices should be in terms of biomass (not numbers) and
representative of the part of the stock vulnerable to the commercial fleets, the
so called exploitable stock biomass (ESB). In many cases, the gear selectivity
of the commercial and scientific fleets do not coincide and thus the stock size
indices have to be corrected to exclude individuals that are not represented in
the commercial fleets.

- Biomass indices are assumed to be snapshots at given points in time.
Therefore, the timing of survey indices `inp$timeI` has to be given as decimal
years reflecting the timing of the survey (e.g. 1995.5 for the middle of the
year). The timing of the survey will be matched to the closest model time which
is dependent on `inp$dteuler` (see below). Commercial CPUE indices should be
associated with the midpoint of the interval of the corresponding catches, i.e.
middle of the year if they are based on yearly aggregated catches and effort.

- If input data spans long periods and there is ecological evidence
for potential changes in productivity, it is possible to fit
productivity regime shifts or long-term gradual changes in
productivity with SPiCT [@tvp].

- The Euler discretisation has to be finer than the data, i.e. for yearly data
the Euler time step `inp$dteuler` has to be smaller than 1. A value of 1 changes
the model to a discrete time production model with different model assumptions.
The default value is 1/16 year, which is sufficient in most cases.


## Checklist for the acceptance of a SPiCT assessment

ICES category 3 stocks can be managed using the official advice rule 3.1.1 and
3.1.2 in @msycat34, which requires an accepted SPiCT assessment. Following
points should be considered for the acceptance of a SPiCT assessment.

1. The assessment **converged** (`fit$opt$convergence` equals $0$).

2. All **variance parameters** of the model parameters are **finite**
(`all(is.finite(fit$sd)`) should be `TRUE`).

3. **No violation of model assumptions** based on one-step-ahead
residuals (bias, auto-correlation, normality). This means, that
p-values are insignificant ($> 0.05$), indicated by green titles in
the graphs of `plotspict.diagnostic(fit)`. Slight violations of these
assumptions do not necessarily invalidate model results.

4. **Consistent patterns in the retrospective analysis** (`fit <-
retro(fit)`). This means that there is no tendency of consistent
under- or overestimation of the relative fishing mortality
($\frac{F}{F_{MSY}}$) and relative biomass ($\frac{B}{B_{MSY}}$) in
successive assessment. The retrospective trajectories of those two
quantities should be inside the credible intervals of the base run.

5. **Realistic production curve**. The shape of the production curve
should not be too skewed ($\frac{B_{MSY}}{K}$ should be between 0.1
and 0.9). Low values of $\frac{B_{MSY}}{K}$ allow for an infinite
population growth rate (`calc.bmsyk(fit)`).

6. **High assessment uncertainty** can indicate a lack of contrast in
the input data or violation of the ecological model assumptions. The
main variance parameters (`logsdb`, `logsdc`, `logsdi`, `logsdf`)
should not be unrealistically high. Credible intervals for
$\frac{B}{B_{MSY}}$ and $\frac{F}{F_{MSY}}$ should not span more than
1 order of magnitude (`calc.om(fit)`).

7. **Initial values do not influence the parameter estimates** (`fit
<- check.ini(fit)`). The estimates should be the same for all initial
values (`fit$check.ini$resmat`). Runs which did not converge should
not be considered in this regard.

## Optional model and assessment tuning

**Disclaimer**: An uncertain assessment is okay; the model should not
be tweaked and parameters should not be fixed unless there is
sufficient information and evidence to make such assumptions. Rather
than avoiding uncertainty, it should be accounted for by means of
stochastic harvest control rules [@wklifeix].

* Increase **iterations of optimisation**. The special error code `false
convergence (8)` could indicate that the optimisation exceeded the maximum
number of iterations `max.iter`. `inp$optimiser.control = list(iter.max = 1e5,
eval.max = 1e5)`

* If the catch time series is longer than the survey(s), **shortening
the catch time series** to cover only the period where there is an
available biomass index may help model convergence. `inp =
shorten.inp(inp, 2005, 2018)`

* **Adjust** the prior for the parameter `logn` determining the **shape of the
  production curve**.

    - Alternative priors for `logn` could be based on meta studies, e.g. for all
      species pooled based on @thorson2012: `inp$priors$logn <- c(log(1.478),
      0.6, 1)`

    - Tighter Schaefer prior for `logn` (or any other prior e.g. from meta
    studies): `inp$priors$logn <- c(log(2), 0.5, 1)`

    - Fixing n to resemble the Schaefer production model (or the meta study,
    alternatively): `inp$ini$logn <- log(2); inp$phases$logn <- -1`

Be aware that these modifications might cause smaller credible bands.

* Use a **prior for the initial depletion level** based on available
information, FOR EXAMPLE:

    - If evidence or expert knowledge allows to infer that there was **low or no**
    exploitation before the beginning of the available data: initial depletion
    level could be assumed to be close to the carrying capacity (**e.g.**
    `inp$priors$logbkfrac <- c(log(0.8),0.5,1)`)

    - If evidence or expert knowledge allows to infer that there was **high**
    exploitation before the beginning of the available data: initial depletion
    level could be assumed to be a certain fraction of the carrying capacity
    (**e.g.** `inp$priors$logbkfrac <- c(log(0.2),0.5,1)`)

* If information on the level of uncertainty in the biomass index and/or the
uncertainty of the catch is available prior distribution could be used for the
observation error term of the indices `logsdi` and catches `logsdc`,
respectively. This requires to remove the priors for the ratios of process to
observation errors `logalpha` and `logbeta`.

## Management based on an accepted SPiCT assessment

SPiCT estimates MSY based reference levels, which can be used to calculate
quantities relevant for fisheries management, such as the Total Allowable Catch
(TAC). The TAC can be calculated based on fishing at $F_{MSY}$ or based more
elaborate harvest control rules (HCR) which include additional components, e.g.
a linearly reduction in F if the biomass is below a certain reference level
(e.g. $B_{trigger}$) as the hockey-stick MSY rule recommended by @msycat34.
Three ICES workshops emphasise the importance of the quantification and
consideration of assessment uncertainty in HCRs [@wklifevii; @wklifeviii;
@wklifeix]. They introduce two modifications of the hockey-stick MSY rule that
account for assessment uncertainty by (i) using a percentile lower than 50 (i.e.
median) for relative biomass $\frac{B}{B_{MSY}}$ and catch $C_{y+1}$ and higher
than 50 for relative fishing mortality $\frac{F}{F_{MSY}}$ [@wklifevii], (ii)
adjusting F to meet the condition that the predicted biomass at the end of the
advice year is equal or smaller than a specified risk level [e.g. 5%, @wklifeviii]. 
After comprehensive simulation testing, @wklifeix recommends to
use the MSY hockey-stick rule with 35th percentiles for all three quantities:

$$ C_{y+1} = q_C(p) $$

$$ F_{y+1} = F_y \frac{min(1, q_B(p))}{q_F (100 - p)} $$

where the advised catch (C) for advice year $y + 1$ corresponds to the
predicted catch given the fishing mortality trajectory in the advice year, and
where $F_{y}$ and $F_{y+1}$ are the fishing mortalities at the beginning and the
end of the advice year, respectively. Other components of the equations are
defined as follows:

| Components |            Definition and purpose                         |
| :--------- | :-------------------------------------------------------- |
|  $q_C$     | Function that takes a percentile and returns the corresponding predicted catch $C_{y+1}$ given the fishing mortality trajectory during the advice year y+1, i.e. $q_C = \Phi_{C_{pred} | F = F_y ... F_{y+1}}^{-1}$|
|  $q_B$     |  Function that takes a percentile and returns the corresponding predicted $\frac{B_{y+1}}{MSY B_{trigger}}$ at the beginning of the advice year and $MSY B_{trigger} = \frac{B_{MSY}}{2}$, i.e. $q_B = \Phi_{2\frac{B_{y+1}}{B_{MSY}}}^{-1}$|
|  $q_F$     | Function that takes a percentile and returns the corresponding predicted $\frac{F_y}{F_{MSY}}$ at the beginning of the advice year y+1, i.e. $q_F = \Phi_{\frac{F_y}{F_{MSY}}}^{-1}$|
|  $p$       |  Specific percentile of the respective distributions, e.g. 35 [WKLIFE IX, @wklifeix].|


Theoretically, any percentiles could be considered for the three quantities
($C_{y+1}, \frac{B_{y+1}}{B_{MSY B_{trigger}}}, \frac{F_{y}}{F_{MSY}}$). In the
simulation testing done during WKLIFE IX, the same percentile (p) was used for
biomass and catch; for fishing mortality the percentile was equal to 100 - p.
The results show that across all tested scenarios and harvest control rules,
management with the MSY-35 rule leads to high levels of relative yield while
retaining the risks at low levels. Higher percentiles than 35th show higher
levels of risk while achieving similar levels of relative yield, while lower
percentiles show a decrease in yield with small change in risk [@wklifeix]. It is
important to mention that informative priors (lower standard deviations)
potentially affect parameter estimates and credible intervals of SPiCT
assessments, and therefore affecting the calculated TAC.

TAC based on a SPiCT assessment can be calculated using the `get.TAC()` function
within the spict package in R, where the function arguments allow to use any of
the above described HCRs or combinations thereof. The hockey-stick MSY rule with
the 35th percentile as recommended by @wklifeix is included as scenario 8 (or
'ices') in the `manage` function, but can also be applied independently of the
manage function by `get.TAC(rep, fractiles = list(catch=0.35, bbmsy=0.35,
ffmsy=0.35), breakpointB=0.5)`.


## Category 1 stocks

During WKMSYSPICT (held online 15--19 February 2021), SPiCT was considered as the 
assessment model for several stocks. Stocks with accepted SPiCT assessments were 
raised to category 1. In accordance with basis for ICES adcvice [@ices] stocks 
and the recommendations from WKLIFEX [@wklifeix], the recommended management HCR 
includes a hockey-stick with $B_{trigger} = B_{MSY} / 2$ and $B_{lim} = B_{MSY} / 3$. 
The TAC advice is the 35th percentile of the short term forecast of the catch 
distribution. This is achieved in `spict` using: 

`rep <- add.man.scenario(rep, "ICES MSY", fractiles = list(catch = 0.35), breakpointB = c(1/3, 1/2))`


## References
