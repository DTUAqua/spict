## ---- echo = FALSE------------------------------------------------------------
knitr::opts_chunk$set(collapse = TRUE, comment = "  ", fig.align = 'center', cache=FALSE)

## ---- eval=FALSE--------------------------------------------------------------
#  devtools::install_github("DTUAqua/spict/spict")

## -----------------------------------------------------------------------------
library(spict)

## -----------------------------------------------------------------------------
data(pol)

## -----------------------------------------------------------------------------
pol$albacore

## -----------------------------------------------------------------------------
inp <- check.inp(pol$albacore)
inp$dtc

## ---- fig.width=5, fig.height=5.5, out.width='0.5\\textwidth', fig.show='hold'----
plotspict.data(pol$albacore)

## ---- fig.width=5, fig.height=5.5, out.width='0.5\\textwidth', fig.show='hold'----
inpshift <- pol$albacore
inpshift$timeC <- inpshift$timeC
inpshift$timeI <- inpshift$timeI + 0.8
plotspict.data(inpshift)

## ---- results='show', message=FALSE, warning=FALSE, fig.width=6.5, fig.height=5.5, fig.show='hold'----
plotspict.ci(pol$albacore)

## ---- results='show', message=FALSE, warning=FALSE, fig.width=6.5, fig.height=5.5, fig.show='hold'----
res <- fit.spict(pol$albacore)

## ---- results='show', message=FALSE, warning=FALSE, fig.width=6.5, fig.height=5.5, fig.show='hold'----
names(res)

## ---- eval=FALSE--------------------------------------------------------------
#  summary(res)

## ---- echo=FALSE,results='show', message=FALSE, warning=FALSE, fig.width=6.5, fig.height=5.5, fig.show='hold'----
capture.output(summary(res))

## ---- results='show', message=FALSE, warning=FALSE, fig.width=6.5, fig.height=5.5, fig.show='hold'----
plot(res)

## ---- echo = FALSE------------------------------------------------------------
knitr::kable(data.frame(## ls(pattern = "plotspict.*", envir = asNamespace("spict")),
  Function = c("**Data**",
               "`plotspict.ci`", "`plotspict.data`",
               "**Estimates**",
               "`plotspict.bbmsy`", "`plotspict.biomass`", "`plotspict.btrend`", "`plotspict.catch`",
               "`plotspict.f`", "`plotspict.fb`", "`plotspict.ffmsy`", "`plotspict.priors`", "`plotspict.production`",
               "`plotspict.season`",
               "**Diagnostics & extras**",
               "`plotspict.diagnostic`", "`plotspict.osar`", "`plotspict.likprof`", "`plotspict.retro`", "`plotspict.retro.fixed`",
               "`plotspict.infl`", "`plotspict.inflsum`", "`plotspict.tc`"),
  Plot = c("",
           "Basic data plotting (see section \\ref{pldat})",
           "Advanced data plotting (see section \\ref{advpldat})",
           "",
           "Relative biomass $B/B_{MSY}$ estimates with uncertainty",
           "Absolute (and relative) biomass estimates with uncertainty",
           "Expected biomass trend",
           "Catch data and estimates",
           "Absolute (and relative) fishing mortality $F$",
           "Kobe plot of relative fishing mortality over biomass estimates",
           "Relative fishing mortality $F/F_{MSY}$",
           "Prior-posterior distribution of all parameters that are estimated using priors",
           "Production over $B/K$",
           "Seasonal pattern of fishing mortality $F$",
           "",
           "OSA residual analysis to evaluate the fit",
           "One-step-ahead residual plots, one for data time-series",
           "Profile likelihood of one or two parameters",
           "Retrospective analysis",
           "Retrospective analysis of fixed effects",
           "Influence statistics of observations",
           "Summary of influence of observations",
            "Time to $B_{MSY}$ under different scenarios about $F$"
           )),
  caption = "Available plotting functions.")

## ----plotspict.biomass, results='show', message=FALSE, warning=FALSE, fig.width=4, fig.height=4, fig.show='hold'----
plotspict.biomass(res)

## ----plotspict.bbmsy, results='show', message=FALSE, warning=FALSE, fig.width=4, fig.height=4, fig.show='hold'----
plotspict.bbmsy(res)

## ----plotspict.f, results='show', message=FALSE, warning=FALSE, fig.width=3, fig.height=3.3, fig.show='hold'----
plotspict.f(res, main='', qlegend=FALSE, rel.axes=FALSE, rel.ci=FALSE)
plotspict.ffmsy(res, main='', qlegend=FALSE)

## ----plotspict.catch, results='show', message=FALSE, warning=FALSE, fig.width=4, fig.height=4, fig.show='hold'----
plotspict.catch(res)

## ----plotspict.fb, results='show', message=FALSE, warning=FALSE, fig.width=4, fig.height=4, fig.show='hold'----
plotspict.fb(res, ylim=c(0, 1.3), xlim=c(0, 300))

## ----plotspict.diagnostic, results='show', message=FALSE, warning=FALSE, fig.width=5.5, fig.height=7, fig.show='hold'----
res <- calc.osa.resid(res)
plotspict.diagnostic(res)

## ----get.par, results='show', message=FALSE, warning=FALSE, fig.width=6.5, fig.height=5.5, fig.show='hold'----
get.par('logBmsy', res)

## ----get.par2, results='show', message=FALSE, warning=FALSE, fig.width=6.5, fig.height=5.5, fig.show='hold'----
get.par('logBmsy', res, exp=TRUE)

## ----list.quantities, results='show', message=FALSE, warning=FALSE, fig.width=6.5, fig.height=5.5, fig.show='hold'----
list.quantities(res)

## ----cov.fixed, results='show', message=FALSE, warning=FALSE------------------
res$cov.fixed

## ----cor.fixed, results='show', message=FALSE, warning=FALSE------------------
cov2cor(res$cov.fixed)

## ----cor.bf, results='show', message=FALSE, warning=FALSE, fig.width=4, fig.height=4, fig.show='hold'----
cov2cor(get.cov(res, 'logBmsy', 'logFmsy'))

