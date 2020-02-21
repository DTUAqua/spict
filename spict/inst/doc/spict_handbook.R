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
inpshift$timeC <- inpshift$timeC + 0.3
inpshift$timeI <- inpshift$timeI + 0.8
plotspict.data(inpshift)

## ---- results='show', message=FALSE, warning=FALSE, fig.width=6.5, fig.height=5.5, fig.show='hold'----
plotspict.ci(pol$albacore)

## ---- results='show', message=FALSE, warning=FALSE, fig.width=6.5, fig.height=5.5, fig.show='hold'----
res <- fit.spict(pol$albacore)

## ---- results='show', message=FALSE, warning=FALSE, fig.width=6.5, fig.height=5.5, fig.show='hold'----
names(res)

## ---- results='show', message=FALSE, warning=FALSE, fig.width=6.5, fig.height=5.5, fig.show='hold'----
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
               "`plotspict.diagnostic`", "`plotspict.osar`", "`plotspict.likprof`", "`plotspict.retro`",
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
           "Influence statistics of observations",
           "Summary of influence of observations",
            "Time to $B_{MSY}$ under different scenarios about $F$"
           )),
  caption = "Available plotting functions.")

## ---- results='show', message=FALSE, warning=FALSE, fig.width=4, fig.height=4, fig.show='hold'----
plotspict.biomass(res)

## ---- results='show', message=FALSE, warning=FALSE, fig.width=4, fig.height=4, fig.show='hold'----
plotspict.bbmsy(res)

## ---- results='show', message=FALSE, warning=FALSE, fig.width=3, fig.height=3.3, fig.show='hold'----
plotspict.f(res, main='', qlegend=FALSE, rel.axes=FALSE, rel.ci=FALSE)
plotspict.ffmsy(res, main='', qlegend=FALSE)

## ---- results='show', message=FALSE, warning=FALSE, fig.width=4, fig.height=4, fig.show='hold'----
plotspict.catch(res)

## ---- results='show', message=FALSE, warning=FALSE, fig.width=4, fig.height=4, fig.show='hold'----
plotspict.fb(res, ylim=c(0, 1.3), xlim=c(0, 300))

## ---- results='show', message=FALSE, warning=FALSE, fig.width=5.5, fig.height=7, fig.show='hold'----
res <- calc.osa.resid(res)
plotspict.diagnostic(res)

## ---- results='show', message=FALSE, warning=FALSE, fig.width=6.5, fig.height=5.5, fig.show='hold'----
get.par('logBmsy', res)

## ---- results='show', message=FALSE, warning=FALSE, fig.width=6.5, fig.height=5.5, fig.show='hold'----
get.par('logBmsy', res, exp=TRUE)

## ---- results='show', message=FALSE, warning=FALSE, fig.width=6.5, fig.height=5.5, fig.show='hold'----
list.quantities(res)

## ---- results='show', message=FALSE, warning=FALSE, fig.width=4, fig.height=4, fig.show='hold'----
res$cov.fixed

## ---- results='show', message=FALSE, warning=FALSE, fig.width=4, fig.height=4, fig.show='hold'----
cov2cor(res$cov.fixed)

## ---- results='show', message=FALSE, warning=FALSE, fig.width=4, fig.height=4, fig.show='hold'----
cov2cor(get.cov(res, 'logBmsy', 'logFmsy'))

## ---- results='show', message=FALSE, warning=FALSE, fig.width=5, fig.height=5, fig.show='hold'----
## res <- fit.spict(pol$albacore)
res <- retro(res, nretroyear = 4)
plotspict.retro(res)

## ---- results='show', message=FALSE, warning=FALSE, fig.width=6, fig.height=5, fig.show='hold'----
set.seed(123)
inp <- list(timeC=pol$albacore$timeC, obsC=pol$albacore$obsC)
inp$timeI <- list(pol$albacore$timeI, pol$albacore$timeI[10:23]+0.25)
inp$obsI <- list()
inp$obsI[[1]] <- pol$albacore$obsI * exp(rnorm(23, sd=0.1)) # Index 1
inp$obsI[[2]] <- 10*pol$albacore$obsI[10:23] # Index 2
res <- fit.spict(inp)
sumspict.parest(res)
plotspict.biomass(res)

## ---- results='show', message=FALSE, warning=FALSE, fig.width=6, fig.height=5, fig.show='hold'----
inpeff <- list(timeC=pol$albacore$timeC, obsC=pol$albacore$obsC,
               timeE=pol$albacore$timeC, obsE=pol$albacore$obsC/pol$albacore$obsI)
repeff <- fit.spict(inpeff)
sumspict.parest(repeff)
par(mfrow=c(2, 2))
plotspict.bbmsy(repeff)
plotspict.ffmsy(repeff, qlegend=FALSE)
plotspict.catch(repeff, qlegend=FALSE)
plotspict.fb(repeff)

## ---- results='show', message=FALSE, warning=FALSE, fig.width=5.5, fig.height=6.5, fig.show='hold'----
inp <- pol$albacore
res1 <- fit.spict(inp)
inp$stdevfacC <- rep(1, length(inp$obsC))
inp$stdevfacC[1:10] <- 5
res2 <- fit.spict(inp)
par(mfrow=c(2, 1))
plotspict.catch(res1, main='No scaling')
plotspict.catch(res2, main='With scaling', qlegend=FALSE)

## ---- results='show', message=FALSE, warning=FALSE, fig.width=6.5, fig.height=5.5, fig.show='hold'----
inp <- check.inp(pol$albacore)
sim <- sim.spict(inp)
plotspict.data(sim)

## ---- results='show', message=FALSE, warning=FALSE, fig.width=6.5, fig.height=5.5, fig.show='hold'----
inp <- list(ini=list(logK=log(100), logm=log(10), logq=log(1)))
sim <- sim.spict(inp, nobs=50)
plotspict.data(sim)

## ---- results='show', message=FALSE, warning=FALSE, fig.width=6.5, fig.height=5.5, fig.show='hold'----
set.seed(31415926)
inp <- list(ini=list(logK=log(100), logm=log(10), logq=log(1),
                     logbkfrac=log(1), logF0=log(0.3), logsdc=log(0.1),
                     logsdf=log(0.3)))
sim <- sim.spict(inp, nobs=30)
res <- fit.spict(sim)
sumspict.parest(res)
par(mfrow=c(2, 2))
plotspict.biomass(res)
plotspict.f(res, qlegend=FALSE)
plotspict.catch(res, qlegend=FALSE)
plotspict.fb(res)

## ---- results='show', message=FALSE, warning=FALSE, fig.width=6, fig.height=5, fig.show='hold'----
set.seed(1234)
inp <- list(nseasons=4, splineorder=3)
inp$timeC <- seq(0, 30-1/inp$nseasons, by=1/inp$nseasons)
inp$timeI <- seq(0, 30-1/inp$nseasons, by=1/inp$nseasons)
inp$ini <- list(logK=log(100), logm=log(20), logq=log(1),
                logbkfrac=log(1), logsdf=log(0.4), logF0=log(0.5),
                logphi=log(c(0.05, 0.1, 1.8)))
seasonsim <- sim.spict(inp)
plotspict.data(seasonsim)

## ---- results='show', message=FALSE, warning=FALSE, fig.width=5.5, fig.height=4.5, fig.show='hold'----
set.seed(432)
inp <- list(nseasons=4, seasontype=2)
inp$timeC <- seq(0, 30-1/inp$nseasons, by=1/inp$nseasons)
inp$timeI <- seq(0, 30-1/inp$nseasons, by=1/inp$nseasons)
inp$ini <- list(logK=log(100), logm=log(20), logq=log(1),
                logbkfrac=log(1), logsdf=log(0.4), logF0=log(0.5))
seasonsim2 <- sim.spict(inp)
plotspict.data(seasonsim2)

## ---- results='show', message=FALSE, warning=FALSE, fig.width=5.5, fig.height=4.5, fig.show='hold'----
seasonres <- fit.spict(seasonsim)
plotspict.biomass(seasonres)
plotspict.f(seasonres, qlegend=FALSE)
plotspict.season(seasonres)

## ---- results='show', message=FALSE, warning=FALSE, fig.width=5.5, fig.height=4.5, fig.show='hold'----
seasonres2 <- fit.spict(seasonsim2)
sumspict.parest(seasonres2)
plotspict.biomass(seasonres2)
plotspict.f(seasonres2, qlegend=FALSE)

## ---- results='show', message=FALSE, warning=FALSE, fig.width=5.5, fig.height=7, fig.show='hold'----
inp2 <- list(obsC=seasonsim2$obsC, obsI=seasonsim2$obsI,
             timeC=seasonsim2$timeC, timeI=seasonsim2$timeI,
             seasontype=1, true=seasonsim2$true)
rep2 <- fit.spict(inp2)
rep2 <- calc.osa.resid(rep2)
plotspict.diagnostic(rep2)

## ---- results='show', message=FALSE, warning=FALSE, fig.width=4, fig.height=4, fig.show='hold'----
inp <- pol$albacore
inp$ini$logK <- log(100)

## ---- results='show', message=FALSE, warning=FALSE, fig.width=4, fig.height=4, fig.show='hold'----
inp <- check.inp(pol$albacore)
inp$ini$logK

## ---- results='show', message=FALSE, warning=FALSE, fig.width=4, fig.height=4, fig.show='hold'----
set.seed(123)
check.ini(pol$albacore, ntrials=4)

## ---- results='show', message=FALSE, warning=FALSE, fig.width=4, fig.height=4, fig.show='hold'----
inp <- pol$albacore
inp$phases$logsdb <- 2
res <- fit.spict(inp)

## ---- results='show', message=FALSE, warning=FALSE, fig.width=4, fig.height=4, fig.show='hold'----
inp <- pol$albacore
inp$phases$logsdb <- -1
inp$ini$logsdb <- log(0.1)
res <- fit.spict(inp)
summary(res)

## ---- results='show', message=FALSE, warning=FALSE, fig.width=4, fig.height=4, fig.show='hold'----
inp <- pol$albacore
inp$priors$logn <- c(1, 1, 0)
inp$priors$logalpha <- c(1, 1, 0)
inp$priors$logbeta <- c(1, 1, 0)
fit.spict(inp)

## ---- results='show', message=FALSE, warning=FALSE, fig.width=4, fig.height=4, fig.show='hold'----
list.possible.priors()

## ---- results='show', message=FALSE, warning=FALSE, fig.width=4, fig.height=4, fig.show='hold'----
inp <- pol$albacore
inp$priors$logK <- c(log(300), 2, 1)
fit.spict(inp)

## ---- results='show', message=FALSE, warning=FALSE, fig.width=7, fig.height=3.5, fig.show='hold'----
inp <- pol$albacore
inp$priors$logB <- c(log(80), 0.1, 1, 1980)
par(mfrow=c(1, 2), mar=c(5, 4.1, 3, 4))
plotspict.biomass(fit.spict(pol$albacore), ylim=c(0, 500))
plotspict.biomass(fit.spict(inp), qlegend=FALSE, ylim=c(0, 500))

## ---- results='show', message=FALSE, warning=FALSE, fig.width=6, fig.height=5, fig.show='hold'----
set.seed(123)
inp <- list(timeC=pol$albacore$timeC, obsC=pol$albacore$obsC)
inp$timeI <- list(pol$albacore$timeI, pol$albacore$timeI[10:23]+0.25)
inp$obsI <- list()
inp$obsI[[1]] <- pol$albacore$obsI * exp(rnorm(23, sd=0.1)) # Index 1
inp$obsI[[2]] <- 10*pol$albacore$obsI[10:23] # Index 2
inp$priors$logsdi <- list(c(1, 1, 0),          # No prior for index 1
                          c(log(0.1), 0.2, 1)) # Set a prior on index 2
res <- fit.spict(inp)
sumspict.priors(res)

## ---- results='show', message=FALSE, warning=FALSE, fig.width=7, fig.height=3.5, fig.show='hold'----
inp <- pol$albacore
inp$priors$logn <- c(log(2), 1e-3)
inp$priors$logalpha <- c(log(1), 1e-3)
inp$priors$logbeta <- c(log(1), 1e-3)
fit.spict(inp)

## ---- results='show', message=FALSE, warning=FALSE, fig.width=7, fig.height=3.5, fig.show='hold'----
inp <- pol$albacore
inp$obsC[10] <- 3*inp$obsC[10]
res1 <- fit.spict(inp)
inp$robflagc <- 1
res2 <- fit.spict(inp)
sumspict.parest(res2)
par(mfrow=c(1, 2))
plotspict.catch(res1, main='Regular fit')
plotspict.catch(res2, qlegend=FALSE, main='Robust fit')

## ---- results='show', message=FALSE, warning=FALSE, fig.width=5.5, fig.height=5, fig.show='hold'----
inp <- pol$albacore
inp$maninterval <- c(1990,1991)
inp$maneval <- 1991
inp$ffac <- 0.75
rep <- fit.spict(inp)

## ---- results='show', message=FALSE, warning=FALSE, fig.width=5.5, fig.height=5, fig.show='hold'----
sumspict.predictions(rep)

## ---- results='show', message=FALSE, warning=FALSE, fig.width=8, fig.height=7, fig.show='hold'----
plot2(rep)

## ---- results='show', message=FALSE, warning=FALSE, fig.width=5.5, fig.height=5, fig.show='hold'----
inp <- pol$albacore
rep <- fit.spict(inp)
rep <- manage(rep)

## ---- results='show', message=FALSE, warning=FALSE, fig.width=5.5, fig.height=5, fig.show='hold'----
sumspict.manage(rep)

## ---- results='show', message=FALSE, warning=FALSE, fig.width=8, fig.height=7, fig.show='hold'----
plot2(rep)

## ---- results='show', message=FALSE, warning=FALSE, fig.width=5.5, fig.height=5, fig.show='hold'----
man.timeline(inp)

## ---- results='show', message=FALSE, warning=FALSE, fig.width=5.5, fig.height=5, fig.show='hold'----
man.timeline(rep)

## ---- results='show', message=FALSE, warning=FALSE, fig.width=5.5, fig.height=5, fig.show='hold'----
inp$maninterval <- c(1991,1992)
man.timeline(inp)

## ---- results='show', message=FALSE, warning=FALSE, fig.width=5.5, fig.height=5, fig.show='hold'----
repIntPer <- manage(rep, scenarios = c(1,2),
                    maninterval = c(1991,1992), maneval = 1992)
plotspict.catch(repIntPer)

## ---- results='show', message=FALSE, warning=FALSE, fig.width=9, fig.height=5, fig.show='hold'----
repIntPerC <- manage(rep, scenarios=c(8), maninterval = c(1991,1992), intermediatePeriodCatch = 20)
par(mfrow=c(1,2))
plotspict.biomass(repIntPerC)
plotspict.catch(repIntPerC)

## ---- results='show', message=FALSE, warning=FALSE, fig.width=5.5, fig.height=5, fig.show='hold'----
repIntPer <- add.man.scenario(repIntPer, ffac = 1.5)
sumspict.manage(repIntPer)
plotspict.f(repIntPer)

## ---- results='show', message=FALSE, warning=FALSE, fig.width=5.5, fig.height=5, fig.show='hold'----
repIntPer <- add.man.scenario(repIntPer, cfac = 0.64,
                              scenarioTitle = "reduced_catch")
names(repIntPer$man)

## ---- results='show', message=FALSE, warning=FALSE, fig.width=5.5, fig.height=5, fig.show='hold'----
repIntPer <- add.man.scenario(repIntPer, scenarioTitle = "ices",
                         breakpointB = 0.5,
                         fractiles = list(catch=0.35, bbmsy=0.35, ffmsy=0.35))
sumspict.manage(repIntPer)

## ---- results='show', message=FALSE, warning=FALSE, fig.width=5.5, fig.height=5, fig.show='hold'----
repVarPer <- manage(rep, scenarios = c(1,4), maninterval = c(1990.5,1991))

## ---- results='show', message=FALSE, warning=FALSE, fig.width=5.5, fig.height=5, fig.show='hold'----
repVarPer <- add.man.scenario(repVarPer, maninterval = c(1991, 1993),
                              breakpointB = 0.5)
plotspict.f(repVarPer)

## ---- results='show', message=FALSE, warning=FALSE, fig.width=5.5, fig.height=5, fig.show='hold'----
sumspict.manage(repVarPer)

## ---- results='show', message=FALSE, warning=FALSE, fig.width=5.5, fig.height=5, fig.show='hold'----
names(rep$man)
rep <- man.select(rep, scenarios = c(1,3))

## ---- results='show', message=FALSE, warning=FALSE, fig.width=5.5, fig.height=5, fig.show='hold'----
man.tac(rep)

## ---- results='show', message=FALSE, warning=FALSE, fig.width=5.5, fig.height=5, fig.show='hold'----
get.TAC(rep, maninterval = c(1990.25, 1991),
        safeguardB = list(limitB = 0.3), verbose = TRUE)

## ---- results='show', message=FALSE, warning=FALSE, fig.width=8, fig.height=7, fig.show='hold'----
## with index
rep_scenario2 <- rep$man[[2]]
## with scenario name
names(rep$man)
rep_cc <- rep$man$currentCatch
## plotted as standard spict object
plot2(rep_cc)

## ---- results='show', message=FALSE, warning=FALSE, fig.width=5.5, fig.height=5, fig.show='hold'----
rep <- man.select(rep, scenarios = NULL)

