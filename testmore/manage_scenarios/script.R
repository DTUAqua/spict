## Testing new management functionality
## T.K. Mildenberger <t.k.mildenberger@gmail.com>
## 11/12/2019

## set seed
set.seed(123)

## load spict
require(spict)

## load testmore functions
source("../funcs.R")


header("1: Testing add.man.scenarios with annual data", append = FALSE)
######################################################

header("1.1: No intermediate year")
#################################

inp <- pol$lobster
inp$maninterval <- c(1991, 1992)
inp$maneval <- 1992
## increase speed
inp$dteuler <- 1/4
fit <- fit.spict(inp)
man.timeline(fit)


test_this("1.1.1: Fishing at Fmsy", {
    round(get.TAC(fit),3)
})
test_this("1.1.2:", {
    round(get.TAC(fit, fractiles = list(catch=0.2, ffmsy = 0.1)),3)
})

test_this("1.1.3: Fishing at Fmsy with PA buffer", {
    round(get.TAC(fit, safeguardB = list(limitB=0.3, prob=0.8)),3)
})
test_this("1.1.4:", {
    round(get.TAC(fit, safeguardB = list(limitB=0.5, prob=0.9),
                                         fractiles = list(catch=0.2, ffmsy = 0.1)),3)
})

test_this("1.1.5: MSY hockey-stick rule", {
    round(get.TAC(fit, breakpointB = 0.5),3)
})
test_this("1.1.6:", {
    round(get.TAC(fit, breakpointB = 0.5, fractiles = list(catch=0.2, ffmsy=0.1, bbmsy=0.1)),3)
})

test_this("1.1.7: MSY hockey-stick rule with PA buffer", {
    round(get.TAC(fit, breakpointB = 0.5, safeguardB = list(limitB=0.3, prob=0.9)),3)
})
test_this("1.1.8:", {
    round(get.TAC(fit, breakpointB = 0.3, safeguardB = list(limitB=0.5, prob=0.9),
                                         fractiles = list(catch=0.2, ffmsy = 0.1)),3)
})


header("1.2: intermediate year with constant F")
#################################

inp$maninterval <- c(1992, 1993)
inp$maneval <- 1993
fit <- fit.spict(inp)

test_this("1.2.1: Fishing at Fmsy", {
    round(get.TAC(fit),3)
})
test_this("1.2.2:", {
    round(get.TAC(fit, fractiles = list(catch=0.2)),3)
})

test_this("1.2.3: Fishing at Fmsy with biomass safeguard", {
    round(get.TAC(fit, safeguardB = list(limitB=0.3, prob=0.9)),3)
})
test_this("1.2.4:", {
    round(get.TAC(fit, safeguardB = list(limitB=0.5, prob=0.6)),3)
})

test_this("1.2.5: MSY hockey-stick rule", {
    round(get.TAC(fit, breakpointB = 0.3),3)
})
test_this("1.2.6:", {
    round(get.TAC(fit, breakpointB = 0.5,
                                        fractiles = list(catch=0.2, ffmsy=0.1, limitB=0.1)),3)
})

test_this("1.2.7: MSY hockey-stick rule with biomass safeguard", {
    round(get.TAC(fit, breakpointB = 0.3, safeguardB = list(limitB=0.3, prob=0.8)),3)
})
test_this("1.2.8:", {
    round(get.TAC(fit, breakpointB = 0.3, safeguardB = list(limitB=0.5, prob=0.9),
                                        fractiles = list(catch=0.2, ffmsy = 0.1)),3)
})


header("1.3: intermediate year with constant catch")
#################################

inp$maninterval <- c(1992, 1993)
inp$maneval <- 1993
fit <- fit.spict(inp)
lastC <- tail(inp$obsC,1)

test_this("1.3.1: Fishing at Fmsy", {
    round(get.TAC(fit, intermediatePeriodCatch = lastC),3)
})
test_this("1.3.2:", {
    round(get.TAC(fit, fractiles = list(catch=0.2), intermediatePeriodCatch = lastC),3)
})

test_this("1.3.3: Fishing at Fmsy with biomass safeguard", {
    round(get.TAC(fit, safeguardB = list(limitB=0.3, prob=0.9), intermediatePeriodCatch = lastC),3)
})
test_this("1.3.4:", {
    round(get.TAC(fit, safeguardB = list(limitB=0.5, prob=0.6), intermediatePeriodCatch = lastC),3)
})

test_this("1.3.5: MSY hockey-stick rule", {
    round(get.TAC(fit, breakpointB = 0.3, intermediatePeriodCatch = lastC),3)
})
test_this("1.3.6:", {
    round(get.TAC(fit, breakpointB = 0.5,
                                        fractiles = list(catch=0.2, ffmsy=0.1, limitB=0.1),
                                        intermediatePeriodCatch = lastC),3)
})

test_this("1.3.7: MSY hockey-stick rule with biomass safeguard", {
    round(get.TAC(fit, breakpointB = 0.3, safeguardB = list(limitB=0.3, prob=0.8),
                                        intermediatePeriodCatch = lastC),3)
})
test_this("1.3.8:", {
    round(get.TAC(fit, breakpointB = 0.3, safeguardB = list(limitB=0.5, prob=0.9),
                                        fractiles = list(catch=0.2, ffmsy = 0.1),
                                        intermediatePeriodCatch = lastC),3)
})


header("2: Testing management scenarios with seasonal data")
#####################################################################

nt <- 50
inp <- list(nseasons = 4, splineorder = 3)
inp$timeC <- seq(0, nt - 1 / inp$nseasons, by = 1 / inp$nseasons)
inp$timeI <- seq(0.1, nt - 1 / inp$nseasons, by = 0.5)
inp$ini <- list(logK = log(1000), logm=log(800), logq = log(1), logn = log(2),
                logbkfrac = log(0.9), logsdf = log(0.3), logF0 = log(0.8),
                logphi = log(c(0.3, 0.5, 1.8)))
inp$dteuler <- 1/4
inpsim <- sim.spict(inp)


header("2.1: standard advice")
###############################

inpsim$maneval <- nt+1
inpsim$maninterval <- c(nt,nt+1)
inp <- check.inp(inpsim)
fit <- fit.spict(inp)

test_this("2.1.1: Fishing at Fmsy", {
    round(get.TAC(fit),3)
})
test_this("2.1.2:", {
    round(get.TAC(fit, fractiles = list(catch=0.2)),3)
})
test_this("2.1.3:", {
    round(get.TAC(fit, fractiles = list(catch=0.2, ffmsy = 0.1)),3)
})

test_this("2.1.4: Fishing at Fmsy with PA buffer", {
    round(get.TAC(fit, safeguardB = list(limitB=0.3, prob=0.8)),3)
})
test_this("2.1.5:", {
    round(get.TAC(fit, safeguardB = list(limitB=0.5, prob=0.65)),3)
})
test_this("2.1.6:", {
    round(get.TAC(fit, safeguardB = list(limitB=0.5, prob=0.9),
                  fractiles = list(catch=0.2,ffmsy = 0.1)),3)
})

test_this("2.1.7: MSY hockey-stick rule", {
    round(get.TAC(fit, breakpointB = 0.3),3)
})
test_this("2.1.8:", {
    round(get.TAC(fit, breakpointB = 0.5, fractiles = list(catch=0.2, ffmsy=0.1, bbmsy=0.1)),3)
})

test_this("2.1.9: MSY hockey-stick rule with PA buffer", {
    round(get.TAC(fit, breakpointB = 0.3, safeguardB = list(limitB=0.3, prob=0.8)),3)
})
test_this("2.1.10:", {
    round(get.TAC(fit, breakpointB = 0.3, safeguardB = list(limitB=0.5, prob=0.9),
                  fractiles = list(catch=0.2, ffmsy = 0.1)),3)
})


header("2.2: in year advice")
###############################

inpsim$maneval <- nt+1.5
inpsim$maninterval <- c(nt+0.5,nt+1.5)
inp <- check.inp(inpsim)
fit <- fit.spict(inp)

test_this("2.2.1: Fishing at Fmsy", {
    round(get.TAC(fit),3)
})
test_this("2.2.2:", {
    round(get.TAC(fit, fractiles = list(catch=0.2, ffmsy = 0.1)),3)
})

test_this("2.2.3: Fishing at Fmsy with PA buffer", {
    round(get.TAC(fit, safeguardB = list(limitB=0.3, prob=0.8)),3)
})
test_this("2.2.4:", {
    round(get.TAC(fit, safeguardB = list(limitB=0.5, prob=0.65)),3)
})

test_this("2.2.5: MSY hockey-stick rule", {
    round(get.TAC(fit, breakpointB = 0.5, fractiles = list(catch=0.2, ffmsy=0.1)),3)
})


test_this("2.2.6: MSY hockey-stick rule with PA buffer", {
    round(get.TAC(fit, breakpointB = 0.3, safeguardB = list(limitB=0.5, prob=0.9),
                  fractiles = list(catch=0.2, ffmsy = 0.1)),3)
})

test_this("2.2.7: Fishing at Fmsy with TAC during assessment year", {
    round(get.TAC(fit, intermediatePeriodCatch = 100,
                  fractiles = list(catch=0.2, ffmsy = 0.1)),3)
})





header("3: Testing management scenarios with mixed data")
#####################################################################

set.seed(123)
nt <- 50
inp <- list(nseasons=4, splineorder=3)
inp$timeC <- c(seq(0, nt/2-1, by=1),seq(nt/2, nt-1/inp$nseasons, by=1/inp$nseasons))
inp$timeI <- seq(0.1, nt-1/inp$nseasons, by=0.5)
inp$ini <- list(logK=log(1000), logm=log(500), logq=log(1), logn=log(2),
                logbkfrac=log(0.9), logsdf=log(0.3), logF0=log(1),
                logphi=log(c(0.3, 0.5, 1.8)))
inp$dteuler <- 1/4
inpsim <- sim.spict(inp)




header("3.1: standard advice")
###############################

inpsim$maneval <- nt+1
inpsim$maninterval <- c(nt,nt+1)
inp <- check.inp(inpsim)
fit <- fit.spict(inp)


test_this("3.1.1: Fishing at Fmsy", {
    round(get.TAC(fit),3)
})
test_this("3.1.2: Using catch and F/Fmsy fractiles", {
    round(get.TAC(fit, fractiles = list(catch=0.2, ffmsy = 0.1)),3)
})

test_this("3.1.3: Fishing at Fmsy with PA buffer", {
    round(get.TAC(fit, safeguardB = list(limitB=0.5, prob=0.65)),3)
})

test_this("3.1.4: MSY hockey-stick rule", {
    round(get.TAC(fit, breakpointB = 0.5, fractiles = list(catch=0.2, ffmsy=0.1)),3)
})

test_this("3.1.5: MSY hockey-stick rule with PA buffer", {
    round(get.TAC(fit, breakpointB = 0.3, safeguardB = list(limitB=0.5, prob=0.9),
                  fractiles = list(catch=0.2, ffmsy = 0.1)),3)
})

test_this("3.1.6: Fishing at Fmsy with TAC during assessment year", {
    round(get.TAC(fit, intermediatePeriodCatch = 40,
                  fractiles = list(catch=0.2, ffmsy = 0.1)),3)
})

header("3.2: in year advice")
###############################
inpsim$maneval <- nt+1.5
inpsim$maninterval <- c(nt+0.5,nt+1.5)
inp <- check.inp(inpsim)
fit <- fit.spict(inp)

test_this("3.2.1: Fishing at Fmsy", {
    round(get.TAC(fit),3)
})
test_this("3.2.2: Using catch and F/Fmsy fractiles", {
    round(get.TAC(fit, fractiles = list(catch=0.2, ffmsy = 0.1)),3)
})

test_this("3.2.3: Fishing at Fmsy with PA buffer", {
    round(get.TAC(fit, safeguardB = list(limitB=0.3, prob=0.9)),3)
})

test_this("3.2.4: MSY hockey-stick rule", {
    round(get.TAC(fit, breakpointB = 0.5, fractiles = list(catch=0.2, ffmsy=0.1)),3)
})


test_this("3.2.5: MSY hockey-stick rule with PA buffer", {
    round(get.TAC(fit, breakpointB = 0.3, safeguardB = list(limitB=0.5, prob=0.9),
                  fractiles = list(catch=0.2, ffmsy = 0.1)),3)
})

test_this("3.2.6: Fishing at Fmsy with TAC during assessment year", {
    round(get.TAC(fit, intermediatePeriodCatch = 40,
                  fractiles = list(catch=0.2, ffmsy = 0.1)),3)
})



header("4: Other special cases")
#####################################################################

set.seed(123)
nt <- 50
inp <- list(nseasons = 1, splineorder = 1)
inp$timeC <- seq(0, nt - 1 / inp$nseasons, by = 1 / inp$nseasons)
inp$timeI <- seq(0.375, nt - 1 / inp$nseasons, by = 1)
inp$ini <- list(logK = log(1000), logm = log(200), logq = log(1), logn = log(2),
                logbkfrac = log(0.9), logsdf = log(0.3), logF0 = log(0.4))
inp$dteuler <- 1/4
inpsim <- sim.spict(inp)
## plotspict.data(inpsim)


header("4.1: Last index observation inside intermediate year")
###############################

inp <- list()
inp$obsC <- inpsim$obsC[-((nt-1):nt)]
inp$timeC <- inpsim$timeC[-((nt-1):nt)]
inp$obsI <- inpsim$obsI
inp$timeI <- inpsim$timeI
inp$maneval <- nt
inp$maninterval <- c(nt-1,nt)
inp <- check.inp(inp)

man.timeline(inp)
tail(inp$timeC + inp$dtc,1)
tail(inp$timeI[[1]],1)

## Fit spict
fit <- fit.spict(inp)


## Management scenarios
fit <- add.man.scenario(fit)
fit <- add.man.scenario(fit, ffac = 0)
fit <- add.man.scenario(fit, ffac = 1)


test_this("4.1.1: Summary without intermediatePeriodCatch", {
    sumspict.manage(fit)
})

plotspict.catch(fit)



## Management scenarios with specified intermediate period catch
fit <- add.man.scenario(fit, intermediatePeriodCatch = mean(inp$obsC))
fit <- add.man.scenario(fit, ffac = 0, intermediatePeriodCatch = mean(inp$obsC))
fit <- add.man.scenario(fit, ffac = 1, intermediatePeriodCatch = mean(inp$obsC))

test_this("4.1.4: Summary with intermediatePeriodCatch", {
    sumspict.manage(fit)
})

plotspict.catch(fit)
