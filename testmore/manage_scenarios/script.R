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

out("1.1: No intermediate year")
#################################
inp <- pol$lobster
inp$maninterval <- c(1991, 1992)
inp$maneval <- 1992
## increase speed
inp$dteuler <- 1/4
fit <- fit.spict(inp)
man.timeline(fit)


out("Fishing at Fmsy")
test_this("1.1.1:", {
    round(get.TAC(fit),3)
})
test_this("1.1.2:", {
    round(get.TAC(fit, fractiles = list(catch=0.2, ffmsy = 0.1)),3)
})

out("Fishing at Fmsy with PA buffer")
test_this("1.1.3:", {
    round(get.TAC(fit, safeguardB = list(limitB=0.3, prob=0.8)),3)
})
test_this("1.1.4:", {
    round(get.TAC(fit, safeguardB = list(limitB=0.5, prob=0.9),
                                         fractiles = list(catch=0.2, ffmsy = 0.1)),3)
})

out("MSY hockey-stick rule")
test_this("1.1.5:", {
    round(get.TAC(fit, breakpointB = 0.5),3)
})
test_this("1.1.6:", {
    round(get.TAC(fit, breakpointB = 0.5, fractiles = list(catch=0.2, ffmsy=0.1, bbmsy=0.1)),3)
})

out("MSY hockey-stick rule with PA buffer")
test_this("1.1.7:", {
    round(get.TAC(fit, breakpointB = 0.5, safeguardB = list(limitB=0.3, prob=0.9)),3)
})
test_this("1.1.8:", {
    round(get.TAC(fit, breakpointB = 0.3, safeguardB = list(limitB=0.5, prob=0.9),
                                         fractiles = list(catch=0.2, ffmsy = 0.1)),3)
})


out("1.2: intermediate year with constant F")
#################################
inp$maninterval <- c(1992, 1993)
inp$maneval <- 1993
fit <- fit.spict(inp)

out("Fishing at Fmsy")
test_this("1.2.1:", {
    round(get.TAC(fit),3)
})
test_this("1.2.2:", {
    round(get.TAC(fit, fractiles = list(catch=0.2)),3)
})

out("Fishing at Fmsy with biomass safeguard")
test_this("1.2.3:", {
    round(get.TAC(fit, safeguardB = list(limitB=0.3, prob=0.9)),3)
})
test_this("1.2.4:", {
    round(get.TAC(fit, safeguardB = list(limitB=0.5, prob=0.6)),3)
})

out("MSY hockey-stick rule")
test_this("1.2.5:", {
    round(get.TAC(fit, breakpointB = 0.3),3)
})
test_this("1.2.6:", {
    round(get.TAC(fit, breakpointB = 0.5,
                                        fractiles = list(catch=0.2, ffmsy=0.1, limitB=0.1)),3)
})

out("MSY hockey-stick rule with biomass safeguard")
test_this("1.2.7:", {
    round(get.TAC(fit, breakpointB = 0.3, safeguardB = list(limitB=0.3, prob=0.8)),3)
})
test_this("1.2.8:", {
    round(get.TAC(fit, breakpointB = 0.3, safeguardB = list(limitB=0.5, prob=0.9),
                                        fractiles = list(catch=0.2, ffmsy = 0.1)),3)
})


out("1.3: intermediate year with constant catch")
#################################
inp$maninterval <- c(1992, 1993)
inp$maneval <- 1993
fit <- fit.spict(inp)
lastC <- tail(inp$obsC,1)

out("Fishing at Fmsy")
test_this("1.3.1:", {
    round(get.TAC(fit, intermediatePeriodCatch = lastC),3)
})
test_this("1.3.2:", {
    round(get.TAC(fit, fractiles = list(catch=0.2), catchintermediateYear = lastC),3)
})

out("Fishing at Fmsy with biomass safeguard")
test_this("1.3.3:", {
    round(get.TAC(fit, safeguardB = list(limitB=0.3, prob=0.9), intermediatePeriodCatch = lastC),3)
})
test_this("1.3.4:", {
    round(get.TAC(fit, safeguardB = list(limitB=0.5, prob=0.6), intermediatePeriodCatch = lastC),3)
})

out("MSY hockey-stick rule")
test_this("1.3.5:", {
    round(get.TAC(fit, breakpointB = 0.3, intermediatePeriodCatch = lastC),3)
})
test_this("1.3.6:", {
    round(get.TAC(fit, breakpointB = 0.5,
                                        fractiles = list(catch=0.2, ffmsy=0.1, limitB=0.1),
                                        intermediatePeriodCatch = lastC),3)
})

out("MSY hockey-stick rule with biomass safeguard")
test_this("1.3.7:", {
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


out("2.1: standard advice")
###############################
inpsim$maneval <- nt+1
inpsim$maninterval <- c(nt,nt+1)
inp <- check.inp(inpsim)
fit <- fit.spict(inp)

out("Fishing at Fmsy")
test_this("2.1.1:", {
    round(get.TAC(fit),3)
})
test_this("2.1.2:", {
    round(get.TAC(fit, fractiles = list(catch=0.2)),3)
})
test_this("2.1.3:", {
    round(get.TAC(fit, fractiles = list(catch=0.2, ffmsy = 0.1)),3)
})

out("Fishing at Fmsy with PA buffer")
test_this("2.1.4:", {
    round(get.TAC(fit, safeguardB = list(limitB=0.3, prob=0.8)),3)
})
test_this("2.1.5:", {
    round(get.TAC(fit, safeguardB = list(limitB=0.5, prob=0.65)),3)
})
test_this("2.1.6:", {
    round(get.TAC(fit, safeguardB = list(limitB=0.5, prob=0.9),
                  fractiles = list(catch=0.2,ffmsy = 0.1)),3)
})

out("MSY hockey-stick rule")
test_this("2.1.7:", {
    round(get.TAC(fit, breakpointB = 0.3),3)
})
test_this("2.1.8:", {
    round(get.TAC(fit, breakpointB = 0.5, fractiles = list(catch=0.2, ffmsy=0.1, bbmsy=0.1)),3)
})

out("MSY hockey-stick rule with PA buffer")
test_this("2.1.9:", {
    round(get.TAC(fit, breakpointB = 0.3, safeguardB = list(limitB=0.3, prob=0.8)),3)
})
test_this("2.1.10:", {
    round(get.TAC(fit, breakpointB = 0.3, safeguardB = list(limitB=0.5, prob=0.9),
                  fractiles = list(catch=0.2, ffmsy = 0.1)),3)
})


out("2.2: in year advice")
###############################
inpsim$maneval <- nt+1.5
inpsim$maninterval <- c(nt+0.5,nt+1.5)
inp <- check.inp(inpsim)
fit <- fit.spict(inp)

out("Fishing at Fmsy")
test_this("2.2.1:", {
    round(get.TAC(fit),3)
})
test_this("2.2.2:", {
    round(get.TAC(fit, fractiles = list(catch=0.2, ffmsy = 0.1)),3)
})

out("Fishing at Fmsy with PA buffer")
test_this("2.2.3:", {
    round(get.TAC(fit, safeguardB = list(limitB=0.3, prob=0.8)),3)
})
test_this("2.2.4:", {
    round(get.TAC(fit, safeguardB = list(limitB=0.5, prob=0.65)),3)
})

out("MSY hockey-stick rule")
test_this("2.2.5:", {
    round(get.TAC(fit, breakpointB = 0.5, fractiles = list(catch=0.2, ffmsy=0.1)),3)
})


out("MSY hockey-stick rule with PA buffer")
test_this("2.2.6:", {
    round(get.TAC(fit, breakpointB = 0.3, safeguardB = list(limitB=0.5, prob=0.9),
                  fractiles = list(catch=0.2, ffmsy = 0.1)),3)
})

out("Fishing at Fmsy with TAC during assessment year")
test_this("2.2.7:", {
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



out("3.1: standard advice")
###############################
inpsim$maneval <- nt+1
inpsim$maninterval <- c(nt,nt+1)
inp <- check.inp(inpsim)
fit <- fit.spict(inp)



out("Fishing at Fmsy")
test_this("3.1.1:", {
    round(get.TAC(fit),3)
})
test_this("3.1.2:", {
    round(get.TAC(fit, fractiles = list(catch=0.2, ffmsy = 0.1)),3)
})

out("Fishing at Fmsy with PA buffer")
test_this("3.1.3:", {
    round(get.TAC(fit, safeguardB = list(limitB=0.5, prob=0.65)),3)
})

out("MSY hockey-stick rule")
test_this("3.1.4:", {
    round(get.TAC(fit, breakpointB = 0.5, fractiles = list(catch=0.2, ffmsy=0.1)),3)
})


out("MSY hockey-stick rule with PA buffer")
test_this("3.1.5:", {
    round(get.TAC(fit, breakpointB = 0.3, safeguardB = list(limitB=0.5, prob=0.9),
                  fractiles = list(catch=0.2, ffmsy = 0.1)),3)
})

out("Fishing at Fmsy with TAC during assessment year")
test_this("4.1.6:", {
    round(get.TAC(fit, intermediatePeriodCatch = 40,
                  fractiles = list(catch=0.2, ffmsy = 0.1)),3)
})



out("3.2: in year advice")
###############################
inpsim$maneval <- nt+1.5
inpsim$maninterval <- c(nt+0.5,nt+1.5)
inp <- check.inp(inpsim)
fit <- fit.spict(inp)

out("Fishing at Fmsy")
test_this("4.2.1:", {
    round(get.TAC(fit),3)
})
test_this("4.2.2:", {
    round(get.TAC(fit, fractiles = list(catch=0.2, ffmsy = 0.1)),3)
})

out("Fishing at Fmsy with PA buffer")
test_this("4.2.3:", {
    round(get.TAC(fit, safeguardB = list(limitB=0.3, prob=0.9)),3)
})

out("MSY hockey-stick rule")
test_this("4.2.4:", {
    round(get.TAC(fit, breakpointB = 0.5, fractiles = list(catch=0.2, ffmsy=0.1)),3)
})


out("MSY hockey-stick rule with PA buffer")
test_this("4.2.5:", {
    round(get.TAC(fit, breakpointB = 0.3, safeguardB = list(limitB=0.5, prob=0.9),
                  fractiles = list(catch=0.2, ffmsy = 0.1)),3)
})

out("Fishing at Fmsy with TAC during assessment year")
test_this("4.2.6:", {
    round(get.TAC(fit, intermediatePeriodCatch = 40,
                  fractiles = list(catch=0.2, ffmsy = 0.1)),3)
})
