## Testing new function "get.TAC"
## Tobias K. Mildenberger <t.k.mildenberger@gmail.com>
## November 2019

## load spict
suppressMessages(require(spict))

## annual data
## -----------------------------------------------------------------------------
inp <- pol$lobster

##########################
## no intermediate year ##
##########################
inp$maninterval <- c(1991, 1992) ## default anyways
inp$maneval <- 1992  ## crucial for PA buffer because timepredi = 1991 by default
fit <- fit.spict(inp)

## Fishing at Fmsy
round(get.TAC(fit)$TAC,3)
round(get.TAC(fit, fractileList = list(catch=0.2))$TAC,3)
round(get.TAC(fit, fractileList = list(catch=0.2, ffmsy = 0.1))$TAC,3)

## Fishing at Fmsy with PA buffer
round(get.TAC(fit, paList = list(bbmsy=0.3, prob=0.8))$TAC,3)
round(get.TAC(fit, paList = list(bbmsy=0.3, prob=0.9))$TAC,3)
round(get.TAC(fit, paList = list(bbmsy=0.3, prob=0.93))$TAC,3)
round(get.TAC(fit, paList = list(bbmsy=0.3, prob=0.95))$TAC,3)
round(get.TAC(fit, paList = list(bbmsy=0.5, prob=0.6))$TAC,3)
round(get.TAC(fit, paList = list(bbmsy=0.5, prob=0.65))$TAC,3)
round(get.TAC(fit, paList = list(bbmsy=0.5, prob=0.7))$TAC,3)
round(get.TAC(fit, paList = list(bbmsy=0.3, prob=0.9), fractileList = list(catch=0.2))$TAC,3)
round(get.TAC(fit, paList = list(bbmsy=0.5, prob=0.9), fractileList = list(catch=0.2, ffmsy = 0.1))$TAC,3)

## MSY hockey-stick rule
round(get.TAC(fit, breakpoint_bbmsy = 0.3)$TAC,3)
round(get.TAC(fit, breakpoint_bbmsy = 0.5)$TAC,3)
round(get.TAC(fit, breakpoint_bbmsy = 0.5, fractileList = list(catch=0.2))$TAC,3)
round(get.TAC(fit, breakpoint_bbmsy = 0.5, fractileList = list(catch=0.2, ffmsy=0.1))$TAC,3)
round(get.TAC(fit, breakpoint_bbmsy = 0.5, fractileList = list(catch=0.2, ffmsy=0.1, bbmsy=0.1))$TAC,3)

## MSY hockey-stick rule with PA buffer
round(get.TAC(fit, breakpoint_bbmsy = 0.3, paList = list(bbmsy=0.3, prob=0.8))$TAC,3)
round(get.TAC(fit, breakpoint_bbmsy = 0.5, paList = list(bbmsy=0.3, prob=0.9))$TAC,3)
round(get.TAC(fit, breakpoint_bbmsy = 0.5, paList = list(bbmsy=0.3, prob=0.93))$TAC,3)
round(get.TAC(fit, breakpoint_bbmsy = 0.3, paList = list(bbmsy=0.3, prob=0.95))$TAC,3)
round(get.TAC(fit, breakpoint_bbmsy = 0.3, paList = list(bbmsy=0.5, prob=0.6))$TAC,3)
round(get.TAC(fit, breakpoint_bbmsy = 0.3, paList = list(bbmsy=0.5, prob=0.65))$TAC,3)
round(get.TAC(fit, breakpoint_bbmsy = 0.5, paList = list(bbmsy=0.5, prob=0.7))$TAC,3)
round(get.TAC(fit, breakpoint_bbmsy = 0.5, paList = list(bbmsy=0.3, prob=0.9), fractileList = list(catch=0.2))$TAC,3)
round(get.TAC(fit, breakpoint_bbmsy = 0.3, paList = list(bbmsy=0.5, prob=0.9), fractileList = list(catch=0.2, ffmsy = 0.1))$TAC,3)

#######################
## intermediate year ##
#######################
inp$maninterval <- c(1992, 1993) ## default anyways
inp$maneval <- 1993  ## crucial for PA buffer because timepredi = 1991 by default
fit <- fit.spict(inp)

## Fishing at Fmsy
round(get.TAC(fit)$TAC,3)
round(get.TAC(fit, fractileList = list(catch=0.2))$TAC,3)

## Fishing at Fmsy with PA buffer
round(get.TAC(fit, paList = list(bbmsy=0.3, prob=0.9))$TAC,3)
round(get.TAC(fit, paList = list(bbmsy=0.5, prob=0.6))$TAC,3)

## MSY hockey-stick rule
round(get.TAC(fit, breakpoint_bbmsy = 0.3)$TAC,3)
round(get.TAC(fit, breakpoint_bbmsy = 0.5, fractileList = list(catch=0.2, ffmsy=0.1, bbmsy=0.1))$TAC,3)

## MSY hockey-stick rule with PA buffer
round(get.TAC(fit, breakpoint_bbmsy = 0.3, paList = list(bbmsy=0.3, prob=0.8))$TAC,3)
round(get.TAC(fit, breakpoint_bbmsy = 0.3, paList = list(bbmsy=0.5, prob=0.9), fractileList = list(catch=0.2, ffmsy = 0.1))$TAC,3)



## seasonal data
## -----------------------------------------------------------------------------
set.seed(1234)
inp <- list(nseasons=4, splineorder=3)
inp$timeC <- seq(0, 30-1/inp$nseasons, by=1/inp$nseasons)
inp$timeI <- seq(0, 30-1/inp$nseasons, by=1/inp$nseasons)
inp$ini <- list(logK=log(100), logm=log(20), logq=log(1),
                logbkfrac=log(1), logsdf=log(0.4), logF0=log(1),
                logphi=log(c(0.05, 0.1, 1.8)))
inpsim <- sim.spict(inp)

#####################
## standard advice ##
#####################
inpsim$maneval <- 31
inpsim$maninterval <- c(30,31)
inp <- check.inp(inpsim)
fit <- fit.spict(inpsim)

## Fishing at Fmsy
round(get.TAC(fit)$TAC,3)
round(get.TAC(fit, fractileList = list(catch=0.2))$TAC,3)
round(get.TAC(fit, fractileList = list(catch=0.2, ffmsy = 0.1))$TAC,3)

## Fishing at Fmsy with PA buffer
round(get.TAC(fit, paList = list(bbmsy=0.3, prob=0.8))$TAC,3)
round(get.TAC(fit, paList = list(bbmsy=0.5, prob=0.65))$TAC,3)
round(get.TAC(fit, paList = list(bbmsy=0.5, prob=0.9), fractileList = list(catch=0.2, ffmsy = 0.1))$TAC,3)

## MSY hockey-stick rule
round(get.TAC(fit, breakpoint_bbmsy = 0.3)$TAC,3)
round(get.TAC(fit, breakpoint_bbmsy = 0.5, fractileList = list(catch=0.2, ffmsy=0.1, bbmsy=0.1))$TAC,3)

## MSY hockey-stick rule with PA buffer
round(get.TAC(fit, breakpoint_bbmsy = 0.3, paList = list(bbmsy=0.3, prob=0.8))$TAC,3)
round(get.TAC(fit, breakpoint_bbmsy = 0.5, paList = list(bbmsy=0.3, prob=0.9))$TAC,3)
round(get.TAC(fit, breakpoint_bbmsy = 0.5, paList = list(bbmsy=0.3, prob=0.9), fractileList = list(catch=0.2))$TAC,3)
round(get.TAC(fit, breakpoint_bbmsy = 0.3, paList = list(bbmsy=0.5, prob=0.9), fractileList = list(catch=0.2, ffmsy = 0.1))$TAC,3)

####################
## in-year advice ##
####################
inpsim$maneval <- 31.5
inpsim$maninterval <- c(30.5,31.5)
inp <- check.inp(inpsim)
fit <- fit.spict(inpsim)

## Fishing at Fmsy
round(get.TAC(fit)$TAC,3)
round(get.TAC(fit, fractileList = list(catch=0.2))$TAC,3)
round(get.TAC(fit, fractileList = list(catch=0.2, ffmsy = 0.1))$TAC,3)

## Fishing at Fmsy with PA buffer
round(get.TAC(fit, paList = list(bbmsy=0.3, prob=0.8))$TAC,3)
round(get.TAC(fit, paList = list(bbmsy=0.3, prob=0.9))$TAC,3)
round(get.TAC(fit, paList = list(bbmsy=0.5, prob=0.6))$TAC,3)
round(get.TAC(fit, paList = list(bbmsy=0.5, prob=0.65))$TAC,3)

## MSY hockey-stick rule
round(get.TAC(fit, breakpoint_bbmsy = 0.3)$TAC,3)
round(get.TAC(fit, breakpoint_bbmsy = 0.5)$TAC,3)
round(get.TAC(fit, breakpoint_bbmsy = 0.5, fractileList = list(catch=0.2, ffmsy=0.1))$TAC,3)


## MSY hockey-stick rule with PA buffer
round(get.TAC(fit, breakpoint_bbmsy = 0.5, paList = list(bbmsy=0.3, prob=0.9))$TAC,3)
round(get.TAC(fit, breakpoint_bbmsy = 0.5, paList = list(bbmsy=0.3, prob=0.9), fractileList = list(catch=0.2))$TAC,3)
round(get.TAC(fit, breakpoint_bbmsy = 0.3, paList = list(bbmsy=0.5, prob=0.9), fractileList = list(catch=0.2, ffmsy = 0.1))$TAC,3)

## Fishing at Fmsy with TAC during assessment year
round(get.TAC(fit, catch_pred = 8)$TAC,3)
round(get.TAC(fit, catch_pred = 4)$TAC,3)
round(get.TAC(fit, catch_pred = 4, fractileList = list(catch=0.2))$TAC,3)
round(get.TAC(fit, catch_pred = 4, fractileList = list(catch=0.2, ffmsy = 0.1))$TAC,3)
