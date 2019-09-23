## Testing some new functions and modifications
## T.K. Mildenberger <t.k.mildenberger@gmail.com>
## 23/09/2019
## ----------------------------------------
library(spict)


## set seed
## ----------------------------------------
set.seed(1234)


## simulate data
## ----------------------------------------
inp <- list(nseasons=4, splineorder=3)
inp$timeC <- seq(0, 30-1/inp$nseasons, by=1/inp$nseasons)
inp$timeI <- seq(0, 30-1/inp$nseasons, by=1/inp$nseasons)
inp$ini <- list(logK=log(100), logm=log(20), logq=log(1),
                logbkfrac=log(1), logsdf=log(0.4), logF0=log(0.5),
                logphi=log(c(0.05, 0.1, 1.8)))
inp <- sim.spict(inp)


## new functions
## ----------------------------------------
## shorten inp
## ----------------------------------------
inp <- shorten.inp(inp, 10, 25)
capture.output(inp, file="res.out")

## remove priors
## ----------------------------------------
inp$priors$logbkfrac <- c(log(0.2),0.2,1)
inp <- remove.priors(inp, priors = c("logn","logbkfrac"))
capture.output(inp$priors, file="res.out", append = TRUE)

## change and check euler
## ----------------------------------------
inp <- pol$albacore
inp <- shorten.inp(inp, 1975)
inpalt <- change.euler(inp, dteuler = 1/4)
repalt <- fit.spict(inpalt)
eulercheck <- check.euler(repalt, dteuler=1/16)
capture.output(eulercheck, file="res.out", append = TRUE)

## Bmsy/K
## ----------------------------------------
bmsyk <- calc.bmsyk(repalt)
capture.output(bmsyk, file="res.out", append = TRUE)

## order of magnitude of CI range
## ----------------------------------------
om <- calc.om(repalt)
capture.output(om, file="res.out", append = TRUE)


## bug fixes
## ----------------------------------------
## inp summary with effort data
inp <- list()
inp$obsC <- seq(100,500,length.out = 10)
inp$timeC <- 1991:2000
inp$obsE <- seq(500,100,length.out = 10)
inp$timeE <- 1991:2000
inp <- check.inp(inp)
capture.output(inp, file="res.out", append = TRUE)
