##################################################
###   Check process residuals functionality    ###
##################################################

suppressMessages(library(spict))
options(device = png)

source("../funcs.R")

header("Check process residuals functionality", append = FALSE)


header("1. Process residuals with two productivity regimes")
## ---------------------------------------------------------

## Simulate data with two regimes
set.seed(12345)
nt <- 70
inp <- list(nseasons = 1)
inp$timeC <- seq(0, nt, by = 1/inp$nseasons)
inp$timeI <- list(seq(0.78, nt, by = 1),
                  seq(0.11, nt, by = 1))
inp$ini <- list(logK = log(2000), logm = c(log(900),log(200)), logq = c(log(0.4),log(0.2)),
                logn = log(2), logqf = 0.47,
                logsdb = log(0.1), logsdc = log(0.02),
                logbkfrac = log(0.8), logsdf = log(0.3),
                logF0 = log(0.6), logsdi = c(log(0.05),log(0.1)))
inp$MSYregime <- factor(c(rep(1,300), rep(2,277)))
inp$dteuler <- 1/8
inpsim <- check.inp(sim.spict(inp))

inpsim$priors$logn <- c(0,0,0)
inpsim$priors$logalpha <- c(0,0,0)
inpsim$priors$logbeta <- c(0,0,0)

## Fit with correct regimes
inpsim$MSYregime <- factor(c(rep(1,300), rep(2,277)))
fitreg1 <- fit.spict(inpsim)
fitreg1 <- calc.osa.resid(fitreg1)
fitreg1 <- calc.process.resid(fitreg1)

## Fit without regimes
inpsim$MSYregime <- factor(c(rep(1,300), rep(1,277)))
fitreg2 <- fit.spict(inpsim)
fitreg2 <- calc.osa.resid(fitreg2)
fitreg2 <- calc.process.resid(fitreg2)


## Estimated process residuals
test_this("1.1. Process residuals with regimes: ", signif(fitreg1$process.resid,2))
test_this("1.2. Process residuals without regimes: ", signif(fitreg2$process.resid,2))


## Tests
tmp <- sumspict.diagnostics(fitreg1)
ind <- sapply(1:ncol(tmp), function(x) is.numeric(tmp[,x]))
test_this("1.3. Process residuals tests with regimes okay: ", cbind(signif(tmp[,ind],2),tmp[,!ind]))
tmp <- sumspict.diagnostics(fitreg2)
ind <- sapply(1:ncol(tmp), function(x) is.numeric(tmp[,x]))
test_this("1.4. Process residuals tests without regimes violated: ", cbind(signif(tmp[,ind],2),tmp[,!ind]))


## Plots
out("1.5: Plots")
plotspict.diagnostic(fitreg1)
plotspict.diagnostic.process(fitreg1)

plotspict.diagnostic(fitreg2)
plotspict.diagnostic.process(fitreg2)
