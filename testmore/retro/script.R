##################################################
### Check retrospective analysis functionality ###
##################################################

suppressMessages(library(spict))
options(device = png)

source("../funcs.R")

header("Check retrospective analysis functionality", append = FALSE)

header("1. Retro and Mohn's rho exclude not converged runs")
inp <- pol$albacore
inp$dteuler <- 1/4
inp <- shorten.inp(pol$albacore, maxtime = 1980)

fit <- fit.spict(inp)
fit <- retro(fit, nretroyear = 7)

test_this("1.1. Convergence of last two retros is not achieved", sapply(fit$retro, function(x) x$opt$convergence))

test_this("1.2. Mohn's rho is calculated correctly", capture.output(mr1 <- mohns_rho(fit), type = "message"))

test_this("1.3. Plot is shown correctly", capture.output(mr2 <- plotspict.retro(fit), type = "message"))

test_this("1.4. The results from plotspict.retro and monh's rho should be the same", all.equal(mr1, mr2))

## Expect error - not a fitted object
header("2. Error expectations")
test_this("2.1. No spictcls object produces an error", retro(pol$albacore))

## Expect error - not converged run
set.seed(1234)
inp <- list(obsC = runif(10, 100, 1000),
            timeC = 1:10,
            obsI = runif(10, 0.5, 1.5),
            timeI = 1:10)
out("\n")

suppressWarnings(fitnc <- fit.spict(inp))
test_this("2.2. Not converged fitted object produces an error", retro(fitnc))


