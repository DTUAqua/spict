##################################################
###      Check hindcasting functionality       ###
##################################################

suppressMessages(library(spict))
options(device = png)

set.seed(12345)

source("../funcs.R")

header("Check hindcasting functionality", append = FALSE)

header("1. Expected basic functionality")
inp <- pol$albacore
inp$dteuler <- 1/4

fit <- fit.spict(inp)
suppressWarnings(fit <- hindcast(fit, nyears = 5))

test_this("1.1. Decreasing length of index observations",
          (nobsi <- sapply(fit$hindcast, function(x) length(x$inp$obsI[[1]])))
          )

test_this("1.2. Decreasing length of index use indicator",
          (nobsiind <- sapply(fit$hindcast, function(x) length(which(x$inp$iuse == TRUE))))
          )

test_this("1.3. Except for first year (ref), difference corresponds to number of indices in this year:",
          (nobsi - nobsiind)
          )

## Last index observation after last catch observation
inp <- pol$albacore
inp$dteuler <- 1/4
inp$timeI <- inp$timeI + 1.8

fit <- fit.spict(inp)
suppressWarnings(fit <- hindcast(fit, nyears = 5))

test_this("1.4. Same pattern when last index observation after last catch observation",
          sapply(fit$hindcast, function(x) length(x$inp$obsI[[1]])) -
          sapply(fit$hindcast, function(x) length(which(x$inp$iuse == TRUE)))
          )


header("2. Multiple indices with different lengths")
nt <- 50
inp <- list(nseasons = 1)
inp$timeC <- seq(0, nt - 1 / inp$nseasons, by = 1)
inp$timeI <- list(seq(0.1, nt, by = 1),
                  seq(0.8, nt-3, by = 1))
inp$ini <- list(logK = log(1000), logm=log(800), logq = log(1), logn = log(2),
                logbkfrac = log(0.9), logsdf = log(0.3), logF0 = log(0.2))
inp$dteuler <- 1/4
inpsim <- sim.spict(inp)
fit <- fit.spict(inpsim)

suppressWarnings(fit <- hindcast(fit, nyears = 7))


test_this("2.1. Decreasing length of index observations",
          (nobsi <- sapply(fit$hindcast, function(x) length(unlist(x$inp$obsI))))
          )

test_this("2.2. Decreasing length of index use indicator",
          (nobsiind <- sapply(fit$hindcast, function(x) length(which(x$inp$iuse == TRUE))))
          )

test_this("2.3. Except for first year (ref), difference corresponds to number of indices in this year:",
          (nobsi - nobsiind)
          )



header("3. Missing index observations at end of time series")
nt <- 50
inp <- list(nseasons = 1)
inp$timeC <- seq(0, nt - 1 / inp$nseasons, by = 1)
inp$timeI <- list(seq(0.1, nt-4, by = 1),
                  seq(0.8, nt-4, by = 1))
inp$ini <- list(logK = log(1000), logm=log(800), logq = log(1), logn = log(2),
                logbkfrac = log(0.9), logsdf = log(0.3), logF0 = log(0.2))
inp$dteuler <- 1/4
inpsim <- sim.spict(inp)
fit <- fit.spict(inpsim)

suppressWarnings(fit <- hindcast(fit, nyears = 7))

test_this("3.1. Decreasing length of index observations",
          (nobsi <- sapply(fit$hindcast, function(x) length(unlist(x$inp$obsI))))
          )

test_this("3.2. Decreasing length of index use indicator",
          (nobsiind <- sapply(fit$hindcast, function(x) length(which(x$inp$iuse == TRUE))))
          )

test_this("3.3. Except for first year (ref), difference corresponds to number of indices in this year:",
          (nobsi - nobsiind)
          )

mase1 <- calc.mase(fit)
mase2 <- plotspict.hindcast(fit)

test_this("3.4. The results from plotspict.hindcast and mase should be the same",
          all.equal(mase1, mase2)
          )


header("4. Missing index observations in middle of time series")
set.seed(12345)
nt <- 50
inp <- list(nseasons = 1)
inp$timeC <- seq(0, nt - 1 / inp$nseasons, by = 1)
inp$timeI <- list(seq(0.1, nt, by = 1))
inp$timeI[[1]] <- inp$timeI[[1]][-(nt-2)]
inp$ini <- list(logK = log(1000), logm=log(800), logq = log(1), logn = log(2),
                logbkfrac = log(0.9), logsdf = log(0.3), logF0 = log(0.4))
inp$dteuler <- 1/4
inpsim <- sim.spict(inp)
fit <- fit.spict(inpsim)

suppressWarnings(fit <- hindcast(fit, nyears = 5))

test_this("4.1. All peels converged:",
          sapply(fit$hindcast, function(x) x$opt$convergence)
          )

test_this("4.2. Decreasing length of index observations",
          (nobsi <- sapply(fit$hindcast, function(x) length(unlist(x$inp$obsI))))
          )

test_this("4.3. Decreasing length of index use indicator",
          (nobsiind <- sapply(fit$hindcast, function(x) length(which(x$inp$iuse == TRUE))))
          )

test_this("4.4. Except for first year (ref), difference corresponds to number of indices in this year:",
          (nobsi - nobsiind)
          )

test_this("4.5. MASE is calculated correctly plus warning about unequal spacing",
          capture.output(mase1 <- calc.mase(fit), type = "message")
          )

test_this("4.6. Plot is shown correctly plus warning about unequal spacing",
          capture.output(mase2 <- plotspict.hindcast(fit), type = "message")
          )

test_this("4.7. The results from plotspict.hindcast and mase should be the same",
          all.equal(mase1, mase2)
          )



header("5. Hindcasts exclude not converged runs")
inp <- pol$albacore
inp$dteuler <- 1/4
inp <- shorten.inp(pol$albacore, maxtime = 1980)

fit <- fit.spict(inp)
suppressWarnings(fit <- hindcast(fit, nyears = 8))

test_this("5.1. Convergence of last 4 peels is not achieved",
          sapply(fit$hindcast, function(x) x$opt$convergence)
          )

test_this("5.2. MASE is calculated correctly",
          capture.output(mase1 <- calc.mase(fit), type = "message")
          )

test_this("5.3. Plot is shown correctly",
          capture.output(mase2 <- plotspict.hindcast(fit), type = "message")
          )

test_this("5.4. The results from plotspict.hindcast and mase should be the same",
          all.equal(mase1, mase2)
          )

msg <- capture.output(plotspict.hindcast(fit), type = "message")
test_this("5.5. Plot shows message with excluded runs", length(msg) != 0)



## Expect error - not a fitted object
header("6. Error expectations")
test_this("6.1. No spictcls object produces an error", hindcast(pol$albacore))

## Expect error - not converged run
set.seed(1234)
inp <- list(obsC = runif(10, 100, 1000),
            timeC = 1:10,
            obsI = runif(10, 0.5, 1.5),
            timeI = 1:10)
out("\n")

suppressWarnings(fitnc <- fit.spict(inp))
test_this("6.2. Not converged fitted object produces an error", hindcast(fitnc))
