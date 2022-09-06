##################################################
###      Check hindcasting functionality       ###
##################################################

suppressMessages(library(spict))

options(device = function() png(width = 1100, height = 700))

set.seed(12345)

source("../funcs.R")

header("Check hindcasting functionality", append = FALSE)

header("1. Expected basic functionality")

header("1.1 Default albacore data set and model")
inp <- pol$albacore
inp$dteuler <- 1/4
fit <- fit.spict(inp)
suppressWarnings(fit <- hindcast(fit, npeels = 5))

test_this("1.1.1 Decreasing length of index observations",
          (nobsi <- sapply(fit$hindcast[-1], function(x) length(x$inp$obsI[[1]])))
          )

test_this("1.1.2 Decreasing length of index use indicator",
          (nobsiind <- sapply(fit$hindcast[-1], function(x) length(which(x$inp$iuse == TRUE))))
          )

test_this("1.1.3 Difference corresponds to number of index obs remove in given year",
          cumsum(nobsi - nobsiind)
          )

test_this("1.1.4 MASE:",
          signif(calc.mase(fit)[,2],4)
          )

test_this("1.1.5 First peel have the same number of catch obs as baserun, then cont. decreasing",
          sapply(fit$hindcast, function(x) length(x$inp$timeC))
          )

plotspict.hindcast(fit)


suppressWarnings(fit2 <- hindcast(fit, npeels = 5, peel.dtc = TRUE))

test_this("1.1.6 For an annual model MASE with peel.dtc = TRUE should be the same",
          signif(calc.mase(fit)[,2],4) == signif(calc.mase(fit2)[,2],4)
          )


header("1.2 Albacore data set assuming midyear index")
inp <- pol$albacore
inp$dteuler <- 1/4
inp$timeI <- inp$timeI + 0.5
fit <- fit.spict(inp)
suppressWarnings(fit <- hindcast(fit, npeels = 5))

test_this("1.2.1 Decreasing length of index observations",
          (nobsi <- sapply(fit$hindcast[-1], function(x) length(x$inp$obsI[[1]])))
          )

test_this("1.2.2 Decreasing length of index use indicator",
          (nobsiind <- sapply(fit$hindcast[-1], function(x) length(which(x$inp$iuse == TRUE))))
          )

test_this("1.2.3 Difference corresponds to number of index obs remove in given year",
          cumsum(nobsi - nobsiind)
          )

test_this("1.2.4 MASE:",
          signif(calc.mase(fit)[,2],4)
          )

test_this("1.2.5 First peel have the same number of catch obs as baserun, then cont. decreasing",
          sapply(fit$hindcast, function(x) length(x$inp$timeC))
          )

plotspict.hindcast(fit)

suppressWarnings(fit2 <- hindcast(fit, npeels = 5, peel.dtc = TRUE))

test_this("1.2.6 For an annual model MASE with peel.dtc = TRUE should be the same",
          signif(calc.mase(fit)[,2],4) == signif(calc.mase(fit2)[,2],4)
          )




header("1.3 Albacore data set assuming midyear index and last catch missing")
inp <- pol$albacore
inp$dteuler <- 1/4
inp$timeI <- inp$timeI + 0.5
inp$timeC <- inp$timeC[-length(inp$timeC)]
inp$obsC <- inp$obsC[-length(inp$obsC)]
fit <- fit.spict(inp)
suppressWarnings(fit <- hindcast(fit, npeels = 5))

test_this("1.3.1 Decreasing length of index observations",
          (nobsi <- sapply(fit$hindcast[-1], function(x) length(x$inp$obsI[[1]])))
          )

test_this("1.3.2 Decreasing length of index use indicator",
          (nobsiind <- sapply(fit$hindcast[-1], function(x) length(which(x$inp$iuse == TRUE))))
          )

test_this("1.3.3 Difference corresponds to number of index obs remove in given year",
          cumsum(nobsi - nobsiind)
          )

test_this("1.3.4 MASE:",
          signif(calc.mase(fit)[,2],4)
          )

test_this("1.3.5 First two peels have the same number of catch obs as baserun, then cont. decreasing",
          sapply(fit$hindcast, function(x) length(x$inp$timeC))
          )

plotspict.hindcast(fit)

suppressWarnings(fit2 <- hindcast(fit, npeels = 5, peel.dtc = TRUE))

test_this("1.3.6 For an annual model MASE with peel.dtc = TRUE should be the same",
          signif(calc.mase(fit)[,2],4) == signif(calc.mase(fit2)[,2],4)
          )




header("1.4 Albacore data set assuming midyear index and last indices missing")
inp <- pol$albacore
inp$dteuler <- 1/4
inp$timeI <- inp$timeI + 0.5
inp$timeI <- inp$timeI[-((length(inp$timeI)-2):length(inp$timeI))]
inp$obsI <- inp$obsI[-((length(inp$timeI)-2):length(inp$timeI))]
fit <- fit.spict(inp)
suppressWarnings(fit <- hindcast(fit, npeels = 7))

test_this("1.4.1 Decreasing length of index observations",
          (nobsi <- sapply(fit$hindcast[-1], function(x) length(x$inp$obsI[[1]])))
          )

test_this("1.4.2 Decreasing length of index use indicator",
          (nobsiind <- sapply(fit$hindcast[-1], function(x) length(which(x$inp$iuse == TRUE))))
          )

test_this("1.4.3 Difference corresponds to number of index obs remove in given year",
          cumsum(nobsi - nobsiind)
          )

test_this("1.4.4 MASE:",
          signif(calc.mase(fit)[,2],4)
          )

test_this("1.4.5 First peel have the same number of catch obs as baserun, then cont. decreasing",
          sapply(fit$hindcast, function(x) length(x$inp$timeC))
          )

plotspict.hindcast(fit)

suppressWarnings(fit2 <- hindcast(fit, npeels = 7, peel.dtc = TRUE))

test_this("1.4.6 For an annual model MASE with peel.dtc = TRUE should be the same",
          signif(calc.mase(fit)[,2],4) == signif(calc.mase(fit2)[,2],4)
          )



header("1.5 Albacore data set assuming midyear index and intermediate indices missing")
inp <- pol$albacore
inp$dteuler <- 1/4
inp$timeI <- inp$timeI + 0.5
inp$timeI <- inp$timeI[-((length(inp$timeI)-4):(length(inp$timeI)-2))]
inp$obsI <- inp$obsI[-((length(inp$timeI)-4):(length(inp$timeI)-2))]
fit <- fit.spict(inp)
suppressWarnings(fit <- hindcast(fit, npeels = 7))

test_this("1.5.1 Decreasing length of index observations",
          (nobsi <- sapply(fit$hindcast[-1], function(x) length(x$inp$obsI[[1]])))
          )

test_this("1.5.2 Decreasing length of index use indicator",
          (nobsiind <- sapply(fit$hindcast[-1], function(x) length(which(x$inp$iuse == TRUE))))
          )

test_this("1.5.3 Difference corresponds to number of index obs remove in given year",
          cumsum(nobsi - nobsiind)
          )

test_this("1.5.4 MASE:",
          signif(calc.mase(fit)[,2],4)
          )

test_this("1.5.5 First peel have the same number of catch obs as baserun, then cont. decreasing",
          sapply(fit$hindcast, function(x) length(x$inp$timeC))
          )

plotspict.hindcast(fit)

suppressWarnings(fit2 <- hindcast(fit, npeels = 7, peel.dtc = TRUE))

test_this("1.5.6 For an annual model MASE with peel.dtc = TRUE should be the same",
          signif(calc.mase(fit)[,2],4) == signif(calc.mase(fit2)[,2],4)
          )



header("1.6 Albacore data set and last index is exactly at end of last catch interval")
inp <- pol$albacore
inp$dteuler <- 1/4
inp <- shorten.inp(pol$albacore, maxtime = 1985)

fit <- fit.spict(inp)
suppressWarnings(fit <- hindcast(fit, npeels = 5))

test_this("1.6.1 Decreasing length of index observations",
          (nobsi <- sapply(fit$hindcast[-1], function(x) length(x$inp$obsI[[1]])))
          )

test_this("1.6.2 Decreasing length of index use indicator",
          (nobsiind <- sapply(fit$hindcast[-1], function(x) length(which(x$inp$iuse == TRUE))))
          )

test_this("1.6.3 Difference corresponds to number of index obs remove in given year",
          cumsum(nobsi - nobsiind)
          )

test_this("1.6.4 MASE:",
          signif(calc.mase(fit)[,2],4)
          )

test_this("1.6.5 First two peels have the same number of catch obs as baserun, then cont. decreasing",
          sapply(fit$hindcast, function(x) length(x$inp$timeC))
          )

plotspict.hindcast(fit)

suppressWarnings(fit2 <- hindcast(fit, npeels = 5, peel.dtc = TRUE))

test_this("1.6.6 For an annual model MASE with peel.dtc = TRUE should be the same",
          signif(calc.mase(fit)[,2],4) == signif(calc.mase(fit2)[,2],4)
          )




header("1.7 Default hake data set (4 peels)")
inp <- pol$hake
inp$dteuler <- 1/4
fit <- fit.spict(inp)
suppressWarnings(fit <- hindcast(fit, npeels = 4))

test_this("1.7.1 Decreasing length of index observations",
          (nobsi <- sapply(fit$hindcast[-1], function(x) length(x$inp$obsI[[1]])))
          )

test_this("1.7.2 Decreasing length of index use indicator",
          (nobsiind <- sapply(fit$hindcast[-1], function(x) length(which(x$inp$iuse == TRUE))))
          )

test_this("1.7.3 Difference corresponds to number of index obs remove in given year",
          cumsum(nobsi - nobsiind)
          )

test_this("1.7.4 MASE:",
          signif(calc.mase(fit)[,2],4)
          )

test_this("1.7.5 First peel have the same number of catch obs as baserun, then cont. decreasing",
          sapply(fit$hindcast, function(x) length(x$inp$timeC))
          )

plotspict.hindcast(fit)


header("1.8 Default hake data set (7 peels)")
inp <- pol$hake
inp$dteuler <- 1/4
fit <- fit.spict(inp)
suppressWarnings(fit <- hindcast(fit, npeels = 7))

test_this("1.8.1 Decreasing length of index observations",
          (nobsi <- sapply(fit$hindcast[-1], function(x) length(x$inp$obsI[[1]])))
          )

test_this("1.8.2 Decreasing length of index use indicator",
          (nobsiind <- sapply(fit$hindcast[-1], function(x) length(which(x$inp$iuse == TRUE))))
          )

test_this("1.8.3 Difference corresponds to number of index obs remove in given year",
          cumsum(nobsi - nobsiind)
          )

test_this("1.8.4 MASE:",
          signif(calc.mase(fit)[,2],4)
          )

test_this("1.8.5 First peel have the same number of catch obs as baserun, then cont. decreasing",
          sapply(fit$hindcast, function(x) length(x$inp$timeC))
          )

plotspict.hindcast(fit)


header("1.9 Hake data set and more stable model (7 peels)")
inp <- pol$hake
inp$dteuler <- 1/4
inp$timeI <- inp$timeI + 0.5
inp$priors$logalpha <- c(0,0,0)
inp$priors$logbeta <- c(0,0,0)
inp$priors$logsdb <- c(log(0.1), 1, 1)
inp$ini$logn <- log(2)
inp$phases$logn <- -1
fit <- fit.spict(inp)
suppressWarnings(fit <- hindcast(fit, npeels = 7))

test_this("1.9.1 Decreasing length of index observations",
          (nobsi <- sapply(fit$hindcast[-1], function(x) length(x$inp$obsI[[1]])))
          )

test_this("1.9.2 Decreasing length of index use indicator",
          (nobsiind <- sapply(fit$hindcast[-1], function(x) length(which(x$inp$iuse == TRUE))))
          )

test_this("1.9.3 Difference corresponds to number of index obs remove in given year",
          cumsum(nobsi - nobsiind)
          )

test_this("1.9.4 MASE:",
          signif(calc.mase(fit)[,2],4)
          )

test_this("1.9.5 First peel have the same number of catch obs as baserun, then cont. decreasing",
          sapply(fit$hindcast, function(x) length(x$inp$timeC))
          )

plotspict.hindcast(fit)



header("2. Multiple indices with different lengths")
set.seed(12345)
nt <- 50
inp <- list(nseasons = 1)
inp$timeC <- seq(0, nt - 1 / inp$nseasons, by = 1)
inp$timeI <- list(seq(0.1, nt, by = 1),
                  seq(0.8, nt-3, by = 1))
inp$ini <- list(logK = log(1000), logm=log(600), logq = log(1), logn = log(2),
                logbkfrac = log(0.9), logsdf = log(0.5), logF0 = log(0.5), logsdi = c(log(0.4),log(0.1)))
inp$dteuler <- 1/4
inpsim <- sim.spict(inp)
fit <- fit.spict(inpsim)

suppressWarnings(fit <- hindcast(fit, npeels = 7))


test_this("2.1 Decreasing length of index observations",
          (nobsi <- sapply(fit$hindcast[-1], function(x) length(unlist(x$inp$obsI))))
          )

test_this("2.2 Decreasing length of index use indicator",
          (nobsiind <- sapply(fit$hindcast[-1], function(x) length(which(x$inp$iuse == TRUE))))
          )

test_this("2.3 Two indices removed each year:",
          cumsum(nobsi - nobsiind)
          )

test_this("2.4 MASE:",
          signif(calc.mase(fit)[,2],4)
          )

plotspict.hindcast(fit, legend.pos = "bottomright")

test_this("2.5 First peel have the same number of catch obs as baserun, then cont. decreasing",
          sapply(fit$hindcast, function(x) length(x$inp$timeC))
          )


header("3. Missing index observations at end of time series")
set.seed(12345)
nt <- 40
inp <- list(nseasons = 1)
inp$timeC <- seq(0, nt - 1 / inp$nseasons, by = 1)
inp$timeI <- list(seq(0.1, nt-4, by = 1),
                  seq(0.8, nt-4, by = 1))
inp$ini <- list(logK = log(1000), logm=log(400), logq = log(1), logn = log(2),
                logbkfrac = log(0.9), logsdf = log(0.3), logF0 = log(0.01))
inp$dteuler <- 1/4
inpsim <- sim.spict(inp)
fit <- fit.spict(inpsim)

suppressWarnings(fit <- hindcast(fit, npeels = 7))

test_this("3.1 Decreasing length of index observations",
          (nobsi <- sapply(fit$hindcast[-1], function(x) length(unlist(x$inp$obsI))))
          )

test_this("3.2 Decreasing length of index use indicator",
          (nobsiind <- sapply(fit$hindcast[-1], function(x) length(which(x$inp$iuse == TRUE))))
          )

test_this("3.3 Last 4 yars, no indices, then two per year removed:",
          cumsum(nobsi - nobsiind)
          )

test_this("3.4 MASE:",
          signif(calc.mase(fit)[,2],4)
          )

plotspict.hindcast(fit)

test_this("3.5 First peel have the same number of catch obs as baserun, then cont. decreasing",
          sapply(fit$hindcast, function(x) length(x$inp$timeC))
          )


header("4. Missing index observations in middle of time series")
set.seed(123456)
nt <- 50
inp <- list(nseasons = 1)
inp$timeC <- seq(0, nt - 1 / inp$nseasons, by = 1)
inp$timeI <- list(seq(0.1, nt, by = 1))
inp$timeI[[1]] <- inp$timeI[[1]][-(nt-c(2:3))]
inp$ini <- list(logK = log(1000), logm=log(800), logq = log(1), logn = log(2),
                logbkfrac = log(0.7), logsdf = log(0.3), logF0 = log(0.2), logsdi = log(0.05))
inp$dteuler <- 1/4
inpsim <- sim.spict(inp)
fit <- fit.spict(inpsim)

suppressWarnings(fit <- hindcast(fit, npeels = 5))

test_this("4.1 All peels converged:",
          sapply(fit$hindcast, function(x) x$opt$convergence)
          )

test_this("4.2 Decreasing length of index observations",
          (nobsi <- sapply(fit$hindcast[-1], function(x) length(unlist(x$inp$obsI))))
          )

test_this("4.3 Decreasing length of index use indicator",
          (nobsiind <- sapply(fit$hindcast[-1], function(x) length(which(x$inp$iuse == TRUE))))
          )

test_this("4.4 Except for first year (ref), difference corresponds to number of indices in this year:",
          cumsum(nobsi - nobsiind)
          )

test_this("4.5 MASE is calculated correctly plus warning about unequal spacing",
          capture.output(signif(calc.mase(fit)[,2],4), type = "message")
          )

test_this("4.6 Plot is shown correctly plus warning about unequal spacing",
          capture.output(plotspict.hindcast(fit), type = "message")
          )

test_this("4.7 First peel have the same number of catch obs as baserun, then cont. decreasing",
          sapply(fit$hindcast, function(x) length(x$inp$timeC))
          )



header("5. Hindcasts exclude not converged runs")
inp <- pol$albacore
inp$dteuler <- 1/4
inp <- shorten.inp(pol$albacore, maxtime = 1978)

fit <- fit.spict(inp)
suppressWarnings(fit <- hindcast(fit, npeels = 6))

test_this("5.1 Convergence of last 3 peels is not achieved",
          sapply(fit$hindcast[-1], function(x) x$opt$convergence)
          )

test_this("5.2 MASE is calculated correctly",
          capture.output(signif(calc.mase(fit)[,2],4), type = "message")
          )

test_this("5.3 Plot is shown correctly",
          capture.output(plotspict.hindcast(fit), type = "message")
          )

msg <- capture.output(plotspict.hindcast(fit), type = "message")
test_this("5.5 Plot shows message with excluded runs", length(msg) != 0)


## Expect error - not a fitted object
header("6. Error expectations")
test_this("6.1 No spictcls object produces an error", hindcast(pol$albacore))

## Expect error - not converged run
set.seed(1234)
inp <- list(obsC = runif(10, 100, 1000),
            timeC = 1:10,
            obsI = runif(10, 0.5, 1.5),
            timeI = 1:10)
out("\n")

suppressWarnings(fitnc <- fit.spict(inp))
test_this("6.2 Not converged fitted object produces an error", hindcast(fitnc))



header("7. Seasonal data with multiple indices with various time series lengths plus catch in last half year is missing")
set.seed(1234)
nt <- 50
inp <- list(nseasons = 4)
inp$timeC <- seq(0, nt - 3 / inp$nseasons , by = 1/inp$nseasons)
inp$timeI <- list(seq(0.78, nt, by = 1),
                  seq(0.35, nt - 1, by = 1),
                  seq(0.56, nt - 20, by = 1),
                  seq(0.11, nt - 2, by = 1))
inp$ini <- list(logK = log(1000), logm = log(400), logq = log(1), logn = log(2),
                logbkfrac = log(0.9), logsdf = log(0.3),
                logF0 = log(0.3), logsdi = c(log(0.05),log(0.1),log(0.1),log(0.4)))
inp$dteuler <- 1/8
inpsim <- sim.spict(inp)
fit <- fit.spict(inpsim)

header("7.1 Annual peels")

suppressWarnings(fit <- hindcast(fit, npeels = 7))

test_this("7.1.1 All peels converged:",
          sapply(fit$hindcast, function(x) x$opt$convergence)
          )

test_this("7.1.2 Decreasing length of index observations",
          (nobsi <- sapply(fit$hindcast[-1], function(x) length(unlist(x$inp$obsI))))
          )

test_this("7.1.3 Decreasing length of index use indicator",
          (nobsiind <- sapply(fit$hindcast[-1], function(x) length(which(x$inp$iuse == TRUE))))
          )

test_this("7.1.4 Except for first year (ref), difference corresponds to number of index obs remove in given year",
          cumsum(nobsi - nobsiind)
          )

test_this("7.1.5 MASE is calculated correctly plus message that one index is outside of peels",
          capture.output(signif(calc.mase(fit)[,2],4), type = "message")
          )

test_this("7.1.6 Plot is shown correctly plus message that one index is outside of peels",
          capture.output(plotspict.hindcast(fit, legend.ncol = 2), type = "message")
          )

test_this("7.1.7 First peel have the same number of catch obs as baserun, then cont. decreasing (first 2 seasons, then 4 seasons)",
          sapply(fit$hindcast, function(x) length(x$inp$timeC))
          )


header("7.2 Seasonal peels")

suppressWarnings(fit <- hindcast(fit, npeels = 20, peel.dtc = TRUE))


test_this("7.2.1 All peels converged:",
          sapply(fit$hindcast, function(x) x$opt$convergence)
          )

test_this("7.2.2 Decreasing length of index observations",
          (nobsi <- sapply(fit$hindcast[-1], function(x) length(unlist(x$inp$obsI))))
          )

test_this("7.2.3 Decreasing length of index use indicator",
          (nobsiind <- sapply(fit$hindcast[-1], function(x) length(which(x$inp$iuse == TRUE))))
          )

test_this("7.2.4 Except for first year (ref), difference corresponds to number of index obs removed in given year",
          cumsum(nobsi - nobsiind)
          )

test_this("7.5 MASE is calculated correctly plus message that one index is outside of peels",
          capture.output(signif(calc.mase(fit)[,2],4), type = "message")
          )

test_this("7.2.6 Plot is shown correctly plus message that one index is outside of peels",
          capture.output(plotspict.hindcast(fit, legend.ncol = 4), type = "message")
          )

test_this("7.2.7 First three peels have the same number of catch obs as baserun, then cont. decreasing ",
          sapply(fit$hindcast, function(x) length(x$inp$timeC))
          )
