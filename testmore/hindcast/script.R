##################################################
###      Check hindcasting functionality       ###
##################################################

suppressMessages(library(spict))
options(device = png)

source("../funcs.R")

header("Check hindcasting functionality", append = FALSE)

header("1. Hindcasts exclude not converged runs")
inp <- pol$albacore
inp$dteuler <- 1/4
inp <- shorten.inp(pol$albacore, maxtime = 1980)

fit <- fit.spict(inp)
suppressWarnings(fit <- hindcast(fit, nyears = 8))

test_this("1.1. Convergence of last three peels is not achieved",
          sapply(fit$hindcast, function(x) x$opt$convergence))

test_this("1.2. MASE is calculated correctly",
          capture.output(mase1 <- calc.mase(fit),
                         type = "message"))

test_this("1.3. Plot is shown correctly",
          capture.output(mase2 <- plotspict.hindcast(fit),
                         type = "message"))

test_this("1.4. The results from plotspict.hindcast and mase should be the same",
          all.equal(mase1, mase2))

msg <- capture.output(plotspict.hindcast(fit), type = "message")
test_this("1.5. Plot shows message with excluded runs", length(msg) != 0)


## Expect error - not a fitted object
header("2. Error expectations")
test_this("2.1. No spictcls object produces an error", hindcast(pol$albacore))

## Expect error - not converged run
set.seed(1234)
inp <- list(obsC = runif(10, 100, 1000),
            timeC = 1:10,
            obsI = runif(10, 0.5, 1.5),
            timeI = 1:10)
out("\n")

suppressWarnings(fitnc <- fit.spict(inp))
test_this("2.2. Not converged fitted object produces an error", hindcast(fitnc))
