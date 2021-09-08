## This script tests the TMB simulation functionality.
## Created: 03/09/2021
##
## Comments: Only simulation of 3 default priors possible: logn, logalpha,
##   logbeta; simulation code for other priors remains to be added. (Some priors
##   require significant more code)


## Load spict
## -------------------
require(spict)


## Test functions
## -------------------
source("../funcs.R")



header("1: Simple tests", append = FALSE)
## -------------------

inp <- pol$hake
inp$dteuler <- 1/4
set.seed(123)
fit1 <- fit.spict(inp)


set.seed(123)
sim1a <- sim.spict(fit1, use.tmb = FALSE)
set.seed(123)
sim1b <- sim.spict(fit1, use.tmb = TRUE)

test_this("Test 1.1: R simulated catch observations", {
  round(sim1a$obsC, 2)
})

test_this("Test 1.2: TMB simulated catch observations", {
  round(sim1b$obsC, 2)
})

test_this("Test 1.3: R simulated index observations", {
  round(sim1a$obsI[[1]], 2)
})

test_this("Test 1.4: TMB simulated index observations", {
  round(sim1b$obsI[[1]], 2)
})




header("2: Testing sim.random.effects", append = TRUE)
## -------------------

inp <- pol$hake
inp$dteuler <- 1/4
inp$sim.random.effects <- FALSE
set.seed(123)
fit2 <- fit.spict(inp)
set.seed(123)
sim2 <- sim.spict(fit2, use.tmb = TRUE)


test_this("Test 2.1: Are simulated logB with and without sim.random.effects different (except first value)?", {
  all(sim1b$true$B[-1] != sim2$true$B[-1])
})

test_this("Test 2.2: Do all B states correspond to estimate Bs with sim.random.effects == TRUE?", {
  all(sim2$true$B == as.vector(get.par("logB",fit1, exp = TRUE)[,2]))
})

test_this("Test 2.3: Are simulated logF with and without sim.random.effects different (except first value)?", {
  all(sim1b$true$F[-1] != sim2$true$F[-1])
})

test_this("Test 2.4: Do all F states correspond to estimate Fs with sim.random.effects == TRUE?", {
  all(sim2$true$F == as.vector(get.par("logF",fit1, exp = TRUE)[,2]))
})




header("3: Testing sim.priors", append = TRUE)
## -------------------

inp <- pol$hake
inp$dteuler <- 1/4
inp$sim.priors <- 1
set.seed(123)
fit3 <- fit.spict(inp)
set.seed(123)
sim3 <- sim.spict(fit3, use.tmb = TRUE)

test_this("Test 3.1: Are simulated logB with and without sim.priors different?", {
  all(sim1b$true$B != sim3$true$B)
})

test_this("Test 3.2: Are simulated logF with and without sim.priors different?", {
  all(sim1b$true$F[-1] != sim3$true$F[-1])
})



header("4: Tests with seasonal data", append = TRUE)
## -------------------

## seasontype = 1
## -------------------------
set.seed(123)
nt <- 40
inp <- list(nseasons = 4)
inp$timeC <- seq(0, nt - 1 / inp$nseasons, by = 1 / inp$nseasons)
inp$timeI <- list(
    seq(1/4, nt - 1 / inp$nseasons, by = 1 / inp$nseasons),
    seq(3/4, nt - 1 / inp$nseasons, by = 1 / inp$nseasons)
)
inp$ini <- list(
    logK = log(1000),
    logm=log(900),
    logq = c(log(0.3),log(1)),
    logn = log(1.4),
    logbkfrac = log(0.9), logsdf = log(0.1),
    logF0 = log(0.01))
inp$dteuler <- 1/4
inp$seasontype <- 1
## simulate
inp <- sim.spict(inp, use.tmb = TRUE)


test_this("Test 4.1: Catch observations simulated with seasontype = 1", {
  round(inp$obsC,2)
})


## seasontype = 2
## -------------------------
set.seed(123)
nt <- 40
inp <- list(nseasons = 4)
inp$timeC <- seq(0, nt - 1 / inp$nseasons, by = 1 / inp$nseasons)
inp$timeI <- list(
    seq(1/4, nt - 1 / inp$nseasons, by = 1 / inp$nseasons),
    seq(3/4, nt - 1 / inp$nseasons, by = 1 / inp$nseasons)
)
inp$ini <- list(
    logK = log(1000),
    logm=log(900),
    logq = c(log(0.3),log(1)),
    logn = log(1.4),
    logbkfrac = log(0.9), logsdf = log(0.1),
    logF0 = log(0.01))
inp$dteuler <- 1/4
inp$seasontype <- 2
## simulate
inp <- sim.spict(inp, use.tmb = TRUE)

test_this("Test 4.2: Catch observations simulated with seasontype = 2", {
  round(inp$obsC,2)
})



## seasontype = 3
## -------------------------
set.seed(123)
nt <- 40
inp <- list(nseasons = 4)
inp$timeC <- seq(0, nt - 1 / inp$nseasons, by = 1 / inp$nseasons)
inp$timeI <- list(
    seq(1/4, nt - 1 / inp$nseasons, by = 1 / inp$nseasons),
    seq(3/4, nt - 1 / inp$nseasons, by = 1 / inp$nseasons)
)
inp$ini <- list(
    logK = log(1000),
    logm=log(900),
    logq = c(log(0.3),log(1)),
    logn = log(1.4),
    logbkfrac = log(0.9), logsdf = log(0.1),
    logF0 = log(0.01))
inp$dteuler <- 1/4
inp$seasontype <- 3
## simulate
inp <- sim.spict(inp, use.tmb = TRUE)

test_this("Test 4.3: Catch observations simulated with seasontype = 3", {
  round(inp$obsC,2)
})


header("5: Simulating effort data", append = TRUE)
## -------------------

set.seed(123)
nt <- 100
inp <- list(nseasons = 1)
inp$timeC <- seq(0, nt - 1 / inp$nseasons, by = 1 / inp$nseasons)
inp$timeE <- seq(0, nt - 1 / inp$nseasons, by = 1 / inp$nseasons)
inp$ini <- list(logK = log(1000), logm=log(100), logqf = log(0.3), logn = log(2),
                 logbkfrac = log(0.9), logsdf = log(0.3), logF0 = log(0.1))
inp$dteuler <- 1/4

sim5a <- sim.spict(inp, use.tmb = FALSE)
sim5b <- sim.spict(inp, use.tmb = TRUE)

test_this("Test 5.1: Simulated effort", {
    round(sim5a$obsE,2)
})

test_this("Test 5.2: Simulated effort (with TMB)", {
    round(sim5b$obsE,2)
})


## sim.comm.cpue
inp$sim.comm.cpue <- TRUE
sim5c <- sim.spict(inp, use.tmb = FALSE)
sim5d <- sim.spict(inp, use.tmb = TRUE)

test_this("Test 5.3: Commercial CPUE observations", {
    round(sim5c$obsI,2)
})

test_this("Test 5.4: Commercial CPUE observations (with TMB)", {
    round(sim5d$obsI,2)
})



header("6: Simulating multiple indices", append = TRUE)
## -------------------
set.seed(123)
nt <- 40
inpI <- list(nseasons = 1)
inpI$timeC <- seq(0, nt - 1 / inpI$nseasons, by = 1 / inpI$nseasons)
inpI$timeI <- list(seq(1/4, nt - 1 / inpI$nseasons, by = 1 / inpI$nseasons),
                   seq(2/4, nt - 1 / inpI$nseasons, by = 1 / inpI$nseasons),
                   seq(3/4, nt - 1 / inpI$nseasons, by = 1 / inpI$nseasons),
                   seq(3.5/4, nt - 1 / inpI$nseasons, by = 1 / inpI$nseasons)
                   )
inpI$ini <- list(logK = log(1000), logm=log(100),
                 logq = c(log(0.3),log(1),log(3),log(0.3)),
                 logn = log(2),
                 logbkfrac = log(0.9), logsdf = log(0.3), logF0 = log(0.1))
inpI$dteuler <- 1/4
sim6a <- sim.spict(inpI, use.tmb = TRUE)

test_this("Test 6.1: Four indices with four catchability coefficients.", {
    round(do.call(cbind, sim6a$obsI),2)
})



## With one q provided for all indices
set.seed(123)
nt <- 40
inpI <- list(nseasons = 1)
inpI$timeC <- seq(0, nt - 1 / inpI$nseasons, by = 1 / inpI$nseasons)
inpI$timeI <- list(seq(1/4, nt - 1 / inpI$nseasons, by = 1 / inpI$nseasons),
                   seq(2/4, nt - 1 / inpI$nseasons, by = 1 / inpI$nseasons),
                   seq(3/4, nt - 1 / inpI$nseasons, by = 1 / inpI$nseasons),
                   seq(3.5/4, nt - 1 / inpI$nseasons, by = 1 / inpI$nseasons)
                   )
inpI$ini <- list(logK = log(1000), logm=log(100),
                 logq = log(1),
                 logn = log(2),
                 logbkfrac = log(0.9), logsdf = log(0.3), logF0 = log(0.1))
inpI$dteuler <- 1/4
sim6b <- sim.spict(inpI, use.tmb = TRUE)

test_this("Test 6.2: Four indices with single catchability coefficients.", {
    round(do.call(cbind, sim6b$obsI),2)
})





header("7: TMB's check consistency", append = TRUE)
## -------------------
nrep <- 100
checkcon <- TMB::checkConsistency(fit1$obj, n = nrep)
tmp <- summary(checkcon, na.rm=TRUE)

test_this("Test 7.1: Summary of join likelihood", {
    lapply(tmp$joint, function(x) round(x,2))
})

test_this("Test 7.2: Summary of marginal likelihood", {
    lapply(tmp$marginal, function(x) round(x,2))
})
