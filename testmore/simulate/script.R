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
inp$sim.random.effects <- TRUE
set.seed(123)
fit1 <- fit.spict(inp)
fit1$inp$sim.random.effects

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

test_this("Test 1.5: R and TMB simulated catch observations the same?", {
  all(round(sim1a$obsC,2) == round(sim1b$obsC,2))
})

test_this("Test 1.6: R and TMB simulated index observations the same?", {
  all(round(sim1a$obsI[[1]],2) == round(sim1b$obsI[[1]],2))
})




header("2: Testing sim.random.effects", append = TRUE)
## -------------------

inp <- pol$hake
inp$dteuler <- 1/4
inp$sim.random.effects <- FALSE ## can be set before or after fitting
set.seed(123)
fit2 <- fit.spict(inp)
set.seed(123)
sim2a <- sim.spict(fit2, use.tmb = FALSE)
set.seed(123)
sim2b <- sim.spict(fit2, use.tmb = TRUE)

test_this("Test 2.0: Do R and TMB give the same values when sim.random.effects = FALSE?", {
    all(round(sim2a$obsC,1) == round(sim2b$obsC,1))
})

test_this("Test 2.1: Are simulated logB with and without sim.random.effects different?", {
  all(sim1a$true$B != sim2a$true$B)
})

test_this("Test 2.2: Are simulated logB with and without sim.random.effects different (TMB)?", {
  all(sim1b$true$B != sim2b$true$B)
})

test_this("Test 2.3: Do all B states correspond to estimate Bs with sim.random.effects == TRUE?", {
  all(sim2a$true$B == as.vector(get.par("logB",fit1, exp = TRUE)[,2]))
})

test_this("Test 2.4: Do all B states correspond to estimate Bs with sim.random.effects == TRUE (TMB)?", {
  all(sim2b$true$B == as.vector(get.par("logB",fit1, exp = TRUE)[,2]))
})

test_this("Test 2.5: Are simulated logF with and without sim.random.effects different?", {
  all(sim1a$true$F != sim2a$true$F)
})

test_this("Test 2.6: Are simulated logF with and without sim.random.effects different (TMB)?", {
  all(sim1b$true$F != sim2b$true$F)
})

test_this("Test 2.7: Do all F states correspond to estimate Fs with sim.random.effects == TRUE?", {
  all(sim2a$true$F == as.vector(get.par("logF",fit1, exp = TRUE)[,2]))
})

test_this("Test 2.8: Do all F states correspond to estimate Fs with sim.random.effects == TRUE (TMB)?", {
  all(sim2b$true$F == as.vector(get.par("logF",fit1, exp = TRUE)[,2]))
})


## settting sim.random.effects after fitting.
fit1$inp$sim.random.effects

set.seed(123)
fit1$inp$sim.random.effects <- FALSE
sim2a2 <- sim.spict(fit1, use.tmb = FALSE)

set.seed(123)
fit1$inp$sim.random.effects <- FALSE
sim2b2 <- sim.spict(fit1, use.tmb = TRUE)

test_this("Test 2.9: Do all F states correspond to estimate Fs with sim.random.effects == TRUE (even when sim.random.effects is used after fitting)?", {
    all(sim2a2$true$F == as.vector(get.par("logF",fit1, exp = TRUE)[,2]))
})

test_this("Test 2.10: Do all F states correspond to estimate Fs with sim.random.effects == TRUE (even when sim.random.effects is used after fitting, TMB)?", {
    all(sim2b2$true$F == as.vector(get.par("logF",fit1, exp = TRUE)[,2]))
})





header("3: Testing sim.priors", append = TRUE)
## -------------------

out("to come")

## inp <- pol$hake
## inp$dteuler <- 1/4
## inp$sim.priors <- 1
## set.seed(123)
## fit3 <- fit.spict(inp)
## set.seed(123)
## sim3 <- sim.spict(fit3, use.tmb = TRUE)

## test_this("Test 3.1: Are simulated logB with and without sim.priors different?", {
##   all(sim1b$true$B != sim3$true$B)
## })

## test_this("Test 3.2: Are simulated logF with and without sim.priors different?", {
##   all(sim1b$true$F[-1] != sim3$true$F[-1])
## })



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
    logphi = log(c(0.05, 0.1, 1.8)),
    logbkfrac = log(0.9), logsdf = log(0.1),
    logF0 = log(0.01))
inp$dteuler <- 1/4
inp$seasontype <- 1
## simulate
set.seed(123)
inp4a1 <- sim.spict(inp, use.tmb = FALSE)
set.seed(123)
inp4b1 <- sim.spict(inp, use.tmb = TRUE)

test_this("Test 4.1: Catch observations simulated with seasontype = 1", {
    round(inp4a1$obsC,2)
})

test_this("Test 4.2: Catch observations simulated with seasontype = 1 (TMB)", {
    round(inp4b1$obsC,2)
})

test_this("Test 4.3: R and TMB catches the same?", {
    all(round(inp4a1$obsC,2) == round(inp4b1$obsC,2))
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
set.seed(123)
inp4a2 <- sim.spict(inp, use.tmb = FALSE)
set.seed(123)
inp4b2 <- sim.spict(inp, use.tmb = TRUE)

test_this("Test 4.3: Catch observations simulated with seasontype = 2", {
  round(inp4a2$obsC,2)
})

test_this("Test 4.4: Catch observations simulated with seasontype = 2 (TMB)", {
  round(inp4b2$obsC,2)
})

## Values are close but not exactly the same.


## seasontype = 3
## -------------------------
set.seed(123)
nt <- 100
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
    logphi = log(c(0.05, 0.1, 1.8)),
    logF0 = log(0.01))
inp$dteuler <- 1/4
inp$seasontype <- 3
## simulate
set.seed(503)
inp4a3 <- sim.spict(inp, use.tmb = FALSE)
set.seed(503)
inp4b3 <- sim.spict(inp, use.tmb = TRUE)

test_this("Test 4.5: Catch observations simulated with seasontype = 3", {
  round(inp4a3$obsC,2)
})

test_this("Test 4.6: Catch observations simulated with seasontype = 3 (TMB)", {
  round(inp4b3$obsC,2)
})


## from seasonal fitted object
fit4 <- fit.spict(inp4b3)


fit4$inp$sim.random.effects <- FALSE
set.seed(123)
inp4a4 <- sim.spict(fit4, use.tmb = FALSE)
set.seed(123)
inp4b4 <- sim.spict(fit4, use.tmb = TRUE)


test_this("Test 4.7: R and TMB catch observations simulated with seasontype = 3 the same (with sim.random.effects = FALSE)?", {
    all(round(inp4a4$obsC,2) == round(inp4b4$obsC,2))
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
set.seed(123)
sim5a <- sim.spict(inp, use.tmb = FALSE)
set.seed(123)
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
sim6a <- sim.spict(inpI, use.tmb = FALSE)
sim6b <- sim.spict(inpI, use.tmb = TRUE)

test_this("Test 6.1: Four indices with four catchability coefficients.", {
    round(do.call(cbind, sim6a$obsI),2)
})

test_this("Test 6.2: Four indices with four catchability coefficients (TMB).", {
    round(do.call(cbind, sim6b$obsI),2)
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
sim6c <- sim.spict(inpI, use.tmb = TRUE)

test_this("Test 6.2: Four indices with single catchability coefficients.", {
    round(do.call(cbind, sim6c$obsI),2)
})





header("7: TMB's check consistency", append = TRUE)
## -------------------
sim6b$priors$logn <- c(0,0,0)
sim6b$priors$logalpha <- c(0,0,0)
sim6b$priors$logbeta <- c(0,0,0)
sim6b$simulate <- FALSE
fit7 <- fit.spict(sim6b)


nrep <- 100
checkcon <- TMB::checkConsistency(fit7$obj, n = nrep)
tmp <- summary(checkcon, na.rm=TRUE)
tmp

test_this("Test 7.1: Pvalue of joint likelihood", {
    round(tmp$joint$p.value,3)
})

test_this("Test 7.2: Biases of joint likelihood", {
    sapply(tmp$joint$bias, function(x) round(x,3))
})

test_this("Test 7.3: Pvalue of marginal likelihood", {
    round(tmp$marginal$p.value,3)
})

test_this("Test 7.4: Biases of marginal likelihood", {
    sapply(tmp$marginal$bias, function(x) round(x,3))
})
