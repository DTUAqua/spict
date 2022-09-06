## Testing management functionality in uncommon situations
## Alex Kokkalis <alko@aqua.dtu.dk>
## 2022-03-07

## set seed
set.seed(123)

## load spict
# library(remotes)
# remotes::install_github("DTUAqua/spict/spict", ref = github_pull("148"))
require(spict)

## load testmore functions
source("../funcs.R")


header("1: Testing add.man.scenarios with mixed data", append = FALSE)


set.seed(1234)
nt <- 50
inp <- list(nseasons=4, splineorder=3)
s1 <- seq(0, nt/2-1, by=1)
s2 <- seq(nt/2, nt-1/inp$nseasons, by=1/inp$nseasons)
inp$timeC <- c(s1,s2)
inp$dtc <- c(rep(1, length(s1)),
             rep(0.25, length(s2)))
inp$timeI <- seq(0.1, nt-1/inp$nseasons, by=0.5)
inp$ini <- list(logK=log(1000), logm=log(123), logq=log(1), logn=log(2),
                logbkfrac=log(0.3), logsdf=log(0.22), logF0=log(0.4),
                logphi=log(c(0.3, 0.5, 1.8)))
inp$dteuler <- 1/8 ## Does not converge with 1/4
inpsim <- sim.spict(inp, verbose = TRUE)

header("1.1: long management period")
###############################

inpsim$maneval <- nt+5
inpsim$maninterval <- c(nt+1,nt+5)
inp <- check.inp(inpsim)

fit <- fit.spict(inp)

fit <- add.man.scenario(fit, "Fmsy")
fit <- add.man.scenario(fit, "C35", fractiles = list(catch = 0.35))

test_this("1.1.1: Fmsy and Fmsy with 35th fractile on catch", {
  sumspict.manage(fit)
})
  # plotspict.f(fit) ## Strange: F/Fmsy seems to be around 1.5 in the prediction
  # plotspict.ffmsy(fit)
  # plotspict.biomass(fit)
  # plotspict.catch(fit)
  # fit <- retro(fit, 3)
  # plotspict.retro(fit)


header("2: Testing add.man.scenarios with mixed data - annual last")


set.seed(1234)
nt <- 50
inp <- list(nseasons=4, splineorder=3)
s1 <- seq(0, 10, by=1)
s2 <- seq(11, nt-3-1/inp$nseasons, by=1/inp$nseasons)
s3 <- seq(nt-3, nt-1, by = 1)
inp$timeC <- c(s1,s2,s3)
inp$dtc <- c(rep(1, length(s1)),
             rep(0.25, length(s2)),
             rep(1, length(s3)))
inp$timeI <- seq(0.1, nt-1/inp$nseasons, by=0.5)
inp$ini <- list(logK=log(1000), logm=log(123), logq=log(1), logn=log(2),
                logbkfrac=log(0.7), logsdf=log(0.123), logF0=log(0.4321),
                logphi=log(c(0.3, 0.5, 1.8)))
inp$dteuler <- 1/8
inpsim <- sim.spict(inp)
## inpsim$dtc <- inp$dtc ## Very important, sim.spict drops it
inpsim$maneval <- nt+5
inpsim$maninterval <- c(nt,nt+5)

header("2.1: long management period")
###############################
inpsim$priors$logbkfrac <- c(log(0.7), 0.001, 1)
inpsim$priors$logn <- c(log(2), 0.001, 1)
inpsim <- check.inp(inpsim)

fit <- fit.spict(inpsim)


fit <- add.man.scenario(fit, "Fmsy")
fit <- add.man.scenario(fit, "C35", fractiles = list(catch = 0.35))
fit <- add.man.scenario(fit, "All35", fractiles = list(catch = 0.35, ffmsy=0.35, bbmsy=0.35))
fit <- add.man.scenario(fit, "FB35_noC", fractiles = list( ffmsy=0.35, bbmsy=0.35))



test_this("2.1.1: Fmsy and Fmsy with 35th fractile on catch", {
  sumspict.manage(fit)
})
