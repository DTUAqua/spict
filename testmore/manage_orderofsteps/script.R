## Testing order of operations
## T.K. Mildenberger <t.k.mildenberger@gmail.com>
## 10/12/2019

## set seed
set.seed(123)

## load spict
require(spict)

## load testmore functions
source("../funcs.R")

## data
inpori <- pol$albacore
## increase speed of testmore
inpori$dteuler <- 1/4

## Default assessment
inp <- check.inp(inpori)
rep <- fit.spict(inp)

header("1: Changing maninterval in 3 ways", append = FALSE)
###################################################################

## Before check.inp
inp1 <- inpori
inp1$maninterval <- c(1991,1992)
inp1 <- check.inp(inp1)
rep1 <- fit.spict(inp1)

## After check.inp
inp2 <- inp
inp2$maninterval <- c(1991,1992)
inp2 <- check.inp(inp2)
rep2 <- fit.spict(inp2)

## With argument of new function on fitted object
rep3 <- check.man.time(rep, maninterval = c(1991,1992))

## All approaches should give same parameter estimates
out(all.equal(sumspict.parest(rep1),
          sumspict.parest(rep2), tolerance = 0.01))

out(all.equal(sumspict.parest(rep1),
          sumspict.parest(rep3), tolerance = 0.01))

b=get.par("logBBmsy",rep)[,2]
b1=get.par("logBBmsy",rep1)[,2]
b2=get.par("logBBmsy",rep2)[,2]
b3=get.par("logBBmsy",rep3)[,2]

## Last 3 reps are one year longer than first rep
out(all.equal(length(b) + 1/rep$inp$dteuler, length(b1), tolerance = 0.01))
out(all.equal(length(b) + 1/rep$inp$dteuler, length(b2), tolerance = 0.01))
out(all.equal(length(b) + 1/rep$inp$dteuler, length(b3), tolerance = 0.01))

## Last 3 reps should give same predictions
out(all.equal(sumspict.parest(rep1),
              sumspict.parest(rep2), tolerance = 0.01))

out(all.equal(sumspict.parest(rep1),
          sumspict.parest(rep3), tolerance = 0.01))



header("2: Changing maneval in 3 ways", append = TRUE)
###################################################################

## Before check.inp
inp1 <- inpori
inp1$maneval <- 1992
inp1 <- check.inp(inp1)
rep1 <- fit.spict(inp1)

## After check.inp
inp2 <- inp
inp2$maneval <- 1992
inp2 <- check.inp(inp2)
rep2 <- fit.spict(inp2)

## With argument of new function on fitted object
rep3 <- check.man.time(rep, maneval = 1992)

## All approaches should give same parameter estimates
out(all.equal(sumspict.parest(rep1),
          sumspict.parest(rep2), tolerance = 0.01))

out(all.equal(sumspict.parest(rep1),
          sumspict.parest(rep3), tolerance = 0.01))

b=get.par("logBBmsy",rep)[,2]
b1=get.par("logBBmsy",rep1)[,2]
b2=get.par("logBBmsy",rep2)[,2]
b3=get.par("logBBmsy",rep3)[,2]

## Last 3 reps are one year longer than first rep
out(all.equal(length(b) + 1/rep$inp$dteuler, length(b1), tolerance = 0.01))
out(all.equal(length(b) + 1/rep$inp$dteuler, length(b2), tolerance = 0.01))
out(all.equal(length(b) + 1/rep$inp$dteuler, length(b3), tolerance = 0.01))

b1=get.par("logBp",rep1)[,2]
b2=get.par("logBp",rep2)[,2]
b3=get.par("logBp",rep3)[,2]

## Bp of last 3 reps shoule be the same
out(all.equal(b1,b2, tolerance = 0.01))
out(all.equal(b2,b3, tolerance = 0.01))

## Last 3 reps should give same predictions
out(all.equal(sumspict.parest(rep1),
              sumspict.parest(rep2), tolerance = 0.01))
out(all.equal(sumspict.parest(rep1),
              sumspict.parest(rep3), tolerance = 0.01))
