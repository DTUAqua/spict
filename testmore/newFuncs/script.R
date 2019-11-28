## Testing some new functions and modifications
## T.K. Mildenberger <t.k.mildenberger@gmail.com>
## 26/11/2019

## set seed
set.seed(123)

## load spict
suppressMessages(library(spict))

## error settings
options(error = expression(NULL))


## 1: annual data
##--------------------
writeLines("Tests with annual data")
inp <- pol$lobster
suppressWarnings(fit1 <- fit.spict(inp))
fit1$opt$convergence

## 1.1: calculate Bmsy/K ratio
writeLines("Test 1.1: ")
round(calc.bmsyk(fit1),3)

## 1.2: calculate order of magnitude
writeLines("Test 1.2: ")
round(calc.om(fit1),3)




## 2: seasonal data
##--------------------
writeLines("\n\nTests with seasonal data")
nt <- 50
inp <- list(nseasons=4, splineorder=3)
inp$timeC <- seq(0, nt-1/inp$nseasons, by=1/inp$nseasons)
inp$timeI <- seq(0.1, nt-1/inp$nseasons, by=0.5)
inp$ini <- list(logK=log(1000), logm=log(800), logq=log(1), logn=log(2),
                logbkfrac=log(0.9), logsdf=log(0.3), logF0=log(0.8),
                logphi=log(c(0.3, 0.5, 1.8)))
inpsim <- sim.spict(inp)
suppressWarnings(fit2 <- fit.spict(inpsim))
fit2$opt$convergence


## 2.1: calculate Bmsy/K ratio
writeLines("Test 2.1: ")
round(calc.bmsyk(fit2),3)

## 2.2: calculate order of magnitude
writeLines("Test 2.2: ")
round(calc.om(fit2),3)




## 3: mixed data
##--------------------
writeLines("\n\nTests with mixed data")
nt <- 50
inp <- list(nseasons=4, splineorder=3)
inp$timeC <- c(seq(0, nt/2-1, by=1),seq(nt/2, nt-1/inp$nseasons, by=1/inp$nseasons))
inp$timeI <- seq(0.1, nt-1/inp$nseasons, by=0.5)
inp$ini <- list(logK=log(1000), logm=log(800), logq=log(1), logn=log(2),
                logbkfrac=log(0.9), logsdf=log(0.3), logF0=log(0.8),
                logphi=log(c(0.3, 0.5, 1.8)))
inpsim <- sim.spict(inp)


inp$timeI <- seq(0.6, 29.6, by=1)
inp$ini <- list(logK=log(100), logm=log(60), logq=log(1),
                logbkfrac=log(1), logsdf=log(0.3), logF0=log(0.5),
                logphi=log(c(0.05, 0.1, 1.8)))
inpsim <- sim.spict(inp)
suppressWarnings(fit3 <- fit.spict(inpsim))
fit3$opt$convergence

## 3.1: calculate Bmsy/K ratio
writeLines("Test 3.1: ")
round(calc.bmsyk(fit3),3)

## 3.2: calculate order of magnitude
writeLines("Test 3.2: ")
round(calc.om(fit3),3)




## 4: incorrect input
##--------------------
writeLines("\n\nTests with incorrect input")

## 4.1: wrong input
## ---------------
## 4.1.1: calculate Bmsy/K ratio
writeLines("Test 4.1.1: ")
calc.bmsyk(inp)

## 4.1.2: calculate order of magnitude
writeLines("Test 4.1.2: ")
calc.om(inp)



## 4.2: with non converged model
## ---------------
inp <- pol$lobster
inp$obsC <- rnorm(46, 2000, 3)
suppressWarnings(fit4 <- fit.spict(inp))
fit4$opt$convergence

## 4.2.1: calculate Bmsy/K ratio
writeLines("Test 4.2.1: ")
round(calc.bmsyk(fit4),3)

## 4.2.2: calculate order of magnitude
writeLines("Test 4.2.2: ")
round(calc.om(fit4),3)
