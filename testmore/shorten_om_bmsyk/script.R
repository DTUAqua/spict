## Testing some new functions and modifications
## T.K. Mildenberger <t.k.mildenberger@gmail.com>
## 26/11/2019

## set seed
set.seed(123)

## load spict
library(spict)

## Test functions
out <- function(..., sep = "", append = TRUE)
  cat(..., "\n", file = "res.out", append = append, sep = sep)
get_nchar <- function(...)
  nchar(paste(as.character(unlist(list(...))), collapse = " "))
header <- function(..., append = TRUE)
  out("\n", ..., "\n", rep("=", get_nchar(...)), "\n", append = append)
test_this <- function(title, expression) {
  out(title)
  tryCatch(out(capture.output(eval(expression)), sep = "\n"),
           error = function(e) out("Error:", e$message))
}


header("1: annual data", append = FALSE)

inp <- pol$lobster
fit1 <- fit.spict(inp)

out(fit1$opt$convergence)

test_this("Test 1.1: calculate Bmsy/K ratio", {
  round(calc.bmsyk(fit1), 3)
})

test_this("Test 1.2: calculate order of magnitude", {
  round(calc.om(fit1), 3)
})

## ## 1.3: shorten input
## writeLines("Test 1.3: shorten.inp")
## inp_full <- check.inp(inp)

## ## 1.3.1: no mintime and maxtime
## inp_full2 <- shorten.inp(inp)
## all.equal(inp_full, inp_full2)

## ## 1.3.2: only mintime
## inp_mintime <- shorten.inp(inp, mintime=1980)
## all.equal(capture.output(inp_full), capture.output(inp_mintime))

## ## 1.3.3: only maxtime
## inp_maxtime <- shorten.inp(inp, maxtime=1982)
## all.equal(inp_full, inp_full2)

## ## 1.3.4: mintime and maxtime
## inp_both <- shorten.inp(inp, mintime=1975.001, maxtime=1981.999)
## all.equal(inp_full, inp_full2)


header("2: seasonal data")

nt <- 50
inp <- list(nseasons = 4, splineorder = 3)
inp$timeC <- seq(0, nt - 1 / inp$nseasons, by = 1 / inp$nseasons)
inp$timeI <- seq(0.1, nt - 1 / inp$nseasons, by = 0.5)
inp$ini <- list(logK = log(1000), logm=log(800), logq = log(1), logn = log(2),
                logbkfrac = log(0.9), logsdf = log(0.3), logF0 = log(0.8),
                logphi = log(c(0.3, 0.5, 1.8)))
inpsim <- sim.spict(inp)
fit2 <- fit.spict(inpsim)

out(fit2$opt$convergence)

test_this("2.1: calculate Bmsy/K ratio", {
  round(calc.bmsyk(fit2), 3)
})


test_this("2.2: calculate order of magnitude", {
  round(calc.om(fit2), 3)
})



header("3: Tests with mixed data")

nt <- 50
inp <- list(nseasons=4, splineorder=3)
inp$timeC <- c(seq(0, nt/2-1, by=1),seq(nt/2, nt-1/inp$nseasons, by=1/inp$nseasons))
inp$timeI <- seq(0.1, nt-1/inp$nseasons, by=0.5)
inp$ini <- list(logK=log(1000), logm=log(800), logq=log(1), logn=log(2),
                logbkfrac=log(0.9), logsdf=log(0.3), logF0=log(0.8),
                logphi=log(c(0.3, 0.5, 1.8)))
inpsim <- sim.spict(inp)


inp$timeI <- seq(0.6, 29.6, by = 1)
inp$ini <- list(logK = log(100), logm = log(60), logq = log(1),
                logbkfrac = log(1), logsdf = log(0.3), logF0 = log(0.5),
                logphi = log(c(0.05, 0.1, 1.8)))
inpsim <- sim.spict(inp)
fit3 <- fit.spict(inpsim)

out(fit3$opt$convergence)

test_this("3.1: calculate Bmsy/K ratio", {
  round(calc.bmsyk(fit3), 3)
})

test_this("3.2: calculate order of magnitude", {
  round(calc.om(fit3), 3)
})



header("4: incorrect input")

out("4.1: wrong input")

test_this("4.1.1: Bmsy/K ratio", {
  calc.bmsyk(inp)
})


test_this("4.1.2: calculate order of magnitude", {
  calc.om(inp)
})


out("4.2: with non converged model")

inp <- pol$lobster
inp$obsC <- rnorm(46, 2000, 3)
fit4 <- fit.spict(inp)

out(fit4$opt$convergence)

test_this("4.2.1: calculate Bmsy/K ratio", {
  round(calc.bmsyk(fit4), 3)
})
test_this("4.2.2: calculate order of magnitude", {
  round(calc.om(fit4), 3)
})
