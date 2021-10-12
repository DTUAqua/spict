## Testing new function to retape fitted spict objects
## T.K. Mildenberger <t.k.mildenberger@gmail.com>
## 12/12/2019

dec.val <- 2
dteuler <- 1/4

## set seed
set.seed(123)

## load spict
require(spict)

## load testmore functions
source("../funcs.R")

## load data
inp <- pol$albacore
inp$dteuler <- dteuler
inp <- check.inp(inp)

## fit spict
rep <- fit.spict(inp)


header("1: Changing prediction horizon", append = FALSE)
#####################################################


rept <- spict:::retape.spict(rep, inp)

test_this("1.1: providing input list without change",{
    round(sumspict.predictions(rept),dec.val)
})

out(all(rep$inp$maninterval == rept$inp$maninterval))

out(all.equal(sumspict.predictions(rep),
              sumspict.predictions(rept),tolerance = 0.01))

inpt <- check.man.time(inp, maninterval = c(1991,1992), verbose = FALSE)
## with retape
rept <- spict:::retape.spict(rep, inpt)
## alternative way
inpx <- inp
inpx$maninterval <- c(1991,1992)
suppressWarnings(repx <- fit.spict(inpx))
test_this("1.2: extending prediction horizon with maninterval",{
    round(sumspict.predictions(rept),dec.val)
})

out(all(rep$inp$maninterval + 1 == rept$inp$maninterval))

out(all.equal(sumspict.predictions(repx),
              sumspict.predictions(rept),tolerance = 0.01))


inpt <- check.man.time(inp, maneval = 1996, verbose = FALSE)
rept <- spict:::retape.spict(rep, inpt)
## alternative way
inpx <- inp
inpx$maneval <- 1996
suppressWarnings(repx <- fit.spict(inpx))
test_this("1.3: extending prediction horizon with maneval",{
    round(sumspict.predictions(rept),dec.val)
})

out(all(rep$inp$maneval + 5 == rept$inp$maneval))

out(all.equal(sumspict.predictions(repx),
              sumspict.predictions(rept),tolerance = 0.01))

inpt$maneval <- 1989
inpt <- check.man.time(inpt, verbose = FALSE)
rept <- spict:::retape.spict(rep, inpt)
## alternative way
inpx <- inp
inpx$maneval <- 1989
repx <- fit.spict(inpx)
test_this("1.3: Shortening prediction horizon",{
    round(sumspict.predictions(rept),dec.val)
})

out(all(rep$inp$maneval - 2 == rept$inp$maneval))

out(all.equal(sumspict.predictions(repx),
              sumspict.predictions(rept),tolerance = 0.01))

inpt <- inp
inpt$maninterval <- c(1990,1990.5)
inpt <- check.man.time(inpt, verbose = FALSE)
rept <- spict:::retape.spict(rep, inpt)
## alternative way
inpx <- inp
inpx$maninterval <- c(1990,1990.5)
repx <- fit.spict(inpx)
test_this("1.4: Half year management interval",{
    round(sumspict.predictions(rept),dec.val)
})

out(get.par("logCp",rept,exp=TRUE)[2] < get.par("logCp",rep,exp=TRUE)[2])

out(all.equal(sumspict.predictions(repx),
              sumspict.predictions(rept),tolerance = 0.01))



header("2: Changing F vector", append = TRUE)
#####################################################

## Testing ffac
inpt <- make.ffacvec(inp, ffac = 0.4)

## old code
plt <- inpt$parlist
datint <- make.datin(inpt)
objt <- make.obj(datint, plt, inpt, phase=1)
objt$retape()
objt$fn(rep$opt$par)
rep1 <- TMB::sdreport(objt)
rep1$inp <- inpt
rep1$obj <- objt
rep1$opt <- rep$opt
class(rep1) <- "spictcls"

## new function
rep2 <- spict:::retape.spict(rep, inpt)

test_this("2.1: ffac == 0.4",{
    all.equal(sumspict.predictions(rep1),
              sumspict.predictions(rep2), tolerance = 0.01)
})


## Testing fcon
inpt <- make.fconvec(inp, fcon = 1.3)

## old code
plt <- inpt$parlist
datint <- make.datin(inpt)
objt <- make.obj(datint, plt, inpt, phase=1)
objt$retape()
objt$fn(rep$opt$par)
rep1 <- TMB::sdreport(objt)
rep1$inp <- inpt
rep1$obj <- objt
rep1$opt <- rep$opt
class(rep1) <- "spictcls"

## new function
rep2 <- spict:::retape.spict(rep, inpt)

test_this("2.2: fcon == 1.3",{
    all.equal(sumspict.predictions(rep1),
              sumspict.predictions(rep2), tolerance = 0.01)
})




header("3: Testing warnings & errors", append = TRUE)
#####################################################

test_this("3.1: wrong input to make.ffacvec",{
    inpt <- make.ffacvec(pol$albacore, ffac=1.2)
})

test_this("3.2: wrong input to make.fconvec",{
    inpt <- make.ffacvec(pol$albacore, ffac=1.2)
})

test_this("3.2: incorrect rep",{
    rept <- spict:::retape.spict(inpt, inpt)
})

test_this("3.3: non meaningful inp",{
    rept <- spict:::retape.spict(rep, list())
})
