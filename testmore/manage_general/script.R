## Testing new management functionality
## T.K. Mildenberger <t.k.mildenberger@gmail.com>
## 16/02/2020

## flag for printing to consol (dev mode)
cnsl <- FALSE

## set seed
set.seed(123)

## load spict
require(spict)

## load testmore functions
source("../funcs.R")

## all tests based on pol$albacore and seasonal data
## annual
inpA <- pol$albacore
inpa <- check.inp(inpA)
## seasonal
nt <- 40
inp <- list(nseasons = 4, splineorder = 3)
inp$timeC <- seq(0, nt - 1 / inp$nseasons, by = 1 / inp$nseasons)
inp$timeI <- seq(0.1, nt - 1 / inp$nseasons, by = 0.5)
inp$ini <- list(logK = log(1000), logm=log(800), logq = log(1), logn = log(2),
                logbkfrac = log(0.9), logsdf = log(0.3), logF0 = log(0.8),
                logphi = log(c(0.3, 0.5, 1.8)))
inpS <- sim.spict(inp)
inps <- check.inp(inpS)

header("1: check.inp with new functionality", append = FALSE)
######################################################

out("1.1: Same defaults of management variables as before (v1.2.8)")
## ---------------------------------------------------
out("1.1.1: Annual data")

## defaults - v1.2.8 - annual
timepredi.v128 <- 1990
timepredc.v128 <- 1990
dtpredc.v128 <- 1
manstart.v128 <- 1990

out("Most management variables have the same default as spict v1.2.8:")
out(inpa$timepredc == timepredc.v128)
out(inpa$dtpredc == dtpredc.v128)
out(inpa$manstart == manstart.v128)
out(inpa$maninterval[1] == manstart.v128)

out("New default for timepredi (at end of maninterval) compared to v1.2.8:")
out(inpa$timepredi == timepredi.v128 + 1)

out("Consistency between new variables:")
out(inpa$maneval == inpa$timepredi)
out(all.equal(inpa$maninterval,c(inpa$timepredc,inpa$timepredc+inpa$dtpredc)))
out(inpa$maninterval[1] == inpa$manstart)


out("1.1.2: Seasonal data")

## defaults - v1.2.8 - seasonal
timepredi.v128 <- 40
timepredc.v128 <- 40
dtpredc.v128 <- 0.25
manstart.v128 <- 40

out("Most management variables have the same default as spict v1.2.8:")
out(inps$timepredc == timepredc.v128)
out(inps$manstart == manstart.v128)
out(inps$maninterval[1] == manstart.v128)

out("New default management interval length for seasonal data (1 year instead of min(dt)) compared to v1.2.8:")
out(inps$dtpredc == dtpredc.v128 * 4)

out("New default for timepredi (at end of maninterval) compared to v1.2.8:")
out(inps$timepredi == timepredi.v128 + 1)

out("Consistency between new variables:")
out(inps$maneval == inps$timepredi)
out(all(inps$maninterval == c(inps$timepredc,inps$timepredc+inps$dtpredc)))
out(inps$maninterval[1] == inps$manstart)


out("1.2: Changing maninterval - changes old variables & reverse")
## ---------------------------------------------------
out("1.2.1: Annual data")
out("1.2.1.1: From new to old")

inp <- inpA
inp$maninterval <- c(1991,1992)
inp$maneval <- 1994
inp <- check.inp(inp)

out(inp$maneval == inp$timepredi)
out(all.equal(inp$maninterval,c(inp$timepredc,inp$timepredc+inp$dtpredc)))
out(inp$maninterval[1] == inp$manstart)

out("1.2.1.2: From old to new")

inp <- inpA
inp$timepredc <- 1992
inp$dtpredc <- 2
inp$timepredi <- 1994
inp$manstart <- 1992
inp <- check.inp(inp) ## all three variables need to be provided

out(inp$maneval == inp$timepredi)
out(all.equal(inp$maninterval,c(inp$timepredc,inp$timepredc+inp$dtpredc)))
out(inp$maninterval[1] == inp$manstart)

out("1.2.2: Seasonal data")
out("1.2.1.1: From new to old")

inp <- inpS
inp$maninterval <- c(41,42)
inp$maneval <- 44
inp <- check.inp(inp)

out(inp$maneval == inp$timepredi)
out(all.equal(inp$maninterval,c(inp$timepredc,inp$timepredc+inp$dtpredc)))
out(inp$maninterval[1] == inp$manstart)

out("1.2.1.2: From old to new")

inp <- inpS
inp$timepredc <- 42
inp$dtpredc <- 2
inp$timepredi <- 44
inp$manstart <- 42
inp <- check.inp(inp)

out(inp$maneval == inp$timepredi)
out(all.equal(inp$maninterval,c(inp$timepredc,inp$timepredc+inp$dtpredc)))
out(inp$maninterval[1] == inp$manstart)

inp <- pol$albacore
inp <- check.inp(inp)
inp$maneval <- 1994
suppressWarnings(inp <- check.inp(inp))
fit <- fit.spict(inp)

test_this("1.3: Adjust all time-dependent variables in check.inp",{
    fit$opt$convergence
})


out("1.4: Management variables with shorten.inp")
inp <- check.inp(inpS)
inp2 <- shorten.inp(inp, mintime = 20)
out(inp$maninterval[1] == inp2$maninterval[1])
out(inp$maninterval[2] == inp2$maninterval[2])

inp2 <- shorten.inp(inp, maxtime = 38)
out(inp$maninterval[1] == inp2$maninterval[1] + 2)
out(inp$maninterval[2] == inp2$maninterval[2] + 2)
out(inp$maneval == inp2$maneval + 2)


header("2: check.man.time")
######################################################

inpb <- inpa
inpb$maninterval <- inpb$maninterval + 1
inpb <- check.man.time(inpb)
out("on inp list")
out(length(inpa$seasonindex2) + 1/inpa$dteuler == length(inpb$seasonindex2))

resa <- fit.spict(inpa)
resb <- resa
resb$inp$maninterval <- resb$inp$maninterval + 1
resb <- check.man.time(resb)
out("on rep list")
out(length(resa$inp$logmcovariatein) + 1/resa$inp$dteuler == length(resb$inp$logmcovariatein))


header("3: timeline")
######################################################

test_this("3.1:", man.timeline(inpa))
test_this("3.2:", man.timeline(inpb))
test_this("3.3:", man.timeline(resa))
test_this("3.4:", man.timeline(resb))


header("4: manage / summary / scenarios / intermediate periods")
######################################################
inp <- inpA
inp$maninterval <- c(1990,1991)
inp <- check.inp(inp)

fit <- fit.spict(inp)

## no intermediate period
mana <- manage(fit, scenarios=c(1,2))

mana1 <- add.man.scenario(mana, maninterval = c(1991.5,1992.5), maneval = 1993)
test_this("4.1: Management periods do not match between scenarios",{
    lapply(man.tac(mana1),round,2)
})

out(man.timeline(mana1))

mana2 <- add.man.scenario(mana, intermediatePeriodCatch = 12)
test_this("4.2: Scenario with continue F vs continue C in intermediate period",{
    lapply(man.tac(mana2),round,2)
})

out(man.timeline(mana2))

## intermediate period
mana3 <- manage(fit, c(1,4), maninterval = c(1993,1995))

test_this("4.3: New summary function",{
    sumspict.manage(mana3)
})

test_this("4.4: Old summary function still works:",{
    mansummary(mana3)
})

out(sumspict.manage(tmp <- man.select(mana3,1))$est)

out(sumspict.manage(tmp <- man.select(mana3,"noF"))$est)

header("5: estimating TAC directly/indirectly")
######################################################

test_this("5.1: TAC from all management scenarios in rep",{
    lapply(man.tac(mana2),round,2)
})

test_this("5.2: TAC for one scenario",{
    round(get.TAC(fit),2)
})


header("6: Provoke warnings & errors & challenge functionality")
######################################################

out("6.1: Message when manstart and maninterval used and not equal")
inp <- inpA
inp$maninterval <- c(1992,1993)
inp$manstart <- 1991
inp <- check.inp(inp)

out("6.2: Message when timepredc and maninterval used and not equal")
inp <- inpA
inp$maninterval <- c(1992,1993)
inp$timepredc <- 1991
inp <- check.inp(inp)

out("6.3: Message when timepredi and maneval used and not equal")
inp <- inpA
inp$maneval <- 1992
inp$timepredi <- 1991
inp <- check.inp(inp)

out("6.4: Maninterval does not match dteuler time steps")
inp <- inpA
inp$maninterval <- c(1991.573838,1992)
inp <- check.inp(inp)

out("6.5: Maninterval smaller than dteuler")
inp <- inpA
inp$maninterval <- c(1991,1991.004)
suppressWarnings(inp <- check.inp(inp))

out("6.6: Warning message when manstart larger than timepredc")
inp <- inpA
inp$timepredc <- 1992
inp$dtpredc <- 2
inp$timepredi <- 1994
inp$manstart <- 1993
suppressWarnings(inp <- check.inp(inp))

out("6.7: manstart smaller than timepredc")
inp1 <- inpA
inp1$timepredi <- 1993
inp1$timepredc <- 1992
inp1$dtpredc <- 1
inp1$manstart <- 1991
inp1$ffac <- 0.2
inp1 <- check.inp(inp1)
set.seed(124)
fit1 <- fit.spict(inp1)

## with Cp right after ffac applied
inp2 <- inpA
inp2$timepredi <- 1993
inp2$timepredc <- 1991
inp2$dtpredc <- 1
inp2$manstart <- 1991
inp2$ffac <- 0.2
inp2 <- check.inp(inp2)
set.seed(124)
fit2 <- fit.spict(inp2)
## Cp differ:
out(fit1$Cp != fit2$Cp)
## But values the same just timing different:
out(round(fit2$Cp,3) == round(get.par("logCpred",fit1,exp=TRUE)[which(inp1$timeCpred == inp2$timepredc),2],3))
## => Thus old functionality and independence between manstart and timepredc is maintained.

## missing one of the three required variables => error
test_this("6.8: Not all variables provided", {
    inp <- inpA
    inp$timepredc <- 1992
    inp$manstart <- 1992
    inp <- check.inp(inp) ## all three variables need to be provided
    inp <- inpA
    inp$timepredc <- 1992
    inp$dtpredc <- 2
    inp <- check.inp(inp) ## all three variables need to be provided
    inp <- inpA
    inp$dtpredc <- 2
    inp$manstart <- 1992
    inp <- check.inp(inp) ## all three variables need to be provided
})

test_this("6.9: wrong input to functions, creating warnings...",{
    man.select(inp)
    manage(inp)
    add.man.scenario(inp)
    get.TAC(inp)
    man.tac(inp)
    spict:::check.man(inp)
    man.select(fit)
    man.tac(fit)
    spict:::check.man(fit)
})


header("7: plotting")
######################################################

if(FALSE){
    ## annual
    ## ----------------
    inp <- pol$albacore
    fit <- fit.spict(inp)

    ## man period = 1, int period = 0
    mana1 <- manage(fit, c(1,2,4), maninterval = c(1990,1991))

    ## man period = 1, int period = 1
    mana2 <- manage(fit, c(1,2,4), maninterval = c(1991,1992))

    ## man period = 2, int period = 1
    mana3 <- manage(fit, c(1,2,4), maninterval = c(1991,1993))

    ## man period = 0.5, int period = 1
    mana4 <- manage(fit, c(1,2,4), maninterval = c(1991,1991.5))

    ## man period = 1, int period = 0.5
    mana5 <- manage(fit, c(1,2,4), maninterval = c(1990.5,1991.5))

    ## man period = 1.5, int period = 0.5
    mana6 <- manage(fit, c(1,2,4), maninterval = c(1990.5,1992))

    ## man period = 0.5, int period = 0.5
    mana7 <- manage(fit, c(1,2,4), maninterval = c(1990.5,1991))

    ## man period = 2.5, int period = 1.5, intercatch
    mana8 <- manage(fit, c(1,2,4), maninterval = c(1991.5,1994), intermediatePeriodCatch = 29)


    ## check plots
    pdf("manscenariosCatchPlot.pdf",width = 10,height=12)
    opar <- par(mfrow=c(4,2))
    for(i in 1:8){
        tmp <- get(paste0("mana",i))
        plotspict.catch(tmp)
    }
    par(opar)
    dev.off()

    ## check plots
    pdf("manscenariosBBmsyPlot.pdf",width = 10,height=12)
    opar <- par(mfrow=c(4,2))
    for(i in 1:8){
        tmp <- get(paste0("mana",i))
        plotspict.bbmsy(tmp)
    }
    par(opar)
    dev.off()

    ## check plots
    pdf("manscenariosFFmsyPlot.pdf",width = 10,height=12)
    opar <- par(mfrow=c(4,2))
    for(i in 1:8){
        tmp <- get(paste0("mana",i))
        plotspict.ffmsy(tmp)
    }
    par(opar)
    dev.off()

    ## seasonal
    ## ----------------
    set.seed(455)
    nt <- 25
    inp <- list(nseasons = 4, splineorder = 3)
    inp$timeC <- seq(0, nt - 1 / inp$nseasons, by = 1 / inp$nseasons)
    inp$timeI <- seq(0.1, nt - 1 / inp$nseasons, by = 0.5)
    inp$ini <- list(logK = log(1000), logm=log(800), logq = log(1), logn = log(2),
                    logbkfrac = log(0.9), logsdf = log(0.3), logF0 = log(0.8),
                    logphi = log(c(0.3, 0.5, 1.8)))
    inpsim <- sim.spict(inp)
    fits <- fit.spict(inpsim)

    ## man period = 1, int period = 0
    mana1s <- manage(fits, c(1,2,4), maninterval = c(25,26))

    ## man period = 1, int period = 1
    mana2s <- manage(fits, c(1,2,4), maninterval = c(26,27))

    ## man period = 2, int period = 1
    mana3s <- manage(fits, c(1,2,4), maninterval = c(26,28))

    ## man period = 0.5, int period = 1
    mana4s <- manage(fits, c(1,2,4), maninterval = c(26,26.5))

    ## man period = 1, int period = 0.5
    mana5s <- manage(fits, c(1,2,4), maninterval = c(25.5,26.5))

    ## man period = 1.5, int period = 0.5
    mana6s <- manage(fits, c(1,2,4), maninterval = c(25.5,27))

    ## man period = 0.5, int period = 0.5
    mana7s <- manage(fits, c(1,2,4), maninterval = c(25.5,26))

    ## man period = 2.5, int period = 1.5, intercatch
    mana8s <- manage(fits, c(1,2,4), maninterval = c(26.5,29), intermediatePeriodCatch = 29)


    pdf("manscenariosCatchPlot_sea.pdf",width = 10,height=12)
    ## check plots
    opar <- par(mfrow=c(4,2))
    for(i in 1:8){
        tmp <- get(paste0("mana",i,"s"))
        plotspict.catch(tmp)
    }
    par(opar)
    dev.off()

    pdf("manscenariosBBmsyPlot_sea.pdf",width = 10,height=12)
    ## check plots
    opar <- par(mfrow=c(4,2))
    for(i in 1:8){
        tmp <- get(paste0("mana",i,"s"))
        plotspict.bbmsy(tmp)
    }
    par(opar)
    dev.off()

    pdf("manscenariosFFmsyPlot_sea.pdf",width = 10,height=12)
    ## check plots
    opar <- par(mfrow=c(4,2))
    for(i in 1:8){
        tmp <- get(paste0("mana",i,"s"))
        plotspict.ffmsy(tmp)
    }
    par(opar)
    dev.off()

}
