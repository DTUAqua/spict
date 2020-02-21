## Testing new management functionality
## T.K. Mildenberger <t.k.mildenberger@gmail.com>
## 16/02/2020

## flag for printing to consol (dev mode)
cnsl = FALSE

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

header("1: check.inp with new functionality", append = FALSE, cnsl=cnsl)
######################################################

out("1.1: Same defaults of management variables as before (v1.2.8)", cnsl=cnsl)
## ---------------------------------------------------
out("1.1.1: Annual data", cnsl=cnsl)

## defaults - v1.2.8 - annual
timepredi.v128 <- 1990
timepredc.v128 <- 1990
dtpredc.v128 <- 1
manstart.v128 <- 1990

out("Most management variables have the same default as spict v1.2.8:")
out(inpa$timepredc == timepredc.v128, cnsl=cnsl)
out(inpa$dtpredc == dtpredc.v128, cnsl=cnsl)
out(inpa$manstart == manstart.v128, cnsl=cnsl)
out(inpa$maninterval[1] == manstart.v128, cnsl=cnsl)

out("New default for timepredi (at end of maninterval) compared to v1.2.8:")
out(inpa$timepredi == timepredi.v128 + 1, cnsl=cnsl)

out("Consistency between new variables:")
out(inpa$maneval == inpa$timepredi, cnsl=cnsl)
out(all.equal(inpa$maninterval,c(inpa$timepredc,inpa$timepredc+inpa$dtpredc)), cnsl=cnsl)
out(inpa$maninterval[1] == inpa$manstart, cnsl=cnsl)


out("1.1.2: Seasonal data", cnsl=cnsl)
## defaults - v1.2.8 - seasonal
timepredi.v128 <- 40
timepredc.v128 <- 40
dtpredc.v128 <- 0.25
manstart.v128 <- 40

out("Most management variables have the same default as spict v1.2.8:")
out(inps$timepredc == timepredc.v128, cnsl=cnsl)
out(inps$manstart == manstart.v128, cnsl=cnsl)
out(inps$maninterval[1] == manstart.v128, cnsl=cnsl)

out("New default management interval length for seasonal data (1 year instead of min(dt)) compared to v1.2.8:")
out(inps$dtpredc == dtpredc.v128 * 4, cnsl=cnsl)

out("New default for timepredi (at end of maninterval) compared to v1.2.8:")
out(inps$timepredi == timepredi.v128 + 1, cnsl=cnsl)

out("Consistency between new variables:")
out(inps$maneval == inps$timepredi, cnsl=cnsl)
out(all(inps$maninterval == c(inps$timepredc,inps$timepredc+inps$dtpredc)), cnsl=cnsl)
out(inps$maninterval[1] == inps$manstart, cnsl=cnsl)


out("1.2: Changing maninterval - changes old variables & reverse", cnsl=cnsl)
## ---------------------------------------------------
out("1.2.1: Annual data", cnsl=cnsl)
out("1.2.1.1: From new to old", cnsl=cnsl)
inp <- inpA
inp$maninterval <- c(1991,1992)
inp$maneval <- 1994
inp <- check.inp(inp)

out(inp$maneval == inp$timepredi, cnsl=cnsl)
out(all.equal(inp$maninterval,c(inp$timepredc,inp$timepredc+inp$dtpredc)), cnsl=cnsl)
out(inp$maninterval[1] == inp$manstart, cnsl=cnsl)

out("1.2.1.2: From old to new", cnsl=cnsl)
inp <- inpA
inp$timepredc <- 1992
inp$dtpredc <- 2
inp$timepredi <- 1994
inp$manstart <- 1992
inp <- check.inp(inp) ## all three variables need to be provided

out(inp$maneval == inp$timepredi, cnsl=cnsl)
out(all.equal(inp$maninterval,c(inp$timepredc,inp$timepredc+inp$dtpredc)), cnsl=cnsl)
out(inp$maninterval[1] == inp$manstart, cnsl=cnsl)

out("1.2.2: Seasonal data", cnsl=cnsl)
out("1.2.1.1: From new to old", cnsl=cnsl)
inp <- inpS
inp$maninterval <- c(41,42)
inp$maneval <- 44
inp <- check.inp(inp)

out(inp$maneval == inp$timepredi, cnsl=cnsl)
out(all.equal(inp$maninterval,c(inp$timepredc,inp$timepredc+inp$dtpredc)), cnsl=cnsl)
out(inp$maninterval[1] == inp$manstart, cnsl=cnsl)

out("1.2.1.2: From old to new", cnsl=cnsl)
inp <- inpS
inp$timepredc <- 42
inp$dtpredc <- 2
inp$timepredi <- 44
inp$manstart <- 42
inp <- check.inp(inp)

out(inp$maneval == inp$timepredi, cnsl=cnsl)
out(all.equal(inp$maninterval,c(inp$timepredc,inp$timepredc+inp$dtpredc)), cnsl=cnsl)
out(inp$maninterval[1] == inp$manstart, cnsl=cnsl)

test_this("1.3: Adjust all time-dependent variables in check.inp",{
    inp <- pol$albacore
    inp <- check.inp(inp)
    inp$maneval <- 1994
    inp <- check.inp(inp)
    fit <- fit.spict(inp)
    fit$opt$convergence
},cnsl=cnsl)


out("1.4: Management variables with shorten.inp",cnsl=cnsl)
inp <- check.inp(inpS)
inp2 <- shorten.inp(inp, mintime = 20)
out(inp$maninterval[1] == inp2$maninterval[1],cnsl=cnsl)
out(inp$maninterval[2] == inp2$maninterval[2],cnsl=cnsl)

inp2 <- shorten.inp(inp, maxtime = 38)
out(inp$maninterval[1] == inp2$maninterval[1] + 2,cnsl=cnsl)
out(inp$maninterval[2] == inp2$maninterval[2] + 2,cnsl=cnsl)
out(inp$maneval == inp2$maneval + 2,cnsl=cnsl)


header("2: check.man.time", cnsl=cnsl)
######################################################

inpb <- inpa
inpb$maninterval <- inpb$maninterval + 1
inpb <- check.man.time(inpb)
out("on inp list", cnsl=cnsl)
out(length(inpa$seasonindex2) + 1/inpa$dteuler == length(inpb$seasonindex2), cnsl=cnsl)

resa <- fit.spict(inpa)
resb <- resa
resb$inp$maninterval <- resb$inp$maninterval + 1
resb <- check.man.time(resb)
out("on rep list", cnsl=cnsl)
out(length(resa$inp$logmcovariatein) + 1/resa$inp$dteuler == length(resb$inp$logmcovariatein), cnsl=cnsl)


header("3: timeline")
######################################################

test_this("3.1:", man.timeline(inpa),cnsl=cnsl)
test_this("3.2:", man.timeline(inpb),cnsl=cnsl)
test_this("3.3:", man.timeline(resa),cnsl=cnsl)
test_this("3.4:", man.timeline(resb),cnsl=cnsl)


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
},cnsl=cnsl)

out(man.timeline(mana1),cnsl=cnsl)

mana2 <- add.man.scenario(mana, intermediatePeriodCatch = 12)
test_this("4.2: Scenario with continue F vs continue C in intermediate period",{
    lapply(man.tac(mana2),round,2)
},cnsl=cnsl)

out(man.timeline(mana2),cnsl=cnsl)

## intermediate period
mana3 <- manage(fit, c(1,4), maninterval = c(1993,1995))

test_this("4.3: New summary function",{
    sumspict.manage(mana3)
},cnsl=cnsl)

test_this("4.4: Old summary function still works:",{
    mansummary(mana3)
},cnsl=cnsl)

out(sumspict.manage(tmp <- man.select(mana3,1)),cnsl=cnsl)

out(sumspict.manage(tmp <- man.select(mana3,"noF")),cnsl=cnsl)

header("5: estimating TAC directly/indirectly")
######################################################

test_this("5.1: TAC from all management scenarios in rep",{
    lapply(man.tac(mana2),round,2)
},cnsl=cnsl)

test_this("5.2: TAC for one scenario",{
    round(get.TAC(fit),2)
},cnsl=cnsl)


header("6: Provoke warnings & errors & challenge functionality")
######################################################

test_this("6.1: Message when manstart and maninterval used and not equal", {
    inp <- inpA
    inp$maninterval <- c(1992,1993)
    inp$manstart <- 1991
    inp <- check.inp(inp)
    TRUE
}, cnsl=cnsl)

test_this("6.2: Message when timepredc and maninterval used and not equal", {
    inp <- inpA
    inp$maninterval <- c(1992,1993)
    inp$timepredc <- 1991
    inp <- check.inp(inp)
    TRUE
}, cnsl=cnsl)

test_this("6.3: Message when timepredi and maneval used and not equal", {
    inp <- inpA
    inp$maneval <- 1992
    inp$timepredi <- 1991
    inp <- check.inp(inp)
    TRUE
}, cnsl=cnsl)

test_this("6.4: Maninterval does not match dteuler time steps", {
    inp <- inpA
    inp$maninterval <- c(1991.573838,1992)
    inp <- check.inp(inp)
    TRUE
}, cnsl=cnsl)

test_this("6.5: Maninterval smaller than dteuler", {
    inp <- inpA
    inp$maninterval <- c(1991,1991.004)
    inp <- check.inp(inp)
    TRUE
}, cnsl=cnsl)

test_this("6.6: Warning message when manstart larger than timepredc", {
    inp <- inpA
    inp$timepredc <- 1992
    inp$dtpredc <- 2
    inp$timepredi <- 1994
    inp$manstart <- 1993
    inp <- check.inp(inp)
    TRUE
}, cnsl=cnsl)

test_this("6.7: manstart smaller than timepredc", {
    inp1 <- inpA
    inp1$timepredi <- 1993
    inp1$timepredc <- 1992
    inp1$dtpredc <- 1
    inp1$manstart <- 1991
    inp1$ffac <- 0.2
    inp1 <- check.inp(inp1)
    set.seed(124)
    fit1 <- fit.spict(inp1)
    TRUE
}, cnsl=cnsl)

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
out(fit1$Cp != fit2$Cp, cnsl=cnsl)
## But values the same just timing different:
out(fit2$Cp == get.par("logCpred",fit1,exp=TRUE)[which(inp1$timeCpred == inp2$timepredc),2],
    cnsl=cnsl)
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
}, cnsl=cnsl)

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
}, cnsl=cnsl)


header("7: plotting",cnsl=cnsl)
######################################################

plotspict.catch(mana)
plotspict.bbmsy(mana)
plot2(mana)
