## Testing new management functionality
## T.K. Mildenberger <t.k.mildenberger@gmail.com>
## 19/01/2020

## flag for printing to consol (dev mode)
cnsl = TRUE

## set seed
set.seed(123)

## load spict
## devtools::install_github("DTUAqua/spict/spict")  ## REMOVE:
require(spict)

## load testmore functions
source("../../funcs.R")

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

reps <- fit.spict(inps)

mans <- manage(reps, scenarios=c(1,6))

plot2(reps)
dev.print(pdf,"plot2.spictclsSeasonal.pdf")


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

out("All management variables have the same default as spict v1.2.8:")
out(inpa$timepredi == timepredi.v128, cnsl=cnsl)
out(inpa$timepredc == timepredc.v128, cnsl=cnsl)
out(inpa$dtpredc == dtpredc.v128, cnsl=cnsl)
out(inpa$manstart == manstart.v128, cnsl=cnsl)
out(inpa$maninterval[1] == manstart.v128, cnsl=cnsl)
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

out("All management variables have the same default as spict v1.2.8:")
out(inps$timepredc == timepredc.v128, cnsl=cnsl)
out(inps$manstart == manstart.v128, cnsl=cnsl)
out(inps$maninterval[1] == manstart.v128, cnsl=cnsl)
out(inps$dtpredc == dtpredc.v128, cnsl=cnsl)
out(inps$timepredi == timepredi.v128, cnsl=cnsl)
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
inp <- check.inp(inp)

out(inp$maneval == inp$timepredi, cnsl=cnsl)
out(all.equal(inp$maninterval,c(inp$timepredc,inp$timepredc+inp$dtpredc)), cnsl=cnsl)
out(inp$maninterval[1] == inp$manstart, cnsl=cnsl)

out("using manstart instead of timepredc")
inp <- inpA
inp$manstart <- 1992
inp <- check.inp(inp)

out(inp$maneval == inp$timepredi, cnsl=cnsl)
out(all.equal(inp$maninterval,c(inp$timepredc,inp$timepredc+inp$dtpredc)), cnsl=cnsl)
out("This does not work for maninterval[1] == manstart:", cnsl=cnsl)
out(inp$maninterval[1] == inp$manstart, cnsl=cnsl)
out("Because maninterval[1] is set by default to timepredc, which makes more sense as maninterval replaces manstart. However, this could cause some changes in scripts which uses manstart and timepredc and set those variables to different values.", cnsl=cnsl)


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
inp <- check.inp(inp)

out(inp$maneval == inp$timepredi, cnsl=cnsl)
out(all.equal(inp$maninterval,c(inp$timepredc,inp$timepredc+inp$dtpredc)), cnsl=cnsl)
out(inp$maninterval[1] == inp$manstart, cnsl=cnsl)

out("using manstart instead of timepredc")
inp <- inpS
inp$manstart <- 42
inp <- check.inp(inp)

out(inp$maneval == inp$timepredi, cnsl=cnsl)
out(all.equal(inp$maninterval,c(inp$timepredc,inp$timepredc+inp$dtpredc)), cnsl=cnsl)
out("This does not work for maninterval[1] == manstart:", cnsl=cnsl)
out(inp$maninterval[1] == inp$manstart, cnsl=cnsl)
out("Because maninterval[1] is set by default to timepredc, which makes more sense as maninterval replaces manstart. However, this could cause some changes in scripts which uses manstart and timepredc and set those variables to different values.", cnsl=cnsl)

header("2: check man time", cnsl=cnsl)
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


header("2: timeline")
######################################################

out(man.timeline(inpa),cnsl=cnsl)
out(man.timeline(inpb),cnsl=cnsl)

out(man.timeline(resa),cnsl=cnsl)
out(man.timeline(resb),cnsl=cnsl)


header("3: manage / summary / timeline")
######################################################

mana <- manage(resb, scenarios=c(1,2))

mana <- manage(resb, scenarios=c("Fmsy","currentCatch"))

plotspict.catch(mana)

inp <- pol$albacore
inp$maninterval <- c(1990,1991)
inp$maneval <- 1991
inp$ffac <- 0.75
rep <- fit.spict(inp)

mana <- manage(rep, scenarios=c("Fmsy","currentCatch"),maninterval = c(1991,1992))

out(man.timeline(mana),cnsl=cnsl)

out("Management periods do not match between scenarios",cnsl=cnsl)
mana1 <- add.man.scenario(mana, maninterval = c(1991.5,1992.5))

out(man.timeline(mana1),cnsl=cnsl)
## BUG: shows management period!

out("Scenario with continue F vs continue C in intermediate period",cnsl=cnsl)
mana2 <- add.man.scenario(mana, intermediatePeriodCatch = 12)

out(man.timeline(mana2),cnsl=cnsl)


rep2 <- manage(rep, c(1,4),maninterval = c(1993,1995))

## BUG: mansumaary error about ypred!!!

## manage
test_this("1.1: run manage all scenarios",{
    repman <- manage(rep)
})

test_this("1.2: print manage summary",{
    sumspict.manage(repman)
})

test_this("1.3: use deprecated summary function",{
    mansummary(repman)
})

test_this("1.4: Make selection of management scenarios",{
    sumspict.manage(man.select(repman,c(1,5,6,8)))
})


test_this("1.5: Change base management scenario (manbase and of rep)" ,{
    names(repman$man)
    repman2 <- manbase.select(repman,"ices")
    all.equal(sumspict.predictions(repman2),
              sumspict.predictions(repman$man[[8]]), tolerance = 0.01)
    all.equal(sumspict.predictions(repman2$manbase),
              sumspict.predictions(repman$man[[8]]), tolerance = 0.01)

})


## testing adding scenarios
test_this("1.6: Adding scenarios",{
    repman1 <- add.man.scenario(repman)
    repman1 <- add.man.scenario(repman1, ffac = 0.5)
    repman1 <- add.man.scenario(repman1, ffac = 0.25)
})

## scenarios should not be overwritten:
out(length(repman1$man))

## names should bemeaningful:
out(names(repman1$man))

## sumspict should still work
out(sumspict.manage(repman1))

## change management interval:
test_this("1.7: Change management interval" ,{
    inp <- pol$albacore
    inp$maninterval <- c(1992,1993)
    res <- fit.spict(inp)
    repman <- manage(res, c(1,5))
    ## should be the same as this:
    repman2 <- manage(res, c(1,5), maninterval = c(1992,1993))

    all.equal(sumspict.predictions(repman),
              sumspict.predictions(repman2), tolerance = 0.01)
})


# TODO: manage()

# TODO: printing/summary()

rep4 <- man.select(rep2)

rep4 <- man.select(rep2, 2)

rep4 <- man.select(rep2, "noF")
names(rep4$man)

## check summary printing (problems with manbase?)
## dF, dB working correctly?

header("3: adding manage scenarios")
######################################################

rep3 <- add.man.scenario(rep2, maninterval=c(1991,1992),
                         maneval=1993, ffac=2,
                         scenarioTitle = "moreFishing")

# TODO: add.man.scenario()
# TODO: printing/summary()


rep

rep2 <- add.man.scenario(rep, "loadsFishing", maninterval = c(1992,1993),
                         catchIntermediatePeriod = 50)
rep2 <- add.man.scenario(rep2, "loadsFishing2", maninterval = c(1992,1993),
                         ffac=2)

rep2$man

plotspict.catch(rep2)
rep2$man[[1]]$inp$timeCpred

inp <- check.inp(pol$albacore)
inp$reportmode <- 1
repx <- fit.spict(inp)
plot(repx)

test <- manage(repx)



## report tac only
test_this("1.8: Calculate TAC for fitted spict object without $man" ,{
    add.man.scenario(res, getFit=FALSE)
})


test_this("1.9: Add scenario to fitted spict object" ,{
    res <- add.man.scenario(res, cfac=0.5, scenarioTitle = "Half catch")
})

## how is the default dealt with?
## what is catchIntermediatePeriod based on as default?

header("4: estimating TAC directly/indirectly")
######################################################

# TODO: calc.tac()

## report tac only
test_this("1.8: Calculate the TAC for all add.man.scenarios" ,{
    man.tac(repman1)
})


header("4: Challenges")
######################################################

## manstart and maninerval specified
## maninterval + timepredc specified
## maneval and timepredi specified
## management time outside of dteuler time steps



rep4 <- man.select(inp)

## wrong input to functions, creating warnings...

## in combination with timevarying growth, MSYregimes, etc.

# TODO: error/ wanrings

header("4: plotting")
######################################################

plotspict.catch(rep2)
plotspict.bbmsy(rep2)



## tests
## TODO: all tests should be numbered
## TODO: reduce scenario testing
## TODO: 2 intermediate years, 2 man years, ...
