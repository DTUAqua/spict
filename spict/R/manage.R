# Stochastic surplus Production model in Continuous-Time (SPiCT)
#    Copyright (C) 2015  Martin Waever Pedersen, mawp@dtu.dk or wpsgodd@gmail.com
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.



#' @name manage
#' @title Calculate predictions under different management scenarios
#' @details Scenarios that are currently implemented include:
#' \itemize{
#'   \item{"1"}{ Take a specific catch. Default catch: MSY.}
#'   \item{"2"}{ Fish at Fmsy.}
#'   \item{"3"}{ No fishing, reduce to 5\% of last F.}
#'   \item{"4"}{ Reduce F by X\%. Default X = 25.}
#'   \item{"5"}{ Increase F by X\%. Default X = 25.}
#' }
#' @param repin Result list from fit.spict().
#' @param scenarios Vector of integers specifying which scenarios to run. Default: 'all'.
#' @param dbg Debug flag, dbg=1 some output, dbg=2 more ourput.
#' @return List containing results of management calculations.
#' @export
#' @examples
#' data(pol)
#' rep <- fit.spict(pol$albacore)
#' repman <- manage(rep)
manage <- function(repin, scenarios='all', dbg=0){
    if(scenarios == 'all') scenarios <- 1:5
    
    inpin <- list()
    inpin$timeC <- repin$inp$timeC
    inpin$obsC <- repin$inp$obsC
    inpin$timeI <- repin$inp$timeI
    inpin$obsI <- repin$inp$obsI
    timelastobs <- repin$inp$time[repin$inp$indlastobs]
    if(!repin$inp$timepredc < timelastobs+1){
        # Always predict at least two years
        inpin$timepredc <- repin$inp$timepredc
        inpin$timepredi <- repin$inp$timepredi

        inp <- list()
        repman <- list() # Output list

        if(1 %in% scenarios){
            #at('Scenario 1: Take specific catch = MSY...')
            # 1. Specify the catch, which will be taken each year in the prediction period
            catch <- get.par('MSY', repin)[2]
            repman[[1]] <- take.c(catch, inpin, repin, dbg=0)
            #at('done!\n')
        }
        if(2 %in% scenarios){
            #at('Scenario 2: Fish at Fmsy...')
            # Fish at Fmsy
            Fmsy <- get.par('logFmsy', repin, exp=TRUE)[2]
            Flast <- get.par('logF', repin, exp=TRUE)[repin$inp$indpred[1], 2]
            fac2 <- Fmsy / Flast
            repman[[2]] <- prop.F(fac2, inpin, repin, dbg=dbg)
            #at('done!\n')
        }
        if(3 %in% scenarios){
            #at('Scenario 3: No fishing (reduce to 5% of last F)...')
            # No fishing, reduce to 5% of last F
            fac3 <- 0.01
            repman[[3]] <- prop.F(fac3, inpin, repin, dbg=dbg)
            #at('done!\n')
        }
        if(4 %in% scenarios){
            #at('Scenario 4: Reduce F by 25%...')
            # Reduce F by X%
            fac4 <- 0.75
            repman[[4]] <- prop.F(fac4, inpin, repin, dbg=dbg)
            #at('done!\n')
        }
        if(5 %in% scenarios){
            #at('Scenario 5: Increase F by 25%...')
            # Increase F by X%
            fac5 <- 1.25
            repman[[5]] <- prop.F(fac5, inpin, repin, dbg=dbg)
            #at('done!\n')
        }
        repin$man <- repman
    } else {
        stop('Error: Could not do management calculations because prediction horizon is too short. Increase inp$timepredc to be at least one year into the future.\n')
    }
    return(repin)
}


#' @name prop.F
#' @title Calculate management for changing F by a given factor.
#' @param fac Factor to multiply current F with.
#' @param inpin Input list.
#' @param repin Results list.
#' @param dbg Debug flag, dbg=1 some output, dbg=2 more ourput.
#' @return List containing results of management calculations.
prop.F <- function(fac, inpin, repin, dbg=0){
    inpt <- check.inp(inpin)
    plt <- repin$obj$env$parList(repin$opt$par)
    datint <- make.datin(inpt, dbg=dbg)
    inds <- inpt$indpred
    ninds <- length(inds)
    ffacvec <- rep(1, ninds)
    # Set F fac
    ffacvec[1] <- fac
    datint$ffacvec[inds] <- ffacvec
    objt <- make.obj(datint, plt, inpt, phase=1)
    objt$fn(repin$opt$par)
    repmant <- sdreport(objt)
    repmant$inp <- inpt
    repmant$obj <- objt
    return(repmant)
}


#' @name take.c
#' @title Calculate management when taking a constant catch (proxy for setting a TAC).
#' @param catch Annual catch to take in the prediction period.
#' @param inpin Input list.
#' @param repin Results list.
#' @param dbg Debug flag, dbg=1 some output, dbg=2 more ourput.
#' @return List containing results of management calculations.
take.c <- function(catch, inpin, repin, dbg=0){
    inpc <- check.inp(inpin)
    inpt <- inpin
    plt <- repin$obj$env$parList(repin$opt$par)
    timecatch <- annual(inpc$time[inpc$indpred], numeric(length(inpc$indpred)))$anntime
    ncatch <- length(timecatch)
    obscatch <- rep(catch, ncatch)
    inpt$timeC <- c(inpt$timeC, timecatch)
    inpt$obsC <- c(inpt$obsC, obscatch)
    inpt <- check.inp(inpt)
    datint <- make.datin(inpt, dbg=dbg)
    objt <- make.obj(datint, plt, inpt, phase=1)
    objt$fn(repin$opt$par)
    repmant <- sdreport(objt)
    repmant$inp <- inpt
    repmant$obj <- objt
    return(repmant)
}


#' @name mansummary
#' @title Print management summary.
#' @param rep Result list as output from manage().
#' @param ypred Show results for ypred years into the future.
#' @return Data frame containing management summary.
#' @export
mansummary <- function(rep, ypred=1){
    repman <- rep$man
    get.pdelta <- function(rep, repman, parname='logB'){
        val <- get.par(parname, rep, exp=TRUE)[indstart, 2]
        val1 <- get.par(parname, repman, exp=TRUE)[indnext, 2]
        return(round((val1 - val)/val*100, 1))
    }
    #ypred <- 2
    indstart <- rep$inp$indpred[1] # Current time
    curtime <- rep$inp$time[indstart] # We are at 1 January this year
    indnext <- which(rep$inp$time == curtime+ypred) # Current time + 1 year
    if(length(indnext)==1){
        indnextC <- which(rep$inp$timeCpred == curtime+ypred-1)
        nsc <- length(repman)
        Cnextyear <- numeric(nsc)
        Bnextyear <- numeric(nsc)
        Fnextyear <- numeric(nsc)
        perc.dB <- numeric(nsc)
        perc.dF <- numeric(nsc)
        for(i in 1:nsc){
            perc.dB[i] <- get.pdelta(rep, repman[[i]], parname='logB')
            perc.dF[i] <- get.pdelta(rep, repman[[i]], parname='logF')
            Cnextyear[i] <- round(get.par('logCpred', repman[[i]], exp=TRUE)[indnextC, 2], 1)
            Bnextyear[i] <- get.par('logB', repman[[i]], exp=TRUE)[indnext, 2]
            Fnextyear[i] <- get.par('logF', repman[[i]], exp=TRUE)[indnext, 2]
        }
        Cn <- paste0('C', curtime+ypred-1)
        Bn <- paste0('B', curtime+ypred)
        Fn <- paste0('F', curtime+ypred)
        BBn <- paste0('B/Bmsy', curtime+ypred)
        FFn <- paste0('F/Fmsy', curtime+ypred)
        df <- list()
        df[[Cn]] <- Cnextyear
        df[[Bn]] <- round(Bnextyear, 1)
        df[[Fn]] <- round(Fnextyear, 3)
        Bmsy <- get.par('logBmsy', rep, exp=TRUE)[2]
        df[[BBn]] <- round(Bnextyear/Bmsy, 2)
        Fmsy <- get.par('logFmsy', rep, exp=TRUE)[2]
        df[[FFn]] <- round(Fnextyear/Fmsy, 2)
        df <- cbind(as.data.frame(df), perc.dB, perc.dF)
        rn <- c('1: Specify catch (MSY)', '2: Fish at Fmsy', '3: No fishing', '4: reduce F X% (25)', '5: increase F X% (25)')
        rownames(df) <- rn
        print(df)
        invisible(df)
    } else {
        cat('Warning: Could not show management results because ypred is larger than the calculated management time frame. Reduce ypred.\n')
    }
}
