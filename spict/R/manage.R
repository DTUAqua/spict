# Stochastic surplus Production model in Continuous-Time (SPiCT)
#    Copyright (C) 2015-2016  Martin W. Pedersen, mawp@dtu.dk, wpsgodd@gmail.com
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
#'   \item{"1"}{ Keep the catch of the current year (i.e. the last observed catch).}
#'   \item{"2"}{ Keep the F of the current year.}
#'   \item{"3"}{ Fish at Fmsy i.e. F=Fmsy.}
#'   \item{"4"}{ No fishing, reduce to 1\% of current F.}
#'   \item{"5"}{ Reduce F by X\%. Default X = 25.}
#'   \item{"6"}{ Increase F by X\%. Default X = 25.}
#' }
#' @param repin Result list from fit.spict().
#' @param scenarios Vector of integers specifying which scenarios to run. Default: 'all'.
#' @param manstart Year that management should be initiated.
#' @param dbg Debug flag, dbg=1 some output, dbg=2 more ourput.
#' @return List containing results of management calculations.
#' @export
#' @examples
#' data(pol)
#' rep <- fit.spict(pol$albacore)
#' repman <- manage(rep)
manage <- function(repin, scenarios='all', manstart=NULL, dbg=0){
    if (scenarios == 'all'){
        scenarios <- 1:6
    }
    if (is.null(manstart)){
        manstart <- repin$inp$manstart
    } else {
        repin$inp$manstart <- manstart
    }
    maninds <- which(repin$inp$time >= manstart)
    # inpin is a list containing only observations (later prediction horizons are added)
    inpin <- list()
    inpin$dteuler <- repin$inp$dteuler
    inpin$timeC <- repin$inp$timeC
    inpin$obsC <- repin$inp$obsC
    inpin$timeI <- repin$inp$timeI
    inpin$obsI <- repin$inp$obsI
    timelastobs <- repin$inp$time[repin$inp$indlastobs]
    if (!repin$inp$timepredc < timelastobs+1){
        # Always predict at least two years
        inpin$timepredc <- repin$inp$timepredc
        inpin$timepredi <- repin$inp$timepredi
        inp <- list()
        repman <- list() # Output list
        if (1 %in% scenarios){
            # 1. Specify the catch, which will be taken each year in the prediction period
            catch <- tail(inpin$obsC, 1)
            repman[[1]] <- take.c(catch, inpin, repin, maninds, dbg=dbg)
        }
        if (2 %in% scenarios){
            # Keep current F
            fac2 <- 1.0
            repman[[2]] <- prop.F(fac2, inpin, repin, maninds, dbg=dbg)
        }
        if (3 %in% scenarios){
            # Fish at Fmsy
            Fmsy <- get.par('logFmsy', repin, exp=TRUE)[2]
            Flast <- get.par('logF', repin, exp=TRUE)[repin$inp$indpred[1], 2]
            fac3 <- Fmsy / Flast
            repman[[3]] <- prop.F(fac3, inpin, repin, maninds, dbg=dbg)
        }
        if (4 %in% scenarios){
            # No fishing, reduce to 0.1% of last F
            fac4 <- 0.001
            repman[[4]] <- prop.F(fac4, inpin, repin, maninds, dbg=dbg)
        }
        if (5 %in% scenarios){
            # Reduce F by X%
            fac5 <- 0.75
            repman[[5]] <- prop.F(fac5, inpin, repin, maninds, dbg=dbg)
        }
        if (6 %in% scenarios){
            # Increase F by X%
            fac6 <- 1.25
            repman[[6]] <- prop.F(fac6, inpin, repin, maninds, dbg=dbg)
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
#' @param maninds Indices of time vector for which to apply management.
#' @param corF Make correction to F process such that the drift (-0.5*sdf^2*dt) is cancelled and F remains constant in projection mode
#' @param dbg Debug flag, dbg=1 some output, dbg=2 more ourput.
#' @return List containing results of management calculations.
prop.F <- function(fac, inpin, repin, maninds, corF=FALSE, dbg=0){
    inpt <- check.inp(inpin)
    plt <- repin$obj$env$parList(repin$opt$par)
    datint <- make.datin(inpt, dbg=dbg)
    maninds <- inpt$indpred
    nmaninds <- length(maninds)
    # Set F fac
    if (corF){
        # This is to compensate for the -0.5*sdf^2*dt term in the logF process
        sdf <- get.par('logsdf', repin, exp=TRUE)[2]
        ffacvec <- exp(0.5 * sdf^2 * inpt$dt[maninds])
    } else {
        ffacvec <- rep(1, nmaninds)
    }
    ffacvec[1] <- ffacvec[1] * fac
    datint$ffacvec[maninds] <- ffacvec
    # Make object
    objt <- make.obj(datint, plt, inpt, phase=1)
    objt$fn(repin$opt$par)
    repmant <- sdreport(objt)
    repmant$inp <- inpt
    repmant$obj <- objt
    if (!is.null(repmant)){
        class(repmant) <- "spictcls"
    }
    return(repmant)
}


#' @name take.c
#' @title Calculate management when taking a constant catch (proxy for setting a TAC).
#' @param catch Annual catch to take in the prediction period.
#' @param inpin Input list.
#' @param repin Results list.
#' @param maninds Indices of time vector for which to apply management.
#' @param dbg Debug flag, dbg=1 some output, dbg=2 more ourput.
#' @return List containing results of management calculations.
take.c <- function(catch, inpin, repin, maninds, dbg=0){
    inpc <- check.inp(inpin)
    inpt <- inpin
    plt <- repin$obj$env$parList(repin$opt$par)
    nmaninds <- length(maninds)
    timecatch <- annual(inpc$time[maninds[-nmaninds]],
                        numeric(length(maninds[-nmaninds])))$anntime
    ncatch <- length(timecatch)
    obscatch <- rep(catch, ncatch)
    inpt$timeC <- c(inpt$timeC, timecatch)
    inpt$obsC <- c(inpt$obsC, obscatch)
    inpt <- check.inp(inpt)
    # Make TMB data and object
    datint <- make.datin(inpt, dbg=dbg)
    objt <- make.obj(datint, plt, inpt, phase=1)
    # Get updated sd report
    objt$fn(repin$opt$par)
    repmant <- sdreport(objt)
    repmant$inp <- inpt
    repmant$obj <- objt
    if (!is.null(repmant)){
        class(repmant) <- "spictcls"
    }
    return(repmant)
}


#' @name mansummary
#' @title Print management summary.
#' @param rep Result list as output from manage().
#' @param ypred Show results for ypred years into the future.
#' @param include.EBinf Include EBinf/Bmsy in the output.
#' @param include.unc Include uncertainty of management quantities.
#' @param verbose Print more details on observed and predicted time intervals.
#' @return Data frame containing management summary.
#' @export
mansummary <- function(rep, ypred=1, include.EBinf=FALSE, include.unc=TRUE, verbose=TRUE){
    repman <- rep$man
    # Calculate percent difference.
    get.pdelta <- function(rep, repman, indstart, indnext, parname='logB'){
        val <- get.par(parname, rep, exp=TRUE)[indstart, 2]
        val1 <- get.par(parname, repman, exp=TRUE)[indnext, 2]
        return(round((val1 - val)/val*100, 1))
    }
    indstart <- which(rep$inp$time == rep$inp$manstart) - 1
    #indstart <- rep$inp$indpred[1]-1 # Current time (last time interval of last year)
    #curtime <- rep$inp$time[indstart+1]
    curtime <- rep$inp$manstart
    indnext <- which(rep$inp$time == curtime+ypred) # Current time + ypred
    if (length(indnext) == 1){
        Cn <- paste0('C')
        Bn <- paste0('B')
        Fn <- paste0('F')
        get.cn <- function(nn){
            nc <- nchar(nn)
            tl <- 7 # Total length
            # Add spaces
            #pad <- ifelse(include.unc, paste0(rep(' ', max(0, tl-nc)), collapse=''), '')
            pad <- ''
            return(c(paste0(nn, '.lo'), paste0(pad, nn), paste0(nn, '.hi')))
        }
        BBn <- paste0('BqBmsy') # Should use / instead of q, but / is not accepted in varnames
        FFn <- paste0('FqFmsy')
        EBinfBn <- paste0('EBinfqBmsy')
        indnextC <- which((rep$inp$timeCpred+rep$inp$dtcp) == curtime+ypred)
        nsc <- length(repman)
        Cnextyear <- matrix(0, nsc, 3)
        colnames(Cnextyear) <- get.cn(Cn)
        Bnextyear <- matrix(0, nsc, 3)
        colnames(Bnextyear) <- get.cn(Bn)
        Fnextyear <- matrix(0, nsc, 3)
        colnames(Fnextyear) <- get.cn(Fn)
        BBnextyear <- matrix(0, nsc, 3)
        colnames(BBnextyear) <- get.cn(BBn)
        FFnextyear <- matrix(0, nsc, 3)
        colnames(FFnextyear) <- get.cn(FFn)
        perc.dB <- numeric(nsc)
        perc.dF <- numeric(nsc)
        EBinf <- numeric(nsc)
        for(i in 1:nsc){
            EBinf[i] <- get.EBinf(repman[[i]])
            perc.dB[i] <- get.pdelta(rep, repman[[i]], indstart, indnext, parname='logB')
            perc.dF[i] <- get.pdelta(rep, repman[[i]], indstart, indnext, parname='logF')
            Cnextyear[i, ] <- round(get.par('logCpred', repman[[i]], exp=TRUE)[indnextC, 1:3], 1)
            Bnextyear[i, ] <- round(get.par('logB', repman[[i]], exp=TRUE)[indnext, 1:3], 1)
            Fnextyear[i, ] <- round(get.par('logF', repman[[i]], exp=TRUE)[indnext, 1:3], 3)
            BBnextyear[i, ] <- round(get.par('logBBmsy', repman[[i]], exp=TRUE)[indnext, 1:3], 3)
            FFnextyear[i, ] <- round(get.par('logFFmsy', repman[[i]], exp=TRUE)[indnext, 1:3], 3)
        }
        FBtime <- fd(curtime+ypred)
        Ctime1 <- fd(rep$inp$timeCpred[indnextC])
        Ctime2 <- fd(rep$inp$timeCpred[indnextC]+rep$inp$dtcp[indnextC])
        if (!verbose){
            Cn <- paste0('C', Ctime1)
            Bn <- paste0('B', FBtime)
            Fn <- paste0('F', FBtime)
        }
        # Data frame with predictions
        df <- cbind(Cnextyear[, 2], Bnextyear[, 2], Fnextyear[, 2], BBnextyear[, 2],
                    FFnextyear[, 2], perc.dB, perc.dF)
        colnames(df)[1:5] <- c(Cn, Bn, Fn, BBn, FFn)
        qinds <- grep('q', colnames(df))
        colnames(df)[qinds] <- sub('q', '/', colnames(df)[qinds]) # Replace q with /
        # Data frame with uncertainties of absolute predictions
        inds <- c(1, 3)
        dfabs <- cbind(Cnextyear[, inds], Bnextyear[, inds], Fnextyear[, inds])
        colnames(dfabs) <- c(colnames(Cnextyear)[inds], colnames(Bnextyear)[inds],
                             colnames(Fnextyear)[inds])
        # Data frame with uncertainties of relateive predictions
        dfrel <- cbind(BBnextyear[, inds], FFnextyear[, inds])
        colnames(dfrel) <- c(colnames(BBnextyear)[inds], colnames(FFnextyear)[inds])
        qinds <- grep('q', colnames(dfrel))
        colnames(dfrel)[qinds] <- sub('q', '/', colnames(dfrel)[qinds]) # Replace q with /
        # Set row names
        rn <- c('1. Keep current catch', '2. Keep current F', '3. Fish at Fmsy',
                '4. No fishing', '5. Reduce F 25%', '6. Increase F 25%')
        rownames(df) <- rn
        rownames(dfrel) <- rn
        rownames(dfabs) <- rn
        #cat('Management summary\n')
        timerangeI <- range(unlist(rep$inp$timeI))
        timerangeC <- range(rep$inp$timeC)
        lastcatchseen <- tail(rep$inp$timeC+rep$inp$dtc, 1)
        # Start printing stuff
        if (verbose){ # Time interval information
            cat(paste0('Observed interval, index:  ',
                       fd(timerangeI[1]),
                       ' - ',
                       fd(timerangeI[2]),
                       '\n'))
            cat(paste0('Observed interval, catch:  ',
                       fd(timerangeC[1]),
                       ' - ',
                       fd(lastcatchseen),
                       '\n\n'))
            cat(paste0('Fishing mortality (F) prediction: ',
                       FBtime, '\n'))
            cat(paste0('Biomass (B) prediction:           ',
                       FBtime, '\n'))
            cat(paste0('Catch (C) prediction interval:    ',
                       Ctime1,
                       ' - ',
                       Ctime2,
                       '\n\n'))
            if (rep$inp$catchunit != ''){
                cat(paste('Catch/biomass unit:', rep$inp$catchunit, '\n\n'))
            }
            cat('Predictions\n')
        }
        print(df)
        if (include.unc){
            cat('\n95% CIs of absolute predictions\n')
            print(dfabs)
            cat('\n95% CIs of relative predictions\n')
            print(dfrel)
        }
        invisible(df)
    } else {
        cat('Warning: Could not show management results because ypred is larger than the calculated management time frame. Reduce ypred.\n')
    }
}
