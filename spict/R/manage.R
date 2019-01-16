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
#' mansummary(repman) # To print projections
manage <- function(repin, scenarios='all', manstart=NULL, dbg=0, catch=NULL, catchList=NULL){
    if (scenarios == 'all'){
        scenarios <- 1:6
    }
    if (is.null(manstart)){
        manstart <- repin$inp$manstart
    } else {
        repin$inp$manstart <- manstart
    }
    maninds <- which(repin$inp$time >= manstart)
    inpin <- repin$inp

    timelastobs <- repin$inp$time[repin$inp$indlastobs]
    if (!repin$inp$timepredc < timelastobs){
        # Add prediction horizons
        inpin$timepredc <- repin$inp$timepredc
        inpin$timepredi <- repin$inp$timepredi
        inpin$manstart <- repin$inp$manstart
        repman <- list() # Output list
        attr(repman, "scenarios") <- scenarios 
        if (1 %in% scenarios){
            # 1. Specify the catch, which will be taken each year in the prediction period
            lastyearidxs <- min( which( cumsum(rev(inpin$dtc))>=1 ) ) ## warning: this will not make sense with subannual/mixed data with missing values
            if(is.null(catch)) catch <- sum(tail(inpin$obsC, lastyearidxs))
            repman[[1]] <- take.c(catch, inpin, repin, dbg=dbg, catchList=catchList)
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
        # Create an baseline F trajectory with constant F and store
        repin$manbase <- prop.F(fac=1, inpin, repin, maninds, dbg=dbg)
    } else {
        stop('Error: Could not do management calculations because prediction horizon is too short. Increase inp$timepredc to be at least one timestep into the future.\n')
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
#' @export
prop.F <- function(fac, inpin, repin, maninds, corF=FALSE, dbg=0){
    inpt <- check.inp(inpin)
    inpt <- make.ffacvec(inpt, fac)
    # Make object
    plt <- repin$obj$env$parList(repin$opt$par)
    datint <- make.datin(inpt, dbg=dbg)
    objt <- make.obj(datint, plt, inpt, phase=1)
    objt$fn(repin$opt$par)
    ## repmant <- sdreport(objt)
    verflag <- as.numeric(gsub('[.]', '', as.character(packageVersion('TMB')))) >= 171
    if (verflag) { 
      repmant <- sdreport(objt,
                          getJointPrecision=repin$inp$getJointPrecision,
                          bias.correct=repin$inp$bias.correct,
                          bias.correct.control=repin$inp$bias.correct.control,
                          getReportCovariance=repin$inp$getReportCovariance)
    } else {
      repmant <-sdreport(objt,
                         getJointPrecision=repin$inp$getJointPrecision,
                         bias.correct=repin$inp$bias.correct,
                         bias.correct.control=repin$inp$bias.correct.control)
    }
    repmant$inp <- inpt
    repmant$obj <- objt
    repmant$opt <- list(convergence=0)
    if (!is.null(repmant)){
        class(repmant) <- "spictcls"
    }
    return(repmant)
}


#' @name take.c
#' @title Calculate management when taking a constant catch (proxy for setting a TAC).
#' @param catch Take this catch 'dtpredc' ahead from manstart time 
#' @param inpin Input list.
#' @param repin Results list.
#' @param dbg Debug flag, dbg=1 some output, dbg=2 more output.
#' @param sdfac Take catch with this 'stdevfacC' (default = 1e-3) 
#' @return List containing results of management calculations.
#' @export
take.c <- function(catch, inpin, repin, dbg=0, sdfac=1e-3, catchList=NULL){
    
    inpt <- inpin
    if(is.null(catchList)){
        tmpTime <- repin$inp$timeCpred  
        maninds <- which(tmpTime >= inpin$manstart)
        inpt$timeC <- c( inpt$timeC, tmpTime[maninds] )
        inpt$obsC <- c( inpt$obsC, rep(catch, length(maninds)) )
        inpt$stdevfacC <- c(inpt$stdevfacC, rep(sdfac, length(maninds)) )  
        inpt$dtc <- c(inpt$dtc, rep(inpt$dtpredc, length(maninds)) )
    } else {
        inpt$timeC <- c( inpt$timeC, catchList$timeC )
        inpt$obsC <- c( inpt$obsC, catchList$obsC )
        if(is.null(catchList$stdevfacC))
            inpt$stdevfacC <- c(inpt$stdevfacC, rep(sdfac, length(catchList$timeC)) )  else
            inpt$stdevfacC <- c(inpt$stdevfacC, catchList$stdevfacC)
        
        inpt$dtc <- c(inpt$dtc, catchList$dtc )
    }

    
    inpt <- check.inp(inpt)
    # Make TMB data and object
    plt <- repin$obj$env$parList(repin$opt$par)
    datint <- make.datin(inpt, dbg=dbg)
    objt <- make.obj(datint, plt, inpt, phase=1)
    # Get updated sd report
    objt$fn(repin$opt$par)
    repmant <- sdreport(objt)
    repmant$inp <- inpt
    repmant$obj <- objt
    repmant$opt <- list(convergence=0)
    if (!is.null(repmant)){
        class(repmant) <- "spictcls"
    }
    return(repmant)
}


#' @name mansummary
#' @title Print management summary.
#' @param repin Result list as output from manage().
#' @param ypred Show results for ypred years from manstart.
#' @param include.EBinf Include EBinf/Bmsy in the output.
#' @param include.unc Include uncertainty of management quantities.
#' @param verbose Print more details on observed and predicted time intervals.
#' @return Data frame containing management summary.
#' @export
mansummary <- function(repin, ypred=1, include.EBinf=FALSE, include.unc=TRUE, verbose=TRUE){
    if (!'man' %in% names(repin)){
        stop('Management calculations not found, run manage() to include them.')
    } else {
        repman <- repin$man
        rep <- repin$manbase
        # Calculate percent difference.
        get.pdelta <- function(rep, repman, indstart, indnext, parname='logB'){
            val <- get.par(parname, rep, exp=TRUE)[indstart, 2]
            val1 <- get.par(parname, repman, exp=TRUE)[indnext, 2]
            return(round((val1 - val)/val*100, 1))
        }
        indstart <- which(rep$inp$time == rep$inp$manstart)
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

            scenarios <- attr(repman,"scenarios")
            nsc <- length( scenarios )
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
                rp <- repman[[ scenarios[i] ]]  ##repman[[i]]
                EBinf[i] <- get.EBinf(rp)
                perc.dB[i] <- get.pdelta(rep, rp, indstart, indnext, parname='logB')
                perc.dF[i] <- get.pdelta(rep, rp, indstart, indnext, parname='logF')
                indnextC <- which((rp$inp$timeCpred + rp$inp$dtcp) == curtime+ypred)
                Cnextyear[i, ] <- round(get.par('logCpred', rp, exp=TRUE)[indnextC, 1:3], 1)
                Bnextyear[i, ] <- round(get.par('logB', rp, exp=TRUE)[indnext, 1:3], 1)
                Fnextyear[i, ] <- round(get.par('logF', rp, exp=TRUE)[indnext, 1:3], 3)
                BBnextyear[i, ] <- round(get.par('logBBmsy', rp, exp=TRUE)[indnext, 1:3], 3)
                FFnextyear[i, ] <- round(get.par('logFFmsy', rp, exp=TRUE)[indnext, 1:3], 3)
            }
            indnextCrep <- which((rep$inp$timeCpred+rep$inp$dtcp) == curtime+ypred)
            FBtime <- fd(curtime+ypred)
            Ctime1 <- fd(rep$inp$timeCpred[indnextCrep])
            Ctime2 <- fd(rep$inp$timeCpred[indnextCrep]+rep$inp$dtcp[indnextCrep])
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
            dfabs <- cbind(Cnextyear[, inds,drop=FALSE], Bnextyear[, inds,drop=FALSE], Fnextyear[, inds,drop=FALSE])
            colnames(dfabs) <- c(colnames(Cnextyear)[inds], colnames(Bnextyear)[inds],
                                 colnames(Fnextyear)[inds])
            # Data frame with uncertainties of relateive predictions
            dfrel <- cbind(BBnextyear[, inds,drop=FALSE], FFnextyear[, inds,drop=FALSE])
            colnames(dfrel) <- c(colnames(BBnextyear)[inds], colnames(FFnextyear)[inds])
            qinds <- grep('q', colnames(dfrel))
            colnames(dfrel)[qinds] <- sub('q', '/', colnames(dfrel)[qinds]) # Replace q with /
            # Set row names
            scenarios <- attr(repman, "scenarios")
            rn <- c('1. Keep current catch', '2. Keep current F', '3. Fish at Fmsy',
                    '4. No fishing', '5. Reduce F 25%', '6. Increase F 25%')[scenarios]
            
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
            cat('Warning: Could not show management results because ypred is larger than the calculated management time frame. Reduce ypred or increase inp$timepredc and run fit.spict() and manage() again.\n')
        }
    }
}


#' @name pred.catch
#' @title Predict the catch of the prediction interval specified in
#'     inp
#' @param rep Result list as output from fit.spict().
#' @param fmsyfac Projection are made using F = fmsyfac * Fmsy.
#' @param MSEmode logical; if TRUE (default) only rel predicted states
#'     and catches are ADreported
#' @param get.sd Get uncertainty of the predicted catch.
#' @param exp If TRUE report exp of log predicted catch.
#' @param dbg Debug flag, dbg=1 some output, dbg=2 more ourput.
#' @return A vector containing predicted catch (possibly with
#'     uncertainty).
#' @export
pred.catch <- function(repin, fmsyfac=1, MSEmode = TRUE, get.sd=FALSE, exp=FALSE, dbg=0){    
    inpin <- list()
    inpin$dteuler <- repin$inp$dteuler
    inpin$timeC <- repin$inp$timeC
    inpin$obsC <- repin$inp$obsC
    inpin$timeI <- repin$inp$timeI
    inpin$obsI <- repin$inp$obsI
    timelastobs <- repin$inp$time[repin$inp$indlastobs]
    # Always predict at least two years
    inpin$timepredc <- repin$inp$timepredc
    inpin$timepredi <- repin$inp$timepredi
    Fmsy <- get.par('logFmsy', repin, exp=TRUE)[2]
    Flast <- get.par('logF', repin, exp=TRUE)[repin$inp$indpred[1], 2]
    fac <- (fmsyfac + 1e-6) * Fmsy / Flast
    inpt <- check.inp(inpin)
    # Set F fac
    inpt <- make.ffacvec(inpt, fac)
    # Make object
    datint <- make.datin(inpt, dbg=dbg)
    datint$MSEmode <- MSEmode 
    plt <- repin$obj$env$parList(repin$opt$par)
    objt <- make.obj(datint, plt, inpt, phase=1)
    objt$fn(repin$opt$par)
    if (get.sd){
        repmant <- sdreport(objt)
        Cp <- get.par('logCp', repmant, exp=exp)
    } else {
        Cp <- c(NA, log(objt$report()$Cp), NA, NA, NA)
        if (exp){
            Cp <- exp(Cp)
        }
        names(Cp) <- names(get.par('logK', repin)) # Just to get the names
    }
    return(Cp)
}


#' @name get.TAC
#' @title Estimate Total Allowable Catch (TAC)
#' @param repin Result list as output from fit.spict().
#' @param reps The number of stochastic samples of the TAC
#'     recommendation (not used for this HCR).
#' @param fractileC The fractile of the catch distribution to be used
#'     for setting the TAC. Default is median (0.5).
#' @param fractileFFmsy The fractile of the distribution of
#'     F/Fmsy. Default is 0.5 (median).
#' @param fractileBBmsy The fractile of the distribution of
#'     B/Bmsy. Default is 0.5 (median).
#' @param pa Logical; indicating if the precautionary approach should
#'     be applied (reduce F if P(B<Blim) < prob). Default is FALSE.
#' @param prob Probability for the precautionary approach (see
#'     argument 'pa', default is 0.95).
#' @param bbmsyfrac Fraction of B/Bmsy for precautionary approach
#' @param stabilityClause Logical; If true F multiplication factor is
#'     bound between two values set in lower and upper. Default:
#'     FALSE.
#' @param lower Lower bound of the stability clause. Default is 0.8,
#'     used if uncertaintyCap = TRUE.
#' @param upper Upper bound of the stability clause. Default is 1.2,
#'     used if uncertaintyCap = TRUE.
#' @param amtint Assessment interval. Default is 1, which indicates
#'     annual assessments.
#' @param npriorSD Standard deviation of logn prior (Default: 2). If
#'     NA, the logn prior is removed
#' @param getFit Logical; if TRUE the fitted results list with
#'     adjusted fsihing mortality value is returned. Default is FALSE.
#' @return A list with estimated TAC based on harvest control rule
#'     settings or the fitted rep list with adjusted fishing mortality
#'     values if getFit = TRUE and a logical value indicating if the
#'     stability clause was hit or not (if in use).
#' @export
#' @examples
#' data(pol)
#' rep <- fit.spict(pol$albacore)
#' get.TAC(rep)
get.TAC  <- function(repin, reps = 1,
                     fractileC = 0.5,
                     fractileFFmsy=0.5,
                     fractileBBmsy=0.5,
                     pa=0,
                     prob=0.95,
                     bbmsyfrac=0.3,       
                     stabilityClause=FALSE,
                     lower=0.8,
                     upper=1.2,
                     amtint = 1,
                     npriorSD = 2,
                     getFit = FALSE){
    ## create inp list (note: put in function)
    inp <- list()    
    inp$dteuler <- repin$inp$dteuler
    inp$timeC <- repin$inp$timeC
    inp$obsC <- repin$inp$obsC
    inp$timeI <- repin$inp$timeI
    inp$obsI <- repin$inp$obsI
    inp$timepredc <- repin$inp$timepredc
    inp$timepredi <- inp$timepredc + amtint
    inp$do.sd.report <- TRUE
    inp$getReportCovariance <- FALSE
    inp$MSEmode <- TRUE
    ## check inp
    inp <- check.inp(inp)
    ## stronger prior
    if(!is.null(npriorSD) && is.finite(npriorSD)){
        inp$priors$logn <- c(log(2),npriorSD,1)
    }else{
        inp$priors$logn <- c(0,0,0)
    }
    ## fit spict
    rep <- try(spict::fit.spict(inp),silent=TRUE)
    ## stop if not converged
    if(is.null(rep) || is(rep, "try-error") || rep$opt$convergence != 0 ||
       any(is.infinite(rep$sd))) return(list(TAC=rep(NA, reps),hitSC=FALSE))
    ## get quantities
    logFpFmsy <- get.par("logFpFmsy", rep)
    logBpBmsy <- get.par("logBpBmsy",rep)
    Fmsy <- get.par('logFmsy', rep, exp=TRUE)[2]
    Flast <- get.par('logFl', rep, exp=TRUE)[2]
    ## second non-convergence stop
    if(any(is.null(c(logFpFmsy[2],logBpBmsy[2],Flast,Fmsy))) ||
       !all(is.finite(c(logFpFmsy[2],logBpBmsy[2],Flast,Fmsy))))
        return(list(TAC=rep(NA, reps),hitSC=FALSE))
    ## F multiplication factor based on uncertainty in F/Fmsy. Default = median        
    fi <- 1-fractileFFmsy
    fm <- exp(qnorm(fi, logFpFmsy[2], logFpFmsy[4]))
    fm5 <- exp(qnorm(0.5, logFpFmsy[2], logFpFmsy[4]))
    ## F multiplication factor based on uncertainty in B/Bmsy. Default = median
    bi <- 2 * exp(qnorm(fractileBBmsy, logBpBmsy[2], logBpBmsy[4]))
    fmult <- fm5 / fm * min(1, bi)         
    fabs <- (fmult + 1e-6) * Fmsy / Flast
    ## precautionary approach
    if(pa == 1){
        repPA <- rep
        inpPA <- make.ffacvec(repPA$inp, fabs)
        repPA$obj$env$data$ffacvec <- inpPA$ffacvec
        repPA$obj$env$data$MSEmode <- 1
        repPA$obj$retape()
        repPA$obj$fn(rep$opt$par)
        sdr <- try(sdreport(repPA$obj),silent=TRUE)
        ## stop if not converged
        if(is.null(sdr) || is(sdr, "try-error")) return(list(TAC=rep(NA, reps),hitSC=FALSE))
        ## get quantities
        logBpBmsyPA <- get.par("logBpBmsy",sdr)
        ll <- qnorm(1-prob,logBpBmsyPA[2],logBpBmsyPA[4])
        bbmsyQ5 <- exp(ll)
        ## stop if not finite
        if(is.null(bbmsyQ5) || !is.finite(bbmsyQ5)) return(list(TAC=rep(NA, reps),hitSC=FALSE))
        ## check if precautionary
        if((bbmsyQ5 - bbmsyfrac) < -1e-3){
            tmp <- try(spict:::get.PAffac(rep, bbmsyfrac=bbmsyfrac, prob=prob))
            if(is.null(tmp) ||
               is(tmp, "try-error") || !is.finite(tmp)) return(list(TAC=rep(NA, reps),hitSC=FALSE))
            ## debugging:
            ## if(tmp > fabs) print(paste0("ffacpa",round(tmp,2)," > ffacmsy",round(fabs,2)))
            fabs <- tmp
            fmult <- fabs * Flast / Fmsy
        }
    }
    if(is.null(fmult) || !is.finite(fmult)) return(list(TAC=rep(NA, reps),hitSC=FALSE))
    ## Stability clause
    if(stabilityClause){
        fmult <- spict:::stabilityClause(fmult, lower, upper)
        if(any(fmult < lower) || any(fmult > upper)) hitSC <- TRUE else hitSC <- FALSE
    }else hitSC <- FALSE
    fabs <- fmult * Fmsy / Flast            
    ## predict catch with fabs
    if(is.null(fabs) || !is.finite(fabs)) return(list(TAC=rep(NA, reps),hitSC=FALSE))    
    TACi <- spict:::get.TACi(rep, fabs, fractileC)
    ## hack for DLMtool
    TAC <- rep(TACi, reps)
    ## get fitted object
    if(getFit){
        inpt <- make.ffacvec(rep$inp, fabs)
        inpt$MSEmode <- FALSE
        fit <- try(fit.spict(inpt),silent=TRUE)
        if(!is(fit,"try-error")) return(fit)
    }
    return(list(TAC=TAC, hitSC=hitSC))
}
