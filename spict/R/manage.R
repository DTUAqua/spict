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
#' @references{ICES. 2017. Report of the Workshop on the Development of the ICES approach to providing MSY advice for category 3 and 4 stocks (WKMSYCat34), 6-10 March 2017, Copenhagen, Denmark. ICES CM 2017/ ACOM:47. 53 pp.}
#' @details Scenarios that are currently implemented include:
#' \itemize{
#'   \item{"1"}{ Keep the catch of the current year (i.e. the last observed catch).}
#'   \item{"2"}{ Keep the F of the current year.}
#'   \item{"3"}{ Fish at Fmsy i.e. F=Fmsy.}
#'   \item{"4"}{ No fishing, reduce to 1\% of current F.}
#'   \item{"5"}{ Reduce F by X\%. Default X = 25.}
#'   \item{"6"}{ Increase F by X\%. Default X = 25.}
#'   \item{"7"}{ Use ICES MSY advice rule.}
#'
#' Scenario 7 implements the ICES MSY advice rule for stocks that are assessed using spict (ICES 2017). MSY B_{trigger} is set equal to B_{MSY} / 2. Then fishing mortality in the short forecast is calculated as:
#'
#' F(y+1) =  F(y) * min{ 1, median[B(y+1) / MSY B_{trigger}] } / median[F(y)/F_{MSY}]
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
        scenarios <- 1:7
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
            repman[[1]] <- take.c(catch, inpin, repin, dbg=dbg, catchList=catchList, sdfac = 1)
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
        if (7 %in% scenarios){
          ## F is equal to Fmsy if B>MSYBtrigger. F is reduced linearly to zero if B<MSYBtrigger
          Fmsy <- get.par('logFmsy', repin, exp=TRUE)[2]
          Bmsy <- get.par('logBmsy', repin, exp=TRUE)[2]
          MSYBtrig <- Bmsy / 2
          Flast <- get.par('logFp', repin, exp=TRUE)[2]
          Blast <- get.par('logBp', repin, exp=TRUE)[2]
          Bcorrection <- min(1, Blast / MSYBtrig)
          fac7 <- Bcorrection * (Fmsy / Flast)
          repman[[7]] <- prop.F(fac7, inpin, repin, maninds, dbg=dbg)

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
        inpt$dtc <- c(inpt$dtc, catchList$dtc)
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
                    '4. No fishing', '5. Reduce F 25%', '6. Increase F 25%', '7. MSY advice rule')[scenarios]

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
#' @title Predict the catch of the prediction interval specified in inp
#' @param rep Result list as output from fit.spict().
#' @param fmsyfac Projection are made using F = fmsyfac * Fmsy.
#' @param get.sd Get uncertainty of the predicted catch.
#' @param exp If TRUE report exp of log predicted catch.
#' @param dbg Debug flag, dbg=1 some output, dbg=2 more ourput.
#' @return A vector containing predicted catch (possibly with uncertainty).
#' @export
pred.catch <- function(repin, fmsyfac=1, get.sd=FALSE, exp=FALSE, dbg=0){
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
#' @title Estimate the Total Allowable Catch (TAC)
#' @param rep Result list from fit.spict().
#' @param fractileList List defining the fractiles of the 3 distributions of
#'     "catch", "bbmsy", and "ffmsy" (see details for more information). By
#'     default median is used for all 3 quantities.
#' @param breakpoint_bbmsy Breakpoint in terms of \eqn{B/B_{MSY}} for the
#'     hockey-stick HCR (see details for more information). By default (0) no
#'     breakpoint is assumed.
#' @param paList List defining an optional precautionary buffer by the fraction
#'     of "bbmsy" (default = 0, i.e. deactivating the PA buffer) and "prob" the
#'     risk aversion probability (default = 0.95).
#' @param catch_pred Catch during assessment year (corresponding to argument
#'     \code{catch} in \code{\link{take.c}}), e.g. last year's TAC (default:
#'     \code{NULL}; see details for more information).
#' @param sdfac Factor for the multiplication of the standard deveiation of the
#'     catch during the assessment year (\code{stdevfacC}; default = 1; see
#'     \code{\link{take.c}}).
#' @param getFit Logical; if \code{TRUE} the function returns the fitted
#'     'spictcls' object with respective HCR (\code{FALSE} by default).
#'
#' @details The combination of the arguments in the "fractileList",
#'     "breakpoint_bbmsy", and "paList" allow a number of different harvest
#'     control rules (HCRs):
#' \itemize{
#'
#' \item{Fishing at F_{MSY} (or any fractile of it in terms of \eqn{F/F_{MSY}}
#' or catch) if \code{breakpoint_bbmsy == 0} and \code{paList$bbmsy == 0}.}
#'
#' \item{MSY hockey-stick rule: Fishing at F_{MSY} above a certain fraction of
#'     B_{MSY} (\code{breakpoint_bbmsy}) and below the fraction of B_{MSY}
#'     fishing is reduced to 0 linearly as suggested in ICES (2017) if
#'     \code{breakpoint_bbmsy != 0} and \code{paList$bbmsy = 0}.}
#'
#' \item{MSY (hockey-stick) rule with additional precautionary buffer: As long
#'     as the probability of the predicted biomass relative to a reference
#'     biomass (e.g. 0.3 B_{MSY}, defined by \code{paList$bbmsy}) is smaller or
#'     equal to a specified risk aversion probability (e.g. 95%, defined by
#'     \code{paList$prob}), fishing at F_{MSY} or according to hockey-stick rule
#'     (if \code{breakpoint !=0 }) otherwise reduce fishing effort to meet
#'     specified risk aversion probability (\code{paList$prob}) as introduced in
#'     ICES (2018).}
#'
#' \item{By ICES (2019) recommended MSY hockey-stick rule with 35th percentiles:
#'     Fishing at 35th percentile of F_{MSY} above the 35th percentile of 0.5
#'     \eqn{B/B_{MSY}} (\code{breakpoint_bbmsy = 0.5}) and 35th percentile of
#'     linearly reduced F_{MSY} below the 35th percentile of 0.5
#'     \eqn{B/B_{MSY}}. TAC corresponds to 35th percentile of predicted catch.
#'     Rule is applied with \code{fractileList = list(catch=0.35, bbmsy=0.35,
#'     ffmsy=0.35)}, \code{breakpoint_bbmsy = 0.5}, and \code{paList =
#'     list(bbmsy = 0, prob = 0.95)}.}}
#'
#' More information about the function arguments controlling the HCRs. The
#' arguments of the "fractileList" are: \itemize{
#'
#' \item{catch - Fractile of the predicted catch distribution. Default: 0.5.}
#'
#' \item{bbmsy - Fractile of the \eqn{B/B_{MSY}} distribution. Default: 0.5.}
#'
#' \item{ffmsy - Fractile of the \eqn{F/F_{MSY}} distribution. Default: 0.5.}}
#'
#' The argument "breakpoint_bbmsy" allows to define the MSY hockey-stick rule,
#' which reduces fishing linearly if the biomass is below specified reference
#' level as specified here relative to B_{MSY}. Theoretically, any value below 1
#' is meaningful, but ICES (2017 to 2019) recommend 50% of B_{MSY}
#' (\code{breakpoint_bbmsy = 0.5}). The argument list "paList" includes:
#' \itemize{
#'
#' \item{bbmsy - Reference level for the evaluation of the predicted biomass
#'   defined as fraction of \eqn{B/B_{MSY}}. By default \code{paList$bbmsy == 0}
#'   the PA buffer is not used. Theoretically, any value smaller than 1 is
#'   meaningful, but an ICES recommended value would be 30% \code{paList$bbmsy
#'   = 0.3} (ICES, 2018).}
#'
#' \item{prob - Risk aversion probability of the predicted biomass relative to
#'   specified reference level (\code{paList$bbmsy}) for all rules with PA
#'   buffer (\code{paList$bbmsy != 0}). Default: 0.95 as recommended by ICES
#'   (2018).}}
#'
#' Dependent on the start of the management period (e.g. advice year), there
#' might be a time lag between the last observation and the start of the
#' management period. If this is the case, an assumption about the intermediate
#' time period (e.g. assessment year) has to be made. Either the fishing
#' mortality is extrapolated for the intermediate time period (\code{catch_pred
#' = NULL}; default), or the argument \code{catch_pred} can be used to set the
#' catch in that period. The argument \code{sdfac} allows to adjust the standard
#' deviation of the catch in the intermediate time period.
#'
#' @return A list with the TAC and management relevant quantities; if
#'     \code{getFit} is \code{TRUE} the fitted object with the respective HCR is
#'     returned.
#'
#' @references
#' ICES. 2017. Report of the Workshop on the Development of the ICES
#' approach to providing MSY advice for category 3 and 4 stocks
#' (WKMSYCat34), 6-10 March 2017, Copenhagen, Denmark. ICES CM 2017/
#' ACOM:47. 53 pp.
#'
#' ICES. 2018. Report of the Eighth Workshop on the Development of
#' Quantitative Assessment Methodologies based on LIFE-history traits,
#' exploitation characteristics, and other relevant parameters for
#' data-limited stocks (WKLIFE VIII), 8-12 October 2018, Lisbon,
#' Portugal. ICES CM 2018/ACOM:40. 172 pp.
#'
#' ICES 2019. Report of the Ninth Workshop on the Development of
#' Quantitative Assessment Methodologies based on LIFE-history traits,
#' exploitation characteristics, and other relevant parameters for
#' data-limited stocks (WKLIFE IX), 30 September-4 October 2019,
#' Lisbon, Portugal.
#'
#' @export
#' @examples
#' rep <- fit.spict(pol$albacore)
#'
#' ## Fishing at Fmsy
#' get.TAC(rep)
#'
#' ## MSY hockey-stick rule
#' get.TAC(rep, breakpoint_bbmsy = 0.5)
#'
#' ## ICES (2019) recommended HCR
#' get.TAC(rep, fractileList = list(catch=0.35, bbmsy=0.35, ffmsy=0.35), breakpoint_bbmsy=0.5)
#'
get.TAC <- function(rep,
                    fractileList = list(catch = 0.5, bbmsy = 0.5, ffmsy = 0.5),
                    breakpoint_bbmsy = 0,
                    paList = list(bbmsy = 0, prob = 0.95),
                    catch_pred = NULL,
                    sdfac = 1,
                    getFit = FALSE){
    reppa <- repin <- rep
    inpin <- repin$inp

    ## Default lists (if user mis-specifies list or changes single element only)
    ## fractile list
    if(is.list(fractileList)) stop("Please provide 'fractileList' with the arguments: 'catch', 'bbmsy', and 'ffmsy'!")
    default_fractileList = list(catch = 0.5, bbmsy = 0.5, ffmsy = 0.5)
    fList <- default_fractileList[which(!names(default_fractileList) %in% names(fractileList))]
    fList <- c(fList,fractileList)
    ## PA list
    if(is.list(paList)) stop("Please provide 'paList' with the arguments: 'bbmsy' and 'prob'!")
    default_paList = list(bbmsy = 0, prob = 0.95)
    pList <- default_paList[which(!names(default_paList) %in% names(paList))]
    pList <- c(pList,paList)

    ## option for assessment year (intermediate year between last data and advice year)
    inttime <- inpin$dtpredcinds[1] - min(inpin$indpred)
    if(inttime > 0 && !is.null(catch_pred) && !is.na(catch_pred) && is.numeric(catch_pred)){
        ## make catchList for projected years (timesteps) before manstart
        catchList <- list()
        catchList$timeC <- inpin$time[min(inpin$indpred)]
        catchList$obsC <- catch_pred
        catchList$stdevfacC <- sdfac
        catchList$dtc <- (inpin$dtpredcinds[1] - min(inpin$indpred)) * inpin$dteuler
        inpin$reportmode <- 1
        repin <- take.c(catch_pred, inpin, repin, catchList = catchList)
        inpin <- repin$inp
    }

    ## quantities
    fmanstart <- get.par('logFm', repin, exp=TRUE)[2]
    fmsy <- get.par('logFmsy', repin, exp=TRUE)[2]
    bmsy <- get.par('logBmsy', repin, exp=TRUE)[2]
    logFpFmsy <- get.par("logFpFmsynotS", repin)
    logBpBmsy <- get.par("logBpBmsy", repin)
    logFmFmsy <- get.par("logFmFmsynotS", repin)
    logBmBmsy <- get.par("logBmBmsy", repin)

    ## HCRs
    ## FFmsy component
    fi <- 1 - fList$ffmsy
    fmfmsyi <- exp(qnorm(fi, logFmFmsy[2], logFmFmsy[4]))
    fmfmsy5 <- exp(qnorm(0.5, logFmFmsy[2], logFmFmsy[4]))
    fred <- fmfmsy5 / fmfmsyi
    ## BBmsy component (hockey stick HCR)
    if(!is.na(breakpoint_bbmsy) && is.numeric(breakpoint_bbmsy) && breakpoint_bbmsy == 0){
        bmbmsyi <- 1/breakpoint_bbmsy * exp(qnorm(fList$bbmsy, logBmBmsy[2], logBmBmsy[4]))
        fred <- fred * min(1, bmbmsyi)
    }
    ## F reduction factor
    ffac <- (fred + 1e-8) * fmsy / fmanstart
    ## PA component
    if(!is.na(pList$bbmsy) && is.numeric(pList$bbmsy) && pList$bbmsy == 0){
        inppa <- make.ffacvec(inpin, ffac)
        reppa$obj$env$data$ffacvec <- inppa$ffacvec
        reppa$obj$env$data$reportmode <- 1
        reppa$obj$retape()
        reppa$obj$fn(repin$opt$par)
        sdr <- try(sdreport(reppa$obj), silent=TRUE)
        if(is(sdr,"try-error")) stop("Model variances could not be estimated with the PA.")
        logBpBmsyPA <- get.par("logBpBmsy", sdr)
        probi <- 1 - pList$prob
        bpbmsyiPA <- exp(qnorm(probi, logBpBmsyPA[2], logBpBmsyPA[4]))
        if((bpbmsyiPA - pList$bbmsy) < -1e-3){
            ffac <- try(get.ffac(reppa, bfrac=pList$bbmsy, prob=pList$prob,
                                 quant="logBpBmsy",
                                 reportmode = 2), silent=TRUE)
            if(is(ffac,"try-error")) stop("F multiplication factor could not be estimated with the PA.")
        }
    }
    ## Catch component
    tac <- try(calc.tac(repin, ffac, fList$catch), silent=TRUE)
    if(is(tac,"try-error")) stop("TAC could not be estimated.")

    ## get fitted object
    if(getFit){
        inpt <- make.ffacvec(repin$inp, ffac)
        repin$obj$env$data$ffacvec <- inpt$ffacvec
        repin$obj$env$data$reportmode <- 0
        repin$obj$retape()
        repin$obj$fn(repin$opt$par)
        fit <- try(sdreport(repin$obj),silent=TRUE)
        if(!is(fit,"try-error")) return(fit) else stop("The model could not be fitted.")
    }

    ## results
    reslist <- list(TAC = tac,
                    maninterval = inpin$maninterval,
                    maneval = inpin$maneval,
                    ffac = ffac,
                    fractileList = fList,
                    breakpoint_bbmsy = breakpoint_bbmsy,
                    paList = pList,
                    catch_pred = catch_pred,
                    sdfac = sdfac)
    ## return
    return(reslist)
}
