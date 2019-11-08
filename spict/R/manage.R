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
#' @param hcr SPiCT specific harvest control rule (HCR). Possible
#'     rules are: \code{"MSY"}, \code{"MSY-HS"}, or \code{"MSY-HS-PA"}
#'     (see details for more information).  below.
#' @param args A list with specific arguments for the respective HCR
#'     (see details for possible arguments).
#' @param curc Optional; catch during assessment year, e.g. last
#'     year's TAC (default: \code{NULL}; see details for more
#'     information).
#' @param sdfac Factor for the multiplication of the standard
#'     deveiation of the catch during the assessment year
#'     (\code{stdevfacC}; default = 1e-3).
#' @param getFit Logical; if \code{TRUE} the function returns the
#'     fitted 'spictcls' object with respective HCR (\code{FALSE} by
#'     default).
#'
#' @details The possible harvest control rules (argument \code{hcr}) are:
#' \itemize{
#'   \item{MSY- MSY rule: Fishing at F_{MSY}.}
#'   \item{MSY-HS - MSY Hockey-Stick rule: Above 0.5B_{MSY} fishing at F_{MSY}, below 0.5B_{MSY} fishing linearly reduced to 0 as suggested in ICES (2017).}
#'    \item{MSY-ICES - MSY ICES rule: Same as the \code{MSY-HS} rule but using the 35th percentiles for predicted catch, \eqn{B/B_{MSY}}, and \eqn{F/F_{MSY}} following ICES (2019).}
#'   \item{MSY-PA - MSY rule with additional precautionary buffer: As long as the probability of the predicted biomass relative to a reference biomass (e.g. 0.3B_{MSY}) is above a specified level (e.g. 5%), rule corresponds to \code{MSY}, otherwise reduce F to meet specified probability as introduced in ICES (2018).}
#'   \item{MSY-HS-PA - MSY Hockey-stick rule with additional precautionary buffer: As long as the probability of the predicted biomass relative to a reference biomass (e.g. 0.3B_{MSY}) is above a specified level (e.g. 5%), rule corresponds to \code{MSY-HS}, otherwise reduce F to meet specified probability as introduced in ICES (2018).}
#' }
#'
#' The possible arguments of the \code{args} list are:
#' \itemize{
#'   \item{fracc - Fractile of the predicted catch distribution. Default: 0.5.}
#'   \item{fracf - Fractile of the \eqn{F/F_{MSY}} distribution. Default: 0.5.}
#'   \item{fracb - Fractile of the \eqn{B/B_{MSY}} distribution. Default: 0.5.}
#'   \item{bfrac - Reference level for the evaluation of the predicted biomass defined as fraction of \eqn{B/B_{MSY}}. Default: 0.3.}
#'   \item{babs - Reference level for the evaluation of the predicted biomass defined in absolute terms. Default: \code{NA}, which means that the reference level is defined in relative terms (\code{bfrac}).}
#'   \item{prop - Risk aversion probability level of the predicted biomass relative to specified reference level for the \code{MSY-HS-PA} rule. Default: 0.95, which corresponds to an accepted risk of 5% (1-0.95).}
#'   \item{reportmode - Determining which objects will be adreported. Default: 2 = only objects relevant for get.TAC are adreported.}
#' }
#'
#' Dependent on the start of the management period (e.g. advice year),
#' there might be a time lag between the last observation and the
#' start of the management period. If this is the case, an assumption
#' about the intermediate time period (e.g. assessment year) has to be
#' made. Either the fishing mortality is extrapolated for the
#' intermediate time period (\code{curc = NULL}; default), or the
#' argument \code{curc} can be used to set the catch in that period.
#'
#' @return A list with absolute and relative reference levels and
#'     states and estimated TAC; if \code{getFit} is \code{TRUE} the
#'     fitted object with the respective HCR is returned.
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
#' get.TAC(rep)
get.TAC <- function(rep,
                    hcr = c("MSY","MSY-HS","MSY-HS-frac","MSY-PA","MSY-HS-PA")[1],
                    args = NULL,
                    curc = NULL,
                    sdfac = 1e-3,
                    getFit = FALSE){
    repin <- rep
    inpin <- repin$inp

    ## all hcrs
    allhcrs <- c("MSY","MSY-HS","MSY-PA","MSY-HS-PA")

    ## elements of args
    if(is.null(args)) args <- list()
    args$fracc <- if(!"fracc" %in% names(args)) 0.5 else args$fracc
    args$fracf <- if(!"fracf" %in% names(args)) 0.5 else args$fracf
    args$fracb <- if(!"fracb" %in% names(args)) 0.5 else args$fracb
    args$bfrac <- if(!"bfrac" %in% names(args)) 0.3 else args$bfrac
    args$babs <- if(!"babs" %in% names(args)) NA else args$babs
    args$prop <- if(!"prop" %in% names(args)) 0.95 else args$prop
    args$reportmode <- if(!"reportmode" %in% names(args)) 1 else args$reportmode

    ## option for assessment year (intermediate year between last data and advice year)
    inttime <- inpin$dtpredcinds[1] - min(inpin$indpred)
    if(inttime > 0 && !is.null(curc) && !is.na(curc) && is.numeric(curc)){
        ## make catchList for projected years (timesteps) before manstart
        catchList <- list()
        catchList$timeC <- inpin$time[min(inpin$indpred)]
        catchList$obsC <- curc
        catchList$stdevfacC <- sdfac
        catchList$dtc <- (inpin$dtpredcinds[1] - min(inpin$indpred)) * inpin$dteuler
        inpin$reportmode <- 1
        repin <- take.c(curc, inpin, repin, catchList = catchList)
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

    ## derived HCRs
    if(hcr == "MSY-ICES"){
        args$fracf <- 0.35
        args$fracb <- 0.35
        args$fracc <- 0.35
        hcr2 <- "MSY-HS"
    }else if(hcr == "MSY-HS-PA"){
        hcr2 <- "MSY-PA"
    }else{
        hcr2 <- hcr
    }

    ## rules
    switch(hcr2,
           "MSY" = {
               ffac <- fmsy / fmanstart
               tac <- calc.tac(repin, ffac, args$fracc)
           },
           "MSY-HS" = {
               fi <- 1 - args$fracf
               fmfmsyi <- exp(qnorm(fi, logFmFmsy[2], logFmFmsy[4]))
               fmfmsy5 <- exp(qnorm(0.5, logFmFmsy[2], logFmFmsy[4]))
               bmbmsyi <- 2 * exp(qnorm(args$fracb, logBmBmsy[2], logBmBmsy[4]))
               fred <- fmfmsy5 / fmfmsyi * min(1, bmbmsyi)
               ffac <- (fred + 1e-8) * fmsy / fmanstart
               tac <- calc.tac(repin, ffac, args$fracc)
           },
           "MSY-PA" = {
               quant <- "logBpBmsy"
               repcop <- repin
               bmbmsyi <- 2 * exp(qnorm(args$fracb, logBmBmsy[2], logBmBmsy[4]))
               ffac <- fmsy / fms * min(1, bmbmsyi)
               if(hcr == "MSY-HS-PA"){
                   fi <- 1 - args$fracf
                   fmfmsyi <- exp(qnorm(fi, logFmFmsy[2], logFmFmsy[4]))
                   fmfmsy5 <- exp(qnorm(0.5, logFmFmsy[2], logFmFmsy[4]))
                   fred <- fmfmsy5 / fmfmsyi * min(1, bmbmsyi)
                   ffac <- (fred + 1e-8) * fmsy / fmanstart
               }
               inpcop <- make.ffacvec(inpin, ffac)
               repcop$obj$env$data$ffacvec <- inpcop$ffacvec
               repcop$obj$env$data$reportmode <- 2
               repcop$obj$retape()
               repcop$obj$fn(repin$opt$par)
               sdr <- try(sdreport(repcop$obj), silent=TRUE)
               logBpBmsycop <- get.par(quant, sdr)
               propi <- 1 - args$prob
               bpbmsyi <- exp(qnorm(propi, logBpBmsycop[2], logBpBmsycop[4]))
               bfrac <- args$bfrac
               if(!is.na(args$babs) && is.numeric(args$babs)){
                   bfrac <- args$babs / bmsy
               }
               if((bpbmsyi - bfrac) < -1e-3){
                   ffac <- try(get.ffac(repcop, bfrac=bfrac, prob=args$prob,
                                        quant=quant, reportmode = args$reportmode), silent = TRUE)

               }
               tac <- calc.tac(repin, ffac, args$fracc)
           },
           stop(paste0("The specified 'hcr' is not known. Please choose between: ",
                       paste0(allhcrs, collapse = "; "))))

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
    reslist <- list(hcr = hcr,
                    manint = inpin$manint,
                    maneval = inpin$maneval,
                    TAC = tac,
                    ffac = ffac,
                    curc = curc,
                    sdfac = sdfac,
                    args = args)
    ## return
    return(reslist)
}
