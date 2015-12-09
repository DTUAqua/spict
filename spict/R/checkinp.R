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


#' @name check.inp
#' @title Check list of input variables
#' @details Fills in defalut values if missing.
#'
#' Required inputs:
#' 
#' \itemize{
#'  \item{"inp$obsC"}{ Vector of catch observations.}
#'  \item{"inp$obsI"}{ List containing vectors of index observations.}
#' }
#' 
#' Optional inputs:
#' 
#' - Data
#' \itemize{
#'  \item{"inp$timeC"}{ Vector of catch times. Default: even time steps starting at 1.}
#'  \item{"inp$timeI"}{ List containing vectors of index times. Default: even time steps starting at 1.}
#'  \item{"inp$dtc"}{ Time interval for catches, e.g. for annual catches inp$dtc=1, for quarterly catches inp$dtc=0.25. Can be given as a scalar, which is then used for all catch observations. Can also be given as a vector specifying the catch interval of each catch observation. Default: min(diff(inp$timeC)). }
#'  \item{"inp$nseasons"}{ Number of within-year seasons in data. If inp$nseasons > 1 then a seasonal pattern is used in F. Valid values of inp$nseasons are 1, 2 or 4. Default: number of unique within-year time points present in data.}
#' }
#' 
#' - Parameters
#' 
#' \itemize{
#'  \item{"inp$ini$logn"}{ Pella-Tomlinson exponent determining shape of production function. Default: log(2) corresponding to the Schaefer formulation.}
#'  \item{"inp$ini$logm"}{ Initial value for logm (log maximum sustainable yield). Default: log(mean(catch)).}
#'  \item{"inp$ini$logK"}{ Initial value for logK (log carrying capacity). Default: log(4*max(catch)).}
#'  \item{"inp$ini$logq"}{ Initial value for logq (log catchability of index). Default: log(max(index)/K).}
#'  \item{"inp$ini$logsdb"}{ Initial value for logsdb (log standard deviation of biomass process). Default: log(0.2).}
#'  \item{"inp$ini$logsdf"}{ Initial value for logsdf (log standard deviation of fishing mortality process). Default: log(0.2).}
#'  \item{"inp$ini$logsdi"}{ Initial value for logsdi (log standard deviation of index observation error). Default: log(0.2).}
#'  \item{"inp$ini$logsdc"}{ Initial value for logsdc (log standard deviation of catch observation error). Default: log(0.2).}
#'  \item{"inp$ini$phi"}{ Vector for cyclic B spline representing within-year seasonal variation. Default: rep(1, inp$nseasons).}
#'  \item{"inp$ini$logsdu"}{ Initial value for logsdu (log standard deviation of log U, the state of the coupled SDE representation of seasonality). Default: log(0.1).}
#'  \item{"inp$ini$loglambda"}{ Initial value for loglambda (log damping parameter of the coupled SDE representation of seasonality). Default: log(0.1).}
#' }
#'
#' - Unobserved states estimated as random effects
#' 
#' \itemize{
#'   \item{"inp$ini$logF"}{ Log fishing mortality. Default: log(0.2*r), with r derived from m and K.}
#'   \item{"inp$ini$logB"}{ Log biomass. Default: log(0.5*K).}
#'   \item{"inp$ini$logU"}{ Log U, the state of the coupled SDE representation of seasonality. Default: log(1).}
#' }
#' 
#' - Priors
#' 
#' Priors on model parameters are assumed Gaussian and specified in a vector of length 2: c(log(mean), stdev in log domain, useflag [optional]). NOTE: if specifying a prior for logB, then a 4th element is required specifying the year the prior should be applied.
#' log(mean): log of the mean of the prior distribution.
#' stdev in log: standard deviation of the prior distribution in log domain.
#' useflag: if 1 then the prior is used, if 0 it is not used. Default is 0.
#' To list parameters to which priors can be applied run list.possible.priors().
#' Example: intrinsic growth rate of 0.8
#'  inp$priors$logr <- c(log(0.8), 0.1)
#'  inp$priors$logr <- c(log(0.8), 0.1, 1) # This includes the optional useflag
#' Example: Biomass prior of 200 in 1985
#'  inp$priors$logB <- c(log(200), 0.2, 1985)
#'  inp$priors$logB <- c(log(200), 0.2, 1, 1985) # This includes the optional useflag
#' 
#' - Settings/Options/Preferences
#' 
#' \itemize{
#'  \item{"inp$dtpredc"}{ Length of catch prediction interval in years. Default: max(inp$dtc). Should be 1 to get annual predictions and 0.25 for quarterly predictions.}
#'  \item{"inp$timepredc"}{ Predict catches in interval lengths given by $dtpredc until this time. Default: Time of last observation. Example: inp$timepredc <- 2012}
#'  \item{"inp$timepredi"}{ Predict index until this time. Default: Time of last observation. Example: inp$timepredi <- 2012}
#'  \item{"inp$do.sd.report"}{ Flag indicating whether SD report (uncertainty of derived quantities) should be calculated. For small values of inp$dteuler this may require a lot of memory. Default: TRUE.}
#'  \item{"inp$reportall"}{ Flag indicating whether quantities derived from state vectors (e.g. B/Bmsy, F/Fmsy etc.) should be calculated by SD report. For small values of inp$dteuler (< 1/32) reporting all may have to be set to FALSE for sdreport to run. Additionally, if only reference points of parameter estimates are of interest one can set to FALSE to gain a speed-up. Default: TRUE.}
#' \item{"inp$robflagc"}{ Flag indicating whether robust estimation should be used for catches (either 0 or 1). Default: 0.}
#'  \item{"inp$robflagi"}{ Flag indicating whether robust estimation should be used for indices (either 0 or 1). Default: 0.}
#'  \item{"inp$ffac"}{ Management scenario represented by a factor to multiply F with when calculating the F of the next time step. ffac=0.8 means a 20\% reduction in F over the next year. The factor is only used when predicting beyond the data set. Default: 1 (0\% reduction).}
#'  \item{"inp$dteuler"}{ Length of Euler time step in years. Default: 1/16 year.}
#'  \item{"inp$phases"}{ Phases can be used to fix/free parameters and estimate in different stages or phases. To fix e.g. logr at inp$ini$logr set inp$phases$logr <- -1. To free logalpha and estimate in phase 1 set inp$phases$logalpha <- 1.}
#'  \item{"inp$osar.method"}{ Method to use in TMB's oneStepPredict function. Valid methods include: "oneStepGaussianOffMode", "fullGaussian", "oneStepGeneric", "oneStepGaussian", "cdf". See TMB help for more information. Default: "none" (i.e. don't run this).}
#'  \item{"inp$osar.trace"}{ If TRUE print OSAR calculation progress to screen. Default: FALSE.}
#'  \item{"inp$osar.parallel"}{ If TRUE parallelise OSAR calculation for speed-up. Default: FALSE.}
#'  \item{"inp$catchunit"}{ Specify unit of catches to be used in plotting legends. Default: ''.}
#'  \item{"inp$stdevfacC"}{ Factors to multiply the observation error standard deviation of each individual catch observation. Can be used if some observations are more uncertain than others. Must be same length as observation vector. Default: 1.}
#'  \item{"inp$stdevfacI"}{ Factors to multiply the observation error standard deviation of each individual index observation. Can be used if some observations are more uncertain than others. A list with vectors of same length as observation vectors. Default: 1.}
#'  \item{"inp$mapsdi"}{ Vector of length equal to the number of index series specifying which indices that should use the same sdi. For example: in case of 3 index series use inp$mapsdi <- c(1, 1, 2) to have series 1 and 2 share sdi and have a separate sdi for series 3. Default: 1:nindex, where nindex is number of index series.}
#' }
#' @param inp List of input variables, see details for required variables.
#' @return An updated list of input variables checked for consistency and with defaults added.
#' @examples
#' data(pol)
#' (inp <- check.inp(pol$albacore))
#' @export
check.inp <- function(inp){
    check.ini <- function(parname, inp, min=NULL, max=NULL){
        if(!parname %in% names(inp$ini)) stop('Please specify an initial value for ', parname, '!')
    }
    # -- DATA --
    # Check catches
    if('obsC' %in% names(inp)){
        if(!'timeC' %in% names(inp)){
            if(!'nseasons' %in% names(inp)){
                inp$timeC <- 0:(length(inp$obsC)-1)
            } else {
                to <- (length(inp$obsC)-1)/inp$nseasons
                inp$timeC <- seq(0, to, length=length(inp$obsC))
            }
        }
        if(any(diff(inp$timeC)<=0)) stop('Catch times are not strictly increasing!')
        if(length(inp$obsC) != length(inp$timeC)) stop('Time and observation vector do not match in length for catch series.')
        if(!'stdevfacC' %in% names(inp)) inp$stdevfacC <- rep(1, length(inp$obsC))
        if(length(inp$obsC) != length(inp$stdevfacC)) stop('stdevfac and observation vector do not match in length for catch series.')
        if(sum(inp$stdevfacC <= 0) > 0) stop('Non-positive values entered in stdevfac for catches.')
        neg <- which(inp$obsC<=0 | is.na(inp$obsC))
        if(length(neg)>0){
            inp$obsC <- inp$obsC[-neg]
            inp$timeC <- inp$timeC[-neg]
            inp$stdevfacC <- inp$stdevfacC[-neg]
            cat(paste('Removing zero, negative, and NAs in catch series\n'))
        }
        inp$nobsC <- length(inp$obsC)
        inp$obsidC <- 1:inp$nobsC
    } else {
        stop('No catch observations included. Please include them as a vector in inp$obsC.')
    }

    # Check indices
    if('obsI' %in% names(inp)){
        if(class(inp$obsI)!='list'){
            tmp <- inp$obsI
            inp$obsI <- list()
            inp$obsI[[1]] <- tmp
        }
        inp$nindex <- length(inp$obsI)
        # Time vector
        if(!'timeI' %in% names(inp)){
            inp$timeI <- list()
            if(!'nseasons' %in% names(inp)){
                inp$timeI[[1]] <- 0:(length(inp$obsI[[1]])-1)
            } else {
                to <- (length(inp$obsI[[1]])-1)/inp$nseasons
                inp$timeI[[1]] <- seq(0, to, length=length(inp$obsI[[1]]))
            }
        } else {
            if(class(inp$timeI)!='list'){
                tmp <- inp$timeI
                inp$timeI <- list()
                inp$timeI[[1]] <- tmp
            }
        }
        # Standard deviation factor
        if(!'stdevfacI' %in% names(inp)){
            inp$stdevfacI <- list()
            for(i in 1:inp$nindex) inp$stdevfacI[[i]] <- rep(1, length(inp$obsI[[i]]))
        } else {
            if(class(inp$stdevfacI)!='list'){
                tmp <- inp$stdevfacI
                inp$stdevfacI <- list()
                inp$stdevfacI[[1]] <- tmp
            }
        }
        if(inp$nindex != length(inp$timeI)) stop('length(inp$timeI) is not equal to length(inp$obsI)!')
        if(inp$nindex != length(inp$stdevfacI)) stop('length(inp$stdevfacI) is not equal to length(inp$obsI)!')
        inp$nobsI <- rep(0, inp$nindex)
        for(i in 1:inp$nindex){
            if(length(inp$obsI[[i]]) != length(inp$timeI[[i]])) stop('Time and observation vector do not match in length for index series ', i)
            if(length(inp$obsI[[i]]) != length(inp$stdevfacI[[i]])) stop('stdevfac and observation vector do not match in length for index series ', i)
            if(sum(inp$stdevfacI[[i]] <= 0) > 0) stop('Non-positive values entered in stdevfac for index ', i)
            if(length(inp$obsI[[i]])>0){
                neg <- which(inp$obsI[[i]]<=0 | is.na(inp$obsI[[i]]))
                if(length(neg)>0){
                    inp$obsI[[i]] <- inp$obsI[[i]][-neg]
                    inp$timeI[[i]] <- inp$timeI[[i]][-neg]
                    inp$stdevfacI[[i]] <- inp$stdevfacI[[i]][-neg]
                    cat(paste('Removing zero, negative and NAs in index series',i,'\n'))
                }
            }
            inp$nobsI[i] <- length(inp$obsI[[i]])
            if(i==1){
                inp$obsidI[[i]] <- (1:inp$nobsI[i]) + inp$nobsC
            } else {
                inp$obsidI[[i]] <- (1:inp$nobsI[i]) + tail(inp$obsidI[[i-1]], 1)
            }
        }
    } else {
        stop('No index observations included. Please include them as a list in inp$obsI.')
    }

    # -- MAX MIN RATIO (used in simulation only) --
    inp$maxminratio <- rep(0, inp$nindex)
    names(inp$maxminratio) <- paste0('I', 1:inp$nindex)
    for(i in 1:inp$nindex){
        if(length(inp$obsI[[i]])>0) inp$maxminratio[i] <- max(inp$obsI[[i]])/min(inp$obsI[[i]])
    }

    # -- MODEL OPTIONS --
    if(!"RE" %in% names(inp)) inp$RE <- c('logF', 'logu', 'logB')
    if(!"scriptname" %in% names(inp)) inp$scriptname <- 'spict'
    if(!"onealpha" %in% names(inp)){
        if(!"onesdi" %in% names(inp)){
            inp$onealpha <- FALSE
        } else {
            inp$onealpha <- inp$onesdi
        }
    }
    if(!"onesdi" %in% names(inp)) inp$onesdi <- inp$onealpha
    if(!"mapsdi" %in% names(inp)){
        inp$mapsdi <- 1:inp$nindex
        if("onealpha" %in% names(inp)) if(inp$onealpha) inp$mapsdi <- rep(1, inp$nindex)
    }
    if("mapsdi" %in% names(inp)) inp$nsdi <- length(unique(inp$mapsdi))
    if(!"mapq" %in% names(inp)) inp$mapq <- 1:inp$nindex
    if("mapq" %in% names(inp)) inp$nq <- length(unique(inp$mapq))
    if(!"catchunit" %in% names(inp)) inp$catchunit <- ''
    if(!"reportall" %in% names(inp)) inp$reportall <- TRUE
    if(!"do.sd.report" %in% names(inp)) inp$do.sd.report <- TRUE
    # Simulation options
    if(!"armalistF" %in% names(inp)) inp$armalistF <- list() # Used for simulating arma noise for F instead of white noise.
    # Optimiser options
    if(!"optimiser" %in% names(inp)) inp$optimiser <- 'nlminb'
    if(!"optimiser.control" %in% names(inp)) inp$optimiser.control <- list()
    # OSAR options
    if(!"osar.method" %in% names(inp)) inp$osar.method <- 'none'
    if(!"osar.trace" %in% names(inp))  inp$osar.trace <- FALSE
    if(!"osar.parallel" %in% names(inp)) inp$osar.parallel <- FALSE
    # Season options
    if(!"seasontype" %in% names(inp)) inp$seasontype <- 1
    if(!"omega" %in% names(inp)) inp$omega <- 2*pi # Annual cycle of noisy oscillator
    # Robust options
    if(!"robflagc" %in% names(inp)) inp$robflagc <- 0
    inp$robflagc <- as.numeric(inp$robflagc)
    if(!"robflagi" %in% names(inp)) inp$robflagi <- 0
    inp$robflagi <- as.numeric(inp$robflagi)
    # ASPIC options 
    if(!"aspic" %in% names(inp)) inp$aspic <- list()
    if(!"mode" %in% names(inp$aspic)) inp$aspic$mode <- 'FIT'
    if(!"verbosity" %in% names(inp$aspic)) inp$aspic$verbosity <- '102'
    if(!"nboot" %in% names(inp$aspic)) inp$aspic$nboot <- 1000
    if(!"ciperc" %in% names(inp$aspic)) inp$aspic$ciperc <- 95

    # Options for simple model
    if(!"simple" %in% names(inp)) inp$simple <- 0
    if(inp$simple==1){ # Set parameters for the simple model (catch assumed known, no F process).
        umodtimeC <- unique(inp$timeC%%1)
        if(length(umodtimeC) != 1) stop('When inp$simple = 1, inp$timeC must have a fixed regular time step of 1 year!')
        if(umodtimeC != 0) inp$timeC <- floor(inp$timeC)
        for(i in 1:inp$nindex){
            umodtimeI <- unique(inp$timeI[[i]]%%1)
            if(length(umodtimeI) != 1) stop('When inp$simple = 1, inp$timeI must have a fixed regular time step of 1 year!')
            if(umodtimeI != 0) inp$timeI[[i]] <- floor(inp$timeI[[i]])
        }
        inp$dteuler <- 1
        # Fix parameters
        inp$phases <- list()
        inp$phases$logF <- -1
        inp$phases$logsdf <- -1
        inp$phases$logbeta <- -1
        # Prediction not possible
        inp$dtpredc <- NULL
        inp$timepredc <- max(inp$timeC)
        inp$timepredi <- max(unlist(inp$timeI))
    }
    # MSY type options
    if(!"msytype" %in% names(inp)){
        inp$msytype <- 's'
    } else {
        if(!inp$msytype %in% c('s', 'd')) stop('inp$msytype must be either "s" (stochastic) or "d" (deterministic!')
    }
    #if("phases" %in% names(inp)){
    #    if('logn' %in% names(inp$phases)){
    #        if(inp$phases$logn > -1){
    #            if(!"msytype" %in% names(inp)){
    #                inp$msytype <- 'd'
    #                cat('Using msytype = "d" because logn is estimated. Force SPiCT to use msytype = "s" by manually specifying it.\n')
    #            } else {
                    # Dangerous to use stochastic msys when estimating n because estimate of n could result in invalid reference points.
                    #if(inp$msytype != 'd') cat('Using msytype = "d" because logn is estimated.\n')
    #                if(inp$msytype != 'd') warning('Dangerous to use stochastic msys when estimating n because estimate of n could result in invalid reference points. Check difference between deterministic and stochastic reference points.')
    #            }
    #        }
    #    }
    #}
    # Catch intervals (dtc)
    if(!"dtc" %in% names(inp)){
        dtc <- diff(inp$timeC)
        if(length(dtc)>0){
            inp$dtc <- min(dtc)
            #cat(paste('Catch interval (dtc) not specified. Assuming an interval of:', inp$dtc, 'year.\n'))
        } else {
            inp$dtc <- 1
            cat(paste('Catch interval (dtc) not specified and length of catch time series shorter than 2. Assuming an interval of 1 year.\n'))
        }
    }
    if(length(inp$dtc)==1) inp$dtc <- rep(inp$dtc, inp$nobsC)
    if(length(inp$dtc) != inp$nobsC) stop('Catch interval vector (inp$dtc) does not match catch observation vector (inp$obsC) in length')

    # - Prediction horizons -
    timeobsall <- sort(c(inp$timeC, inp$timeC + inp$dtc, unlist(inp$timeI)))
    # Catch prediction time step (dtpredc)
    if(!"dtpredc" %in% names(inp)){
        if(length(inp$dtc)>0){
            inp$dtpredc <- max(inp$dtc)
        } else {
            inp$dtpredc <- 1
            cat('Assuming a 1 year prediction interval for catch.\n')
        }
    }
    # Time point to predict catches until
    if(!"timepredc" %in% names(inp)){
        inp$timepredc <- max(timeobsall)
    } else {
        if(inp$timepredc < max(inp$timeC)) cat('inp$timepredc:', inp$timepredc, ' must be equal to or later than last catch observation: ', max(inp$timeC), '!')
    }
    # Time point to predict indices until
    if(!"timepredi" %in% names(inp)){
        inp$timepredi <- max(timeobsall)
    } else {
        if(inp$timepredi < max(unlist(inp$timeI))) stop('inp$timepredi must be equal to or later than last index observation!')
    }

    # Numerical Euler discretisation time used by SDE solver
    #timesteps <- diff(sort(timeobs))
    #timesteps <- timesteps[timesteps > 0]
    #if("dteuler" %in% names(inp)) if(inp$dteuler > min(timesteps)){
    #    cat('inp$dteuler is', inp$dteuler, 'while the minimum time step of observations is', min(timesteps), 'inp$dteuler will be changed!')
    #    inp$dteuler <- NULL
    #}
    #if(!"dteuler" %in% names(inp)) inp$dteuler <- 0.5*min(timesteps) # half because often a time step finer than obs is required
    #if(!"dteuler" %in% names(inp)) inp$dteuler <- min(1/16, 0.5*min(timesteps)) # half because often a time step finer than obs is required
    
    # Euler time step
    if(!"dteuler" %in% names(inp)) inp$dteuler <- 1/16
    if("dteuler" %in% names(inp)){
        alloweddteuler <- 1/2^(6:0)
        if(!inp$dteuler %in% alloweddteuler){ # Check if dteuler is among the alloweddteuler
            ind <- cut(inp$dteuler, alloweddteuler, right=FALSE, labels=FALSE)
            if(is.na(ind)){
                if(inp$dteuler > max(alloweddteuler)) inp$dteuler <- max(alloweddteuler)
                if(inp$dteuler < min(alloweddteuler)) inp$dteuler <- min(alloweddteuler)
            } else {
                inp$dteuler <- alloweddteuler[ind]
            }
            cat('The dteuler used is not allowed! using inp$dteuler:', inp$dteuler, '\n')
        }
    }
    # Euler types:
    # hard: time discretisation is equidistant with step length = dteuler. Observations are assigned to intervals
    # soft: time discretisation is equidistant with step length = dteuler, but with time points of observations inserted such that they can be assigned accurately to a time point instead of an interval.
    # Include dtc because a catch observation at time t includes information in the interval [t; t+dtc[
    if(!"eulertype" %in% names(inp)){
        inp$eulertype <- 'hard'
    }
    if("eulertype" %in% names(inp)){
        if(inp$eulertype == 'hard'){
            # Hard Euler discretisation
            time <- seq(floor(min(timeobsall)), max(inp$timepredi, inp$timepredc+inp$dtpredc), by=inp$dteuler)
            inp$time <- time
        }
        if(inp$eulertype == 'soft'){
            # Include times of observations (including dtc)
            time <- seq(ceiling(min(timeobsall)), max(inp$timepredi, inp$timepredc+inp$dtpredc), by=inp$dteuler)
            inp$time <- sort(unique(c(timeobsall, time)))
        }
        if(!inp$eulertype %in% c('soft', 'hard'))
            stop('inp$eulertype must be either "soft" or "hard"!')
    }
    # Calculate time steps
    inp$dt <- c(diff(inp$time), inp$dteuler)
    #inp$dt <- rep(inp$dteuler, inp$ns)
    inp$ns <- length(inp$time)

    # -- DERIVED VARIABLES --
    #timeobsnodtc <- c(inp$timeC, unlist(inp$timeI))
    
    #inp$indlastobs <- which(inp$time == max(c(inp$timeC, unlist(inp$timeI))))
    inp$timerange <- range(timeobsall)
    inp$indlastobs <- cut(max(c(inp$timeC, unlist(inp$timeI))), inp$time, right=FALSE, labels=FALSE)
    inp$indest <- which(inp$time <= inp$timerange[2])
    inp$indpred <- which(inp$time >= inp$timerange[2])
    # Management
    if(!"ffac" %in% names(inp)) inp$ffac <- 1
    if("ffac" %in% names(inp)){
        if(inp$ffac < 0){
            cat('Warning: ffac < 0, which is not allowed, setting ffac = 0.')
            inp$ffac <- 0
        }
    }
    if(!"fcon" %in% names(inp)) inp$fcon <- 0
    if("fcon" %in% names(inp)){
        if(inp$fcon < 0){
            cat('Warning: fcon < 0, which is not allowed, setting fcon = 0.')
            inp$fcon <- 0
        }
    }
    if(!"ffacvec" %in% names(inp)){
        inp$ffacvec <- numeric(inp$ns) + 1
        # -1 in indpred because 1 is for plotting
        inp$ffacvec[inp$indpred[-1]] <- inp$ffac + 1e-8 # Add small to avoid taking log of 0
    }
    if(!"fconvec" %in% names(inp)){
        inp$fconvec <- numeric(inp$ns)
        # -1 in indpred because 1 is for plotting
        inp$fconvec[inp$indpred[-1]] <- inp$fcon + 1e-8 # Add small to avoid taking log of 0
    }
    inp$ffaceuler <- inp$ffac^inp$dteuler
    # Seasons
    if(!"nseasons" %in% names(inp)){
        expnseasons <- 1/min(inp$dtc)
        if(expnseasons >= 4){
            inp$nseasons <- 4
        } else {
            inp$nseasons <- 1
        }
        #inp$nseasons <- length(unique(timeobs %% 1))
        #if(inp$seasons>4) inp$seasons <- 4
    }
    if("nseasons" %in% names(inp)){
       if(!inp$nseasons %in% c(1, 2, 4)) stop('inp$nseasons (=', inp$nseasons, ') must be either 1, 2 or 4.')
    }
    if(inp$nseasons == 1) inp$seasontype <- 0 # seasontype = 0 means seasons are disabled.
    # Calculate seasonal spline
    if("splineorder" %in% names(inp)){
        if(inp$nseasons<4 & inp$splineorder>2) inp$splineorder <- 2
    } else {
        inp$splineorder <- ifelse(inp$nseasons<4, 2, 3)
    }
    inp$splinemat <- make.splinemat(inp$nseasons, inp$splineorder, dtfine=inp$dteuler)
    inp$splinematfine <- make.splinemat(inp$nseasons, inp$splineorder, dtfine=1/100)
    inp$seasonindex <- 1/inp$dteuler*(inp$time %% 1)
    inp$seasons <- rep(0, inp$ns)
    for(i in 1:inp$nseasons){
        frac <- 1/inp$nseasons
        modtime <- inp$time %% 1
        inds <- which(modtime>=((i-1)*frac) & modtime<(i*frac))
        inp$seasons[inds] <- i
    }
    # ic is the indices of inp$time to which catch observations correspond
    if(length(inp$dtc)>0){
        dtcpred <- min(inp$dtc)
    } else {
        dtcpred <- 1
    }
    inp$timeCpred <- unique(c(inp$timeC, (seq(tail(inp$timeC,1), inp$timepredc, by=dtcpred))))
    #inp$timeCp <- tail(inp$timeCpred, 1)
    inp$nobsCp <- length(inp$timeCpred)
    inp$dtcp <- c(inp$dtc, rep(dtcpred, inp$nobsCp-inp$nobsC))
    #inp$ic <- match(inp$timeCpred, inp$time)
    inp$ic <- cut(inp$timeCpred, inp$time, right=FALSE, labels=FALSE)
    # nc is number of states to integrate a catch observation over
    inp$nc <- rep(0, inp$nobsCp)
    for(i in 1:inp$nobsCp) inp$nc[i] <- sum(inp$time >= inp$timeCpred[i] & inp$time < (inp$timeCpred[i]+inp$dtcp[i]))
    if(any(inp$nc == 0)) stop('Current inp$dteuler is too large to accommodate some catch intervals. Make inp$dteuler smaller!')
    # ii is the indices of inp$time to which index observations correspond
    inp$ii <- list()
    #for(i in 1:inp$nindex) inp$ii[[i]] <- match(inp$timeI[[i]], inp$time)
    for(i in 1:inp$nindex) inp$ii[[i]] <- cut(inp$timeI[[i]], inp$time, right=FALSE, labels=FALSE)
    # Translate index observations from a list to a vector
    inp$obsIin <- unlist(inp$obsI)
    inp$stdevfacIin <- unlist(inp$stdevfacI)
    inp$iiin <- unlist(inp$ii)
    #inp$iqin <- rep(1:inp$nindex, times=inp$nobsI)
    inp$iqin <- rep(inp$mapq, times=inp$nobsI)
    inp$isdiin <- rep(inp$mapsdi, times=inp$nobsI)
    # Add helper variable such that predicted catch can be calculated using small euler steps
    # Need to include timerange[2] and exclude timerange[2]+dtpred because the catch at t is acummulated over t to t+dtc.
    #inp$dtpredcinds <- which(inp$time >= inp$timerange[2] & inp$time < (inp$timerange[2]+inp$dtpredc))
    inp$dtpredcinds <- which(inp$time >= inp$timepredc & inp$time < (inp$timepredc+inp$dtpredc))
    inp$dtpredcnsteps <- length(inp$dtpredcinds)
    #inp$dtprediind <- which(inp$time == inp$timepredi)
    inp$dtprediind <- cut(inp$timepredi, inp$time, right=FALSE, labels=FALSE)
    

    # - Sort observations in time and store in one vector -
    timeobsseen <- c(inp$timeC+inp$dtc-1e-4, unlist(inp$timeI)) # Add dtc to timeC because the total catch is first "seen" at the end of the given catch interval (typically year or quarter)
    srt <- sort(timeobsseen, index=TRUE)
    timeobs <- c(inp$timeC, unlist(inp$timeI))
    timeobssrt <- timeobs[srt$ix]
    obs <- log(c(inp$obsC, unlist(inp$obsI)))
    obsid <- c(inp$obsidC, unlist(inp$obsidI))
    inp$obssrt <- obs[srt$ix]
    inp$timeobssrt <- timeobs[srt$ix]
    inp$obsidsrt <- obsid[srt$ix]
    inp$isc <- match(1:inp$nobsC, srt$ix)
    inp$isi <- match((inp$nobsC+1):(inp$nobsC+sum(inp$nobsI)), srt$ix)
    if(sum(inp$nobsI) != length(inp$isi)){
        cat('Warning: Mismatch between length(inp$isi)', length(inp$isi), 'and sum(inp$nobsI)', sum(inp$nobsI), '\n')
    }
    inp$osar.conditional <- which(inp$timeobssrt < inp$time[1]+1) # Condition on the first year of data.
    inp$osar.subset <- setdiff(1:length(inp$obssrt), inp$osar.conditional)

    # -- PRIORS --
    # Priors are assumed Gaussian and specified in a vector of length 3: c(log(mean), stdev in log, useflag).
    # log(mean): log of the mean of the prior distribution.
    # stdev in log: standard deviation of the prior distribution in log domain.
    # useflag: if 1 then the prior is used, if 0 it is not used. Default is 0.
    check.prior <- function(priors, priorname){
        priorvec <- priors[[priorname]]
        if(priorname %in% repriors){ # RE priors
            if(length(priorvec) < 3){
                priorvec <- rep(0, 5)
                warning('Invalid prior length specified for', priorname, ', must be 3 (without useflag or 4 (with useflag). Not using this prior.')
            }
            if(length(priorvec) == 3){
                #warning('Length of ', priorname, ' is 3. Proceeding assuming useflag has not been specified.')
                priorvec <- c(priorvec[1:2], 1, priorvec[3])
            }
            if(length(priorvec) == 4){
                ib <- match(priorvec[4], inp$time)
                if(is.na(ib)){
                    ib <- 0
                    priorvec[3] <- 0
                    warning('Year for prior on ', priorname, ' (', priorvec[3], ') did not match times where this RE is estimated. Not using this prior. To fix this use a year where an observation is available.')
                }
                priorvec <- c(priorvec, ib) # Add index in time vec to which this year corresponds
            }
        } else { # FE priors
            if(!length(priorvec) %in% 2:3){
                priorvec <- rep(0, 3)
                warning('Invalid prior length specified for', priorname, ', must be 2 (without useflag or 3 (with useflag). Not using this prior.')
            }
            if(length(priorvec) == 2){
                #warning('Length of ', priorname, ' is 2. Proceeding assuming useflag has not been specified.')
                priorvec <- c(priorvec, 1)
            }
        }
        if(priorvec[3] == 1){
            # Check st dev
            if(priorvec[2] <= 0){
                warning('Invalid standard deviation specified in prior for', priorname, '(must be > 0). Not using this prior.')
                priorvec[3] <- 0
            }
        }
        return(priorvec)
    }
    possiblepriors <- c('logn', 'logalpha', 'logbeta', 'logr', 'logK', 'logm', 'logq', 'logbkfrac', 'logB', 'logF', 'logBBmsy', 'logFFmsy', 'logsdb', 'logsdf', 'logsdi', 'logsdc')
    repriors <- c('logB', 'logF', 'logBBmsy', 'logFFmsy')
    npossiblepriors <- length(possiblepriors)
    if(!"priors" %in% names(inp)){
        inp$priors <- list()
    }
    # Default priors
    lognflag <- TRUE
    logalphaflag <- TRUE
    logbetaflag <- TRUE
    logn <- log(2)
    logalpha <- log(1)
    logbeta <- log(1)
    small <- 1e-3
    wide <- 2 # Value suggested by Casper 
    lognsd <- wide
    logalphasd <- wide
    logbetasd <- wide
    # Default phases are not set at this point, so if phases exist the user has chosen them manually
    # If a phase is set > -1 then the parameter is estimated without prior.
    # If a phase is set <= -1 then a narrow priors is used for this parameter fixing it.
    if('phases' %in% names(inp)){ # Only use default priors if not assigned a phase
        if('logn' %in% names(inp$phases)){
            if(inp$phases$logn > -1){
                lognflag <- FALSE # Estimate logn
            } else {
                lognsd <- small # Fix n
            }
        }
        if('logalpha' %in% names(inp$phases)){
            if(inp$phases$logalpha > -1){
                logalphaflag <- FALSE # Estimate logalpha i.e. sdb and sdi untied
            } else {
                logalphasd <- small # Fix alpha
            }
        }
        if('logbeta' %in% names(inp$phases)){  
            if(inp$phases$logbeta > -1){
                logbetaflag <- FALSE # Estimate logalpha i.e. sdf and sdc untied
            } else {
                logbetasd <- small # Fix beta
            }
        }
    }
    # Default ini values have not been set yet at this point
    # If an ini value is available the user has set it manually
    # If an ini value is available it is used as mean of the prior for that parameter.
    # The sd of the prior depends on whether a phase has been specified.
    if('ini' %in% names(inp)){ 
        # Log n
        if('logn' %in% names(inp$ini)){
            logn <- inp$ini$logn
        }
        # Log alpha
        if(!'logalpha' %in% names(inp$ini) & 'logsdb' %in% names(inp$ini) & 'logsdi' %in% names(inp$ini)){
            inp$ini$logalpha <- inp$ini$logsdi - inp$ini$logsdb
        }
        if('logalpha' %in% names(inp$ini)){
            logalpha <- inp$ini$logalpha
        }
        # Log beta
        if(!'logbeta' %in% names(inp$ini) & 'logsdf' %in% names(inp$ini) & 'logsdc' %in% names(inp$ini)){
            inp$ini$logbeta <- inp$ini$logsdc - inp$ini$logsdf
        }
        if('logbeta' %in% names(inp$ini)){
            logbeta <- inp$ini$logbeta
        }
    }
    if(!'logn' %in% names(inp$priors) & lognflag) inp$priors$logn <- c(logn, lognsd)
    if(!'logalpha' %in% names(inp$priors) & logalphaflag) inp$priors$logalpha <- c(logalpha, logalphasd) 
    if(!'logbeta' %in% names(inp$priors) & logbetaflag) inp$priors$logbeta <- c(logbeta, logbetasd)
    if("priors" %in% names(inp)){
        # Remove wrong priors names
        nms <- names(inp$priors)
        inds <- which(is.na(match(nms, possiblepriors)))
        if(length(inds)>0){
            warning('Wrong prior names specified: ', nms[inds])
            inp$priors[inds] <- NULL
        }
        # Check priors
        for(nm in possiblepriors){
            if(!nm %in% names(inp$priors)){
                # Set default prior values. These will not be used
                if(nm %in% repriors){
                    inp$priors[[nm]] <- c(log(1e-4), 0.2, 0, 2000, 0) # RE prior
                } else {
                    inp$priors[[nm]] <- c(log(1e-4), 0.2, 0) # FE priors
                }
            } else {
                inp$priors[[nm]] <- check.prior(inp$priors, nm)
            }
        }
    }
    npriors <- length(inp$priors)
    inp$priorsuseflags <- numeric(npriors)
    for(i in 1:npriors) inp$priorsuseflags[i] <- inp$priors[[i]][3]
           
    # -- MODEL PARAMETERS --
    # logn
    if(!"logn" %in% names(inp$ini)){
        inp$ini$logn <- logn
    } else {
        if(inp$ini$logn==0) stop('Initial value for logn == 0, that is not valid!')
    }
    # Calculate gamma from n
    n <- exp(inp$ini$logn)
    inp$ini$gamma <- calc.gamma(n)
    # logK
    if(!'logK' %in% names(inp$ini)) inp$ini$logK <- log(4*max(inp$obsC))
    # logr
    if('logr' %in% names(inp$ini) & 'logm' %in% names(inp$ini)) inp$ini$logr <- NULL # If both r and m are specified use m and discard r
    if(!'logr' %in% names(inp$ini)){
        if(!'logm' %in% names(inp$ini)){
            inp$ini$logm <- log(guess.m(inp))
        }
        r <- inp$ini$gamma * exp(inp$ini$logm) / exp(inp$ini$logK) # n > 1
        if(n>1){
            inp$ini$logr <- log(r)
        } else {
            inp$ini$logr <- log(-r)
        }
    }
    if('logr' %in% names(inp$ini)){
        nr <- length(inp$ini$logr)
        if(!'ir' %in% names(inp) | nr==1){
            inp$ir <- rep(0, inp$ns)
            for(i in 1:nr){
                frac <- 1/nr
                modtime <- inp$time %% 1
                inds <- which(modtime>=((i-1)*frac) & modtime<(i*frac))
                inp$ir[inds] <- i
            }
        } else {
            if(length(unique(inp$ir)) != nr) stop('Mismatch between specified inp$ir and inp$ini$logr!')
            nir <- length(inp$ir)
            if(nir != inp$ns){
                if(nir == inp$nobsC){ # Assume that inp$ir fits with inp$timeC
                    ir <- rep(0, inp$ns)
                    for(i in 1:nir){
                        inds <- which(inp$time >= inp$timeC[i] & inp$time < (inp$timeC[i]+inp$dtc[i]))
                        ir[inds] <- inp$ir[i]
                    }
                    inds <- which(inp$time >= inp$timeC[nir])
                    ir[inds] <- inp$ir[nir]
                    inp$ir <- ir
                } else {
                    if(nir == inp$nobsI[[1]]){ # Assume that inp$ir fits with inp$timeI[[1]]
                        ir <- rep(0, inp$ns)
                        for(i in 2:nir){
                            inds <- which(inp$time >= inp$timeI[[1]][i-1] & inp$time < inp$timeI[[1]][i])
                            ir[inds] <- inp$ir[i]
                        }
                        inds <- which(inp$time >= inp$timeI[[1]][nir])
                        ir[inds] <- inp$ir[nir]
                        inp$ir <- ir
                    } else {
                        stop('inp$ir is misspecified, only advanced users should specify this manually!')
                    }
                }
            }
        }
    }
    if(!'logq' %in% names(inp$ini)) inp$ini$logq <- log(max(inp$obsI[[1]])) - inp$ini$logK
    if('logq' %in% names(inp$ini)){
        if(length(inp$ini$logq) != inp$nq){ # nq is given by mapq
            if(length(inp$ini$logq) == 1){
                inp$ini$logq <- rep(inp$ini$logq, inp$nq)
            } else {
                stop('The length of inp$ini$logq (', length(inp$ini$logq), ') does not fit with the number of qs to be estimated (', inp$nq, ')')
            }
        }
    }
    if(!'logsdf' %in% names(inp$ini)) inp$ini$logsdf <- log(0.2)
    if(!'logsdu' %in% names(inp$ini)) inp$ini$logsdu <- log(0.1)
    if(!'logsdb' %in% names(inp$ini)) inp$ini$logsdb <- log(0.2)
    if(!'logsdc' %in% names(inp$ini)) inp$ini$logsdc <- log(0.2)
    if(!'logsdi' %in% names(inp$ini)) inp$ini$logsdi <- log(0.2)
    if('logsdi' %in% names(inp$ini)){
        if(length(inp$ini$logsdi) != inp$nsdi){ # nsdi is given by mapsdi
            if(length(inp$ini$logsdi) == 1){
                inp$ini$logsdi <- rep(inp$ini$logsdi, inp$nsdi)
            } else {
                stop('The length of inp$ini$logsdi (', length(inp$ini$logsdi), ') does not fit with the number of sdis to be estimated (', inp$nsdi, ')')
            }
        }
    }
    if(!"logalpha" %in% names(inp$ini)) inp$ini$logalpha <- logalpha
    if('logalpha' %in% names(inp$ini)){
        if(length(inp$ini$logalpha) != inp$nsdi){ # nsdi is given by mapsdi
            if(length(inp$ini$logalpha) == 1){
                inp$ini$logalpha <- rep(inp$ini$logalpha, inp$nsdi)
            } else {
                stop('The length of inp$ini$logalpha (', length(inp$ini$logalpha), ') does not fit with the number of sdis to be estimated (', inp$nsdi, ')')
            }
        }
    }
    if(!"logbeta" %in% names(inp$ini))  inp$ini$logbeta <- logbeta

    #if('logsdi' %in% names(inp$ini) & 'logalpha' %in% names(inp$ini)){
    #    if(inp$onealpha | inp$onesdi){
    #           if(length(inp$ini$logsdi) != 1) inp$ini$logsdi <- log(0.2)
    #           if(length(inp$ini$logalpha) != 1) inp$ini$logalpha <- log(1)
    #    } else {
    #        if(length(inp$ini$logsdi) != inp$nindex){
    #            if(length(inp$ini$logsdi) == 1){
    #                inp$ini$logsdi <- rep(inp$ini$logsdi, inp$nindex)
    #            } else {
    #                stop('The length of inp$ini$logsdi (', length(inp$ini$logsdi), ') does not fit with the number of index series (', inp$nindex, ')')
    #            }
    #        }
    #        if(length(inp$ini$logalpha) != inp$nindex){
    #            if(length(inp$ini$logalpha) == 1){
    #                inp$ini$logalpha <- rep(inp$ini$logalpha, inp$nindex)
    #            } else {
    #                stop('The length of inp$ini$logalpha (', length(inp$ini$logalpha), ') does not fit with the number of index series (', inp$nindex, ')')
    #            }
    #        }
    #    }
    #}

    if(!"logm" %in% names(inp$ini)){
        gamma <- inp$ini$gamma
        r <- exp(inp$ini$logr)
        K <- exp(inp$ini$logK)
        m <- r * K / gamma # n > 1
        if(n>1){
            inp$ini$logm <- log(m)
        } else {
            inp$ini$logm <- log(-m)
        }
    }
    # Fill in unspecified (more rarely user defined) model parameter values
    if(!"loglambda" %in% names(inp$ini)) inp$ini$loglambda <- log(0.1)
    #if("lambda" %in% names(inp)) if(inp$lambda <= 0) cat('Error: lambda must be positive!')
    if("logphi" %in% names(inp$ini)){
        if(length(inp$ini$logphi)+1 != dim(inp$splinemat)[2]){
            cat('Mismatch between length of ini$logphi and number of columns of splinemat! removing prespecified ini$logphi and setting default.\n')
            inp$ini$logphi <- NULL
        }
    }
    if(!"logphi" %in% names(inp$ini)) inp$ini$logphi <- rep(0, inp$nseasons-1)
    if(!"logitpp" %in% names(inp$ini)) inp$ini$logitpp <- log(0.95/(1-0.95))
    if(!"logp1robfac" %in% names(inp$ini)) inp$ini$logp1robfac <- log(15-1)
    if(!"logbkfrac" %in% names(inp$ini)) inp$ini$logbkfrac <- log(0.8)
    #if(!"logp" %in% names(inp$ini)) inp$ini$logp <- log(1.0)
    #if(!"logF0" %in% names(inp$ini)) inp$ini$logF0 <- log(0.2*exp(inp$ini$logr[1]))
    if(!"logF" %in% names(inp$ini)){
        #inp$ini$logF <- rep(inp$ini$logF0, inp$ns)
        inp$ini$logF <- rep(log(0.2) + inp$ini$logr[1], inp$ns)
    } else {
        if(length(inp$ini$logF) != inp$ns){
            cat('Wrong length of inp$ini$logF:', length(inp$ini$logF), ' Should be equal to inp$ns:', inp$ns, ' Setting length of logF equal to inp$ns (removing beyond inp$ns).\n')
            inp$ini$logF <- inp$ini$logF[1:inp$ns]
        }
    }
    if("logu" %in% names(inp$ini)){
        if(dim(inp$ini$logu)[1] != 2*length(inp$ini$logsdu) & dim(inp$ini$logu)[2] != inp$ns){
            cat('Wrong dimension of inp$ini$logu:', dim(inp$ini$logu)[1], 'x', dim(inp$ini$logu)[2], ' should be equal to 2*length(inp$ini$logsdu) x inp$ns: ', 2*length(inp$ini$logsdu), 'x', inp$ns,', Filling with log(1).\n')
            inp$ini$logu <- NULL
        }
    }
    if(!"logu" %in% names(inp$ini)){
        inp$ini$logu <- matrix(log(1)+1e-3, 2*length(inp$ini$logsdu), inp$ns)
    }
    if(!"logB" %in% names(inp$ini)){
        inp$ini$logB <- rep(inp$ini$logK + log(0.5), inp$ns)
    } else {
        if(length(inp$ini$logB) != inp$ns){
            cat('Wrong length of inp$ini$logB:', length(inp$ini$logB), ' Should be equal to inp$ns:', inp$ns, ' Setting length of logF equal to inp$ns (removing beyond inp$ns).\n')
            inp$ini$logB <- inp$ini$logB[1:inp$ns]
        }
    }

    # Reorder parameter list
    inp$parlist <- list(logm=inp$ini$logm,
                        logK=inp$ini$logK,
                        logq=inp$ini$logq,
                        logn=inp$ini$logn,
                        logsdb=inp$ini$logsdb,
                        logsdu=inp$ini$logsdu,
                        logsdf=inp$ini$logsdf,
                        logsdi=inp$ini$logsdi,
                        logsdc=inp$ini$logsdc,
                        #logalpha=inp$ini$logalpha,
                        #logbeta=inp$ini$logbeta,
                        logphi=inp$ini$logphi,
                        loglambda=inp$ini$loglambda,
                        logitpp=inp$ini$logitpp,
                        logp1robfac=inp$ini$logp1robfac,
                        logF=inp$ini$logF,
                        logu=inp$ini$logu,
                        logB=inp$ini$logB)
    
    # Determine phases and fixed parameters
    fixpars <- c('logalpha', 'logbeta', 'logn', 'logitpp', 'logp1robfac') # These are fixed unless otherwise specified
    forcefixpars <- c() # Parameters that are forced to be fixed.
    if(inp$nseasons==1){
        forcefixpars <- c('logphi', 'logu', 'logsdu', 'loglambda', forcefixpars)
    } else {
        if(inp$seasontype==1){ # Use spline
            forcefixpars <- c('logu', 'logsdu', 'loglambda', forcefixpars)
        }
        if(inp$seasontype==2){ # Use coupled SDEs
            forcefixpars <- c('logphi', forcefixpars)
        }
    }
    if(!"phases" %in% names(inp)){
        inp$phases <- list()
    } else {
        if("logr" %in% names(inp$phases)){
            inp$phases$logm <- inp$phases$logr
            inp$phases$logr <- NULL
        }
        nms <- names(inp$phases)
        inds <- which(is.na(match(nms, c(names(inp$ini), 'logalpha', 'logbeta'))))
        if(length(inds)>0){
            stop('phase specified for invalid parameter(s): ', paste(nms[inds], collapse=', '))
        }
    }
    # If robust flags are set to 1 then set phases for robust parameters to 1
    if(inp$robflagc==1 | inp$robflagi==1){
        if(!"logitpp" %in% names(inp$phases)) inp$phases$logitpp <- 1
        if(!"logp1robfac" %in% names(inp$phases)) inp$phases$logp1robfac <- 1
    }
    if("phases" %in% names(inp)){
        for(i in 1:length(forcefixpars)) inp$phases[[forcefixpars[i]]] <- -1
        for(i in 1:length(fixpars)) if(!fixpars[i] %in% names(inp$phases)) inp$phases[[fixpars[i]]] <- -1
    }
    # Assign phase 1 to parameters without a phase
    nms <- names(inp$parlist)
    nnms <- length(nms)
    for(nm in nms){
        if(!nm %in% names(inp$phases)){
            inp$phases[[nm]] <- 1
        }
    }
    # If a parameter is assigned a prior then estimate this parameter
    inds <- which(inp$priorsuseflags == 1)
    priornmsuse <- names(inp$priors)[inds]
    for(nm in priornmsuse) if(nm %in% names(inp$phases)) inp$phases[[nm]] <- 1
    
    nphasepars <- length(inp$phases)
    phasevec <- unlist(inp$phases)
    phases <- unique(unname(phasevec))
    phases <- phases[phases>0] # Don't include phase -1 (fixed value)
    inp$nphases <- length(phases)
    inp$map <- list()
    for(j in 1:inp$nphases){
        inp$map[[j]] <- list() # Map for phase j
        inds <- which(phasevec > j | phasevec == -1)
        for(i in inds){
            parnam <- names(phasevec)[i]
            if(parnam %in% names(inp$parlist)){
                phase <- inp$phases[[parnam]]
                inp$map[[j]][[parnam]] <- factor(rep(NA, length(inp$parlist[[parnam]])))
            } else {
                warning('Phase specified for an invalid parameter:', parnam)
            }
        }
    }
    if(!is.null(inp)) class(inp) <- "spictcls"
    return(inp)
}


#' @name list.possible.priors
#' @title List parameters to which priors can be added
#' @return Prints parameters to which priors can be added.
#' @export
list.possible.priors <- function(){
    data(pol)
    inp <- check.inp(pol$albacore)
    print(names(inp$priors))
}
