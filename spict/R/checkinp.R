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


#' @name check.inp
#' @title Check list of input variables
#'
#' @param inp List of input variables, see details for required variables.
#' @param verbose Should detailed outputs be provided (default: TRUE).
#' @param mancheck Should the time-dependent objects in \code{inp} be
#'     checked against the management time and corrected if necessary? (Default: TRUE)
#'
#'
#' @details Fills in default values if missing.
#'
#' Required inputs:
#'
#' \itemize{
#'  \item{"inp$obsC"}{ Vector of catch observations.}
#'  \item{"inp$obsI and/or inp$obsE"}{ List containing vectors of index observations and/or a vector of effort information.}
#' }
#'
#' Optional inputs:
#'
#' - Data
#' \itemize{
#'  \item{"inp$timeC"}{ Vector of catch times. Default: even time steps starting at 1.}
#'  \item{"inp$timeI"}{ List containing vectors of index times. Default: even time steps starting at 1.}
#'  \item{"inp$timeE"}{ Vector of effort times. Default: even time steps starting at 1.}
#'  \item{"inp$dtc"}{ Time interval for catches, e.g. for annual catches inp$dtc=1, for quarterly catches inp$dtc=0.25. Can be given as a scalar, which is then used for all catch observations. Can also be given as a vector specifying the catch interval of each catch observation. Default: min(diff(inp$timeC)). }
#'  \item{"inp$dte"}{ Time interval for effort observations. For annual effort inp$dte=1, for quarterly effort inp$dte=0.25. Default: min(diff(inp$timeE)). }
#'  \item{"inp$nseasons"}{ Number of within-year seasons in data. If inp$nseasons > 1 then a seasonal pattern is used in F. Valid values of inp$nseasons are 1, 2 or 4. Default: number of unique within-year time points present in data.}
#' }
#'
#' - Initial parameter values
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
#' - Initial values for unobserved states estimated as random effects
#'
#' \itemize{
#'   \item{"inp$ini$logF"}{ Log fishing mortality. Default: log(0.2*r), with r derived from m and K.}
#'   \item{"inp$ini$logB"}{ Log biomass. Default: log(0.5*K).}
#'   \item{"inp$ini$logU"}{ Log U, the state of the coupled SDE representation of seasonality. Default: log(1).}
#' }
#'
#' - Priors
#'
#' Priors on model parameters are assumed generally assumed Gaussian and specified in a vector of length 2: c(log(mean), stdev in log domain, useflag [optional]). NOTE: if specifying a prior for a value in a temporal vector e.g. logB, then a fourth element is required specifying the year the prior should be applied.
#' log(mean): log of the mean of the prior distribution.
#' stdev in log: standard deviation of the prior distribution in log domain.
#' useflag: if 1 then the prior is used, if 0 it is not used. Default is 1.
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
#' \item{"inp$maninterval"}{ Start and end time of management period. Default: One year interval starting at the beginning of the new year after the last observation. Example: inp$maninterval <- c(2020.25,2021.25)}
#' \item{"inp$maneval"}{ Time for the estimation of predicted model states (biomass and fishing mortality), which can be used to evaluate the implications of management scenarios. Default: At the end of the management interval \code{inp$maninterval[2]}. Example: inp$maneval <- 2021.25}
#'  \item{"inp$timepredc"}{ Deprecated: Predict accummulated catch in the interval starting at $timepredc and $dtpredc into the future. Default depends on \code{inp$maninterval}.}
#'  \item{"inp$dtpredc"}{ Deprecated: Length of catch prediction interval in years. Default depends on \code{inp$maninterval}.}
#'  \item{"inp$timepredi"}{ Deprecated: Predict index until this time. Default depends on \code{inp$maneval}.}
#' \item{"inp$manstart"}{ Deprecated: Start of the management period. Updated argument \code{inp$maninterval}. Default depends on \code{inp$maninterval}.}
#'  \item{"inp$do.sd.report"}{ Flag indicating whether SD report (uncertainty of derived quantities) should be calculated. For small values of inp$dteuler this may require a lot of memory. Default: TRUE.}
#'  \item{"inp$reportall"}{ Flag indicating whether quantities derived from state vectors (e.g. B/Bmsy, F/Fmsy etc.) should be calculated by SD report. For small values of inp$dteuler (< 1/32) reporting all may have to be set to FALSE for sdreport to run. Additionally, if only reference points of parameter estimates are of interest one can set to FALSE to gain a speed-up. Default: TRUE.}
#' \item{"inp$reportmode"}{ Integer between 0 and 2 determining which objects will be adreported. Default: 0 = all quantities are adreported. Example: inp$reportmode <- 1}
#' \item{"inp$reportRel"}{ Flag indicating whether mean 1 standardized states (i.e. B/mean(B), F/mean(F) etc.) should be calculated by SD report. Default: FALSE.}
#' \item{"inp$robflagc"}{ Flag indicating whether robust estimation should be used for catches (either 0 or 1). Default: 0.}
#'  \item{"inp$robflagi"}{ Vector of flags indicating whether robust estimation should be used for indices (either 0 or 1). Default: 0.}
#'  \item{"inp$ffac"}{ Management scenario represented by a factor to multiply F with when calculating the F of the next time step. ffac=0.8 means a 20\% reduction in F over the next year. The factor is only used when predicting beyond the data set. Default: 1 (0\% reduction).}
#'  \item{"inp$dteuler"}{ Length of Euler time step in years. Default: 1/16 year.}
#'  \item{"inp$phases"}{ Phases can be used to fix/free parameters and estimate in different stages or phases. To fix e.g. logr at inp$ini$logr set inp$phases$logr <- -1. To free logalpha and estimate in phase 1 set inp$phases$logalpha <- 1.}
#'  \item{"inp$osar.method"}{ Method to use in TMB's oneStepPredict function. Valid methods include: "oneStepGaussianOffMode", "fullGaussian", "oneStepGeneric", "oneStepGaussian", "cdf". See TMB help for more information. Default: "none" (i.e. don't run this).}
#'  \item{"inp$osar.trace"}{ If TRUE print OSAR calculation progress to screen. Default: FALSE.}
#'  \item{"inp$osar.parallel"}{ If TRUE parallelise OSAR calculation for speed-up. Default: FALSE.}
#'  \item{"inp$catchunit"}{ Specify unit of catches to be used in plotting legends. Default: ''.}
#'  \item{"inp$stdevfacC"}{ Factors to multiply the observation error standard deviation of each individual catch observation. Can be used if some observations are more uncertain than others. Must be same length as observation vector. Default: 1.}
#'  \item{"inp$stdevfacI"}{ Factors to multiply the observation error standard deviation of each individual index observation. Can be used if some observations are more uncertain than others. A list with vectors of same length as observation vectors. Default: 1.}
#'  \item{"inp$stdevfacE"}{ Factors to multiply the observation error standard deviation of each individual effort observation. Can be used if some observations are more uncertain than others. A list with vectors of same length as observation vectors. Default: 1.}
#'  \item{"inp$mapsdi"}{ Vector of length equal to the number of index series specifying which indices that should use the same sdi. For example: in case of 3 index series use \code{inp$mapsdi <- c(1, 1, 2)} to have series 1 and 2 share sdi and have a separate sdi for series 3. Default: 1:nindex, where nindex is number of index series.}
#'  \item{"inp$seasontype"}{ If set to 1 use the spline-based representation of seasonality. If set to 2 use the oscillatory SDE system (this is more unstable and difficult to fit, but also more flexible).}
#'  \item{"inp$sim.random.effects"}{Should random effects (logB, logF, etc.) be simulated (default) or the same random effects be used (as specified in \code{inp$ini} or in an fitted spict object)?}
#'  \item{"inp$sim.fit"}{Should the estimated parameters from the last fit of a fitted spict object be used for simulation (\code{env$last.par}, default) or the inital values (specified in \code{inp$ini})?. Note, that this only works if a fitted spict object is provided as input to \code{sim.spict}.}
#' }
#'
#' @return An updated list of input variables checked for consistency and with defaults added.
#'
#' @examples
#' data(pol)
#' (inp <- check.inp(pol$albacore))
#'
#' @export
check.inp <- function(inp, verbose = TRUE, mancheck = TRUE){
    ## # Check management settings if inp is 'checked' list
    ## isChecked <- ifelse(inherits(inp, "spictcls"), 1, 0)
    ## if(isChecked && mancheck){
    ##     inp <- check.man.time(inp, printTimeline = FALSE, verbose = verbose)
    ## }

    set.default <- function(inpin, key, val){
        if (!key %in% names(inpin)){
            inpin[[key]] <- val
        }
        return(inpin)
    }
    check.ini <- function(parname, inp, min=NULL, max=NULL){
        if (!parname %in% names(inp$ini)){
            stop('Please specify an initial value for ', parname, '!')
        }
    }
    rm.neg <- function(tmpin, nam, nms, j=''){
        tmpout <- tmpin
        if (!is.null(tmpin[[nms[1]]])){
            nnms <- length(nms)
            neg <- which(tmpin[[nms[1]]]<=0 | is.na(tmpin[[nms[1]]]))
            if (length(neg)>0){
                for (i in 1:nnms){
                    tmpout[[nms[i]]] <- tmpin[[nms[i]]][-neg]
                }
                if(verbose) cat('Removing zero, negative, and NAs in ', nam, ' series ', j, ' \n')
            }
            return(tmpout)
        }
    }
    remove.neg <- function(inp, nam, extrakeys=NULL){
        nms <- c(paste0(c('obs', 'time', 'stdevfac'), nam), extrakeys)
        nnms <- length(nms)
        flag <- class(inp[[nms[1]]]) == 'list'
        if (flag){
            nna <- length(inp[[nms[1]]])
            if (nna>0){
                for (j in 1:nna){
                    tmp <- list()
                    for (i in 1:nnms){
                        tmp[[nms[i]]] <- inp[[nms[i]]][[j]]
                    }
                    tmpup <- rm.neg(tmp, nam, nms, j)
                    for (i in 1:nnms){
                        inp[[nms[i]]][[j]] <- tmpup[[nms[i]]]
                    }
                }
            }
        } else {
            tmp <- list()
            if (length(inp[[nms[1]]] > 0)){ # nms[1] is the obs name
                for (i in 1:nnms){
                    tmp[[nms[i]]] <- inp[[nms[i]]]
                }
                tmpup <- rm.neg(tmp, nam, nms)
                for (i in 1:nnms){
                    inp[[nms[i]]] <- tmpup[[nms[i]]]
                }
            }
        }
        return(inp)
    }
    base.checks <- function(obs, time, stdevfac, nam){
        if (length(obs) != length(time)){
            stop('Time and observation vector do not match in length for  ', nam, ' series.')
        }
        if (length(obs) != length(stdevfac)){
            stop('stdevfac and observation vector do not match in length for ', nam, ' series.')
        }
        if (sum(stdevfac <= 0) > 0){
            stop('Non-positive values entered in stdevfac for ', nam, 'series.')
        }
    }
    make.list <- function(invar){
        if (class(invar)!='list'){
            outvar <- list()
            outvar[[1]] <- invar
        } else {
            outvar <- invar
        }
        return(outvar)
    }
    make.time <- function(inp, nam){
        time <- paste0('time', nam)
        obs <- paste0('obs', nam)
        nobs <- length(unlist(inp[[obs]]))
        if (!time %in% names(inp) & nobs > 0){
            if (!'nseasons' %in% names(inp)){
                inp[[time]] <- 0:(nobs-1)
            } else {
                to <- (nobs-1)/inp$nseasons
                inp[[time]] <- seq(0, to, length=nobs)
            }
        }
        return(inp)
    }
    match.times <- function(times, eulertimes){
        return(apply(abs(outer(times, eulertimes, FUN='-')), 1, which.min))
    }

    checkandadd <- function(name, default, cls) {
        if (!name %in% names(inp)) {
            out <- default
        } else if (inherits(inp[[name]], cls)){
            out <- inp[[name]]
        } else{
            stop("The element 'inp$", name, "' has to be a ", cls,
                 " including any of the following elements: ",
                 paste(names(default), collapse = ", "),
                 "!", call. = FALSE)
        }
        inp[[name]] <- default[which(!names(default) %in% names(out))]
        inp[[name]] <<- c(inp[[name]], out)
    }

    # -- DATA --
    # Check catch observations
    if ('obsC' %in% names(inp)){
        inp <- make.time(inp, 'C')
        if (any(diff(inp$timeC)<=0)){
            stop('Catch times are not strictly increasing!')
        }
        # Catch intervals (dtc)
        if (!"dtc" %in% names(inp)){
            dtc <- diff(inp$timeC)
            if (length(dtc) > 0){
                inp$dtc <- min(dtc)
            } else {
                inp$dtc <- 1
                if(verbose) cat(paste('Catch interval (dtc) not specified and length of catch time series shorter than 2. Assuming an interval of 1 year.\n'))
            }
        }
        if (length(inp$dtc) == 1){
            inp$dtc <- rep(inp$dtc, length(inp$obsC))
        }
        if (!'stdevfacC' %in% names(inp)){
            inp$stdevfacC <- rep(1, length(inp$obsC))
        }
        base.checks(inp$obsC, inp$timeC, inp$stdevfacC, 'C')
        inp <- remove.neg(inp, 'C', extrakeys='dtc')
        inp$nobsC <- length(inp$obsC)
        inp$obsidC <- 1:inp$nobsC
    } else {
        stop('No catch observations included. Please include them as a vector in inp$obsC.')
    }
    if (length(inp$dtc) != inp$nobsC){
        stop('Catch interval vector (inp$dtc) does not match catch observation vector (inp$obsC) in length')
    }

    if (!any(c('obsI', 'obsE') %in% names(inp))){
        stop('No effort or index observations. Please include index observations as a vector in inp$obsI and effort observations in inp$obsE.')
    }

    # Check index observations
    inp$obsI <- make.list(inp$obsI)
    inp$nindex <- length(inp$obsI)
    if (inp$nindex > 0){
        inp$nindexseq <- 1:inp$nindex
    }
    # Time vector
    inp <- make.time(inp, 'I')
    inp$timeI <- make.list(inp$timeI)
    # Standard deviation factor
    if (!'stdevfacI' %in% names(inp)){
        inp$stdevfacI <- list()
        for (i in inp$nindexseq){
            inp$stdevfacI[[i]] <- rep(1, length(inp$obsI[[i]]))
        }
    }
    inp$stdevfacI <- make.list(inp$stdevfacI)
    if (inp$nindex != length(inp$timeI)){
        stop('length(inp$timeI) is not equal to length(inp$obsI)!')
    }
    if (inp$nindex != length(inp$stdevfacI)){
        stop('length(inp$stdevfacI) is not equal to length(inp$obsI)!')
    }
    for (i in inp$nindexseq){
        base.checks(inp$obsI[[i]], inp$timeI[[i]], inp$stdevfacI[[i]], paste0('I', i))
    }
    inp <- remove.neg(inp, 'I')
    inp$nobsI <- rep(0, inp$nindex) # Need to be after negative have been removed
    for (i in inp$nindexseq){
        inp$nobsI[i] <- length(inp$obsI[[i]])
    }
    inp$obsidI <- list()
    for (i in inp$nindexseq){
        if (i == 1){
            inp$obsidI[[i]] <- (1:inp$nobsI[i]) + inp$nobsC
        } else {
            inp$obsidI[[i]] <- (1:inp$nobsI[i]) + tail(inp$obsidI[[i-1]], 1)
        }
    }

    # Check effort observations
    inp <- make.time(inp, 'E')
    if (!is.null(inp$timeE)){
        if (class(inp$timeE) != 'numeric' & class(inp$timeE) != 'integer'){
            stop('class(inp$timeE) is not numeric!')
        }
    }
    if (!is.null(inp$obsE)){
        if (class(inp$obsE) != 'numeric' & class(inp$obsE) != 'integer'){
            stop('class(inp$obsE) is not numeric!')
        }
    }
    if (any(diff(inp$timeE) <= 0)){
        stop('Effort times are not strictly increasing!')
    }
    # Effort intervals (dte)
    if (!"dte" %in% names(inp) & length(inp$obsE) > 0){
        dte <- diff(inp$timeE)
        if (length(dte) > 0){
            inp$dte <- min(dte)
        } else {
            inp$dte <- 1
            if(verbose) cat(paste('Effort interval (dte) not specified and length of effort time series shorter than 2. Assuming an interval of 1 year.\n'))
        }
    }
    if (length(inp$dte) == 1){
        inp$dte <- rep(inp$dte, length(inp$obsE))
    }
    if (!'stdevfacE' %in% names(inp)){
        inp$stdevfacE <- rep(1, length(inp$obsE))
    }
    base.checks(inp$obsE, inp$timeE, inp$stdevfacE, 'E')
    inp <- remove.neg(inp, 'E', extrakeys='dte')
    inp$nobsE <- length(inp$obsE)
    if (inp$nobsE > 0){
        inp$obsidE <- (1:inp$nobsE) + inp$nobsC + sum(inp$nobsI)
    } else {
        inp$obsidE <- numeric()
    }
    if (length(inp$dte) != inp$nobsE){
        stop('Effort interval vector (inp$dte, ', length(inp$dte),
             ') does not match effort observation vector (inp$obsE, ',
             inp$nobsE, ') in length')
    }
    inp$nseries <- 1 + inp$nindex + as.numeric(inp$nobsE > 0)

    # -- MODEL OPTIONS --
    if (!"RE" %in% names(inp)) inp$RE <- c('logF', 'logu', 'logB', 'logmre','SARvec')
    if (!"scriptname" %in% names(inp)) inp$scriptname <- 'spict'
    # Index related
    if (!"onealpha" %in% names(inp)){
        if (!"onesdi" %in% names(inp)){
            inp$onealpha <- FALSE
        } else {
            inp$onealpha <- inp$onesdi
        }
    }
    if (!"onesdi" %in% names(inp)) inp$onesdi <- inp$onealpha
    if (!"mapsdi" %in% names(inp)){
        inp$mapsdi <- inp$nindexseq
        if ("onealpha" %in% names(inp)) if (inp$onealpha) inp$mapsdi <- rep(1, inp$nindex)
    }
    if ("mapsdi" %in% names(inp)) inp$nsdi <- length(unique(inp$mapsdi))
    if (!"mapq" %in% names(inp)) inp$mapq <- inp$nindexseq
    if ("mapq" %in% names(inp)) inp$nq <- length(unique(inp$mapq))
    # Effort related
    if (!"efforttype" %in% names(inp)) inp$efforttype <- 1
    if (!"effortmodel" %in% names(inp)) inp$effortmodel <- 'RW'
    if (!"mapsde" %in% names(inp)) inp$mapsde <- inp$neffortseq
    if ("mapsde" %in% names(inp)) inp$nsde <- length(unique(inp$mapsde))
    if (!"mapqf" %in% names(inp)) inp$mapqf <- inp$neffortseq
    if ("mapqf" %in% names(inp)) inp$nqf <- length(unique(inp$mapqf))
    # Catch related
    if (!"catchunit" %in% names(inp)) inp$catchunit <- ''
    # Reporting
    if (!"reportall" %in% names(inp)) inp$reportall <- TRUE
    if (!"reportRel" %in% names(inp)) inp$reportRel <- FALSE
    if (!"do.sd.report" %in% names(inp)) inp$do.sd.report <- TRUE
    if (!"bias.correct" %in% names(inp)) inp$bias.correct <- FALSE # This is time consuming
    if (!"bias.correct.control" %in% names(inp)) inp$bias.correct.control <- list(sd=FALSE) # This is time consuming
    if (!"getReportCovariance" %in% names(inp)) inp$getReportCovariance <- TRUE # Set to FALSE to save memory
    if (!"getJointPrecision" %in% names(inp)){
        inp$getJointPrecision <- FALSE
    }
    # Simulation options
    #if (!"armalistF" %in% names(inp)) inp$armalistF <- list() # Used for simulating arma noise for F instead of white noise.
    # When simulating using sim.comm.cpue == TRUE, the simulation calculates
    # commercial CPUE by dividing catch and effort.
    # This is output as an biomass index.
    if (!"sim.comm.cpue" %in% names(inp)) inp$sim.comm.cpue <- FALSE
    # Optimiser options
    if (!"optimiser" %in% names(inp)) inp$optimiser <- 'nlminb'
    if (!"optimiser.control" %in% names(inp)) inp$optimiser.control <- list(iter.max = 1e4, eval.max = 1e4)
    if (!"optim.method" %in% names(inp)) inp$optim.method <- 'BFGS'
    if (!"stabilise" %in% names(inp)) inp$stabilise <- 1 # If 1 wide uninformative priors are imposed on some parameters to stabilise optimisation (this happens inside the cpp file)
    # OSAR options
    if (!"osar.method" %in% names(inp)) inp$osar.method <- 'none'
    if (!"osar.trace" %in% names(inp))  inp$osar.trace <- FALSE
    if (!"osar.parallel" %in% names(inp)) inp$osar.parallel <- FALSE
    # Season options
    if (!"seasontype" %in% names(inp)) inp$seasontype <- 1
    if (!"omega" %in% names(inp)) inp$omega <- 2*pi # Annual cycle of noisy oscillator
    # Robust options
    if (!"robflagc" %in% names(inp)) inp$robflagc <- 0
    inp$robflagc <- as.numeric(inp$robflagc)
    if (!"robflagi" %in% names(inp)) inp$robflagi <- rep(0, length(inp$nindexseq))
    inp$robflagi <- as.numeric(inp$robflagi)
    if (!"robflage" %in% names(inp)) inp$robflage <- 0
    inp$robflage <- as.numeric(inp$robflage)
    # Time-varying growth option
    inp <- set.default(inp, 'timevaryinggrowth', FALSE)
    # ASPIC options
    if (!"aspic" %in% names(inp)) inp$aspic <- list()
    if (!"mode" %in% names(inp$aspic)) inp$aspic$mode <- 'FIT'
    if (!"verbosity" %in% names(inp$aspic)) inp$aspic$verbosity <- '102'
    if (!"nboot" %in% names(inp$aspic)) inp$aspic$nboot <- 1000
    if (!"ciperc" %in% names(inp$aspic)) inp$aspic$ciperc <- 95
    if (!"bootlimtime" %in% names(inp$aspic)) inp$aspic$bootlimtime <- 100 # Seconds
    # Meyer & Millar model options
    if (!"meyermillar" %in% names(inp)) inp$meyermillar <- list()
    if (!"n.iter" %in% names(inp$meyermillar)) inp$meyermillar$n.iter <- 225000
    if (!"n.chains" %in% names(inp$meyermillar)) inp$meyermillar$n.chains <- 2
    if (!"burnin" %in% names(inp$meyermillar)) inp$meyermillar$burnin <- 25000
    if (!"thin" %in% names(inp$meyermillar)) inp$meyermillar$thin <- 25
    if (!"cleanup" %in% names(inp$meyermillar)) inp$meyermillar$cleanup <- TRUE
    if (!"bugfn" %in% names(inp$meyermillar)) inp$meyermillar$bugfn <- 'sp.bug'

    # Options for simple model
    if (!"simple" %in% names(inp)) inp$simple <- 0
    if (inp$simple == 1){ # Set parameters for the simple model (catch assumed known, no F process).
        if (!"btype" %in% names(inp)) inp$btype <- 'naive'
        umodtimeC <- unique(inp$timeC %% 1)
        if (length(umodtimeC) != 1){
            stop('When inp$simple = 1, inp$timeC must have a fixed regular time step of 1 year!')
        }
        if (umodtimeC != 0){
            inp$timeC <- floor(inp$timeC)
        }
        for (i in inp$nindexseq){
            umodtimeI <- unique(inp$timeI[[i]] %% 1)
            if (length(umodtimeI) != 1){
                stop('When inp$simple = 1, inp$timeI must have a fixed regular time step of 1 year!')
            }
            if (umodtimeI != 0){
                inp$timeI[[i]] <- floor(inp$timeI[[i]])
            }
        }
        inp$dteuler <- 1
        # Fix parameters that are not used
        inp$phases$logF <- -1
        inp$phases$logsdf <- -1
        #inp$phases$logbeta <- -1
        inp$phases$logsdc <- -1
        # Disable priors
        inp$priors$logbeta <- c(log(1), 1, 0)
        # Prediction not possible
        inp$dtpredc <- NULL
        inp$timepredc <- max(inp$timeC)
        inp$timepredi <- max(unlist(inp$timeI))
    }
    # Biomass related
    if (!"btype" %in% names(inp)) inp$btype <- 'lamperti'
    # MSY type options
    if (!"msytype" %in% names(inp)){
        inp$msytype <- 's'
    } else {
        if (!inp$msytype %in% c('s', 'd')){
            stop('inp$msytype must be either "s" (stochastic) or "d" (deterministic!')
        }
    }

    # Euler time step
    if (!"dteuler" %in% names(inp)){
        inp$dteuler <- 1/16
    }
    if ("dteuler" %in% names(inp)){
        if (inp$dteuler > 1){
            inp$dteuler <- 1
            if(verbose) cat('The dteuler used is not allowed! using inp$dteuler:', inp$dteuler, '\n')
        }
    }

    # - Prediction horizons -
    timeobsall <- sort(c(inp$timeC, inp$timeC + inp$dtc,
                         unlist(inp$timeI),
                         inp$timeE, inp$timeE + inp$dte))

    # - Management variables -
    # Management period, Interval for catch prediction and manstart
    if (any(names(inp) == "maninterval")){
        manflag <- FALSE
        if (length(inp$maninterval) < 2){
            manflag <- TRUE
            if(verbose) warning("Only one of the two times of the management interval specified! The default management interval will be used!\n")
        }else if (inp$maninterval[1] == inp$maninterval[2]){
            manflag <- TRUE
            if(verbose) warning("The times of the specified management interval are equal! The default management interval will be used!\n")
        }else if (inp$maninterval[1] > inp$maninterval[2]){
            inp$maninterval <- sort(inp$maninterval)
            if(verbose) warning("The specified management interval is not increasing! 'inp$maninterval' =",
                                paste0("[",inp$maninterval[1],",",inp$maninterval[2],"]"), "will be used!\n")
        }else{
        if (inp$maninterval[1] < max(timeobsall)){
            if(mancheck){  ## necessary to be able to switch off for HCRs with intermediate periods
                manflag <- TRUE
                if(verbose) warning("The specified management interval starts before the time of the last observation: ", max(timeobsall),"! The default management interval will be used!\n")
            }
        }
        if (abs(diff(inp$maninterval)) < inp$dteuler){
            manflag <- TRUE
            if(verbose) warning("The specified management interval is smaller than the Euler discretisation time step:", inp$dteuler,"! The default management interval will be used!\n")
        }
        if (verbose && any(names(inp) == "timepredc") && inp$timepredc != min(inp$maninterval) && any(names(inp) == "manstart") && inp$manstart != min(inp$maninterval)){
            cat("The arguments 'inp$maninterval', 'inp$timepredc', and 'inp$manstart' are specified and differ. Only 'inp$maninterval' =", paste0("[",inp$maninterval[1],",",inp$maninterval[2],"]"), "will be used! \n")
        }else if (verbose && any(names(inp) == "timepredc") && inp$timepredc != min(inp$maninterval)){
            cat("Both arguments 'inp$maninterval' and 'inp$timepredc' are specified and differ. Only 'inp$maninterval' =", paste0("[",inp$maninterval[1],",",inp$maninterval[2],"]"), "will be used! \n")
        }else if (verbose && any(names(inp) == "manstart") && inp$manstart != min(inp$maninterval))
            cat("Both arguments 'inp$maninterval' and 'inp$manstart' are specified and differ. Only 'inp$maninterval' =", paste0("[",inp$maninterval[1],",",inp$maninterval[2],"]"), "will be used! \n")
        }
        if(manflag){
            manstart <- ceiling(max(timeobsall))
            inp$maninterval <- c(manstart, manstart + 1)
        }
        inp$timepredc <- min(inp$maninterval)
        inp$dtpredc <- abs(diff(inp$maninterval))
        inp$manstart <- min(inp$maninterval)
    }else{
        if(!any(names(inp) == "timepredc") && !any(names(inp) == "dtpredc") && !any(names(inp) == "manstart")){
            ## if no old variables set -> default of maninterval as start of new year after last observation
            manstart <- ceiling(max(timeobsall))
            inp$maninterval <- c(manstart, manstart + 1)
            inp$timepredc <- min(inp$maninterval)
            inp$dtpredc <- abs(diff(inp$maninterval))
            inp$manstart <- manstart
        }else if(any(names(inp) == "timepredc") && any(names(inp) == "dtpredc") && any(names(inp) == "manstart")){
            ## if ALL old variables set -> use old ones and maninterval dependent on them
            ## stop if manstart after timepredc or any before time of last observation
            if(inp$manstart < max(timeobsall)) stop("inp$manstart is set before time of last observation!")
            if(inp$timepredc < max(timeobsall)) stop("inp$timepredc is set before time of last observation!")
            if(verbose && inp$timepredc < inp$manstart) warning("inp$manstart is set after inp$timepredc. This can have unpredictable effects on the management scenarios.")
            inp$maninterval <- c(inp$timepredc, inp$timepredc+inp$dtpredc)
        }else
            ## if any of the old variables missing -> informative error message that all have to be specified
            stop("Specify all three variables: inp$timepredc, inp$dtpredc, and inp$manstart OR use inp$maninterval!")
    }

    # Time point to evaluate model states for management (p states)
    if (any(names(inp) == "maneval")){
        if(verbose && any(names(inp) == "timepredi") && inp$timepredi != inp$maneval)
            cat("Both arguments 'inp$maneval' and 'inp$timepredi' are specified. Only 'inp$maneval' =",
                inp$maneval, "will be used! \n")
        inp$timepredi <- inp$maneval
    }else{
        if(!any(names(inp) == "timepredi")){
            ## default of maneval as end of the management interval
            inp$maneval <- inp$maninterval[2]
            inp$timepredi <- inp$maneval
        }else{
            inp$maneval <- inp$timepredi
        }
    }

    # Interval for effort prediction
    if (!"timeprede" %in% names(inp)){
        if (inp$nobsE > 0){
            inp$timeprede <- ceiling(max(timeobsall))
        } else {
            inp$timeprede <- numeric(0)
        }
    } else {
        if (inp$nobsE > 0){
            if (inp$timeprede < max(inp$timeE + inp$dte)){
                if(verbose) cat('inp$timeprede:', inp$timeprede,
                    ' must be equal to or later than last effort observation: ', max(inp$timeE + inp$dte), '!  \n')
            }
        }
    }
    if (!"dtprede" %in% names(inp)){
        if (inp$nobsE > 0){
            if (length(inp$dte)>0){
                inp$dtprede <- max(inp$dte)
            } else {
                inp$dtprede <- 1
                if(verbose) cat('Assuming a 1 year prediction interval for effort.\n')
            }
        } else {
            inp$dtprede <- numeric(0)
        }
    }
    # This may give a problem if effort data has later time points than catches or index
    if (inp$nobsE > 0 && sum(inp$nobsI) > 0){
        if (max(inp$timeE + inp$dte) > max(unlist(inp$timeI), inp$timeC + inp$dtc)){
            stop('Effort data must overlap temporally with index or catches')
        }
    }

    # Numerical Euler discretisation time used by SDE solver
    # Euler time step
    if (!"dteuler" %in% names(inp)){
        inp$dteuler <- 1/16
    }
    if ("dteuler" %in% names(inp)){
        if (inp$dteuler > 1){
            inp$dteuler <- 1
            if(verbose) cat('The dteuler used is not allowed! using inp$dteuler:', inp$dteuler, '\n')
        }
    }
    if (FALSE){
        alloweddteuler <- 1/2^(6:0)
        if (!inp$dteuler %in% alloweddteuler){ # Check if dteuler is among the alloweddteuler
            ind <- cut(inp$dteuler, alloweddteuler, right=FALSE, labels=FALSE)
            if (is.na(ind)){
                if (inp$dteuler > max(alloweddteuler)){
                    inp$dteuler <- max(alloweddteuler)
                }
                if (inp$dteuler < min(alloweddteuler)){
                    inp$dteuler <- min(alloweddteuler)
                }
            } else {
                inp$dteuler <- alloweddteuler[ind]
            }
            if(verbose) cat('The dteuler used is not allowed! using inp$dteuler:', inp$dteuler, '\n')
        }
    }
    # Euler types:
    # hard: time discretisation is equidistant with step length = dteuler. Observations are assigned to intervals
    # soft: time discretisation is equidistant with step length = dteuler, but with time points of observations inserted such that they can be assigned accurately to a time point instead of an interval.
    # Include dtc because a catch observation at time t includes information in the interval [t; t+dtc[
    if (!"eulertype" %in% names(inp)){
        inp$eulertype <- 'hard'
    }
    if (!"start.in.first.data.point" %in% names(inp)){
        inp$start.in.first.data.point <- TRUE
    }
    if ("eulertype" %in% names(inp)){
        if (inp$eulertype == 'hard'){
            # Hard Euler discretisation
            if (inp$start.in.first.data.point){
                time <- seq(min(timeobsall), max(inp$timepredi, inp$timepredc+inp$dtpredc),
                            by=inp$dteuler)
            } else {
                # Here we take floor of time of first data point
                # This can sometimes cause problems when estimating the random effect
                # prior to the first data point, because of no data support
                time <- seq(floor(min(timeobsall)), max(inp$timepredi, inp$timepredc+inp$dtpredc),
                            by=inp$dteuler)
            }
            inp$time <- time
        }
        if (inp$eulertype == 'soft'){
            # Include times of observations (including dtc)
            time <- seq(ceiling(min(timeobsall)), max(inp$timepredi, inp$timepredc+inp$dtpredc), by=inp$dteuler)
            inp$time <- sort(unique(c(timeobsall, time)))
        }
        if (!inp$eulertype %in% c('soft', 'hard'))
            stop('inp$eulertype must be either "soft" or "hard"!')
    }
    # Calculate time steps
    inp$dt <- c(diff(inp$time), inp$dteuler)
    inp$ns <- length(inp$time)

    # -- DERIVED VARIABLES --
    inp$timerange <- range(timeobsall)
    inp$indlastobs <- cut(max(c(inp$timeC + inp$dtc - inp$dteuler, unlist(inp$timeI),
                                inp$timeE + inp$dte - inp$dteuler)),
                          inp$time, right=FALSE, labels=FALSE)
    inp$indest <- which(inp$time <= inp$timerange[2])
    inp$indpred <- which(inp$time >= inp$timerange[2])
    inp$indCpred <- which(inp$time >= max(inp$timeC + inp$dtc))
    # Redefine management variables in case they do not match any dteuler time steps
    if(!any(inp$time == inp$maninterval[1]) || !any(inp$time == inp$maninterval[2])){
        indmanint <- match.times(inp$maninterval, inp$time)
        inp$maninterval <- inp$time[indmanint]
        inp$manstart <- min(inp$maninterval)
        inp$timepredc <- inp$manstart
        inp$dtpredc <- diff(inp$maninterval)
        if(verbose) cat("Specified management interval does not match to any time step. Management variables are overwritten and 'inp$maninterval' =",paste0("[",inp$maninterval[1],", ",inp$maninterval[2],"]"), "will be used!\n")
    }
    if(!any(inp$time == inp$maneval)){
        indmaneval <- match.times(inp$maneval, inp$time)
        inp$maneval <- inp$time[indmaneval]
        inp$timepredi <- inp$maneval
        if(verbose) cat("Specified management evaluation time does not match to any time step. 'inp$maneval' =",inp$maneval, "will be used!\n")
    }

    if (!"ffac" %in% names(inp)) inp$ffac <- 1
    if ("ffac" %in% names(inp)){
        if (!is.numeric(inp$ffac)){
            if(verbose) cat('Warning: ffac not numeric, which is not allowed, setting ffac = 1. \n')
            inp$ffac <- 1
        }
        if (inp$ffac < 0){
            if(verbose) cat('Warning: ffac < 0, which is not allowed, setting ffac = 0. \n')
            inp$ffac <- 0
        }
    }
    if (!"fcon" %in% names(inp)) inp$fcon <- 0
    if ("fcon" %in% names(inp)){
        if (!is.numeric(inp$fcon)){
            if(verbose) cat('Warning: fcon not numeric, which is not allowed, setting fcon = 0. \n')
            inp$fcon <- 0
        }
        if (inp$fcon < 0){
            if(verbose) cat('Warning: fcon < 0, which is not allowed, setting fcon = 0. \n')
            inp$fcon <- 0
        }
    }
    if (!"ffacvec" %in% names(inp)){
        inp <- make.ffacvec(inp, inp$ffac)
    }
    if (length(inp$ffacvec) != inp$ns){
        if(verbose) warning('Wrong length of inp$ffacvec: ', length(inp$ffacvec),
                            '. Should be equal to inp$ns: ', inp$ns,
                            '. Resetting ffacvec.')
        inp <- make.ffacvec(inp, inp$ffac)
    }
    if (!"fconvec" %in% names(inp)){
        inp$fconvec <- numeric(inp$ns)
        # -1 in indpred because 1 is for plotting
        inp$fconvec[inp$indpred[-1]] <- inp$fcon + 1e-8 # Add small to avoid taking log of 0
    }
    if (length(inp$fconvec) != inp$ns){
        if(verbose) warning('Wrong length of inp$fconvec: ', length(inp$fconvec),
                            '. Should be equal to inp$ns: ', inp$ns,
                            '. Resetting fconvec.')
        inp$fconvec[inp$indpred[-1]] <- inp$fcon + 1e-8 # Add small to avoid taking log of 0
    }

    # Seasons
    if (!"nseasons" %in% names(inp)){
        expnseasons <- 1/min(inp$dtc)
        if (expnseasons >= 4){
            inp$nseasons <- 4
        } else {
            inp$nseasons <- 1
        }
    }
    if ("nseasons" %in% names(inp)){
       if (!inp$nseasons %in% c(1, 2, 4)){
           stop('inp$nseasons (=', inp$nseasons, ') must be either 1, 2 or 4.')
       }
    }
    if (inp$nseasons == 1) inp$seasontype <- 0 # seasontype = 0 means seasons are disabled.
    # Calculate seasonal spline
    if ("splineorder" %in% names(inp)){
        if (inp$nseasons < 4 & inp$splineorder > 2){
            inp$splineorder <- 2
        }
    } else {
        inp$splineorder <- ifelse(inp$nseasons < 4, 2, 3)
    }
    inp$splinemat <- make.splinemat(inp$nseasons, inp$splineorder, dtfine=inp$dteuler)
    inp$splinematfine <- make.splinemat(inp$nseasons, inp$splineorder, dtfine=1/100)
    inp$seasonindex <- 1/inp$dteuler*(inp$time %% 1)
    inp$seasons <- rep(0, inp$ns)
    inp$seasonindex2 <- rep(1:inp$ns,each=inp$nseasons,length.out=inp$ns)
    for (i in 1:inp$nseasons){
        frac <- 1/inp$nseasons
        modtime <- inp$time %% 1
        inds <- which(modtime>=((i-1)*frac) & modtime<(i*frac))
        inp$seasons[inds] <- i
    }
    # ic is the indices of inp$time to which catch observations correspond
    #if (length(inp$dtc) > 0){
    #    dtcpred <- min(inp$dtc)
    #} else {
    #    dtcpred <- 1
    #}
    dtcpred <- inp$dtpredc
    if(inp$timepredc > max(inp$timeC)){
        inp$timeCpred <- unique(c(inp$timeC, (seq(tail(inp$timeC,1), inp$timepredc, by=tail(inp$dtc,1))), inp$timepredc))
    }else{ ## for constant catch scenarios which add obs to inp
        inp$timeCpred <- unique(inp$timeC)
    }
    inp$nobsCp <- length(inp$timeCpred)
    if( inp$nobsCp > inp$nobsC ){
        ## inp$dtcp <- c(inp$dtc, rep(tail(inp$dtc,1), inp$nobsCp-inp$nobsC-1),dtcpred) ##wrong if int period 0.5yr
        inp$dtcp <- c(diff(inp$timeCpred), dtcpred)
    }else inp$dtcp <- inp$dtc

    inp$ic <- cut(inp$timeCpred, inp$time, right=FALSE, labels=FALSE)

    # nc is number of states to integrate a catch observation over
    inp$nc <- rep(0, inp$nobsCp)
    for (i in 1:inp$nobsCp){
        inp$nc[i] <- sum(inp$time >= inp$timeCpred[i] & inp$time < (inp$timeCpred[i]+inp$dtcp[i]))
    }
    if (any(inp$nc == 0)){
        stop('Current inp$dteuler is too large to accommodate some catch intervals. Make inp$dteuler smaller!')
    }
    # ie is the indices of inp$time to which effort observations correspond
    if (length(inp$dte) > 0){
        dtepred <- min(inp$dte)
    } else {
        dtepred <- 1
    }
    if (inp$nobsE > 0){
        inp$timeEpred <- unique(c(inp$timeE, (seq(tail(inp$timeE,1), inp$timeprede, by=dtepred))))
    } else {
        inp$timeEpred <- numeric(0)
    }
    inp$nobsEp <- length(inp$timeEpred)
    inp$dtep <- c(inp$dte, rep(dtepred, inp$nobsEp-inp$nobsE))
    inp$ie <- cut(inp$timeEpred, inp$time, right=FALSE, labels=FALSE)
    # ne is number of states to integrate an effort observation over
    inp$ne <- rep(0, inp$nobsEp)
    if (inp$nobsE > 0){
        for (i in 1:inp$nobsEp){
            inp$ne[i] <- sum(inp$time >= inp$timeEpred[i]
                             & inp$time < (inp$timeEpred[i]+inp$dtep[i]))
        }
    }
    if (any(inp$ne == 0)){
        stop('Current inp$dteuler is too large to accommodate some effort intervals. Make inp$dteuler smaller!')
    }
    # ii is the indices of inp$time to which index observations correspond
    inp$ii <- list()
    for (i in inp$nindexseq){
        # Find indices in inp$time to which observations correspond
        # cut finds indices based on intervals in inp$time
        #inp$ii[[i]] <- cut(inp$timeI[[i]], inp$time, right=FALSE, labels=FALSE)
        # This approach uses the index of inp$time with the smallest temporal difference
        # to the observations, this difference will typically be zero.
        inp$ii[[i]] <- match.times(inp$timeI[[i]], inp$time)
    }
    # Translate index observations from a list to a vector
    inp$obsIin <- unlist(inp$obsI)
    inp$stdevfacIin <- unlist(inp$stdevfacI)
    if (is.null(inp$stdevfacIin)) inp$stdevfacIin <- numeric(0) # If zero index observations
    inp$iiin <- unlist(inp$ii)
    if (is.null(inp$iiin)) inp$iiin <- numeric(0)
    inp$iqin <- rep(inp$mapq, times=inp$nobsI)
    if (is.null(inp$iqin)) inp$iqin <- numeric(0)
    inp$isdiin <- rep(inp$mapsdi, times=inp$nobsI)
    if (is.null(inp$isdiin)) inp$isdiin <- numeric(0)
    # Add helper variable such that predicted catch can be calculated using small euler steps
    # Need to include timerange[2] and exclude timerange[2]+dtpred because the catch at t is acummulated over t to t+dtc.
    inp$dtpredcinds <- which(inp$time >= inp$timepredc & inp$time < (inp$timepredc+inp$dtpredc))
    inp$dtpredcnsteps <- length(inp$dtpredcinds)
    #inp$dtprediind <- cut(inp$timepredi, inp$time, right=FALSE, labels=FALSE)
    inp$dtprediind <- match.times(inp$timepredi, inp$time)
    inp$dtpredeinds <- which(inp$time >= inp$timeprede & inp$time < (inp$timeprede+inp$dtprede))
    inp$dtpredensteps <- length(inp$dtpredeinds)

    # - Sort observations in time and store in one vector -
    timeobsseen <- c(inp$timeC + inp$dtc - 1e-4, unlist(inp$timeI), unlist(inp$timeE) + unlist(inp$dte) - 1e-5) # Add dtc to timeC because the total catch is first "seen" at the end of the given catch interval (typically year or quarter), similarly for quarter. By arbitrary convention catches are assumed to be "seen" before effort although they should be observed at the same time.
    srt <- sort(timeobsseen, index=TRUE)
    timeobs <- c(inp$timeC, unlist(inp$timeI), unlist(inp$timeE))
    timeobssrt <- timeobs[srt$ix]
    obs <- log(c(inp$obsC, unlist(inp$obsI), unlist(inp$obsE)))
    obsid <- c(inp$obsidC, unlist(inp$obsidI), unlist(inp$obsidE))
    inp$obssrt <- obs[srt$ix]
    inp$timeobssrt <- timeobs[srt$ix]
    inp$obsidsrt <- obsid[srt$ix]
    inp$isc <- match(1:inp$nobsC, srt$ix)
    if (sum(inp$nobsI)>0){
        inp$isi <- match((inp$nobsC+1):(inp$nobsC+sum(inp$nobsI)), srt$ix)
    } else {
        inp$isi <- numeric(0)
    }
    if (sum(inp$nobsE)>0){
        inp$ise <- match((inp$nobsC+sum(inp$nobsI)+1):(inp$nobsC+sum(inp$nobsI)+sum(inp$nobsE)), srt$ix)
    } else {
        inp$ise <- numeric(0)
    }
    if (sum(inp$nobsI) != length(inp$isi)){
        warning('Mismatch between length(inp$isi) ', length(inp$isi),
                ' and sum(inp$nobsI) ', sum(inp$nobsI), '.')
    }
    if (sum(inp$nobsE) != length(inp$ise)){
        warning('Mismatch between length(inp$ise) ', length(inp$ise),
                ' and sum(inp$nobsE) ', sum(inp$nobsE), '.')
    }
    # Condition on the first year of data.
    inp$osar.conditional <- which(inp$timeobssrt < (inp$time[1] + 1))
    inp$osar.subset <- setdiff(1:length(inp$obssrt), inp$osar.conditional)

    # -- COVARIATES --
    inp$logmcovflag <- FALSE
    inp <- set.default(inp, 'logmcovspar', 0.5)
    if ('logmcovariate' %in% names(inp)){
        if (!'logmcovariatetime' %in% names(inp)){
            stop('inp$logmcovariatetime unspecified but required!')
        }
        if ('logmcovariatetime' %in% names(inp)){
            if (length(inp$logmcovariatetime) != length(inp$logmcovariate)){
                stop('length(inp$logmcovariatetime) != length(inp$logmcovariate), cannot continue!')
            }
        }
        # Check for NAs and remove
        nainds <- which(is.na(inp$logmcovariate))
        if (length(nainds) > 0){
            inp$logmcovariate <- inp$logmcovariate[-nainds]
            inp$logmcovariatetime <- inp$logmcovariatetime[-nainds]
        }
        dat <- data.frame(x=inp$logmcovariatetime, y=inp$logmcovariate)
        smoocov <- smooth.spline(dat$x, dat$y, spar=inp$logmcovspar)
        covpred <- predict(smoocov, x=inp$time)
        #plot(covpred$x, covpred$y, typ='l')
        #points(inp$logmcovariatetime, inp$logmcovariate)
        inp$logmcovariatein <- covpred$y
        inp$logmcovflag <- TRUE
    }
    # Fill in dummy defaults if unspecified
    #inp <- set.default(inp,'logmcovariate', rep(0, inp$nobsC))
    #inp <- set.default(inp, 'logmcovariatetime', 1:inp$nobsC)
    inp <- set.default(inp, 'logmcovariatein', rep(0, inp$ns))
    if(length(inp$logmcovariatein) != inp$ns){
        if(verbose) warning('Wrong length of inp$logmcovariatein: ', length(inp$logmcovariatein),
                            '. Should be equal to inp$ns: ', inp$ns,
                            '. Resetting logmcovariatein.')
        inp$logmcovariatein <- NULL
        inp <- set.default(inp, 'logmcovariatein', rep(0, inp$ns))
    }


    # -- MODEL PARAMETERS --
    # Default values
    lognflag <- TRUE
    logalphaflag <- TRUE
    logbetaflag <- TRUE
    logn <- log(2)
    logalpha <- log(1)
    logbeta <- log(1)
    if (!"logn" %in% names(inp$ini)){
        inp$ini$logn <- logn
    } else {
        if (inp$ini$logn == 0){
            stop('Initial value for logn == 0, that is not valid!')
        }
    }
    # Calculate gamma from n
    n <- exp(inp$ini$logn)
    inp$ini$gamma <- calc.gamma(n)
    # logK
    if (!'logK' %in% names(inp$ini)){
        inp$ini$logK <- log(4*max(inp$obsC))
    }
    # logr
    if ('logr' %in% names(inp$ini) & 'logm' %in% names(inp$ini)){
        inp$ini$logr <- NULL # If both r and m are specified use m and discard r
    }

    # find number of regimes for 'm'
    if(!'MSYregime' %in% names(inp)){
        inp$MSYregime <- factor(rep(1,length(inp$time)))
    } else if(length(inp$MSYregime) < length(inp$time)){ # manage changes number of time steps!
        if(verbose) warning('Wrong length of inp$MSYregime: ', length(inp$MSYregime),
                            '. Should be equal to inp$ns: ', inp$ns,
                            '. Resetting MSYregime.')
        inp$MSYregime <- factor(c(inp$MSYregime,
                                  rep( tail(inp$MSYregime,1), length(inp$time)-length(inp$MSYregime))))
    }
    if(inp$timevaryinggrowth && nlevels(inp$MSYregime)>1)
        stop("'timevaryinggrowth' and multiple MSYregimes cannot be used at the same time")
    inp$noms<-nlevels(inp$MSYregime)
    inp$ir<-as.numeric(inp$MSYregime)

    if (!'logitSARphi' %in% names(inp$ini)) inp$ini$logitSARphi <- 0
    if (!'logSdSAR' %in% names(inp$ini)) inp$ini$logSdSAR <- -2


    if (!'logr' %in% names(inp$ini)){
        if (!'logm' %in% names(inp$ini) | length(inp$ini$logm)!=inp$noms){
            inp$ini$logm <- rep(unname(log(guess.m(inp))), inp$noms)
        }
        r <- exp(inp$ini$logm)/exp(inp$ini$logK) * n^(n/(n-1))
        inp$ini$logr <- log(r)
        #r <- inp$ini$gamma * exp(inp$ini$logm) / exp(inp$ini$logK) # n > 1 # rold
        #if (n>1){
        #    inp$ini$logr <- log(r)
        #} else {
        #    inp$ini$logr <- log(-r)
        #}
    }
    if(length(inp$ini$logm)!=inp$noms) inp$ini$logm <- rep(inp$ini$logm,noms)


    if ('logr' %in% names(inp$ini)){
        nr <- length(inp$ini$logr)
        if (!'ir' %in% names(inp) | nr == 1){
            inp$ir <- rep(0, inp$ns)
            for (i in 1:nr){
                frac <- 1/nr
                modtime <- inp$time %% 1
                inds <- which(modtime>=((i-1)*frac) & modtime<(i*frac))
                inp$ir[inds] <- i
            }
        } else {
            if (length(unique(inp$ir)) != nr){
                stop('Mismatch between specified inp$ir and inp$ini$logr!')
            }
            nir <- length(inp$ir)
            if (nir != inp$ns){
                if (nir == inp$nobsC){ # Assume that inp$ir fits with inp$timeC
                    ir <- rep(0, inp$ns)
                    for (i in 1:nir){
                        inds <- which(inp$time >= inp$timeC[i] & inp$time < (inp$timeC[i]+inp$dtc[i]))
                        ir[inds] <- inp$ir[i]
                    }
                    inds <- which(inp$time >= inp$timeC[nir])
                    ir[inds] <- inp$ir[nir]
                    inp$ir <- ir
                } else {
                    if (nir == inp$nobsI[[1]]){ # Assume that inp$ir fits with inp$timeI[[1]]
                        ir <- rep(0, inp$ns)
                        for (i in 2:nir){
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

    check.mapped.ini <- function(inp, nam, nnam){
        if (nam %in% names(inp$ini)){
            if (length(inp$ini[[nam]]) != inp[[nnam]]){ # nq is given by mapq
                if (length(inp$ini[[nam]]) == 1){
                    inp$ini[[nam]] <- rep(inp$ini[[nam]], inp[[nnam]])
                } else {
                    stop('The length of ', nam, ' in inp$ini (', length(inp$ini[[nam]]),
                         ') does not fit with the number of parameters to be estimated (',
                         inp[[nam]], ').')
                }
            }
        }
        return(inp)
    }

    if (length(inp$obsI) == 0){
        logmaxI <- 0
    } else {
        logmaxI <- log(max(inp$obsI[[1]]))
    }
    #logmaxI <- ifelse(length(inp$obsI)==0, 0, log(max(inp$obsI[[1]])))
    if (!'logq' %in% names(inp$ini)) inp$ini$logq <- logmaxI - inp$ini$logK
    if (sum(inp$nobsI)>0) inp <- check.mapped.ini(inp, 'logq', 'nq')
    if (inp$nobsE > 0){
        logmaxE <- log(max(inp$obsE))
    } else {
        logmaxE <- 0
    }
    #logmaxE <- ifelse(length(inp$obsE)==0, 0, log(max(inp$obsE[[1]])))
    if (!'logqf' %in% names(inp$ini)) inp$ini$logqf <- inp$ini$logr - logmaxE
    #if (sum(inp$nobsE)>0) inp <- check.mapped.ini(inp, 'logqf', 'nqf')
    inp$isdf <- rep(1, inp$ns)
    if (!is.null(inp$sdfsplityear)){
        inds <- which(inp$time >= inp$sdfsplityear)
        inp$isdf[inds] <- 2
        inp$nsdf <- length(unique(inp$isdf))
    } else {
        inp$nsdf <- 1
    }

    if (!'logsdf' %in% names(inp$ini)){
        inp$ini$logsdf <- log(rep(0.2, inp$nsdf))
    }
    if (!'logsdu' %in% names(inp$ini)) inp$ini$logsdu <- log(0.1)
    if (!'logsdb' %in% names(inp$ini)) inp$ini$logsdb <- log(0.2)
    inp$ini <- set.default(inp$ini, 'logsdm', log(1e-8))
    if (!'logsdc' %in% names(inp$ini)) inp$ini$logsdc <- log(0.2)
    if (!'logsdi' %in% names(inp$ini)) inp$ini$logsdi <- log(0.2)
    if (sum(inp$nobsI)>0) inp <- check.mapped.ini(inp, 'logsdi', 'nsdi')
    if (!'logsde' %in% names(inp$ini)) inp$ini$logsde <- log(0.2)
    #if (sum(inp$nobsE)>0) inp <- check.mapped.ini(inp, 'logsde', 'nsde')
    if (!"logalpha" %in% names(inp$ini)) inp$ini$logalpha <- logalpha
    if (sum(inp$nobsI)>0) inp <- check.mapped.ini(inp, 'logalpha', 'nsdi')
    if (!"logbeta" %in% names(inp$ini))  inp$ini$logbeta <- logbeta

    if (!"logm" %in% names(inp$ini)){
        gamma <- inp$ini$gamma
        r <- exp(inp$ini$logr)
        K <- exp(inp$ini$logK)
        n <- exp(inp$ini$logn)
        m <- rep( r * K / (n^(n/(n-1))), inp$noms)
    }
    logm <- inp$ini$logm # Store this to be able to set logmre later
    if (inp$timevaryinggrowth){
        inp$ini$logm <- log(1)
    }
    # Fill in unspecified (more rarely user defined) model parameter values
    inp$ini <- set.default(inp$ini, 'logpsi', log(1e-8))
    inp$ini <- set.default(inp$ini, 'mu', 0)
    if (!"loglambda" %in% names(inp$ini)) inp$ini$loglambda <- log(0.1)
    if (!"logdelta" %in% names(inp$ini)) inp$ini$logdelta <- log(1e-8) # Strength of mean reversion of OU for F
    if (!"logeta" %in% names(inp$ini)) inp$ini$logeta <- log(0.2) # Mean of OU for F
    if ("logphi" %in% names(inp$ini)){
        if (length(inp$ini$logphi)+1 != dim(inp$splinemat)[2]){
            if(verbose) cat('Mismatch between length of ini$logphi and number of columns of splinemat! removing prespecified ini$logphi and setting default.\n')
            inp$ini$logphi <- NULL
        }
    }
    if (!"logphi" %in% names(inp$ini)) inp$ini$logphi <- rep(0, inp$nseasons-1)
    if (!"logitpp" %in% names(inp$ini)) inp$ini$logitpp <- log(0.95/(1-0.95))
    if (!"logp1robfac" %in% names(inp$ini)) inp$ini$logp1robfac <- log(15-1)
    if (!"logbkfrac" %in% names(inp$ini)) inp$ini$logbkfrac <- log(0.8)
    if (!"logF" %in% names(inp$ini)){
        inp$ini$logF <- rep(log(0.2) + inp$ini$logr[1], inp$ns)
    } else if (length(inp$ini$logF) > inp$ns){
        if(verbose) warning('Wrong length of inp$ini$logF: ', length(inp$ini$logF),
                            '. Should be equal to inp$ns: ', inp$ns,
                            '. Setting length of logF equal to inp$ns (removing beyond inp$ns).')
        inp$ini$logF <- inp$ini$logF[1:inp$ns]
    }else if (length(inp$ini$logF) < inp$ns){
        if(verbose) warning('Wrong length of inp$ini$logF: ', length(inp$ini$logF),
                            '. Should be equal to inp$ns: ', inp$ns,
                            '. Resetting logF.')
        inp$ini$logF <- rep(log(0.2) + inp$ini$logr[1], inp$ns)
    }
    if ("logu" %in% names(inp$ini)){
        if (dim(inp$ini$logu)[1] != 2*length(inp$ini$logsdu) || dim(inp$ini$logu)[2] != inp$ns){
            if(verbose) warning('Wrong dimension of inp$ini$logu: ', dim(inp$ini$logu)[1], 'x',
                                dim(inp$ini$logu)[2], '. Should be equal to 2*length(inp$ini$logsdu) x inp$ns: ',
                                2*length(inp$ini$logsdu), 'x', inp$ns,'. Filling with log(1).')
            inp$ini$logu <- NULL
        }
    }
    if (!"logu" %in% names(inp$ini)){
        inp$ini$logu <- matrix(log(1)+1e-3, 2*length(inp$ini$logsdu), inp$ns)
    }
    if (!"logB" %in% names(inp$ini)){
        inp$ini$logB <- rep(inp$ini$logK + log(0.5), inp$ns)
    } else if (length(inp$ini$logB) > inp$ns){
        if(verbose) warning('Wrong length of inp$ini$logB: ', length(inp$ini$logB),
                            '. Should be equal to inp$ns: ', inp$ns,
                            '. Setting length of logB equal to inp$ns (removing beyond inp$ns).')
        inp$ini$logB <- inp$ini$logB[1:inp$ns]
    } else if (length(inp$ini$logB) < inp$ns){
        if(verbose) warning('Wrong length of inp$ini$logB: ', length(inp$ini$logB),
                            '. Should be equal to inp$ns: ', inp$ns,
                            '. Resetting logB.')
        inp$ini$logB <- rep(inp$ini$logK + log(0.5), inp$ns)
    }
    if (!"logmre" %in% names(inp$ini)){
        inp$ini$logmre <- rep(log(1), inp$ns)
    } else if (length(inp$ini$logmre) > inp$ns){
        if(verbose) warning('Wrong length of inp$ini$logmre: ', length(inp$ini$logmre),
                            '. Should be equal to inp$ns: ', inp$ns,
                            '. Setting length of logmre equal to inp$ns (removing beyond inp$ns).')
        inp$ini$logmre <- inp$ini$logmre[1:inp$ns]
    } else if (length(inp$ini$logmre) < inp$ns){
        if(verbose) warning('Wrong length of inp$ini$logmre: ', length(inp$ini$logmre),
                            '. Should be equal to inp$ns: ', inp$ns,
                            '. Resetting logmre.')
        inp$ini$logmre <- rep(log(1), inp$ns)
    }

    #if ("logmre" %in% names(inp$ini)){
    #    inp$ini$logmre <- check.mat(inp$ini$logmre, c(inp$nstocks, inp$ns), 'inp$ini$logmre')
    #}
    inp$ini$SARvec <- rep(0, max(inp$seasonindex2))

    ## reporting
    if(!"reportmode" %in% names(inp)) inp$reportmode <- 0

    ## timerange of original observations (required for intermediate catch)
    if(!"lastCatchObs" %in% names(inp)) inp$lastCatchObs <- max(inp$timeC + inp$dtc)
    if(!"timerangeObs" %in% names(inp)) inp$timerangeObs <- inp$timerange

    ## fractiles used for management
    checkandadd("manFractiles", list(catch = 0.5, bbmsy = 0.5, ffmsy = 0.5), "list")

    ## biomass breakpoint used for management
    if(!"manBreakpointB" %in% names(inp)) inp$manBreakpointB <- 0

    ## biomass safeguard used for management
    checkandadd("manSafeguardB", list(limitB = 0, prob = 0.95), "list")

    ## cfac, ffac, bfac
    checkandadd("manfacs", list(cfac = NULL, ffac = NULL, bfac = NULL), "list")

    ## Simulate random effects?
    if(!"sim.random.effects" %in% names(inp)) inp$sim.random.effects <- TRUE

    ## Simulate using fitted object (env$last.par)?
    if(!"sim.fit" %in% names(inp)) inp$sim.fit <- TRUE

    ##
    if(!"iuse" %in% names(inp) || length(inp$iuse) != length(unlist(inp$obsI)))
        inp$iuse <- rep(TRUE, length(unlist(inp$obsI)))



    # Reorder parameter list
    inp$parlist <- list(logm=inp$ini$logm,
                        mu=inp$ini$mu,
                        logK=inp$ini$logK,
                        logq=inp$ini$logq,
                        logqf=inp$ini$logqf,
                        logn=inp$ini$logn,
                        logsdb=inp$ini$logsdb,
                        logsdu=inp$ini$logsdu,
                        logsdf=inp$ini$logsdf,
                        logsdi=inp$ini$logsdi,
                        logsde=inp$ini$logsde,
                        logsdc=inp$ini$logsdc,
                        logsdm=inp$ini$logsdm,
                        logpsi=inp$ini$logpsi,
                        logphi=inp$ini$logphi,
                        loglambda=inp$ini$loglambda,
                        logdelta=inp$ini$logdelta,
                        logeta=inp$ini$logeta,
                        logitpp=inp$ini$logitpp,
                        logp1robfac=inp$ini$logp1robfac,
                        logF=inp$ini$logF,
                        logu=inp$ini$logu,
                        logB=inp$ini$logB,
                        logmre=inp$ini$logmre,
                        SARvec=inp$ini$SARvec,
                        logitSARphi=inp$ini$logitSARphi,
                        logSdSAR=inp$ini$logSdSAR)


    # -- PRIORS --
    # Priors are assumed Gaussian and specified in a vector of length 3:
    # c(log(mean), stdev in log, useflag).
    # log(mean): log of the mean of the prior distribution.
    # stdev in log: standard deviation of the prior distribution in log domain.
    # useflag: if 1 then the prior is used, if 0 it is not used. Default is 0.
    check.prior <- function(priorvec, priorname){
        #priorvec <- priors[[priorname]]
        if (priorname %in% repriors){ # RE priors
            if (length(priorvec) < 3){
                priorvec <- rep(0, 5)
                warning('Invalid prior length specified for', priorname,
                        ', must be 3 (without useflag or 4 (with useflag). Not using this prior.')
            }
            if (length(priorvec) == 3){
                #warning('Length of ', priorname, ' is 3. Proceeding assuming useflag has not been specified.')
                priorvec <- c(priorvec[1:2], 1, priorvec[3])
            }
            if (length(priorvec) == 4){
                ib <- match(priorvec[4], inp$time)
                if (is.na(ib)){
                    ib <- 0
                    priorvec[3] <- 0
                    warning('Year for prior on ', priorname, ' (', priorvec[3], ') did not match times where this RE is estimated. Not using this prior. To fix this use a year where an observation is available.')
                }
                priorvec <- c(priorvec, ib) # Add index in time vec to which this year corresponds
            }
        } else { # FE priors
            if (!length(priorvec) %in% 2:3){
                priorvec <- rep(0, 3)
                warning('Invalid prior length specified for', priorname,
                        ', must be 2 (without useflag or 3 (with useflag). Not using this prior.')
            }
            if (length(priorvec) == 2){
                #warning('Length of ', priorname, ' is 2. Proceeding assuming useflag has not been specified.')
                priorvec <- c(priorvec, 1)
            }
        }
        if (priorvec[3] == 1){
            # Check st dev
            if (priorvec[2] <= 0){
                warning('Invalid standard deviation specified in prior for', priorname,
                        '(must be > 0). Not using this prior.')
                priorvec[3] <- 0
            }
        }
        return(priorvec)
    }
    possiblepriors <- c('logn', 'logalpha', 'logbeta', 'logr', 'logK', 'logm', 'logq',
                        'logqf', 'logbkfrac', 'logB', 'logF', 'logBBmsy',
                        'logFFmsy', 'logsdb', 'logsdf',
                        'logsdi', 'logsde','logsdc',
                        'logsdm', 'logpsi', 'mu', 'BmsyB0','logngamma')
    repriors <- c('logB', 'logF', 'logBBmsy', 'logFFmsy')
    matrixpriors <- c('logsdi','logq')
    npossiblepriors <- length(possiblepriors)
    if (!"priors" %in% names(inp)){
        inp$priors <- list()
    }

    # Default prior uncertainties
    small <- 1e-3
    wide <- 2 # Value suggested by Casper
    lognsd <- wide
    logalphasd <- wide
    logbetasd <- wide
    if (!'logn' %in% names(inp$priors)){
        inp$priors$logn <- c(logn, lognsd)
    }
    if (!'logalpha' %in% names(inp$priors)){
        inp$priors$logalpha <- c(logalpha, logalphasd)
    }
    if (!'logbeta' %in% names(inp$priors)){
        inp$priors$logbeta <- c(logbeta, logbetasd)
    }
    if (inp$timevaryinggrowth){
        inp$priors <- set.default(inp$priors, 'logsdm', c(log(0.2), wide))
        inp$priors <- set.default(inp$priors, 'logpsi', c(log(0.01), wide))
    }
    # Remaining priors, set to something, but will not be used
    if ("priors" %in% names(inp)){
        # Remove wrong priors names
        nms <- names(inp$priors)
        inds <- which(is.na(match(nms, possiblepriors)))
        if (length(inds)>0){
            warning('Wrong prior names specified: ', nms[inds])
            inp$priors[inds] <- NULL
        }
        # Check priors
        for (nm in possiblepriors){
            if (!nm %in% names(inp$priors)){
                # Set default prior values. These will not be used
                if (nm %in% repriors){
                    inp$priors[[nm]] <- c(log(1e-4), 0.2, 0, 2000, 0) # RE prior
                } else {
                    if (nm %in% matrixpriors){
                        npar <- length(inp$ini[[nm]])
                        inp$priors[[nm]] <- list()
                        for (i in 1:npar){
                            inp$priors[[nm]][[i]] <- c(log(1e-4), 0.2, 0) # FE priors
                        }
                    } else {
                        inp$priors[[nm]] <- c(log(1e-4), 0.2, 0) # FE priors
                    }
                }
            } else {
                if (nm %in% matrixpriors){
                    if (is(inp$priors[[nm]], 'list')){
                        npar <- length(inp$ini[[nm]])
                        if (length(inp$priors[[nm]]) == npar){
                            for (i in 1:npar){
                                inp$priors[[nm]][[i]] <- check.prior(inp$priors[[nm]][[i]], nm)
                            }
                        } else {
                            stop('Wrong number of priors specified for ', nm)
                        }
                    } else {
                        if (is(inp$priors[[nm]], 'numeric')){
                            warning('The prior specified for ', nm, ' is replicated for each time series.')
                            thisprior <- check.prior(inp$priors[[nm]], nm)
                            npar <- length(inp$ini[[nm]])
                            inp$priors[[nm]] <- rep(list(thisprior), each=npar)
                        } else {
                            stop('Prior specified for ', nm, ' should be a list of length ', npar, ' or a vector.')
                        }
                    }
                } else {
                    inp$priors[[nm]] <- check.prior(inp$priors[[nm]], nm)
                }
            }
        }
    }
    # Translate matrixpriors from lists to matrices
    inp$matrixpriors <- list()
    for (nm in matrixpriors){
        inp$matrixpriors[[nm]] <- do.call(rbind, inp$priors[[nm]])
        rownames(inp$matrixpriors[[nm]]) <- paste0(nm, 1:length(inp$priors[[nm]]))
    }
    npriors <- length(inp$priors)
    inp$priorsuseflags <- numeric(npriors)
    for (i in 1:npriors){
        nm <- names(inp$priors)[i]
        if (nm %in% matrixpriors){
            inp$priorsuseflags[i] <- max(inp$matrixpriors[[nm]][, 3])
        } else {
            inp$priorsuseflags[i] <- inp$priors[[i]][3]
        }
    }


    # ------------------------------
    # Determine fixed parameters
    forcefixpars <- c() # Parameters that are forced to be fixed.
    if (inp$nseasons == 1){
        forcefixpars <- c('logphi', 'logu', 'logsdu', 'loglambda', 'SARvec','logitSARphi','logSdSAR',forcefixpars)
    } else {
        if (inp$seasontype == 1){ # Use spline
            forcefixpars <- c('logu', 'logsdu', 'loglambda','SARvec','logitSARphi','logSdSAR', forcefixpars)
        }
        if (inp$seasontype == 2){ # Use coupled SDEs
            forcefixpars <- c('logphi','SARvec','logitSARphi','logSdSAR', forcefixpars)
        }
        if(inp$seasontype == 3){ # Use spline + AR
            forcefixpars <- c('logu', 'logsdu', 'loglambda', forcefixpars)
        }
    }
    if (inp$robflagc == 0 & all(inp$robflagi == 0) & inp$robflage == 0){
        forcefixpars <- c('logitpp', 'logp1robfac', forcefixpars)
    }
    if (inp$nobsE == 0){
        forcefixpars <- c('logqf', 'logsde', forcefixpars)
    }
    if (sum(inp$nobsI) == 0){
        forcefixpars <- c('logq', 'logsdi', forcefixpars)
    }
    if (inp$effortmodel == 'RW'){
        forcefixpars <- c('logdelta', 'logeta', forcefixpars)
    }
    if (!inp$logmcovflag){
        forcefixpars <- c('mu', forcefixpars)
    }
    if (!inp$timevaryinggrowth){
        forcefixpars <- c('logmre', 'logsdm', 'logpsi', forcefixpars)
    }

    # Determine phases
    if (!"phases" %in% names(inp)){
        inp$phases <- list()
    } else {
        if ("logr" %in% names(inp$phases)){
            inp$phases$logm <- inp$phases$logr
            inp$phases$logr <- NULL
        }
        nms <- names(inp$phases)
        inds <- which(is.na(match(nms, c(names(inp$ini), 'logalpha', 'logbeta'))))
        if (length(inds) > 0){
            stop('phase specified for invalid parameter(s): ', paste(nms[inds], collapse=', '))
        }
    }
    if ("phases" %in% names(inp)){
        for (i in 1:length(forcefixpars)){
            inp$phases[[forcefixpars[i]]] <- -1
        }
    }
    # Assign phase 1 to parameters without a phase
    nms <- names(inp$parlist)
    nnms <- length(nms)
    for (nm in nms){
        if (!nm %in% names(inp$phases)){
            inp$phases[[nm]] <- 1
        }
    }

    nphasepars <- length(inp$phases)
    phasevec <- unlist(inp$phases)
    phases <- unique(unname(phasevec))
    phases <- phases[phases>0] # Don't include phase -1 (fixed value)
    nphases <- length(phases)
    inp$map <- list()
    for (j in 1:nphases){
        inp$map[[j]] <- list() # Map for phase j
        inds <- which(phasevec > j | phasevec == -1)
        for (i in inds){
            parnam <- names(phasevec)[i]
            if (parnam %in% names(inp$parlist)){
                phase <- inp$phases[[parnam]]
                inp$map[[j]][[parnam]] <- factor(rep(NA, length(inp$parlist[[parnam]])))
            } else {
                warning('Phase specified for an invalid parameter: ', parnam)
                inp$phases[[parnam]] <- NULL # Remove invalid parameter
            }
        }
    }
    inpphases <- unlist(inp$phases)
    inpphases <- inpphases[inpphases > 0] # Don't include phase -1 (fixed value)
    inp$nphases <- length(unique(inpphases))

    # SET DEFAULT PARAMETER RANGES USED IN ROBUSTNESS CHECK
    if (!"minfac" %in% names(inp)){
        inp$minfac <- 0.1
    }
    if (!"maxfac" %in% names(inp)){
        inp$maxfac <- 1 / inp$minfac
    }
    logfacs <- log(c(inp$minfac, inp$maxfac))
    defaultranges <- c('logn', 'logK', 'logm', 'logq', 'logqf',
                       'logsdb', 'logsdf', 'logsdi', 'logsdc', 'logsde')
    inp$ranges <- list()
    # Set default ranges (that aren't set)
    for (nm in defaultranges){
        if (!nm %in% names(inp$ranges)){
            inp$ranges[[nm]] <- outer(inp$ini[[nm]], logfacs, FUN='+')
        }
    }

    # Set class
    if (!is.null(inp)){
        class(inp) <- "spictcls"
    }
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
