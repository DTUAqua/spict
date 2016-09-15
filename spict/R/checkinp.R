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
#' @details Fills in defalut values if missing.
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
#'  \item{"inp$dtef"}{ Time interval for effort observations. For annual effort inp$dtef=1, for quarterly effort inp$dtef=0.25. Default: min(diff(inp$timeE)). }
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
#' useflag: if 1 then the prior is used, if 0 it is not used. Default is 0.
#' To list parameters to which priors can be applied run list.possible.priors().
#' Example: intrinsic growth rate of 0.8
#'  inp$priors$logr <- c(log(0.8), 0.1)
#'  inp$priors$logr <- c(log(0.8), 0.1, 1) # This includes the optional useflag
#' Example: Biomass prior of 200 in 1985
#'  inp$priors$logB <- c(log(200), 0.2, 1985)
#'  inp$priors$logB <- c(log(200), 0.2, 1, 1985) # This includes the optional useflag
#' Example: Inverse gamma prior on sdb^2:
#'  inp$priors$isdb2gamma <- meanvar2shaperate(1/exp(inp$ini$logsdb)^2, 150^2)
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
#'  \item{"inp$stdevfacE"}{ Factors to multiply the observation error standard deviation of each individual effort observation. Can be used if some observations are more uncertain than others. A list with vectors of same length as observation vectors. Default: 1.}
#'  \item{"inp$mapsdi"}{ Vector of length equal to the number of index series specifying which indices that should use the same sdi. For example: in case of 3 index series use inp$mapsdi <- c(1, 1, 2) to have series 1 and 2 share sdi and have a separate sdi for series 3. Default: 1:nindex, where nindex is number of index series.}
#'  \item{"inp$seasontype"}{ If set to 1 use the spline-based representation of seasonality. If set to 2 use the oscillatory SDE system (this is more unstable and difficult to fit, but also more flexible).}
#' }
#' @param inp List of input variables, see details for required variables.
#' @return An updated list of input variables checked for consistency and with defaults added.
#' @examples
#' data(pol)
#' (inp <- check.inp(pol$albacore))
#' @export
check.inp <- function(inp){

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
    rm.neg <- function(tmpin, nam, nms, j='', sii=''){
        tmpout <- tmpin
        if (!is.null(tmpin[[nms[1]]])){
            nnms <- length(nms)
            neg <- which(tmpin[[nms[1]]] <= 0 | is.na(tmpin[[nms[1]]]))
            if (length(neg) > 0){
                for (i in 1:nnms){
                    tmpout[[nms[i]]] <- tmpin[[nms[i]]][-neg]
                }
                cat('Removing zero, negative, and NAs in ', nam, ' series ', j, ' stock ', sii, '\n')
            }
            return(tmpout)
        }
    }
    remove.neg2 <- function(inp, nam, sii, extrakeys=NULL){
        # sii is stock index
        nms <- c(paste0(c('obs', 'time', 'stdevfac'), nam), extrakeys)
        nnms <- length(nms)
        nna <- length(inp[[nms[1]]][[sii]]) # Number of index series
        if (nna > 0){
            for (j in 1:nna){
                tmp <- list()
                for (i in 1:nnms){
                    tmp[[nms[i]]] <- inp[[nms[i]]][[sii]][[j]]
                }
                tmpup <- rm.neg(tmp, nam, nms, j, sii)
                for (i in 1:nnms){
                    inp[[nms[i]]][[sii]][[j]] <- tmpup[[nms[i]]]
                }
            }
        }
        return(inp)
    }
    remove.neg <- function(inp, nam, extrakeys=NULL){
        nms <- c(paste0(c('obs', 'time', 'stdevfac'), nam), extrakeys)
        nnms <- length(nms)
        flag <- class(inp[[nms[1]]]) == 'list'
        if (flag){
            nna <- length(inp[[nms[1]]])
            if (nna > 0){
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
        if (class(invar) != 'list'){
            outvar <- list()
            outvar[[1]] <- invar
        } else {
            outvar <- invar
        }
        return(outvar)
    }
    make.time <- function(inp, nam, sii=NULL){
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
    make.time2 <- function(obs, nseasons=1){
        if (is.null(nseasons)){
            nseasons <- 1
        }
        nobs <- length(obs)
        to <- (nobs-1)/nseasons
        return(seq(0, to, length=nobs))
    }

    # -- DATA --
    if (!any(c('obsI', 'obsE') %in% names(inp))){
        stop('No effort or index observations. Please include index observations as a vector in inp$obsI and effort observations in inp$obsE.')
    }
    
    # CHECK INDEX OBSERVATIONS
    # One stock, one index included as a vector
    inp$obsI <- make.list(inp$obsI)
    inp$nstocks <- length(inp$obsI)
    if (!'stocknames' %in% names(inp)){
        inp$stocknames <- paste0('stock', 1:inp$nstocks)
    }
    if (length(inp$stocknames) != inp$nstocks){
        stop('Wrong length of inp$stocknames. Got ', length(inp$stocknames),
             ' but expected ', inp$nstocks)
    }
    # Stock as list, index as vector
    inp$nindex <- numeric(inp$nstocks)
    inp$nindexseq <- list()
    for (si in 1:inp$nstocks){
        inp$obsI[[si]] <- make.list(inp$obsI[[si]])
        inp$nindex[si] <- length(inp$obsI[[si]])
        if (inp$nindex[si] > 0){
            inp$nindexseq[[si]] <- 1:inp$nindex[si]
        }
    }
    # Stock as list, index as list
    # Time vector (create if doesn't exist)
    if (!'timeI' %in% names(inp)){
        inp$timeI <- vector('list', inp$nstocks)
    }
    if ('timeI' %in% names(inp)){
        inp$timeI <- make.list(inp$timeI)
        if (length(inp$timeI) != inp$nstocks){
            stop('Time vector specified for some index but not all. Needs specification for all or none (in which case assumptions about time are made).')
        }
        for (si in 1:inp$nstocks){
            inp$timeI[[si]] <- make.list(inp$timeI[[si]])
            #if (length(inp$timeI[[si]]) != inp$nindex[si]){
            #    stop('Lacking time vectors for stock ', si, ' got ', length(inp$timeI[[si]]),
            #         ' but need ', inp$nindex[si])
            #}
        }
    }
    # Standard deviation factor (create if doesn't exist)
    if (!'stdevfacI' %in% names(inp)){
        inp$stdevfacI <- vector('list', inp$nstocks)
    } else {
        inp$stdevfacI <- make.list(inp$stdevfacI)
    }
    inp$nobsI <- list()
    for (si in 1:inp$nstocks){
        # Time vector
        # Create if doesn't exist
        if (is.null(inp$timeI[[si]]) | length(inp$timeI[[si]]) == 0){
            inp$timeI[[si]] <- list()
            for (i in inp$nindexseq[[si]]){
                #inp$timeI[[si]][[i]] <- 1:length(inp$obsI[[si]][[i]])
                inp$timeI[[si]][[i]] <- make.time2(inp$obsI[[si]][[i]], nseasons=inp$nseasons)
            }
        } 
        inp$timeI[[si]] <- make.list(inp$timeI[[si]])
        # Standard deviation factor
        # Create if doesn't exist
        if (is.null(inp$stdevfacI[[si]])){
            inp$stdevfacI[[si]] <- list()
            for (i in inp$nindexseq[[si]]){
                inp$stdevfacI[[si]][[i]] <- rep(1, length(inp$obsI[[si]][[i]]))
            }
        } 
        inp$stdevfacI[[si]] <- make.list(inp$stdevfacI[[si]])
        if (inp$nindex[si] != length(inp$timeI[[si]])){
            stop('length(inp$timeI) is not equal to length(inp$obsI) for stock ', si)
        }
        if (inp$nindex[si] != length(inp$stdevfacI[[si]])){
            stop('length(inp$stdevfacI) is not equal to length(inp$obsI) for stock ', si)
        }
        for (i in inp$nindexseq[[si]]){
            base.checks(inp$obsI[[si]][[i]], inp$timeI[[si]][[i]],
                        inp$stdevfacI[[si]][[i]], paste0('I', i))
        }
        inp <- remove.neg2(inp, 'I', si)
        inp$nobsI[[si]] <- rep(0, inp$nindex[si]) # Need to be after negative have been removed
        for (i in inp$nindexseq[[si]]){
            inp$nobsI[[si]][i] <- length(inp$obsI[[si]][[i]])
        }
        if (length(inp$nobsI[[si]]) == 0){
            inp$nobsI[[si]] <- 0
        }
    }

    # CHECK EFFORT OBSERVATIONS
    # One fleet included as a vector
    inp$obsE <- make.list(inp$obsE)
    # neffort is number of observed effort series (some can contain no observations)
    inp$neffort <- length(inp$obsE)
    if (inp$neffort > 0){
        inp$neffortseq <- 1:inp$neffort
    } else {
        inp$neffortseq <- numeric(0)
    }
    if (!"effortobs2fleet" %in% names(inp)){
        inp$effortobs2fleet <- inp$neffortseq
    }
    if ('effortobs2fleet' %in% names(inp)){
        if (inp$neffort != length(inp$effortobs2fleet)){
            stop('Length of inp$effortobs2fleet does not match number of effort series input!')
        }
    }

    for (i in inp$neffortseq){
        if (is.null(inp$obsE[[i]])){
            inp$obsE[[i]] <- numeric(0)
        }
    }
    if (is.null(inp$timeE)){
        inp$timeE <- list()
    }
    for (i in inp$neffortseq){
        if (is.null(inp$timeE[[i]])){
            inp$timeE[[i]] <- make.time2(inp$obsE[[i]], nseasons=inp$nseasons)
        }
    }
    inp$timeE <- make.list(inp$timeE)
    if (length(inp$timeE) != inp$neffort){
        stop('Lacking effort time vectors, got ', length(inp$timeE),
             ' but need ', inp$neffort, '.')
    }
    # Effort intervals (dtef)
    if (!"dtef" %in% names(inp) & inp$neffort > 0){
        inp$dtef <- vector('list', inp$neffort)
    } else {
        if (length(inp$dtef) != inp$neffort){
            stop('Lacking dtef vectors, got ', length(inp$dtef),
                 ' but need ', inp$neffort, '.')
        }
    }
    if (!"stdevfacE" %in% names(inp) & inp$neffort > 0){
        inp$stdevfacE <- vector('list', inp$neffort)
    } else {
        if (length(inp$stdevfacE) != inp$neffort){
            stop('Lacking stdevfacE vectors, got ', length(inp$stdevfacE),
                 ' but need ', inp$neffort, '.')
        }
    }
    inp$nobsE <- numeric(inp$neffort)
    for (i in inp$neffortseq){
        if (any(diff(inp$timeE[[i]]) <= 0)){
            stop('Effort times for effort series', i, 'are not strictly increasing!')
        }
        if (is.null(inp$dtef[[i]])){
            if (!is.null(inp$timeE[[i]])){
                dtef <- diff(inp$timeE[[i]])
                if (length(dtef) > 0){
                    inp$dtef[[i]] <- min(dtef)
                } else {
                    inp$dtef[[i]] <- 1
                    if (length(inp$obsE[[i]]) != 0){
                        cat(paste('Effort interval (dtef) not specified and length of effort time series',
                                  i, 'shorter than 2. Assuming an interval of 1 year.\n'))
                    }
                }
            }
        }
        if (length(inp$dtef[[i]]) == 1){
            inp$dtef[[i]] <- rep(inp$dtef[[i]], length(inp$obsE[[i]]))
        }
        if (is.null(inp$stdevfacE[[i]])){
            inp$stdevfacE[[i]] <- rep(1, length(inp$obsE[[i]]))
        }
        base.checks(inp$obsE[[i]], inp$timeE[[i]], inp$stdevfacE[[i]], 'E')
    }
    inp <- remove.neg(inp, 'E', extrakeys='dtef')
    for (i in inp$neffortseq){
        inp$nobsE[i] <- length(inp$obsE[[i]])
        if (length(inp$dtef[[i]]) != inp$nobsE[i]){
            stop('Effort interval vector (inp$dtef, ', length(inp$dtef),
                 ') does not match effort observation vector (inp$obsE, ',
                 inp$nobsE, ') in length for effort series ', i)
        }
    }
    if (length(inp$nobsE) == 0){
        inp$nobsE <- 0
    }

    # CHECK CATCH OBSERVATIONS
    if (!'obsC' %in% names(inp)){
        stop('No catch observations included. Please include them as inp$obsC.')
    }
    inp$obsC <- make.list(inp$obsC)
    inp$nfisheries <- length(inp$obsC)
    if (inp$nfisheries > 0){
        inp$nfisheriesseq <- 1:inp$nfisheries
    }
    if (is.null(inp$timeC)){
        inp$timeC <- list()
        for (i in inp$nfisheriesseq){
            inp$timeC[[i]] <- make.time2(inp$obsC[[i]], nseasons=inp$nseasons)
        }
    } 
    inp$timeC <- make.list(inp$timeC)
    if (length(inp$timeC) != inp$nfisheries){
        stop('Lacking catch time vectors, got ', length(inp$timeC),
             ' but need ', inp$nfisheries, '.')
    }
    # Catch intervals (dtc)
    if (!"dtc" %in% names(inp)){
        inp$dtc <- vector('list', inp$nfisheries)
    } else {
        if (length(inp$dtc) != inp$nfisheries){
            stop('Lacking dtc vectors, got ', length(inp$dtc), ' but need ', inp$nfisheries, '.')
        }
    }
    if (!'stdevfacC' %in% names(inp)){
        inp$stdevfacC <- vector('list', inp$nfisheries)
    } else {
        if (length(inp$stdevfacC) != inp$nfisheries){
            stop('Lacking stdevfacC vectors, got ', length(inp$stdevfacC),
                 ' but need ', inp$nfisheries, '.')
        }
    }
    inp$nobsC <- numeric(inp$nfisheries)
    for (i in inp$nfisheriesseq){
        if (any(diff(inp$timeC[[i]]) <= 0)){
            stop('Catch times in series ', i, ' are not strictly increasing!')
        }
        # Catch intervals (dtc)
        if (is.null(inp$dtc[[i]])){
            dtc <- diff(inp$timeC[[i]])
            if (length(dtc) > 0){
                inp$dtc[[i]] <- min(dtc)
            } else {
                inp$dtc[[i]] <- 1
                cat(paste('Catch interval (dtc) not specified and length of catch time series', i, 'shorter than 2. Assuming an interval of 1 year.\n'))
            }
        }
        if (length(inp$dtc[[i]]) == 1){
            inp$dtc[[i]] <- rep(inp$dtc[[i]], length(inp$obsC[[i]]))
        }
        if (is.null(inp$stdevfacC[[i]])){
            inp$stdevfacC[[i]] <- rep(1, length(inp$obsC[[i]]))
        }
        base.checks(inp$obsC[[i]], inp$timeC[[i]], inp$stdevfacC[[i]], 'C')
    }
    inp <- remove.neg(inp, 'C', extrakeys='dtc')
    for (i in inp$nfisheriesseq){    
        inp$nobsC[i] <- length(inp$obsC[[i]])
        if (length(inp$dtc[[i]]) != inp$nobsC[i]){
            stop('Catch interval vector (inp$dtc) does not match catch observation vector (inp$obsC) in length in series ', i)
        }
    }

    # CHECK TARGET MATRIX
    # One row for each stock, one column for each fleet,
    # 0 means no stock-fleet interaction, integer > 0 indicates which catch series relate to that interaction
    if (!'target' %in% names(inp)){
        if ((inp$nstocks * inp$nfisheries) > 1){ # Multiple stocks and/or fleets
            stop('Data for running a multi-fleet/multi-stock model have been input without specifying inp$target. Cannot continue. See ?check.inp for help.')
        }
        inp$target <- matrix(1, 1, 1)
    }
    if (class(inp$target) == 'numeric'){
        if (length(inp$target) == 1){
            inp$target <- matrix(inp$target, 1, 1)
        }
    }
    if (nrow(inp$target) != inp$nstocks){
        stop('nrow(inp$target) (', nrow(inp$target), ') != inp$nstocks (', inp$nstocks,
             ') fix index input data or inp$target.')
    }
    #if (ncol(inp$target) != inp$nfleets){
    #    stop('ncol(inp$target) (', ncol(inp$target), ') != inp$nfleets (', inp$nfleets,
    #         ') fix effort input data or inp$target.')
    #}
    #if (any(is.na(match(as.numeric(inp$target), 0:1)))){ # Check whether only 0 and 1
    #    stop('Values other than 0 and 1 specified in inp$target.')
    #}
    if (all(inp$target == 0)){
        stop('All values of inp$target are 0, there should be at least one entry > 0.')
    }
    # nfleets is number of fleets fishing (some may be unobserved, i.e. no effort data)
    inp$nfleets <- ncol(inp$target)
    if (!'fleetnames' %in% names(inp)){
        inp$fleetnames <- paste0('fleet', 1:inp$nfleets)
    }
    if (length(inp$fleetnames) != inp$nfleets){
        stop('Wrong length of inp$fleetnames. Got ', length(inp$fleetnames),
             ' but expected ', inp$nfleets)
    }
    #inp$nfisheries <- sum(inp$target > 0)
    #if (inp$nfisheries > 0){
    #    inp$nfisheriesseq <- 1:inp$nfisheries
    #}
    # inp$target seems valid

    # Create target map
    ind2sub <- function(M, ind){
        return(c(row(M)[ind], col(M)[ind] ))
    }
    inp$targetmap <- matrix(0, inp$nfisheries, 3)
    rownames(inp$targetmap) <- paste0('F', inp$nfisheriesseq)
    colnames(inp$targetmap) <- c('stock', 'fleet', 'qf')
    nztarget <- as.numeric(inp$target[inp$target != 0])
    for (i in nztarget){
        j <- which(inp$target == i)
        inp$targetmap[i, 1:2] <- ind2sub(inp$target, j) # Get row (stock) and column (fleet)
    }
    # Find nqf, i.e. the number of observed fisheries
    inp$observedfisheries <- which(inp$targetmap[, 2] %in% inp$effortobs2fleet)
    inp$nqf <- length(inp$observedfisheries)
    if (inp$nqf > 0){
        inp$nqfseq <- 1:inp$nqf
    } else {
        inp$nqfseq <- numeric(0)
    }
    inds <- which(inp$targetmap[, 2] %in% inp$effortobs2fleet)
    inp$targetmap[inds, 3] <- inp$nqfseq
    if (nrow(inp$targetmap) != inp$nfisheries){
        stop('Number of catch time series (', inp$nfisheries,
             ') differs from the number of fleet stock interactions given in inp$target (',
             nrow(inp$targetmap), ').')
    }
    # Check whether all required indices are present in target matrix
    if (any(sort(inp$target[inp$target > 0]) != 1:length(inp$obsC))){
        stop('Mismatch between indices specified in inp$target and the number of catch series provided in inp$obsC!')
    }
    rownames(inp$target) <- inp$stocknames
    colnames(inp$target) <- inp$fleetnames

    # Find the fisheries that affect stock si
    inp$stock2f <- list()
    inp$stock2qf <- list()
    for (si in 1:inp$nstocks){
        inp$stock2f[[si]] <- which(inp$targetmap[, 1] == si)
        qfinds <- inp$targetmap[inp$targetmap[, 1] == si, 3]
        inp$stock2qf[[si]] <- qfinds[qfinds > 0]
    }
    
    # Vector containing all observation times
    timeobsall <- sort(c(unlist(inp$timeC), unlist(inp$timeC) + unlist(inp$dtc),
                         unlist(inp$timeI),
                         unlist(inp$timeE), unlist(inp$timeE) + unlist(inp$dtef)))
    dtimeobsall <- diff(timeobsall)
    if (any(dtimeobsall > 50)){
        stop('At least one gap over 50 years exists in data! cannot work with such data.')
    }
    if (any(dtimeobsall > 20)){
        warning('At least one gap over 20 years exists in data!')
    }
    
    # Observation IDs
    # Observation index in the observation input vector containing all observations
    # Initialise
    inp$obsidI <- vector('list', inp$nstocks)
    for (si in 1:inp$nstocks){
        inp$obsidI[[si]] <- vector('list', inp$nindex[si])
    }
    inp$obsidE <- vector('list', inp$neffort)
    inp$obsidC <- vector('list', inp$nfisheries)
    cur <- 0
    # Catch obs ID
    for (i in inp$nfisheriesseq){
        nextcur <- cur + inp$nobsC[i]
        inp$obsidC[[i]] <- (cur+1):(nextcur)
        cur <- nextcur
    }
    # Index obs ID
    for (si in 1:inp$nstocks){
        for (i in inp$nindexseq[[si]]){
            nextcur <- cur + inp$nobsI[[si]][i]
            inp$obsidI[[si]][[i]] <- (cur+1):(nextcur)
            cur <- nextcur
        }
    }
    # Effort obs ID
    for (i in inp$neffortseq){
        nextcur <- cur + inp$nobsE[i]
        if (inp$nobsE[i] > 0){
            ids <- (cur+1):(nextcur)
        } else {
            ids <- numeric(0)
        }
        inp$obsidE[[i]] <- ids
        cur <- nextcur
    }
    
    #inp$nseries <- 1 + inp$nindex + as.numeric(inp$nobsE > 0)
    inp$nseries <- inp$nfisheries + sum(inp$nindex) + inp$neffort

    # -- MODEL OPTIONS --
    if (!"RE" %in% names(inp)){
        inp$RE <- c('logE', 'logu', 'logB', 'logmre')
    }
    if (!"scriptname" %in% names(inp)) inp$scriptname <- 'spict'
    # Index related
    if (!"mapsdi" %in% names(inp)){
        inp$mapsdi <- list()
        sdiend <- 0
        for (si in 1:inp$nstocks){
            sdinext <- sdiend + inp$nindex[si]
            inp$mapsdi[[si]] <- (sdiend+1) : (sdinext)
            sdiend <- sdinext
        }
    }
    # TODO: Perhaps create a check of mapsdi here
    if ("mapsdi" %in% names(inp)){
        inp$nsdi <- numeric(inp$nstocks)
        for (si in 1:inp$nstocks){
            inp$nsdi[si] <- length(unique(inp$mapsdi[[si]]))
        }
    }
    tmp <- list()
    for (si in 1:inp$nstocks){
        tmp[[si]] <- rep(si, inp$nindex[si])
        #inp$index2sdb[inp$mapsdi[[si]]] <- si
    }
    if (!"mapq" %in% names(inp)){
        inp$mapq <- list()
        qend <- 0
        for (si in 1:inp$nstocks){
            qnext <- qend + inp$nindex[si]
            inp$mapq[[si]] <- (qend+1) : (qnext)
            qend <- qnext
        }
    }
    inp$index2sdb <- unlist(tmp)
    inp$index2sdi <- unlist(inp$mapsdi)
    inp$index2q <- unlist(inp$mapq)
    # TODO: Perhaps create a check of mapq here
    if ("mapq" %in% names(inp)){
        inp$nq <- numeric(inp$nstocks)
        for (si in 1:inp$nstocks){
            inp$nq[si] <- length(unique(inp$mapq[[si]]))
        }
    }
    # Effort related
    if (!"efforttype" %in% names(inp)){
        inp$efforttype <- 1
    }
    if (!"mapsde" %in% names(inp)){
        inp$mapsde <- inp$neffortseq
    }
    if ("mapsde" %in% names(inp)){
        inp$nsde <- length(unique(inp$mapsde))
    }
    if (!"mapqf" %in% names(inp)){
        inp$mapqf <- inp$nqfseq
    }
    #if ("mapqf" %in% names(inp)) inp$nqf <- length(unique(inp$mapqf))
    # Catch related
    if (!"catchunit" %in% names(inp)){
        inp$catchunit <- ''
    }
    if (!"mapsdf" %in% names(inp)){
        inp$mapsdf <- inp$nfisheriesseq
    }
    if ("mapsdf" %in% names(inp)){
        inp$nsdf <- length(unique(inp$mapsdf))
    }
    # Reporting
    if (!"reportall" %in% names(inp)){
        inp$reportall <- TRUE
    }
    if (!"do.sd.report" %in% names(inp)){
        inp$do.sd.report <- TRUE
    }
    if (!"bias.correct" %in% names(inp)){
        inp$bias.correct <- FALSE # This is time consuming
    }
    if (!"bias.correct.control" %in% names(inp)){
        inp$bias.correct.control <- list(sd=FALSE) # This is time consuming
    }
    if (!"getReportCovariance" %in% names(inp)){
        inp$getReportCovariance <- TRUE # Set to FALSE to save memory
    }
    # Simulation options
    #if (!"armalistF" %in% names(inp)) inp$armalistF <- list() # Used for simulating arma noise for F instead of white noise.
    # When simulating using sim.comm.cpue == TRUE, the simulation calculates
    # commercial CPUE by dividing catch and effort.
    # This is output as an biomass index.
    if (!"sim.comm.cpue" %in% names(inp)) inp$sim.comm.cpue <- FALSE
    # Optimiser options
    if (!"optimiser" %in% names(inp)) inp$optimiser <- 'nlminb'
    if (!"optimiser.control" %in% names(inp)) inp$optimiser.control <- list()
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
    if (!"robflagi" %in% names(inp)) inp$robflagi <- 0
    inp$robflagi <- as.numeric(inp$robflagi)
    if (!"robflage" %in% names(inp)) inp$robflage <- 0
    inp$robflage <- as.numeric(inp$robflage)
    # Time-varying growth option
    inp <- set.default(inp, 'timevaryinggrowth', inp$nqf > 0)
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

    # Options for simple model (this is probably broken when extending to multi fleet)
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
        inp$phases$logE <- -1
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

    # - Prediction horizons -
    # Catch prediction time step (dtpredc)
    if (!"dtpredc" %in% names(inp)){
        if (length(unlist(inp$dtc)) > 0){
            inp$dtpredc <- max(unlist(inp$dtc))
        } else {
            inp$dtpredc <- 1
            cat('Assuming a 1 year prediction interval for catch.\n')
        }
    }
    # Time point to predict catches until
    if (!"timepredc" %in% names(inp)){
        inp$timepredc <- max(timeobsall)
    } else {
        if (inp$timepredc < max(unlist(inp$timeC))){
            stop('inp$timepredc: ', inp$timepredc,
                ' must be equal to or later than last catch observation: ',
                max(inp$timeC), '!')
        }
    }
    # Time point to predict indices until
    if (!"timepredi" %in% names(inp)){
        inp$timepredi <- max(timeobsall)
    } else {
        if (sum(unlist(inp$nobsI)) > 0){
            if (inp$timepredi < max(unlist(inp$timeI))){
                stop('inp$timepredi must be equal to or later than last index observation!')
            }
        }
    }
    # Effort prediction time step (dtprede)
    if (!"dtprede" %in% names(inp)){
        if (sum(inp$nobsE) > 0){
            if (length(unlist(inp$dtef)) > 0){
                inp$dtprede <- max(unlist(inp$dtef))
            } else {
                inp$dtprede <- 1
                cat('Assuming a 1 year prediction interval for effort.\n')
            }
        } else {
            inp$dtprede <- numeric(0)
        }
    }
    # Time point to predict effort until
    if (!"timeprede" %in% names(inp)){
        if (sum(inp$nobsE) > 0){
            inp$timeprede <- max(timeobsall)
        } else {
            inp$timeprede <- numeric(0)
        }
    } else {
        if (sum(inp$nobsE) > 0){
            if (inp$timeprede < max(unlist(inp$timeE))){
                cat('inp$timeprede:', inp$timeprede,
                    ' must be equal to or later than last effort observation: ',
                    max(unlist(inp$timeE)), '!')
            }
        }
    }

    # This may give a problem if effort data has later time points than catches or index
    if (sum(inp$nobsE) > 0 & sum(unlist(inp$nobsI)) > 0){
        if (max(unlist(inp$timeE)) > max(unlist(inp$timeI), unlist(inp$timeC))){
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
            cat('The dteuler used is not allowed! using inp$dteuler:', inp$dteuler, '\n')            
        }
    }
    if (FALSE){ # The restriction on dteuler could possibly be removed
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
            cat('The dteuler used is not allowed! using inp$dteuler:', inp$dteuler, '\n')
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
    from <- min(timeobsall) # First time point of inp$time
    if (!inp$start.in.first.data.point){
        # Here we take floor of time of first data point
        # This can sometimes cause problems when estimating the random effect
        # prior to the first data point, because of no data support
        from <- floor(from)
    }
    # Include times of observations (including dtc)
    inp$time <- seq(from=from,
                    to=max(inp$timepredi,
                        inp$timepredc + inp$dtpredc,
                        inp$timeprede + inp$dtprede),
                    by=inp$dteuler)
    # inp$time now contains the "hard" Euler discretisation
    if (inp$eulertype == 'soft'){
        # Include observed time points in addition to Euler discretisation
        inp$time <- sort(unique(c(timeobsall, inp$time)))
    }
    if (!inp$eulertype %in% c('soft', 'hard')){
        stop('inp$eulertype must be either "soft" or "hard"!')
    }
    # Calculate time steps
    inp$dt <- c(diff(inp$time), inp$dteuler)
    inp$ns <- length(inp$time)
    # Find indices of inp$time relevant to each stock
    inp$indsstock <- list()
    for (si in 1:inp$nstocks){
        cinds <- which(inp$targetmap[, 1] == si)
        mintime <- min(c(unlist(inp$timeI[[si]]), unlist(inp$timeC[cinds])))
        inp$indsstock[[si]] <- which(inp$time >= mintime)
    }
    
    # -- DERIVED VARIABLES --
    inp$timerange <- range(timeobsall)
    inp$indlastobs <- cut(max(c(unlist(inp$timeC), unlist(inp$timeI), unlist(inp$timeE))),
                          inp$time, right=FALSE, labels=FALSE)
    inp$indest <- which(inp$time <= inp$timerange[2])
    inp$indpred <- which(inp$time >= inp$timerange[2])
    inp$indCpred <- list()
    for (i in inp$nfisheriesseq){
        inp$indCpred[[i]] <- which(inp$time >= max(inp$timeC[[i]] + inp$dtc[[i]]))
    }
    # Management
    inp$manstart <- ceiling(inp$time[inp$indpred[1]])
    if (!"ffac" %in% names(inp)) inp$ffac <- 1
    if ("ffac" %in% names(inp)){
        if (inp$ffac < 0){
            cat('Warning: ffac < 0, which is not allowed, setting ffac = 0.')
            inp$ffac <- 0
        }
    }
    if (!"fcon" %in% names(inp)) inp$fcon <- 0
    if ("fcon" %in% names(inp)){
        if (inp$fcon < 0){
            cat('Warning: fcon < 0, which is not allowed, setting fcon = 0.')
            inp$fcon <- 0
        }
    }
    if (!"ffacvec" %in% names(inp)){
        inp$ffacvec <- numeric(inp$ns) + 1
        # -1 in indpred because 1 is for plotting
        #inp$ffacvec[inp$indpred[-1]] <- inp$ffac + 1e-8 # Add small to avoid taking log of 0
        # Start in indpred[2] because indpred[1] is mainly for plotting
        inp$ffacvec[inp$indpred[2]] <- inp$ffac + 1e-8 # Add small to avoid taking log of 0
    }
    if (!"fconvec" %in% names(inp)){
        inp$fconvec <- numeric(inp$ns)
        # -1 in indpred because 1 is for plotting
        inp$fconvec[inp$indpred[-1]] <- inp$fcon + 1e-8 # Add small to avoid taking log of 0
    }
    inp$ffaceuler <- inp$ffac^inp$dteuler
    # Seasons
    if (!"nseasons" %in% names(inp)){
        expnseasons <- 1/min(unlist(inp$dtc))
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
    inp$seasonindex <- 1 / inp$dteuler * (inp$time %% 1)
    inp$seasons <- rep(0, inp$ns)
    for (i in 1:inp$nseasons){
        frac <- 1 / inp$nseasons
        modtime <- inp$time %% 1
        inds <- which(modtime >= ((i-1)*frac) & modtime < (i*frac))
        inp$seasons[inds] <- i
    }
    # ic is the indices of inp$time to which catch observations correspond
    inp$timeCpred <- list()
    inp$nobsCp <- numeric(inp$nfisheries)
    inp$idCpred <- list() # To which fishery does a catch prediction belong
    inp$dtcp <- list()
    inp$ic <- list()
    inp$icpred <- list() # Indices of timeCpred corresponding to observations
    idend <- 0
    for (i in inp$nfisheriesseq){
        if (length(inp$dtc[[i]]) > 0){
            dtcpred <- min(inp$dtc[[i]])
        } else {
            dtcpred <- 1
        }
        # Find fisheries related to the stock of the current fishery (i)
        ffinds <- which(inp$targetmap[, 1] == inp$targetmap[i, 1])
        # Predict catches for each fishery over the span of all fisheries related to a stock
        # to be able to plot predicted catch even when data are missing
        timeCstock <- unlist(inp$timeC[ffinds])
        timeC <- inp$timeC[[i]]
        inp$timeCpred[[i]] <- unique(c(timeCstock,
                                       seq(tail(timeCstock, 1), inp$timepredc, by=dtcpred)))
        inp$nobsCp[i] <- length(inp$timeCpred[[i]])
        inp$dtcp[[i]] <- numeric(inp$nobsCp[i])
        indsthis <- match(timeC, timeCstock) # This fishery
        indsother <- setdiff(1:inp$nobsCp[i], indsthis) # Other fisheries on this stock
        inp$dtcp[[i]][indsthis] <- inp$dtc[[i]]
        inp$dtcp[[i]][indsother] <- dtcpred
        #inp$dtcp[[i]] <- c(inp$dtc[[i]], rep(dtcpred, inp$nobsCp[i]-inp$nobsC[i]))
        #inp$dtcp[[i]] <- dtcpred
        inp$idCpred[[i]] <- idend + (1:inp$nobsCp[i])
        inp$icpred[[i]] <- inp$idCpred[[i]][indsthis]
        idend <- tail(inp$idCpred[[i]], 1)
        inp$ic[[i]] <- cut(inp$timeCpred[[i]], inp$time, right=FALSE, labels=FALSE)
    }
    # nc is number of states to integrate a catch observation over
    inp$nc <- list()
    inp$iff <- list()
    for (i in inp$nfisheriesseq){
        inp$nc[[i]] <- rep(0, inp$nobsCp[i])
        inp$iff[[i]] <- rep(i, inp$nobsCp[i])
        for (j in 1:inp$nobsCp[i]){
            inp$nc[[i]][j] <- sum(inp$time >= inp$timeCpred[[i]][j]
                                  & inp$time < (inp$timeCpred[[i]][j] + inp$dtcp[[i]][j]))
        }
        if (any(inp$nc[[i]] == 0)){
            stop('Current inp$dteuler is too large to accommodate some catch intervals. Make inp$dteuler smaller!')
        }
    }

    # ie is the indices of inp$time to which effort observations correspond
    if (length(unlist(inp$dtef)) > 0){
        dtefpred <- min(unlist(inp$dtef))
    } else {
        dtefpred <- 1
    }
    inp$timeEpred <- list()
    inp$nobsEp <- numeric(inp$neffort)
    inp$dtefp <- list()
    inp$ie <- list()
    inp$ifleet <- list()
    for (i in inp$neffortseq){
        if (inp$nobsE[i] > 0){
            inp$timeEpred[[i]] <- unique(c(inp$timeE[[i]],
                                           (seq(tail(inp$timeE[[i]], 1),
                                                inp$timeprede, by=dtefpred))))
        } else {
            inp$timeEpred[[i]] <- numeric(0)
        }
        inp$nobsEp[i] <- length(inp$timeEpred[[i]])
        inp$dtefp[[i]] <- c(inp$dtef[[i]], rep(dtefpred, inp$nobsEp[i]-inp$nobsE[i]))
        inp$ie[[i]] <- cut(inp$timeEpred[[i]], inp$time, right=FALSE, labels=FALSE)
        inp$ifleet[[i]] <- rep(inp$effortobs2fleet[i], inp$nobsE[i])
    }
    # ne is number of states to integrate an effort observation over
    inp$ne <- list()
    for (i in inp$neffortseq){
        inp$ne[[i]] <- rep(0, inp$nobsEp[i])
        if (inp$nobsE[i] > 0){
            for (j in 1:inp$nobsEp[i]){
                inp$ne[[i]][j] <- sum(inp$time >= inp$timeEpred[[i]][j]
                                 & inp$time < (inp$timeEpred[[i]][j] + inp$dtefp[[i]][j]))
            }
        }
        if (any(inp$ne[[i]] == 0)){
            stop('Current inp$dteuler is too large to accommodate some effort intervals. Make inp$dteuler smaller!')
        }
    }
    
    # ii is the indices of inp$time to which index observations correspond
    inp$ii <- list()
    inp$ib <- list()
    for (si in 1:inp$nstocks){
        inp$ii[[si]] <- list()
        inp$ib[[si]] <- list()
        for (i in inp$nindexseq[[si]]){
            inp$ii[[si]][[i]] <- cut(inp$timeI[[si]][[i]], inp$time, right=FALSE, labels=FALSE)
            inp$ib[[si]][[i]] <- rep(si, inp$nobsI[[si]][i])
        }
    }

    # Translate observations from a list to a vector
    list2vec <- function(inlist){
        vec <- unlist(inlist)
        if (is.null(vec)){
            vec <- numeric(0)
        }
        return(vec)
    }
    replist <- function(map, nobs){
        out <- rep(map, times=nobs)
        if (is.null(out)){
            out <- numeric(0)
        }
        return(out)
    }
    # Catch related
    inp$icin <- list2vec(inp$ic)
    inp$ncin <- list2vec(inp$nc)
    inp$icpredin <- list2vec(inp$icpred)
    inp$iffin <- list2vec(inp$iff)
    inp$stdevfacCin <- list2vec(inp$stdevfacC)
    inp$isdcin <- rep(inp$nfisheriesseq, inp$nobsC)
    # Index related
    inp$obsIin <- list2vec(inp$obsI)
    inp$stdevfacIin <- list2vec(inp$stdevfacI)
    inp$iiin <- list2vec(inp$ii)
    inp$ibin <- list2vec(inp$ib)
    inp$iqin <- numeric(0)
    for (si in 1:inp$nstocks){
        inp$iqin <- c(inp$iqin, replist(inp$mapq[[si]], inp$nobsI[[si]]))
    }
    inp$isdiin <- numeric(0)
    for (si in 1:inp$nstocks){
        inp$isdiin <- c(inp$isdiin, replist(inp$mapsdi[[si]], inp$nobsI[[si]]))
    }
    # Effort related
    inp$iein <- list2vec(inp$ie)
    inp$nein <- list2vec(inp$ne)
    inp$ifleetin <- list2vec(inp$ifleet)
    inp$stdevfacEin <- list2vec(inp$stdevfacE)
    inp$isdein <- replist(inp$neffortseq, inp$nobsE)
    
    # Add helper variable such that predicted catch can be calculated using small euler steps
    # Need to include timerange[2] and exclude timerange[2]+dtpred because the catch at t is acummulated over t to t+dtc.
    # dtpredcinds is the indices of time of the catch prediction
    inp$dtpredcinds <- which(inp$time >= inp$timepredc & inp$time < (inp$timepredc + inp$dtpredc))
    inp$dtpredcnsteps <- length(inp$dtpredcinds)
    inp$dtprediind <- cut(inp$timepredi, inp$time, right=FALSE, labels=FALSE)
    inp$dtpredeinds <- which(inp$time >= inp$timeprede & inp$time < (inp$timeprede + inp$dtprede))
    inp$dtpredensteps <- length(inp$dtpredeinds)

    # - Sort observations in time and store in one vector -
    timeobsseen <- c(unlist(inp$timeC) + unlist(inp$dtc) - 1e-4,
                     unlist(inp$timeI),
                     unlist(inp$timeE) + unlist(inp$dtef) - 1e-5)
    # Add dtc to timeC because the total catch is first "seen" at the end of the given catch interval (typically year or quarter), similarly for quarter. By arbitrary convention catches are assumed to be "seen" before effort although they should be observed at the same time.
    srt <- sort(timeobsseen, index=TRUE)
    timeobs <- c(unlist(inp$timeC), unlist(inp$timeI), unlist(inp$timeE))
    timeobssrt <- timeobs[srt$ix]
    obs <- log(c(unlist(inp$obsC), unlist(inp$obsI), unlist(inp$obsE)))
    obsid <- c(unlist(inp$obsidC), unlist(inp$obsidI), unlist(inp$obsidE))
    # TODO: Could add a check here that obsid is equal to 1:length(obs)
    inp$obssrt <- obs[srt$ix]
    inp$timeobssrt <- timeobs[srt$ix]
    inp$obsidsrt <- obsid[srt$ix]
    # Find indices of catch, index and effort in sorted observation vector 
    #inp$isc <- match(1:inp$nobsC, srt$ix)
    inp$isc <- match(unlist(inp$obsidC), srt$ix)
    if (sum(unlist(inp$nobsI)) > 0){
        #inp$isi <- match((inp$nobsC+1):(inp$nobsC+sum(inp$nobsI)), srt$ix)
        inp$isi <- match(unlist(inp$obsidI), srt$ix)
    } else {
        inp$isi <- numeric(0)
    }
    if (sum(inp$nobsE) > 0){
        #inp$ise <- match((inp$nobsC+sum(inp$nobsI)+1):(inp$nobsC+sum(inp$nobsI)+sum(inp$nobsE)), srt$ix)
        inp$ise <- match(unlist(inp$obsidE), srt$ix)
    } else {
        inp$ise <- numeric(0)
    }
    if (sum(inp$nobsC) != length(inp$isc)){
        warning('Mismatch between length(inp$isc) ', length(inp$isc),
                ' and sum(inp$nobsC) ', sum(inp$nobsC), '.')
    }
    if (sum(unlist(inp$nobsI)) != length(inp$isi)){
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

    # -- PRIORS --
    # Priors are assumed Gaussian and specified in a vector of length 3:
    # c(log(mean), stdev in log, useflag).
    # log(mean): log of the mean of the prior distribution.
    # stdev in log: standard deviation of the prior distribution in log domain.
    # useflag: if 1 then the prior is used, if 0 it is not used. Default is 0.
    check.prior <- function(priors, priorname){
        priorvec <- priors[[priorname]]
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
                    warning('Year for prior on ', priorname, ' (', priorvec[3],
                            ') did not match times where this RE is estimated. Not using this prior. To fix this use a year where an observation is available.')
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
                        'iqgamma', 'logqf', 'logbkfrac', 'logB', 'logF', 'logBBmsy',
                        'logFFmsy', 'logsdb', 'isdb2gamma', 'logsdf', 'isdf2gamma',
                        'logsdi', 'isdi2gamma', 'logsde', 'isde2gamma', 'logsdc', 'isdc2gamma')
    repriors <- c('logB', 'logF', 'logBBmsy', 'logFFmsy')
    npossiblepriors <- length(possiblepriors)
    if (!"priors" %in% names(inp)){
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
    if (!'logn' %in% names(inp$priors)){
        inp$priors$logn <- c(logn, lognsd)
    }
    if (!'logalpha' %in% names(inp$priors)){
        inp$priors$logalpha <- c(logalpha, logalphasd)
    }
    if (!'logbeta' %in% names(inp$priors)){
        inp$priors$logbeta <- c(logbeta, logbetasd)
    }

    # Remaining priors, set to something, but will not be used
    if ("priors" %in% names(inp)){
        # Remove wrong priors names
        nms <- names(inp$priors)
        inds <- which(is.na(match(nms, possiblepriors)))
        if (length(inds) > 0){
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
                    inp$priors[[nm]] <- c(log(1e-4), 0.2, 0) # FE priors
                }
            } else {
                inp$priors[[nm]] <- check.prior(inp$priors, nm)
            }
        }
    }
    npriors <- length(inp$priors)
    inp$priorsuseflags <- numeric(npriors)
    for (i in 1:npriors){
        inp$priorsuseflags[i] <- inp$priors[[i]][3]
    }
    
    # -- MODEL PARAMETERS --
    check.ini <- function(inp, parnm, expectedlength){
        lngt <- length(inp$ini[[parnm]])
        if (lngt != expectedlength){
            stop('Specified initial value for ', parnm, ' has length ', lngt, ' but ',
                 expectedlength, ' is expected!')
        }
    }
    rm.zeros <- function(vec){
        return(vec[vec > 0])
    }
    # logn
    if (!"logn" %in% names(inp$ini)){
        inp$ini$logn <- rep(logn, inp$nstocks)
    } else {
        if (any(inp$ini$logn == 0)){
            stop('Initial value(s) for logn == 0. That is not valid!')
        }
    }
    check.ini(inp, 'logn', inp$nstocks)
    # Calculate gamma from n
    n <- exp(inp$ini$logn)
    inp$ini$gamma <- calc.gamma(n)
   
    # logK
    if (!'logK' %in% names(inp$ini)){
        inp$ini$logK <- numeric(inp$nstocks)
        for (si in 1:inp$nstocks){
            inp$ini$logK[si] <- log(4 * max(unlist(inp$obsC[rm.zeros(inp$target[si, ])])))
        }
    }
    check.ini(inp, 'logK', inp$nstocks)
    # logr
    if ('logr' %in% names(inp$ini) & 'logm' %in% names(inp$ini)){
        inp$ini$logr <- NULL # If both r and m are specified use m and discard r
    }
    if (!'logr' %in% names(inp$ini)){
        if (!'logm' %in% names(inp$ini)){
            inp$ini$logm <- numeric(inp$nstocks)
            for (si in 1:inp$nstocks){
                obsC <- unlist(inp$obsC[rm.zeros(inp$target[si, ])[1]])
                obsI <- inp$obsI[[si]][[1]]
                inp$ini$logm[si] <- unname(log(guess.m2(obsC, obsI)))
            }
        }
        n <- exp(inp$ini$logn)
        r <- exp(inp$ini$logm) / exp(inp$ini$logK) * n^(n/(n-1))
        inp$ini$logr <- log(r)
    }
    # logm
    if (!"logm" %in% names(inp$ini)){
        n <- exp(inp$ini$logn)
        m <- exp(inp$ini$logr) * exp(inp$ini$logK) / (n^(n/(n-1)))
        inp$ini$logm <- log(m)
    }
    logm <- inp$ini$logm # Store this to be able to set logmre later
    if (inp$timevaryinggrowth){
        inp$ini$logm <- log(1)
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

    check.mapped.ini.qsdi <- function(inp, nam, map){
        if (nam %in% names(inp$ini)){
            if (length(inp$ini[[nam]]) != sum(inp$nindex)){ # nq is given by mapq
                if (length(inp$ini[[nam]]) == 1){
                    val <- inp$ini[[nam]]
                    for (si in 1:inp$nstocks){
                        for (i in inp$nindexseq[[si]]){
                            j <- map[[si]][i]
                            inp$ini[[nam]][j] <- val
                        }
                    }
                } else {
                    stop('The length of ', nam, ' in inp$ini (', length(inp$ini[[nam]]),
                         ') does not fit with the number of parameters to be estimated (',
                         sum(inp$nindex), ').')
                }
            }
        }
        return(inp)
    }

    # logq
    if (!'logq' %in% names(inp$ini)){
        inp$ini$logq <- numeric(sum(inp$nindex))
        for (si in 1:inp$nstocks){
            for (i in inp$nindexseq[[si]]){
                j <- inp$mapq[[si]][i]
                logmaxI <- log(max(inp$obsI[[si]][[i]]))
                inp$ini$logq[j] <- logmaxI - inp$ini$logK[si]
            }
        }
    }
    inp <- check.mapped.ini.qsdi(inp, 'logq', inp$mapq)
    
    #if ('logq' %in% names(inp$ini)){
    #    if (length(inp$ini$logq) != sum(inp$nindex)){
    #        if (length(inp$ini$logq) == 1){
    #            logq <- inp$ini$logq
    #            for (si in 1:inp$nstocks){
    #                for (i in inp$nindexseq[[si]]){
    #                    j <- inp$mapq[[si]][i]
    #                    inp$ini$logq[j] <- logq
    #                }
    #            }
    #        } else {
    #            stop('The length of logq in inp$ini (', length(inp$ini$logq),
    #                 ') does not fit with the number of parameters to be estimated (',
    #                 sum(inp$nindex), ').')
    #        }
    #    }
    #}
    
    # logsdi
    if (!'logsdi' %in% names(inp$ini)){
        inp$ini$logsdi <- log(0.2)
    }
    if (sum(unlist(inp$nobsI)) > 0){
        inp <- check.mapped.ini.qsdi(inp, 'logsdi', inp$mapsdi)
    }

    # logsdb
    if (!'logsdb' %in% names(inp$ini)){
        inp$ini$logsdb <- log(0.2)
    }
    inp <- check.mapped.ini(inp, 'logsdb', 'nstocks')
    # logsdm
    inp$ini <- set.default(inp$ini, 'logsdm', log(0.2))
    inp <- check.mapped.ini(inp, 'logsdm', 'nstocks')
    # logqf
    if (!'logqf' %in% names(inp$ini)){
        inp$ini$logqf <- numeric(inp$nqf)
        for (i in inp$nqfseq){
            #j <- which(inp$target == i)
            #rc <- ind2sub(inp$target, j) # Get row (stock) and column (fleet)
            si <- inp$targetmap[inp$observedfisheries[i], 1]
            logmaxE <- log(max(inp$obsE[[i]]))
            inp$ini$logqf[i] <- inp$ini$logr[si] - logmaxE
        }
    }
    if (sum(inp$nobsE) > 0){
        inp <- check.mapped.ini(inp, 'logqf', 'nqf')
    }
    # logsdf (process stdev)
    if (!'logsdf' %in% names(inp$ini)){
        inp$ini$logsdf <- log(0.2)
    }
    inp <- check.mapped.ini(inp, 'logsdf', 'nfisheries')
    # logsdu
    if (!'logsdu' %in% names(inp$ini)){
        inp$ini$logsdu <- log(0.1)
    }
    inp <- check.mapped.ini(inp, 'logsdu', 'nfisheries')
    # logsdc
    if (!'logsdc' %in% names(inp$ini)){
        inp$ini$logsdc <- log(0.2)
    }
    inp <- check.mapped.ini(inp, 'logsdc', 'nfisheries')
    
    # NOTE: sde and sdf may need rethinking depending on whether random effects will be E or F
    # RE: E will be RE but sdf will still be process noise and sde will be observation noise.
    if (!'logsde' %in% names(inp$ini)){
        inp$ini$logsde <- log(0.2)
    }
    if (sum(inp$nobsE) > 0){
        inp <- check.mapped.ini(inp, 'logsde', 'neffort')
    }
    
    #if (sum(inp$nobsE)>0) inp <- check.mapped.ini(inp, 'logsde', 'nsde')
    #if (!"logalpha" %in% names(inp$ini)) inp$ini$logalpha <- logalpha
    #if (sum(inp$nobsI)>0) inp <- check.mapped.ini(inp, 'logalpha', 'nsdi')
    #if (!"logbeta" %in% names(inp$ini))  inp$ini$logbeta <- logbeta

    # Fill in unspecified (more rarely user defined) model parameter values
    if (!"loglambda" %in% names(inp$ini)){
        inp$ini$loglambda <- rep(log(0.1), inp$nfisheries)
    }
    if ("logphi" %in% names(inp$ini)){
        if (length(inp$ini$logphi)+1 != dim(inp$splinemat)[2]){
            cat('Mismatch between length of ini$logphi and number of columns of splinemat! removing prespecified ini$logphi and setting default.\n')
            inp$ini$logphi <- NULL
        }
    }
    if (!"logphi" %in% names(inp$ini)){
        inp$ini$logphi <- rep(0, inp$nseasons-1)
    }
    if (!"logitpp" %in% names(inp$ini)){
        inp$ini$logitpp <- log(0.95/(1-0.95))
    }
    if (!"logp1robfac" %in% names(inp$ini)){
        inp$ini$logp1robfac <- log(15-1)
    }
    if (!"logbkfrac" %in% names(inp$ini)){
        inp$ini$logbkfrac <- log(0.8)
    }
    check.mat <- function(mat, dims, nam){
        if (any(is.na(mat))){
            stop('NAs in ', nam, ', cannot continue!')
        }
        if (is.null(dim(mat))){
            mat <- matrix(mat, 1, 1)
        }
        if (dim(mat)[1] != dims[1] | dim(mat)[2] != dims[2]){
            warning('Wrong dimension of ', nam, ': ', dim(mat)[1], ' x ',
                    dim(mat)[2], ' should be equal to ', dims[1], ' x ', dims[2],
                    '. Redefining to correct size.')
            mat <- matrix(mat[1, 1], dims[1], dims[2])
        }
        return(mat)
    }
    if (!"logE" %in% names(inp$ini)){
        inp$ini$logE <- matrix(log(0.2) + inp$ini$logr[1], inp$nfleets, inp$ns)
    }
    if ("logE" %in% names(inp$ini)){    
        inp$ini$logE <- check.mat(inp$ini$logE, c(inp$nfleets, inp$ns), 'inp$ini$logE')
    }
    if (!"logu" %in% names(inp$ini)){
        inp$ini$logu <- matrix(log(1) + 1e-3, 2*length(inp$ini$logsdu), inp$ns)
    }
    if ("logF" %in% names(inp$ini)){    
        inp$ini$logu <- check.mat(inp$ini$logu, c(2*length(inp$ini$logsdu), inp$ns), 'inp$ini$logu')
    }
    if (!"logB" %in% names(inp$ini)){
        inp$ini$logB <- matrix(inp$ini$logK + log(0.5), inp$nstocks, inp$ns)
    }
    if ("logB" %in% names(inp$ini)){    
        inp$ini$logB <- check.mat(inp$ini$logB, c(inp$nstocks, inp$ns), 'inp$ini$logB')
    }
    if (!"logmre" %in% names(inp$ini)){
        if (inp$timevaryinggrowth){
            inp$ini$logmre <- matrix(logm, inp$nstocks, inp$ns)
        } else {
            inp$ini$logmre <- matrix(log(1), inp$nstocks, inp$ns)
        }
    }
    if ("logmre" %in% names(inp$ini)){    
        inp$ini$logmre <- check.mat(inp$ini$logmre, c(inp$nstocks, inp$ns), 'inp$ini$logmre')
    }
    
    # ids to extract logB, logF etc
    #inp$idstock <- unlist(lapply(1:inp$nstocks, function(x) rep(x, inp$ns)))
    #inp$idfleet <- unlist(lapply(1:inp$nfleets, function(x) rep(x, inp$ns)))
    #inp$idqf <- unlist(lapply(1:inp$nqf, function(x) rep(x, inp$ns)))
    inp$idstock <- rep(1:inp$nstocks, inp$ns)
    inp$idfleet <- rep(1:inp$nfleets, inp$ns)
    inp$idfishery <- rep(inp$nfisheriesseq, inp$ns)
    
    # Reorder parameter list
    inp$parlist <- list(logm = inp$ini$logm,
                        logK = inp$ini$logK,
                        logq = inp$ini$logq,
                        logqf = inp$ini$logqf,
                        logn = inp$ini$logn,
                        logsdb = inp$ini$logsdb,
                        logsdu = inp$ini$logsdu,
                        logsdf = inp$ini$logsdf,
                        logsdi = inp$ini$logsdi,
                        logsde = inp$ini$logsde,
                        logsdc = inp$ini$logsdc,
                        logsdm = inp$ini$logsdm,
                        logphi = inp$ini$logphi,
                        loglambda = inp$ini$loglambda,
                        logitpp = inp$ini$logitpp,
                        logp1robfac = inp$ini$logp1robfac,
                        logE = inp$ini$logE,
                        logu = inp$ini$logu,
                        logB = inp$ini$logB,
                        logmre = inp$ini$logmre)
    
    # Determine fixed parameters
    forcefixpars <- c() # Parameters that are forced to be fixed.
    if (inp$nseasons == 1){
        forcefixpars <- c('logphi', 'logu', 'logsdu', 'loglambda', forcefixpars)
    } else {
        if (inp$seasontype == 1){ # Use spline
            forcefixpars <- c('logu', 'logsdu', 'loglambda', forcefixpars)
        }
        if (inp$seasontype == 2){ # Use coupled SDEs
            forcefixpars <- c('logphi', forcefixpars)
        }
    }
    if (inp$robflagc == 0 & inp$robflagi == 0 & inp$robflage == 0){
        forcefixpars <- c('logitpp', 'logp1robfac', forcefixpars)
    }
    if (sum(inp$nobsE) == 0){
        forcefixpars <- c('logqf', 'logsde', forcefixpars)
    }
    if (sum(unlist(inp$nobsI)) == 0){
        forcefixpars <- c('logq', 'logsdi', forcefixpars)
    }
    if (inp$timevaryinggrowth){
        forcefixpars <- c('logm', forcefixpars)
    } else {
        forcefixpars <- c('logmre', 'logsdm', forcefixpars)
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
    for (nm in names(inp$parlist)){
        if (!nm %in% names(inp$phases)){
            inp$phases[[nm]] <- 1
        }
    }
    
    nphasepars <- length(inp$phases)
    phasevec <- unlist(inp$phases)
    phases <- unique(unname(phasevec))
    phases <- phases[phases > 0] # Don't include phase -1 (fixed value)
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
