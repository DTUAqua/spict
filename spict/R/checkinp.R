#' @name check.inp
#' @title Check list of input variables
#' @details Fills in defalut values if missing.
#' Required inputs:
#' \itemize{
#'  \item{"inp$obsC"}{ Vector of catch observations.}
#'  \item{"inp$obsI"}{ List containing vectors of index observations.}
#'  \item{"inp$ini$logr"}{ Initial value(s) for logr (log intrinsic growth rate). Can be specified as a scalar, or a vector of length 2 or 4. If length(logr)=2 different r values are estimated for first and second half of the year, if length(logr)=4 different r values are estimated for the four quarters of the year.}
#'  \item{"inp$ini$logK"}{ Initial value for logK (log carrying capacity).}
#'  \item{"inp$ini$logq"}{ Initial value for logq (log catchability of index).}
#'  \item{"inp$ini$logsdb"}{ Initial value for logsdb (log standard deviation of log biomass).}
#'  \item{"inp$ini$logsdf"}{ Initial value for logsdf (log standard deviation of log fishing mortality).}
#' }
#' Optional inputs:
#' - Data
#' \itemize{
#'  \item{"inp$timeC"}{ Vector of catch times. Default: even time steps starting at 1.}
#'  \item{"inp$timeI"}{ List containing vectors of index times. Default: even time steps starting at 1.}
#'  \item{"inp$dtc"}{ Time interval for catches, e.g. for annual catches inp$dtc=1, for quarterly catches inp$dtc=0.25. Can be given as a scalar, which is then used for all catch observations. Can also be given as a vector specifying the catch interval of each catch observation. Default: min(diff(inp$timeC)). }
#'  \item{"inp$nseasons"}{ Number of within-year seasons in data. If inp$nseasons > 1 then a cyclic B spline is used to impose a seasonal pattern in F. The parameters of the spline are phi and are estimated. Valid values are 1, 2 or 4. Default: number of unique within-year time points present in data.}
#' }
#' - Parameters
#' \itemize{
#'  \item{"inp$ini$logn"}{ Pella-Tomlinson exponent. Default: log(2) corresponding to the Schaefer formulation.}
#'  \item{"inp$ini$phi"}{ Vector for cyclic B spline representing within-year seasonal variation. Default: rep(1, inp$nseasons).}
#'  \item{"inp$ini$logalpha"}{ sdi = alpha*sdb. Default: log(1).}
#'  \item{"inp$ini$logbeta"}{ sdc = beta*sdf. Default: log(1).}
#'  \item{"inp$ini$logF"}{ Default: log(0.2*r).}
#'  \item{"inp$ini$logB"}{ Default: log(0.5*K).}
#' }
#' - Priors
#' Priors on model parameters are assumed Gaussian and specified in a vector of length 3: c(log(mean), stdev in log, useflag).
#' log(mean): log of the mean of the prior distribution.
#' stdev in log: standard deviation of the prior distribution in log domain.
#' useflag: if 1 then the prior is used, if 0 it is not used. Default is 0.
#' Example: inp$priors$logr <- c(log(0.8), 0.1, 1)
#' - Settings
#' \itemize{
#'  \item{"inp$dtpredc"}{ Length of catch prediction interval in years. Default: max(inp$dtc).}
#'  \item{"inp$dtpredi"}{ Length of index prediction interval in years. Default: 0.}
#'  \item{"inp$do.sd.report"}{ Flag indicating whether SD report (uncertainty of derived quantities) should be calculated. For small values of inp$dteuler this may require a lot of memory. Default: TRUE.}
#'  \item{"inp$reportall"}{ Flag indicating whether quantities derived from state vectors (e.g. B/Bmsy, F/Fmsy etc.) should be calculated by SD report. For small values of inp$dteuler (< 1/32) reporting all may have to be set to FALSE for sdreport to run. Additionally, if only reference points of parameter estimates are of interest one can set to FALSE to gain a speed-up. Default: TRUE.}
#' \item{"inp$robflagc"}{ Flag indicating whether robust estimation should be used for catches (either 0 or 1). Default: 0.}
#'  \item{"inp$robflagi"}{ Flag indicating whether robust estimation should be used for indices (either 0 or 1). Default: 0.}
#'  \item{"inp$ffac"}{ Management scenario represented by a factor to multiply F with when calculating the F of the next time step. ffac=0.8 means a 20\% reduction in F over the next year. The factor is only used when predicting beyond the data set. Default: 1 (0\% reduction).}
#'  \item{"inp$dteuler"}{ Length of Euler time step in years. Default: min(inp$dtc).}
#'  \item{"inp$phases"}{ Phases can be used to fix/free parameters and estimate in different stages or phases. To fix e.g. logr at inp$ini$logr set inp$phases$logr <- -1. To free logalpha and estimate in phase 1 set inp$phases$logalpha <- 1.}
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
        if(length(inp$obsC) != length(inp$timeC)) stop('Time and observation vector do not match in length for catch series')
        neg <- which(inp$obsC<0 | is.na(inp$obsC))
        if(length(neg)>0){
            inp$obsC <- inp$obsC[-neg]
            inp$timeC <- inp$timeC[-neg]
            cat(paste('Removing negative and NAs in catch series\n'))
        }
        inp$nobsC <- length(inp$obsC)
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
        inp$nindex <- length(inp$obsI)
        if(inp$nindex != length(inp$timeI)){
            stop('The length(inp$timeI) does not match length(inp$obsI)!')
        }
        inp$nobsI <- rep(0, inp$nindex)
        for(i in 1:inp$nindex){
            if(length(inp$obsI[[i]]) != length(inp$timeI[[i]])) stop('Time and observation vector do not match in length for index series ',i)
            if(length(inp$obsI[[i]])>0){
                neg <- which(inp$obsI[[i]]<0 | is.na(inp$obsI[[i]]))
                if(length(neg)>0){
                    inp$obsI[[i]] <- inp$obsI[[i]][-neg]
                    inp$timeI[[i]] <- inp$timeI[[i]][-neg]
                    cat(paste('Removing negative and NAs in index series',i,'\n'))
                }
            }
            inp$nobsI[i] <- length(inp$obsI[[i]])
        }
    } else {
        stop('No index observations included. Please include them as a list in inp$obsI.')
    }

    # -- MAX MIN RATIO --
    inp$maxminratio <- rep(0, inp$nindex)
    names(inp$maxminratio) <- paste0('I', 1:inp$nindex)
    for(i in 1:inp$nindex){
        if(length(inp$obsI[[i]])>0) inp$maxminratio[i] <- max(inp$obsI[[i]])/min(inp$obsI[[i]])
    }
    # -- MODEL OPTIONS --
    if(!"msytype" %in% names(inp)){
        inp$msytype <- 's'
    } else {
        if(!inp$msytype %in% c('s', 'd')) stop('inp$msytype must be either "s" (stochastic) or "d" (deterministic!')
    }
    if(!"reportall" %in% names(inp)) inp$reportall <- TRUE
    if(!"do.sd.report" %in% names(inp)) inp$do.sd.report <- TRUE
    if(!"robflagc" %in% names(inp)) inp$robflagc <- 0
    inp$robflagc <- as.numeric(inp$robflagc)
    if(!"robflagi" %in% names(inp)) inp$robflagi <- 0
    inp$robflagi <- as.numeric(inp$robflagi)
    if(!"ffac" %in% names(inp)) inp$ffac <- 1
    if(!"lamperti" %in% names(inp)) inp$lamperti <- 1
    inp$lamperti <- 1 # Since Pella-Tomlinson form was implemented only Lamperti is allowed
    if(!"euler" %in% names(inp)) inp$euler <- 1
    inp$euler <- 1 # Since Pella-Tomlinson form was implemented only Euler is allowed
    if(!"dtc" %in% names(inp)){
        dtc <- diff(inp$timeC)
        if(length(dtc)>0){
            inp$dtc <- min(dtc)
            #cat(paste('Catch interval (dtc) not specified. Assuming an interval of:', inp$dtc, 'year.\n'))
        } else {
            inp$dtc <- 1
            cat(paste('Catch interval (dtc) not specified and length of catch time series shorter than 1. Assuming an interval of 1 year.\n'))
        }
    }
    if(length(inp$dtc)==1) inp$dtc <- rep(inp$dtc, inp$nobsC)
    if(length(inp$dtc) != inp$nobsC) stop('Catch interval vector (inp$dtc) does not match catch observation vector (inp$obsC) in length')
    if(!"dteuler" %in% names(inp)){
        if(length(inp$dtc)>0){
            inp$dteuler <- min(inp$dtc)
        } else {
            inp$dteuler <- diff(inp$timeI[[1]])
        }
    }
    inp$ffaceuler <- inp$ffac^inp$dteuler
    timeobs <- c(inp$timeC, inp$timeC + inp$dtc, unlist(inp$timeI))
    if(!"timepredc" %in% names(inp)){
        inp$timepredc <- max(timeobs)
    } else {
        if(inp$timepredc < max(timeobs)) stop('inp$timepredc must be later than last observation!')
    }
    if(!"dtpredc" %in% names(inp)){
        if("dtpred" %in% names(inp)){
            inp$dtpredc <- inp$dtpred
        } else {
            if(length(inp$dtc)>0){
                inp$dtpredc <- max(inp$dtc)
            } else {
                inp$dtpredc <- 1
                cat('Assuming a 1 year prediction interval for catch.\n')
            }
        }
    }
    if(!"dtpredi" %in% names(inp)){
        if("dtpred" %in% names(inp)){
            inp$dtpredi <- inp$dtpred
        } else {
            #if(length(inp$dtc)>0){
            #    inp$dtpredi <- max(inp$dtc)
            #} else {
                inp$dtpredi <- 0
                #cat('Assuming a 0 year prediction interval for index.\n')
            #}
        }
    }
    dtpredmax <- max(c(inp$dtpredc, inp$dtpredi))
    if(!"RE" %in% names(inp)) inp$RE <- c('logF', 'logB')
    #if(!"RE" %in% names(inp)) inp$RE <- c('logF', 'logB', 'logbkfrac')
    #if(!"RE" %in% names(inp)) inp$RE <- c('logF', 'logB', 'Cpredcum')
    if(!"scriptname" %in% names(inp)) inp$scriptname <- 'spict'

    # -- ASPIC SETTINGS --
    if(!"aspic" %in% names(inp)) inp$aspic <- list()
    if(!"mode" %in% names(inp$aspic)) inp$aspic$mode <- 'FIT'
    if(!"verbosity" %in% names(inp$aspic)) inp$aspic$verbosity <- '102'
    if(!"nboot" %in% names(inp$aspic)) inp$aspic$nboot <- 1000
    if(!"ciperc" %in% names(inp$aspic)) inp$aspic$ciperc <- 95

    # -- SIMULATION SETTINGS --
    if(!"armalistF" %in% names(inp)) inp$armalistF <- list() # Used for simulating arma noise for F instead of white noise.

    # -- DERIVED VARIABLES --
    timeobsnodtc <- c(inp$timeC, unlist(inp$timeI))
    # Include dtc because a catch observation at time t includes information in the interval [t; t+dtc[ 
    inp$timerange <- range(timeobs)
    #time <- seq(min(timeobs), max(timeobs)+dtpredmax, by=inp$dteuler)
    time <- seq(min(timeobs), max(max(timeobs)+inp$dtpredi, inp$timepredc+inp$dtpredc), by=inp$dteuler)
    # Remove duplicate time points and store time in inp list
    inp$time <- sort(unique(c(timeobs, time)))
    inp$ns <- length(inp$time)
    inp$indlastobs <- which(inp$time == max(c(inp$timeC, unlist(inp$timeI))))
    inp$indest <- which(inp$time <= inp$timerange[2])
    inp$indpred <- which(inp$time >= inp$timerange[2])
    # Seasons
    if(!"nseasons" %in% names(inp)){
        inp$nseasons <- length(unique(timeobs %% 1))
        #if(inp$seasons>4) inp$seasons <- 4
    }
    if("nseasons" %in% names(inp)){
       if(!inp$nseasons %in% c(1, 2, 4)) stop('inp$nseasons (=', inp$nseasons, ') must be either 1, 2 or 4.')
    }
    # Calculate seasonal spline
    if("splineorder" %in% names(inp)){
        if(inp$nseasons<4 & inp$splineorder>2) inp$splineorder <- 2
    } else {
        inp$splineorder <- ifelse(inp$nseasons<4, 2, 3)
    }
    inp$splinemat <- make.splinemat(inp$nseasons, inp$splineorder, dtfine=inp$dteuler)
    #if(inp$dteuler < 1){
        #knots <- seq(0, 1, by=1/inp$nseasons)
        #x <- seq(0, 1-inp$dteuler, by=inp$dteuler)
        #inp$splinemat <- cSplineDes(x, knots, ord=inp$splineorder)
    #} else {
     #   inp$splinemat <- matrix(1, 1, 1)
    #}
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
    inp$ic <- match(inp$timeCpred, inp$time)
    # nc is number of states to integrate a catch observation over
    inp$nc <- rep(0, inp$nobsCp)
    for(i in 1:inp$nobsCp) inp$nc[i] <- sum(inp$time >= inp$timeCpred[i] & inp$time < (inp$timeCpred[i]+inp$dtcp[i]))
    # ii is the indices of inp$time to which index observations correspond
    inp$ii <- list()
    for(i in 1:inp$nindex) inp$ii[[i]] <- match(inp$timeI[[i]], inp$time)
    # Translate index observations from a list to a vector
    inp$obsIin <- unlist(inp$obsI)
    inp$iiin <- unlist(inp$ii)
    inp$iqin <- rep(1:inp$nindex, times=inp$nobsI)
    # Calculate time steps
    inp$dt <- diff(inp$time)
    # Add helper variable such that predicted catch can be calculated using small euler steps
    # Need to include timerange[2] and exclude timerange[2]+dtpred because the catch at t is acummulated over t to t+dtc.
    #inp$dtpredcinds <- which(inp$time >= inp$timerange[2] & inp$time < (inp$timerange[2]+inp$dtpredc))
    inp$dtpredcinds <- which(inp$time >= inp$timepredc & inp$time < (inp$timepredc+inp$dtpredc))
    inp$dtpredcnsteps <- length(inp$dtpredcinds)
    inp$dtprediind <- which(inp$time == (inp$timerange[2]+inp$dtpredi))

    # -- PRIORS --
    # Priors are assumed Gaussian and specified in a vector of length 3: c(log(mean), stdev in log, useflag).
    # log(mean): log of the mean of the prior distribution.
    # stdev in log: standard deviation of the prior distribution in log domain.
    # useflag: if 1 then the prior is used, if 0 it is not used. Default is 0.
    if(!"priors" %in% names(inp)){
        inp$priors <- list()
    }
    if(!"logr" %in% names(inp$priors)){
        inp$priors$logr <- c(log(0.8), 0.1, 0)
    } else {
        if(inp$priors$logr[3] == 1){ # If prior is to be used
            if(inp$priors$logr[2] <= 0) cat(paste('WARNING: invalid standard deviation specified in prior for r (must be > 0). Not using this prior.\n'))
        }
    }
    
    # -- MODEL PARAMETERS --
    # logn
    if(!"logn" %in% names(inp$ini)){
        inp$ini$logn <- log(2.0)
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
    #print(inp$ini$logr)
    if('logr' %in% names(inp$ini)){
        nr <- length(inp$ini$logr)
        #if(!nr %in% c(1, 2, 4)){
        #    stop('logr must be a vector of length either 1, 2 or 4.')
        #}
        #min <- log(0.3)-3
        #max <- log(0.3)+3
        #for(i in 1:nr) if(inp$ini$logr[i] < min | inp$ini$logr[i] > max) stop('Please specify a value for logr(', i, ') in the interval: [', min, ';', max, ']')
        # ir is the mapping from time to the r vector
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
        if(length(inp$ini$logq) != inp$nindex){
            if(length(inp$ini$logq) == 1){
                inp$ini$logq <- rep(inp$ini$logq, inp$nindex)
            } else {
                stop('The length of inp$ini$logq (', length(inp$ini$logq), ') does not fit with the number of index series (', inp$nindex, ')')
            }
        }
    }
    if(!'logsdb' %in% names(inp$ini)) inp$ini$logsdb <- log(0.2)
    if(!'logsdf' %in% names(inp$ini)) inp$ini$logsdf <- log(0.2)
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
    if(!"logphi" %in% names(inp$ini)) inp$ini$logphi <- rep(0, inp$nseasons-1)
    if(!"logitpp" %in% names(inp$ini)) inp$ini$logitpp <- log(0.95/(1-0.95))
    if(!"logp1robfac" %in% names(inp$ini)) inp$ini$logp1robfac <- log(20-1)
    if(!"logalpha" %in% names(inp$ini)) inp$ini$logalpha <- log(1)
    if(!"logbeta" %in% names(inp$ini)) inp$ini$logbeta <- log(1)
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
    if(!"logB" %in% names(inp$ini)){
        #inp$ini$logB <- log(rep(exp(inp$ini$logbkfrac)*exp(inp$ini$logK), inp$ns))
        inp$ini$logB <- rep(inp$ini$logK + log(0.5), inp$ns)
    } else {
        if(length(inp$ini$logB) != inp$ns){
            cat('Wrong length of inp$ini$logB:', length(inp$ini$logB), ' Should be equal to inp$ns:', inp$ns, ' Setting length of logF equal to inp$ns (removing beyond inp$ns).\n')
            inp$ini$logB <- inp$ini$logB[1:inp$ns]
        }
    }
    #if(!"Cpredcum" %in% names(inp$ini)){
    #    inp$ini$Cpredcum <- rep(sum(inp$obsC)/inp$ns, inp$ns-1)
    #}

    # Reorder parameter list
    inp$parlist <- list(logphi=inp$ini$logphi,
                        logitpp=inp$ini$logitpp,
                        logp1robfac=inp$ini$logp1robfac,
                        logalpha=inp$ini$logalpha,
                        logbeta=inp$ini$logbeta,
                        #logbkfrac=inp$ini$logbkfrac,
                        #logF0=inp$ini$logF0,
                        logm=inp$ini$logm,
                        logK=inp$ini$logK,
                        logq=inp$ini$logq,
                        #logqf=inp$ini$logqf,
                        logn=inp$ini$logn,
                        logsdf=inp$ini$logsdf,
                        logsdb=inp$ini$logsdb,
                        logF=inp$ini$logF,
                        logB=inp$ini$logB,
                        #Cpredcum=inp$ini$Cpredcum,
                        dum=0.0)
    # Determine phases and fixed parameters
    if(inp$nseasons==1){
        fixpars <- c('logphi', 'logalpha', 'logbeta', 'logn', 'logitpp', 'logp1robfac', 'dum') # These are fixed unless specified
    } else {
        fixpars <- c('logalpha', 'logbeta', 'logn', 'logitpp', 'logp1robfac', 'dum') # These are fixed unless otherwise specified
    }
    if(!"phases" %in% names(inp)){
        inp$phases <- list()
        for(i in 1:length(fixpars)) inp$phases[[fixpars[i]]] <- -1
    } else {
        if("logr" %in% names(inp$phases)){
            inp$phases$logm <- inp$phases$logr
            inp$phases$logr <- NULL
        }
        for(i in 1:length(fixpars)) if(!fixpars[i] %in% names(inp$phases)) inp$phases[[fixpars[i]]] <- -1
    }
    # If robust flags are set to 1 then set phases for robust parameters to 1
    if(inp$robflagc==1 | inp$robflagi==1){
        inp$phases$logitpp <- 1
        inp$phases$logp1robfac <- 1
    }
    # Assign phase 1 to parameters without a phase
    nms <- names(inp$parlist)
    nnms <- length(nms)
    for(i in 1:nnms){
        if(!nms[i] %in% names(inp$phases)){
            inp$phases[[nms[i]]] <- 1
        }
    }
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
                cat(paste('WARNING: Phase specified for an invalid parameter:', parnam, '\n'))
            }
        }
    }
    if(!is.null(inp)) class(inp) <- "spictcls"
    return(inp)
}
