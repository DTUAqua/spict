#' @exportClass spictcls
setClass("spictcls")

#' @name test.spict
#' @title Example of a spict analysis.
#' @details Loads a data set, fits the model, calculates one-step-ahead residuals, plots the results.
#' @return A result report as given by fit.spict().
#' @examples
#' rep <- test.spict()
#' @export
test.spict <- function(dataset='albacore'){
    # Load data
    data(pol)
    inp <- pol[[dataset]]
    if(dataset=='albacore'){
        inp$dtpred <- 6
        inp$ffac <- 0.9
    }
    # Fit model
    rep <- fit.spict(inp)
    # Calculate one-step-ahead residuals
    rep <- calc.osa.resid(rep)
    # Plot results
    graphics.off()
    dev.new(width=10, height=10)
    plot(rep)
    summary(rep)
    return(rep)
}


#' @name pol
#' @title Fisheries data for south Atlantic albacore tuna 1967-1989.
#' @details This data were included in Polacheck et al. (1993).
#' @docType data
#' @keywords datasets
#' @usage data(pol)
#' @source Polacheck et al. (1993), Canadian Journal of Fisheries and Aquatic Science, vol 50, pp. 2597-2607.
#' @examples
#' data(pol)
#' rep <- fit.spict(inp=pol$albacore)
#' @format pol$albacore is a list containing the data and initial values for estimation formatted to be used as an input to fit.spict().
NULL


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
#' \itemize{
#'  \item{"inp$timeC"}{ Vector of catch times. Default: even time steps starting at 1.}
#'  \item{"inp$timeI"}{ List containing vectors of index times. Default: even time steps starting at 1.}
#'  \item{"inp$ini$logp"}{ Pella-Tomlinson exponent. Default: log(1) corresponding to the Schaefer formulation.}
#'  \item{"inp$ini$logbkfrac"}{ Default: log(0.3).}
#'  \item{"inp$ini$logF0"}{ Default: log(0.2*r).}
#'  \item{"inp$ini$phi1"}{ Default: 1.}
#'  \item{"inp$ini$phi2"}{ Default: 0.}
#'  \item{"inp$ini$logalpha"}{ Default: 0.}
#'  \item{"inp$ini$logbeta"}{ Default: 0.}
#'  \item{"inp$ini$logF"}{ Default: logF0.}
#'  \item{"inp$ini$logB"}{ Default: bkfrac*K.}
#'  \item{"inp$ini$robflagc"}{ Flag indicating whether robust estimation should be used for catches (either 0 or 1). Default: 0.}
#'  \item{"inp$ini$robflagi"}{ Flag indicating whether robust estimation should be used for indices (either 0 or 1). Default: 0.}
#'  \item{"inp$ffac"}{ Management scenario represented by a factor to multiply F with when calculating the F of the next time step. ffac=0.8 means a 20\% reduction in F over the next year. The factor is only used when predicting beyond the data set. Default: 1 (0\% reduction).}
#'  \item{"inp$lamperti"}{ Logical indicating whether to use Lamperti transformed equations (recommended). Default: TRUE.}
#'  \item{"inp$euler"}{ Logical indicating whether to use Euler time discretisation (recommended). Default: TRUE.}
#'  \item{"inp$dtc"}{ Time interval for catches, e.g. for annual catches inp$dtc=1, for quarterly catches inp$dtc=0.25. Can be given as a scalar, which is then used for all catch observations. Can also be given as a vector specifying the catch interval of each catch observation. Default: min(diff(inp$timeC)). }
#'  \item{"inp$dteuler"}{ Length of Euler time step in years. Default: min(inp$dtc).}
#'  \item{"inp$delay"}{ Include a delayed effect in the the F-process: F_i = phi1*F_i-1 + phi2*F_i-delay. Note the delay is counted in the number of Euler time steps (dteuler). So if dteuler is refined delay also needs to change. Default: 1.}
#'  \item{"inp$phases"}{ Phases can be used to fix/free parameters and estimate in different stages or phases. To fix e.g. logr at inp$ini$logr set inp$phases$logr <- -1. To free logalpha and estimate in phase 1 set inp$phases$logalpha <- 1.}
#' }
#' @param inp List of input variables, see details for required variables.
#' @return An updated list of input variables checked for consistency and with defaults added.
#' @examples
#' data(pol)
#' (inp <- check.inp(pol$albacore))
#' @export
check.inp <- function(inp){
    # -- DATA --
    # Check catches
    if('obsC' %in% names(inp)){
        if(!'timeC' %in% names(inp)) inp$timeC <- 1:length(inp$obsC)
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
            inp$timeI[[1]] <- 1:length(inp$obsI[[1]])
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
            neg <- which(inp$obsI[[i]]<0 | is.na(inp$obsI[[i]]))
            if(length(neg)>0){
                inp$obsI[[i]] <- inp$obsI[[i]][-neg]
                inp$timeI[[i]] <- inp$timeI[[i]][-neg]
                cat(paste('Removing negative and NAs in index series',i,'\n'))
            }
            inp$nobsI[i] <- length(inp$obsI[[i]])
        }
    } else {
        stop('No index observations included. Please include them as a list in inp$obsI.')
    }

    # -- MODEL OPTIONS --
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
    if(!"delay" %in% names(inp)) inp$delay <- 1
    if(!"dtpredc" %in% names(inp)){
        if("dtpred" %in% names(inp)){
            inp$dtpredc <- inp$dtpred
        } else {
            if(length(inp$dtc)>0){
                inp$dtpredc <- min(inp$dtc)
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
            if(length(inp$dtc)>0){
                inp$dtpredi <- min(inp$dtc)
            } else {
                inp$dtpredi <- 1
                cat('Assuming a 1 year prediction interval for index.\n')
            }
        }
    }
    dtpredmax <- max(c(inp$dtpredc, inp$dtpredi))
    if(!"RE" %in% names(inp)) inp$RE <- c('logF', 'logB')
    if(!"scriptname" %in% names(inp)) inp$scriptname <- 'spict'

    # -- ASPIC SETTINGS --
    if(!"aspic" %in% names(inp)) inp$aspic <- list()
    if(!"mode" %in% names(inp$aspic)) inp$aspic$mode <- 'FIT'
    if(!"verbosity" %in% names(inp$aspic)) inp$aspic$verbosity <- '102'
    if(!"nboot" %in% names(inp$aspic)) inp$aspic$nboot <- 1000
    if(!"ciperc" %in% names(inp$aspic)) inp$aspic$ciperc <- 95

    # -- DERIVED VARIABLES --
    timeobs <- inp$timeC
    timeobs <- c(timeobs, inp$timeC + inp$dtc)
    for(i in 1:inp$nindex) timeobs <- c(timeobs, inp$timeI[[i]])
    inp$timerange <- range(timeobs)
    time <- seq(inp$timerange[1], inp$timerange[2]+dtpredmax, by=inp$dteuler)
    # Remove duplicate time points and store time in inp list
    inp$time <- sort(unique(c(timeobs, time)))
    inp$ns <- length(inp$time)
    inp$indlastobs <- which(inp$time == inp$timerange[2])
    inp$indest <- which(inp$time <= inp$timerange[2])
    inp$indpred <- which(inp$time >= inp$timerange[2])
    # ic is the indices of inp$time to which catch observations correspond
    if(length(inp$dtc)>0){
        dtcpred <- min(inp$dtc)
    } else {
        dtcpred <- 1
    }
    inp$timeCp <- unique(c(inp$timeC, (tail(inp$timeC,1) + seq(0, inp$dtpredc, by=dtcpred))))
    inp$nobsCp <- length(inp$timeCp)
    inp$dtcp <- c(inp$dtc, rep(dtcpred, inp$nobsCp-inp$nobsC))
    inp$ic <- match(inp$timeCp, inp$time)
    # nc is number of states to integrate a catch observation over
    inp$nc <- rep(0, inp$nobsCp)
    for(i in 1:inp$nobsCp) inp$nc[i] <- sum(inp$time >= inp$timeCp[i] & inp$time < (inp$timeCp[i]+inp$dtcp[i]))
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
    # Need to include timerange[2] and exclude timerange[2]+dtpred because the catch at t is acummulated over t to t+dtc.
    inp$dtpredcinds <- which(inp$time >= inp$timerange[2] & inp$time < (inp$timerange[2]+inp$dtpredc))
    inp$dtpredcnsteps <- length(inp$dtpredcinds)
    inp$dtprediind <- which(inp$time == (inp$timerange[2]+inp$dtpredi))
    
    # -- MODEL PARAMETERS --
    # Check that required initial values are specified
    if(!'logr' %in% names(inp$ini)){
        stop('Please specify initial value(s) for logr!')
    } else {
        nr <- length(inp$ini$logr)
        if(!nr %in% c(1, 2, 4)){
            stop('logr must be a vector of length either 1, 2 or 4.')
        }
        #min <- log(0.3)-3
        #max <- log(0.3)+3
        #for(i in 1:nr) if(inp$ini$logr[i] < min | inp$ini$logr[i] > max) stop('Please specify a value for logr(', i, ') in the interval: [', min, ';', max, ']')
        # ir is the mapping from time to the r vector
        inp$ir <- rep(0, inp$ns)
        for(i in 1:nr){
            frac <- 1/nr
            modtime <- inp$time %% 1
            inds <- which(modtime>=((i-1)*frac) & modtime<(i*frac))
            inp$ir[inds] <- i
        }
    }
    check.ini('logK', inp)
    check.ini('logq', inp)
    if(!'logq' %in% names(inp$ini)){
        stop('Please specify initial value(s) for logq!')
    } else {
        if(length(inp$ini$logq) != inp$nindex){
            if(length(inp$ini$logq) == 1){
                inp$ini$logq <- rep(inp$ini$logq, inp$nindex)
            } else {
                stop('The length of inp$ini$logq (', length(inp$ini$logq), ') does not fit with the number of index series (', inp$nindex, ')')
            }
        }
    }
    check.ini('logsdb', inp, min=log(0.03), max=log(5))
    check.ini('logsdf', inp, min=log(0.03), max=log(5))
    # Fill in unspecified less important model parameter values
    if(!"phi1" %in% names(inp$ini)) inp$ini$phi1 <- 1
    if(!"phi2" %in% names(inp$ini)) inp$ini$phi2 <- 0
    if(!"logitpp" %in% names(inp$ini)) inp$ini$logitpp <- log(0.95/(1-0.95))
    if(!"logp1robfac" %in% names(inp$ini)) inp$ini$logp1robfac <- log(20-1)
    if(!"logalpha" %in% names(inp$ini)) inp$ini$logalpha <- log(1)
    if(!"logbeta" %in% names(inp$ini)) inp$ini$logbeta <- log(1)
    if(!"logbkfrac" %in% names(inp$ini)) inp$ini$logbkfrac <- log(0.3)
    if(!"logp" %in% names(inp$ini)) inp$ini$logp <- log(1.0)
    if(!"logF0" %in% names(inp$ini)) inp$ini$logF0 <- log(0.2*exp(inp$ini$logr[1]))
    if(!"logF" %in% names(inp$ini)){
        inp$ini$logF <- rep(inp$ini$logF0, inp$ns)
    } else {
        if(length(inp$ini$logF) != inp$ns){
            cat('Wrong length of inp$ini$logF: ', length(inp$ini$logF), ' Should be equal to inp$ns: ', inp$ns, '. Setting length of logF equal to inp$ns (removing beyond inp$ns).\n')
            inp$ini$logF <- inp$ini$logF[1:inp$ns]
        }
    }
    if(!"logB" %in% names(inp$ini)){
        inp$ini$logB <- log(rep(exp(inp$ini$logbkfrac)*exp(inp$ini$logK), inp$ns))
    } else {
        if(length(inp$ini$logB) != inp$ns){
            cat('Wrong length of inp$ini$logB: ', length(inp$ini$logB), ' Should be equal to inp$ns: ', inp$ns, '. Setting length of logF equal to inp$ns (removing beyond inp$ns).\n')
            inp$ini$logB <- inp$ini$logB[1:inp$ns]
        }
    }

    # Reorder parameter list
    inp$ini <- list(phi1=inp$ini$phi1,
                    phi2=inp$ini$phi2,
                    logitpp=inp$ini$logitpp,
                    logp1robfac=inp$ini$logp1robfac,
                    logalpha=inp$ini$logalpha,
                    logbeta=inp$ini$logbeta,
                    #loggamma=inp$ini$loggamma,
                    logbkfrac=inp$ini$logbkfrac,
                    logF0=inp$ini$logF0,
                    logr=inp$ini$logr,
                    logK=inp$ini$logK,
                    logq=inp$ini$logq,
                    logp=inp$ini$logp,
                    logsdf=inp$ini$logsdf,
                    logsdb=inp$ini$logsdb,
                    logF=inp$ini$logF,
                    logB=inp$ini$logB)
    # Determine phases and fixed parameters
    if(inp$delay==1){
        fixpars <- c('phi1', 'phi2', 'logalpha', 'logbeta', 'logp', 'logitpp', 'logp1robfac') # These are fixed unless specified
    } else {
        fixpars <- c('logalpha', 'logbeta', 'logp', 'logitpp', 'logp1robfac') # These are fixed unless otherwise specified
    }
    if(!"phases" %in% names(inp)){
        inp$phases <- list()
        for(i in 1:length(fixpars)) inp$phases[[fixpars[i]]] <- -1
    } else {
        for(i in 1:length(fixpars)) if(!fixpars[i] %in% names(inp$phases)) inp$phases[[fixpars[i]]] <- -1
    }
    # If robust flags are set to 1 then set phases for robust parameters to 1
    if(inp$robflagc==1 | inp$robflagi==1){
        inp$phases$logitpp <- 1
        inp$phases$logp1robfac <- 1
    }
    # Assign phase 1 to parameters without a phase
    nms <- names(inp$ini)
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
            if(parnam %in% names(inp$ini)){
                phase <- inp$phases[[parnam]]
                inp$map[[j]][[parnam]] <- factor(rep(NA, length(inp$ini[[parnam]])))
                #cat('WARNING: Phases not yet implemented! will estimate everything in one phase.\n')
            } else {
                cat(paste('WARNING: Phase specified for an invalid parameter:', parnam, '\n'))
            }
        }
    }
        
    #inp$checked <- TRUE
    return(inp)
}


#' @name check.ini
#' @title Check initial values of required model parameters.
#' @details Checks that specified initial values are valid (i.e. within a certain required range). If they are not the program stops and the user is notified with an error.
#' @param parname Parameter name to be checked.
#' @param inp List of input variables as output by check.inp.
#' @param min Minimum allowed value of parameter.
#' @param max Maximum allowed value of parameter.
#' @return If the check passes nothing is returned, otherwise an error is returned.
check.ini <- function(parname, inp, min=NULL, max=NULL){
    if(!parname %in% names(inp$ini)) stop('Please specify an initial value for ', parname, '!')
    #if(!is.null(min) & !is.null(max)){
    #    if(inp$ini[[parname]] < min | inp$ini[[parname]] > max) stop('Please specify a value for ', parname, ' in the interval: [', min, ';', max, ']')
    #}
}


#' @name get.predmap
#' @title Get the map used to calculate one-step-ahead predictions.
#' @details When calculating osa predictions using TMB spict is run by sequentially including one data point at a time while keeping the model parameters fixed at their ML estimates. This is obtained using a map, which contains the names of the parameters that need to be fixed.
#' @param guess The guess used when running the TMB estimation.
#' @param RE The names of the random effects of the model.
#' @return A map that can be input to TMB to fix all model parameters except the random effects.
get.predmap <- function(guess, RE){
    FE <- setdiff(names(guess), RE)
    predmap <- rep(factor(NA), length(FE))
    names(predmap) <- FE
    predmap <- as.list(predmap)
    return(predmap)
}


#' @name fit.spict
#' @title Fit a continuous-time surplus production model to data.
#' @details Fits the model using the TMB package and returns a result report containing estimates of model parameters, random effects (biomass and fishing mortality), reference points (Fmsy, Bmsy, MSY) including uncertainties given as standard deviations.
#'
#' Fixed effects:
#' \itemize{
#'   \item{"logbkfrac"}{ Log of B0/K.}
#'   \item{"logF0"}{ Log of initial value (F0) of the fishing mortality process.}
#'   \item{"logr"}{ Log of intrinsic growth rate.}
#'   \item{"logK"}{ Log of carrying capacity.}
#'   \item{"logq"}{ Log of catchability vector.}
#'   \item{"logsdf"}{ Log of standard deviation of fishing mortality process noise.}
#'   \item{"logsdb"}{ Log of standard deviation of biomass process noise.}
#' }
#'
#' Optional fixed effects (which are normally not estimated):
#' \itemize{
#'   \item{"phi1"}{ Used in the delayed F-process: F_i = phi1*F_i-1 + phi2*F_i-delay. (normally set to phi1=1)}
#'   \item{"phi2"}{ Used in the delayed F-process: F_i = phi1*F_i-1 + phi2*F_i-delay. (normally set to phi2=1)}
#'   \item{"logalpha"}{ Proportionality factor for the observation noise of the indices and the biomass process noise: sdi = exp(logalpha)*sdb. (normally set to logalpha=0)}
#'   \item{"logbeta"}{ Proportionality factor for the observation noise of the catches and the fishing mortality process noise: sdc = exp(logbeta)*sdf. (this is often difficult to estimate and can result in divergence of the optimisation. Normally set to logbeta=0)}
#' }
#'
#' Random effects:
#' \itemize{
#'   \item{"logB"}{ Log of the biomass process given by the continuous-time stochastic Schaefer equation: dB_t = r*B_t*(1-B_t/K)*dt + sdb*dW_t, where dW_t is Brownian motion.}
#'   \item{"logF"}{ Log of the fishing mortality process given by: dF_t = sdf*dV, where dV is Brownian motion.}
#' }
#'
#' Derived parameters:
#' \itemize{
#'   \item{"logBmsy"}{ Log of the equilibrium biomass (Bmsy) when fished at Fmsy.}
#'   \item{"logFmsy"}{ Log of the fishing mortality (Fmsy) leading to the maximum sustainable yield.}
#'   \item{"MSY"}{ The yield when the biomass is at Bmsy and the fishing mortality is at Fmsy, i.e. the maximum sustainable yield.}
#' }
#'
#' The above parameter values can be extracted from the fit.spict() results using get.par().
#' @param inp List of input variables as output by check.inp.
#' @param dbg Debugging option. Will print out runtime information useful for debugging if set to 1. Will print even more if set to 2.
#' @return A result report containing estimates of model parameters, random effects (biomass and fishing mortality), reference points (Fmsy, Bmsy, MSY) including uncertainties given as standard deviations.
#' @export
#' @examples
#' data(pol)
#' rep <- fit.spict(pol$albacore)
#' summary(rep)
#' plot(rep)
#' @import TMB
fit.spict <- function(inp, dbg=0){
    # Check input list
    #if(!'checked' %in% names(inp)) inp <- check.inp(inp)
    #if(!inp$checked)
    inp <- check.inp(inp)
    datin <- list(delay=inp$delay, dt=inp$dt, dtpredcinds=inp$dtpredcinds, dtpredcnsteps=inp$dtpredcnsteps, dtprediind=inp$dtprediind, obsC=inp$obsC, ic=inp$ic, nc=inp$nc, I=inp$obsIin, ii=inp$iiin, iq=inp$iqin, ir=inp$ir, ffac=inp$ffaceuler, indpred=inp$indpred, robflagc=inp$robflagc, robflagi=inp$robflagi, lamperti=inp$lamperti, euler=inp$euler, dbg=dbg)
    pl <- inp$ini
    for(i in 1:inp$nphases){
        if(inp$nphases>1) cat(paste('Estimating - phase',i,'\n'))
        obj <- TMB::MakeADFun(data=datin, parameters=pl, random=inp$RE, DLL=inp$scriptname, hessian=TRUE, map=inp$map[[i]])
        config(trace.optimize=0, DLL=inp$scriptname)
        verbose <- FALSE
        obj$env$tracemgc <- verbose
        obj$env$inner.control$trace <- verbose
        obj$env$silent <- ! verbose
        obj$fn(obj$par)
        if(dbg==0){
            opt <- try(nlminb(obj$par, obj$fn, obj$gr))
            pl <- obj$env$parList(opt$par)
        }
    }
    rep <- NULL
    if(dbg==0){
        opt <- try(nlminb(obj$par, obj$fn, obj$gr))
        if(class(opt)!='try-error'){
            # Results
            rep <- try(TMB::sdreport(obj))
            if(class(rep)=='try-error'){
                rep <- list()
                rep$sderr <- 1
                rep$par.fixed <- opt$par
                rep$cov.fixed <- matrix(NA, length(opt$par), length(opt$par))
                cat('WARNING: Could not calculate sdreport.\n')
            }
            rep$pl <- obj$env$parList(opt$par)
            obj$fn()
            rep$Cp <- obj$report()$Cp
            rep$inp <- inp
            rep$opt <- opt
            #  - Calculate statistics -
            rep$stats <- list()
            K <- get.par('logK', rep, exp=TRUE)[2]
            r <- get.par('logr', rep, exp=TRUE)[2]
            p <- get.par('logp', rep, exp=TRUE)[2]
            Best <- get.par('logB', rep, exp=TRUE)
            Bests <- Best[rep$inp$indest, 2]
            # R-square
            Pest <- get.par('P', rep)
            B <- Best[-inp$ns, 2]
            inds <- unique(c(inp$ic[1:inp$nobsC], unlist(inp$ii)))
            gr <- r*(1-(B[inds]/K)^p) # Pella-Tomlinson intrinsic growth
            grobs <- Pest[inds, 2]/inp$dt[inds]/B[inds]
            ssqobs <- sum((grobs - mean(grobs))^2)
            ssqres <- sum((grobs - gr)^2)
            rep$stats$pseudoRsq <- 1 - ssqres/ssqobs
            # Prager's nearness
            Bmsy <- get.par('logBmsy', rep, exp=TRUE)[2]
            Bdiff <- Bests - Bmsy
            if(any(diff(sign(Bdiff))!=0)){
                rep$stats$nearness <- 1
            } else {
                ind <- which.min(abs(Bdiff))
                Bstar <- Bests[ind]
                if(Bstar > Bmsy){
                    rep$stats$nearness <- (K-Bstar)/(K-Bmsy)
                } else {
                    rep$stats$nearness <- Bstar/Bmsy
                }
            }
            # Prager's coverage
            rep$stats$coverage <- diff(range(Bests))/K
        } else {
            stop('Could not fit model, try changing the initial parameter guess in inp$ini. Error msg:', opt)
        }
    }
    if(!is.null(rep)) class(rep) <- "spictcls"
    return(rep)
}


#' @name calc.osa.resid
#' @title Calculate one-step-ahead residuals.
#' @details In TMB one-step-ahead residuals are calculated by sequentially including one data point at a time while keeping the model parameters fixed at their ML estimates. The calculated residuals are tested for independence in lag 1 using the Ljung-Box test (see Box.test).
#' @param rep A result report as generated by running fit.spict.
#' @param dbg Debugging option. Will print out runtime information useful for debugging if set to 1. Will print even more if set to 2.
#' @return An updated result report, which contains one-step-ahead residuals stored in the $osar variable.
#' @export
#' @examples
#' data(pol)
#' rep <- fit.spict(pol$albacore)
#' rep <- calc.osa.resid2(rep)
#' plotspict.osar(rep)
#' @import TMB
calc.osa.resid <- function(rep, dbg=0){
    inp <- rep$inp
    inp$ffac <- 1
    if("logF" %in% names(inp$ini)) inp$ini$logF <- NULL
    if("logB" %in% names(inp$ini)) inp$ini$logB <- NULL
    if("ns" %in% names(inp)) inp$ns <- NULL
    predmap <- get.predmap(rep$pl, inp$RE)
    plnew <- rep$pl
    logCpred <- rep(0, inp$nobsC-1)
    logIpred <- list()

    delay <- inp$delay*inp$dteuler # Delay in years
    timeobs <- inp$timeC
    for(i in 1:inp$nindex){
        logIpred[[i]] <- rep(0, inp$nobsI[i]-1)
        timeobs <- c(timeobs, inp$timeI[[i]])
    }
    timeobs <- sort(unique(timeobs))
    inds <- which(timeobs >= timeobs[1]+delay)
    timepred <- timeobs[inds]
    npred <- length(timepred)
    nser <- inp$nindex+1
    obsmat <- matrix(0, npred, nser)
    haveobs <- matrix(0, npred, nser)
    rm.na <- function(rr) rr[!is.na(rr)]
    icm <- rm.na(match(inp$timeC, timepred))
    ic <- rm.na(match(timepred, inp$timeC))
    obsmat[icm, 1] <- inp$obsC[ic]
    haveobs[icm, 1] <- 1
    for(i in 1:inp$nindex){
        iim <- rm.na(match(inp$timeI[[i]], timepred))
        ii <- rm.na(match(timepred, inp$timeI[[i]]))
        obsmat[iim, i+1] <- inp$obsI[[i]][ii]
        haveobs[iim, i+1] <- 1
    }
    predmat <- matrix(0, npred, nser)

    for(j in 1:npred){
        inp2 <- inp
        endtimes <- rep(-1, nser)
        # Catch
        cind <- which(inp$timeC < timepred[j])
        inp2$timeC <- inp$timeC[cind]
        if(length(cind)>0) endtimes[1] <- tail(inp2$timeC,1)
        inp2$obsC <- inp$obsC[cind]
        inp2$dtc <- inp$dtc[cind]
        cindm <- which(inp$timeC == timepred[j])
        if(length(cindm)==1) obsmat[j, 1] <- inp$obsC[cindm]
        # Indices
        for(i in 1:inp$nindex){
            iind <- which(inp$timeI[[i]] < timepred[j])
            inp2$timeI[[i]] <- inp$timeI[[i]][iind]
            if(length(iind)>0) endtimes[i+1] <- tail(inp2$timeI[[i]],1)
            inp2$obsI[[i]] <- inp$obsI[[i]][iind]
            iindm <- which(inp$timeI[[i]] == timepred[j])
            if(length(iindm)==1) obsmat[j, i+1] <- inp$obsI[[i]][iindm]
        }
        inp2$dtpredc <- timepred[j] - max(endtimes)
        inp2$dtpredi <- timepred[j] - max(endtimes)
        if(haveobs[j, 1] == 1) if(inp2$dtpredc < inp$dtc[which(inp$timeC == timepred[j])]) stop('Cannot calculate OSAR because index has a finer time step than catch. This needs to be implemented!')
        inp2 <- check.inp(inp2)
        datnew <- list(delay=inp2$delay, dt=inp2$dt, dtpredcinds=inp2$dtpredcinds, dtpredcnsteps=inp2$dtpredcnsteps, dtprediind=inp2$dtprediind, obsC=inp2$obsC, ic=inp2$ic, nc=inp2$nc, I=inp2$obsIin, ii=inp2$iiin, iq=inp2$iqin, ir=inp$ir, ffac=inp$ffac, indpred=inp2$indpred, robflagc=inp2$robflagc, robflagi=inp2$robflagi, lamperti=inp2$lamperti, euler=inp2$euler, dbg=dbg)
        for(k in 1:length(inp2$RE)) plnew[[inp2$RE[k]]] <- rep$pl[[inp2$RE[k]]][1:inp2$ns]
        objpred <- TMB::MakeADFun(data=datnew, parameters=plnew, map=predmap, random=inp2$RE, DLL=inp2$scriptname, hessian=TRUE, tracemgc=FALSE)
        verbose <- FALSE
        objpred$env$tracemgc <- verbose
        objpred$env$inner.control$trace <- verbose
        objpred$env$silent <- ! verbose
        objpred$fn()
        if(haveobs[j, 1] == 1) predmat[j, 1] <- objpred$report()$Cp
        for(i in 1:inp$nindex) if(haveobs[j, i+1] == 1) predmat[j, i+1] <- exp(objpred$report()$logIp[i])
    }
    # Catches
    inds <- which(haveobs[, 1]==1)
    timeC <- timepred[inds]
    logCpred <- log(predmat[inds, 1])
    logCpres <- log(obsmat[inds, 1]) - logCpred
    logCpboxtest <- Box.test(logCpres, lag=1, type='Ljung-Box')
    rep$stats$osarCpval <- logCpboxtest$p.value
    # Indices
    timeI <- list()
    logIpred <- list()
    logIpres <- list()
    logIpboxtest <- list()
    pvals <- rep(0, inp$nindex)
    names(pvals) <- paste0('I', 1:inp$nindex)
    for(i in 1:inp$nindex){
        inds <- which(haveobs[, i+1]==1)
        timeI[[i]] <- timepred[inds]
        logIpred[[i]] <- log(predmat[inds, i+1])
        logIpres[[i]] <- log(obsmat[inds, i+1]) - logIpred[[i]]
        logIpboxtest[[i]] <- Box.test(logIpres[[i]], lag=1, type='Ljung-Box')
        pvals[i] <- logIpboxtest[[i]]$p.value
    }
    rep$stats$osarIpval <- pvals
    rep$osar <- list(logCpres=logCpres, logCpred=logCpred, timeC=timeC, logCpboxtest=logCpboxtest, timeI=timeI, logIpres=logIpres, logIpred=logIpred, logIpboxtest=logIpboxtest)
    return(rep)
}


#' @name get.par
#' @title Extract parameters from a result report as generated by fit.spict.
#' @details Helper function for extracting the value and uncertainty of a specific model parameter, random effect or derived quantity.
#' @param parname Character string containing the name of the variable of interest.
#' @param rep A result report as generated by running fit.spict.
#' @param exp Take exp of the variable? TRUE/FALSE.
#' @param random DUMMY not used anymore. (Is the variable a random effect? TRUE/FALSE.)
#' @param fixed DUMMY not used anymore. (Is the variable a fixed effect? TRUE/FALSE.)
#' @return A matrix with four columns containing respectively: 1) the lower 95% confidence limit; 2) the parameter estimate; 3) the upper 95% confidence limit; 4) the parameter standard deviation in the domain it was estimated (log or non-log).
#' @export
#' @examples
#' data(pol)
#' rep <- fit.spict(pol$albacore)
#' Bmsy <- get.par('logBmsy', rep, exp=TRUE)
#' Best <- get.par('logB', rep, exp=TRUE)
#' K <- get.par('logK', rep, exp=TRUE)
get.par <- function(parname, rep=rep, exp=FALSE, random=FALSE, fixed=FALSE){
    if(!'sderr' %in% names(rep)){
        indran <- which(names(rep$par.random)==parname)
        indfix <- which(names(rep$par.fixed)==parname)
        indsdr <- which(names(rep$value)==parname)
        indopt <- which(names(rep$opt$par)==parname)
        est <- NULL
        if(length(indran)>0){
            est <- rep$par.random[indran]
            sd <- sqrt(rep$diag.cov.random[indran])
            ll <- est - 1.96*sd
            ul <- est + 1.96*sd
        }
        if(length(indfix)>0){
            est <- rep$par.fixed[indfix]
            sd <- sqrt(diag(rep$cov.fixed))[indfix]
            ll <- est - 1.96*sd
            ul <- est + 1.96*sd
        }
        if(length(indsdr)>0){
            est <- rep$value[indsdr]
            sd <- rep$sd[indsdr]
            ll <- est - 1.96*sd
            ul <- est + 1.96*sd
        }
        if(length(est)==0){
            ll <- NA
            ul <- NA
            sd <- NA
            est <- NA
            if(length(indopt)>0){
                est <- rep$opt$par[indopt]
            } else {
                if('phases' %in% names(rep$inp)){
                    if(rep$inp$phases[[parname]] == -1){
                        cat(paste('WARNING: did not estimate', parname, 'extracting fixed initial value.\n'))
                        est <- rep$inp$ini[[parname]]
                        ll <- est
                        ul <- est
                    }
                } else {
                    cat(paste('WARNING: could not extract', parname, '\n'))
                }
            }
        }
        if(exp){
            return(cbind(ll=exp(ll), est=exp(est), ul=exp(ul), sd))
        } else {
            return(cbind(ll, est, ul, sd))
        }
    }
}


#' @name make.ellipse
#' @title Calculate confidence ellipsis.
#' @details Calculates the confidence ellipsis of two reported model parameters. This is particularly useful as a detailed view of the uncertainty of two correlated parameters.
#' @param inds Indices of the two reported model parameters.
#' @param rep A result report as generated by running fit.spict.
#' @return A matrix with two columns containing the x and y coordinates of the ellipsis.
make.ellipse <- function(inds, rep){
    require(ellipse)
    covBF <- rep$cov[inds,inds]
    corBF <- cov2cor(covBF)
    parBF <- rep$value[inds]
    return(ellipse(corBF[1,2], scale=sqrt(diag(covBF)), centre=parBF))
}


#' @name arrow.line
#' @title Draw a line with arrow heads.
#' @details Add to an existing plot a continuous line with arrow heads showing the direction between each data point
#' @param x X coordinates.
#' @param y Y coordinates.
#' @param length See documentation for arrows.
#' @param angle See documentation for arrows.
#' @param code See documentation for arrows.
#' @param col See documentation for arrows.
#' @param lty See documentation for arrows.
#' @param lwd See documentation for arrows.
#' @param ... See documentation for arrows.
#' @return Nothing, but an arrow line is added to the current plot.
arrow.line <- function(x, y, length = 0.25, angle = 30, code = 2, col = par("fg"), lty = par("lty"), lwd = par("lwd"), ...){
    n <- length(x)
    for(i in 2:n) arrows(x[i-1], y[i-1], x[i], y[i], length, angle, code, col, lty, lwd, ...)
}


#' @name annual
#' @title Convert from quarterly (or other sub-annual) data to annual means.
#' @param inp An inp list.
#' @param vec The vector of values to convert to annual means
#' @return A list containing the annual means and a corresponding time vector.
annual <- function(intime, vec){
    anntime <- intime[which(intime %% 1 ==0)]
    nanntime <- length(anntime)
    nstepvec <- rep(0, nanntime)
    floortime <- floor(intime)
    for(i in 1:nanntime) nstepvec[i] <- sum(anntime[i]==floortime)
    nsteps <- max(nstepvec)
    # Remove years that are not full
    anntime <- anntime[which(nstepvec==max(nstepvec))]
    nanntime <- length(anntime)
    annvec <- rep(0, nanntime)
    for(i in 1:nanntime){
        inds <- which(anntime[i]==floortime)
        annvec[i] <- mean(vec[inds])
    }
    return(list(anntime=anntime, annvec=annvec))
}


#' @name plotspict.biomass
#' @title Plot estimated biomass.
#' @details Plots estimated biomass, Bmsy with confidence limits, and Binf.
#' @param rep A result report as generated by running fit.spict.
#' @param logax Take log of y-axis? default: FALSE
#' @return Nothing.
#' @examples
#' data(pol)
#' rep <- fit.spict(pol$albacore)
#' plotspict.biomass(rep)
#' @export
plotspict.biomass <- function(rep, logax=FALSE){
    if(!'sderr' %in% names(rep)){
        log <- ifelse(logax, 'y', '')
        inp <- rep$inp
        # Biomass plot
        Best <- get.par('logB', rep, exp=TRUE, random=TRUE)
        ns <- dim(Best)[1]
        Kest <- get.par('logK', rep, exp=TRUE, fixed=TRUE)
        Bmsy <- get.par('logBmsy', rep, exp=TRUE)
        Binf <- get.par('logBinf', rep, exp=TRUE)
        qest <- get.par('logq', rep, fixed=TRUE, exp=TRUE)
        BB <- get.par('logBBmsy', rep, exp=TRUE)
        inds <- which(is.na(Binf) | Binf<0)
        Binf[inds] <- 1e-12
        annlist <- annual(inp$time, Binf[, 2])
        Binftime <- annlist$anntime
        Binfs <- annlist$annvec
        Bp <- get.par('logBp', rep, exp=TRUE)
        scal <- 1
        cicol <- 'lightgray'
        obsI <- list()
        for(i in 1:inp$nindex) obsI[[i]] <- inp$obsI[[i]]/qest[i, 2]
        par(mar=c(5,4,4,4))
        ylim <- range(Best[, 1:3], Bp[1:3], unlist(obsI))/scal
        plot(inp$time, Best[,2]/scal, typ='n', xlab='Time', ylab=paste('Biomass'), main=paste('- Bmsy:',round(Bmsy[2]),' K:',round(Kest[2])), ylim=ylim, xlim=range(c(inp$time, tail(inp$time,1)+1)), log=log)
        axis(4, labels=pretty(ylim/Bmsy[2]), at=pretty(ylim/Bmsy[2])*Bmsy[2])
        mtext("B/Bmsy", side=4, las=0, line=2, cex=par('cex'))
        polygon(c(inp$time[1]-5,tail(inp$time,1)+5,tail(inp$time,1)+5,inp$time[1]-5), c(Bmsy[1],Bmsy[1],Bmsy[3],Bmsy[3]), col=cicol, border=cicol)
        cicol2 <- rgb(0, 0, 1, 0.1)
        polygon(c(inp$time, rev(inp$time)), c(BB[,1], rev(BB[,3]))/scal*Bmsy[2], col=cicol2, border=cicol2)
        abline(v=tail(inp$time[inp$indest],1), col='gray')
        for(i in 1:inp$nindex) points(inp$timeI[[i]], inp$obsI[[i]]/qest[i, 2], pch=i, cex=0.7)
        # Highlight influential index observations
        if('infl' %in% names(rep)){
            infl <- rep$infl$infl
            indslast <- inp$nobsC # Start after catch observations
            for(i in 1:inp$nindex){
                iinds <- indslast + 1:inp$nobsI[i]
                infl2 <- infl[iinds, ]
                cols <- apply(!is.na(infl2), 1, sum)
                ncols <- length(unique(cols))
                inds <- which(cols>0)
                points(inp$timeI[[i]][inds], inp$obsI[[i]][inds]/qest[i, 2], pch=21, cex=0.9, bg=cols[inds])
            }
        }
        if('true' %in% names(inp)){
            lines(inp$time, inp$true$B/scal, col='orange') # Plot true
            abline(h=inp$true$Bmsy, col='orange', lty=2)
        }
        lines(inp$time[inp$indest], Best[inp$indest,2]/scal, col='blue', lwd=1.5)
        lines(inp$time[inp$indpred], Best[inp$indpred,2]/scal, col='blue', lty=3)
        abline(h=Bmsy[2]/scal, col='black')
        # B CI
        lines(inp$time[inp$indest], Best[inp$indest,1]/scal, col=4, lty=2, lwd=1.5)
        lines(inp$time[inp$indest], Best[inp$indest,3]/scal, col=4, lty=2, lwd=1.5)
        lines(inp$time[inp$indpred], Best[inp$indpred,1]/scal, col=4, lty=2)
        lines(inp$time[inp$indpred], Best[inp$indpred,3]/scal, col=4, lty=2)
        # B/Bmsy CI
        cicol3 <- rgb(0, 0, 1, 0.2)
        lines(inp$time[inp$indest], BB[inp$indest,1]/scal*Bmsy[2], col=cicol3, lty=1, lwd=1)
        lines(inp$time[inp$indest], BB[inp$indest,3]/scal*Bmsy[2], col=cicol3, lty=1, lwd=1)
        lines(inp$time[inp$indpred], BB[inp$indpred,1]/scal*Bmsy[2], col=cicol3, lty=1, lwd=1)
        lines(inp$time[inp$indpred], BB[inp$indpred,3]/scal*Bmsy[2], col=cicol3, lty=1, lwd=1)
        lines(Binftime, Binfs/scal, col='green', lty=1)
        tp <- tail(inp$time,1)
        points(tp, tail(Best[,2],1)/scal, pch=21, bg='yellow')
        legend('topright', legend=c('Equilibrium',paste(tp,'pred.')), lty=c(1,NA), pch=c(NA,21), col=c('green',1), pt.bg=c(NA,'yellow'), bg='white')
        box()
    }
}


#' @name plotspict.bbmsy
#' @title Plot estimated B/Bmsy.
#' @details Plots estimated B/Bmsy.
#' @param rep A result report as generated by running fit.spict.
#' @param logax Take log of y-axis? default: FALSE
#' @return Nothing.
#' @examples
#' data(pol)
#' rep <- fit.spict(pol$albacore)
#' plotspict.bbmsy(rep)
#' @export
plotspict.bbmsy <- function(rep, logax=FALSE){
    if(!'sderr' %in% names(rep)){
        log <- ifelse(logax, 'y', '')
        inp <- rep$inp
        # Biomass plot
        Kest <- get.par('logK', rep, exp=TRUE, fixed=TRUE)
        Bmsy <- get.par('logBmsy', rep, exp=TRUE)
        Binf <- get.par('logBinf', rep, exp=TRUE)
        qest <- get.par('logq', rep, fixed=TRUE, exp=TRUE)
        BB <- get.par('logBBmsy', rep, exp=TRUE)
        ns <- dim(BB)[1]
        inds <- which(is.na(Binf) | Binf<0)
        cicol <- 'lightgray'
        obsI <- list()
        for(i in 1:inp$nindex) obsI[[i]] <- inp$obsI[[i]]/qest[i, 2]/Bmsy[2]
        par(mar=c(5,4,4,4))
        ylim <- range(c(BB[, 1:3], unlist(obsI)))
        plot(inp$time, BB[,2], typ='n', xlab='Time', ylab=paste('B/Bmsy'), ylim=ylim, xlim=range(c(inp$time, tail(inp$time,1)+1)), log=log, main='Relative biomass')
        cicol2 <- rgb(0, 0, 1, 0.1)
        polygon(c(inp$time, rev(inp$time)), c(BB[,1], rev(BB[,3])), col=cicol2, border=cicol2)
        abline(v=tail(inp$time[inp$indest],1), col='gray')
        for(i in 1:inp$nindex) points(inp$timeI[[i]], obsI[[i]], pch=i, cex=0.7)
        # Highlight influential index observations
        if('infl' %in% names(rep)){
            infl <- rep$infl$infl
            indslast <- inp$nobsC # Start after catch observations
            for(i in 1:inp$nindex){
                iinds <- indslast + 1:inp$nobsI[i]
                infl2 <- infl[iinds, ]
                cols <- apply(!is.na(infl2), 1, sum)
                ncols <- length(unique(cols))
                inds <- which(cols>0)
                points(inp$timeI[[i]][inds], inp$obsI[[i]][inds]/qest[i, 2]/Bmsy[2], pch=21, cex=0.9, bg=cols[inds])
            }
        }
        if('true' %in% names(inp)){
            lines(inp$time, inp$true$B/inp$true$Bmsy, col='orange') # Plot true
        }
        lines(inp$time[inp$indest], BB[inp$indest,2], col='blue', lwd=1.5)
        lines(inp$time[inp$indpred], BB[inp$indpred,2], col='blue', lty=3)
        cicol3 <- rgb(0, 0, 1, 0.2)
        lines(inp$time[inp$indest], BB[inp$indest,1], col=cicol3, lty=1, lwd=1)
        lines(inp$time[inp$indest], BB[inp$indest,3], col=cicol3, lty=1, lwd=1)
        lines(inp$time[inp$indpred], BB[inp$indpred,1], col=cicol3, lty=1, lwd=1)
        lines(inp$time[inp$indpred], BB[inp$indpred,3], col=cicol3, lty=1, lwd=1)
        abline(h=1)
        box()
    }
}


#' @name plotspict.osar
#' @title Plot one-step-ahead residuals
#' @details Plots observed versus predicted catches.
#' @param rep A result report as generated by running fit.spict.
#' @return Nothing.
#' @examples
#' data(pol)
#' rep <- fit.spict(pol$albacore)
#' rep <- calc.osa.resid(rep)
#' plotspict.osar(rep)
#' @export
plotspict.osar <- function(rep){
    inp <- rep$inp
    Cscal <- 1
    Cpred <- rep$osar$logCpred
    plot(rep$osar$logCpres, main=paste('Catch LB p-val:',round(rep$osar$logCpboxtest$p.value,4)), xlab='Time', ylab='OSA catch res.')
    abline(h=0, lty=3)
    #plot(inp$timeC, log(inp$obsC), typ='p', ylim=range(c(Cpred,log(inp$obsC),1.13*c(Cpred,log(inp$obsC)))), main=paste('Ljung-Box test p-value:',round(rep$osar$logCpboxtest$p.value,5)), ylab=paste('log Catch, Cscal:',Cscal), xlim=range(c(inp$timeC,inp$timeC[inp$nobsC]+1)), xlab='Time', col=1)
    #clr <- 'blue'
    #lines(rep$osar$timeC, Cpred, col=clr)
    #points(rep$osar$timeC, Cpred, pch=20, cex=0.7, col=clr)
    #legend('topright', c('Observations', 'OSA pred.'), lty=c(NA,1), col=c(1,clr), pch=c(1,20), pt.cex=c(1,0.7))
}



#' @name plotspict.f
#' @title Plot estimated fishing mortality.
#' @details Plots estimated fishing mortality with Fmsy and associated confidence interval.
#' @param rep A result report as generated by running fit.spict.
#' @param logax Take log of y-axis? default: FALSE
#' @return Nothing.
#' @examples
#' data(pol)
#' rep <- fit.spict(pol$albacore)
#' plotspict.f(rep)
#' @export
plotspict.f <- function(rep, logax=FALSE){
    if(!'sderr' %in% names(rep)){
        log <- ifelse(logax, 'y', '')
        inp <- rep$inp
        cicol <- 'lightgray'
        Fest <- get.par('logF', rep, exp=TRUE, random=TRUE)
        FF <- get.par('logFFmsy', rep, exp=TRUE)
        Fmsy <- get.par('logFmsy', rep, exp=TRUE)
        rest <- get.par('logr', rep, exp=TRUE, fixed=TRUE)
        rest <- apply(rest, 2, mean)

        if(min(inp$dtc) < 1){
        #if(TRUE){
            # Annual    
            al1 <- annual(inp$time, Fest[, 1])
            al2 <- annual(inp$time, Fest[, 2])
            al3 <- annual(inp$time, Fest[, 3])
            inds <- which(!is.na(al1$annvec) & al2$anntime<=tail(inp$time[inp$indest],1))
            time <- al1$anntime[inds]
            cl <- al1$annvec[inds]
            F <- al2$annvec[inds]
            cu <- al3$annvec[inds]
            inds <- which(!is.na(al1$annvec) & al2$anntime>=tail(inp$time[inp$indest],1))
            timep <- al1$anntime[inds]
            clp <- al1$annvec[inds]
            Fp <- al2$annvec[inds]
            cup <- al3$annvec[inds]
            al1f <- annual(inp$time, FF[, 1])
            al2f <- annual(inp$time, FF[, 2])
            al3f <- annual(inp$time, FF[, 3])
            inds <- which(!is.na(al1f$annvec))
            timef <- al1f$anntime[inds]
            clf <- al1f$annvec[inds]*Fmsy[2]
            Ff <- al2f$annvec[inds]*Fmsy[2]
            cuf <- al3f$annvec[inds]*Fmsy[2]
        } else {
            time <- inp$time[inp$indest]
            cl <- Fest[inp$indest, 1]
            F <- Fest[inp$indest, 2]
            cu <- Fest[inp$indest, 3]
            timep <- inp$time[inp$indpred]
            clp <- Fest[inp$indpred, 1]
            Fp <- Fest[inp$indpred, 2]
            cup <- Fest[inp$indpred, 3]
            timef <- inp$time
            clf <- FF[, 1]*Fmsy[2]
            Ff <- FF[, 2]*Fmsy[2]
            cuf <- FF[, 3]*Fmsy[2]
        }

        ylim <- range(c(cl, cu))
        plot(timef, Ff, typ='n', main=paste('Fmsy:',round(Fmsy[2],3),' ffac:',inp$ffac), ylim=ylim, col='blue', ylab='F', xlab='Time')
        axis(4, labels=pretty(ylim/Fmsy[2]), at=pretty(ylim/Fmsy[2])*Fmsy[2])
        mtext("F/Fmsy", side=4, las=0, line=2, cex=par('cex'))
        polygon(c(inp$time[1]-5,tail(inp$time,1)+5,tail(inp$time,1)+5,inp$time[1]-5), c(Fmsy[1],Fmsy[1],Fmsy[3],Fmsy[3]), col=cicol, border=cicol)
        cicol2 <- rgb(0, 0, 1, 0.1)
        polygon(c(timef, rev(timef)), c(clf, rev(cuf)), col=cicol2, border=cicol2)
        if(min(inp$dtc) < 1){
            lines(inp$time, Fest[, 2], col=rgb(0, 0, 0, 0.2))
        }
        abline(v=tail(inp$time[inp$indest],1), col='gray')
        if('true' %in% names(inp)){
            lines(inp$time, inp$true$F, col='orange') # Plot true
            abline(h=inp$true$Fmsy, col='orange', lty=2)
        }
        maincol <- 'blue'
        lines(time, cl, col=maincol, lwd=1.5, lty=2)
        lines(time, F, col=maincol, lwd=1.5)
        lines(time, cu, col=maincol, lwd=1.5, lty=2)
        lines(timep, clp, col=maincol, lty=2)
        lines(timep, Fp, col=maincol, lty=3)
        lines(timep, cup, col=maincol, lty=2)
        lines(timef, clf, col=rgb(0, 0, 1, 0.2))
        lines(timef, cuf, col=rgb(0, 0, 1, 0.2))
        points(tail(inp$time,1), tail(Fest[, 2],1), pch=21, bg='yellow')
        abline(h=Fmsy[2], col='black')
        abline(h=rest[2], col='red')
        legend('topleft',c(paste(tail(inp$time,1),'prediction')), pch=21, pt.bg=c('yellow'), bg='white')
        box()
    }
}


#' @name plotspict.fb
#' @title Plot fishing mortality versus biomass.
#' @details Plots estimated fishing mortality as a function of biomass together with reference points and the prediction for next year given a constant F. The equilibrium biomass for F fixed to the current value is also plotted.
#' @param rep A result report as generated by running fit.spict.
#' @param logax Take log of x and y-axes? default: FALSE
#' @return Nothing.
#' @examples
#' data(pol)
#' rep <- fit.spict(pol$albacore)
#' plotspict.fb(rep)
#' @export
plotspict.fb <- function(rep, logax=FALSE){
    if(!'sderr' %in% names(rep)){
        log <- ifelse(logax, 'xy', '')
        inp <- rep$inp
        scal <- 1
        cicol <- 'lightgray'
        Bp <- get.par('logBp', rep, exp=TRUE)
        Best <- get.par('logB', rep, exp=TRUE, random=TRUE)
        Bmsy <- get.par('logBmsy', rep, exp=TRUE)
        Binf <- get.par('logBinf', rep, exp=TRUE)
        ns <- dim(Best)[1]
        qest <- get.par('logq', rep, exp=TRUE, fixed=TRUE)
        rest <- get.par('logr', rep, exp=TRUE, fixed=TRUE)
        rest <- apply(rest, 2, mean)
        Fest <- get.par('logF', rep, exp=TRUE, random=TRUE)
        Fmsy <- get.par('logFmsy', rep, exp=TRUE)
        Fp <- Fest[ns,]
        inds <- c(which(names(rep$value)=='logBmsy'), which(names(rep$value)=='logFmsy'))
        cl <- make.ellipse(inds, rep)
        xlim <- range(c(tail(Binf[!is.na(Binf[,2]),2],1),exp(cl[,1]),Best[,2])/scal)
        ylim <- range(c(exp(cl[,2]),Fest[,2]))
        main <- paste('logalpha:',rep$pl$logalpha,'r:',round(rest[2],3),'q:',round(qest[2],3))
        par(mar=c(5,4,4,4))
        plot(Bmsy[2]/scal, Fmsy[2], typ='p', xlim=xlim, xlab='Biomass', ylab='F',  pch=24, bg='blue', ylim=ylim, log=log)
        axis(3, labels=pretty(xlim/Bmsy[2]), at=pretty(xlim/Bmsy[2])*Bmsy[2])
        mtext("B/Bmsy", side=3, las=0, line=2, cex=par('cex'))
        axis(4, labels=pretty(ylim/Fmsy[2]), at=pretty(ylim/Fmsy[2])*Fmsy[2])
        mtext("F/Fmsy", side=4, las=0, line=2, cex=par('cex'))
        alpha <- 0.15
        polygon(c(Bmsy[2], Bmsy[2], xlim[2], xlim[2]), c(Fmsy[2], 0, 0, Fmsy[2]), col=rgb(0,0.6,0,alpha), border=NA)
        polygon(c(Bmsy[2], Bmsy[2], 0, 0), c(Fmsy[2], 0, 0, Fmsy[2]), col=rgb(1,1,0,alpha), border=NA)
        polygon(c(Bmsy[2], Bmsy[2], xlim[2], xlim[2]), c(Fmsy[2], ylim[2], ylim[2], Fmsy[2]), col=rgb(1,1,0,alpha), border=NA)
        polygon(c(Bmsy[2], Bmsy[2], 0, 0), c(Fmsy[2], ylim[2], ylim[2], Fmsy[2]), col=rgb(0.6,0,0,alpha), border=NA)
        cicol2 <- 'gray'
        polygon(exp(cl[,1])/scal, exp(cl[,2]), col=cicol, border=cicol2)
        abline(h=Fmsy[2], lty=3)
        abline(v=Bmsy[2], lty=3)
        #arrow.line(Best[,2]/scal, Fest[,2], length=0.05, col='blue')
        if('true' %in% names(inp)){
            lines(inp$true$B/scal, inp$true$F, col='orange') # Plot true
            points(inp$true$Bmsy, inp$true$Fmsy, pch=24, bg='orange')
        }
        maincol <- rgb(0,0,1,0.8)
        if(min(inp$dtc) < 1){
        #if(FALSE){
            alf <- annual(inp$time, Fest[, 2])
            alb <- annual(inp$time, Best[, 2])
            lines(alb$annvec/scal, alf$annvec, col=maincol, lwd=1.5)
        } else {
            lines(Best[inp$indest,2]/scal, Fest[inp$indest,2], col=maincol, lwd=1.5)
            lines(Best[inp$indpred,2]/scal, Fest[inp$indpred,2], col=maincol, lty=3)
        }
        points(Bmsy[2]/scal, Fmsy[2], pch=24, bg='black')
        #lines(tail(Best[,2],1)/scal, tail(Fest[,2],1), col='blue', lty=3)
        #lines(c(tail(Best[,2],1), Bp[2])/scal, rep(Fp[2],2), col='blue', lty=3)
        points(tail(Best[,2],1)/scal, tail(Fest[,2],1), pch=21, bg='yellow')
        #points(Bp[2]/scal, Fp[2], pch=21, bg='yellow')
        points(tail(Binf[,2],1)/scal, Fp[2], pch=22, bg='green', cex=2)
        arrow.line(c(tail(Best[,2],1), tail(Binf[,2],1))/scal, rep(Fp[2],2), col='black', length=0.05)
        legend('topright', c('Estimated MSY',paste(tail(inp$time,1),'prediction'),'Equilibrium'), pch=c(24,21,22), pt.bg=c('black','yellow','green'), bg='white')
    }
}


#' @name plotspict.catch
#' @title Plot observed catch and predictions.
#' @details Plots observed catch and predictions using the current F and Fmsy. The plot also contains the equilibrium catch if the current F is maintained.
#' @param rep A result report as generated by running fit.spict.
#' @return Nothing.
#' @examples
#' data(pol)
#' rep <- fit.spict(pol$albacore)
#' plotspict.catch(rep)
#' @export
plotspict.catch <- function(rep){
    if(!'sderr' %in% names(rep)){
        inp <- rep$inp
        Cscal <- 1
        cicol <- 'lightgray'
        MSY <- get.par('MSY', rep, exp=FALSE)
        Cpredsub <- get.par('Cpredsub', rep)
        Pest <- get.par('P', rep)
        Cpredest <- get.par('logCpred', rep, exp=TRUE)
        Cpredest[Cpredest<0] <- 0
        rep$Cp[rep$Cp<0] <- 0
        indest <- which(inp$timeCp <= tail(inp$timeC,1))
        indpred <- which(inp$timeCp >= tail(inp$timeC,1))
        dtc <- inp$dtcp
        if(min(inp$dtc) < 1){
        #if(FALSE){
            alo <- annual(inp$timeC, inp$obsC/inp$dtc)
            timeo <- alo$anntime
            obs <- alo$annvec
            al2 <- annual(inp$timeCp[indest], Cpredest[indest, 2]/dtc[indest])
            inds <- which(!is.na(al2$annvec))
            time <- al2$anntime[inds]
            c <- al2$annvec[inds]
            al2p <- annual(inp$timeCp[indpred], Cpredest[indpred, 2]/dtc[indpred])
            inds <- which(!is.na(al2p$annvec))
            timep <- al2p$anntime[inds]
            cp <- al2p$annvec[inds]
            al1f <- annual(inp$timeCp, Cpredest[, 1]/dtc)
            al2f <- annual(inp$timeCp, Cpredest[, 2]/dtc)
            al3f <- annual(inp$timeCp, Cpredest[, 3]/dtc)
            inds <- which(!is.na(al2f$annvec))
            timef <- al2f$anntime[inds]
            clf <- al1f$annvec[inds]
            cf <- al2f$annvec[inds]
            cuf <- al3f$annvec[inds]
        } else {
            timeo <- inp$timeC
            obs <- inp$obsC/inp$dtc
            time <- inp$timeCp[indest]
            c <- Cpredest[indest, 2]/dtc[indest]
            timep <- inp$timeCp[indpred]
            cp <- Cpredest[indpred, 2]/dtc[indpred]
            timef <- inp$timeCp
            clf <- Cpredest[, 1]
            cf <- Cpredest[, 2]/dtc
            cuf <- Cpredest[, 3]
        }
        plot(time, c, typ='n', main=paste('MSY:',round(MSY[2]/Cscal)), xlab='Time', ylab=paste('Catch'), xlim=range(c(inp$time, tail(inp$time,1))), ylim=range(c(clf, cuf))/Cscal)
        polygon(c(inp$time[1]-5,tail(inp$time,1)+5,tail(inp$time,1)+5,inp$time[1]-5), c(MSY[1],MSY[1],MSY[3],MSY[3])/Cscal, col=cicol, border=cicol)
        cicol2 <- rgb(0, 0, 1, 0.1)
        polygon(c(timef, rev(timef)), c(clf, rev(cuf)), col=cicol2, border=cicol2)

        abline(v=tail(inp$timeC,1), col='gray')
        points(timeo, obs/Cscal)
        # Highlight influential index observations
        if('infl' %in% names(rep) & min(inp$dtc) == 1){
            infl <- rep$infl$infl[1:inp$nobsC, ]
            cols <- apply(!is.na(infl), 1, sum)
            ncols <- length(unique(cols))
            inds <- which(cols>0)
            points(inp$timeC[inds], inp$obsC[inds]/Cscal, pch=21, cex=0.9, bg=cols[inds])
        }
        if('true' %in% names(inp)) abline(h=inp$true$MSY, col='orange', lty=2)
        abline(h=MSY[2]/Cscal)
        lines(time, c, col=4, lwd=1.5)
        lines(timep, cp, col=4, lty=3)
        points(tail(inp$timeCp,1), tail(Cpredest[,2],1)/Cscal, pch=21, bg='yellow')
        legend('topleft',c(paste(tail(inp$timeCp,1),'prediction')), pch=21, pt.bg=c('yellow'), bg='white')
        box()
    }
}


#' @name plotspict.production
#' @title Plot theoretical production curve and estimates.
#' @details Plots the theoretical production curve (production as a function of biomass) as calculated from the estimated model parameters. Overlaid is the estimated production/biomass trajectory.
#' @param rep A result report as generated by running fit.spict.
#' @return Nothing.
#' @examples
#' data(pol)
#' rep <- fit.spict(pol$albacore)
#' plotspict.production(rep)
#' @export
plotspict.production <- function(rep){
    if(!'sderr' %in% names(rep)){
        inp <- rep$inp
        Best <- get.par('logB', rep, exp=TRUE, random=TRUE)
        Kest <- get.par('logK', rep, exp=TRUE, fixed=TRUE)
        rest <- get.par('logr', rep, exp=TRUE, fixed=TRUE)
        rest <- apply(rest, 2, mean)
        p <- get.par('logp', rep, exp=TRUE)
        Bmsy <- get.par('logBmsy', rep, exp=TRUE)
        Pest <- get.par('P', rep)
        nBplot <- 200
        Bplot <- seq(0.5*min(Best[, 2]), 1*max(Best[, 2]), length=nBplot)
        pfun <- function(r, K, p, B) r*B*(1 - (B/K)^p)
        Pst <- pfun(rest[2], Kest[2], p[2], Bplot)
        xlim <- range(Bplot/Kest[2])
        Bvec <- Best[-1, 2]
        dt <- inp$dt[-1]
        inde <- inp$indest[-length(inp$indest)]
        indp <- inp$indpred[-1]-1
        plot(Bvec[inde]/Kest[2], Pest[inde, 2]/dt[inde]/Bmsy[2], typ='l', ylim=range(Pest[,2]/inp$dt/Bmsy[2], Pst/Bmsy[2]), xlim=xlim, xlab='B/K', ylab='Production/Bmsy', col=4, main='Production curve')
        lines(Bvec[indp]/Kest[2], Pest[indp, 2]/dt[indp]/Bmsy[2], col=4, lty=3)
        lines(Bplot/Kest[2], Pst/Bmsy[2], col=1)
        #arrow.line(Best[-1, 2]/Bmsy[2], Pest[,2]/inp$dt/Bmsy[2], length=0.05)
        m <- p[2]+1
        mx <- (1/m)^(1/(m-1))
        abline(v=mx, lty=3)
        abline(h=0, lty=3)
    }
}


#' @name plotspict.growthrate
#' @title Plot intrinsic growth rate as a function of biomass.
#' @param rep A result report as generated by running fit.spict.
#' @return Nothing.
#' @examples
#' data(pol)
#' rep <- fit.spict(pol$albacore)
#' plotspict.growthrate(rep)
#' @export
plotspict.growthrate <- function(rep){
    if(!'sderr' %in% names(rep)){
        inp <- rep$inp
        Best <- get.par('logB', rep, exp=TRUE)
        Kest <- get.par('logK', rep, exp=TRUE)
        rest <- get.par('logr', rep, exp=TRUE)
        pest <- get.par('logp', rep, exp=TRUE)
        Pest <- get.par('P', rep)
        inds <- unique(c(inp$ic, unlist(inp$ii)))
        B <- Best[-inp$ns,2]
        r <- rest[2]
        p <- pest[2]
        K <- Kest[2]
        gr <- r*(1-(B[inds]/K)^p)
        grobs <- Pest[inds, 2]/inp$dt[inds]/B[inds]
        plot(B[inds], grobs, typ='p', xlab='B', ylab='Intrinsic growth rate', col='blue', main=bquote(pseudoR^2 == .(round(rep$stats$pseudoRsq, 4))))
        lines(B[inds], gr)
        abline(h=0, lty=3)
        legend('topright', legend=c('Observed growth', expression(r*(1-(B/K)^p))), col=c(4,1), pch=c(1,NA), lty=c(NA,1))
    }
}


#' @name plotspict.tc
#' @title Plot time constant.
#' @details Plots the time required for the biomass to reach a certain proportion of Bmsy. The time required to reach 95% of Bmsy is highlighted.
#' @param rep A result report as generated by running fit.spict.
#' @return Nothing.
#' @examples
#' data(pol)
#' rep <- fit.spict(pol$albacore)
#' plotspict.tc(rep)
#' @export
plotspict.tc <- function(rep){
    if(!'sderr' %in% names(rep)){
        inp <- rep$inp
        Best <- get.par('logB', rep, exp=TRUE)
        Fest <- get.par('logF', rep, exp=TRUE)
        Kest <- get.par('logK', rep, exp=TRUE)
        rest <- get.par('logr', rep, exp=TRUE)
        p <- get.par('logp', rep, exp=TRUE)
        sdbest <- get.par('logsdb', rep, exp=TRUE)
        Fmsy <- get.par('logFmsy', rep, exp=TRUE)
        Bmsy <- get.par('logBmsy', rep, exp=TRUE)
        B0cur <- Best[inp$indlastobs, 2]
        if(B0cur < Bmsy[2]) do.flag <- ifelse(B0cur/Bmsy[2]>0.95, FALSE, TRUE)
        if(B0cur > Bmsy[2]) do.flag <- ifelse(Bmsy[2]/B0cur>0.95, FALSE, TRUE)
        if(do.flag){
            if(B0cur < Bmsy[2]) facvec <- c(0, 0.75, 0.95, 1)
            if(B0cur > Bmsy[2]) facvec <- c(2, 1.25, 1.05, 1)
            Fvec <- round(facvec*Fmsy[2], digits=4)
            nFvec <- length(Fvec)
            g <- function(F, K, r, p, sdb, B0, dt, lamperti){
                if(lamperti){
                    return(exp(log(B0) + (r - r*(B0/K)^p - F - 0.5*sdb^2)*dt))
                } else {
                    return(B0 + B0*r*(1 - (B0/K)^p - F)*dt)
                }
            }
            simdt <- 0.01
            nt <- 5000
            Bsim <- matrix(0, nFvec, nt)
            time <- matrix(0, nFvec, nt)
            for(i in 1:nFvec){
                time[i, ] <- seq(0, simdt*(nt-1), by=simdt)
                Bsim[i, ] <- rep(0, nt)
                Bsim[i, 1] <- B0cur
                for(j in 2:nt){
                    Bsim[i, j] <- g(Fvec[i], Kest[2], rest[2], p[2], sdbest[2], Bsim[i, j-1], simdt, inp$lamperti)
                }
            }
            Bsim <- Bsim/Bmsy[2]
            frac <- 0.95
            if(B0cur < Bmsy[2]) inds <- which(Bsim[nFvec, ]<0.99)
            if(B0cur > Bmsy[2]) inds <- which(Bsim[nFvec, ]>(1/0.99))
            ylim <- range(Bsim[nFvec, ])
            #print(Bsim[nFvec, ])
            #print(ylim)
            plot(time[1, ], Bsim[1, ], typ='l', xlim=range(time[nFvec, inds]), ylim=ylim, col=3, ylab='Proportion of Bmsy', xlab='Years to Bmsy', main='Time to Bmsy')
            abline(h=c(frac, 1/frac), lty=1, col='lightgray')
            abline(h=1, lty=3)
            for(i in 2:nFvec) lines(time[i, ], Bsim[i, ], col=i+2)
            vt <- rep(0, nFvec)
            if(B0cur < Bmsy[2]) for(i in 1:nFvec) vt[i] <- time[i, max(which(Bsim[i, ]<frac))]
            if(B0cur > Bmsy[2]) for(i in 1:nFvec) vt[i] <- time[i, max(which(1/Bsim[i, ]<frac))]
            for(i in 1:nFvec) abline(v=vt[i], col=i+2, lty=2)
            lgnplace <- 'bottomright'
            if(B0cur > Bmsy[2]) lgnplace <- 'topright'
            legend(lgnplace, legend=paste('F =',facvec,'x Fmsy'), lty=1, col=2+(1:nFvec), lwd=rep(1,nFvec), bg='white')
            points(vt, rep(par('usr')[3], nFvec), col=3:(nFvec+2), pch=4)
        }
    }
}



#' @name plot.spictcls
#' @title 3x3 plot illustrating spict results.
#' @details Create a 3x3 plot containing the following:
#' \itemize{
#'  \item{1. Biomass using plotspict.biomass().}
#'  \item{2. One-step-ahead residuals, only if calculated, using plotspict.osar().}
#'  \item{3. One-step-ahead auto-correlation function (only if calculated).}
#'  \item{4. Estimated F versus estimated B using plotspict.fb().}
#'  \item{5. Estimated fishing mortality using plotspict.f().}
#'  \item{6. Observed versus predicted catches using plotspict.catch().}
#'  \item{7. Observed versus theoretical production using plotspict.production().}
#'  \item{8. Observed versus theoretical growth rate using plotspict.growthrate().}
#'  \item{9. Calculated time-constant using plotspict.tc().}
#' }
#' @param rep A result report as generated by running fit.spict.
#' @param logax Take log of relevant axes? default: FALSE
#' @return Nothing.
#' @examples
#' data(pol)
#' rep <- fit.spict(pol$albacore)
#' plot(rep)
#' @export
plot.spictcls <- function(rep, logax=FALSE){
    inp <- rep$inp
    #dev.new(width=10, height=10)
    #if('osar' %in% names(rep)){
        par(mfrow=c(4, 3))
    #} else {
    #    par(mfrow=c(2, 3))
    #}
    # Biomass
    plotspict.biomass(rep, logax=logax)
    # B/Bmsy
    plotspict.bbmsy(rep, logax=logax)
    # F
    plotspict.f(rep, logax=logax)
    # F versus B
    plotspict.fb(rep, logax=logax)
    # Catch
    plotspict.catch(rep)
    # Production curve
    plotspict.production(rep)
    # Intrinsic growth rate
    plotspict.growthrate(rep)
    # Time constant
    plotspict.tc(rep)
    if('osar' %in% names(rep)){
        # One-step-ahead catch residuals
        plotspict.osar(rep)
        acf(rep$osar$logCpres, main='ACF of catch OSAR')
    }
    if('infl' %in% names(rep)){
        # Plot influence summary
        plotspict.inflsum(rep)
    }
}


#' @name summary.spictcls
#' @title Output a summary of a fit.spict() run.
#' @details The output includes, the convergence message from the optimiser, the likelihood value of the parameters, the parameter estimates with 95% confidence intervals, estimates of derived parameters (BMsy, Fmsy, MSY) with 95% confidence intervals, and predictions of biomass, fishing mortality, and catch for the value of inp$dtpred.
#' @param object A result report as generated by running fit.spict.
#' @param numdigits Present values with this number of digits after the dot.
#' @return Nothing.
#' @examples
#' data(pol)
#' rep <- fit.spict(pol$albacore)
#' summary(rep)
#' @export
summary.spictcls <- function(object, numdigits=8){
    #options("scipen"=-100, "digits"=numdigits)
    options("scipen"=-3)
    rep <- object
    cat(paste('Convergence: ', rep$opt$convergence, '  MSG: ', rep$opt$message, '\n', sep=''))
    if(rep$opt$convergence>0) cat('WARNING: Model did not obtain proper convergence! Estimates and uncertainties are most likely invalid and should not be used.\n')
    if('sderr' %in% names(rep)) cat('WARNING: Could not calculate standard deviations. The optimum found may be invalid. Proceed with caution.\n')
    cat(paste('Negative log likelihood: ', round(rep$opt$objective, numdigits), '\n', sep=''))
    cat('\nFit statistics\n')
    statout <- unlist(rep$stats)
    cat('', paste(capture.output(statout),' \n'))
    cat('\nModel parameter estimates w 95% CI \n')
    sd <- sqrt(diag(rep$cov.fixed))
    nms <- names(rep$par.fixed)
    loginds <- grep('log', nms)
    logp1inds <- grep('logp1',nms)
    logitinds <- grep('logit',nms)
    loginds <- setdiff(loginds, c(logp1inds, logitinds))
    est <- rep$par.fixed
    est[loginds] <- exp(est[loginds])
    est[logitinds] <- invlogit(est[logitinds])
    est[logp1inds] <- invlogp1(est[logp1inds])
    cilow <- rep$par.fixed-1.96*sd
    cilow[loginds] <- exp(cilow[loginds])
    cilow[logitinds] <- invlogit(cilow[logitinds])
    cilow[logp1inds] <- invlogp1(cilow[logp1inds])
    ciupp <- rep$par.fixed+1.96*sd
    ciupp[loginds] <- exp(ciupp[loginds])
    ciupp[logitinds] <- invlogit(ciupp[logitinds])
    ciupp[logp1inds] <- invlogp1(ciupp[logp1inds])
    if('true' %in% names(rep$inp)){
        npar <- length(nms)
        truepar <- rep(0, npar)
        for(i in 1:npar) truepar[i] <- rep$inp$true[[nms[i]]]
        truepar[loginds] <- exp(truepar[loginds])
        truepar[logitinds] <- invlogit(truepar[logitinds])
        truepar[logp1inds] <- invlogp1(truepar[logp1inds])
        ci <- rep(0, npar)
        for(i in 1:npar) ci[i] <- as.numeric(truepar[i] > cilow[i] & truepar[i] < ciupp[i])
        resout <- cbind(estimate=round(est,numdigits), true=round(truepar,numdigits), cilow=round(cilow,numdigits), ciupp=round(ciupp,numdigits), true.in.ci=ci, est.in.log=round(rep$par.fixed,numdigits))
    } else {
        resout <- cbind(estimate=round(est,numdigits), cilow=round(cilow,numdigits), ciupp=round(ciupp,numdigits), est.in.log=round(rep$par.fixed,numdigits))
    }
    #nms[loginds] <- substr(names(rep$par.fixed[loginds]),4,60)
    nms[loginds] <- sub('log', '', names(rep$par.fixed[loginds]))
    nms[logitinds] <- sub('logit', '', names(rep$par.fixed[logitinds]))
    nms[logp1inds] <- sub('logp1', '', names(rep$par.fixed[logp1inds]))
    rownames(resout) <- nms
    cat('', paste(capture.output(resout),' \n'))
    if(!'sderr' %in% names(rep)){
        # Derived estimates
        cat('\nDerived estimates w 95% CI\n')
        cat(' Deterministic\n')
        derout <- rbind(get.par(parname='logBmsyd', rep, exp=TRUE)[c(2,1,3,2)],
                        get.par(parname='logFmsyd', rep, exp=TRUE)[c(2,1,3,2)],
                        get.par(parname='MSYd', rep)[c(2,1,3,2)])
        derout[, 4] <- log(derout[, 4])
        derout <- round(derout, numdigits)
        colnames(derout) <- c('estimate', 'cilow', 'ciupp', 'est.in.log')
        rownames(derout) <- c('Bmsyd', 'Fmsyd', 'MSYd')
        if('true' %in% names(rep$inp)){
            trueder <- c(rep$inp$true$Bmsyd, rep$inp$true$Fmsyd, rep$inp$true$MSYd)
            cider <- rep(0, 3)
            for(i in 1:3) cider[i] <- as.numeric(trueder[i] > derout[i, 2] & trueder[i] < derout[i, 3])
            derout <- cbind(estimate=derout[, 1], true=round(trueder,numdigits), derout[, 2:3], true.in.ci=cider, est.in.log=derout[, 4])
        }
        cat('', paste(capture.output(derout),' \n'))
        # Stochastic derived estimates
        cat(' Stochastic\n')
        derout <- rbind(get.par(parname='logBmsy', rep, exp=TRUE)[c(2,1,3,2)],
                        get.par(parname='logFmsy', rep, exp=TRUE)[c(2,1,3,2)],
                        get.par(parname='MSY', rep)[c(2,1,3,2)])
        derout[, 4] <- log(derout[, 4])
        derout <- round(derout, numdigits)
        colnames(derout) <- c('estimate', 'cilow', 'ciupp', 'est.in.log')
        rownames(derout) <- c('Bmsy', 'Fmsy', 'MSY')
        if('true' %in% names(rep$inp)){
            trueder <- c(rep$inp$true$Bmsy, rep$inp$true$Fmsy, rep$inp$true$MSY)
            cider <- rep(0, 3)
            for(i in 1:3) cider[i] <- as.numeric(trueder[i] > derout[i, 2] & trueder[i] < derout[i, 3])
            derout <- cbind(estimate=derout[, 1], true=round(trueder,numdigits), derout[, 2:3], true.in.ci=cider, est.in.log=derout[, 4])
        }
        cat('', paste(capture.output(derout),' \n'))
        # States
        cat('\nStates w 95% CI \n')
        stateout <- rbind(
            get.par(parname='logBl', rep, exp=TRUE)[c(2,1,3,2)],
            get.par(parname='logFl', rep, exp=TRUE)[c(2,1,3,2)],
            get.par(parname='logBlBmsy', rep, exp=TRUE)[c(2,1,3,2)],
            get.par(parname='logFlFmsy', rep, exp=TRUE)[c(2,1,3,2)])
        stateout[, 4] <- log(stateout[, 4])
        stateout <- round(stateout, numdigits)
        colnames(stateout) <- c('state est.', 'cilow', 'ciupp', 'est.in.log')
        et <- tail(rep$inp$time[rep$inp$indest],1)
        rownames(stateout) <- c(paste0('B_',et), paste0('F_',et), paste0('B_',et,'/Bmsy'), paste0('F_',et,'/Fmsy'))
        cat('', paste(capture.output(stateout),' \n'))
        # Predictions
        cat('\nPredictions w 95% CI \n')
        predout <- rbind(
            get.par(parname='logBp', rep, exp=TRUE)[c(2,1,3,2)],
            get.par(parname='logFp', rep, exp=TRUE)[c(2,1,3,2)],
            get.par(parname='logBpBmsy', rep, exp=TRUE)[c(2,1,3,2)],
            get.par(parname='logFpFmsy', rep, exp=TRUE)[c(2,1,3,2)],
            tail(get.par(parname='logCpred', rep, exp=TRUE),1)[c(2,1,3,2)])
        #        get.par(parname='logCp', rep, exp=TRUE)[c(2,1,3,2)])
        predout[, 4] <- log(predout[, 4])
        predout <- round(predout, numdigits)
        colnames(predout) <- c('prediction', 'cilow', 'ciupp', 'est.in.log')
        et <- tail(rep$inp$time,1)
        rownames(predout) <- c(paste0('B_',et), paste0('F_',et), paste0('B_',et,'/Bmsy'), paste0('F_',et,'/Fmsy'), paste0('Catch_',tail(rep$inp$timeCp,1)))
        cat('', paste(capture.output(predout),' \n'))
    }
    #if('osar' %in% names(rep)){
    #    cat('\nOne-step-ahead residuals - Ljung-box test for independence of lag 1\n')
    #    cat(paste(' Catches p-value:', round(rep$osar$logCpboxtest$p.value, numdigits), '\n'))
    #    for(i in 1:rep$inp$nindex) cat(paste(' Index', i, 'p-value:', round(rep$osar$logIpboxtest[[i]]$p.value, numdigits), '\n'))
    #}
}



#' @name read.aspic
#' @title Reads ASPIC input file.
#' @details Reads an input file following the ASPIC 7 format described in the ASPIC manual (found here http://www.mhprager.com/aspic.html).
#' @param filename Path of the ASPIC input file.
#' @return A list of input variables that can be used as input to fit.spict().
#' @examples
#' \dontrun{
#' filename <- 'YFT-SSE.a7inp' # or some other ASPIC 7 input file
#' inp <- read.aspic(filename)
#' rep <- fit.spict(inp)
#' summary(rep)
#' plot(rep)
#' }
#' @export
read.aspic <- function(filename){
    rawdat <- readLines(filename)
    # Read meta data
    dat <- list()
    c <- 0
    i <- 0
    iq <- 0
    found <- FALSE
    while(c<16){
        i <- i+1
        found <- !substr(rawdat[i],1,1)=='#'
        c <- c + found
        if(c==1 & found) dat$version <- rawdat[i]
        if(c==2 & found) dat$title <- rawdat[i]
        # Number of observations
        if(c==5 & found) dat$nobs <- as.numeric(read.table(filename, skip=i-1, nrows=1)[1])
        # Number of data series
        if(c==5 & found){
            dat$nseries <- as.numeric(read.table(filename, skip=i-1, nrows=1)[2])
            dat$qini <- rep(0, dat$nseries)
        }
        # Initial parameters
        if(c==10 & found) dat$B1Kini <- as.numeric(read.table(filename, skip=i-1, nrows=1)[2])
        if(c==11 & found) dat$MSYini <- as.numeric(read.table(filename, skip=i-1, nrows=1)[2])
        if(c==12 & found) dat$Fmsyini <- as.numeric(read.table(filename, skip=i-1, nrows=1)[2])
        if('nseries' %in% names(dat)){
            if(c==(13+iq) & found & iq<dat$nseries){
                iq <- iq+1
                dat$qini[iq] <- as.numeric(read.table(filename, skip=i-1, nrows=1)[2])
            }
        }
        found <- FALSE
    }
    # Read actual observations
    obsdat <- list()
    idat <- which(rawdat=='DATA')
    type <- rep('', dat$nseries)
    name <- rep('', dat$nseries)
    for(j in 1:dat$nseries){
        c <- 0
        i <- 0
        while(c<2){
            i <- i+1
            c <- c + !substr(rawdat[idat+i],1,1)=='#'
            if(c==1) name[j] <- rawdat[idat+i]
            if(c==2) type[j] <- rawdat[idat+i]
        }
        datstart <- idat+i
        obsdat[[j]] <- read.table(filename, skip=datstart, sep='', nrows=dat$nobs)
        obsdat[[j]][obsdat[[j]]<0] <- NA
        if(type[j]=='CE'){
            obsdat[[j]] <- obsdat[[j]][, 1:3]
            colnames(obsdat[[j]]) <- c('time', 'effort', 'catch')
            obsdat[[j]] <- cbind(obsdat[[j]], cpue=obsdat[[j]]$catch/obsdat[[j]]$effort)
        }
        if(type[j]=='CC'){
            obsdat[[j]] <- obsdat[[j]][, 1:3]
            colnames(obsdat[[j]]) <- c('time', 'cpue', 'catch')
        }
        # ASPIC differentiates between these index types, but SPiCT does not (yet). See the ASPIC manual for details.
        if(type[j] %in% c('I0','I1','I2','B0','B1','B2')){ 
            obsdat[[j]] <- obsdat[[j]][, 1:2]
            colnames(obsdat[[j]]) <- c('time', 'index')
        }
        idat <- datstart + dat$nobs
    }
    # Convert to SPiCT data by storing as an inp list
    inp <- list()
    # Insert catches
    ind <- grep('C', type)
    inp$timeC <- obsdat[[ind]]$time[!is.na(obsdat[[ind]]$catch)]
    inp$obsC <- obsdat[[ind]]$catch[!is.na(obsdat[[ind]]$catch)]
    # Insertindices
    inp$timeI <- list()
    inp$obsI <- list()
    for(j in 1:dat$nseries){
        if(type[j] %in% c('CC', 'CE')){
            inp$timeI[[j]] <- obsdat[[j]]$time[!is.na(obsdat[[j]]$cpue)]
            inp$obsI[[j]] <- obsdat[[j]]$cpue[!is.na(obsdat[[j]]$cpue)]
        }
        if(type[j] %in% c('I0','I1','I2','B0','B1','B2')){
            inp$timeI[[j]] <- obsdat[[j]]$time[!is.na(obsdat[[j]]$index)]
            inp$obsI[[j]] <- obsdat[[j]]$index[!is.na(obsdat[[j]]$index)]
        }
    }
    # Insert initial values
    inp$ini <- list()
    inp$ini$logbkfrac <- log(dat$B1Kini)
    inp$ini$logr <- log(2*dat$Fmsyini)
    inp$ini$logK <- log(2*dat$MSYini/dat$Fmsyini)
    inp$ini$logq <- log(dat$qini)
    inp$ini$logsdf <- log(1)
    inp$ini$logsdb <- log(1)
    inp$lamperti <- 1
    inp$euler <- 1
    inp$dteuler <- 1
    return(inp)
}


#' @useDynLib spict

.onUnload <- function (lib) {
  library.dynam.unload("spict", lib)
}


#' @name calc.rate
#' @title Helper function for sim.spict().
#' @param r Intrinsic growth rate.
#' @param F Fishing mortality.
#' @param sdb Standard deviation of biomass process.
#' @return Effective growth rate.
calc.rate <- function(r, F, sdb=0) r - F - 0.5*sdb^2


#' @name calc.binf
#' @title Helper function for sim.spict().
#' @param K Carrying capacity.
#' @param F Fishing mortality.
#' @param r Intrinsic growth rate.
#' @param p Pella-Tomlinson exponent.
#' @param sdb Standard deviation of biomass process.
#' @param lamperti Optional logical, use lamperti transformation?
#' @return B_infinity.
calc.binf <- function(K, F, r, p, sdb=0, lamperti=FALSE){
    if(!lamperti) sdb <- 0
    Binf <- K * (1 - F/r)^(1/p) * (1 - (p+1)/2 / (1-(1.0-p*r+p*F)^2) * sdb^2);
    return(Binf)
}


#' @name predict.b
#' @title Helper function for sim.spict().
#' @param B0 Initial biomass.
#' @param Binf Equilibrium biomass.
#' @param F Fishing mortality.
#' @param r Intrinsic growth rate.
#' @param K Carrying capacity.
#' @param dt Time step.
#' @param sdb Standard deviation of biomass process.
#' @param lamperti Optional logical, use lamperti transformation?
#' @param euler Optional logical, use Euler approximation?
#' @return Predicted biomass at the end of dt.
predict.b <- function(B0, Binf, F, r, K, dt, sdb=0, lamperti=FALSE, euler=FALSE){
    if(euler) lamperti <- 1
    if(!lamperti) sdb <- 0
    rate <- calc.rate(r, F, sdb)
    if(euler){
        return(exp( log(B0) + (rate - r/K*B0)*dt )) # Euler
    } else {
        return(1 / ( 1/Binf + (1/B0 - 1/Binf) * exp(-rate*dt) )) # Approximative analytical
    }
}


#' @name predict.c
#' @title Helper function for sim.spict().
#' @param F Fishing mortality.
#' @param K Carrying capacity.
#' @param r Intrinsic growth rate.
#' @param B0 Initial biomass.
#' @param Binf Equilibrium biomass.
#' @param dt Time step.
#' @param sdb Standard deviation of biomass process.
#' @param lamperti Optional logical, use lamperti transformation?
#' @param euler Optional logical, use Euler approximation?
#' @return Predicted catch during dt.
predict.c <- function(F, K, r, B0, Binf, dt, sdb=0, lamperti=FALSE, euler=FALSE){
    if(euler) lamperti <- 1
    if(!lamperti) sdb <- 0
    rate <- calc.rate(r, F, sdb)
    if(euler){
        return(F*B0*dt)
    } else {
        return(K/r*F * log( 1 - B0/Binf * (1 - exp(rate*dt))))
    }
}


#' @name sim.spict
#' @title Simulate data from Pella-Tomlinson model
#' @details Simulates data using either manually specified parameters values or parameters estimated by fit.spict().
#'
#' Manual specification:
#' To specify parameters manually use the inp$ini format similar to when specifying initial values for running fit.spict(). Observations can be simulated at specific times using inp$timeC and inp$timeI. If these are not specified then the length of inp$obsC or inp$obsI is used to determine the number of observations of catches and indices respectively. If none of these are specified then nobs observations of catch and index will be simulated evenly distributed in time.
#'
#' Estimated parameters:
#' Simply take the output from a fit.spict() run and use as input to sim.spict().
#' 
#' @param input Either an inp list with an ini key (see ?check.inp) or a rep list where rep is the output of running fit.spict().
#' @param nobs Optional specification of the number of simulated observations.
#' @return A list containing the simulated data.
#' @examples
#' data(pol)
#' rep <- fit.spict(pol$albacore)
#' sim <- sim.spict(rep)
#' repsim <- fit.spict(sim)
#' summary(repsim) # Note true values are listed in the summary
#'
#' inp <- pol$albacore
#' inp$obsC <- NULL
#' inp$obsI <- NULL
#' sim2 <- sim.spict(inp, nobs=150)
#' @export
sim.spict <- function(input, nobs=100){
    # Check if input is a inp (initial values) or rep (results).
    if('par.fixed' %in% names(input)){
        cat('Detected input as a SPiCT result, proceeding...\n')
        rep <- input
        inp <- rep$inp
        pl <- rep$pl
        plin <- inp$ini
    } else {
        if('ini' %in% names(input)){
            cat('Detected input as a SPiCT inp, proceeding...\n')
            inp <- input
            nm <- names(inp)
            if(!'timeC' %in% nm){
                if(!'obsC' %in% nm){
                    inp$nobsC <- nobs
                } else {
                    inp$nobsC <- length(inp$obsC)
                }
                inp$timeC <- 1:inp$nobsC
            } else {
                inp$nobsC <- length(inp$timeC)
            }
            if(!'obsC' %in% nm) inp$obsC <- rep(10, inp$nobsC) # Insert dummy, required by check.inp().
            if(!'logq' %in% names(inp$ini)) stop('logq not specified in inp$ini!')
            inp$nindex <- length(inp$ini$logq)
            if(!'timeI' %in% nm){
                if(!'obsI' %in% nm){
                    inp$nobsI <- nobs
                } else {
                    if(class(inp$obsI)!='list'){
                        tmp <- inp$obsI
                        inp$obsI <- list()
                        inp$obsI[[1]] <- tmp
                    }
                    inp$nobsI <- rep(0, inp$nindex)
                    for(i in 1:inp$nindex) inp$nobsI[i] <- length(inp$obsI[[i]])
                }
                inp$timeI <- list()
                for(i in 1:inp$nindex) inp$timeI[[i]] <- 1:inp$nobsI[i]
            } else {
                inp$nobsI <- rep(0, inp$nindex)
                for(i in 1:inp$nindex) inp$nobsI[i] <- length(inp$timeI[[i]])
            }
            if(!'obsI' %in% nm){
                inp$obsI <- list()
                for(i in 1:inp$nindex) inp$obsI[[i]] <- rep(10, inp$nobsI[i]) # Insert dummy
            }
            plin <- inp$ini
            inp <- check.inp(inp)
            pl <- inp$ini
        } else {
            stop('Invalid input! use either an inp list or a fit.spict() result.')
        }
    }
    time <- inp$time
    euler <- inp$euler
    lamperti <- inp$lamperti

    # Calculate derived variables
    dt <- inp$dteuler
    nt <- length(time)
    B0 <- exp(pl$logbkfrac)*exp(pl$logK)
    F0 <- exp(pl$logF0)
    r <- exp(pl$logr)
    K <- exp(pl$logK)
    q <- exp(pl$logq)
    p <- exp(pl$logp)
    R <- mean(r*p/(p+1)) # Take mean in the case r is a vector
    sdb <- exp(pl$logsdb)
    sdf <- exp(pl$logsdf)
    alpha <- exp(pl$logalpha)
    beta <- exp(pl$logbeta)
    sdi <- alpha * sdb
    sdc <- beta * sdf

    # B[t] is biomass at the beginning of the time interval starting at time t
    # I[t] is an index of biomass (e.g. CPUE) at the beginning of the time interval starting at time t
    # P[t] is the accumulated biomass production over the interval starting at time t
    # F[t] is the constant fishing mortality during the interval starting at time t
    # C[t] is the catch removed during the interval starting at time t. 
    # obsC[j] is the catch removed over the interval starting at time j, this will typically be the accumulated catch over the year.

    # - Fishing mortality -
    flag <- TRUE
    while(flag){
        F <- c(F0, exp(log(F0) + cumsum(rnorm(nt-1, 0, sdf*sqrt(dt))))) # Fishing mortality
        flag <- any(F >= 1.5*max(r)) # Do this to avoid crazy F values
    }
    # - B_infinity
    Binf <- rep(0,nt)
    for(t in 1:nt) Binf[t] <- calc.binf(K, F[t], r[inp$ir[t]], p, sdb, lamperti)
    # - Biomass -
    B <- rep(0,nt)
    B[1] <- B0
    e <- exp(rnorm(nt-1, 0, sdb*sqrt(dt)))
    for(t in 2:nt) B[t] <- predict.b(B[t-1], Binf[t-1], F[t-1], r[inp$ir[t]], K, dt, sdb, lamperti, euler) * e[t-1]
    # - Catch -
    Csub <- rep(0,nt)
    for(t in 1:nt) Csub[t] <- predict.c(F[t], K, r[inp$ir[t]], B[t], Binf[t], dt, sdb, lamperti, euler)
    # - Production -
    Psub <- rep(0,nt)
    for(t in 2:nt) Psub[t-1] <- B[t] - B[t-1] + Csub[t-1]
    # - Catch observations -
    C <- rep(0, inp$nobsC)
    obsC <- rep(0, inp$nobsC)
    for(i in 1:inp$nobsC){
        for(j in 1:inp$nc[i]){
            ind <- inp$ic[i] + j-1
            C[i] <- C[i] + Csub[ind];
        }
        obsC[i] <- exp(log(C[i]) + rnorm(1, 0, sdc))
    }
    if('outliers' %in% names(inp)){
        if('noutC' %in% names(inp$outliers)){
            #if(!'facoutC' %in% names(inp$outliers)) inp$outliers$facoutC <- 20
            fac <- invlogp1(inp$ini$logp1robfac)
            inp$outliers$orgobsC <- obsC
            inp$outliers$indsoutC <- sample(1:inp$nobsC, inp$outliers$noutC)
            obsC[inp$outliers$indsoutC] <- exp(log(obsC[inp$outliers$indsoutC]) + rnorm(inp$outliers$noutC, 0, fac*sdc))
        }
    }
    # - Index observations -
    obsI <- list()
    errI <- list()
    for(I in 1:inp$nindex){
        obsI[[I]] <- rep(0, inp$nobsI[I])
        errI[[I]] <- rep(0, inp$nobsI[I])
        for(i in 1:inp$nobsI[I]){
            errI[[I]][i] <- rnorm(1, 0, sdi)
            logItrue <- log(q[I]) + log(B[inp$ii[[I]][i]])
            obsI[[I]][i] <- exp(logItrue + errI[[I]][i])
        }
    }
    if('outliers' %in% names(inp)){
        if('noutI' %in% names(inp$outliers)){
            if(length(inp$outliers$noutI)==1) inp$outliers$noutI <- rep(inp$outliers$noutI, inp$nindex)
            #if(!'facoutI' %in% names(inp$outliers)) inp$outliers$facoutI <- 20
            fac <- invlogp1(inp$ini$logp1robfac)
            #if(length(inp$outliers$facoutI)==1) inp$outliers$facoutI <- rep(inp$outliers$facoutI, inp$nindex)
            inp$outliers$orgobsI <- obsI
            inp$outliers$indsoutI <- list()
            for(i in 1:inp$nindex){
                inp$outliers$indsoutI[[i]] <- sample(1:inp$nobsI[i], inp$outliers$noutI[i])
                obsI[[i]][inp$outliers$indsoutI[[i]]] <- exp(log(obsI[[i]][inp$outliers$indsoutI[[i]]]) + rnorm(inp$outliers$noutI[i], 0, fac*sdi))
            }
        }
    }
    
    sim <- list()
    sim$obsC <- obsC
    sim$timeC <- inp$timeC
    sim$obsI <- obsI
    sim$timeI <- inp$timeI
    sim$ini <- plin #list(logr=log(r), logK=log(K), logq=log(q), logsdf=log(sdf), logsdb=log(sdb))
    sim$dteuler <- inp$dteuler
    sim$euler <- inp$euler
    sim$lamperti <- inp$lamperti
    sim$phases <- inp$phases
    sim$outliers <- inp$outliers
    sim$true <- pl
    sim$true$B <- B
    sim$true$F <- F
    sim$true$R <- R
    sim$true$Bmsyd <- K/((p+1)^(1/p))
    sim$true$Fmsyd <- R
    sim$true$MSYd <- sim$true$Bmsy * sim$true$Fmsy
    # From Bordet & Rivest (2014)
    sim$true$Bmsy <- K/(p+1)^(1/p) * (1- (1+R*(p-1)/2)/(R*(2-R)^2)*sdb^2)
    sim$true$Fmsy <- R - p*(1-R)*sdb^2/((2-R)^2)
    sim$true$MSY <- K*R/((p+1)^(1/p)) * (1 - (p+1)/2*sdb^2/(1-(1-R)^2))
    sim$true$errI <- errI
    sim$true$logB <- NULL
    sim$true$logF <- NULL
    return(sim)
}


#' @name validate.spict
#' @title Simulate data and reestimate parameters
#' @details Given input parameters simulate a number of data sets. Then estimate the parameters from the simulated data and compare with the true values. Specifically, the one-step-ahead residuals are checked for autocorrelation and the confidence intervals of the estimated Fmsy and Bmsy are checked for consistency.
#'
#' WARNING: One should simulate at least 50 data sets and preferably more than 100 to obtain reliable results. This will take some time (potentially hours).
#' @param inp An inp list with an ini key (see ?check.inp). If you want to use estimated parameters for the simulation create the inp$ini from the pl key of a result of fit.spict().
#' @param nsim Number of simulated data sets in each batch.
#' @param nobsvec Vector containing the number of simulated observations of each data set in each batch.
#' @param estinp The estimation uses the true parameters as starting guess. Other initial values to be used for estimation can be specified in estinp$ini.
#' @param backup Since this procedure can be slow a filename can be specified in backup where the most recent results will be available.
#' @return A list containing the results of the validation with the following keys:
#' \itemize{
#'  \item{"osarpvals"}{ P-values of the Ljun-Box test for uncorrelated one-step-ahead residuals.}
#'  \item{"*msyci"}{Logical. TRUE if the true value of B/Fmsy was inside the 95\% confidence interval for the estimate, otherwise FALSE}
#'  \item{"*msyciw"}{ Width of the 95\% confidence interval of the estimate of B/Fmsy.}
#' }
#' @examples
#' data(pol)
#' rep0 <- fit.spict(pol$albacore)
#' inp <- list()
#' inp$ini <- rep0$pl
#' set.seed(1234)
#' validate.spict(inp, nsim=10, nobsvec=c(30, 60), backup='validate.RData')
#' @export
validate.spict <- function(inp, nsim=50, nobsvec=c(15, 60, 240), estinp=NULL, backup=NULL){
    nnobsvec <- length(nobsvec)
    if('logF' %in% names(inp$ini)) inp$ini$logF <- NULL
    if('logB' %in% names(inp$ini)) inp$ini$logB <- NULL
    ss <- list()
    
    for(i in 1:nnobsvec){
        ss[[i]] <- list()
        for(j in 1:nsim){
            cat(paste(Sys.time(), '- validating:  i:', i, 'j:', j, '\n'))
            sim <- sim.spict(inp, nobs=nobsvec[i])
            if(!is.null(estinp)) sim$ini <- estinp$ini
            rep <- try(fit.spict(sim))
            if(!class(rep)=='try-error'){
                rep <- calc.osa.resid(rep)
                ss[[i]][[j]] <- extract.simstats(rep)
            }
            if(!is.null(backup)) save(ss, file=backup)
        }
    }
    return(ss)
}


#' @name extract.simstats
#' @title Extracts relevant statistics from the estimation of a simulated data set.
#' @details TBA
#' @param rep A result report as generated by running fit.spict.
#' @return A list containing the relevant statistics.
#' @examples
#' data(pol)
#' sim <- sim.spict(pol$albacore)
#' rep <- fit.spict(sim)
#' extract.simstats(rep)
#' @export
extract.simstats <- function(rep){
    if('true' %in% names(rep$inp)){
        ss <- list()
        # Convergence
        ss$conv <- rep$opt$convergence
        # Fit stats
        ss$stats <- rep$stats
        # OSA residuals p-values
        if('osar' %in% names(rep)) ss$osarpvalC <- rep$osar$logCpboxtest$p.value
        if('osar' %in% names(rep)) ss$osarpvalI <- rep$osar$logIpboxtest[[1]]$p.value
        # Fmsy estimate
        Fmsy <- get.par('logFmsy', rep, exp=TRUE)
        ss$Fmsyci <- rep$inp$true$Fmsy > Fmsy[1] & rep$inp$true$Fmsy < Fmsy[3]
        ss$Fmsyciw <- Fmsy[3] - Fmsy[1]
        # Bmsy estimate
        Bmsy <- get.par('logBmsy', rep, exp=TRUE)
        ss$Bmsyci <- rep$inp$true$Bmsy > Bmsy[1] & rep$inp$true$Bmsy < Bmsy[3]
        ss$Bmsyciw <- Bmsy[3] - Bmsy[1]
        # MSY estimate
        MSY <- get.par('MSY', rep)
        ss$MSYci <- rep$inp$true$MSY > MSY[1] & rep$inp$true$MSY < MSY[3]
        ss$MSYciw <- MSY[3] - MSY[1]
        # Final biomass estimate
        Best <- get.par('logB', rep, exp=TRUE, random=TRUE)
        ind <- which(tail(rep$inp$timeC, 1)==rep$inp$time)
        ss$Bci <- rep$inp$true$B[ind] > Best[ind, 1] & rep$inp$true$B[ind] < Best[ind, 3]
        ss$Bciw <- Best[ind, 3] - Best[ind, 1]
        # Final B/Bmsy estimate
        BB <- get.par('logBBmsy', rep, exp=TRUE)
        ind <- which(tail(rep$inp$timeC, 1)==rep$inp$time)
        ss$BBci <- rep$inp$true$B[ind]/rep$inp$true$Bmsy > BB[ind, 1] & rep$inp$true$B[ind]/rep$inp$true$Bmsy < BB[ind, 3]
        ss$BBciw <- BB[ind, 3] - BB[ind, 1]
        # Final F/Fmsy estimate
        FF <- get.par('logFFmsy', rep, exp=TRUE)
        ind <- which(tail(rep$inp$timeC, 1)==rep$inp$time)
        ss$FFci <- rep$inp$true$F[ind]/rep$inp$true$Fmsy > FF[ind, 1] & rep$inp$true$F[ind]/rep$inp$true$Fmsy < FF[ind, 3]
        ss$FFciw <- FF[ind, 3] - FF[ind, 1]
        return(ss)
    } else {
        stop('These results do not come from the estimation of a simulated data set!')
    }
}


#' @name write.aspic
#' @title Takes a SPiCT input list and writes it as an Aspic input file.
#' @details TBA
#' @param input List of input variables or the output of a simulation using sim.spict().
#' @param filename Name of the file to write.
#' @return Noting.
#' @examples
#' data(pol)
#' sim <- (pol$albacore)
#' write.aspic(sim)
#' @export
write.aspic <- function(input, filename='spictout.a7inp'){
    inp <- check.inp(input)
    timeobs <- sort(unique(c(inp$timeC, inp$timeI[[1]])))
    nobs <- length(timeobs)
    dat <- cbind(timeobs, rep(-1, nobs), rep(-1, nobs))
    inds <- match(inp$timeI[[1]], timeobs)
    dat[inds, 2] <- inp$obsI[[1]]
    inds <- match(inp$timeC, timeobs)
    dat[inds, 3] <- inp$obsC
    cat('Writing input to:', filename, '\n')
    cat('ASPIC-V7\n', file=filename)
    cat(paste('# File generated by SPiCT function write.aspic at', Sys.time(), '\n'), file=filename, append=TRUE)
    cat('"Unknown stock"\n', file=filename, append=TRUE)
    cat('# Program mode (FIT/BOT), verbosity, N bootstraps, [opt] user percentile:\n', file=filename, append=TRUE)
    cat(paste0(inp$aspic$mode, '  ', inp$aspic$verbosity, '  ', inp$aspic$nboot, '  ', inp$aspic$ciperc, '\n'), file=filename, append=TRUE)
    #cat('BOT  102  1000  95\n', file=filename, append=TRUE)
    cat('# Model shape, conditioning (YLD/EFT), obj. fn. (SSE/LAV/MLE/MAP):\n', file=filename, append=TRUE)
    cat('LOGISTIC  YLD  SSE\n', file=filename, append=TRUE)
    cat('# N years, N series:\n', file=filename, append=TRUE)
    cat(paste0(nobs, '  ', inp$nindex, '\n'), file=filename, append=TRUE)
    cat('# Monte Carlo mode (0/1/2), N trials:\n', file=filename, append=TRUE)
    cat('0  30000\n', file=filename, append=TRUE)
    cat('# Convergence criteria (3 values):\n', file=filename, append=TRUE)
    cat('1.00E-08  3.00E-08  1.00E-04\n', file=filename, append=TRUE)
    cat('# Maximum F, N restarts, [gen. model] N steps/yr:\n', file=filename, append=TRUE)
    cat('8.00E+00  6  24\n', file=filename, append=TRUE)
    cat('# Random seed (large integer):\n', file=filename, append=TRUE)
    cat('1234\n', file=filename, append=TRUE)
    cat('# Initial guesses and bounds follow:\n', file=filename, append=TRUE)
    estbkfrac <- 1
    if(!is.null(inp$phases$logbkfrac)) if(inp$phases$logbkfrac==-1) estbkfrac <- 0
    #estbkfrac <- 0
    cat(sprintf('B1K   %3.2E  %i  %3.2E  %3.2E  penalty  %3.2E\n', exp(inp$ini$logbkfrac), estbkfrac, 0.01*exp(inp$ini$logbkfrac), 100*exp(inp$ini$logbkfrac), 0), file=filename, append=TRUE)
    MSY <- exp(inp$ini$logK + inp$ini$logr)/4
    cat(sprintf('MSY   %3.2E  1  %3.2E  %3.2E\n', MSY, 0.03*MSY, 5000*MSY), file=filename, append=TRUE)
    Fmsy <- exp(inp$ini$logr)/2
    cat(sprintf('Fmsy  %3.2E  1  %3.2E  %3.2E\n', Fmsy, 0.01*Fmsy, 100*Fmsy), file=filename, append=TRUE)
    for(i in 1:inp$nindex) cat(sprintf('q     %3.2E  1  %3.2E  %3.2E  %3.2E\n', exp(inp$ini$logq[i]), 1, 0.001*exp(inp$ini$logq[i]), 100*exp(inp$ini$logq[i])), file=filename, append=TRUE)
    cat('DATA\n', file=filename, append=TRUE)
    cat('"Combined-Fleet Index, Total Landings"\n', file=filename, append=TRUE)
    cat('CC\n', file=filename, append=TRUE)
    for(i in 1:nobs) cat(sprintf('  %4i    % 6.4E    % 6.4E\n', timeobs[i], dat[i, 2], dat[i, 3]), file=filename, append=TRUE)
    if(inp$nindex>1){
        for(I in 2:inp$nindex){
            cat(paste0('"Index',I-1,'"\n'), file=filename, append=TRUE)
            cat('I1\n', file=filename, append=TRUE)
            dat <- cbind(timeobs, rep(-1, nobs))
            inds <- match(inp$timeI[[I]], timeobs)
            dat[inds, 2] <- inp$obsI[[I]]
            for(i in 1:nobs) cat(sprintf('  %4i    % 6.4E\n', dat[i ,1], dat[i, 2]), file=filename, append=TRUE)
        }
    }
}


#' @name read.aspic.res
#' @title Reads the parameter estimates of an Aspic result file.
#' @details TBA
#' @param filename Name of the Aspic result file to read
#' @return Vector containing the parameter estimates.
#' @export
read.aspic.res <- function(filename){
    out <- list()
    aspicres <- readLines(filename)
    ind <- grep('Number of years analyzed', aspicres)
    ind2 <- gregexpr('[0-9]+', aspicres[ind])
    nobs <- as.numeric(unlist(regmatches(aspicres[ind], ind2))[1])
    get.num <- function(str, aspicres){
        ind <- grep(str, aspicres)
        ind2 <- gregexpr('[0-9.]*E[+-][0-9]*', aspicres[ind])
        as.numeric(unlist(regmatches(aspicres[ind], ind2)))[1]
    }
    # Parameters
    bkfrac <- get.num('^B1/K', aspicres)
    MSY <- get.num('^MSY', aspicres)
    Fmsy <- get.num('^Fmsy', aspicres)
    q <- get.num('^q', aspicres)
    Bmsy <- MSY/Fmsy
    K <- 2*Bmsy
    r <- 2*Fmsy
    # States
    ind <- grep('ESTIMATED POPULATION TRAJECTORY', aspicres)
    out$states <- read.table(filename, skip=ind+6, sep='', nrows=nobs, strip.white=TRUE)
    out$states <- read.table(filename, skip=ind+6, sep='', nrows=nobs, strip.white=TRUE)
    colnames(out$states) <- c('obs', 'time', 'Fest', 'B0est', 'Best', 'Catch', 'Cest', 'Pest', 'FFmsy', 'BBmsy')
    out$pars <- c(r=r, K=K, q=q, Fmsy=Fmsy, Bmsy=Bmsy, MSY=MSY)
    return(out)
}


#' @name get.osar.pvals
#' @title Gets the p-values of one-step-ahead residuals for all data series.
#' @details TBA
#' @param rep A valid result from fit.spict().
#' @return A vector containing the p-values.
get.osar.pvals <- function(rep){
    osarI <- unlist(rep$osar$logIpboxtest)
    inds <- which(names(osarI)=='p.value')
    return(c(rep$osar$logCpboxtest$p.value, as.numeric(osarI[inds])))
}


#' @name calc.influence
#' @title Calculates influence statistics of observations.
#' @details TBA
#' @param rep A valid result from fit.spict().
#' @return A list equal to the input with the added key "infl" containing influence statistics.
#' @export
calc.influence <- function(rep){
    #parnams <- c('logFmsy', 'MSY', 'logCp', 'logBlBmsy')
    if(!'osar' %in% names(rep)) stop('Need to calculate one-step-ahead residuals before calculating influence statistics. Use the calc.osa.resid() function.')
    inp <- rep$inp
    #parnams <- c('logFmsy', 'MSY')
    parnams <- c('logr', 'logK', 'logq', 'logsdf', 'logsdb', 'MSY')
    parnamsnolog <- sub('log', '', parnams)
    np <- length(parnams)
    doexp <- rep(0, np)
    doexp[grep('log', parnams)] <- 1
    parmat <- matrix(0, np, 4)
    rownames(parmat) <- parnamsnolog
    for(k in 1:np) parmat[k, ] <- get.par(parnams[k], rep)
    detcov <- unlist(determinant(rep$cov.fixed, log=TRUE))[1]
    likval <- rep$opt$objective
    osarpvals <- get.osar.pvals(rep)
    rwnms <- paste0('C_', inp$timeC)
    for(i in 1:inp$nindex) rwnms <- c(rwnms, paste0('I', i, '_', inp$timeI[[i]]))
    nser <- inp$nindex+1
    nobs <- inp$nobsC + sum(inp$nobsI)
    pararr <- array(0, dim=c(np, nobs, 4))
    ddetcov <- rep(0, nobs)
    names(ddetcov) <- rwnms
    dosarpvals <- matrix(0, nobs, length(osarpvals))
    rownames(dosarpvals) <- rwnms
    sernames <- c('C', paste0('I', 1:inp$nindex))
    colnames(dosarpvals) <- sernames
    cat('Calculating influence statistics for', nobs, 'observations:\n')
    inflfun <- function(c, inp, parnams){
        #cat(c, '.. ')
        inp2 <- inp
        if(c <= inp$nobsC){
            # Loop over catches
            inp2$obsC <- inp$obsC[-c]
            inp2$timeC <- inp$timeC[-c]
            inp2$dtc <- inp$dtc[-c]
        } else {
            # Loop over indices
            breaks <- inp$nobsC + cumsum(c(0, inp$nobsI[-inp$nindex]))
            j <- sum(c > breaks)
            i <- c - breaks[j]
            if(class(inp$obsI)=='numeric'){
                inp2$obsI <- inp$obsI[-i]
                inp2$timeI <- inp$timeI[-i]
            } else {
                inp2$obsI[[j]] <- inp$obsI[[j]][-i]
                inp2$timeI[[j]] <- inp$timeI[[j]][-i]
            }
        }
        rep2 <- fit.spict(inp2)
        rep2 <- calc.osa.resid(rep2)
        # Calculate diagnostics
        np <- length(parnams)
        parmat <- matrix(0, np, 4)
        for(k in 1:np) parmat[k, ] <- get.par(parnams[k], rep2)
        detcov <- unlist(determinant(rep2$cov.fixed, log=TRUE))[1]
        dosarpvals <- get.osar.pvals(rep2)
        return(list(parmat=parmat, detcov=detcov, dosarpvals=dosarpvals))
    }
    # Calculate influence
    #partry <- try(library(parallel))
    partry <- try(library(multicore))
    if(class(partry)=='try-error'){
        single.inflfun <- function(c, inp, parnams, nobs){
            res <- inflfun(c, inp, parnams)
            ## Send progress update
            cat(sprintf("Progress: %.2f%%\r", c/nobs * 100))
            return(res)
        }

        res <- lapply(1:nobs, single.inflfun, inp, parnams, nobs)
    } else {
        multi.inflfun <- function(c, inp, parnams, nobs, progfile){
            res <- inflfun(c, inp, parnams)
            ## Send progress update
            if(class(progfile)[1]=='fifo') writeBin(1/nobs, progfile)
            return(res)
        }
        ## Open fifo progress file
        progfile <- fifo(tempfile(), open="w+b", blocking=T)
        if(inherits(fork(), "masterProcess")) {
        #if(inherits(parallel:::mcfork(), "masterProcess")) {
            ## Child
            progress <- 0.0
            while (progress < 1 && !isIncomplete(progfile)) {
                msg <- readBin(progfile, "double")
                progress <- progress + as.numeric(msg)
                cat(sprintf("Progress (multicore): %.2f%%\r", progress * 100))
            } 
            exit()
        }
        res <- mclapply(1:nobs, multi.inflfun, inp, parnams, nobs, progfile, mc.cores=4)
        ## Close progress file
        close(progfile)
    }
    cat('\n')
    # Gather results
    for(i in 1:nobs){
        ddetcov[i] <- detcov - res[[i]]$detcov
        dosarpvals[i, ] <- res[[i]]$dosarpvals
        for(k in 1:np) pararr[k, i, ] <- res[[i]]$parmat[k, ]
    }
    dfbeta <- matrix(0, nobs, np)
    colnames(dfbeta) <- parnamsnolog
    rownames(dfbeta) <- rwnms
    dpar <- dfbeta
    for(k in 1:np){
        if(doexp[k]==0){
            dpar[, k] <- parmat[k, 2]-pararr[k, , 2]
        } else {
            dpar[, k] <- exp(parmat[k, 2])-exp(pararr[k, , 2])
        }
        dfbeta[, k] <- (parmat[k, 2]-pararr[k, , 2])/pararr[k, , 4]
    }
    # - Calculate influence matrix -
    infl <- matrix(NA, nobs, 1+nser+np)
    rownames(infl) <- rwnms
    tnames <- c('CR', paste0('p', sernames), paste0('dfb', parnamsnolog))
    colnames(infl) <- tnames
    # COVRATIO
    dl <- 3*length(rep$par.fixed)/nobs
    covratio <- exp(ddetcov)
    ac <- abs(covratio-1)
    inds <- which(ac > dl)
    infl[inds, 1] <- 1
    # OSAR
    osarpvals <- get.osar.pvals(rep)
    alpha <- 0.05
    orgres <- osarpvals < alpha
    newres <- dosarpvals < alpha
    chgmat <- matrix(0, nobs, nser)
    for(i in 1:nser) chgmat[, i] <- newres[, i]-orgres[i]
    tmp <- matrix(NA, nobs, nser)
    inds <- which(chgmat!=0)
    tmp[inds] <- 1
    infl[, 2:(1+nser)] <- tmp
    # DFBETA
    tmp <- matrix(NA, nobs, np)
    adfbeta <- abs(dfbeta)
    al <- 2/sqrt(nobs)
    inds <- which(adfbeta>al)
    tmp[inds] <- 1
    infl[, (2+nser):(1+nser+np)] <- tmp
    # Compile influence
    ninfl <- dim(infl)[2]
    for(i in 1:ninfl){
        inds <- which(infl[, i]==1)
        infl[inds, i] <- i
    }
    rep$infl <- list(dfbeta=dfbeta, dpar=dpar, ddetcov=ddetcov, dosarpvals=dosarpvals, infl=infl)
    # Add to stats
    rep$stats$inflperc <- sum(apply(is.na(infl), 1, sum)<dim(infl)[2])/dim(infl)[1]
    return(rep)
}

#' @name put.ax
#' @title Adds the x-axis to influence plots
#' @details TBA
#' @param rep A valid result from calc.influence().
#' @return Nothing.
put.xax <- function(rep){
    inp <- rep$inp
    sernames <- colnames(rep$infl$dosarpvals)
    xs <- c(inp$nobsC, inp$nobsI)
    xat <- cumsum(xs[-length(xs)]+0.5)
    xmid <- c(0, xat) + xs/2
    axis(1, at=xat, labels='')
    axis(1, at=xmid, labels=sernames, tick=FALSE)
    abline(v=xat, lty=2, col='gray')
}


#' @name plotspict.infl
#' @title Plots influence statistics of observations.
#' @details TBA
#' @param rep A valid result from calc.influence().
#' @return Nothing.
#' @export
plotspict.infl <- function(rep){
    inp <- rep$inp
    #sernames <- c('C', paste0('I', 1:inp$nindex))
    dfbeta <- rep$infl$dfbeta
    dpar <- rep$infl$dpar
    ddetcov <- rep$infl$ddetcov
    dosarpvals <- rep$infl$dosarpvals
    sernames <- colnames(dosarpvals)
    nser <- inp$nindex+1
    nobs <- inp$nobsC + sum(inp$nobsI)
    parnams <- colnames(dfbeta)
    np <- length(parnams)
    rwnms <- rownames(dfbeta)
    infl <- rep$infl$infl
    par(mfrow=c(2,2))
    # Plot covratio
    dl <- 3*length(rep$par.fixed)/nobs
    covratio <- exp(ddetcov)
    ac <- abs(covratio-1)
    inds <- which(ac > dl)
    plot(ac, ylab='abs(covratio-1)', ylim=c(0, 1.05*max(ac)), xlab='', xaxt='n', main='COVRATIO')
    abline(h=dl)
    text(inds, ac[inds], rwnms[inds], pos=3, cex=0.7)
    put.xax(rep)
    # Plot OSAR p-values
    osarpvals <- get.osar.pvals(rep)
    alpha <- 0.05
    orgres <- osarpvals < alpha
    newres <- dosarpvals < alpha
    chgmat <- matrix(0, nobs, nser)
    for(i in 1:nser) chgmat[, i] <- newres[, i]-orgres[i]
    inds <- which(chgmat!=0)
    rc <- arrayInd(inds, dim(dosarpvals))
    nms <- rwnms[rc[, 1]]
    plot(dosarpvals[, 1], ylab='p-value', ylim=0:1, xlab='', xaxt='n', main='OSAR p-values')
    for(i in 2:nser) points(dosarpvals[, i], col=i)
    abline(h=osarpvals, col=1:nser)
    for(i in 1:length(inds)){
        text(rc[i, 1], dosarpvals[rc[i, 1], rc[i, 2]], nms[i], pos=3, cex=0.7, col=rc[i, 2])
    }
    abline(h=0.05, lwd=2, lty=3)
    legend('topleft', legend=sernames, pch=1, col=1:nser)
    put.xax(rep)
    # Plot dfbeta
    adfbeta <- abs(dfbeta)
    al <- 2/sqrt(nobs)
    inds <- which(adfbeta>al)
    rc <- arrayInd(inds, dim(adfbeta))
    nms <- rwnms[rc[, 1]]
    plot(adfbeta[, 1], ylim=c(0, 1.05*max(adfbeta)), ylab='abs dfbeta', xlab='', xaxt='n', main='DFBETA')
    for(i in 2:np) points(adfbeta[, i], col=i)
    abline(h=al, lwd=2, lty=3)
    for(i in 1:length(inds)){
        text(rc[i, 1], adfbeta[rc[i, 1], rc[i, 2]], nms[i], pos=3, cex=0.7, col=rc[i, 2])
        text(rc[i, 1], adfbeta[rc[i, 1], rc[i, 2]], round(dpar[rc[i, 1], rc[i, 2]], 3), pos=1, cex=0.7, col=rc[i, 2])
    }
    legend('topleft', legend=parnams, pch=1, col=1:np)
    put.xax(rep)
    # Plot of influence
    plotspict.inflsum(rep)
}


#' @name plotspict.inflsum
#' @title Plots summary of influence statistics of observations.
#' @details TBA
#' @param rep A valid result from calc.influence().
#' @return Nothing.
#' @export
plotspict.inflsum <- function(rep){
    infl <- rep$infl$infl
    nobs <- dim(infl)[1]
    ninfl <- dim(infl)[2]
    matplot(infl, pch=1, col=1, yaxt='n', ylab='', xlab='', xaxt='n', main='Overall influence', typ='n')
    axis(2, at=1:ninfl, labels=colnames(infl))
    cols <- apply(!is.na(infl), 1, sum)
    for(i in 1:nobs) abline(v=i, col=cols[i], lwd=(1+cols[i]/5))
    matplot(infl, pch=1, col=1, add=TRUE)
    put.xax(rep)
}

#' @name lprof.spict
#' @title Create profile likelihood
#' @details The "lprof" list must containg the following keys:
#' \itemize{
#'   \item{"pars"}{ A character vector of length equal 1 or 2 containing the name(s) of the parameters to calculate the profile likelihood for.}
#'   \item{"parrange"}{ A vector containing the parameter range(s) to profile over: parrange = c(min(par1), max(par1), min(par2), max(par2)).}
#' }
#' Optional:
#' \itemize{
#'   \item{"nogridpoints"}{ Number of grid points to evaluate the profile likelihood for each parameter. Default: 9. Note: with two parameters the calculation time increases quadratically when increasing the number of gridpoints.}
#' }
#' @param input A list containing observations and initial values for non profiled parameters (essentially an inp list) with the additional key "lprof" (see details for required keys). A valid result from fit.spict() containing an "inp" key with the described properties is also accepted.
#' @return Nothing.
#' @export
lprof.spict <- function(input){
    if('par.fixed' %in% names(input)){
        cat('Detected input as a SPiCT result, proceeding...\n')
        rep <- input
        inp <- rep$inp
        pl <- rep$pl
        repflag <- TRUE
    } else {
        cat('Assuming input is a SPiCT inp, checking for validity...\n')
        inp <- input
        inp <- check.inp(inp)
        repflag <- FALSE
    }
    lprof <- inp$lprof
    # Check lprof input
    if(!'pars' %in% names(lprof)) stop('"pars" key of lprof is unspecified!')
    np <- length(lprof$pars)
    if(!np %in% 1:2) stop('Length of pars vector is ', np, ' but should be 1 or 2!')
    if(!'parrange' %in% names(lprof)) stop('No "parrange" key specified in lprof!')
    npp <- length(lprof$parrange)
    if(npp != (2*np)) stop('Length of "parrange" is ', npp, ' but should be ', 2*np)
    if(lprof$parrange[2] <= lprof$parrange[1]) stop('Upper parameter bound is small than lower bound!')
    if(np==2) if(lprof$parrange[4] <= lprof$parrange[3]) stop('Upper parameter bound is small than lower bound!')
    prange <- matrix(lprof$parrange, np, 2, byrow=TRUE)
    if(!'nogridpoints' %in% names(lprof)) lprof$nogridpoints <- 9
    if(!'phases' %in% names(inp)) inp$phases <- list()
    # Determine good initial values for other parameters
    if(!repflag){
        inpp <- inp
        inpp$phases <- list()
        pp <- rep(0, np)
        for(i in 1:np){
            pp[i] <- mean(prange[i, ])
            inpp$ini[[inp$lprof$pars[i]]] <- pp[i]
            inpp$phases[[inp$lprof$pars[i]]] <- -1
        }
        repp <- fit.spict(inpp)
        pl <- repp$pl
    }
    parvals <- matrix(0, lprof$nogridpoints, np)
    for(i in 1:np) parvals[, i] <- seq(prange[i, 1], prange[i, 2], length=lprof$nogridpoints)
    if(np==1){
        pv <- parvals[, 1, drop=FALSE]
    } else {
        pv <- expand.grid(parvals[, 1], parvals[, 2])
    }
    ngrid <- dim(pv)[1]
    do.grid <- function(i){
        asd2 <- inp
        asd2$ini <- pl
        for(j in 1:np){
            asd2$ini[[lprof$pars[j]]] <- pv[i, j]
            asd2$phases[[lprof$pars[j]]] <- -1
        }
        rep2 <- fit.spict(asd2)
        #cat(i, '..')
        return(rep2$opt$objective)
    }
    single.do.grid <- function(i, ngrid){
        lv <- do.grid(i)
        cat(sprintf("Profiling: %.2f%%\r", i/ngrid * 100))
        return(lv)
    }
    multi.do.grid <- function(i, ngrid, progfile){
        lv <- do.grid(i)
        ## Send progress update
        if(class(progfile)[1]=='fifo') writeBin(1/ngrid, progfile)
        return(lv)
    }
    partry <- try(library(multicore))
    #partry <- try(library(parallel))
    cat('Profiling pars:', paste(lprof$pars, collapse=' and '), 'with', ngrid, 'trials.\n')
    if(class(partry)=='try-error'){
    #if(TRUE){
        likvals <- lapply(1:ngrid, single.do.grid, ngrid)
    } else {
        ## Open fifo progress file
        progfile <- fifo(tempfile(), open="w+b", blocking=T)
        if(inherits(fork(), "masterProcess")) {
        #if(inherits(parallel:::mcfork(), "masterProcess")) {
            ## Child
            progress <- 0.0
            while (progress < 1 && !isIncomplete(progfile)) {
                msg <- readBin(progfile, "double")
                progress <- progress + as.numeric(msg)
                cat(sprintf("Multicore profiling: %.2f%%\r", progress * 100))
            } 
            exit()
        }
        likvals <- mclapply(1:ngrid, multi.do.grid, ngrid, progfile, mc.cores=4)
        ## Close progress file
        close(progfile)
    }
    if(np==1){
        lprof$likvals <- unlist(likvals)
    } else {
        lprof$likvals <- matrix(unlist(likvals), lprof$nogridpoints, lprof$nogridpoints)
    }
    lprof$parvals <- parvals
    if(repflag){
        rep$inp$lprof <- lprof
        return(rep)
    } else {
        inp$lprof <- lprof
        return(inp)
    }
}


#' @name plotspict.lprof
#' @title Plots result of likelihood profiling.
#' @details TBA
#' @param input Result of running lprof.spict().
#' @param logpar If TRUE log of parameters are shown.
#' @return Nothing but shows a plot.
#' @export
plotspict.lprof <- function(input, logpar=FALSE){
    repflag <- 'par.fixed' %in% names(input)
    if(repflag){
        rep <- input
        inp <- rep$inp
        nll <- rep$opt$objective
    } else {
        inp <- input
        nll <- min(inp$lprof$likvals)
    }
    lprof <- inp$lprof
    np <- length(lprof$pars)
    pv <- lprof$parvals
    pars <- lprof$pars
    loginds <- grep('log', pars)
    expinds <- setdiff(1:np, loginds)
    if(!logpar){
        pv[, loginds] <- exp(pv[, loginds])
        pars[loginds] <- gsub('log', '', pars[loginds])
    } else {
        pv[, expinds] <- log(pv[, expinds])
        pars[expinds] <- paste0('log', pars[expinds])
    }
    if(np==1){
        plot(pv[, 1], lprof$likvals, typ='l', xlab=pars[1], ylab='negative log lik')
        lrlim <- 0.5*qchisq(0.95, 1) + nll
        abline(h=lrlim, lty=2)
    } else {
        pvals <- pchisq(2*(lprof$likvals - nll), np)
        contour(pv[, 1], pv[, 2], pvals, xlab=pars[1], ylab=pars[2], levels=c(0.5, 0.8, 0.95))
    }
}


#' @name invlogit
#' @title Inverse logit transform.
#' @param a Value to take inverse logit of.
#' @return Inverse logit.
invlogit <- function(a) 1/(1+exp(-a))


#' @name invlogp1
#' @title Inverse log "plus one" transform
#' @details If a = log(b-1), then the inverse transform is b = 1 + exp(a). Useful for values with lower bound at 1.
#' @param a Value to take inverse logp1 of.
#' @return Inverse logp1.
invlogp1 <- function(a) 1 + exp(a)
