#' @exportClass spictcls
setClass("spictcls")

#' @name test.spict
#' @title Example of a spict analysis.
#' @details Loads a data set, fits the model, calculates one-step-ahead residuals, plots the results.
#' @return A result report as given by fit.spict().
#' @examples
#' rep <- test.spict()
#' @export
test.spict <- function(){
    # Load data
    data(pol)
    # Fit model
    rep <- fit.spict(inp=pol$albacore)
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
#'  \item{"inp$timeC"}{ Vector of catch times.}
#'  \item{"inp$obsI"}{ List containing vectors of index observations.}
#'  \item{"inp$timeI"}{ List containing vectors of index times.}
#'  \item{"inp$logr"}{ Initial value for logr (log intrinsic growth rate).}
#'  \item{"inp$logK"}{ Initial value for logK (log carrying capacity).}
#'  \item{"inp$logq"}{ Initial value for logq (log catchability of index).}
#'  \item{"inp$logsdb"}{ Initial value for logsdb (log standard deviation of log biomass).}
#'  \item{"inp$logsdf"}{ Initial value for logsdf (log standard deviation of log fishing mortality).}
#' }
#' Optional inputs:
#' \itemize{
#' \item{"FILL IN"}{}
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
    if(!"lamperti" %in% names(inp)) inp$lamperti <- 1
    if(!"euler" %in% names(inp)) inp$euler <- 1
    if(!"dtc" %in% names(inp)){
        inp$dtc <- min(diff(inp$timeC))
        cat(paste('Catch interval (dtc) not specified. Assuming an interval of:', inp$dtc, 'year.\n'))
    }
    if(length(inp$dtc)==1) inp$dtc <- rep(inp$dtc, inp$nobsC)
    if(length(inp$dtc) != inp$nobsC) stop('Catch interval vector (inp$dtc) does not match catch observation vector (inp$obsC) in length')
    if(!"dteuler" %in% names(inp)) inp$dteuler <- min(inp$dtc)
    if(!"delay" %in% names(inp)) inp$delay <- 1
    if(!"dtpred" %in% names(inp)) inp$dtpred <- min(inp$dtc)
    if(!"RE" %in% names(inp)) inp$RE <- c('logF', 'logB')
    if(!"scriptname" %in% names(inp)) inp$scriptname <- 'spict'

    # -- DERIVED VARIABLES --
    alltimes <- inp$timeC
    for(i in 1:inp$nindex) alltimes <- c(alltimes, inp$timeI[[i]])
    #for(i in 1:1) alltimes <- c(alltimes, inp$timeI[[i]]) # CHANGE THIS WHEN USING ALL INDICES
    inp$timerange <- range(alltimes)
    # Add two dtc intervals, one for the final catch observation, and one for the predicted catch
    #mindtc <- min(inp$dtc)
    lastdtc <- tail(inp$dtc,1)
    timepad <- lastdtc + inp$dtpred
    inp$time <- seq(inp$timerange[1], inp$timerange[2]+timepad, by=inp$dteuler)
    inp$ns <- length(inp$time)
    inp$indlastobs <- which(inp$time == inp$timerange[2])
    # ic is the indices of inp$time to which catch observations correspond
    inp$ic <- match(inp$timeC, inp$time)
    # nc is number of states to integrate a catch observation over
    inp$nc <- rep(0, inp$nobsC)
    for(i in 1:inp$nobsC) inp$nc[i] <- sum(inp$time >= inp$timeC[i] & inp$time < (inp$timeC[i]+inp$dtc[i]))
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
    inp$dtpredinds <- which(inp$time >= (inp$timerange[2]+lastdtc) & inp$time < (inp$timerange[2]+timepad))
    inp$dtprednsteps <- length(inp$dtpredinds)
    #if(tail(inp$ic,1) == inp$ns) inp$dt <- c(inp$dt, tail(inp$dtc,1))
    
    # -- MODEL PARAMETERS --
    # Check that required initial values are specified
    if(!'logr' %in% names(inp$ini)){
        stop('Please specify initial value(s) for logr!')
    } else {
        nr <- length(inp$ini$logr)
        if(!nr %in% c(1, 2, 4)){
            stop('logr must be a vector of length either 1, 2 or 4.')
        }
        min <- log(0.3)-3
        max <- log(0.3)+3
        for(i in 1:nr) if(inp$ini$logr[i] < min | inp$ini$logr[i] > max) stop('Please specify a value for logr(', i, ') in the interval: [', min, ';', max, ']')
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
    if(!"logalpha" %in% names(inp$ini)) inp$ini$logalpha <- log(1)
    if(!"logbeta" %in% names(inp$ini)) inp$ini$logbeta <- log(1)
    #if(!"loggamma" %in% names(inp$ini)) inp$ini$loggamma <- log(1)
    if(!"logbkfrac" %in% names(inp$ini)) inp$ini$logbkfrac <- log(0.3)
    if(!"logF0" %in% names(inp$ini)) inp$ini$logF0 <- log(0.2*exp(inp$ini$logr[1]))
    if(!"logF" %in% names(inp$ini)){
        inp$ini$logF <- rep(inp$ini$logF0, inp$ns)
    } else {
        if(length(inp$ini$logF) != inp$ns) stop('Wrong length of inp$ini$logF: ', length(inp$ini$logF), ' Should be equal to inp$ns: ', inp$ns, ' remember logF is predicted beyond observations.')
    }
    if(!"logB" %in% names(inp$ini)){
        inp$ini$logB <- log(rep(exp(inp$ini$logbkfrac)*exp(inp$ini$logK), inp$ns))
    } else {
        if(length(inp$ini$logB) != inp$ns) stop('Wrong length of inp$ini$logB: ', length(inp$ini$logB), ' Should be equal to inp$ns: ', inp$ns, ' remember logB is predicted beyond observations.')
    }

    # Reorder parameter list
    inp$ini <- list(phi1=inp$ini$phi1,
                    phi2=inp$ini$phi2,
                    logalpha=inp$ini$logalpha,
                    logbeta=inp$ini$logbeta,
                    #loggamma=inp$ini$loggamma,
                    logbkfrac=inp$ini$logbkfrac,
                    logF0=inp$ini$logF0,
                    logr=inp$ini$logr,
                    logK=inp$ini$logK,
                    logq=inp$ini$logq,
                    logsdf=inp$ini$logsdf,
                    logsdb=inp$ini$logsdb,
                    logF=inp$ini$logF,
                    logB=inp$ini$logB)
    # Determine phases and fixed parameters
    if(inp$delay==1){
        fixpars <- c('phi1', 'phi2', 'logalpha', 'logbeta') # These are fixed unless specified
    } else {
        fixpars <- c('logalpha', 'logbeta') # These are fixed unless otherwise specified
    }
    if(!"phases" %in% names(inp)){
        inp$phases <- list()
        for(i in 1:length(fixpars)) inp$phases[[fixpars[i]]] <- -1
    } else {
        for(i in 1:length(fixpars)) if(!fixpars[i] %in% names(inp$phases)) inp$phases[[fixpars[i]]] <- -1
    }
    nphasepars <- length(inp$phases)
    inp$map <- list()
    for(i in 1:nphasepars){
        parnam <- names(inp$phases)[i]
        if(parnam %in% names(inp$ini)){
            phase <- inp$phases[[parnam]]
            if(phase<0){
                inp$map[[parnam]] <- factor(rep(NA, length(inp$ini[[parnam]])))
            } else {
                cat('WARNING: Phases not yet implemented! will estimate everything in one phase.\n')
            }
        } else {
            cat(paste('WARNING: Phase specified for an invalid parameter:', parnam, '\n'))
        }
    }

    inp$checked <- TRUE
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
    if(!is.null(min) & !is.null(max)){
        if(inp$ini[[parname]] < min | inp$ini[[parname]] > max) stop('Please specify a value for ', parname, ' in the interval: [', min, ';', max, ']')
    }
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
    if(!'checked' %in% names(inp)) inp <- check.inp(inp)
    if(!inp$checked) inp <- check.inp(inp)
    # Currently only able to use one index.
    datin <- list(delay=inp$delay, dt=inp$dt, dtpred=inp$dtpred, dtpredinds=inp$dtpredinds, dtprednsteps=inp$dtprednsteps, obsC=inp$obsC, ic=inp$ic, nc=inp$nc, I=inp$obsIin, ii=inp$iiin, iq=inp$iqin, ir=inp$ir, lamperti=inp$lamperti, euler=inp$euler, dbg=dbg)
    obj <- TMB::MakeADFun(data=datin, parameters=inp$ini, random=inp$RE, DLL=inp$scriptname, hessian=TRUE, map=inp$map)
    config(trace.optimize=0, DLL=inp$scriptname)
    verbose <- FALSE
    obj$env$tracemgc <- verbose
    obj$env$inner.control$trace <- verbose
    obj$env$silent <- ! verbose
    obj$fn(obj$par)
    rep <- NULL
    if(dbg==0){
        opt <- nlminb(obj$par, obj$fn, obj$gr)
        # Results
        rep <- TMB::sdreport(obj)
        rep$pl <- obj$env$parList(opt$par)
        obj$fn()
        rep$Cp <- obj$report()$Cp
        rep$inp <- inp
        rep$opt <- opt
    }
    class(rep) <- "spictcls"
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
#' rep <- calc.osa.resid(rep)
#' plotspict.osar(rep)
#' @import TMB
calc.osa.resid <- function(rep, dbg=0){
    inp <- rep$inp
    if("logF" %in% names(inp$ini)) inp$ini$logF <- NULL
    if("logB" %in% names(inp$ini)) inp$ini$logB <- NULL
    if("ns" %in% names(inp)) inp$ns <- NULL
    predmap <- get.predmap(rep$pl, inp$RE)
    plnew <- rep$pl
    logCpred <- rep(0, inp$nobsC-1)
    for(nadj in inp$delay:(inp$nobsC-1)){
        inp2 <- inp
        inp2$timeC <- inp$timeC[1:nadj]
        inp2$obsC <- inp$obsC[1:nadj]
        inp2$dtc <- inp$dtc[1:nadj]
        endtime <- inp$timeC[nadj]
        for(i in 1:inp$nindex){
            Iind <- which(inp$timeI[[i]] < endtime)
            inp2$timeI[[i]] <- inp$timeI[[i]][Iind] 
            inp2$obsI[[i]] <- inp$obsI[[i]][Iind]
        }
        inp2 <- check.inp(inp2)
        datnew <- list(delay=inp2$delay, dt=inp2$dt, dtpred=inp2$dtpred, dtpredinds=inp2$dtpredinds, dtprednsteps=inp2$dtprednsteps, obsC=inp2$obsC, ic=inp2$ic, nc=inp2$nc, I=inp2$obsIin, ii=inp2$iiin, iq=inp2$iqin, ir=inp$ir, lamperti=inp2$lamperti, euler=inp2$euler, dbg=dbg)
        for(k in 1:length(inp2$RE)) plnew[[inp2$RE[k]]] <- rep$pl[[inp2$RE[k]]][1:inp2$ns]
        objpred <- TMB::MakeADFun(data=datnew, parameters=plnew, map=predmap, random=inp2$RE, DLL=inp2$scriptname, hessian=TRUE, tracemgc=FALSE)
        verbose <- FALSE
        objpred$env$tracemgc <- verbose
        objpred$env$inner.control$trace <- verbose
        objpred$env$silent <- ! verbose
        objpred$fn()
        logCpred[nadj] <- log(objpred$report()$Cp)
    }
    if(inp$delay>1) logCpred <- logCpred[-(1:(inp$delay-1))]
    timeC <- inp$timeC[-(1:inp$delay)]
    logCpres <- log(inp$obsC[-(1:inp$delay)]) - logCpred
    # Test for independence of residuals (one-step-ahead predictions versus observations)
    logCpboxtest <- Box.test(logCpres, lag=1, type='Ljung-Box')
    rep$osar <- list(logCpres=logCpres, logCpred=logCpred, timeC=timeC, logCpboxtest=logCpboxtest)
    return(rep)
}


#' @name get.par
#' @title Extract parameters from a result report as generated by fit.spict.
#' @details Helper function for extracting the value and uncertainty of a specific model parameter, random effect or derived quantity.
#' @param parname Character string containing the name of the variable of interest.
#' @param rep A result report as generated by running fit.spict.
#' @param exp Take exp of the variable? TRUE/FALSE.
#' @param random Is the variable a random effect? TRUE/FALSE.
#' @param fixed Is the variable a fixed effect? TRUE/FALSE.
#' @return A matrix with four columns containing respectively: 1) the lower 95% confidence limit; 2) the parameter estimate; 3) the upper 95% confidence limit; 4) the parameter standard deviation in the domain it was estimated (log or non-log).
#' @export
#' @examples
#' data(pol)
#' rep <- fit.spict(pol$albacore)
#' Bmsy <- get.par('logBmsy', rep, exp=TRUE)
#' Best <- get.par('logB', rep, exp=TRUE, random=TRUE)
#' K <- get.par('logK', rep, exp=TRUE, fixed=TRUE)
get.par <- function(parname, rep=rep, exp=FALSE, random=FALSE, fixed=FALSE){
    if(random){
        ind <- which(names(rep$par.random)==parname)
        est <- rep$par.random[ind]
        sd <- sqrt(rep$diag.cov.random[ind])
        ll <- est - 1.96*sd
        ul <- est + 1.96*sd
    } else {
        if(fixed){
            ind <- which(names(rep$par.fixed)==parname)
            est <- rep$par.fixed[ind]
            sd <- sqrt(diag(rep$cov.fixed))[ind]
            ll <- est - 1.96*sd
            ul <- est + 1.96*sd
        } else {
            ind <- which(names(rep$value)==parname)
            est <- rep$value[ind]
            sd <- rep$sd[ind]
            ll <- est - 1.96*sd
            ul <- est + 1.96*sd
        }
    }
    if(length(est)==0){
        cat(paste('WARNING: could not extract', parname, '\n'))
    } else {
        if(exp){
            cbind(ll=exp(ll), est=exp(est), ul=exp(ul), sd)
        } else {
            cbind(ll, est, ul, sd)
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


#' @name tc.fun2
#' @title Calculate time constant (proportion of B_infinity reached).
#' @details Calculate the time required to reach a certain proportion of the equilibrium biomass (B_infinity) from a given starting point.
#' @param F Fishing mortality.
#' @param K Carrying capacity.
#' @param r Intrinsic growth rate.
#' @param sdb Standard deviation of the log biomass.
#' @param B0 Starting biomass.
#' @param p Proportion of B_infinity.
#' @param lamperti Return the stochastic (TRUE) or deterministic (FALSE) estimate.
#' @return The calculated time to reach the given proportion of B_infinity.
tc.fun2 <- function(F, K, r, sdb, B0, p, lamperti){
    Binf <- calc.binf(K, F, r, sdb, lamperti)
    if(lamperti){
        rate <- F + 0.5*sdb^2 - r
    } else {
        rate <- F - r
    }
    Brat <- Binf/B0
    Brat[Brat<1] <- 1/Brat[Brat<1] # Invert if B0 > Binf
    return( 1/(rate) * log( (1-p) / (p*(Brat-1))) )
}

#' @name annual
#' @title Convert from quarterly (or other sub-annual) data to annual means.
#' @param inp An inp list.
#' @param vec The vector of values to convert to annual means
#' @return A list containing the annual means and a corresponding time vector.
annual <- function(inp, vec){
    anntime <- inp$time[which(inp$time %% 1 ==0)]
    nanntime <- length(anntime)
    annvec <- rep(0, nanntime)
    floortime <- floor(inp$time)
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
#' @export
plotspict.biomass <- function(rep, logax=FALSE){
    log <- ifelse(logax, 'y', '')
    inp <- rep$inp
    # Biomass plot
    Best <- get.par('logB', rep, exp=TRUE, random=TRUE)
    ns <- dim(Best)[1]
    Kest <- get.par('logK', rep, exp=TRUE, fixed=TRUE)
    Bmsy <- get.par('logBmsy', rep, exp=TRUE)
    Binf <- get.par('logBinf', rep, exp=TRUE)
    qest <- get.par('logq', rep, fixed=TRUE, exp=TRUE)
    inds <- which(is.na(Binf) | Binf<0)
    Binf[inds] <- 1e-12
    annlist <- annual(inp, Binf[, 2])
    Binftime <- annlist$anntime
    Binfs <- annlist$annvec
    Bp <- get.par('logBp', rep, exp=TRUE)
    scal <- 1
    cicol <- 'lightgray'
    par(mar=c(5,4,4,4))
    ylim <- range(Best[, 1:3], Bp[1:3])/scal
    plot(inp$time, Best[,2]/scal, typ='n', xlab='Time', ylab=paste('Biomass'), main=paste('- Bmsy:',round(Bmsy[2]),' K:',round(Kest[2])), ylim=ylim, xlim=range(c(inp$time, tail(inp$time,1)+1)), log=log)
    axis(4, labels=pretty(ylim/Bmsy[2]), at=pretty(ylim/Bmsy[2])*Bmsy[2])
    mtext("B/Bmsy", side=4, las=0, line=2)
    polygon(c(inp$time[1]-5,tail(inp$time,1)+5,tail(inp$time,1)+5,inp$time[1]-5), c(Bmsy[1],Bmsy[1],Bmsy[3],Bmsy[3]), col=cicol, border=cicol)
    for(i in 1:inp$nindex) points(inp$timeI[[i]], inp$obsI[[i]]/qest[i, 2], pch=i, cex=0.7)
    if('true' %in% names(inp)){
        lines(inp$time, inp$true$B/scal, col='orange') # Plot true
        abline(h=inp$true$Bmsy, col='orange', lty=2)
    }
    lines(inp$time, Best[,2]/scal, col='blue')
    abline(h=Bmsy[2]/scal, col='black')
    lines(inp$time, Best[,1]/scal, col=4, lty=2)
    lines(inp$time, Best[,3]/scal, col=4, lty=2)
    lines(Binftime, Binfs/scal, col='green', lty=1)
    et <- tail(inp$time,1)
    tp <- et + inp$dtpred
    lines(c(et, tp), c(tail(Best[,2],1), Bp[2])/scal, col='blue', lty=3)
    points(tp, Bp[2], pch=21, bg='yellow')
    lines(c(et, tp), c(tail(Best[,1],1), Bp[1])/scal, col='blue', lty=3)
    lines(c(et, tp), c(tail(Best[,3],1), Bp[3])/scal, col='blue', lty=3)
    legend('topright', legend=c('Equilibrium',paste(tp,'pred.')), lty=c(1,NA), pch=c(NA,21), col=c('green',1), pt.bg=c(NA,'yellow'), bg='white')
    box()
}


#' @name plotspict.osar
#' @title Plot one-step-ahead residuals
#' @details Plots observed versus predicted catches.
#' @param rep A result report as generated by running fit.spict.
#' @return Nothing.
#' @export
plotspict.osar <- function(rep){
    inp <- rep$inp
    Cscal <- 1
    Cpred <- rep$osar$logCpred
    plot(inp$timeC, log(inp$obsC), typ='p', ylim=range(c(Cpred,log(inp$obsC),1.13*c(Cpred,log(inp$obsC)))), main=paste('Ljung-Box test p-value:',round(rep$osar$logCpboxtest$p.value,5)), ylab=paste('log Catch, Cscal:',Cscal), xlim=range(c(inp$timeC,inp$timeC[inp$nobsC]+1)), xlab='Time', col=1)
    clr <- 'blue'
    lines(rep$osar$timeC, Cpred, col=clr)
    points(rep$osar$timeC, Cpred, pch=20, cex=0.7, col=clr)
    legend('topright', c('Observations', 'OSA pred.'), lty=c(NA,1), col=c(1,clr), pch=c(1,20), pt.cex=c(1,0.7))
}


#' @name plotspict.fb
#' @title Plot fishing mortality versus biomass.
#' @details Plots estimated fishing mortality as a function of biomass together with reference points and the prediction for next year given a constant F. The equilibrium biomass for F fixed to the current value is also plotted.
#' @param rep A result report as generated by running fit.spict.
#' @param logax Take log of x and y-axes? default: FALSE
#' @return Nothing.
#' @export
plotspict.fb <- function(rep, logax=FALSE){
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
    Fest <- get.par('logF', rep, exp=TRUE, random=TRUE)
    Fmsy <- get.par('logFmsy', rep, exp=TRUE)
    Fp <- Fest[ns,]
    inds <- c(which(names(rep$value)=='logBmsy'), which(names(rep$value)=='logFmsy'))
    cl <- make.ellipse(inds, rep)
    xlim <- range(c(exp(cl[,1]),Bp[1:3],Best[,2])/scal)
    ylim <- range(c(exp(cl[,2]),Fest[,2]))
    main <- paste('logalpha:',rep$pl$logalpha,'r:',round(rest[2],3),'q:',round(qest[2],3))
    par(mar=c(5,4,4,4))
    plot(Bmsy[2]/scal, Fmsy[2], typ='p', xlim=xlim, xlab='Biomass', ylab='F',  pch=24, bg='blue', ylim=ylim, log=log)
    axis(3, labels=pretty(xlim/Bmsy[2]), at=pretty(xlim/Bmsy[2])*Bmsy[2])
    mtext("B/Bmsy", side=3, las=0, line=2)
    axis(4, labels=pretty(ylim/Fmsy[2]), at=pretty(ylim/Fmsy[2])*Fmsy[2])
    mtext("F/Fmsy", side=4, las=0, line=2)
    polygon(exp(cl[,1])/scal, exp(cl[,2]), col=cicol, border=cicol)
    abline(h=Fmsy[2], lty=3)
    abline(v=Bmsy[2], lty=3)
    #arrow.line(Best[,2]/scal, Fest[,2], length=0.05, col='blue')
    if('true' %in% names(inp)){
        lines(inp$true$B/scal, inp$true$F, col='orange') # Plot true
        points(inp$true$Bmsy, inp$true$Fmsy, pch=24, bg='orange')
    }
    lines(Best[,2]/scal, Fest[,2], col='blue')
    points(Bmsy[2]/scal, Fmsy[2], pch=24, bg='black')
    lines(c(tail(Best[,2],1), Bp[2])/scal, rep(Fp[2],2), col='blue', lty=3)
    points(Bp[2]/scal, Fp[2], pch=21, bg='yellow')
    points(tail(Binf[,2],1)/scal, Fp[2], pch=22, bg='green', cex=2)
    arrow.line(c(Bp[2], tail(Binf[,2],1))/scal, rep(Fp[2],2), col='black', length=0.05)
    legend('topright', c('Estimated MSY',paste(tail(inp$time,1)+1,'prediction'),'Equilibrium'), pch=c(24,21,22), pt.bg=c('black','yellow','green'), bg='white')
}


#' @name plotspict.f
#' @title Plot estimated fishing mortality.
#' @details Plots estimated fishing mortality with Fmsy and associated confidence interval.
#' @param rep A result report as generated by running fit.spict.
#' @param logax Take log of y-axis? default: FALSE
#' @return Nothing.
#' @export
plotspict.f <- function(rep, logax=FALSE){
    log <- ifelse(logax, 'y', '')
    inp <- rep$inp
    cicol <- 'lightgray'
    Fest <- get.par('logF', rep, exp=TRUE, random=TRUE)
    Fmsy <- get.par('logFmsy', rep, exp=TRUE)
    rest <- get.par('logr', rep, exp=TRUE, fixed=TRUE)
    ylim <- range(Fest[,1:3])
    plot(inp$time, Fest[, 2], typ='n', main=paste('Fmsy:',round(Fmsy[2],3)), ylim=ylim, col='blue', ylab='F', xlim=range(inp$time,tail(inp$time+1,2)), xlab='Time', log=log)
    axis(4, labels=pretty(ylim/Fmsy[2]), at=pretty(ylim/Fmsy[2])*Fmsy[2])
    mtext("F/Fmsy", side=4, las=0, line=2)
    polygon(c(inp$time[1]-5,tail(inp$time,1)+5,tail(inp$time,1)+5,inp$time[1]-5), c(Fmsy[1],Fmsy[1],Fmsy[3],Fmsy[3]), col=cicol, border=cicol)
    #abline(h=Fmax/Fmsy, col='red')
    #lines(inp$time, Fest/Fmsy, col='black')
    if('true' %in% names(inp)){
        lines(inp$time, inp$true$F, col='orange') # Plot true
        abline(h=inp$true$Fmsy, col='orange', lty=2)
    }
    maincol <- 'blue'
    if(min(inp$dtc) < 1){
        al1 <- annual(inp, Fest[, 1])
        al2 <- annual(inp, Fest[, 2])
        al3 <- annual(inp, Fest[, 3])
        lines(al1$anntime, al1$annvec, col=maincol, lwd=2, lty=2)
        lines(al2$anntime, al2$annvec, col=maincol, lwd=2)
        lines(al3$anntime, al3$annvec, col=maincol, lwd=2, lty=2)
        maincol <- 'cyan'
    }
    lines(inp$time, Fest[, 2], col=maincol)
    #points(inp$time, C/Best/Fmsy)
    abline(h=Fmsy[2], col='black')
    abline(h=rest[2], col='red')
    #abline(h=Fmsyci, col=3, lty=2)
    #lines(inp$time, Festul/Fmsy, col='forestgreen', lty=2)
    #lines(inp$time, Festll/Fmsy, col='forestgreen', lty=2)
    lines(inp$time, Fest[, 1], col=maincol, lty=2)
    lines(inp$time, Fest[, 3], col=maincol, lty=2)
    box()
}


#' @name plotspict.catch
#' @title Plot observed catch and predictions.
#' @details Plots observed catch and predictions using the current F and Fmsy. The plot also contains the equilibrium catch if the current F is maintained.
#' @param rep A result report as generated by running fit.spict.
#' @return Nothing.
#' @export
plotspict.catch <- function(rep){
    inp <- rep$inp
    Cscal <- 1
    cicol <- 'lightgray'
    MSY <- get.par('MSY', rep, exp=FALSE)
    Cpmsy <- get.par('Cpmsy', rep)
    Cpmsy[Cpmsy<0] <- 0
    Cinfp <- get.par('Cinfp', rep)
    Cinfp[Cinfp<0] <- 0
    Cpredsub <- get.par('Cpredsub', rep)
    Pest <- get.par('P', rep)
    Cpredest <- get.par('logCpred', rep, exp=TRUE)
    Cpredest[Cpredest<0] <- 0
    rep$Cp[rep$Cp<0] <- 0
    plot(inp$timeC, inp$obsC/Cscal, typ='n', main=paste('MSY:',round(MSY[2]/Cscal)), xlab='Time', ylab=paste('Catch'), xlim=range(c(inp$time, tail(inp$time,1)+1)), ylim=range(c(1.3*inp$obsC, Cpredest[,1:3], 0.8*inp$obsC, Cinfp[2], Cpmsy[2], rep$Cp))/Cscal)
    polygon(c(inp$time[1]-5,tail(inp$time,1)+5,tail(inp$time,1)+5,inp$time[1]-5), c(MSY[1],MSY[1],MSY[3],MSY[3])/Cscal, col=cicol, border=cicol)
    points(inp$timeC, inp$obsC/Cscal)
    points(tail(inp$timeC,1)+inp$dtpred, rep$Cp/Cscal, pch=21, bg='yellow')
    #points(tail(inp$timeC,1)+1, Cpmsy[2]/Cscal, pch=21, bg='black')
    #points(tail(inp$timeC,1)+1, Cinfp[2]/Cscal, pch=21, bg='green')
    if('true' %in% names(inp)) abline(h=inp$true$MSY, col='orange', lty=2)
    abline(h=MSY[2]/Cscal)
    lines(inp$timeC, Cpredest[, 2]/Cscal, col=4)
    lines(inp$timeC, Cpredest[, 1]/Cscal, col=4, lty=2)
    lines(inp$timeC, Cpredest[, 3]/Cscal, col=4, lty=2)
    #legend('topleft',c(paste(tail(inp$time,1)+inp$dtpred,'prediction'),'Pred w Fmsy','Equilibrium'),pch=21,pt.bg=c('yellow','black','green'), bg='white')
    legend('topleft',c(paste(tail(inp$time,1)+inp$dtpred,'prediction')), pch=21, pt.bg=c('yellow'), bg='white')
    box()
}


#' @name plotspict.production
#' @title Plot theoretical production curve and estimates.
#' @details Plots the theoretical production curve (production as a function of biomass) as calculated from the estimated model parameters. Overlaid is the estimated production/biomass trajectory.
#' @param rep A result report as generated by running fit.spict.
#' @return Nothing.
#' @export
plotspict.production <- function(rep){
    inp <- rep$inp
    Best <- get.par('logB', rep, exp=TRUE, random=TRUE)
    Kest <- get.par('logK', rep, exp=TRUE, fixed=TRUE)
    rest <- get.par('logr', rep, exp=TRUE, fixed=TRUE)
    Bmsy <- get.par('logBmsy', rep, exp=TRUE)
    Pest <- get.par('P', rep)
    nBplot <- 200
    Bplot <- seq(0.5*min(Best[, 2]), 1*max(Best[, 2]), length=nBplot)
    pfun <- function(r, K, B) r*B*(1 - B/K)
    Pst <- pfun(rest[2], Kest[2], Bplot)
    xlim <- range(Bplot/Bmsy[2])
    plot(Best[-1, 2]/Bmsy[2], Pest[, 2]/inp$dt/Bmsy[2], typ='l', ylim=range(Pest[,2]/inp$dt/Bmsy[2], Pst/Bmsy[2]), xlim=xlim, xlab='B/Bmsy', ylab='Production/Bmsy', col=4, main='Production curve')
    lines(Bplot/Bmsy[2], Pst/Bmsy[2], col=1)
    #arrow.line(Best[-1, 2]/Bmsy[2], Pest[,2]/inp$dt/Bmsy[2], length=0.05)
    abline(v=1, lty=3)
}


#' @name plotspict.tc
#' @title Plot time constant.
#' @details Plots the time required for the biomass to reach a certain proportion of Bmsy. The time required to reach 95% of Bmsy is highlighted.
#' @param rep A result report as generated by running fit.spict.
#' @return Nothing.
#' @export
plotspict.tc <- function(rep){
    inp <- rep$inp
    Best <- get.par('logB', rep, exp=TRUE, random=TRUE)
    Fest <- get.par('logF', rep, exp=TRUE, random=TRUE)
    Kest <- get.par('logK', rep, exp=TRUE, fixed=TRUE)
    rest <- get.par('logr', rep, exp=TRUE, fixed=TRUE)
    sdbest <- get.par('logsdb', rep, exp=TRUE, fixed=TRUE)
    Fmsy <- get.par('logFmsy', rep, exp=TRUE)
    Bmsy <- get.par('logBmsy', rep, exp=TRUE)
    B0cur <- Best[inp$indlastobs, 2]
    if(B0cur < Bmsy[2]) facvec <- c(0, 0.75, 0.95, 1)
    if(B0cur > Bmsy[2]) facvec <- c(2, 1.25, 1.05, 1)
    Fvec <- round(facvec*Fmsy[2], digits=4)
    nFvec <- length(Fvec)
    g <- function(F, K, r, sdb, B0, dt, lamperti){
        if(lamperti){
            return(exp(log(B0) + (r - r/K*B0 - F - 0.5*sdb^2)*dt))
        } else {
            return(B0 + B0*r*(1 - B0/K - F)*dt)
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
            Bsim[i, j] <- g(Fvec[i], Kest[2], rest[2], sdbest[2], Bsim[i, j-1], simdt, inp$lamperti)
        }
    }
    Bsim <- Bsim/Bmsy[2]
    frac <- 0.95
    if(B0cur < Bmsy[2]) inds <- which(Bsim[nFvec, ]<0.99)
    if(B0cur > Bmsy[2]) inds <- which(Bsim[nFvec, ]>(1/0.99))
    plot(time[1, ], Bsim[1, ], typ='l', xlim=range(time[nFvec, inds]), ylim=range(Bsim[nFvec, ]), col=3, ylab='Proportion of Bmsy', xlab='Years to Bmsy')
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


#' @name plot.spictcls
#' @title 3x3 plot illustrating spict results.
#' @details Create a 3x3 plot containing the following:
#' \itemize{
#'  \item{1. Biomass.}
#'  \item{2. One-step-ahead residuals (only if calculated).}
#'  \item{3. One-step-ahead auto-correlation function (only if calculated).}
#'  \item{4. Estimated F versus estimated B.}
#'  \item{5. Estimated fishing mortality.}
#'  \item{6. Observed versus predicted catches.}
#'  \item{7. Observed versus theoretical production.}
#'  \item{8. Calculated time-constant.}
#' }
#' @param rep A result report as generated by running fit.spict.
#' @param logax Take log of relevant axes? default: FALSE
#' @return Nothing.
#' @export
plot.spictcls <- function(rep, logax=FALSE){
    inp <- rep$inp
    #dev.new(width=10, height=10)
    if('osar' %in% names(rep)){
        par(mfrow=c(3, 3))
    } else {
        par(mfrow=c(2, 3))
    }
    # Biomass plot
    plotspict.biomass(rep, logax=logax)
    if('osar' %in% names(rep)){
        # One-step-ahead catch residuals
        plotspict.osar(rep)
        # OSAR ACF
        acf(rep$osar$logCpres, main='ACF of OSA residuals')
    }
    # F versus B
    plotspict.fb(rep, logax=logax)
    # F
    plotspict.f(rep, logax=logax)
    # Catch
    plotspict.catch(rep)
    # Production curve
    plotspict.production(rep)
    # Time constant
    plotspict.tc(rep)
}


#' @name summary.spictcls
#' @title Output a summary of a fit.spict() run.
#' @details The output includes, the convergence message from the optimiser, the likelihood value of the parameters, the parameter estimates with 95% confidence intervals, estimates of derived parameters (BMsy, Fmsy, MSY) with 95% confidence intervals, and predictions of biomass, fishing mortality, and catch for the value of inp$dtpred.
#' @param object A result report as generated by running fit.spict.
#' @param numdigits Present values with this number of digits after the dot.
#' @return Nothing.
#' @export
summary.spictcls <- function(object, numdigits=4){
    rep <- object
    cat(paste('Convergence: ', rep$opt$convergence, '  MSG: ', rep$opt$message, '\n', sep=''))
    if(rep$opt$convergence>0) cat('WARNING: Model did not obtain proper convergence! Estimates and uncertainties are most likely invalid and should not be used.\n')
    cat(paste('Negative log likelihood: ', round(rep$opt$objective, numdigits), '\n', sep=''))
    cat('\nModel parameter estimates \n')
    sd <- sqrt(diag(rep$cov.fixed))
    nms <- names(rep$par.fixed)
    loginds <- grep('log', nms)
    est <- rep$par.fixed
    est[loginds] <- exp(est[loginds])
    cilow <- rep$par.fixed-1.96*sd
    cilow[loginds] <- exp(cilow[loginds])
    ciupp <- rep$par.fixed+1.96*sd
    ciupp[loginds] <- exp(ciupp[loginds])
    if('true' %in% names(rep$inp)){
        npar <- length(nms)
        truepar <- rep(0, npar)
        for(i in 1:npar) truepar[i] <- rep$inp$true[[nms[i]]]
        truepar[loginds] <- exp(truepar[loginds])
        ci <- rep(0, npar)
        for(i in 1:npar) ci[i] <- as.numeric(truepar[i] > cilow[i] & truepar[i] < ciupp[i])
        resout <- cbind(estimate=round(est,numdigits), true=round(truepar,numdigits), cilow=round(cilow,numdigits), ciupp=round(ciupp,numdigits), true.in.ci=ci, est.in.log=round(rep$par.fixed,numdigits))
    } else {
        resout <- cbind(estimate=round(est,numdigits), cilow=round(cilow,numdigits), ciupp=round(ciupp,numdigits), est.in.log=round(rep$par.fixed,numdigits))
    }
    nms[loginds] <- substr(names(rep$par.fixed[loginds]),4,60)
    rownames(resout) <- nms
    cat(paste(capture.output(resout),' \n'))
    cat('\nDerived estimates \n')
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
    cat(paste(capture.output(derout),' \n'))
    cat('\nPredictions \n')
    predout <- rbind(
        get.par(parname='logBp2', rep, exp=TRUE)[c(2,1,3,2)],
        get.par(parname='logFp2', rep, exp=TRUE)[c(2,1,3,2)],
        get.par(parname='logCp2', rep, exp=TRUE)[c(2,1,3,2)])
    predout[, 4] <- log(predout[, 4])
    predout <- round(predout, numdigits)
    colnames(predout) <- c('prediction', 'cilow', 'ciupp', 'est.in.log')
    rownames(predout) <- c('B', 'F', 'Catch')
    cat(paste(capture.output(predout),' \n'))
    if('osar' %in% names(rep)){
        cat('\nOne-step-ahead residuals \n')
        cat(paste('Ljung-box test for independence of lag 1, p-value:', round(rep$osar$logCpboxtest$p.value, numdigits), '\n'))
    }
}



#' @name read.aspic
#' @title Reads ASPIC input file.
#' @details Reads an input file following the ASPIC 7 format described in the ASPIC manual (found here http://www.mhprager.com/aspic.html).
#' @param filename Path of the ASPIC input file.
#' @return A list of input variables that can be used as input to fit.spict().
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
#' @param sdb Standard deviation of biomass process.
#' @param lamperti Optional logical, use lamperti transformation?
#' @return B_infinity.
calc.binf <- function(K, F, r, sdb=0, lamperti=FALSE){
    if(!lamperti) sdb <- 0
    rate <- calc.rate(r, F, sdb)
    return(K/r*rate)
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
#' @title Simulate data from Schaefer model
#' @details Simulates data using either manually specified parameters values or parameters estimated by fit.spict().
#'
#' Manual specification:
#' To specify parameters manually use the inp$ini format similar to when specifying initial values for running fit.spict(). Observations can be simulated at specific times using inp$timeC and inp$timeI. If these are not specified then nobs observations will be simulated evenly distributed in time.
#'
#' Estimated parameters:
#' Simply take the output from a fit.spict() run and use as input to sim.spict().
#' 
#' @param input Either an inp list with an ini key (see ?check.inp) or a rep list where rep is the output of running fit.spict().
#' @param nobs Optional specification of the number of simulated observations.
#' @return A list containing the simulated data.
#' @export
sim.spict <- function(input, nobs=100){
    # Check if input is a inp (initial values) or rep (results).
    if('par.fixed' %in% names(input)){
        cat('Detected input as a SPiCT result, proceeding...\n')
        rep <- input
        inp <- rep$inp
        p <- rep$pl
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
            inp <- check.inp(inp)
            p <- inp$ini
            
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
    B0 <- exp(p$logbkfrac)*exp(p$logK)
    F0 <- exp(p$logF0)
    r <- exp(p$logr)
    K <- exp(p$logK)
    q <- exp(p$logq)
    sdb <- exp(p$logsdb)
    sdf <- exp(p$logsdf)
    alpha <- exp(p$logalpha)
    beta <- exp(p$logbeta)
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
        flag <- any(F >= 1.5*r) # Do this to avoid crazy F values
    }
    # - B_infinity
    Binf <- rep(0,nt)
    for(t in 1:nt) Binf[t] <- calc.binf(K, F[t], r, sdb, lamperti)
    # - Biomass -
    B <- rep(0,nt)
    B[1] <- B0
    e <- exp(rnorm(nt-1, 0, sdb*sqrt(dt)))
    for(t in 2:nt) B[t] <- predict.b(B[t-1], Binf[t-1], F[t-1], r, K, dt, sdb, lamperti, euler) * e[t-1]
    # - Catch -
    Csub <- rep(0,nt)
    for(t in 1:nt) Csub[t] <- predict.c(F[t], K, r, B[t], Binf[t], dt, sdb, lamperti, euler)
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
    # - Index observations -
    obsI <- list()
    for(I in 1:inp$nindex){
        obsI[[I]] <- rep(0, inp$nobsI[I])
        for(i in 1:inp$nobsI[I]) obsI[[I]][i] <- exp(log(q[I] * B[inp$ii[[I]][i]]) + rnorm(1, 0, sdi))
    }
    
    sim <- list()
    sim$obsC <- obsC
    sim$timeC <- inp$timeC
    sim$obsI <- obsI
    sim$timeI <- inp$timeI
    sim$ini <- list(logr=log(r), logK=log(K), logq=log(q), logsdf=log(sdf), logsdb=log(sdb))
    sim$dteuler <- inp$dteuler
    sim$euler <- inp$euler
    sim$lamperti <- inp$lamperti
    sim$phases <- inp$phases
    sim$true <- p
    sim$true$B <- B
    sim$true$F <- F
    sim$true$Bmsy <- K/2
    sim$true$Fmsy <- ifelse(inp$lamperti, r/2 - sdb^2/2, r/2)
    sim$true$MSY <- sim$true$Bmsy * sim$true$Fmsy
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
    Fmsytrue <- exp(inp$ini$logr)/2 - exp(inp$ini$logsdb)^2/2
    Bmsytrue <- exp(inp$ini$logK)/2
    matstan <- matrix(0, nsim, nnobsvec)
    colnames(matstan) <- nobsvec
    conv <- matrix(NA, nsim, nnobsvec)
    colnames(conv) <- nobsvec
    osarpvals <- matstan
    Fmsyci <- matstan
    Fmsyciw <- matstan
    Bmsyci <- matstan
    Bmsyciw <- matstan
    sim <- list()
    rep <- list()

    for(i in 1:nnobsvec){
        sim[[i]] <- list()
        rep[[i]] <- list()
        for(j in 1:nsim){
            sim[[i]][[j]] <- sim.spict(inp, nobs=nobsvec[i])
            if(!is.null(estinp)) sim[[i]][[j]]$ini <- estinp$ini
            fit <- try(fit.spict(sim[[i]][[j]]))
            if(class(fit)=='try-error'){
                rep[[i]][[j]] <- NA
            } else {
                rep[[i]][[j]] <- fit
                rep[[i]][[j]] <- calc.osa.resid(rep[[i]][[j]])
                conv[j, i] <- rep[[i]][[j]]$opt$convergence
                osarpvals[j, i] <- rep[[i]][[j]]$osar$logCpboxtest$p.value
                Fmsy <- get.par('logFmsy', rep[[i]][[j]], exp=TRUE)
                Fmsyci[j, i] <- Fmsytrue > Fmsy[1] & Fmsytrue < Fmsy[3]
                Fmsyciw[j, i] <- Fmsy[3] - Fmsy[1]
                Bmsy <- get.par('logBmsy', rep[[i]][[j]], exp=TRUE)
                Bmsyci[j, i] <- Bmsytrue > Bmsy[1] & Bmsytrue < Bmsy[3]
                Bmsyciw[j, i] <- Bmsy[3] - Bmsy[1]
            }
            cat(paste(Sys.time(), '- validating:  i:', i, 'j:', j, '\n'))
            if(!is.null(backup)) save(Fmsyci, Fmsyciw, Bmsyci, Bmsyciw, osarpvals, conv, file=backup)
        }
    }
    return(list(osarpvals=osarpvals, Fmsyci=Fmsyci, Fmsyciw=Fmsyciw, Bmsyci=Bmsyci, Bmsyciw=Bmsyciw, conv=conv))
}
