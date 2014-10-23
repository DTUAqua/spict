
#' @name test.spict
#' @title Example of a spict analysis.
#' @details Loads a data set, fits the model, calculates one-step-ahead residuals, plots the results.
#' @return Nothing.
#' @export
test.spict <- function(){
    # Load data
    data(pol.albacore)
    # Fit model
    rep <- fit.spict(inp)
    # Calculate one-step-ahead residuals
    rep <- calc.osa.resid(rep)
    # Plot results
    graphics.off()
    dev.new(width=10, height=10)
    plotspict(rep)
}


#' @name pol.albacore
#' @title Fisheries data for south Atlantic albacore tuna 1967-1989.
#' @details This data were included in Polacheck et al. (1993), Canadian Journal of Fisheries and Aquatic Science, vol 50, pp. 2597-2607.
#' @docType data
#' @keywords datasets
#' @usage data(pol.albacore)
#' @format An inp list containing the data and initial values for estimation.
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
#' @export
check.inp <- function(inp){
    # -- DATA --
    # Check catches
    if('obsC' %in% names(inp) & 'timeC' %in% names(inp)){
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
    if('obsI' %in% names(inp) & 'timeI' %in% names(inp)){
        if(class(inp$obsI)!='list'){
            tmp <- inp$obsI
            inp$obsI <- list()
            inp$obsI[[1]] <- tmp
        }
        inp$nindex <- length(inp$obsI)
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
        inp$dtc <- 1
        cat('Annual catch observations assumed\n')
    }
    if(length(inp$dtc)==1) inp$dtc <- rep(inp$dtc, inp$nobsC)
    if(length(inp$dtc) != inp$nobsC) stop('Catch interval vector (inp$dtc) does not match catch observation vector (inp$obsC) in length')
    if(!"timefrac" %in% names(inp)) inp$timefrac <- min(inp$dtc)
    if(!"map" %in% names(inp)) inp$map <- list(phi1=factor(NA), phi2=factor(NA), alpha=factor(NA), beta=factor(NA), loggamma=factor(NA))
    if(!"delay" %in% names(inp)) inp$delay <- 1
    if(!"dtpred" %in% names(inp)) inp$dtpred <- min(inp$dtc)
    if(!"RE" %in% names(inp)) inp$RE <- c('logF', 'logB')
    if(!"scriptname" %in% names(inp)) inp$scriptname <- 'spict'

    # -- DERIVED VARIABLES --
    alltimes <- inp$timeC
    for(i in 1:inp$nindex) alltimes <- c(alltimes, inp$timeI[[i]])
    timerange <- range(alltimes)
    inp$time <- seq(timerange[1], timerange[2]+1, by=inp$timefrac)
    inp$ns <- length(inp$time)
    # ic is the indices of inp$time to which catch observations correspond
    inp$ic <- match(inp$timeC, inp$time)
    # nc is number of states to integrate a catch observation over
    inp$nc <- rep(0, inp$nobsC)
    for(i in 1:inp$nobsC) inp$nc[i] <- sum(inp$time >= inp$timeC[i] & inp$time < (inp$timeC[i]+inp$dtc[i]))
    # ii is the indices of inp$time to which index observations correspond
    inp$ii <- list()
    for(i in 1:inp$nindex) inp$ii[[i]] <- match(inp$timeI[[i]], inp$time)
    inp$dt <- diff(inp$time)
    #if(tail(inp$ic,1) == inp$ns) inp$dt <- c(inp$dt, tail(inp$dtc,1))
    
    # -- MODEL PARAMETERS --
    # Check that required initial values are specified
    check.ini('logr', inp, min=log(0.3)-3, max=log(0.3)+1)
    check.ini('logK', inp)
    check.ini('logq', inp)
    check.ini('logsdb', inp, min=log(0.1), max=log(5))
    check.ini('logsdf', inp, min=log(0.1), max=log(5))
    # Fill in unspecified less important model parameter values
    if(!"phi1" %in% names(inp$ini)) inp$ini$phi1 <- 1
    if(!"phi2" %in% names(inp$ini)) inp$ini$phi2 <- 0
    if(!"alpha" %in% names(inp$ini)) inp$ini$alpha <- 1
    if(!"beta" %in% names(inp$ini)) inp$ini$beta <- 1
    if(!"loggamma" %in% names(inp$ini)) inp$ini$loggamma <- log(1)
    if(!"logF" %in% names(inp$ini)) inp$ini$logF <- log(rep(0.2*exp(inp$ini$logr), inp$ns))
    if(!"logB" %in% names(inp$ini)) inp$ini$logB <- log(rep(0.3*exp(inp$ini$logK), inp$ns))
    # Reorder parameter list
    inp$ini <- list(phi1=inp$ini$phi1,
                    phi2=inp$ini$phi2,
                    alpha=inp$ini$alpha,
                    beta=inp$ini$beta,
                    loggamma=inp$ini$loggamma,
                    logr=inp$ini$logr,
                    logK=inp$ini$logK,
                    logq=inp$ini$logq,
                    logsdf=inp$ini$logsdf,
                    logsdb=inp$ini$logsdb,
                    logF=inp$ini$logF,
                    logB=inp$ini$logB)

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


#' @name calc.osa.resid
#' @title Calculate one-step-ahead residuals.
#' @details In TMB one-step-ahead residuals are calculated by sequentially including one data point at a time while keeping the model parameters fixed at their ML estimates. The calculated residuals are tested for independence in lag 1 using the Ljung-Box test (see Box.test).
#' @param rep A result report as generated by running fit.spict.
#' @return An updated result report, which contains one-step-ahead residuals stored in the $osar variable.
#' @export
calc.osa.resid <- function(rep){
    inp <- rep$inp
    ns <- length(time)
    predmap <- get.predmap(rep$pl, inp$RE)
    plnew <- rep$pl
    logCpred <- rep(0, inp$nobsC-1)
    for(nadj in inp$delay:(inp$nobsC-1)){
        endtime <- inp$timeC[nadj]
        ptime <- inp$time[inp$time <= endtime]
        pns <- length(ptime)
        pdt <- diff(ptime)
        Cind <- which(inp$timeC < endtime)
        ptimeC <- inp$timeC[Cind]
        pCobs <- inp$obsC[Cind]
        pnCobs <- length(pCobs)
        Iind <- which(inp$timeI[[1]] < endtime) # Must change as timeI is a list
        ptimeI <- inp$timeI[[1]][Iind] # Must change as timeI is a list
        pI <- inp$obsI[[1]][Iind] # Must change as timeI is a list
        pic <- match(ptimeC, ptime)
        pdtc <- c(diff(inp$timeC), inp$dtpred)
        # nc is number of states to integrate a catch observation over
        pnc <- rep(0, pnCobs) # Initialise
        for(i in 1:pnCobs) pnc[i] <- sum(ptime >= ptimeC[i] & ptime < (ptimeC[i]+pdtc[i]))
        pii <- match(ptimeI, ptime)
        datnew <- list(delay=inp$delay, dt=pdt, dtpred=inp$dtpred, Cobs=pCobs, ic=pic, nc=pnc, I=pI, ii=pii, isum=rep(0,inp$ns), lamperti=inp$lamperti, euler=inp$euler, dbg=0)
        for(k in 1:length(inp$RE)) plnew[[inp$RE[k]]] <- rep$pl[[inp$RE[k]]][1:length(ptime)]
        objpred <- MakeADFun(data=datnew, parameters=plnew, map=predmap, random=inp$RE, DLL=inp$scriptname, hessian=TRUE, tracemgc=FALSE)
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


#' @name fit.spict
#' @title Fit a continuous-time surplus production model to data.
#' @details Fits the model using the TMB package and returns a result report containing estimates of model parameters, random effects (biomass and fishing mortality), reference points (Fmsy, Bmsy, MSY) including uncertainties given as standard deviations.
#' @param inp List of input variables as output by check.inp.
#' @param dbg Debugging option. Will print out runtime information useful for debugging if set to 1. Will print even more if set to 2.
#' @return A result report containing estimates of model parameters, random effects (biomass and fishing mortality), reference points (Fmsy, Bmsy, MSY) including uncertainties given as standard deviations.
#' @export
fit.spict <- function(inp, dbg=0){
    # Check input list
    if(!'checked' %in% names(inp)) inp <- check.inp(inp)
    if(!inp$checked) inp <- check.inp(inp)
    #require(TMB)
    scriptname <- 'spict'
    #compile(paste(scriptname,'.cpp',sep=''))
    #dyn.load(paste(scriptname,'.so',sep=''))

    datin <- list(delay=inp$delay, dt=inp$dt, dtpred=inp$dtpred, Cobs=inp$obsC, ic=inp$ic, nc=inp$nc, I=inp$obsI[[1]], ii=inp$ii[[1]], isum=rep(0,inp$ns), lamperti=inp$lamperti, euler=inp$euler, dbg=dbg)
    obj <- MakeADFun(data=datin, parameters=inp$ini, random=inp$RE, DLL=inp$scriptname, hessian=TRUE, tracemgc=FALSE, map=inp$map)
    config(trace.optimize=0, DLL=inp$scriptname)
    obj$env$tracemgc <- FALSE # Make TMB even more quiet
    obj$fn(obj$par)
    rep <- NULL
    if(dbg==0){
        opt <- nlminb(obj$par, obj$fn, obj$gr)
        # Results
        rep <- sdreport(obj)
        rep$pl <- obj$env$parList(opt$par)
        obj$fn()
        rep$Cp <- obj$report()$Cp
        rep$inp <- inp
        #osap <- one.step.ahead.pred(pl, RE, time, Cobs, timeC, I, timeI, dtpred, delay, scriptname, lamperti=lamperti, euler=euler)
    }
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
    if(exp){
        cbind(ll=exp(ll), est=exp(est), ul=exp(ul), sd)
    } else {
        cbind(ll, est, ul, sd)
    }
}


#' @name make.ellipse
#' @title Calculate confidence ellipsis.
#' @details Calculates the confidence ellipsis of two reported model parameters. This is particularly useful as a detailed view of the uncertainty of two correlated parameters.
#' @param inds Indices of the two reported model parameters.
#' @param rep A result report as generated by running fit.spict.
#' @return A matrix with two columns containing the x and y coordinates of the ellipsis.
#' @export
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


#' @name calc.Binf
#' @title Calculate B_infinity.
#' @details Calculate from the model parameters the equilibrium value of the biomass, B_infinity.
#' @param K Carrying capacity.
#' @param r Intrinsic growth rate.
#' @param sdb Standard deviation of the log biomass.
#' @param F Fishing mortality.
#' @param lamperti Return the stochastic (TRUE) or deterministic (FALSE) estimate.
#' @return The calculated estimate of B_infinity.
calc.Binf <- function(K, r, sdb, F, lamperti){
  if(lamperti){
    return( K * (1 - F/r - 0.5*sdb^2/r))
  } else {
    return( K * (1 - F/r))
  }
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
    Binf <- calc.Binf(K, r, sdb, F, lamperti)
    if(lamperti){
        rate <- F + 0.5*sdb^2 - r
    } else {
        rate <- F - r
    }
    Brat <- Binf/B0
    Brat[Brat<1] <- 1/Brat[Brat<1] # Invert if B0 > Binf
    return( 1/(rate) * log( (1-p) / (p*(Brat-1))) )
}

#' @name plotspict
#' @title 3x3 plot illustrating spict results.
#' @details Create a 3x3 plot containing the following:
#' \itemize{
#'  \item{1. Biomass.}
#'  \item{2. One-step-ahead residuals.}
#'  \item{3. One-step-ahead auto-correlation function.}
#'  \item{4. Estimated F versus estimated B.}
#'  \item{5. Estimated fishing mortality.}
#'  \item{6. Observed versus predicted catches.}
#'  \item{7. Observed versus theoretical production.}
#'  \item{8. Calculated time-constant.}
#' }
#' @param rep A result report as generated by running fit.spict.
#' @return Nothing.
#' @export
plotspict <- function(rep){
    inp <- rep$inp
    #dev.new(width=10, height=10)
    par(mfrow=c(3,3))

    # Biomass plot
    Best <- get.par('logB', rep, exp=TRUE, random=TRUE)
    ns <- dim(Best)[1]
    Kest <- get.par('logK', rep, exp=TRUE, fixed=TRUE)
    Bmsy <- get.par('logBmsy', rep, exp=TRUE)
    Binf <- get.par('logBinf', rep, exp=TRUE)
    Bp <- get.par('logBp', rep, exp=TRUE)
    scal <- 1
    cicol <- 'lightgray'
    par(mar=c(5,4,4,4))
    ylim <- range(Best[, 1:3], Bp[1:3])/scal
    plot(inp$time, Best[,2]/scal, typ='n', xlab='Time', ylab=paste('Biomass'), main=paste('- Bmsy:',round(Bmsy[2]),' K:',round(Kest[2])), ylim=ylim, xlim=range(c(inp$time, tail(inp$time,1)+1)))
    axis(4, labels=pretty(ylim/Bmsy[2]), at=pretty(ylim/Bmsy[2])*Bmsy[2])
    mtext("B/Bmsy", side=4, las=0, line=2)
    polygon(c(inp$time[1]-5,tail(inp$time,1)+5,tail(inp$time,1)+5,inp$time[1]-5), c(Bmsy[1],Bmsy[1],Bmsy[3],Bmsy[3]), col=cicol, border=cicol)
    lines(inp$time, Best[,2]/scal, col='blue')
    abline(h=Bmsy[2]/scal, col='black')
    lines(inp$time, Best[,1]/scal, col=4, lty=2)
    lines(inp$time, Best[,3]/scal, col=4, lty=2)
    lines(inp$time, Binf[,2]/scal, col=3, lty=1)
    et <- tail(inp$time,1)
    tp <- et + inp$dtpred
    lines(c(et, tp), c(tail(Best[,2],1), Bp[2])/scal, col='blue', lty=3)
    points(tp, Bp[2], pch=21, bg='yellow')
    lines(c(et, tp), c(tail(Best[,1],1), Bp[1])/scal, col='blue', lty=3)
    lines(c(et, tp), c(tail(Best[,3],1), Bp[3])/scal, col='blue', lty=3)
    legend('topright', legend=c('Binf','Bpred'), lty=c(1,NA), pch=c(NA,21), col=c(3,1), pt.bg=c(NA,'yellow'))
    box()

    # One-step-ahead catch residuals
    Cscal <- 1
    Cpred <- rep$osar$logCpred
    plot(inp$timeC, log(inp$obsC), typ='p', ylim=range(c(Cpred,log(inp$obsC),1.08*c(Cpred,log(inp$obsC)))), main=paste('Ljung-Box test p-value:',round(rep$osar$logCpboxtest$p.value,3)), ylab=paste('Catch, Cscal:',Cscal), xlim=range(c(inp$timeC,inp$timeC[inp$nobsC]+1)))
    clr <- 'black'
    lines(rep$osar$timeC, Cpred, col=clr)
    points(rep$osar$timeC, Cpred, pch=20, cex=0.7, col=clr)
    legend('topright', 'One-step pred.', lty=1, col=clr, pch=20, pt.cex=0.7)

    # OSAR ACF
    acf(rep$osar$logCpres, main='ACF of OSA residuals')

    # F versus B
    qest <- get.par('logq', rep, exp=TRUE, fixed=TRUE)
    rest <- get.par('logr', rep, exp=TRUE, fixed=TRUE)
    Fest <- get.par('logF', rep, exp=TRUE, random=TRUE)
    Fmsy <- get.par('logFmsy', rep, exp=TRUE)
    Fp <- Fest[ns,]
    inds <- c(which(names(rep$value)=='logBmsy'), which(names(rep$value)=='logFmsy'))
    cl <- make.ellipse(inds, rep)
    xlim <- range(c(exp(cl[,1]),Bp[1:3],Best[,2])/scal)
    ylim <- range(c(exp(cl[,2]),Fest[,2]))
    main <- paste('alpha:',rep$pl$alpha,'r:',round(rest[2],3),'q:',round(qest[2],3))
    par(mar=c(5,4,4,4))
    plot(Bmsy[2]/scal, Fmsy[2], typ='p', xlim=xlim, xlab='Biomass', ylab='F',  pch=24, bg='blue', ylim=ylim)
    axis(3, labels=pretty(xlim/Bmsy[2]), at=pretty(xlim/Bmsy[2])*Bmsy[2])
    mtext("B/Bmsy", side=3, las=0, line=2)
    axis(4, labels=pretty(ylim/Fmsy[2]), at=pretty(ylim/Fmsy[2])*Fmsy[2])
    mtext("F/Fmsy", side=4, las=0, line=2)
    polygon(exp(cl[,1])/scal, exp(cl[,2]), col=cicol, border=cicol)
    abline(h=Fmsy[2], lty=3)
    abline(v=Bmsy[2], lty=3)
    #arrow.line(Best[,2]/scal, Fest[,2], length=0.05, col='blue')
    lines(Best[,2]/scal, Fest[,2], col='blue')
    points(Bmsy[2]/scal, Fmsy[2], pch=24, bg='blue')
    lines(c(tail(Best[,2],1), Bp[2])/scal, rep(Fp[2],2), col='blue', lty=3)
    points(Bp[2]/scal, Fp[2], pch=21, bg='yellow')
    points(tail(Binf[,2],1)/scal, Fp[2], pch=22, bg='orange', cex=2)
    arrow.line(c(Bp[2], tail(Binf[,2],1))/scal, rep(Fp[2],2), col='black', length=0.05)
    legend('topright', c('Estimated MSY',paste(tail(inp$time,1)+1,'prediction'),'Equilibrium'), pch=c(24,21,22), pt.bg=c('blue','yellow','orange'))

    # F
    ylim <- range(Fest[,1:3])
    plot(inp$time, Fest[, 2], typ='n', main=paste('Fmsy:',round(Fmsy[2],3)), ylim=ylim, col='blue', ylab='F', xlim=range(inp$time,tail(inp$time+1,2)))
    axis(4, labels=pretty(ylim/Fmsy[2]), at=pretty(ylim/Fmsy[2])*Fmsy[2])
    mtext("F/Fmsy", side=4, las=0, line=2)
    polygon(c(inp$time[1]-5,tail(inp$time,1)+5,tail(inp$time,1)+5,inp$time[1]-5), c(Fmsy[1],Fmsy[1],Fmsy[3],Fmsy[3]), col=cicol, border=cicol)
    #abline(h=Fmax/Fmsy, col='red')
    #lines(inp$time, Fest/Fmsy, col='black')
    lines(inp$time, Fest[, 2], col='blue')
    #points(inp$time, C/Best/Fmsy)
    abline(h=Fmsy[2], col='black')
    abline(h=rest[2], col='red')
    #abline(h=Fmsyci, col=3, lty=2)
    #lines(inp$time, Festul/Fmsy, col='forestgreen', lty=2)
    #lines(inp$time, Festll/Fmsy, col='forestgreen', lty=2)
    lines(inp$time, Fest[, 1], col='blue', lty=2)
    lines(inp$time, Fest[, 3], col='blue', lty=2)
    box()

    # Catch
    Cscal <- 1
    MSY <- get.par('MSY', rep, exp=FALSE)
    Cpmsy <- get.par('Cpmsy', rep)
    Cinfp <- get.par('Cinfp', rep)
    Cpredsub <- get.par('Cpredsub', rep)
    Pest <- get.par('P', rep)
    Cpredest <- get.par('logCpred', rep, exp=TRUE)
    plot(inp$timeC, inp$obsC/Cscal, typ='n', main=paste('MSY:',round(MSY[2]/Cscal)), xlab='Time', ylab=paste('Catch'), xlim=range(c(inp$time, tail(inp$time,1)+1)), ylim=range(c(1.3*inp$obsC, Cpredest[,1:3], 0.8*inp$obsC, Cinfp[2], Cpmsy[2], rep$Cp))/Cscal)
    polygon(c(inp$time[1]-5,tail(inp$time,1)+5,tail(inp$time,1)+5,inp$time[1]-5), c(MSY[1],MSY[1],MSY[3],MSY[3])/Cscal, col=cicol, border=cicol)
    points(inp$timeC, inp$obsC/Cscal)
    points(tail(inp$timeC,1)+inp$dtpred, rep$Cp/Cscal, pch=21, bg='yellow')
    #points(tail(inp$time,1)+1, Cpmax/Cscal, pch=21, bg='red')
    points(tail(inp$timeC,1)+1, Cpmsy[2]/Cscal, pch=21, bg='green')
    points(tail(inp$timeC,1)+1, Cinfp[2]/Cscal, pch=21, bg='orange')
    #lines(inp$time[-1], Pest[, 2]/Cscal/dt, col=3)
    abline(h=MSY[2]/Cscal)
    #lines(inp$time[-ns], Cpredsub[-ns, 2]/Cscal/dt, col=4)
    lines(inp$timeC, Cpredest[, 2]/Cscal, col=4)
    lines(inp$timeC, Cpredest[, 1]/Cscal, col=4, lty=2)
    lines(inp$timeC, Cpredest[, 3]/Cscal, col=4, lty=2)
    #legend('topleft',c('Pred@Fmax','Pred@Fcur','Pred@Fmsy'),pch=21,pt.bg=c('red','yellow','green'))
    legend('topleft',c('Pred w Fcur','Pred w Fmsy','Equilibrium'),pch=21,pt.bg=c('yellow','green','orange'), bg='white')
    #legend('topleft',c('Pred@Fcur','Pred@Fmsy','Equilibrium','Production'),pch=c(21,21,21,NA),pt.bg=c('yellow','green','orange',NA),lty=c(NA,NA,NA,1), col=c(1,1,1,3))
    box()

    # Production curve
    nBplot <- 200
    Bplot <- seq(0.5*min(Best[, 2]), 1*max(Best[, 2]), length=nBplot)
    pfun <- function(r, K, B) r*B*(1 - B/K)
    Pst <- pfun(rest[2], Kest[2], Bplot)
    xlim <- range(Bplot/Bmsy[2])
    plot(Best[-1, 2]/Bmsy[2], Pest[, 2]/inp$dt/Bmsy[2], typ='l', ylim=range(Pest[,2]/inp$dt/Bmsy[2], Pst/Bmsy[2]), xlim=xlim, xlab='B/Bmsy', ylab='Production/Bmsy', col=4)
    lines(Bplot/Bmsy[2], Pst/Bmsy[2], col=1)
    #arrow.line(Best[-1, 2]/Bmsy[2], Pest[,2]/inp$dt/Bmsy[2], length=0.05)
    abline(v=1, lty=3)

    # Time constant
    sdbest <- get.par('logsdb', rep, exp=TRUE, fixed=TRUE)
    p <- 1-exp(-1) # this equals 63.2%, which is the usual definition of the time constant for linear time invariant systems
    tc <- tc.fun2(Fest[, 2], Kest[2], rest[2], sdbest[2], Best[, 2], p=0.99, lamperti=inp$lamperti)
    k <- 6
    pvec <- (0:99)/100
    npvec <- length(pvec)
    tcfix <- tc.fun2(Fest[k, 2], Kest[2], rest[2], sdbest[2], Best[k, 2], pvec, lamperti=inp$lamperti)
    Bcur <- Best[ns, 2]
    B0vec <- round(c(0.5*Bcur, 2*Bcur))
    nB0vec <- length(B0vec)
    tcmsy <- matrix(0, nB0vec, npvec)
    for(i in 1:nB0vec) tcmsy[i, ] <- tc.fun2(Fmsy[2], Kest[2], rest[2], sdbest[2], B0vec[i], pvec, lamperti=inp$lamperti)
    tccur <- tc.fun2(Fmsy[2], Kest[2], rest[2], sdbest[2], Best[ns, 2], pvec, lamperti=inp$lamperti)
    #tcmsy <- tc.fun2(Fmsy[2], Kest[2], rest[2], sdbest[2], Bmsy[2], pvec, lamperti=inp$lamperti)
    plot(tccur[tccur>0], pvec[tccur>0], typ='l', col=4, lwd=1.5, ylab='Proportion of Bmsy', xlab='Years to Bmsy when fished at Fmsy', ylim=range(c(1.3*B0vec/Bmsy[2],0.7*B0vec/Bmsy[2])), xlim=c(0, max(tcmsy)))
    abline(h=c(0.95, 1/0.95), lty=1, col='lightgray')
    abline(h=1, lty=3)
    for(i in 1:nB0vec){
        if(B0vec[i]>Bmsy[2]){
            lines(tcmsy[i,tcmsy[i,]>0], 1/pvec[tcmsy[i,]>0], typ='l', col=i+4)
        } else {
            lines(tcmsy[i,tcmsy[i,]>0], pvec[tcmsy[i,]>0], typ='l', col=i+4)
        }
    }
    abline(v=tcmsy[,which(pvec==0.95)], lty=3, col=c(5,6))
    abline(v=tccur[which(pvec==0.95)], lty=3, lwd=1.5, col=4)
    legend('topright', legend=c(paste('Bcur =',round(Best[ns,2])),paste('B =',B0vec)), lty=1, col=c(4,5,6), lwd=c(1.5,rep(1,nB0vec)), bg='white')
}


#' @import TMB
.onLoad <- function(libname, pkgname){
    #path <- paste(libname, pkgname, "TMB", "", sep="/")
    path <- paste0(system.file("TMB", package="spict"),'/')
    if(!file.exists(TMB::dynlib(paste0(path,"spict"))))
        TMB::compile(paste0(path,"spict.cpp"))
    dyn.load(TMB::dynlib(paste0(path,"spict")))
}
