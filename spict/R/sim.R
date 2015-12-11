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


#' @name predict.b
#' @title Helper function for sim.spict().
#' @param B0 Initial biomass.
#' @param F0 Fishing mortality.
#' @param gamma gamma parameter in Fletcher's Pella-Tomlinson formulation.
#' @param m m parameter in Fletcher's Pella-Tomlinson formulation.
#' @param K Carrying capacity.
#' @param n Pella-Tomlinson exponent.
#' @param dt Time step.
#' @param sdb Standard deviation of biomass process.
#' @return Predicted biomass at the end of dt.
predict.b <- function(B0, F0, gamma, m, K, n, dt, sdb){
    exp( log(B0) + (gamma*m/K - gamma*m/K*(B0/K)^(n-1.0) - F0 - 0.5*sdb^2)*dt )
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
#' # Simulate a specific number of observations
#' inp <- pol$albacore
#' inp$obsC <- NULL
#' inp$timeC <- NULL
#' inp$obsI <- NULL
#' inp$timeI <- NULL
#' set.seed(1)
#' sim <- sim.spict(inp, nobs=150)
#' repsim <- fit.spict(sim)
#' summary(repsim) # Note true values are listed in the summary
#' dev.new(width=10, height=10)
#' plot(repsim) # Note true states are shown with orange colour
#'
#' # Simulate data with seasonal F
#' inp <- list()
#' inp$dteuler <- 1/4
#' inp$nseasons <- 2
#' inp$splineorder <- 1
#' inp$obsC <- 1:80
#' inp$obsI <- 1:80
#' inp$ini <- pol$albacore$ini
#' inp$ini$logphi <- log(2) # Seasonality introduced here
#' inp <- check.inp(inp)
#' sim2 <- sim.spict(inp)
#' par(mfrow=c(2, 1))
#' plot(sim2$obsC, typ='l')
#' plot(sim2$obsI[[1]], typ='l')
#' @export
sim.spict <- function(input, nobs=100){
    # Check if input is a inp (initial values) or rep (results).
    if('par.fixed' %in% names(input)){
        #cat('Detected input as a SPiCT result, proceeding...\n')
        rep <- input
        inp <- rep$inp
        pl <- rep$pl
        plin <- inp$ini
    } else {
        if('ini' %in% names(input)){
            #cat('Detected input as a SPiCT inp, proceeding...\n')
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
            pl <- inp$parlist
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
    if('logbkfrac' %in% names(inp$ini)){
        B0 <- exp(inp$ini$logbkfrac)*exp(pl$logK)
    } else {
        B0 <- 0.5*exp(pl$logK)
    }
    if('logF0' %in% names(inp$ini)){
        F0 <- exp(inp$ini$logF0)
    } else {
        F0 <- 0.2*exp(inp$ini$logr)
    }
    m <- exp(pl$logm)
    n <- exp(pl$logn)
    gamma <- calc.gamma(n)
    K <- exp(pl$logK)
    q <- exp(pl$logq)
    logphi <- pl$logphi
    sdb <- exp(pl$logsdb)
    sdu <- exp(pl$logsdu)
    sdf <- exp(pl$logsdf)
    #alpha <- exp(pl$logalpha)
    #beta <- exp(pl$logbeta)
    sdi <- exp(pl$logsdi)
    sdc <- exp(pl$logsdc)
    #sdi <- alpha * sdb
    #sdc <- beta * sdf
    #A <- inp$A
    lambda <- exp(pl$loglambda)
    omega <- inp$omega
    
    # B[t] is biomass at the beginning of the time interval starting at time t
    # I[t] is an index of biomass (e.g. CPUE) at time t
    # P[t] is the accumulated biomass production over the interval starting at time t
    # F[t] is the constant fishing mortality during the interval starting at time t
    # C[t] is the catch removed during the interval starting at time t. 
    # obsC[j] is the catch removed over the interval starting at time j, this will typically be the accumulated catch over the year.

    if(inp$seasontype==1){ # Use spline to generate season
        seasonspline <- get.spline(pl$logphi, order=inp$splineorder, dtfine=dt)
        nseasonspline <- length(seasonspline)
    }
    flag <- TRUE
    recount <- 0
    while(flag){
        # - Fishing mortality -
        ef <- arima.sim(inp$armalistF, nt-1) * sdf*sqrt(dt) # Used to simulate other than white noise in F
        logFbase <- c(log(F0), log(F0) + cumsum(ef)) # Fishing mortality
        # Impose seasons
        season <- numeric(length(logFbase))
        if(inp$seasontype==1){ # Spline-based seasonality
            season <- seasonspline[inp$seasonindex+1]
        }
        if(inp$seasontype==2){ # This one should not be used yet!
            # These expressions are based on analytical results
            # Derivations etc. can be found in Uffe's SDE notes on p. 132
            expmosc <- function(lambda, omega, t) exp(-lambda*t) * matrix(c(cos(omega*t), -sin(omega*t), sin(omega*t), cos(omega*t)), 2, 2, byrow=TRUE)
            nu <- length(sdu)
            for(j in 1:nu){
                u <- matrix(0, 2, inp$ns)
                sduana <- sqrt(sdu[j]^2/(2*lambda)*(1-exp(-2*lambda*dt)))
                u[, 1] <- c(0.1, 0) # Set an initial state different from zero to get some action!
                omegain <- omega*j
                for(i in 2:inp$ns) u[, i] <- rnorm(2, expmosc(lambda, omegain, dt) %*% u[, i-1], sduana)
                season <- season + u[1, ]
            }
        }
        F <- exp(logFbase + season)
        # - Biomass -
        B <- rep(0,nt)
        B[1] <- B0
        e <- exp(rnorm(nt-1, 0, sdb*sqrt(dt)))
        for(t in 2:nt) B[t] <- predict.b(B[t-1], F[t-1], gamma, m[inp$ir[t]], K, n, dt, sdb) * e[t-1]
        flag <- any(B <= 0) # Negative biomass not allowed
        recount <- recount+1
        if(recount > 10) stop('Having problems simulating data where B > 0, check parameter values!')
    }
    if(any(B > 1e15)) warning('Some simulated biomass values are larger than 1e15.')
    # - Catch -
    Csub <- rep(0,nt)
    for(t in 1:nt) Csub[t] <- F[t]*B[t]*dt
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
            fac <- invlogp1(inp$ini$logp1robfac)
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
    sim$ini <- plin
    sim$do.sd.report <- inp$do.sd.report
    sim$reportall <- inp$reportall
    sim$dteuler <- inp$dteuler
    sim$splineorder <- inp$splineorder
    sim$euler <- inp$euler
    sim$lamperti <- inp$lamperti
    sim$phases <- inp$phases
    sim$priors <- inp$priors
    sim$outliers <- inp$outliers
    sim$recount <- recount
    sim$nseasons <- inp$nseasons
    sim$seasontype <- inp$seasontype
    sim$true <- pl
    sim$true$logalpha <- sim$true$logsdi - sim$true$logsdb
    sim$true$logbeta <- sim$true$logsdc - sim$true$logsdf
    sim$true$dteuler <- inp$dteuler
    sim$true$splineorder <- inp$splineorder
    sim$true$time <- time
    sim$true$B <- B
    sim$true$F <- exp(logFbase)
    sim$true$Fs <- F
    sim$true$gamma <- gamma
    sim$true$seasontype <- inp$seasontype
    
    sign <- 1
    R <- (n-1)/n*gamma*mean(m[inp$ir])/K
    p <- n-1
    sim$true$R <- R
    sim$true$logr <- log(2*R)
    # Deterministic reference points
    sim$true$Bmsyd <- K/(n^(1/(n-1)))
    sim$true$MSYd <- mean(m[inp$ir])
    sim$true$Fmsyd <- sim$true$MSYd/sim$true$Bmsyd
    # From Bordet & Rivest (2014)
    sim$true$Bmsy <- K/(p+1)^(1/p) * (1- (1+R*(p-1)/2)/(R*(2-R)^2)*sdb^2)
    sim$true$Fmsy <- R - p*(1-R)*sdb^2/((2-R)^2)
    sim$true$MSY <- K*R/((p+1)^(1/p)) * (1 - (p+1)/2*sdb^2/(1-(1-R)^2))
    sim$true$BBmsy <- B/sim$true$Bmsy
    sim$true$FFmsy <- F/sim$true$Fmsy
    #sim$true$Bmsy <- sim$true$Bmsyd * (1 - (1+R*(p-1)/2)*sdb^2 / (R*(2-R)^2))
    #sim$true$Fmsy <- sim$true$Fmsyd - (p*(1-R)*sdb^2) / ((2-R)^2)
    #sim$true$MSY <- sim$true$MSYd * (1 - ((p+1)/2*sdb^2) / (1-(1-R)^2))
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
#' @param df.out Output data frame instead of list.
#' @param summ.ex.file Save a summary example to this file (to check that parameters have correct priors or are fixed).
#' @return A list containing the results of the validation with the following keys:
#' \itemize{
#'  \item{"osarpvals"}{ P-values of the Ljung-Box test for uncorrelated one-step-ahead residuals.}
#'  \item{"*msyci"}{Logical. TRUE if the true value of B/Fmsy was inside the 95\% confidence interval for the estimate, otherwise FALSE}
#'  \item{"*msyciw"}{ Width of the 95\% confidence interval of the estimate of Bmsy/Fmsy.}
#' }
#' @examples
#' data(pol)
#' rep0 <- fit.spict(pol$albacore)
#' inp <- list()
#' inp$ini <- rep0$pl
#' set.seed(1234)
#' validate.spict(inp, nsim=10, nobsvec=c(30, 60), backup='validate.RData')
#' @export
validate.spict <- function(inp, nsim=50, nobsvec=c(15, 60, 240), estinp=NULL, backup=NULL, df.out=FALSE, summ.ex.file=NULL){
    nnobsvec <- length(nobsvec)
    inp$timeC <- NULL
    inp$timeI <- NULL
    inp$obsC <- NULL
    inp$obsI <- NULL
    if('logF' %in% names(inp$ini)) inp$ini$logF <- NULL
    if('logB' %in% names(inp$ini)) inp$ini$logB <- NULL
    ss <- list()
    require(parallel)
    #nobs <- nobsvec[1]
    fun <- function(i, inp, nobs, estinp, backup){
        cat(paste(Sys.time(), '- validating:  i:', i, 'nobs:', nobs, '\n'))
        sim <- sim.spict(inp, nobs)
        #if(!is.null(estinp)) sim$ini <- estinp$ini
        if(!is.null(estinp)) for(nm in names(estinp$priors)) sim$priors[[nm]] <- estinp$priors[[nm]]
        rep <- try(fit.spict(sim))
        s <- NA
        if(!class(rep)=='try-error'){
            s <- extract.simstats(rep, inp)
            if(!is.null(summ.ex.file)) capture.output(summary(rep), file=summ.ex.file)
        }
        return(s)
    }
    #asd <- fun(1, inp, nobs, estinp, backup)
    for(j in 1:nnobsvec){
        nobs <- nobsvec[j]
        cat(paste(Sys.time(), '- validating nobs:', nobs, '\n'))
        ss[[j]] <- mclapply(1:nsim, fun, inp, nobs, estinp, backup, mc.cores=8)
        if(!is.null(backup)) save(ss, file=backup)
    }
    if(df.out) ss <- validation.data.frame(ss)
    return(ss)
}


#' @name extract.simstats
#' @title Extracts relevant statistics from the estimation of a simulated data set.
#' @details TBA
#' @param rep A result report as generated by running fit.spict.
#' @param inp The input list used as input to the validation.spict function.
#' @return A list containing the relevant statistics.
#' @examples
#' data(pol)
#' sim <- sim.spict(pol$albacore)
#' rep <- fit.spict(sim)
#' extract.simstats(rep)
#' @export
extract.simstats <- function(rep, inp){
    if('true' %in% names(rep$inp)){
        ss <- list()
        ss$nobs <- c(nobsc=rep$inp$nobsC, nobsI=rep$inp$nobsI)
        # Convergence
        ss$conv <- rep$opt$convergence
        # Fit stats
        ss$stats <- rep$stats
        # Max min ratio
        ss$maxminratio <- rep$inp$maxminratio
        # OSA residuals p-values
        if('osar' %in% names(rep)) ss$osarpvalC <- rep$osar$logCpboxtest$p.value
        if('osar' %in% names(rep)) ss$osarpvalI <- rep$osar$logIpboxtest[[1]]$p.value
        # Estimates
        calc.simstats <- function(parname, rep, exp=TRUE, true, ind=NULL){
            par <- get.par(parname, rep, exp)
            if(!is.null(ind)){
                par <- par[ind, ]
            }
            if(is.null(dim(par))){ # par is not a matrix
                ci <- unname(true > par[1] & true < par[3])
                ciw <- unname(par[3] - par[1])
                cv <- par[5]
            } else {
                np <- dim(par)[1]
                ci <- numeric(np)
                ciw <- numeric(np)
                cv <- numeric(np)
                if(length(true)==1) true <- rep(true, np)
                for(i in 1:np){
                    ci[i] <- unname(true[i] > par[i, 1] & true[i] < par[i, 3])
                    ciw[i] <- unname(par[i, 3] - par[i, 1])
                    cv[i] <- par[i, 5]
                }
            }
            return(list(ci=ci, ciw=ciw, cv=cv))
        }
        # sdu estimate
        if('logsdu' %in% names(inp$ini) & rep$inp$phases$logsdu > 0) ss$sdu <- calc.simstats('logsdu', rep, exp=FALSE, rep$inp$true$logsdu)
        # Fmsy estimate
        ss$Fmsy <- calc.simstats('logFmsy', rep, exp=TRUE, rep$inp$true$Fmsy)
        # Bmsy estimate
        ss$Bmsy <- calc.simstats('logBmsy', rep, exp=TRUE, rep$inp$true$Bmsy)
        # MSY estimate
        ss$MSY <- calc.simstats('MSY', rep, exp=FALSE, rep$inp$true$MSY)
        # Final biomass estimate
        ind <- rep$inp$indlastobs
        ss$B <- calc.simstats('logBl', rep, exp=TRUE, rep$inp$true$B[ind])
        # Final B/Bmsy estimate
        ss$BB <- calc.simstats('logBlBmsy', rep, exp=TRUE, rep$inp$true$B[ind]/rep$inp$true$Bmsy)
        # Final F/Fmsy estimate
        ss$FF <- calc.simstats('logFlFmsy', rep, exp=TRUE, rep$inp$true$F[ind]/rep$inp$true$Fmsy)
        # Biomass process noise
        ss$sdb <- calc.simstats('logsdb', rep, exp=TRUE, exp(rep$inp$true$logsdb))
        # Biomass process noise
        ss$sdi <- calc.simstats('logsdi', rep, exp=TRUE, exp(rep$inp$true$logsdi))
        # Convergence for all values
        uss <- unlist(ss)
        ss$convall <- (any(is.na(uss) | !is.finite(uss)) | ss$conv>0)
        return(ss)
    } else {
        stop('These results do not come from the estimation of a simulated data set!')
    }
}



#' @name validation.data.frame.old
#' @title Collect results from the output of running validate.spict (this function is outdated).
#' @param ss Output from validation.spict.
#' @return A data frame containing the formatted validation results.
#' @export
validation.data.frame.old <- function(ss){
    nna <- length(ss)
    nsima <- length(ss[[1]])
    i <- 0
    flag <- TRUE
    while(flag){
        i <- i+1        
        flag <- is.null(ss[[1]][[i]])
        if(!flag){
            usso <- unlist(ss[[1]][[i]])
            nms <- names(usso)
            flag <- length(nms)<12
        }
    }
    usso <- unlist(ss[[1]][[i]])
    nms <- names(usso)
    nnms <- length(nms)
    v <- list()
    for(k in 1:nnms){
        v[[nms[k]]] <- matrix(NA, nsima, nna)
    }
    for(i in 1:nna){
        nsimt <- length(ss[[i]])
        for(j in 1:nsimt){
            uss <- unlist(ss[[i]][[j]])
            if(!is.null(uss)){
                for(k in 1:nnms){
                    v[[nms[k]]][j, i] <- uss[k]
                }
            }
        }
    }
    v[['convall']][is.na(v[['convall']])] <- 1
    return(data.frame(v))
}



#' @name validation.data.frame
#' @title Collect results from the output of running validate.spict.
#' @param ss Output from validation.spict.
#' @return A data frame containing the formatted validation results.
#' @export
validation.data.frame <- function(ss){
    uss <- unlist(ss)
    allnms <- names(uss)
    nms <- unique(allnms)
    inds <- which(nms=='')
    if(length(inds)>0) nms <- nms[-inds]
    nnms <- length(nms)
    nna <- length(ss)
    nobsvec <- numeric(nna)
    df <- data.frame(matrix(0, 1, nnms+1))
    colnames(df) <- c(nms, 'nobs')
    for(i in 1:nna){ # nna is length of nobsvec
        nsim <- length(ss[[i]])
        # While loop to avoid taking nobs from a NA run.
        c <- 1
        flag <- TRUE
        while(flag){
            flag <- all(is.na(ss[[i]][[c]])) & length(ss[[i]][[c]])==1
            c <- c+1
        }
        nobsvec[i] <- ss[[i]][[c]]$nobs[1]
        nobs <- rep(nobsvec[i], nsim)
        mat <- cbind(matrix(NA, nsim, nnms), nobs)
        colnames(mat) <- c(nms, 'nobs')
        for(j in 1:nsim){ # nsim is number of simulations for each nobs
            vals <- unlist(ss[[i]][[j]])
            namvals <- names(vals)
            for(k in 1:nnms){ # nnms is number of names
                ind <- which(namvals==nms[k])
                if(length(ind)==1) mat[j, k] <- vals[ind]
            }
        }
        df <- rbind(df, mat)
    }
    df <- df[-1, ] # Remove dummy first line
    return(df)
}
