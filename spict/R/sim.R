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
#' @param btype If 'lamperti' use Lamperti transformed equation, if 'naive' use naive formulation.
#' @return Predicted biomass at the end of dt.
predict.b <- function(B0, F0, gamma, m, K, n, dt, sdb, btype){
    if (btype == 'lamperti'){
        B <- exp( log(B0) + (gamma*m/K - gamma*m/K*(B0/K)^(n-1.0) - F0 - 0.5*sdb^2)*dt )
    }
    if (btype == 'naive'){
        # Include hack (max) to avoid negative values of biomass
        B <- max(B0 + (gamma*m*B0/K - gamma*m*(B0/K)^n - F0*B0) * dt, 1e-3)
    }
    return(B)
}


#' @name predict.logf
#' @title Helper function for sim.spict().
#' @param logF0 Fishing mortality.
#' @param dt Time step.
#' @param sdf Standard deviation of F process.
#' @param efforttype If 1 use diffusion on logF, if 2 use diffusion of F with state dependent noise (this induces the drift term -0.5*sdf^2 in log domain)
#' @return Predicted F at the end of dt.
predict.logf <- function(logF0, dt, sdf, efforttype){
    if (efforttype == 1){
        return(logF0)
    }
    # This is the Lamperti transformed F process with state dependent noise.
    if (efforttype == 2){
        return(logF0 - 0.5*sdf^2*dt)
    }
}


#' @name predict.loge
#' @title Helper function for sim.spict().
#' @param logE0 Effort
#' @param dt Time step.
#' @param sdf Standard deviation of E process.
#' @param efforttype If 1 use diffusion on logE, if 2 use diffusion of E with state dependent noise (this induces the drift term -0.5*sdf^2 in log domain)
#' @return Predicted F at the end of dt.
predict.loge <- function(logE0, dt, sdf, efforttype){
    if (efforttype == 1){
        return(logE0)
    }
    # This is the Lamperti transformed F process with state dependent noise.
    if (efforttype == 2){
        return(logE0 - 0.5*sdf^2*dt)
    }
}


#' @name predict.logmre
#' @title Helper function for sim.spict().
#' @param logmre0 Initial value
#' @param dt Time step.
#' @param sdm Standard deviation of mre process.
#' @return Predicted mre at the end of dt.
predict.logmre <- function(logmre0, dt, sdm){
    return(logmre0)
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
#' repin <- fit.spict(pol$albacore)
#' # Simulate a specific number of observations
#' inp <- list()
#' inp$dteuler <- 1/4 # To reduce calculation time
#' inp$ini <- repin$inp$ini
#' inp$ini$logF <- NULL
#' inp$ini$logB <- NULL
#' set.seed(1)
#' sim <- sim.spict(inp, nobs=150)
#' repsim <- fit.spict(sim)
#' summary(repsim) # Note true values are listed in the summary
#' plot(repsim) # Note true states are shown with orange colour
#'
#' # Simulate data with seasonal F
#' inp <- list()
#' inp$dteuler <- 1/4
#' inp$nseasons <- 2
#' inp$splineorder <- 1
#' inp$obsC <- 1:80
#' inp$obsI <- 1:80
#' inp$ini <- repin$inp$ini
#' inp$ini$logF <- NULL
#' inp$ini$logB <- NULL
#' inp$ini$logphi <- log(2) # Seasonality introduced here
#' inp <- check.inp(inp)
#' sim2 <- sim.spict(inp)
#' par(mfrow=c(2, 1))
#' plot(sim2$obsC, typ='l')
#' plot(sim2$obsI[[1]], typ='l')
#' @export
sim.spict <- function(input, nobs=100){
    # Check if input is a inp (initial values) or rep (results).
    use.effort.flag <- TRUE
    use.index.flag <- TRUE
    check.effort <- function(inp){
        nm <- names(inp)
        if (!'timeE' %in% nm){
            if (!'obsE' %in% nm){
                #inp$nobsE <- inp$nobsC
                #inp$timeE <- inp$timeC
                inp$nobsE <- 0
                inp$timeE <- numeric(0)
                inp$stdevfacE <- NULL
                inp$timeprede <- NULL
                use.effort.flag <<- FALSE
            } else {
                inp$nobsE <- length(inp$obsE)
                inp$timeE <- 1:inp$nobsE
            }
        } else {
            inp$nobsE <- length(inp$timeE)
        }
        if (!'obsE' %in% nm){
            inp$obsE <- numeric(inp$nobsE) # Insert dummy
        }
        return(inp)
    }
    check.index <- function(inp){
        if (!'logq' %in% names(inp$ini)){
            stop('logq not specified in inp$ini!')
        }
        inp$nindex <- length(inp$ini$logq)
        nm <- names(inp)
        if (!'timeI' %in% nm){
            if (!'obsI' %in% nm){
                inp$nobsI <- inp$nobsC
                inp$timeI <- inp$timeC
                inp$stdevfacI <- NULL
                inp$timepredi <- NULL
                use.index.flag <<- FALSE
                inp$nobsI <- nobs
            } else {
                if (class(inp$obsI) != 'list'){
                    tmp <- inp$obsI
                    inp$obsI <- list()
                    inp$obsI[[1]] <- tmp
                }
                inp$nobsI <- rep(0, inp$nindex)
                for (i in 1:inp$nindex){
                    inp$nobsI[i] <- length(inp$obsI[[i]])
                }
            }
            inp$timeI <- list()
            for (i in 1:inp$nindex){
                inp$timeI[[i]] <- 1:inp$nobsI[i]
            }
        } else {
            if (class(inp$timeI) != 'list'){
                tmp <- inp$timeI
                inp$timeI <- list()
                inp$timeI[[1]] <- tmp
            }
            inp$nobsI <- rep(0, inp$nindex)
            for (i in 1:inp$nindex){
                inp$nobsI[i] <- length(inp$timeI[[i]])
            }
        }
        if (!'obsI' %in% nm){
            inp$obsI <- list()
            for (i in 1:inp$nindex){
                inp$obsI[[i]] <- rep(10, inp$nobsI[i]) # Insert dummy
            }
        }
        return(inp)
    }
    zero.omit <- function(inobj){
        return(inobj[inobj != 0])
    }
    if ('par.fixed' %in% names(input)){
        #cat('Detected input as a SPiCT result, proceeding...\n')
        rep <- input
        inp <- rep$inp
        inp <- check.effort(inp)
        inp <- check.index(inp)
        inp <- check.inp(inp)
        pl <- rep$pl
        plin <- inp$ini
    } else {
        if ('ini' %in% names(input)){
            #cat('Detected input as a SPiCT inp, proceeding...\n')
            inp <- input
            nm <- names(inp)
            # Catch
            if (!'timeC' %in% nm){
                if (!'obsC' %in% nm){
                    inp$nobsC <- nobs
                } else {
                    inp$nobsC <- length(inp$obsC)
                }
                inp$timeC <- 1:inp$nobsC
            } else {
                inp$nobsC <- length(inp$timeC)
            }
            if (!'obsC' %in% nm){
                inp$obsC <- rep(10, inp$nobsC) # Insert dummy, req by check.inp().
            }
            # Effort
            inp <- check.effort(inp)
            # Index
            inp <- check.index(inp)
            # Make parameter lists and check input
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
    if (!use.effort.flag){
        use.index.flag <- TRUE
    }
    if (inp$sim.comm.cpue){
        use.index.flag <- FALSE
        use.effort.flag <- FALSE
    }

    # Calculate derived variables
    dt <- inp$dteuler
    nt <- length(time)
    if ('logbkfrac' %in% names(inp$ini)){
        B0 <- exp(inp$ini$logbkfrac) * exp(pl$logK)
    } else {
        B0 <- 0.5 * exp(pl$logK)
    }
    # Initial effort
    if ('logE0' %in% names(inp$ini)){
        E0 <- exp(inp$ini$logE0)
    } else {
        E0 <- 0.2
    }
    if (length(E0) < inp$nfleets){
        E0 <- rep(E0[1], inp$nfleets)
    }
    # Initial F (not used?)
    if ('logF0' %in% names(inp$ini)){
        F0 <- exp(inp$ini$logF0)
    } else {
        F0 <- 0.2 * exp(inp$ini$logr)
    }
    m <- exp(pl$logm)
    mre <- exp(pl$logmre)
    n <- exp(pl$logn)
    gamma <- calc.gamma(n)
    K <- exp(pl$logK)
    q <- exp(pl$logq)
    qf <- exp(pl$logqf)
    if (length(qf) == 0){
        qf <- rep(1, inp$nfisheries)
    }
    if (length(qf) != inp$nfisheries){
        qf <- rep(qf[1], inp$nfisheries)
    }
    logphi <- pl$logphi
    sdb <- exp(pl$logsdb)
    sdu <- exp(pl$logsdu)
    sdf <- exp(pl$logsdf)
    sdi <- exp(pl$logsdi)
    sdc <- exp(pl$logsdc)
    sde <- exp(pl$logsde)
    sdm <- exp(pl$logsdm)
    lambda <- exp(pl$loglambda)
    omega <- inp$omega
    
    # B[t] is biomass at the beginning of the time interval starting at time t
    # I[t] is an index of biomass (e.g. CPUE) at time t
    # P[t] is the accumulated biomass production over the interval starting at time t
    # F[t] is the constant fishing mortality during the interval starting at time t
    # C[t] is the catch removed during the interval starting at time t. 
    # obsC[j] is the catch removed over the interval starting at time j, this will typically be the accumulated catch over the year.

    if (inp$seasontype == 1){ # Use spline to generate season
        seasonspline <- get.spline(pl$logphi, order=inp$splineorder, dtfine=dt)
        nseasonspline <- length(seasonspline)
    }
    flag <- TRUE
    recount <- 0
    while(flag){
        # - Effort -
        logEbase <- matrix(0, inp$nfleets, nt)
        season <- matrix(0, inp$nfleets, nt)
        e.f <- matrix(0, inp$nfleets, nt-1)
        logEbase[, 1] <- log(E0)
        for (jj in 1:inp$nfleets){
            e.f[jj, ] <- rnorm(nt-1, 0, sdf[jj]*sqrt(dt))
            for (t in 2:nt){
                logEbase[jj, t] <- predict.loge(logEbase[jj, t-1], dt, sdf[jj],
                                                inp$efforttype) + e.f[jj, t-1]
            }
            # Impose seasons
            if (inp$seasontype == 1){ # Spline-based seasonality
                season[jj, ] <- seasonspline[inp$seasonindex + 1]
            }
            if (inp$seasontype == 2){ # This one should not be used yet!
                # These expressions are based on analytical results
                # Derivations etc. can be found in Uffe's SDE notes on p. 132
                expmosc <- function(lambda, omega, t){
                    exp(-lambda*t) * matrix(c(cos(omega*t), -sin(omega*t), sin(omega*t),
                                              cos(omega*t)), 2, 2, byrow=TRUE)
                }
                nu <- length(sdu)
                for (j in 1:nu){
                    u <- matrix(0, 2, inp$ns)
                    sduana <- sqrt(sdu[j]^2/(2*lambda) * (1-exp(-2*lambda*dt)))
                    u[, 1] <- c(0.1, 0) # Set an ini state different from zero to get some action!
                    omegain <- omega*j
                    for (i in 2:inp$ns){
                        u[, i] <- rnorm(2, expmosc(lambda, omegain, dt) %*% u[, i-1], sduana)
                    }
                    season[jj, ] <- season[jj, ] + u[1, ]
                }
            }
        }
        E <- exp(logEbase + season)
        # - Fishing mortality -
        F <- matrix(0, inp$nfisheries, nt)
        for (k in inp$nfisheriesseq){
            #qfind <- inp$targetmap[k, 3]
            flind <- inp$targetmap[k, 2]
            # If qf is not set use qf = 1
            #thisqf <- ifelse(qfind > 0, qf[qfind], 1)
            #F[k, ] <- thisqf * E[flind, ]
            F[k, ] <- qf[k] * E[flind, ]
        }
        # - Growth (time-varying via RW) -
        # Always run this even when timevaryinggrowth == FALSE to obtain same random numbers
        e.m <- matrix(0, inp$nstocks, nt-1)
        logmre <- matrix(0, inp$nstocks, nt)
        for (si in 1:inp$nstocks){
            logmre[si, 1] <- log(mre[si, 1])
            e.m[si, ] <- rnorm(nt-1, 0, sdm[si]*sqrt(dt))
            for (t in 2:nt){
                logmre[si, t] <- predict.logmre(logmre[si, t-1], dt, sdm[si]) + e.m[si, t-1]
            }
        }
        if (inp$timevaryinggrowth){
            mre <- exp(logmre) # To be used in mres below
        }
        # - Biomass -
        B <- matrix(0, inp$nstocks, nt)
        Fstock <- matrix(0, inp$nstocks, nt)
        e.b <- matrix(0, inp$nstocks, nt-1)
        #B <- numeric(nt)
        B[, 1] <- B0
        for (si in 1:inp$nstocks){
            e.b[si, ] <- exp(rnorm(nt-1, 0, sdb[si]*sqrt(dt)))
            finds <- zero.omit(inp$target[si, ])
            Fstock[si, ] <- apply(F[finds, , drop=FALSE], 2, sum)
            for (t in 2:nt){
                mres <- m[si] * mre[si, t-1]
                B[si, t] <- predict.b(B[si, t-1], Fstock[si, t-1], gamma[si], mres, K[si], n[si],
                                      dt, sdb[si], inp$btype) * e.b[si, t-1]
            }
        }
        flag <- any(B <= 0) # Negative biomass not allowed
        recount <- recount + 1
        if (recount > 10){
            stop('Having problems simulating data where B > 0, check parameter values!')
        }
    }
    if (any(B > 1e15)){
        warning('Some simulated biomass values are larger than 1e15.')
    }
    # - Catch -
    Csub <- matrix(0, inp$nfisheries, nt)
    for (k in inp$nfisheriesseq){
        bind <- inp$targetmap[k, 1]
        for (t in 1:nt){
            Csub[k, t] <- F[k, t] * B[bind, t] * dt
        }
    }
    # - Production -
    Psub <- matrix(0, inp$nstocks, nt)
    for (si in 1:inp$nstocks){
        finds <- zero.omit(inp$target[si, ])
        Csubstock <- apply(Csub[finds, , drop=FALSE], 2, sum)
        for (t in 2:nt){
            Psub[t-1] <- B[si, t] - B[si, t-1] + Csubstock[t-1]
        }
    }
    # - Catch observations -
    C <- list()
    obsC <- list()
    e.c <- list()
    for (k in inp$nfisheriesseq){
        C[[k]] <- rep(0, inp$nobsC[k])
        obsC[[k]] <- rep(0, inp$nobsC[k])
        e.c[[k]] <- exp(rnorm(inp$nobsC[k], 0, sdc[k]))
        if (inp$nobsC[k] > 0){
            # Find the indices of timeCpred that was observed
            obsinds <- match(inp$timeC[[k]], inp$timeCpred[[k]])
            ic <- inp$ic[[k]][obsinds]
            nc <- inp$nc[[k]][obsinds]
            for (i in 1:inp$nobsC[k]){
                inds <- ic[i]:(ic[i] + nc[i]-1)
                C[[k]][i] <- sum(Csub[k, inds])
                obsC[[k]][i] <- C[[k]][i] * e.c[[k]][i]
            }
        }
    }
    if ('outliers' %in% names(inp)){ # This probably doesn't work (02.08.2016)
        if ('noutC' %in% names(inp$outliers)){
            fac <- invlogp1(inp$ini$logp1robfac)
            inp$outliers$orgobsC <- obsC
            inp$outliers$indsoutC <- list()
            for (k in inp$nfisheriesseq){
                inp$outliers$indsoutC[[k]] <- sample(1:inp$nobsC[k], inp$outliers$noutC)
                obsC[[k]][inp$outliers$indsoutC[[k]]] <- exp(log(obsC[[k]][inp$outliers$indsoutC[[k]]]) + rnorm(inp$outliers$noutC, 0, fac*sdc[k]))
            }
        }
    }
    # - Effort observations -
    #Esub <- F / qf * dt
    Esub <- E * dt
    Etrue <- list()
    obsE <- list()
    e.e <- list()
    if (inp$neffort > 0){
        for (i in inp$neffortseq){
            # To which fleet does the observed effort belong?
            flind <- inp$targetmap[which(inp$targetmap[, 3] == i), 2]
            Etrue[[i]] <- numeric(inp$nobsE[i])
            obsE[[i]] <- numeric(inp$nobsE[i])
            e.e[[i]] <- exp(rnorm(inp$nobsE[i], 0, sde[i]))
            if (inp$nobsE[i] > 0){
                for (j in 1:inp$nobsE[i]){
                    inds <- inp$ie[[i]][j]:(inp$ie[[i]][j] + inp$ne[[i]][j]-1)
                    Etrue[[i]][j] <- sum(Esub[flind, inds])
                    obsE[[i]][j] <- Etrue[[i]][j] * e.e[[i]][j]
                }
            }
        }
        if ('outliers' %in% names(inp)){
            if ('noutE' %in% names(inp$outliers)){
                fac <- invlogp1(inp$ini$logp1robfae)
                inp$outliers$orgobsE <- obsE
                inp$outliers$indsoutE <- list()
                for (i in inp$neffortseq){
                    inp$outliers$indsoutE[[i]] <- sample(1:inp$nobsE[i], inp$outliers$noutE)
                    obsE[[i]][inp$outliers$indsoutE[[i]]] <- exp(log(obsE[[i]][inp$outliers$indsoutE[[i]]]) + rnorm(inp$outliers$noutE, 0, fac*sde[i]))
                }
            }
        }
    }
    # - Index observations -
    obsI <- list()
    errI <- list()
    logItrue <- list()
    Itrue <- list()
    e.i <- list()
    for (si in 1:inp$nstocks){
        obsI[[si]] <- list()
        errI[[si]] <- list()
        logItrue[[si]] <- list()
        Itrue[[si]] <- list()
        e.i[[si]] <- list()
        if (inp$nindex[si] > 0){
            for (I in inp$nindexseq[[si]]){
                obsI[[si]][[I]] <- rep(0, inp$nobsI[[si]][I])
                errI[[si]][[I]] <- rep(0, inp$nobsI[[si]][I])
                logItrue[[si]][[I]] <- numeric(inp$nobsI[[si]][I])
                sdiind <- inp$mapsdi[[si]][I]
                qind <- inp$mapq[[si]][I]
                e.i[[si]][[I]] <- exp(rnorm(inp$nobsI[[si]][I], 0, sdi[sdiind]))
                if (inp$nobsI[[si]][I] > 0){
                    for (i in 1:inp$nobsI[[si]][I]){
                        errI[[si]][[I]][i] <- rnorm(1, 0, sdi[sdiind])
                        logItrue[[si]][[I]][i] <- log(q[qind]) + log(B[inp$ii[[si]][[I]][i]])
                        obsI[[si]][[I]][i] <- exp(logItrue[[si]][[I]][i]) * e.i[[si]][[I]][i]
                    }
                    Itrue[[si]][[I]] <- exp(logItrue[[si]][[I]])
                }
            }
        }
        if ('outliers' %in% names(inp)){ # This is not done! (02.08.2016)
            if ('noutI' %in% names(inp$outliers)){
                if (length(inp$outliers$noutI) == 1){
                    inp$outliers$noutI <- rep(inp$outliers$noutI, inp$nindex)
                }
                fac <- invlogp1(inp$ini$logp1robfac)
                inp$outliers$orgobsI <- obsI
                inp$outliers$indsoutI <- list()
                if (inp$nobsI[[I]] > 0){
                    for (i in 1:inp$nindex){
                        inp$outliers$indsoutI[[i]] <- sample(1:inp$nobsI[i], inp$outliers$noutI[i])
                        obsI[[i]][inp$outliers$indsoutI[[i]]] <- exp(log(obsI[[i]][inp$outliers$indsoutI[[i]]]) +
                                                                     rnorm(inp$outliers$noutI[i], 0, fac*sdi))
                    }
                }
            }
        }
    }
    
    sim <- list()
    sim$obsC <- obsC
    sim$timeC <- inp$timeC
    if (use.index.flag){
        sim$obsI <- obsI
        sim$timeI <- inp$timeI
    } else {
        sim$otherobs$obsI <- obsI
        sim$otherobs$timeI <- inp$timeI
    }
    if (use.effort.flag){
        sim$obsE <- obsE
        sim$timeE <- inp$timeE
    } else {
        sim$otherobs$obsE <- obsE
        sim$otherobs$timeE <- inp$timeE
    }
    # Commercial CPUE observations
    if (inp$sim.comm.cpue){
        sim$obsI <- list()
        sim$timeI <- list()
    } else {
        sim$otherobs$obsIcommcpue <- list()
        sim$otherobs$timeIcommcpue <- list()
    }
    for (k in inp$nfisheriesseq){
        flind <- inp$targetmap[k, 3]
        if (flind > 0){
            if (inp$nobsC[k] == inp$nobsE[flind]){
                obsItmp <- obsC[[k]] / obsE[[flind]]
                timeItmp <- inp$timeC[[k]]
                if (inp$sim.comm.cpue){
                    sim$obsI[[k]] <- obsItmp
                    sim$timeI[[k]] <- timeItmp
                } else {
                    sim$otherobs$obsIcommcpue[[k]] <- obsItmp
                    sim$otherobs$timeIcommcpue[[k]] <- timeItmp
                }
            }
        }
    }
    # Parameters
    sim$ini <- plin
    for (nm in inp$RE){
        sim$ini[[nm]] <- NULL
    }
    # Copy these
    tocopy <- c('do.sd.report', 'reportall', 'dteuler', 'splineorder', 'euler',
                'lamperti', 'phases', 'priors', 'target', 'stocknames', 'fleetnames',
                'outliers', 'nseasons', 'seasontype', 'sim.comm.cpue', 'meyermillar',
                'aspic', 'timevaryinggrowth')
    for (nm in tocopy){
        sim[[nm]] <- inp[[nm]]
    }

    sim$recount <- recount
    # True
    sim$true <- pl
    sim$true$logalpha <- numeric(sum(inp$nindex))
    for (i in 1:sum(inp$nindex)){
        bind <- inp$index2sdb
        iind <- inp$index2sdi
        sim$true$logalpha[i] <- sim$true$logsdi[iind] - sim$true$logsdb[bind]
    }
    sim$true$logbeta <- sim$true$logsdc - sim$true$logsdf
    sim$true$dteuler <- inp$dteuler
    sim$true$splineorder <- inp$splineorder
    sim$true$time <- time
    sim$true$C <- C
    sim$true$E <- E
    sim$true$I <- Itrue
    sim$true$B <- B
    sim$true$F <- F
    sim$true$Fs <- F
    sim$true$Fstock <- Fstock
    sim$true$mre <- mre
    sim$true$gamma <- gamma
    sim$true$seasontype <- inp$seasontype
    sim$true$e.c <- e.c
    sim$true$e.e <- e.e
    sim$true$e.i <- e.i
    sim$true$e.b <- e.b
    sim$true$e.f <- e.f
    sim$true$e.m <- e.m
    
    sign <- 1
    R <- (n-1)/n * gamma * m / K
    p <- n-1
    sim$true$R <- R
    sim$true$logrold <- log(abs(gamma * m / K))
    sim$true$logr <- log(m / K * n^(n/(n-1.0)))
    sim$true$logrc <- log(2 * R)
    # Deterministic reference points
    sim$true$Bmsyd <- K/(n^(1/(n-1)))
    sim$true$MSYd <- m
    sim$true$Fmsyd <- sim$true$MSYd / sim$true$Bmsyd
    # Stochastic reference points from Bordet & Rivest (2014)
    sim$true$Bmsys <- K/(p+1)^(1/p) * (1- (1+R*(p-1)/2)/(R*(2-R)^2)*sdb^2)
    sim$true$Fmsys <- R - p*(1-R)*sdb^2/((2-R)^2)
    sim$true$MSYs <- K*R/((p+1)^(1/p)) * (1 - (p+1)/2*sdb^2/(1-(1-R)^2))
    if (inp$msytype == 's'){
        sim$true$Bmsy <- sim$true$Bmsys
        sim$true$Fmsy <- sim$true$Fmsys
        sim$true$MSY <- sim$true$MSYs
    } else {
        sim$true$Bmsy <- sim$true$Bmsyd
        sim$true$Fmsy <- sim$true$Fmsyd
        sim$true$MSY <- sim$true$MSYd
    }
    # Calculate relative B and F
    sim$true$BBmsy <- B / sim$true$Bmsy
    sim$true$FFmsy <- F / sim$true$Fmsy
    # include the log of some quantities
    lognames <- c('B', 'F', 'Bmsy', 'Fmsy', 'MSY', 'FFmsy', 'BBmsy')
    for (pn in lognames){
        sim$true[[paste0('log', pn)]] <- log(sim$true[[pn]])
    }
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
#' @param invec Vector containing the number of simulated observations of each data set in each batch.
#' @param estinp The estimation uses the true parameters as starting guess. Other initial values to be used for estimation can be specified in estinp$ini.
#' @param backup Since this procedure can be slow a filename can be specified in backup where the most recent results will be available.
#' @param df.out Output data frame instead of list.
#' @param summ.ex.file Save a summary example to this file (to check that parameters have correct priors or are fixed).
#' @param type Specify what type of information is contained in invec. If type == 'nobs' then invec is assumed to be a vector containing the number of simulated observations of each data set in each batch. If type == 'logsdc' then invec is assumed to be a vector containing values of logsdc over which to loop.
#' @param parnames Vector of parameter names to extract stats for.
#' @param exp Should exp be taken of parameters?
#' @param mc.cores Number of cores to use.
#' @param model If 'spict' estimate using SPiCT. If 'meyermillar' estimate using the model of Meyer & Millar (1999), this requires rjags and coda packages.
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
#' validate.spict(inp, nsim=10, invec=c(30, 60), backup='validate.RData')
#' @export
validate.spict <- function(inp, nsim=50, invec=c(15, 60, 240), estinp=NULL, backup=NULL,
                           df.out=FALSE, summ.ex.file=NULL, type='nobs', parnames=NULL, exp=NULL,
                           mc.cores=8, model='spict'){
    if (is.null(parnames)){
        parnames <- c('logFmsy', 'logBmsy', 'MSY', 'logBl', 'logBlBmsy',
                      'logFlFmsy', 'logsdb', 'logsdi')
    }
    ss <- list()
    #require(parallel)
    #nobs <- invec[1]
    fun <- function(i, inp, nobs, estinp, backup, type, val, parnames, exp){
        #cat(paste(Sys.time(), '- validating:  i:', i, ' type:', type, ' val:', round(val, 4)))
        sim <- sim.spict(inp, nobs)
        if (model == 'spict'){
            #if (!is.null(estinp)) sim$ini <- estinp$ini
            if (!is.null(estinp)){
                for (nm in names(estinp$priors)){
                    sim$priors[[nm]] <- estinp$priors[[nm]]
                }
            }
            rep <- try(fit.spict(sim))
            s <- NA
            str <- paste(Sys.time(), '- validating:  i:', i, ' type:', type, ' val:', round(val, 4))
            if (!class(rep)=='try-error'){
                s <- extract.simstats(rep, inp, exp=exp, parnames=parnames)
                if (!is.null(summ.ex.file)){
                    # This line causes problems when running simulation2.R, the problem is
                    # that log cannot be taken of the derout variable of the summary.
                    capture.output(summary(rep), file=summ.ex.file) 
                }
                s$type <- type
                s[[type]] <- val
                cat(str, ' convall:', as.numeric(s$convall), '\n')
            } else {
                cat(str, ' Error in fit.spict()\n')
            }
        }
        if (model == 'meyermillar'){
            sim <- check.inp(sim)
            sim$meyermillar$bugfn <- paste0('sp', i, '.bug')
            # Fit Meyer & Millar model
            res <- fit.meyermillar(sim)
            
            # Extract relevant values
            resmat <- res$resmat
            nms <- rownames(resmat)
            parnms <- c('Fmsy', 'Bmsy', 'MSY', 'BBlast', 'Blast')
            inds <- match(parnms, nms)
            resmatc <- resmat[inds, ]
            cicov <- numeric(length(parnms))
            names(cicov) <- parnms
            err <- cicov
            # Fmsy
            ind <- which(parnms == 'Fmsy')
            cicov[ind] <- resmatc[ind, 3] < sim$true$Fmsy & resmatc[ind, 4] > sim$true$Fmsy
            err[ind] <- (resmatc[ind, 5] - sim$true$Fmsy) / sim$true$Fmsy
            # Bmsy
            ind <- which(parnms == 'Bmsy')
            cicov[ind] <- resmatc[ind, 3] < sim$true$Bmsy & resmatc[ind, 4] > sim$true$Bmsy
            err[ind] <- (resmatc[ind, 5] - sim$true$Bmsy) / sim$true$Bmsy
            # MSY
            ind <- which(parnms == 'MSY')
            cicov[ind] <- resmatc[ind, 3] < sim$true$MSY & resmatc[ind, 4] > sim$true$MSY
            err[ind] <- (resmatc[ind, 5] - sim$true$MSY) / sim$true$MSY
            # BBlast
            indt <- sim$indlastobs
            true <- sim$true$B[indt] / sim$true$Bmsy
            ind <- which(parnms == 'BBlast')
            cicov[ind] <- resmatc[ind, 3] < true & resmatc[ind, 4] > true
            err[ind] <- (resmatc[ind, 5] - true) / true
            # BBlast
            true <- sim$true$B[indt]
            ind <- which(parnms == 'Blast')
            cicov[ind] <- resmatc[ind, 3] < true & resmatc[ind, 4] > true
            err[ind] <- (resmatc[ind, 5] - true) / true
            # Save results
            s <- list(ci=cicov, cv=resmatc[, 6], err=err, convall=res$diag$convall, nobs=nobs)
            # Capture output
            if (!is.null(summ.ex.file)){
                capture.output(res$resmat, file=summ.ex.file)
            }
            str <- paste(Sys.time(), '- validating:  i:', i, ' type:', type, ' val:', round(val, 4))
            cat(str, ' convall:', as.numeric(res$diag$convall), '\n')
        }
        if (model == 'aspic'){
            sim <- check.inp(sim)
            filebase <- paste0('aspic', i)
            savefile <- paste0(filebase, '.RData')
            res <- try(fit.aspic(sim, do.boot=TRUE, verbose=FALSE, filebase=filebase,
                                 savefile=savefile))
            if (class(res) == 'try-error'){
                # Model failed entirely
                s <- res
                s$errtype <- 'fiterr'
                convall <- 1
                s$convall <- convall
            } else {
                if (res$errorcode != '0'){
                    # Model ran but did not converge according to ASPIC
                    s <- res
                    s$errtype <- 'fiterr'
                    convall <- 1
                    s$convall <- convall
                } else {
                    # Model converged according to ASPIC
                    if ('boot' %in% names(res)){
                        # Bootstrap was performed
                        a <- res$boot[c(2, 3, 7, 9, 10), c(1, 4, 5)]
                        parnms <- c('MSY', 'Fmsy', 'Bmsy', 'BBlast', 'FFlast')
                        cicov <- numeric(length(parnms))
                        names(cicov) <- parnms
                        err <- c(cicov, Blast=0)
                        # MSY
                        cicov[1] <- a[1, 2] < sim$true$MSY & a[1, 3] > sim$true$MSY
                        err[1] <- (res$pars[6] - sim$true$MSY) / sim$true$MSY
                        # Fmsy
                        cicov[2] <- a[2, 2] < sim$true$Fmsy & a[2, 3] > sim$true$Fmsy
                        err[2] <- (res$pars[4] - sim$true$Fmsy) / sim$true$Fmsy
                        # Bmsy
                        cicov[3] <- a[3, 2] < sim$true$Bmsy & a[3, 3] > sim$true$Bmsy
                        err[3] <- (res$pars[5] - sim$true$Bmsy) / sim$true$Bmsy
                        # BBlast
                        indt <- sim$indlastobs
                        true <- sim$true$B[indt] / sim$true$Bmsy
                        cicov[4] <- a[4, 2] < true & a[4, 3] > true
                        err[4] <- (res$states$BBmsy[dim(res$states)[1]] - true) / true
                        # FFlast
                        true <- sim$true$F[indt] / sim$true$Fmsy
                        cicov[5] <- a[5, 2] < true & a[5, 3] > true
                        err[5] <- (res$states$FFmsy[dim(res$states)[1]] - true) / true
                        # Blast (this one doesn't have uncertainty, only an estimated value)
                        true <- sim$true$B[indt]
                        err[6] <- (res$states$B0est[dim(res$states)[1]] - true) / true
                        # Other summaries
                        cv <- (a[, 3]-a[, 2])/(2*1.96)/a[, 1]
                        names(cv) <- parnms
                        convall <- 0
                        s <- list(ci=cicov, cv=cv, err=err, convall=convall, nobs=nobs)
                    } else {
                        # Bootstrap was not performed because too slow
                        s <- res
                        s$errtype <- 'booterr'
                        convall <- 2
                        s$convall <- convall
                    }
                }
            }
            if (!is.null(summ.ex.file)){
                capture.output(res, file=summ.ex.file)
            }
            # Clean up temporary files
            unlink(paste0(filebase, '.*'))
            # Print information to screen
            str <- paste(Sys.time(), '- validating:  i:', i, ' type:', type, ' val:', round(val, 4))
            cat(str, ' convall:', convall, '\n')
        }
        if (model == 'simple'){
            # Not implemented yet
            s <- NULL
        }
        return(s)
    }
    #asd <- fun(1, inp, nobs, estinp, backup, type, val=nobs)
    ninvec <- length(invec)
    if (type == 'nobs'){
        inp$timeC <- NULL
        inp$timeI <- NULL
        inp$obsC <- NULL
        inp$obsI <- NULL
        if ('logF' %in% names(inp$ini)) inp$ini$logF <- NULL
        if ('logB' %in% names(inp$ini)) inp$ini$logB <- NULL
    }

    for (j in 1:ninvec){
        if (type == 'logsdc'){
            inp$ini$logsdc <- invec[j]
        }
        if (type == 'nobs'){
            nobs <- invec[j]
        }
        # Run simulations
        cat(paste0(Sys.time(), ' - validating w mc.cores: ', mc.cores, ' ',
                   type, ': ', round(invec[j], 4), '\n'))
        ss[[j]] <- parallel::mclapply(1:nsim, fun, inp, nobs, estinp, backup,
                                      type, invec[j], parnames, exp, mc.cores=mc.cores)
        if (!is.null(backup)){
            save(ss, file=backup)
        }
    }
    
    if (df.out & model == 'spict'){
        ss <- validation.data.frame(ss)
    }
    return(ss)
}


#' @name extract.simstats
#' @title Extracts relevant statistics from the estimation of a simulated data set.
#' @details TBA
#' @param rep A result report as generated by running fit.spict.
#' @param inp The input list used as input to the validation.spict function.
#' @param exp Should exp be taken of parameters?
#' @param parnames Vector of parameter names to extract stats for.
#' @return A list containing the relevant statistics.
#' @examples
#' data(pol)
#' repin <- fit.spict(pol$albacore)
#' sim <- sim.spict(repin)
#' rep <- fit.spict(sim)
#' extract.simstats(rep)
#' @export
extract.simstats <- function(rep, inp=NULL, exp=NULL, parnames=NULL){
    if ('true' %in% names(rep$inp)){
        if (is.null(parnames)){
            parnames <- c('logFmsy', 'logBmsy', 'MSY', 'logBl', 'logBlBmsy',
                          'logFlFmsy', 'logsdb', 'logsdi')
        }
        ss <- list()
        #ss$nobs <- c(nobsc=rep$inp$nobsC, nobsI=rep$inp$nobsI)
        ss$nobs <- list(nobsc=rep$inp$nobsC, nobsI=rep$inp$nobsI, nobsE=rep$inp$nobsE)
        # Convergence
        ss$conv <- rep$opt$convergence
        # SDreport error
        ss$sderr <- ifelse(is.null(rep$sderr), 0, 1)
        # Estimates
        calc.simstats <- function(parname, rep, exp=TRUE, true, ind=NULL){
            par <- get.par(parname, rep, exp)
            if (!is.null(ind)){
                par <- par[ind, ]
            }
            if (is.null(dim(par))){ # par is not a matrix
                ci <- unname(true > par[1] & true < par[3])
                ciw <- unname(par[3] - par[1])
                cv <- par[5]
                err <- (par[2] - true) / true
            } else {
                np <- dim(par)[1]
                ci <- numeric(np)
                ciw <- numeric(np)
                cv <- numeric(np)
                err <- numeric(np)
                if (length(true) == 1){
                    true <- rep(true, np)
                }
                for (i in 1:np){
                    ci[i] <- unname(true[i] > par[i, 1] & true[i] < par[i, 3])
                    ciw[i] <- unname(par[i, 3] - par[i, 1])
                    cv[i] <- par[i, 5]
                    err[i] <- (par[i, 2] - true[i]) / true[i]
                }
            }
            return(list(ci=ci, ciw=ciw, cv=cv, err=err, exp=exp))
        }
        # Calculate simulation statistics
        nparnames <- length(parnames)
        if (length(exp) != nparnames){
            if (length(exp) == 1){
                exp <- rep(exp, nparnames)
            } else {
                exp <- NULL
            }
        }
        if (is.null(exp)){
            inds <- grep('^log', parnames)
            exp <- logical(nparnames)
            exp[inds] <- TRUE
        }
        for (i in 1:nparnames){
            pn <- parnames[i]
            true <- rep$inp$true[[pn]]
            if (pn == 'logBl'){
                ind <- rep$inp$indlastobs
                true <- log(rep$inp$true$B[ind])
            }
            if (pn == 'logFl'){
                ind <- rep$inp$indlastobs
                true <- log(rep$inp$true$F[ind])
            }
            if (pn == 'logBlBmsy'){
                ind <- rep$inp$indlastobs
                true <- log(rep$inp$true$B[ind] / rep$inp$true$Bmsy)
            }
            if (pn == 'logFlFmsy'){
                ind <- rep$inp$indlastobs
                true <- log(rep$inp$true$F[ind] / rep$inp$true$Fmsy)
            }
            true <- ifelse(exp[i], exp(true), true)
            ss[[pn]] <- calc.simstats(pn, rep, exp=exp[i], true)
        }
        
        # Convergence for all values
        uss <- unlist(ss)
        ss$convall <- (any(is.na(uss) | !is.finite(uss)) | ss$conv > 0)
        # Residual diagnostics
        ss$diagn <- rep$diagn
        return(ss)
    } else {
        stop('These results do not come from the estimation of a simulated data set!')
    }
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
    if (length(inds) > 0){
        nms <- nms[-inds]
    }
    nnms <- length(nms)
    nna <- length(ss)
    # Initialise df (remember to remove dummy line)
    iniflag <- TRUE
    c <- 0
    while (iniflag){
        c <- c + 1
        lngt <- length(unlist(ss[[1]][[c]])) # Create with dummy first line
        iniflag <- lngt != nnms
    }
    df <- as.data.frame(ss[[1]][[c]])
    for (i in 1:nna){ # nna is length of nobsvec
        nsim <- length(ss[[i]])
        # Initialise mat (remember to remove dummy line)
        iniflag <- TRUE
        c <- 0
        while (iniflag){
            c <- c + 1
            lngt <- length(unlist(ss[[i]][[c]])) # Create with dummy first line
            iniflag <- lngt != nnms
        }
        mat <- as.data.frame(ss[[i]][[c]])
        for (j in 1:nsim){ # nsim is number of simulations for each nobs
            vals <- unlist(ss[[i]][[j]])
            if (length(vals) == nnms){
                newrow <- try(as.data.frame(ss[[i]][[j]]))
                if (class(newrow) != 'try-error' & !any(is.na(newrow))){
                    mat <- rbind(mat, newrow)
                }
            }
        }
        mat <- mat[-1, ] # Remove dummy first line create when initialising
        df <- rbind(df, mat)
    }
    df <- df[-1, ] # Remove dummy first line create when initialising
    return(df)
}

