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


#' @name summary.spictcls
#' @title Output a summary of a fit.spict() run.
#' @details The output includes the parameter estimates with 95\% confidence
#'     intervals, estimates of derived parameters (Bmsy, Fmsy, MSY) with 95\%
#'     confidence intervals, and predictions of biomass, fishing mortality, and
#'     catch.
#' @param object A result report as generated by running fit.spict.
#' @param CI Confidence intervals to be calculated, e.g. 0.9 for the 90\%
#'     confidence intervals. By default (CI = 0.95), the 95\% confidence
#'     intervals are estimated.
#' @param ... additional arguments affecting the summary produced.
#'
#' @return Nothing. Prints a summary to the screen.
#'
#' @examples
#' data(pol)
#' rep <- fit.spict(pol$albacore)
#' summary(rep)
#'
#' @export
summary.spictcls <- function(object, CI = 0.95, ...){
    if (!"digits" %in% names(list(...))) digits <- 7
    ndigits <- digits # Present values with this number of digits after the dot.
    rep <- object
    cat(paste('Convergence: ', rep$opt$convergence, '  MSG: ', rep$opt$message, '\n', sep=''))
    if (rep$opt$convergence > 0){
        cat('WARNING: Model did not obtain proper convergence! Estimates and uncertainties are most likely invalid and cannot be trusted.\n')
    }
    if (rep$opt$convergence > 0 | rep$inp$optim.method == 'SANN'){
        grad <- rep$obj$gr()
        names(grad) <- names(rep$par.fixed)
        cat('Gradient at current parameter vector\n')
        cat('', paste(capture.output(grad),' \n'))
        cat('\n')
    }
    if('sderr' %in% names(rep)) cat('WARNING: Could not calculate standard deviations. The optimum found may be invalid. Proceed with caution.\n')
    if (rep$opt$convergence > 0){
        txtobj <- 'Objective function: '
    } else {
        txtobj <- 'Objective function at optimum: '
    }
    cat(paste0(txtobj, round(rep$obj$fn(), ndigits), '\n'))
    #cat(paste0('Computing time (seconds): ', round(rep$computing.time, 3), '\n'))
    cat(paste0('Euler time step (years):  1/', round(1/rep$inp$dteuler, 2),
               ' or ', round(rep$inp$dteuler, 5), '\n'))
    str <- paste0('Nobs C: ', rep$inp$nobsC)
    if (rep$inp$nindex > 0){
        str <- paste0(str, paste0(paste0(',  Nobs I', 1:rep$inp$nindex),
                                                      ': ', rep$inp$nobsI, collapse=''))
    }
    if (rep$inp$nobsE > 0) str <- paste0(str, paste0(',  Nobs E: ', rep$inp$nobsE))
    cat(paste0(str, '\n'))
    # -- Catch/biomass unit --
    if(rep$inp$catchunit != ''){
        cat(paste('Catch/biomass unit:', rep$inp$catchunit, '\n'))
    }
    # -- Residual diagnostics --
    if (length(rep$diagn) > 0){
        cat('\nResidual diagnostics (p-values)\n')
        #cat('.: 0.1>p>0.05, *: 0.05>p>0.1, **: 0.01>p>0.001, ***: 0.001>p\n')
        diagout <- sumspict.diagnostics(rep, ndigits=4)
        cat('', paste(capture.output(diagout),' \n'))
    }
    # -- Priors --
    resout <- sumspict.priors(rep, ndigits)

    # -- Fixed parameters --
    resout <- sumspict.fixedpars(rep, ndigits=ndigits)
    if(!is.null(resout)){
        cat('\nFixed parameters\n')
        cat('', paste(capture.output(resout),' \n'))
    }
    # -- Model parameters --
    cat('\nModel parameter estimates w ', CI * 100, '% CI \n', sep = "")
    resout <- sumspict.parest(rep, ndigits=ndigits, CI = CI)
    cat('', paste(capture.output(resout),' \n'), '\n')
    if(rep$inp$do.sd.report & !'sderr' %in% names(rep)){
        # Deterministic ref points
        cat('Deterministic reference points (Drp)\n')
        derout <- sumspict.drefpoints(rep, ndigits=ndigits, CI = CI)
        cat('', paste(capture.output(derout),' \n'))
        # Stochastic derived estimates
        cat('Stochastic reference points (Srp)\n')
        derout <- sumspict.srefpoints(rep, ndigits=ndigits, CI = CI)
        cat('', paste(capture.output(derout),' \n'))
        # States
        cat('\nStates w ', CI * 100, '% CI (inp$msytype: ', rep$inp$msytype, ')\n', sep = "")
        stateout <- sumspict.states(rep, ndigits=ndigits, CI = CI)
        cat('', paste(capture.output(stateout),' \n'))
        # Predictions
        if(rep$inp$reportall){
            cat('\nPredictions w ', CI * 100, '% CI (inp$msytype: ', rep$inp$msytype, ')\n', sep = "")
            predout <- sumspict.predictions(rep, ndigits=ndigits, CI = CI)
            cat('', paste(capture.output(predout),' \n'))
        } else {
            cat(paste0('\nPredictions omitted because inp$reportall = FALSE\n'))
        }
    }
}


#' @name get.order
#' @title Get order of printed quantities.
#' @return Vector containing indices of printed quantities.
get.order <- function() return(c(2, 1, 3, 2))


#' @name get.colnms
#' @title Get column names for data.frames.
#' @return Vector containing column names of data frames.
get.colnms <- function() return(c('estimate', 'cilow', 'ciupp', 'log.est'))


#' @name sumspict.parest
#' @title Parameter estimates of a fit.spict() run.
#' @param rep A result report as generated by running fit.spict.
#' @param ndigits Present values with this number of digits after the dot.
#' @param CI Confidence intervals to be calculated, e.g. 0.9 for the 90\%
#'     confidence intervals. By default (CI = 0.95), the 95\% confidence
#'     intervals are estimated.
#'
#' @return data.frame containing parameter estimates.
#'
#' @export
sumspict.parest <- function(rep, ndigits=8, CI = 0.95){
    if(CI > 1 || CI < 0) stop("CI has to be between 0 and 1!")
    zscore <- qnorm(CI + (1 - CI)/2)
    if(rep$inp$do.sd.report){
        order <- get.order()
        colnms <- get.colnms()
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
        cilow <- rep$par.fixed-zscore*sd
        cilow[loginds] <- exp(cilow[loginds])
        cilow[logitinds] <- invlogit(cilow[logitinds])
        cilow[logp1inds] <- invlogp1(cilow[logp1inds])
        ciupp <- rep$par.fixed+zscore*sd
        ciupp[loginds] <- exp(ciupp[loginds])
        ciupp[logitinds] <- invlogit(ciupp[logitinds])
        ciupp[logp1inds] <- invlogp1(ciupp[logp1inds])
        if('true' %in% names(rep$inp)){
            npar <- length(nms)
            unms <- unique(nms)
            nupar <- length(unms)
            truepar <- NULL
            parnotest <- NULL
            for(i in 1:nupar){
                tp <- rep$inp$true[[unms[i]]]
                nestpar <- sum(names(est) == unms[i])
                truepar <- c(truepar, tp[1:nestpar])
                if(nestpar < length(tp)){
                    inds <- (nestpar+1):length(tp)
                    parnotest <- c(parnotest, tp[inds])
                    names(parnotest) <- c(names(parnotest), paste0(unms[i], inds))
                }
            }
            truepar[loginds] <- exp(truepar[loginds])
            truepar[logitinds] <- invlogit(truepar[logitinds])
            truepar[logp1inds] <- invlogp1(truepar[logp1inds])
            ci <- rep(0, npar)
            for(i in 1:npar){
                ci[i] <- as.numeric(truepar[i] > cilow[i] & truepar[i] < ciupp[i])
            }
            resout <- cbind(estimate=round(est,ndigits),
                            true=round(truepar,ndigits),
                            cilow=round(cilow,ndigits),
                            ciupp=round(ciupp,ndigits),
                            true.in.ci=ci,
                            log.est=round(rep$par.fixed,ndigits))
        } else {
            resout <- cbind(estimate=round(est,ndigits),
                            cilow=round(cilow,ndigits),
                            ciupp=round(ciupp,ndigits),
                            log.est=round(rep$par.fixed,ndigits))
        }
        nms[loginds] <- sub('log', '', names(rep$par.fixed[loginds]))
        nms[logitinds] <- sub('logit', '', names(rep$par.fixed[logitinds]))
        nms[logp1inds] <- sub('logp1', '', names(rep$par.fixed[logp1inds]))
        unms <- unique(nms)
        for(inm in unms){
            nn <- sum(inm==nms)
            if(nn>1){
                newnms <- paste0(inm, 1:nn)
                inds <- which(inm==nms)
                nms[inds] <- newnms
            }
        }
        rownames(resout) <- nms
        # Derived variables
        if(!'sderr' %in% names(rep)){
            nalpha <- sum(names(rep$par.fixed) == 'logsdi')
            derout <- rbind(get.par(parname='logbeta', rep, exp=TRUE, CI = CI)[, order],
                            get.par(parname='logr', rep, exp=TRUE, CI = CI)[, order],
                            get.par(parname='logrc', rep, exp=TRUE, CI = CI)[, order],
                            get.par(parname='logrold', rep, exp=TRUE, CI = CI)[, order])
            if (nalpha > 0){
                derout <- rbind(get.par(parname='logalpha', rep, exp=TRUE, CI = CI)[1:nalpha, order],
                                derout)
            }
            derout[, 4] <- log(derout[, 4])
            derout <- round(derout, ndigits)
            nr <- dim(derout)[1]
            if('true' %in% names(rep$inp)){
                dertrue <- exp(c(rep$inp$true$logbeta, rep$inp$true$logr))
                if (nalpha > 0){
                    dertrue <- c(rep(exp(rep$inp$true$logalpha), nalpha), dertrue)
                }
                ndertrue <- length(dertrue)
                if(ndertrue == nr){
                    cider <- numeric(ndertrue)
                    for(i in 1:ndertrue) cider[i] <- as.numeric(dertrue[i] > derout[i, 2] & dertrue[i] < derout[i, 3])
                } else {
                    dertrue <- rep(-9, nr)
                    cider <- rep(-9, nr)
                }
                derout <- cbind(est=derout[, 1], true=dertrue, ll=derout[, 2], ul=derout[, 3], tic=cider, eil=derout[, 4])
            }
            if(nr>1 & 'yearsepgrowth' %in% names(rep$inp)){
                rnms <- c('r   ', paste0('r', rep$inp$yearsepgrowth))
                roldnms <- c('rold   ', paste0('rold', rep$inp$yearsepgrowth))
                rcnms <- c('rc   ', paste0('rc', rep$inp$yearsepgrowth))
            } else {
                rnms <- rep('r  ',length(rep$inp$ini$logr))
                roldnms <- rep('rold  ', length(rep$inp$ini$logr))
                rcnms <- rep('rc  ', length(rep$inp$ini$logr))
            }
            if (nalpha > 0){
                if(nalpha > 1){
                    alphanms <- paste0('alpha', 1:nalpha)
                } else {
                    alphanms <- 'alpha'
                }
            } else {
                alphanms <- NULL
            }
            rownames(derout) <- c(alphanms, 'beta', rnms, rcnms, roldnms)
            resout <- rbind(derout, resout)
        }
        if('true' %in% names(rep$inp)){
            colnames(resout) <- c(colnms[1], 'true', colnms[2:3], 'true.in.ci', colnms[4])
        } else {
            colnames(resout) <- colnms
        }
    } else {
        if('opt' %in% names(rep)) resout <- data.frame(estimate=rep$opt$par)
    }
    return(resout)
}


#' @name sumspict.drefpoints
#' @title Deternistic reference points of a fit.spict() run.
#' @param rep A result report as generated by running fit.spict.
#' @param ndigits Present values with this number of digits after the dot.
#' @param CI
#' @return data.frame containing deterministic reference points.
#' @export
sumspict.drefpoints <- function(rep, ndigits=8, CI = 0.95){
    order <- get.order()
    colnms <- get.colnms()
    derout <- rbind(get.par(parname='logBmsyd', rep, exp=TRUE, CI = CI)[,order],
                    get.par(parname='logFmsyd', rep, exp=TRUE, CI = CI)[,order],
                    get.par(parname='logMSYd', rep, exp=TRUE, CI = CI)[,order])
    derout[, 4] <- log(derout[, 4])
    derout <- round(derout, ndigits)
    colnames(derout) <- colnms
    nr <- length(rep$inp$ini$logr)
    if(nr > 1){
        rownames(derout) <- c(t(outer(c('Bmsyd', 'Fmsyd', 'MSYd'), 1:2, paste0)))
    } else {
        rownames(derout) <- c('Bmsyd', 'Fmsyd', 'MSYd')
    }
    if('true' %in% names(rep$inp)){
        trueder <- c(rep$inp$true$Bmsyd, rep$inp$true$Fmsyd, rep$inp$true$MSYd)
        cider <- numeric(3)
        for(i in 1:3) cider[i] <- as.numeric(trueder[i] > derout[i, 2] & trueder[i] < derout[i, 3])
        derout <- cbind(derout[, 1], round(trueder,ndigits), derout[, 2:3], cider, derout[, 4])
        colnames(derout) <- c(colnms[1], 'true', colnms[2:3], 'true.in.ci', colnms[4])
    }
    return(derout)
}


#' @name sumspict.srefpoints
#' @title Stochastic reference points of a fit.spict() run.
#' @param rep A result report as generated by running fit.spict.
#' @param ndigits Present values with this number of digits after the dot.
#' @param CI Confidence intervals to be calculated, e.g. 0.9 for the 90\%
#'     confidence intervals. By default (CI = 0.95), the 95\% confidence
#'     intervals are estimated.
#' @return data.frame containing stochastic reference points.
#' @export
sumspict.srefpoints <- function(rep, ndigits=8, CI = 0.95){
    order <- get.order()
    colnms <- get.colnms()
    derout <- rbind(get.par(parname='logBmsys', rep, exp=TRUE, CI = CI)[,order],
                    get.par(parname='logFmsys', rep, exp=TRUE, CI = CI)[,order],
                    get.par(parname='logMSYs', rep, exp=TRUE, CI = CI)[,order])
    derout[, 4] <- log(derout[, 4])
    derout <- round(derout, ndigits)
    colnames(derout) <- colnms
    nr <- length(rep$inp$ini$logr)
    if(nr > 1){
        rownames(derout) <- c(t(outer(c('Bmsys', 'Fmsys', 'MSYs'), 1:2, paste0)))
    } else {
        rownames(derout) <- c('Bmsys', 'Fmsys', 'MSYs')
    }
    if('true' %in% names(rep$inp)){
        trueder <- c(rep$inp$true$Bmsy, rep$inp$true$Fmsy, rep$inp$true$MSY)
        cider <- rep(0, 3)
        for(i in 1:3) cider[i] <- as.numeric(trueder[i] > derout[i, 2] & trueder[i] < derout[i, 3])
        derout <- cbind(derout[, 1], round(trueder,ndigits), derout[, 2:3], cider, derout[, 4])
        colnames(derout) <- c(colnms[1], 'true', colnms[2:3], 'true.in.ci', colnms[4])
    }
    Drp <- c(get.par(parname='logBmsyd', rep, exp=TRUE, CI = CI)[, 2],
             get.par(parname='logFmsyd', rep, exp=TRUE, CI = CI)[, 2],
             get.par(parname='logMSYd', rep, exp=TRUE, CI = CI)[, 2])
    rel.diff.Drp <- (derout[, 1] - Drp)/derout[, 1]
    if(length(rel.diff.Drp) == dim(derout)[1]) derout <- cbind(derout, rel.diff.Drp)
    return(derout)
}


#' @name sumspict.states
#' @title State estimates of a fit.spict() run.
#' @param rep A result report as generated by running fit.spict.
#' @param ndigits Present values with this number of digits after the dot.
#' @return data.frame containing state estimates.
#' @export
sumspict.states <- function(rep, ndigits=8, CI = 0.95){
    order <- get.order()
    colnms <- get.colnms()
    idx <- rep$obj$env$data$indlastobs
    stateout <- rbind(
        get.par(parname='logBl', rep, exp=TRUE, CI = CI)[order],
        get.par(parname='logFnotS', rep, exp=TRUE, CI = CI)[idx,order],
        get.par(parname='logBlBmsy', rep, exp=TRUE, CI = CI)[order],
        get.par(parname='logFFmsynotS', rep, exp=TRUE, CI = CI)[idx,order])
    stateout[, 4] <- log(stateout[, 4])
    stateout <- round(stateout, ndigits)
    colnames(stateout) <- colnms
    indl <- rep$inp$indlastobs
    et <- fd(rep$inp$time[indl])
    rownames(stateout) <- c(paste0('B_',et), paste0('F_',et), paste0('B_',et,'/Bmsy'), paste0('F_',et,'/Fmsy'))
    if('true' %in% names(rep$inp)){
        truest <- c(rep$inp$true$B[indl], rep$inp$true$F[indl], rep$inp$true$B[indl]/rep$inp$true$Bmsy,
                    rep$inp$true$F[indl]/rep$inp$true$Fmsy)
        nst <- length(truest)
        cist <- numeric(nst)
        for (i in 1:nst){
            cist <- as.numeric(truest[i] > stateout[i, 2] & truest[i] < stateout[i, 3])
        }
        stateout <- cbind(stateout[, 1], round(truest, ndigits), stateout[, 2:3],
                          cist, stateout[, 4])
        colnames(stateout) <- c(colnms[1], 'true', colnms[2:3], 'true.in.ci', colnms[4])
    }
    return(stateout)
}


#' @name sumspict.predictions
#' @title Predictions of a fit.spict() run.
#' @param rep A result report as generated by running fit.spict.
#' @param ndigits Present values with this number of digits after the dot.
#' @return data.frame containing predictions.
#' @export
sumspict.predictions <- function(rep, ndigits=8, CI = 0.95){
    order <- get.order()
    colnms <- get.colnms()
    EBinf <- get.EBinf(rep)
    idx <- rep$obj$env$data$dtprediind
    predout <- rbind(
        get.par(parname='logBp', rep, exp=TRUE, CI = CI)[order],
        get.par(parname='logFnotS', rep, exp=TRUE, CI = CI)[idx,order],
        get.par(parname='logBpBmsy', rep, exp=TRUE, CI = CI)[order],
        get.par(parname='logFFmsynotS', rep, exp=TRUE, CI = CI)[idx,order],
        tail(get.par(parname='logCpred', rep, exp=TRUE, CI = CI),1)[order],
        c(EBinf, NA, NA, EBinf))
    inds <- predout[, 4] <= 0
    predout[inds, 4] <- NA
    inds <- predout[, 4] > 0 & !is.na(predout[, 4])
    predout[inds, 4] <- log(predout[inds, 4])
    predout <- round(predout, ndigits)
    colnames(predout) <- c('prediction', colnms[2:4])
    indp <- rep$inp$dtprediind
    et <- fd(rep$inp$time[indp])
    rownames(predout) <- c(paste0('B_',et), paste0('F_',et), paste0('B_',et,'/Bmsy'),
                           paste0('F_',et,'/Fmsy'), paste0('Catch_', fd(tail(rep$inp$timeCpred,1))), 'E(B_inf)')
    if(rep$inp$dtpredc == 0){
        predout <- predout[-dim(predout)[1], ]
    }
    return(predout)
}


#' @name print.spictcls
#' @title Output a summary of a fit.spict() run.
#' @param x A result report as generated by running fit.spict.
#' @param ... additional arguments affecting the summary produced.
#' @return Nothing.
#' @export
print.spictcls <- function(x, ...){
    if('timeC' %in% names(x)){
        cat('Catch observations:\n')
        print(x$timeC)
        print(x$obsC)
        cat('Index observations:\n')
        print(x$timeI)
        print(x$obsI)
    }
    if('par.fixed' %in% names(x)) summary(x)
}


#' @name sumspict.fixedpars
#' @title Fixed parameters table.
#' @param rep A result report as generated by running fit.spict.
#' @param ndigits Present values with this number of digits after the dot.
#' @return data.frame containing fixed parameter information.
#' @export
sumspict.fixedpars <- function(rep, ndigits=8){
    inds <- which(unlist(rep$inp$phases) < 0)
    nms <- names(rep$inp$phases)[inds]
    # Remove random effects
    reinds <- which(nms %in% rep$inp$RE)
    if(length(reinds)>0) nms <- nms[-reinds]
    # Were index observations used?
    if (sum(rep$inp$nobsI) == 0){
        nms <- nms[-match(c('logsdi', 'logq'), nms)]
    }
    # Were effort observations used?
    if (rep$inp$nobsE == 0){
        nms <- nms[-match(c('logsde', 'logqf'), nms)]
    }
    # Are robust options used? if not remove
    if(!any(rep$inp$robflagi | rep$inp$robflagc | rep$inp$robflage)){
        nms <- nms[-match(c('logitpp', 'logp1robfac'), nms)]
    }
    # Are seasonal spline not used? if yes remove
    if(! rep$inp$seasontype %in% c(1, 3)){
        nms <- nms[-match('logphi', nms)]
    }
    # Are seasonal AR + spline used? if not remove
    if(rep$inp$seasontype != 3){
        nms <- nms[-match(c('logitSARphi','logSdSAR'), nms)]
    }
    # Are seasonal SDE used? if not remove
    if(rep$inp$seasontype != 2){
        nms <- nms[-match(c('logsdu', 'loglambda'), nms)]
    }
    # Is RW effort model used?
    if (rep$inp$effortmodel == 'RW'){
        nms <- nms[-match(c('logeta', 'logdelta'), nms)]
    }
    # Is covariate information used for logm?
    if (!rep$inp$logmcovflag){
        nms <- nms[-match('mu', nms)]
    }
    # Is growth time varying
    if (!rep$inp$timevaryinggrowth){
        nms <- nms[-match(c('logsdm', 'logpsi'),  nms)]
    }
    nnms <- length(nms)
    if(nnms > 0){
        vals <- numeric(0)
        valnms <- character(0)
        for(i in 1:nnms){
            val <- get.par(parname=nms[i], rep)[2]
            vals <- c(vals, val)
            nval <- length(val)
            if(nval>1){
                valnms <- c(valnms, paste0(nms[i], 1:nval))
            } else {
                valnms <- c(valnms, nms[i])
            }
        }
        names(vals) <- valnms
        vals <- trans2real(vals, nms)
        df <- data.frame(fixed.value=vals)
        df <- round(df, ndigits)

        if('true' %in% names(rep$inp)){
            alltrue <- unlist(rep$inp$true)
            inds <- match(nms, names(alltrue))
            truevals <- alltrue[inds]
            truenms <- names(truevals)
            truevals <- trans2real(truevals, truenms, chgnms=FALSE)
            df <- cbind(df, true=truevals)
        }
        return(df)
    } else {
        return(NULL)
    }
}


#' @name trans2real
#' @title Get real parameter values from transformed ones.
#' @param vals Parameters in transformed domain.
#' @param nms Names of transformed parameters (including log etc.)
#' @param chgnms Remove transformation indication from the parameter names (e.g. remove log from logK).
#' @return Parameter values in the natural domain.
#' @export
trans2real <- function(vals, nms, chgnms=TRUE){
    loginds <- grep('log', nms)
    logp1inds <- grep('logp1',nms)
    logitinds <- grep('logit',nms)
    loginds <- setdiff(loginds, c(logp1inds, logitinds))
    vals[loginds] <- exp(vals[loginds])
    vals[logitinds] <- invlogit(vals[logitinds])
    vals[logp1inds] <- invlogp1(vals[logp1inds])
    if(chgnms){
        valnms <- names(vals)
        valnms[logitinds] <- gsub('logit', '', valnms[logitinds])
        valnms[logp1inds] <- gsub('logp1', '', valnms[logp1inds])
        valnms[loginds] <- gsub('log', '', valnms[loginds])
        names(vals) <- valnms
    }
    return(vals)
}


#' @name sumspict.priors
#' @title Fixed parameters table.
#' @param rep A result report as generated by running fit.spict.
#' @param ndigits Present values with this number of digits after the dot.
#' @return data.frame containing fixed parameter information.
#' @export
sumspict.priors <- function(rep, ndigits=8){
  indso <- which(rep$inp$priorsuseflag == 1)
  if(length(indso) > 0){
    cat('\nPriors\n')
    priors <- rep$inp$priors[indso]
    usepriors <- names(priors)

    usepriors <- gsub('gamma', '', usepriors) # Strip gamma-text away
    npriors <- length(usepriors)
    # RE priors
    repriors <- c('logB', 'logF', 'logBBmsy', 'logFFmsy')
    if(any(repriors %in% usepriors)){
      inds <- na.omit(match(repriors, usepriors))
      for(i in 1:length(inds)){
        names(priors)[inds[i]] <- paste0(usepriors[inds[i]], fd(priors[[inds[i]]][4]))
        priors[[inds[i]]] <- priors[[inds[i]]][1:3] # Remove additional columns
      }
    }
    # Matrix priors
    matpriors <- rep$inp$matrixpriors
    nmmatpriors <- names(matpriors)
    priorsmat <- do.call(rbind, priors[! names(priors) %in% nmmatpriors])
    priorsmat <- rbind( priorsmat, do.call(rbind, matpriors))
    storage.mode(priorsmat)<-"double"
    priorsmat <- priorsmat[ priorsmat[,3] == 1 , , drop=FALSE]
    gammainds <- grep('gamma', rownames(priorsmat))
    npriors <- nrow(priorsmat)
    usepriors <- rownames(priorsmat)
    str <- character(npriors)
    maxchar <- max(nchar(usepriors))
    for(i in 1:npriors){
      if (i %in% gammainds){
        shape <- priorsmat[i,1]
        rate <- priorsmat[i,2]
        vec <- shaperate2meanvar(shape, rate)
        str[i] <- paste0('~  dgamma[', round(shape, 3),
                         ', ', round(rate, 3), '] (mean=', round(vec[1], 3), ', sd=', round(vec[3], 3), ')')
      } else {
        str[i] <- paste0('~  dnorm[log(', round(exp(priorsmat[i, 1]), 3),
                         '), ', round(priorsmat[i, 2], 3), '^2]',
                         ifelse(priorsmat[i, 2] <= 1e-3, ' (fixed)', ''))
      }
      usepriors[i] <- formatC(usepriors[i], width = maxchar, flag = 0)
      cat(paste0(' ', usepriors[i], '  ', str[i], '\n'))
    }
  } else {
    cat('\nNo priors are used\n')
  }
}


#' @name sumspict.diagnostics
#' @title Diagnostics table
#' @param rep A result report as generated by running fit.spict.
#' @param ndigits Present values with this number of digits after the dot.
#' @return data.frame containing diagnostics information.
#' @export
sumspict.diagnostics <- function(rep, ndigits=8) {
    if(! is(rep, "spictcls")) {
      stop("rep should be a spictcls object as it is returned by fit.spict")
    }
    if (is.null(rep$diagn)) {
      stop("rep does not contain diagnostics, calculate them using calc.osa.resid and/or calc.process.resid.")
    }
    # Diagnostics matrix
    diagn <- rep$diagn
    nms <- names(diagn)
    if (FALSE){
        # Exclude Ljung-Box test for now, testing phase
        inds <- grep('LB', nms)
        for (i in inds){
            diagn[[nms[i]]] <- NULL
        }
    }
    # Continue with LB removed
    nms <- names(diagn)
    tests <- unique(substr(nms, 1, 3))
    ntests <- length(tests)
    if(any(names(rep) == "process.resid") && any(names(rep) == "osar")){
        nseries <- rep$inp$nseries + 2
        testnames <- sub('C.p', '', nms[grep("C.p", nms)])
    }else if(any(names(rep) == "process.resid")){
        nseries <- 2
        testnames <- sub('F.p', '', nms[grep("F.p", nms)])
    }else{
        nseries <- rep$inp$nseries
        testnames <- sub('C.p', '', nms[grep("C.p", nms)])
    }
    seriesnames <- sub('.p', '', sub(testnames[1], '', grep(testnames[1], nms, value=TRUE)))
    diagnmat <- matrix(round(unlist(diagn), ndigits), nseries, ntests, byrow=TRUE)
    colnames(diagnmat) <- testnames
    rownames(diagnmat) <- seriesnames
    # Stars
    stars <- array('-', dim=dim(diagnmat))
    colnames(stars) <- testnames
    stars[diagnmat < 0.1] <- '.'
    stars[diagnmat < 0.05] <- '*'
    stars[diagnmat < 0.01] <- '**'
    stars[diagnmat < 0.001] <- '***'
    df <- cbind(as.data.frame(diagnmat), as.data.frame(stars))
    return(df)
}

#' @name sumspict.ini
#' @title Sensitivity to the initial parameter values
#' @param rep A result report as generated by running fit.spict.
#' @param numdigits Present values with this number of digits after the dot.
#' @return list containing diagnostics information.
#' @export
sumspict.ini <- function(rep, numdigits) {
  if (! is(rep, "spictcls")) {
    stop("rep should be a spictcls object as it is returned by fit.spict")
  }
  if (is.null(rep$check.ini)) {
    stop("rep does not contain results of the sensitivity analysis to inital parameter values, calculate them using check.ini.")
  }
  return(rep$check.ini)
}
