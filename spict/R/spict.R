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


#' @name fit.spict
#' @title Fit a continuous-time surplus production model to data.
#' @details Fits the model using the TMB package and returns a result report containing estimates of model parameters, random effects (biomass and fishing mortality), reference points (Fmsy, Bmsy, MSY) including uncertainties given as standard deviations.
#'
#' Model parameters using the formulation of Fletcher (1978):
#' \itemize{
#'   \item{"logn"}{ Parameter determining the shape of the production curve as in the generalised form of Pella & Tomlinson (1969).}
#'   \item{"logm"}{ Log of maximum sustainable yield.}
#'   \item{"logK"}{ Log of carrying capacity.}
#'   \item{"logq"}{ Log of catchability vector.}
#'   \item{"logsdb"}{ Log of standard deviation of biomass process error.}
#'   \item{"logsdf"}{ Log of standard deviation of fishing mortality process error.}
#'   \item{"logsdi"}{ Log of standard deviation of index observation error.}
#'   \item{"logsdc"}{ Log of standard deviation of catch observation error.}
#' }
#'
#' Unobserved states estimated as random effects:
#' \itemize{
#'   \item{"logB"}{ Log of the biomass process given by the stochastic differential equation: dB_t = r*B_t*(1-(B_t/K)^n)*dt + sdb*dW_t, where dW_t is Brownian motion.}
#'   \item{"logF"}{ Log of the fishing mortality process given by: dlog(F_t) = f(t, sdf), where the function f depends on the choice of seasonal model.}
#' }
#'
#' Other parameters (which are only needed in certain cases):
#' \itemize{
#'   \item{"logphi"}{ Log of parameters used to specify the cyclic B spline representing seasonal variation. Used when inp$nseasons > 1 and inp$seasontype = 1.}
#'   \item{"logU"}{ Log of the state of the coupled SDE system used to represent seasonal variation, i.e. when inp$nseasons > 1 and inp$seasontype = 2.}
#'   \item{"loglambda"}{ Log of damping parameter when using the coupled SDE system to represent seasonal variation, i.e. when inp$nseasons > 1 and inp$seasontype = 2.}
#'   \item{"logsdu"}{ Log of standard deviation of process error of U_t (the state of the coupled SDE system) used to represent seasonal variation, i.e. when inp$nseasons > 1 and inp$seasontype = 2.}
#'   \item{"logsde"}{ Log of standard deviation of observation error of effort data. Only used if effort data is part of input.}
#'   \item{"logp1robfac"}{ Log plus one of the coefficient to the standard deviation of the observation error when using a mixture distribution robust toward outliers, i.e. when either inp$robflag = 1 and/or inp$robflagi = 1.}
#'   \item{"logitpp"}{Logit of the proportion of narrow distribution when using a mixture distribution robust toward outliers, i.e. when either inp$robflag = 1 and/or inp$robflagi = 1.}
#' }
#'
#' Parameters that can be derived from model parameters:
#' \itemize{
#'   \item{"logr"}{ Log of intrinsic growth rate (r = 4m/K).}
#'   \item{"logalpha"}{ Proportionality factor for the observation noise of the indices and the biomass process noise: sdi = exp(logalpha)*sdb. (normally set to logalpha=0)}
#'   \item{"logbeta"}{ Proportionality factor for the observation noise of the catches and the fishing mortality process noise: sdc = exp(logbeta)*sdf. (this is often difficult to estimate and can result in divergence of the optimisation. Normally set to logbeta=0)}
#'   \item{"logBmsy"}{ Log of the equilibrium biomass (Bmsy) when fished at Fmsy.}
#'   \item{"logFmsy"}{ Log of the fishing mortality (Fmsy) leading to the maximum sustainable yield.}
#'   \item{"MSY"}{ The yield when the biomass is at Bmsy and the fishing mortality is at Fmsy, i.e. the maximum sustainable yield.}
#' }
#'
#' The above parameter values can be extracted from the fit.spict() results using get.par().
#'
#' Model assumptions
#' \itemize{
#'   \item{"1"}{The intrinsic growth rate (r) represents a combination of natural mortality, growth, and recruitment.}
#'   \item{"2"}{The biomass B_t refers to the exploitable part of the stock. Estimates in absolute numbers (K, Bmsy, etc.) should be interpreted in light of this.}
#'   \item{"3"}{The stock is closed to migration.}
#'   \item{"4"}{Age and size-distribution are stable in time.}
#'   \item{"5"}{Constant catchability of the gear used to gather information for the biomass index.}
#' }
#' @param inp List of input variables as output by check.inp.
#' @param verbose Should detailed outputs be provided (default: TRUE).
#' @param dbg Debugging option. Will print out runtime information useful for debugging if set to 1. Will print even more if set to 2.
#' @return A result report containing estimates of model parameters, random effects (biomass and fishing mortality), reference points (Fmsy, Bmsy, MSY) including uncertainties given as standard deviations.
#' @export
#' @examples
#' data(pol)
#' rep <- fit.spict(pol$albacore)
#' Bmsy <- get.par('logBmsy', rep, exp=TRUE)
#' summary(rep)
#' plot(rep)
#' @import TMB
fit.spict <- function(inp, verbose=TRUE, dbg=0){
    rep <- NULL
    # Check input list
    inp <- check.inp(inp, verbose = verbose)
    datin <- make.datin(inp, dbg)
    pl <- inp$parlist
    tic <- Sys.time()
    # Cycle through phases
    for (i in 1:inp$nphases){
        if (inp$nphases > 1) cat(paste('Estimating - phase', i, '\n'))
        # Create TMB object
        obj <- make.obj(datin, pl, inp, phase=i)
        if (dbg<1){
            # Do estimation
            if (inp$optimiser == 'nlminb'){
                opt <- try(nlminb(obj$par, obj$fn, obj$gr, control=inp$optimiser.control))
                if (class(opt)!='try-error'){
                    pl <- obj$env$parList(opt$par)
                }
            }
            if (inp$optimiser == 'optim'){
                if (inp$optim.method == 'SANN'){
                    tmpgr <- obj$gr
                    obj$gr <- NULL # SANN doesn't use gradient
                }
                opt <- try(optim(par=obj$par, fn=obj$fn, gr=obj$gr, method=inp$optim.method,
                                 control=inp$optimiser.control))
                if (inp$optim.method == 'SANN' & class(opt) != 'try-error'){
                    obj$gr <- tmpgr # Restore gradient function
                    cat('SANN optimisation done, switching to BFGS to refine...\n')
                    opt <- try(optim(par=obj$par, fn=obj$fn, gr=obj$gr, method='BFGS',
                                 control=inp$optimiser.control))
                    #obj$fn(opt$par)
                    #obj$gr(opt$par)
                }
                if (class(opt) != 'try-error'){
                    opt$objective <- opt$value
                    pl <- obj$env$parList(opt$par)
                }

            }
        }
    }
    if (dbg<1){
        optfailflag <- class(opt)=='try-error'
        sdfailflag <- FALSE
        if (optfailflag){ # Optimisation failed
            cat('obj$par:\n')
            print(obj$par)
            cat('obj$fn:\n')
            print(obj$fn())
            cat('obj$gr:\n')
            print(obj$gr())
            stop('Could not fit model. Error msg:', opt)
        } else {
            if (inp$do.sd.report){
                # Calculate SD report
                rep <- try(TMB::sdreport(obj,
                                         getJointPrecision=inp$getJointPrecision,
                                         bias.correct=inp$bias.correct,
                                         bias.correct.control=inp$bias.correct.control,
                                         getReportCovariance = inp$getReportCovariance))
                sdfailflag <- class(rep) == 'try-error'
                if (sdfailflag){
                    warning('Could not calculate sdreport.\n')
                    rep <- NULL
                }
            }
            if (is.null(rep)){ # If sdreport failed or was not calculated
                rep <- list()
                if (sdfailflag){
                    rep <- list()
                    rep$sderr <- 1
                    rep$par.fixed <- opt$par
                    rep$cov.fixed <- matrix(NA, length(opt$par), length(opt$par))
                }
            }
            rep$inp <- inp
            rep$obj <- obj
            rep$opt <- opt
            rep$opt$gr <- rep$obj$gr(rep$opt$par)
            rep$pl <- obj$env$parList(opt$par)
            obj$fn()
            rep$Cp <- obj$report()$Cp
            rep$report <- obj$report()
            if (!sdfailflag & inp$reportall){
                #  - Calculate Prager's statistics -
                #rep <- calc.prager.stats(rep)
                # - Built-in OSAR -
                if (!inp$osar.method == 'none'){
                    reposar <- try(calc.osa.resid(rep))
                    if (class(reposar) != 'try-error'){
                        rep <- reposar
                    }
                }
            }
        }
    }
    toc <- Sys.time()
    if (!is.null(rep)){
        rep$computing.time <- as.numeric(toc - tic)
        class(rep) <- "spictcls"
    }
    return(rep)
}



calc.prager.stats <- function(rep){
    if (!'stats' %in% names(rep)){
        rep$stats <- list()
    }
    K <- get.par('logK', rep, exp=TRUE)[2]
    Bests <- get.par('logB', rep, exp=TRUE)[rep$inp$indest, 2]
    Bmsy <- get.par('logBmsy', rep, exp=TRUE)[2]
    if (!any(is.na(Bests)) & !is.na(Bmsy)){
        Bdiff <- Bmsy - Bests
        # Prager's nearness
        if (any(diff(sign(Bdiff))!=0)){
            rep$stats$nearness <- 1
        } else {
            rep$stats$nearness <- 1 - min(abs(Bdiff))/Bmsy
        }
        # Prager's coverage
        rep$stats$coverage <- min(c(2, (min(c(K, max(Bests))) - min(Bests))/Bmsy))
    }
    return(rep)
}

#' @name make.datin
#' @title Create data list used as input to TMB::MakeADFun.
#' @param inp List of input variables as output by check.inp.
#' @param dbg Debugging option. Will print out runtime information useful for debugging if set to 1.
#' @return List to be used as data input to TMB::MakeADFun.
#' @export
make.datin <- function(inp, dbg=0){
    datin <- list(reportall=as.numeric(inp$reportall),
                  reportRel=as.numeric(inp$reportRel),
                  dt=inp$dt,
                  dtpredcinds=inp$dtpredcinds,
                  dtpredcnsteps=inp$dtpredcnsteps,
                  dtprediind=inp$dtprediind,
                  dtpredeinds=inp$dtpredeinds,
                  dtpredensteps=inp$dtpredensteps,
                  indlastobs=inp$indlastobs,
                  obssrt=inp$obssrt,
                  stdevfacc=inp$stdevfacC,
                  stdevfaci=inp$stdevfacIin,
                  stdevface=inp$stdevfacE,
                  isc=inp$isc,
                  isi=inp$isi,
                  ise=inp$ise,
                  nobsC=inp$nobsC,
                  nobsI=sum(inp$nobsI),
                  nobsE=inp$nobsE,
                  ic=inp$ic,
                  nc=inp$nc,
                  ie=inp$ie,
                  ne=inp$ne,
                  ii=inp$iiin,
                  iq=inp$iqin,
                  isdi=inp$isdiin,
                  ir=inp$ir,
                  isdf=inp$isdf,
                  logmcov=inp$logmcovariatein,
                  seasons=inp$seasons,
                  seasonindex=inp$seasonindex,
                  nseasons=inp$nseasons,
                  seasonindex2=inp$seasonindex2,
                  splinemat=inp$splinemat,
                  splinematfine=inp$splinematfine,
                  omega=inp$omega,
                  seasontype=inp$seasontype,
                  efforttype=inp$efforttype,
                  timevaryinggrowth=as.numeric(inp$timevaryinggrowth),
                  logmcovflag=as.numeric(inp$logmcovflag),
                  ffacvec=inp$ffacvec,
                  fconvec=inp$fconvec,
                  indpred=inp$indpred,
                  robflagc=inp$robflagc,
                  robflagi=inp$robflagi,
                  robflage=inp$robflage,
                  stochmsy=ifelse(inp$msytype=='s', 1, 0),
                  stabilise=inp$stabilise,
                  MSYregime=inp$MSYregime,
                  iuse=as.numeric(inp$iuse),

                  priorn=inp$priors$logn,
                  priorngamma=inp$priors$logngamma,
                  priorr=inp$priors$logr,
                  priorK=inp$priors$logK,
                  priorm=inp$priors$logm,
                  priormu=inp$priors$mu,
                  priorq=inp$matrixpriors$logq,
                  priorqf=inp$priors$logqf,
                  priorbkfrac=inp$priors$logbkfrac,
                  priorsdb=inp$priors$logsdb,
                  priorsdm=inp$priors$logsdm,
                  priorsdf=inp$priors$logsdf,
                  priorsdi=inp$matrixpriors$logsdi,
                  priorsde=inp$priors$logsde,
                  priorsdc=inp$priors$logsdc,
                  prioralpha=inp$priors$logalpha,
                  priorbeta=inp$priors$logbeta,
                  priorpsi=inp$priors$logpsi,
                  priorB=inp$priors$logB,
                  priorF=inp$priors$logF,
                  priorBBmsy=inp$priors$logBBmsy,
                  priorFFmsy=inp$priors$logFFmsy,
                  priorBmsyB0=inp$priors$BmsyB0,

                  simple=inp$simple,
                  reportmode=inp$reportmode,
                  simRandomEffects=inp$sim.random.effects,
                  dbg=dbg)
    return(datin)
}


#' @name make.obj
#' @title Create TMB obj using TMB::MakeADFun and squelch screen printing.
#' @param datin Data list.
#' @param pl Parameter list.
#' @param inp List of input variables as output by check.inp.
#' @param phase Estimation phase, integer.
#' @return List to be used as data input to TMB.
#' @export
#' @import TMB
make.obj <- function(datin, pl, inp, phase=1){
    obj <- TMB::MakeADFun(data=datin, parameters=pl, random=inp$RE, DLL=inp$scriptname,
                          hessian=TRUE, map=inp$map[[phase]])
    TMB:::config(trace.optimize=0, DLL=inp$scriptname)
    verbose <- FALSE
    obj$env$tracemgc <- verbose
    obj$env$inner.control$trace <- verbose
    obj$env$silent <- ! verbose
    obj$fn(obj$par)
    return(obj)
}
