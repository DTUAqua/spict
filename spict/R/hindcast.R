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


#' @name hindcast
#' @title Conduct hindcasting analysis
#'
#' @param rep rep Result of \code{\link{fit.spict}}.
#' @param npeels Number of years/seasons (dependent on dtc) of data (catch and
#'     effort) to remove (this is also the total number of model runs).
#' @param reduce.output.size logical, if \code{TRUE} (default) hindcasting is
#'     run with \code{getReportCovariance} and \code{getJointPrecision} set as
#'     \code{FALSE}
#' @param mc.cores Number of cores for \code{parallel::mclapply} function. By
#'     default 1.
#' @param peel.dtc Peel according to catch seasons (dtc) rather than years? It
#'     only differs if the data includes seasonal catches (Default: FALSE)
#'
#' @details This method creates a number of subsets (or peels) specified with
#'     the argument \code{npeels} by sequentially omitting all observations from
#'     the most recent time step (by default corresponding to a year). Then,
#'     spict is fitted to each subset while excluding all index observations in
#'     the most recent year with data. The resulting fitted spict objects are
#'     attached to the base spict object as a list element labeled
#'     \code{hindcast} and resulted by this method.
#'
#'     The prediction of the excluded index observations relative to the actual
#'     excluded index can then be compared to a 'naive' prediction using the
#'     preceding index observation. This allows to estimate the Mean Absolute
#'     Scaled Error (MASE) for each index \code{\link{calc.mase}}.
#'
#'     The predicted indices can be visualised with the function
#'     \code{\link{plotspict.hindcast}}.
#'
#' @return A \code{spictcls} list with the added element \code{hindcast}
#'     containing the results of the hindcasting analysis. Use
#'     \code{\link{plotspict.hindcast}} to plot these results.
#'
#' @importFrom parallel mclapply
#'
#' @examples
#' data(pol)
#' inp <- pol$albacore
#' rep <- fit.spict(inp)
#' rep <- hindcast(rep, npeels = 5)
#' plotspict.hindcast(rep)
#' @export
#'
#' @references
#' Carvalho, F., Winker, H., Courtney, D., Kapur, M., Kell, L., Cardinale, M.,
#' Schirripa, M., Kitakado, T., Yemane, D., Piner, K.R. Maunder, M.N., Taylor,
#' I., Wetzel, C.R., Doering, K., Johnsonm, K.F., Methot, R. D. (2021). A
#' cookbook for using model diagnostics in integrated stock assessments.
#' Fisheries Research, 240, 105959.
#'
#' Kell, L. T., Kimoto, A., & Kitakado, T. (2016). Evaluation of the prediction
#' skill of stock assessment using hindcasting. Fisheries research, 183,
#' 119-127.
#'
#' Kell, L. T., Sharma, R., Kitakado, T., Winker, H., Mosqueira, I., Cardinale,
#' M., & Fu, D. (2021). Validation of stock assessment methods: is it me or my
#' model talking?. ICES Journal of Marine Science, 78(6), 2244-2255.
#'
#' Winker, H., Carvalho, F., & Kapur, M. (2018). JABBA: just another Bayesian
#' biomass assessment. Fisheries Research, 204, 275-288.
#'
hindcast <- function(rep, npeels = 7, reduce.output.size = TRUE, mc.cores = 1, peel.dtc = FALSE){

    check.rep(rep)
    if(rep$opt$convergence != 0) stop("The fitted object did not converge.")

    inpin <- rep$inp

    ## Accounting for the cast that the last index obs is at end of catch
    ## interval: In this case the index obs is interpreted to belong to
    ## lastyear rather than end of lastyear-1, while catch obs belongs to
    ## lastyear-1
    if(max(unlist(inpin$timeI)) == inpin$lastCatchObs){
        peeling <- 0:(npeels-1)
    }else{
        peeling <- 1:npeels
    }
    ## Accounting for sub-annual catch data: This allows to remove seasons and
    ## thus indices subsequently rather than all at once.
    if(peel.dtc){
        tmp <- seq(floor(inpin$timerange[1]), ceiling(inpin$timerange[2] + 100), tail(inpin$dtc,1))
        if(inpin$timerangeObs[2] == inpin$lastCatchObs){
            hindcastTimes <- tmp[as.integer(cut(inpin$timerangeObs[2],tmp,right = FALSE))] -
                c(0,cumsum(rev(inpin$dtc)))[peeling+1]
        }else{
            hindcastTimes <- tmp[as.integer(cut(inpin$timerangeObs[2],tmp,right = FALSE)) + 1] -
                c(0,cumsum(rev(inpin$dtc)))[peeling+1]
        }
        extraTime <- rev(inpin$dtc)[1:npeels]
    }else{
        hindcastTimes <- ceiling(inpin$timerangeObs[2]) - peeling
        extraTime <- rep(1, npeels)
    }

    inpall <- list()
    for (i in 1:npeels) {
        inpall[[i]] <- inpin
        ## Catch
        indsC <- which(inpin$timeC + inpin$dtc <= hindcastTimes[i] + extraTime[i])
        inpall[[i]]$obsC <- inpin$obsC[indsC]
        inpall[[i]]$timeC <- inpin$timeC[indsC]
        inpall[[i]]$stdevfacC <- inpin$stdevfacC[indsC]
        inpall[[i]]$dtc <- inpin$dtc[indsC]
        ## Effort
        indsE <- which(inpin$timeE + inpin$dte <= hindcastTimes[i] + extraTime[i])
        inpall[[i]]$obsE <- inpin$obsE[indsE]
        inpall[[i]]$timeE <- inpin$timeE[indsE]
        inpall[[i]]$stdevfacE <- inpin$stdevfacE[indsE]
        inpall[[i]]$dte <- inpin$dte[indsE]
        ## Index
        inpall[[i]]$obsI <- list()
        inpall[[i]]$timeI <- list()
        for (j in seq_len(inpin$nindex)) {
            if(inpin$timerangeObs[2] > inpin$lastCatchObs){
                indsI <- which(inpin$timeI[[j]] <= hindcastTimes[i] + extraTime[i])
            }else{
                indsI <- which(inpin$timeI[[j]] < hindcastTimes[i] + extraTime[i])
            }
            inpall[[i]]$obsI[[j]] <- inpin$obsI[[j]][indsI]
            inpall[[i]]$timeI[[j]] <- inpin$timeI[[j]][indsI]
            inpall[[i]]$stdevfacI[[j]] <- inpin$stdevfacI[[j]][indsI]
        }
        inpall[[i]]$getReportCovariance <- !reduce.output.size
        inpall[[i]]$getJointPrecision <- !reduce.output.size
        ## Check
        inpall[[i]] <- check.inp(inpall[[i]])
        ## Flag out last index observations
        inpall[[i]]$iuse[cumsum(inpall[[i]]$nobsI)] <-
            sapply(inpall[[i]]$timeI, function(x) !(max(x) >= hindcastTimes[i]))
    }

    asd <- try(parallel::mclapply(inpall, fit.spict, mc.cores = mc.cores))
    if(inherits(asd, "try-error")){
        rep$hindcast <- lapply(inpall, fit.spict)
    }else{
        rep$hindcast <- asd
    }

    ## Add the base run into the retro list
    rep$hindcast <- c(list(prune.baserun(rep)), rep$hindcast)
    return(rep)
}



#' @name extract.hindcast.info
#' @title Extract hindcast info from a fitted spict object
#'
#' @param rep Result of \code{\link{fit.spict}} that contains hindcasted runs
#'     added by \code{\link{hindcast}}.
#' @param CI Confidence intervals to be calculated, e.g. 0.9 for the 90%
#'     confidence intervals. By default (CI = 0.95), the 95% confidence
#'     intervals are estimated.
#' @param verbose Should detailed outputs be provided (default: TRUE).
#'
#' @details Note that a difference in the timing of the index observations of
#'     less than a month are considered acceptable for the estimation of the
#'     naive prediction residuals and no warning is printed. If the variability
#'     exceeds a month, the predictions are still calculated, but a warning is
#'     printed.
#'
#' @return A list with hindcast information.
#'
extract.hindcast.info <- function(rep, CI = 0.95, verbose = TRUE) {

    ## Checks
    check.rep(rep)
    if (!"hindcast" %in% names(rep)) stop("No results of the hindcast function found. Please run hindcasting using the `hindcast` function.")

    inpin <- rep$inp
    hindcast <- rep$hindcast
    npeels <- length(hindcast) - 1
    peels <- 1:npeels
    ## Account for interpretation of cont time
    if(max(unlist(inpin$timeI)) == inpin$lastCatchObs){
        peeling <- 0:(npeels-1)
    }else{
        peeling <- 1:npeels
    }
    ## If latest obs is midyear index than first diff could be < 1
    peel.dtc <- ifelse(length(which(abs(diff(
        sapply(hindcast[-1], function(x) max(c(x$inp$timeC + x$inp$dtc,
                                               x$inp$timeE + x$inp$dte))))) < 1)) >= 2, 1, 0)
    if(peel.dtc){
        tmp <- seq(floor(inpin$timerange[1]), ceiling(inpin$timerange[2] + 100), tail(inpin$dtc,1))
        hindcastTimes <- tmp[as.integer(cut(inpin$timerangeObs[2],tmp,right = FALSE)) + 1] -
            cumsum(rev(inpin$dtc))[peeling]
    }else{
        hindcastTimes <- ceiling(inpin$timerangeObs[2]) - peeling
    }

    ## Indices
    nind <- length(rep$inp$timeI)
    ## Check that surveys overlap with hindcast years
    valid <- NULL
    for(i in 1:nind){
        valid <- c(valid, ifelse(max(rep$inp$timeI[[i]]) >= min(hindcastTimes), TRUE, FALSE))
    }
    nind.val <- length(which(valid))

    ## Error if no valid survey
    if(nind.val < 1) stop("Not a single index overlaps with the peels chosen for hindcasting. Thus, cannot calcute or plot metrics based on hindcasting cross-validation. Choose more peels and re-run the hindcasting.")

    ## convergence
    conv <- sapply(hindcast[-1], function(x) x$opt$convergence)
    conv <- rev(ifelse(conv == 0, TRUE, FALSE))

    ## for extracting logIpred
    nobsI <- lapply(hindcast, function(x) x$inp$nobsI)
    indI <- list()
    for(i in 1:length(nobsI)){
        indI[[i]] <- list()
        for(j in 1:length(nobsI[[i]])){
            if(j == 1){
                indI[[i]][[j]] <- seq_len(nobsI[[i]][j])
            }else{
                indI[[i]][[j]] <- seq(cumsum(nobsI[[i]])[j-1]+1, cumsum(nobsI[[i]])[j], 1)
            }
        }
    }

    mase <- NULL
    res.index <- vector("list", nind.val)
    ival <- 0
    for(i in 1:nind){
        if(valid[i]){
            ival <- ival + 1
            ## observations and predictions
            dat.list <- vector("list", npeels + 1)
            for(j in 1:(npeels+1)){
                obs <- log(hindcast[[j]]$inp$obsI[[i]])
                iuse <- hindcast[[j]]$inp$iuse[indI[[j]][[i]]]
                pred <- unname(get.par("logIpred",
                                       hindcast[[j]], exp = FALSE)[indI[[j]][[i]],1:3])
                colnames(pred) <- c("lc","pred","uc")
                time <- hindcast[[j]]$inp$timeI[[i]]
                peel <- rep(j-1, nobsI[[j]][i])
                dat.list[[j]] <- data.frame(obs, iuse, pred, time, peel)
            }
            dat <- do.call(rbind, dat.list)
            dat0 <- dat[dat$peel == min(dat$peel),]

            ## Predictions
            pred <- data.frame(matrix(NA,npeels,4))
            for(j in 1:npeels){
                dati <- dat[dat$peel == peels[j],]
                ## only if index was peeled in that peel
                if(!tail(dati$iuse,1)){
                    pred[j,1] <- tail(dati$time,1)
                    pred[j,2] <- tail(dati$pred,1) - tail(dati$obs,1)
                    pred[j,3] <- tail(dat0$obs[dat0$time < tail(dati$time,1)],1) - tail(dati$obs,1)
                    pred[j,4] <- tail(dati$time,1) - tail(dat0$time[dat0$time < tail(dati$time,1)],1)
                }
            }
            colnames(pred) <- c("time","pred","naive","spacing")
            ## remove non-converged peels
            pred <- pred[rev(conv),]

            if(nrow(pred) < 1){
                stop("Not enough converged naive predictions available. Increase the number of peels and re-run hindcasting!")
            }
            if(verbose && !all(abs(pred$spacing[!is.na(pred$spacing)] - 1) * 12 <= 1)){
                cat(paste0("Warning: The observations of Index ", i, " are not constant in time and vary by more than a month, because of variable timing of index observations, missing observations in intermediate years, or non-converged peels. This might affect the interpretation of MASE. \n"))
            }

            ## MASE
            mase <- rbind(mase, data.frame(Index = i,
                                           MASE = mean(abs(pred$pred), na.rm = TRUE) /
                                               mean(abs(pred$naive), na.rm = TRUE),
                                           npeels = nrow(pred)))
            res.index[[ival]] <- list(dat = dat, pred = pred)
        }else{
            cat(paste0("No observations in evaluation years to compute prediction residuals for Index ",i),"\n")
        }
    }

    ## Message about non-convergence
    nnotconv <- sum(!conv)
    if(verbose && nnotconv > 0){
        message(nnotconv, " peels corresponding to the year(s) ",
                paste(hindcastTimes[which(!rev(conv))], collapse = ", "),
                " did not converge.")
    }

    ## Output
    res <- list()
    res$index <- res.index
    res$mase <- mase

    return(res)
}





#' @name calc.mase
#' @title Calculate Mean Absolute Scaled Error (MASE)
#'
#' @param rep Result of \code{\link{fit.spict}} that contains hindcasted runs added by \code{\link{hindcast}}.
#' @param verbose Should detailed outputs be provided (default: TRUE).
#'
#' @details This function calculates the Mean Absolute Scaled Error (MASE) for
#'     each index time series by hindcasting. Thus, the application of this
#'     method requires a fitted spict object with the results of the hindcasting
#'     analysis \code{\link{hindcast}}.
#'
#'     The smaller MASE, the higher the predictive power of the spict model
#'     regarding the prediction of index observations. In contrast, a MASE above
#'     1 suggests that the naive prediction of the index observations assuming
#'     the preceding index observations have a higher predictive power than the
#'     spict model. Note, however, that the absolute MASE value depends on
#'     multiple factors such as the number of peels, assumed priors, etc.
#'
#'     Note that a difference in the timing of the index observations of less
#'     than a month are considered acceptable for the estimation of the naive
#'     prediction residuals and no warning is printed. If the variability
#'     exceeds a month, the predictions are still calculated, but a warning is
#'     printed.
#'
#' @return A data frame with estimate the MASE and the number of runs used for
#'     the estimation for each index.
#'
#' @examples
#' data(pol)
#' inp <- pol$albacore
#' rep <- fit.spict(inp)
#' rep <- hindcast(rep)
#' calc.mase(rep)
#' @export
#'
#' @references
#' Carvalho, F., Winker, H., Courtney, D., Kapur, M., Kell, L., Cardinale, M.,
#' Schirripa, M., Kitakado, T., Yemane, D., Piner, K.R. Maunder, M.N., Taylor,
#' I., Wetzel, C.R., Doering, K., Johnsonm, K.F., Methot, R. D. (2021). A
#' cookbook for using model diagnostics in integrated stock assessments.
#' Fisheries Research, 240, 105959.
#'
#' Kell, L. T., Kimoto, A., & Kitakado, T. (2016). Evaluation of the prediction
#' skill of stock assessment using hindcasting. Fisheries research, 183,
#' 119-127.
#'
#' Kell, L. T., Sharma, R., Kitakado, T., Winker, H., Mosqueira, I., Cardinale,
#' M., & Fu, D. (2021). Validation of stock assessment methods: is it me or my
#' model talking?. ICES Journal of Marine Science, 78(6), 2244-2255.
#'
#' Winker, H., Carvalho, F., & Kapur, M. (2018). JABBA: just another Bayesian
#' biomass assessment. Fisheries Research, 204, 275-288.
#'
calc.mase <- function(rep, verbose = TRUE){
    hcInfo <- extract.hindcast.info(rep, verbose = verbose)
    return(hcInfo$mase)
}
