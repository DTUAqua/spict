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
#' @details This hindcasting analysis fits spict to the data set while
#'     sequentially exluding the index observations in the last years.
#' @param rep A valid result from fit.spict.
#' @param nyears Number of years of data (catch and effort) to remove (this is
#'     also the total number of model runs).
#' @param reduce.output.size logical, if \code{TRUE} (default) hindcasting is
#'     run with \code{getReportCovariance} and \code{getJointPrecision} set as
#'     \code{FALSE}
#' @param mc.cores Number of cores for \code{parallel::mclapply} function. By
#'     default 1.
#'
#' @return A spictcls list with the added element 'hindcast' containing the
#'     results of the hindcasting analysis. Use plotspict.hindcast() to plot these
#'     results.
#'
#' @importFrom parallel mclapply
#'
#' @examples
#' data(pol)
#' inp <- pol$albacore
#' rep <- fit.spict(inp)
#' rep <- hindcast(rep, nyears = 5)
#' plotspict.hindcast(rep)
#' @export
hindcast <- function(rep, nyears = 7, reduce.output.size = TRUE, mc.cores = 1){

    if (!"spictcls" %in% class(rep)) stop("This function only works with a fitted spict object (class 'spictcls'). Please run `fit.spict` first.")
    if (rep$opt$convergence != 0) stop("The fitted object did not converged.")

    inpin <- rep$inp
    lastyears <- max(inpin$timerangeObs) - 1:nyears

    inpall <- list()
    for (i in 1:nyears) {
        inpall[[i]] <- inpin
        ## Catch
        indsC <- which(inpin$timeC <= lastyears[i]) ## inpin$timeC[inpin$nobsC] - i + 1)
        inpall[[i]]$obsC <- inpin$obsC[indsC]
        inpall[[i]]$timeC <- inpin$timeC[indsC]
        inpall[[i]]$stdevfacC <- inpin$stdevfacC[indsC]
        inpall[[i]]$dtc <- inpin$dtc[indsC]
        ## Effort
        indsE <- which(inpin$timeE <= lastyears[i]) ## inpin$timeE[inpin$nobsE] - i + 1)
        inpall[[i]]$obsE <- inpin$obsE[indsE]
        inpall[[i]]$timeE <- inpin$timeE[indsE]
        inpall[[i]]$stdevfacE <- inpin$stdevfacE[indsE]
        inpall[[i]]$dte <- inpin$dte[indsE]
        ## Index
        inpall[[i]]$obsI <- list()
        inpall[[i]]$timeI <- list()
        for (j in seq_len(inpin$nindex)) {
            indsI <- which(inpin$timeI[[j]] <= lastyears[i] + 1) ## inpin$timeI[[j]][inpin$nobsI[j]] - i + 1)
            inpall[[i]]$obsI[[j]] <- inpin$obsI[[j]][indsI]
            inpall[[i]]$timeI[[j]] <- inpin$timeI[[j]][indsI]
            inpall[[i]]$stdevfacI[[j]] <- inpin$stdevfacI[[j]][indsI]
        }
        inpall[[i]]$getReportCovariance <- !reduce.output.size
        inpall[[i]]$getJointPrecision <- !reduce.output.size
        ## Check
        inpall[[i]] <- check.inp(inpall[[i]])
        ## Flag out last index observations
        inpall[[i]]$iuse[cumsum(inpall[[i]]$nobsI)] <- !ceiling(sapply(inpall[[i]]$timeI, max)) >= lastyears[i]
    }

    asd <- try(parallel::mclapply(inpall, fit.spict, mc.cores = mc.cores))
    if (class(asd) == "try-error") {
        rep$hindcast <- lapply(inpall, fit.spict)
    } else {
        rep$hindcast <- asd
    }
    ## Add the base run into the retro list
    baserun <- rep
    baserun$hindcast <- NULL
    baserun$cov <- NA
    baserun$man <- NULL
    rep$hindcast <- c(list(baserun), rep$hindcast)
    return(rep)
}







#' @name calc.mase
#' @title Calculate Mean Absolute Scaled Error (MASE)
#' @details A function that calculates the Mean Absolute Scaled Error (MASE) for
#'     the indices.
#' @param rep A valid result from fit.spict
#' @param verbose Should detailed outputs be provided (default: TRUE).
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
calc.mase <- function(rep, verbose = TRUE) {

    if (!"spictcls" %in% class(rep)) stop("This function only works with a fitted spict object (class 'spictcls'). Please run `fit.spict` first.")
    if (!"hindcast" %in% names(rep)) stop("No results of the retro function found. Please run hindcasting using the `hindcast` function.")


    ## Hindcast
    hindcast <- rep$hindcast
    nhindcast <- length(hindcast) - 1
    runs <- c(0, seq_len(nhindcast-1))
    runs <- seq_len(nhindcast)

    ## Time
    timeRangeSurv <- range(unlist(rep$inp$timeI))
    survYears <- seq(floor(timeRangeSurv[1]), floor(timeRangeSurv[2]),1)
    hindcastyears <- survYears[(length(survYears)-nhindcast+1):length(survYears)]
    revhindcastyears <- rev(hindcastyears)

    ## Indices
    nind <- length(rep$inp$timeI)
    indices <- seq_len(nind)
    ## Check that surveys overlap with hindcast years
    valid <- NULL
    for(i in 1:nind){
        valid <- c(valid, ifelse(any(hindcastyears %in% floor(rep$inp$timeI[[i]])), TRUE, FALSE))
    }
    indices <- indices[valid]
    nind <- length(indices)

    ## Error if no valid survey
    if(nind < 1) stop("The survey(s) do(es) not overlap with the hindcast years! Cannot perform hindcasting cross-validation based on the hindcasted runs!")

    ## convergence
    conv <- sapply(hindcast[-1], function(x) x$opt$convergence)
    conv <- ifelse(conv == 0, TRUE, FALSE)
    conv <- rev(conv)

    ## for extracting logIpred
    nobsI <- lapply(hindcast, function(x) x$inp$nobsI)
    indI <- list()
    for(i in 1:length(nobsI)){
        indI[[i]] <- list()
        for(j in 1:length(nobsI[[i]])){
            if(j == 1){
                indI[[i]][[j]] <- seq_len(nobsI[[i]][j])
            }else{
                indI[[i]][[j]] <- seq(nobsI[[i]][j-1]+1, nobsI[[i]][j-1] + nobsI[[i]][j], 1)
            }
        }
    }

    mase <- NULL
    for(i in 1:nind){
        if(valid[i]){
            ## Extract survey obs and preds from all hindcast runs
            dat <- as.data.frame(
                do.call(rbind,
                        lapply(1:(nhindcast+1),
                               function(x)
                                   cbind(obs = hindcast[[x]]$inp$obsI[[i]],
                                         pred = unname(get.par("logIpred",
                                                               hindcast[[x]], exp = TRUE)[indI[[x]][[i]],2]),
                                         lc = unname(get.par("logIpred",
                                                             hindcast[[x]], exp = TRUE)[indI[[x]][[i]],1]),
                                         uc = unname(get.par("logIpred",
                                                             hindcast[[x]], exp = TRUE)[indI[[x]][[i]],3]),
                                         year = floor(hindcast[[x]]$inp$timeI[[i]]),
                                         run = x-1))))

            ## Variables
            years <- sort(unique(dat$year))
            ind <- dat$run == min(dat$run)
            py <- dat$year[ind]
            obs <- dat$obs[ind]
            pred <- dat$pred[ind]
            lc <- dat$lc[ind]
            uc <- dat$uc[ind]

            ## Naive diff
            ind <- which(hindcastyears %in% years & conv)
            npe <- length(ind)
            ## Missing values on the bounds okay, but if in middle of vector -> unequal spacing! -> Warning
            if(verbose && !all(seq(min(ind), max(ind), 1) %in% ind))
                cat('Warning: Unequal spacing of naive predictions residuals may influence the interpretation of MASE. \n')
            obsi <- rep(NA, length(hindcastyears))
            obsi[hindcastyears %in% years] <- dat$obs[dat$run == min(dat$run)][years %in% hindcastyears]
            obsi <- obsi[conv]
            isNAobs <- is.na(obsi)
            naive <- log(obsi[!isNAobs][-length(obsi[!isNAobs])]) - log(obsi[!isNAobs][-1])
            isNAnaive <- is.na(naive)

            ## Predictions
            pred <- NULL
            for(j in 1:nhindcast){
                if(revhindcastyears[j] %in% dat$year && rev(conv)[j]){
                    x <- min(py):max(hindcastyears)
                    x <- x[1:(length(x) - runs[j] + 1)]
                    x <- x[x %in% years]
                    n <- length(x)
                    y <- dat[dat$run == runs[j] & dat$year %in% x,]$pred
                    pred <- c(pred, log(y[n]) - log(obs[n]))
                }
            }

            ## MASE
            masei <- mean(abs(pred))/mean(abs(naive))
            mase <- rbind(mase, data.frame(Index = i,
                                           MASE = masei,
                                           nruns = npe))
        }else{
            cat(paste0("\n","No observations in evaluation years to compute prediction residuals for Index ",i),"\n")
            mase <- rbind(mase, data.frame(Index = i,
                                           MASE = NA,
                                           nruns = 0))
        }
    }

    nnotconv <- sum(!conv)
    if(nnotconv > 0){
        message("Excluded ", nnotconv, " hindcasted runs that ",
                if (nnotconv == 1) "was" else "were" , " not converged: ",
                paste(which(!rev(conv)), collapse = ", "))
    }

    return(mase)
}
