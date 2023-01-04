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


#' @name retro
#' @title Conduct retrospective analysis
#' @details A retrospective analysis consists of estimating the model with later
#'     data points removed sequentially one year at a time.
#' @param rep A valid result from fit.spict.
#' @param nretroyear Number of years of data to remove (this is also the total
#'     number of model runs).
#' @param reduce_output_size logical, if \code{TRUE} (default) the retro is run
#'     with \code{getReportCovariance} and \code{getJointPrecision} set as
#'     \code{FALSE}
#' @param mc.cores Number of cores for \code{parallel::mclapply} function. By
#'     default 1.
#'
#' @importFrom parallel mclapply
#'
#' @return A rep list with the added key retro containing the results of the
#'     retrospective analysis. Use plotspict.retro() to plot these results.
#' @examples
#' data(pol)
#' inp <- pol$albacore
#' rep <- fit.spict(inp)
#' rep <- retro(rep, nretroyear=6)
#' plotspict.retro(rep)
#' @export
retro <- function(rep, nretroyear=5, reduce_output_size = TRUE, mc.cores = 1){
    check.rep(rep)
    if(rep$opt$convergence != 0) stop("The fitted object did not converge.")
    inp1 <- rep$inp

    lastyears <- max(inp1$timerangeObs) - 1:nretroyear

    inpall <- list()
    for (i in 1:nretroyear) {
        inpall[[i]] <- inp1
        ## catch
        indsC <- which(inp1$timeC < lastyears[i]) ## which(inp1$timeC <= inp1$timeC[inp1$nobsC] - i)
        inpall[[i]]$obsC <- inp1$obsC[indsC]
        inpall[[i]]$timeC <- inp1$timeC[indsC]
        inpall[[i]]$stdevfacC <- inp1$stdevfacC[indsC]
        inpall[[i]]$dtc <- inp1$dtc[indsC]
        ## effort
        indsE <- which(inp1$timeE < lastyears[i]) ## which(inp1$timeE <= inp1$timeE[inp1$nobsE] - i)
        inpall[[i]]$obsE <- inp1$obsE[indsE]
        inpall[[i]]$timeE <- inp1$timeE[indsE]
        inpall[[i]]$stdevfacE <- inp1$stdevfacE[indsE]
        inpall[[i]]$dte <- inp1$dte[indsE]
        ## indices
        inpall[[i]]$obsI <- list()
        inpall[[i]]$timeI <- list()
        for (j in seq_len(inp1$nindex)) {
            indsI <- which(inp1$timeI[[j]] < lastyears[i] + 1) ## which(inp1$timeI[[j]] <= inp1$timeI[[j]][inp1$nobsI[j]] - i)
            inpall[[i]]$obsI[[j]] <- inp1$obsI[[j]][indsI]
            inpall[[i]]$timeI[[j]] <- inp1$timeI[[j]][indsI]
            inpall[[i]]$stdevfacI[[j]] <- inp1$stdevfacI[[j]][indsI]
        }
        inpall[[i]]$getReportCovariance <- !reduce_output_size
        inpall[[i]]$getJointPrecision <- !reduce_output_size
        inpall[[i]] <- check.inp(inpall[[i]])
    }
    retroruns <- try(parallel::mclapply(inpall, fit.spict,
                                        mc.cores = mc.cores))
    if (class(retroruns) == "try-error") {
        rep$retro <- lapply(inpall, fit.spict)
    }
    else {
        rep$retro <- retroruns
    }
    ## Add the base run into the retro list
    rep$retro <- c(list(prune.baserun(rep)), rep$retro)
    names(rep$retro) <- c("Baseline", paste("-", seq(nretroyear)))

    return(rep)
}


#' @name mohns_rho
#' @title Calculate Mohn's rho for different estimates
#' @details A function that calculates Mohn's rho for selected estimated quantities. The function allows the user to define the method of aggrgating from the subannual time steps (1/dteuler) into annual values; the default is to take the mean.
#' @param rep A valid result from fit.spict
#' @param what character vector specifying the quantities
#' @param annualfunc function used to convert subannual data into annual
#' @return A named vector with the Monh's rho value for each quantity.
#' @examples
#' data(pol)
#' inp <- pol$albacore
#' rep <- fit.spict(inp)
#' rep <- retro(rep, nretroyear = 4)
#' mohns_rho(rep)
#' @export
mohns_rho <- function(rep, what = c("FFmsy", "BBmsy"), annualfunc = mean) {
  if (!"spictcls" %in% class(rep)) stop("This function only works with a fitted spict object (class 'spictcls'). Please run `fit.spict` first.")
  if (!"retro" %in% names(rep)) stop("No results of the retro function found. Please run the retrospective analysis using the `retro` function.")
  if ("Ipred" %in% what && rep$inp$nindex > 1) warning("Mohn's rho will be calculated only for index 1.")

  getFullYearEstimates <- function(x, what = c("FFmsy", "BBmsy"), annualfunc = mean) {
    res <- lapply(what, function(ww) {
      if (ww == "Ipred" || ww == "Ipred1") {
        if (x$inp$nindex == 1) {
          get.par("logIpred", x, exp = TRUE)[, 2]
        } else {
          get.par("logIpred", x, exp = TRUE)[seq(length(x$inp$timeI[[1]])), 2]
        }
      } else if (startsWith(ww, prefix = "Ipred") || startsWith(ww, "logIpred")) {
        idx <- as.numeric(sub("log", "", sub("Ipred", "", ww)))
        if (idx > x$inp$nindex) stop("There are only ", x$inp$nindex, " indices in the model.")
        start <- sum(sapply(x$inp$timeI[[seq.int(idx - 1)]], length)) + 1
        len <- length(x$inp$timeI[[idx]])
        get.par("logIpred", x, exp = TRUE)[seq.int(start, start + len - 1), 2]
      } else {
        par <- paste0("log", sub("log", "", ww))
        if (!par %in% list.quantities(x)) stop(ww, " is not a reported quantity, use `list.quantities(rep)` to list all available ones.")
        indest <- x$inp$indest
        time <- x$inp$time[indest]
        res <- annual(time, get.par(par, x, exp = TRUE)[, 2], annualfunc)
        setNames(as.array(res$annvec), res$anntime)
      }
    })
    res <- do.call(cbind.data.frame, res)
    names(res) <- what
    setNames(do.call(cbind.data.frame, res), what)
  }
  ## Exclude not converged runs
  conv <- sapply(rep$retro, function(x) x$opt$convergence == 0)
  nnotconv <- sum(!conv)
  if (nnotconv > 0) {
    message("Exluded ", nnotconv, " retrospective runs that ",
            if (nnotconv == 1) "was" else "were" , " not converged: ",
            paste(which(!conv) - 1, collapse = ", "))
  }
  rep$retro <- rep$retro[conv]
  ## Adapted from fishfollower/SAM/stockassessment package
  ret <- lapply(rep$retro, getFullYearEstimates, what = what, annualfunc = annualfunc)
  ref <- ret[[1]]
  bias <- lapply(ret[-1], function(x) {
    y <- rownames(x)[nrow(x)]
    (x[rownames(x) == y, ] - ref[rownames(ref) == y, ]) / ref[rownames(ref) == y, ]
  })
  colMeans(do.call(rbind, bias))
}
