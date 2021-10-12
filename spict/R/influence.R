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


#' @name calc.influence
#' @title Calculates influence statistics of observations.
#' @details TBA
#' @param rep A valid result from fit.spict().
#' @param mc.cores Number of cores for \code{parallel::mclapply} function. By
#'     default 1.
#' @return A list equal to the input with the added key "infl" containing
#'     influence statistics.
#' @export
calc.influence <- function(rep, mc.cores = 1){
    #parnams <- c('logFmsy', 'MSY', 'logCp', 'logBlBmsy')
    if (!'osar' %in% names(rep)){
        stop('Need to calculate one-step-ahead residuals before calculating influence statistics. Use the calc.osa.resid() function.')
    }
    inp <- rep$inp
    #parnams <- c('logFmsy', 'MSY')
    parnams <- c('logm', 'logBmsy', 'logFmsy', 'logsdf', 'logsdb')
    parnamsnolog <- sub('log', '', parnams)
    np <- length(parnams)
    doexp <- rep(0, np)
    doexp[grep('log', parnams)] <- 1
    parmat <- matrix(0, np, 5)
    rownames(parmat) <- parnamsnolog
    for (k in 1:np){
        parmat[k, ] <- get.par(parnams[k], rep)
    }
    detcov <- unlist(determinant(rep$cov.fixed, logarithm=TRUE))[1]
    likval <- rep$opt$objective
    osarpvals <- get.osar.pvals(rep)
    rwnms <- paste0('C_', inp$timeC)
    for (i in 1:inp$nindex){
        rwnms <- c(rwnms, paste0('I', i, '_', inp$timeI[[i]]))
    }
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
        inp2$osar.trace <- FALSE
        if (c <= inp$nobsC){
            # Loop over catches
            inp2$obsC <- inp$obsC[-c]
            inp2$timeC <- inp$timeC[-c]
            inp2$dtc <- inp$dtc[-c]
        } else {
            # Loop over indices
            breaks <- inp$nobsC + cumsum(c(0, inp$nobsI[-inp$nindex]))
            j <- sum(c > breaks)
            i <- c - breaks[j]
            if (class(inp$obsI)=='numeric'){
                inp2$obsI <- inp$obsI[-i]
                inp2$timeI <- inp$timeI[-i]
            } else {
                inp2$obsI[[j]] <- inp$obsI[[j]][-i]
                inp2$timeI[[j]] <- inp$timeI[[j]][-i]
            }
        }
        rep2 <- fit.spict(inp2)
        #rep2 <- calc.osa.resid(rep2)
        # Calculate diagnostics
        np <- length(parnams)
        parmat <- matrix(0, np, 5)
        for (k in 1:np){
            parmat[k, ] <- get.par(parnams[k], rep2)
        }
        detcov <- unlist(determinant(rep2$cov.fixed, logarithm=TRUE))[1]
        dosarpvals <- get.osar.pvals(rep2)
        return(list(parmat=parmat, detcov=detcov, dosarpvals=dosarpvals))
    }
    # Calculate influence
    if (Sys.info()['sysname']!='Linux'){
    #if (class(partry)=='try-error'){
        single.inflfun <- function(c, inp, parnams, nobs){
            res <- inflfun(c, inp, parnams)
            ## Send progress update
            cat(sprintf("Progress: %.2f%%\r", c/nobs * 100))
            return(res)
        }

        res <- lapply(1:nobs, single.inflfun, inp, parnams, nobs)
    } else {
        #require(parallel)
        multi.inflfun <- function(c, inp, parnams, nobs, progfile){
            res <- inflfun(c, inp, parnams)
            ## Send progress update
            if (class(progfile)[1]=='fifo') writeBin(1/nobs, progfile)
            return(res)
        }
        ## Open fifo progress file
        progfile <- fifo(tempfile(), open="w+b", blocking=T)
        #if (inherits(fork(), "masterProcess")) {
        if (inherits(parallel:::mcfork(), "masterProcess")) {
            ## Child
            progress <- 0.0
            while (progress < 1 && !isIncomplete(progfile)) {
                msg <- readBin(progfile, "double")
                progress <- progress + as.numeric(msg)
                cat(sprintf("Progress (multicore): %.2f%%\r", progress * 100))
            }
            #exit()
        }
        res <- parallel::mclapply(1:nobs, multi.inflfun, inp, parnams, nobs, progfile, mc.cores=mc.cores)
        ## Close progress file
        close(progfile)
    }
    cat('\n')
    # Gather results
    for (i in 1:nobs){
        ddetcov[i] <- detcov - res[[i]]$detcov
        dosarpvals[i, ] <- res[[i]]$dosarpvals
        for (k in 1:np){
            pararr[k, i, ] <- res[[i]]$parmat[k, ]
        }
    }
    dfbeta <- matrix(0, nobs, np)
    colnames(dfbeta) <- parnamsnolog
    rownames(dfbeta) <- rwnms
    dpar <- dfbeta
    for (k in 1:np){
        if (doexp[k]==0){
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
    for (i in 1:nser){
        chgmat[, i] <- newres[, i]-orgres[i]
    }
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
    for (i in 1:ninfl){
        inds <- which(infl[, i]==1)
        infl[inds, i] <- i
    }
    rep$infl <- list(dfbeta=dfbeta, dpar=dpar, ddetcov=ddetcov, dosarpvals=dosarpvals, infl=infl)
    # Add to stats
    rep$stats$inflperc <- sum(apply(is.na(infl), 1, sum)<dim(infl)[2])/dim(infl)[1]
    return(rep)
}
