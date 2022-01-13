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


#' @name check.ini
#' @title Check sensitivity of fit to initial parameter values
#' @param input Either an inp list passing check.inp(), or a rep list where rep is the output of running fit.spict().
#' @param ntrials The number of trials with different starting values to run.
#' @param verbose If true write information to screen.
#' @param numdigits Number of digits in reported results.
#' @return List containing results of sensitivity check and associated initial values.
#' @export
check.ini <- function(input, ntrials=10, verbose=TRUE, numdigits=2){
    nd <- numdigits
    if ('par.fixed' %in% names(input)){
        rep <- input
        inp <- rep$inp
    } else {
        if ('obsC' %in% names(input)){
            inp <- check.inp(input)
            rep <- fit.spict(inp)
        } else {
            stop('Invalid input! use either an inp list or a fit.spict() result.')
        }
    }
    calc.dist <- function(vec1, vec2){
        return(sqrt(sum((as.numeric(vec1) - as.numeric(vec2))^2)))
    }
    # Start sensitivity check by drawing random initial parameters
    if ('par.fixed' %in% names(rep)){
        nms <- intersect( names(inp$ranges), names(rep$par.fixed) )
        inibasevec <- unlist(inp$ini[nms])
        resbasevec <- trans2real(rep$opt$par, names(rep$opt$par))
        resmat <- matrix(nrow=ntrials, ncol=length(rep$par.fixed)+1)
        colnames(resmat) <- c(names(rep$par.fixed),"objective_function_value")
        inimat <- matrix(nrow=ntrials, ncol=length(rep$par.fixed))
        colnames(inimat) <- names(rep$par.fixed) 
        propchng <- inimat
        rownames(propchng) <- paste('Trial', 1:ntrials)
        resdist <- numeric(ntrials)
        inidist <- numeric(ntrials)
        if(verbose){
            cat('Checking sensitivity of fit to initial parameter values...\n')
        }
        for (i in 1:ntrials){
            if(verbose){
                cat(' Trial', i, '...')
            }
            inpsens <- rep$inp
            inpsens$do.sd.report <- FALSE # Only check point estimates
            inpsens$ini$logB <- NULL
            inpsens$ini$logF <- NULL
            inpsens$ini$logu <- NULL
            # Draw random initial values
            j <- 1
            for (nm in nms){
                for(kk in 1:nrow(inp$ranges[[nm]])) {
                    randval <- runif(1, min = inp$ranges[[nm]][kk,1], 
                                     max = inp$ranges[[nm]][kk,2])
                    inpsens$ini[[nm]][kk] <- randval
                    inimat[i, j] <- randval
                    propchng[i, j] <- (randval - inp$ini[[nm]][kk])/inp$ini[[nm]][kk]
                    j <- j + 1
                }
            }
            inidist[i] <- calc.dist(inimat[i, ], inibasevec)
            inpsens$osar.method <- "none"
            repsens <- try(fit.spict(inpsens))
            if (class(repsens) != 'try-error' & 'opt' %in% names(repsens)){
                if (repsens$opt$convergence == 0){
                    resmat[i, 1:(ncol(resmat)-1)] <- trans2real(repsens$opt$par, names(repsens$opt$par))
                    resmat[i,ncol(resmat)]        <- repsens$opt$objective
                    resdist[i] <- calc.dist(resmat[i, ], resbasevec)
                    if(verbose){
                        cat(' model fitted!\n')
                    }
                } else {
                    if(verbose){
                        cat(' convergence not obtained!\n')
                    }
                }
            } else {
                if(verbose){
                    cat(' fit failed!\n')
                }
            }
        }
        resmat <- cbind(resdist, resmat)
        resmat <- rbind(c(0, resbasevec), resmat)
        rownames(resmat) <- c('Basevec', paste('Trial', 1:ntrials))
        colnames(resmat)[1] <- 'Distance'
        inimat <- cbind(inidist, inimat)
        inimat <- rbind(c(0, inibasevec), inimat)
        rownames(inimat) <- c('Basevec', paste('Trial', 1:ntrials))
        colnames(inimat)[1] <- 'Distance'
        input$check.ini <- list(propchng=round(propchng, nd), inimat=round(inimat, nd),
                                resmat=round(resmat, nd))
        if (verbose) {
          print(input$check.ini)
        }
        return(input)
    }
}
