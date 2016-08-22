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
#' @return List containing results of sensitivity check and associated initial values.
#' @export
check.ini <- function(input, ntrials=10, verbose=TRUE){
    if ('par.fixed' %in% names(input)){
        rep <- input
        inp <- rep$inp
    } else {
        if ('obsC' %in% names(input)){
            inp <- check.inp(input)
        } else {
            stop('Invalid input! use either an inp list or a fit.spict() result.')
        }
    }
    if (!exists('rep')){
        rep <- fit.spict(inp)
    }

    # Start sensitivity check by drawing random initial parameters
    if ('par.fixed' %in% names(rep)){
        nms <- names(inp$ranges)
        inibasevec <- unlist(inp$ini[names(inp$ranges)])
        resbasevec <- trans2real(rep$opt$par, names(rep$opt$par))
        resmat <- matrix(nrow=ntrials, ncol=length(rep$par.fixed))
        colnames(resmat) <- names(resbasevec)
        inimat <- matrix(nrow=ntrials, ncol=length(nms))
        colnames(inimat) <- nms
        perchange <- inimat
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
                randval <- runif(1, min=inp$ranges[[nm]][1], max=inp$ranges[[nm]][2])
                inpsens$ini[[nm]] <- randval
                inimat[i, j] <- randval
                perchange[i, j] <- randval - inp$ini[[nm]]
                j <- j + 1
            }
            repsens <- try(fit.spict(inpsens))
            if (class(repsens) != 'try-error' & 'opt' %in% names(repsens)){
                if (repsens$opt$convergence == 0){
                    resmat[i, ] <- trans2real(repsens$opt$par, names(repsens$opt$par))
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
        return(list(inibasevec=inibasevec, perchange=perchange, inimat=inimat,
                    resbasevec=resbasevec, resmat=resmat))
    }
}
