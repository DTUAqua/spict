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


#' @name read.aspic
#' @title Reads ASPIC input file.
#' @details Reads an input file following the ASPIC 7 format described in the ASPIC manual (found here http://www.mhprager.com/aspic.html).
#' @param filename Path of the ASPIC input file.
#' @return A list of input variables that can be used as input to fit.spict().
#' @examples
#' \dontrun{
#' filename <- 'YFT-SSE.a7inp' # or some other ASPIC 7 input file
#' inp <- read.aspic(filename)
#' rep <- fit.spict(inp)
#' summary(rep)
#' plot(rep)
#' }
#' @export
read.aspic <- function(filename){
    rawdat <- readLines(filename)
    # Read meta data
    dat <- list()
    c <- 0
    i <- 0
    iq <- 0
    found <- FALSE
    while(c < 16){
        i <- i+1
        found <- !substr(rawdat[i],1,1)=='#'
        c <- c + found
        if (c==1 & found) dat$version <- rawdat[i]
        if (c==2 & found) dat$title <- rawdat[i]
        # Number of observations
        if (c==5 & found) dat$nobs <- as.numeric(read.table(filename, skip=i-1, nrows=1)[1])
        # Number of data series
        if (c==5 & found){
            dat$nseries <- as.numeric(read.table(filename, skip=i-1, nrows=1)[2])
            dat$qini <- rep(0, dat$nseries)
        }
        # Initial parameters
        if (c==10 & found) dat$B1Kini <- as.numeric(read.table(filename, skip=i-1, nrows=1)[2])
        if (c==11 & found) dat$MSYini <- as.numeric(read.table(filename, skip=i-1, nrows=1)[2])
        if (c==12 & found) dat$Fmsyini <- as.numeric(read.table(filename, skip=i-1, nrows=1)[2])
        if ('nseries' %in% names(dat)){
            if (c==(13+iq) & found & iq < dat$nseries){
                iq <- iq+1
                dat$qini[iq] <- as.numeric(read.table(filename, skip=i-1, nrows=1)[2])
            }
        }
        found <- FALSE
    }
    # Read actual observations
    obsdat <- list()
    idat <- which(rawdat=='DATA')
    type <- rep('', dat$nseries)
    name <- rep('', dat$nseries)
    for (j in 1:dat$nseries){
        c <- 0
        i <- 0
        while(c < 2){
            i <- i+1
            c <- c + !substr(rawdat[idat+i],1,1)=='#'
            if (c==1) name[j] <- rawdat[idat+i]
            if (c==2) type[j] <- rawdat[idat+i]
        }
        datstart <- idat+i
        obsdat[[j]] <- read.table(filename, skip=datstart, sep='', nrows=dat$nobs)
        obsdat[[j]][obsdat[[j]]<0] <- NA
        if (type[j] == 'CE'){
            obsdat[[j]] <- obsdat[[j]][, 1:3]
            colnames(obsdat[[j]]) <- c('time', 'effort', 'catch')
            obsdat[[j]] <- cbind(obsdat[[j]], cpue=obsdat[[j]]$catch/obsdat[[j]]$effort)
        }
        if (type[j] == 'CC'){
            obsdat[[j]] <- obsdat[[j]][, 1:3]
            colnames(obsdat[[j]]) <- c('time', 'cpue', 'catch')
        }
        # ASPIC differentiates between these index types, but SPiCT does not (yet). See the ASPIC manual for details.
        if (type[j] %in% c('I0','I1','I2','B0','B1','B2')){ 
            obsdat[[j]] <- obsdat[[j]][, 1:2]
            colnames(obsdat[[j]]) <- c('time', 'index')
        }
        idat <- datstart + dat$nobs
    }
    # Convert to SPiCT data by storing as an inp list
    inp <- list()
    # Insert catches
    ind <- grep('C', type)
    inp$timeC <- obsdat[[ind]]$time[!is.na(obsdat[[ind]]$catch)]
    inp$obsC <- obsdat[[ind]]$catch[!is.na(obsdat[[ind]]$catch)]
    # Insertindices
    inp$timeI <- list()
    inp$obsI <- list()
    for (j in 1:dat$nseries){
        if (type[j] %in% c('CC', 'CE')){
            inp$timeI[[j]] <- obsdat[[j]]$time[!is.na(obsdat[[j]]$cpue)]
            inp$obsI[[j]] <- obsdat[[j]]$cpue[!is.na(obsdat[[j]]$cpue)]
        }
        if (type[j] %in% c('I0','I1','I2','B0','B1','B2')){
            inp$timeI[[j]] <- obsdat[[j]]$time[!is.na(obsdat[[j]]$index)]
            inp$obsI[[j]] <- obsdat[[j]]$index[!is.na(obsdat[[j]]$index)]
        }
    }
    # Insert initial values
    inp$ini <- list()
    inp$ini$logr <- log(2*dat$Fmsyini)
    inp$ini$logK <- log(2*dat$MSYini/dat$Fmsyini)
    inp$ini$logq <- log(dat$qini)
    return(inp)
}


#' @name write.aspic
#' @title Takes a SPiCT input list and writes it as an Aspic input file.
#' @details TBA
#' @param input List of input variables or the output of a simulation using sim.spict().
#' @param filename Name of the file to write.
#' @param verbose If true write information to screen.
#' @return Noting.
#' @examples
#' data(pol)
#' sim <- (pol$albacore)
#' write.aspic(sim)
#' @export
write.aspic <- function(input, filename='spictout.a7inp', verbose=FALSE){
    inp <- check.inp(input)
    if (any(c(inp$timeC, unlist(inp$timeI)) %% 1 != 0)){
        cat('Warning (write.aspic): Observation times were rounded down (floored) to integers. Consider providing only integer observation times.\n')
    }
    inp$timeC <- floor(inp$timeC)
    for (i in 1:inp$nindex){
        inp$timeI[[i]] <- floor(inp$timeI[[i]])
    }
    timeobs <- sort(unique(c(inp$timeC, inp$timeI[[1]])))
    nobs <- length(timeobs)
    dat <- cbind(timeobs, rep(-1, nobs), rep(-1, nobs))
    inds <- match(inp$timeI[[1]], timeobs)
    dat[inds, 2] <- inp$obsI[[1]]
    inds <- match(inp$timeC, timeobs)
    dat[inds, 3] <- inp$obsC
    if (verbose){
        cat('Writing input to:', filename, '\n')
    }
    cat('ASPIC-V7\n', file=filename)
    cat(paste('# File generated by SPiCT function write.aspic at', Sys.time(), '\n'),
        file=filename, append=TRUE)
    cat('"Unknown stock"\n', file=filename, append=TRUE)
    cat('# Program mode (FIT/BOT), verbosity, N bootstraps, [opt] user percentile:\n',
        file=filename, append=TRUE)
    cat(paste0(inp$aspic$mode, '  ', inp$aspic$verbosity, '  ', inp$aspic$nboot, '  ',
               inp$aspic$ciperc, '\n'), file=filename, append=TRUE)
    #cat('BOT  102  1000  95\n', file=filename, append=TRUE)
    cat('# Model shape, conditioning (YLD/EFT), obj. fn. (SSE/LAV/MLE/MAP):\n',
        file=filename, append=TRUE)
    cat('LOGISTIC  YLD  SSE\n', file=filename, append=TRUE)
    cat('# N years, N series:\n', file=filename, append=TRUE)
    cat(paste0(nobs, '  ', inp$nindex, '\n'), file=filename, append=TRUE)
    cat('# Monte Carlo mode (0/1/2), N trials:\n', file=filename, append=TRUE)
    cat('0  30000\n', file=filename, append=TRUE)
    cat('# Convergence criteria (3 values):\n', file=filename, append=TRUE)
    cat('1.00E-08  3.00E-08  1.00E-04\n', file=filename, append=TRUE)
    cat('# Maximum F, N restarts, [gen. model] N steps/yr:\n', file=filename, append=TRUE)
    cat('8.00E+00  6  24\n', file=filename, append=TRUE)
    cat('# Random seed (large integer):\n', file=filename, append=TRUE)
    cat('1234\n', file=filename, append=TRUE)
    cat('# Initial guesses and bounds follow:\n', file=filename, append=TRUE)
    estbkfrac <- 1
    if (!is.null(inp$phases$logbkfrac)){
        if (inp$phases$logbkfrac==-1){
            estbkfrac <- 0
        }
    }
    #estbkfrac <- 0
    if (!'logbkfrac' %in% names(inp$ini)){
        inp$ini$logbkfrac <- log(0.8)
    }
    cat(sprintf('B1K   %3.2E  %i  %3.2E  %3.2E  penalty  %3.2E\n',
                exp(inp$ini$logbkfrac), estbkfrac, 0.01*exp(inp$ini$logbkfrac),
                100*exp(inp$ini$logbkfrac), 0), file=filename, append=TRUE)
    MSY <- exp(inp$ini$logK + inp$ini$logr)/4
    cat(sprintf('MSY   %3.2E  1  %3.2E  %3.2E\n', MSY, 0.03*MSY, 5000*MSY),
        file=filename, append=TRUE)
    Fmsy <- exp(inp$ini$logr)/2
    cat(sprintf('Fmsy  %3.2E  1  %3.2E  %3.2E\n', Fmsy, 0.01*Fmsy, 100*Fmsy),
        file=filename, append=TRUE)
    for (i in 1:inp$nindex){
        cat(sprintf('q     %3.2E  1  %3.2E  %3.2E  %3.2E\n',
                    exp(inp$ini$logq[i]), 1, 0.001*exp(inp$ini$logq[i]),
                    100*exp(inp$ini$logq[i])), file=filename, append=TRUE)
    }
    cat('DATA\n', file=filename, append=TRUE)
    cat('"Combined-Fleet Index, Total Landings"\n', file=filename, append=TRUE)
    cat('CC\n', file=filename, append=TRUE) # CC refers to CPUE annual average, catch annual total.
    for (i in 1:nobs){
        cat(sprintf('  %4i    % 6.4E    % 6.4E\n', timeobs[i], dat[i, 2],
                    dat[i, 3]), file=filename, append=TRUE)
    }
    if (inp$nindex>1){
        for (I in 2:inp$nindex){
            cat(paste0('"Index',I-1,'"\n'), file=filename, append=TRUE)
            cat('I1\n', file=filename, append=TRUE)
            dat <- cbind(timeobs, rep(-1, nobs))
            inds <- match(inp$timeI[[I]], timeobs)
            dat[inds, 2] <- inp$obsI[[I]]
            for (i in 1:nobs){
                cat(sprintf('  %4i    % 6.4E\n', dat[i ,1], dat[i, 2]),
                    file=filename, append=TRUE)
            }
        }
    }
}


#' @name read.aspic.res
#' @title Reads the parameter estimates of an Aspic result file.
#' @details TBA
#' @param filename Name of the Aspic result file to read
#' @return Vector containing the parameter estimates.
#' @export
read.aspic.res <- function(filename){
    out <- list()
    aspicres <- readLines(filename)
    ind <- grep('Number of years analyzed', aspicres)
    ind2 <- gregexpr('[0-9]+', aspicres[ind])
    nobs <- as.numeric(unlist(regmatches(aspicres[ind], ind2))[1])
    get.num <- function(str, aspicres){
        ind <- grep(str, aspicres)
        ind2 <- gregexpr('[0-9.]*E[+-][0-9]*', aspicres[ind])
        as.numeric(unlist(regmatches(aspicres[ind], ind2)))[1]
    }
    # Error code
    errind <- grep('error code', aspicres)[1]
    out$errorcode <- regmatches(aspicres[errind], regexpr('[0-9]+$', aspicres[errind]))
    out$convmsg <- aspicres[errind + 2]
    norestartind <- grep('Number of restarts', aspicres)
    out$no.restarts <- as.numeric(regmatches(aspicres[norestartind],
                                             regexpr('[0-9]+$', aspicres[norestartind])))
    # Parameters
    bkfrac <- get.num('^B1/K', aspicres)
    MSY <- get.num('^MSY', aspicres)
    Fmsy <- get.num('^Fmsy', aspicres)
    q <- get.num('^q', aspicres)
    Bmsy <- MSY/Fmsy
    K <- 2*Bmsy
    r <- 2*Fmsy
    # States
    ind <- grep('ESTIMATED POPULATION TRAJECTORY', aspicres)
    out$states <- read.table(filename, skip=ind+6, sep='', nrows=nobs, strip.white=TRUE)
    out$states <- read.table(filename, skip=ind+6, sep='', nrows=nobs, strip.white=TRUE)
    colnames(out$states) <- c('obs', 'time', 'Fest', 'B0est', 'Best', 'Catch', 'Cest', 'Pest', 'FFmsy', 'BBmsy')
    out$pars <- c(r=r, K=K, q=q, Fmsy=Fmsy, Bmsy=Bmsy, MSY=MSY)
    # Time spent
    string <- strsplit(aspicres[length(aspicres)], ' ')
    et <- as.numeric(c(string[[1]][3], string[[1]][5], string[[1]][8]))
    timespent <- et[1]*60*60 + et[2]*60 + et[3]
    #names(timespent) <- c('hr', 'min', 'sec')
    out$timespent <- timespent
    return(out)
}



#' @name fit.aspic
#' @title Fits aspic to the data contained in the input file
#' @details Only works on Linux. This furthermore requires that wine is installed and that aspic7 is installed and available to the PATH.
#' @param input A spict input list containing observations.
#' @param do.boot Do bootstrap to get uncertainties of estimates?
#' @param nboot Number of bootstrap runs (only used if do.boot=TRUE). Prager suggests in the ASPIC manual p. 13 to use nboot > 1000 if ciperc > 80.
#' @param ciperc Coverage percentage (integer between 0 and 100) of bootstrapped confidence intervals.
#' @param verbose If TRUE write information to screen.
#' @param filebase Basename of all generated aspic files.
#' @param savefile Save results to this file.
#' @return List containing aspic results.
#' @export
fit.aspic <- function(input, do.boot=FALSE, nboot=NULL, ciperc=NULL, verbose=FALSE, filebase='tmp', savefile=NULL){
    # Write aspic input file
    fn <- paste0(filebase, '.a7inp')
    write.aspic(input, filename=fn, verbose=verbose)

    # Fit aspic
    if (file.exists(paste0(filebase, '.fit'))){
        file.remove(paste0(filebase, '.fit')) # Remove old fit
    }
    system(paste0('wine aspic7 ', fn), intern=!verbose)
    #system(paste0('aspic7.exe ', fn)) # Alternative system call

    # Read aspic results
    res <- read.aspic.res(paste0(filebase, '.fit'))

    #if (do.boot & res$timespent < input$aspic$bootlimtime){
    if (do.boot & res$errorcode == '0'){
        # Bootstrap aspic to get uncertainties
        inpaspic <- input
        inpaspic$aspic$mode <- 'BOT'
        if (!is.null(nboot)){
            inpaspic$aspic$nboot <- nboot
        }
        if (!is.null(ciperc)){
            stopifnot(ciperc > 0, ciperc < 100)
            inpaspic$aspic$ciperc <- ciperc
        }

        # Run 10 bootstrap samples to check calculation time
        inpaspictest <- inpaspic
        inpaspictest$aspic$nboot <- 10
        write.aspic(inpaspictest, filename=fn, verbose=verbose)
        if (file.exists(paste0(filebase, '.bot'))){
            file.remove(paste0(filebase, '.bot'))
        }
        system(paste0('wine aspic7 ', filebase, '.a7inp'), intern=!verbose)
        filename <- paste0(filebase, '.bot')
        aspicres <- readLines(filename)
        string <- strsplit(aspicres[length(aspicres)], ' ')
        string <- subset(string[[1]], string[[1]] != '') # Remove blanks
        et <- as.numeric(c(string[3], string[5], string[7]))
        boottime <- et[1]*60*60 + et[2]*60 + et[3]

        if (boottime < input$aspic$bootlimtime){
            # Run full bootstrap with all samples if preliminary run was quick enough
            write.aspic(inpaspic, filename=fn, verbose=verbose)
            if (file.exists(paste0(filebase, '.bot'))){
                file.remove(paste0(filebase, '.bot'))
            }
            #system(paste0('aspic7.exe ', filebase, '.a7inp'), intern=!verbose)
            system(paste0('wine aspic7 ', filebase, '.a7inp'), intern=!verbose)
            
            # Read aspic bootstrap results
            filename <- paste0(filebase, '.bot')
            aspicres <- readLines(filename)
            ind <- grep('ESTIMATES FROM BOOTSTRAP ANALYSIS', aspicres)
            bootres <- read.table(filename, skip=ind+5, sep='', nrows=11, strip.white=TRUE)
            rnms <- as.character(bootres[, 1])
            bootres <- bootres[, -1]
            rownames(bootres) <- rnms
            colnames(bootres) <- c('estimate', 'low80', 'upp80', 'low95', 'upp95', 'IQR', 'rIQR')
            res$boot <- bootres
        } else {
            res$booterrorcode <- paste('Bootstrap of', inpaspictest$aspic$nboot, 'trials took', boottime, 'seconds, which was longer than allowed:', input$aspic$bootlimtime)
        }
    }

    # Save results
    if (!is.null(savefile)){
        save('res', file=savefile)
    }
    
    # Clean temporary files
    
    return(res)
}
