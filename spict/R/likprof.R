
#' @name likprof.spict
#' @title Create profile likelihood
#' @details The "likprof" list must containg the following keys:
#' \itemize{
#'   \item{"pars"}{ A character vector of length equal 1 or 2 containing the name(s) of the parameters to calculate the profile likelihood for.}
#'   \item{"parrange"}{ A vector containing the parameter range(s) to profile over: parrange = c(min(par1), max(par1), min(par2), max(par2)).}
#' }
#' Optional:
#' \itemize{
#'   \item{"nogridpoints"}{ Number of grid points to evaluate the profile likelihood for each parameter. Default: 9. Note: with two parameters the calculation time increases quadratically when increasing the number of gridpoints.}
#' }
#' @param input A list containing observations and initial values for non profiled parameters (essentially an inp list) with the additional key "likprof" (see details for required keys). A valid result from fit.spict() containing an "inp" key with the described properties is also accepted.
#' @return The output is the input with the likelihood profile information added to the likprof key of either inp or rep$inp.
#' @examples
#' data(pol)
#' inp <- pol$albacore
#' inp$likprof <- list()
#' inp$likprof$pars <- 'logsdb'
#' inp$likprof$parrange <- c(log(0.05), log(0.4))
#' inp$likprof$nogridpoints <- 15
#' rep <- fit.spict(inp)
#' rep <- likprof.spict(rep)
#' plotspict.likprof(rep, logpar=TRUE)
#' @export
likprof.spict <- function(input){
    if('par.fixed' %in% names(input)){
        cat('Detected input as a SPiCT result, proceeding...\n')
        rep <- input
        inp <- rep$inp
        pl <- rep$pl
        repflag <- TRUE
    } else {
        cat('Assuming input is a SPiCT inp, checking for validity...\n')
        inp <- input
        inp <- check.inp(inp)
        repflag <- FALSE
    }
    likprof <- inp$likprof
    # Check likprof input
    if(!'pars' %in% names(likprof)) stop('"pars" key of likprof is unspecified!')
    np <- length(likprof$pars)
    if(!np %in% 1:2) stop('Length of pars vector is ', np, ' but should be 1 or 2!')
    if(!'parrange' %in% names(likprof)) stop('No "parrange" key specified in likprof!')
    npp <- length(likprof$parrange)
    if(npp != (2*np)) stop('Length of "parrange" is ', npp, ' but should be ', 2*np)
    if(likprof$parrange[2] <= likprof$parrange[1]) stop('Upper parameter bound is small than lower bound!')
    if(np==2) if(likprof$parrange[4] <= likprof$parrange[3]) stop('Upper parameter bound is small than lower bound!')
    prange <- matrix(likprof$parrange, np, 2, byrow=TRUE)
    if(!'nogridpoints' %in% names(likprof)) likprof$nogridpoints <- 9
    if(!'phases' %in% names(inp)) inp$phases <- list()
    # Determine good initial values for other parameters
    if(!repflag){
        inpp <- inp
        inpp$phases <- list()
        pp <- rep(0, np)
        for(i in 1:np){
            pp[i] <- mean(prange[i, ])
            inpp$ini[[inp$likprof$pars[i]]] <- pp[i]
            inpp$phases[[inp$likprof$pars[i]]] <- -1
        }
        repp <- fit.spict(inpp)
        pl <- repp$pl
    }
    parvals <- matrix(0, likprof$nogridpoints, np)
    for(i in 1:np) parvals[, i] <- seq(prange[i, 1], prange[i, 2], length=likprof$nogridpoints)
    if(np==1){
        pv <- parvals[, 1, drop=FALSE]
    } else {
        pv <- expand.grid(parvals[, 1], parvals[, 2])
    }
    ngrid <- dim(pv)[1]
    do.grid <- function(i){
        inptmp <- inp
        inptmp$ini <- pl
        for(j in 1:np){
            inptmp$ini[[likprof$pars[j]]] <- pv[i, j]
            inptmp$phases[[likprof$pars[j]]] <- -1
        }
        reptmp <- try(fit.spict(inptmp))
        if(class(reptmp)=='try-error'){
            cat('likprof WARNING: par[i], i=', i, 'failed to fit!\n')
            rep2 <- list()
            rep2$opt <- list()
            rep2$opt$objective <- NA
        } else {
            rep2 <- reptmp
        }
        #cat(i, '..')
        return(rep2$opt$objective)
    }
    single.do.grid <- function(i, ngrid){
        lv <- do.grid(i)
        cat(sprintf("Profiling: %.2f%%\r", i/ngrid * 100))
        return(lv)
    }
    multi.do.grid <- function(i, ngrid, progfile){
        lv <- do.grid(i)
        ## Send progress update
        if(class(progfile)[1]=='fifo') writeBin(1/ngrid, progfile)
        return(lv)
    }
    #partry <- try(library(multicore))
    #partry <- try(library(parallel))
    cat('Profiling pars:', paste(likprof$pars, collapse=' and '), 'with', ngrid, 'trials.\n')
    if(Sys.info()['sysname']!='Linux'){
    #if(class(partry)=='try-error'){
    #if(TRUE){
        likvals <- lapply(1:ngrid, single.do.grid, ngrid)
    } else {
        require(parallel)
        ## Open fifo progress file
        progfile <- fifo(tempfile(), open="w+b", blocking=T)
        #if(inherits(fork(), "masterProcess")) {
        if(inherits(parallel:::mcfork(), "masterProcess")) {
            ## Child
            progress <- 0.0
            while (progress < 1 && !isIncomplete(progfile)) {
                msg <- readBin(progfile, "double")
                progress <- progress + as.numeric(msg)
                cat(sprintf("Multicore profiling: %.2f%%\r", progress * 100))
            } 
            #exit()
        }
        likvals <- mclapply(1:ngrid, multi.do.grid, ngrid, progfile, mc.cores=4)
        ## Close progress file
        close(progfile)
    }
    if(np==1){
        likprof$likvals <- unlist(likvals)
    } else {
        likprof$likvals <- matrix(unlist(likvals), likprof$nogridpoints, likprof$nogridpoints)
    }
    likprof$parvals <- parvals
    if(repflag){
        rep$inp$likprof <- likprof
        return(rep)
    } else {
        inp$likprof <- likprof
        return(inp)
    }
}