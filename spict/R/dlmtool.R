#' Get function to estimate TAC for the DLMtool package
#' 
#' This function creates harvest control rules (HCRs) which can be incorporated into a
#' management strategy evaluation framework (DLMtool package). HCRs are saved with a
#' generic name to the global environment and the names of the HCRs are returned if results of the
#' function are assigned to an object. HCR runs a SPiCT assessment using catch and
#' relative biomass index observations. Stock status estimates are used to set the TAC
#' for the next year. TAC can be based on the distribution of predicted catches (percentileC)
#' and/or the distribution of the Fmsy reference level (percentileFmsy).
#' Additionally, a cap can be applied to account for low biomass levels (below Bmsy).
#' Arguments of returned function are 'x' - the position in a data-limited mehods data object,
#' 'Data' - the data-limited methods data object (see DLMtool), and 'reps' - the number of
#' stochastic samples of the TAC recommendation (not used for this HCR).
#' One or several arguments of the function can be provided as vectors to generate several
#' HCRs at once (several vectors have to have same length).
#'
#' @param fractileC The fractile of the catch distribution to be used for setting TAC. Default
#'   is median (0.5).
#' @param fractileFFmsy The fractile of the distribution of F/Fmsy. Default is 0.5 (median).
#' @param pa Logical; indicating if the precautionary approach should be applied (reduce F if P(B<Blim) < prob). Default is FALSE.
#' @param prob Probability for the precautionary approach (see argument 'pa', default is 0.95).
#' @param fractileBBmsy The fractile of the distribution of B/Bmsy. Default is 0.5 (median).
#' @param uncertaintyCap Logical; If true TAC is bound between two values set in lower and upper. Default: FALSE.
#' @param lower lower bound of the uncertainty cap. Default is 0.8, used if uncertaintyCap = TRUE.
#' @param upper upper bound of the uncertainty cap. Default is 1.2, used if uncertaintyCap = TRUE.
#' @param interval Assessment interval. Default is 1, which indicates annual assessments.
#' @param env environment where the harvest control rule function(s) are assigned to.
#' @return A function which can estimate TAC recommendations based on SPiCT assessment,
#'   taking assessment uncertainty into account.
#' 
#' @export
#'
#' @examples
#' \dontrun{
#' ## Put together an operating model from the available DLM toolkit examples
#' StockEx <- Herring
#' FleetEx <- Generic_IncE
#' ObsEx <- Precise_Unbiased
#' ## Remove changes in life history parameters
#' StockEx@Mgrad <- c(0,0)
#' StockEx@Kgrad <- c(0,0)
#' StockEx@Linfgrad <- c(0,0)
#' StockEx@Size_area_1 <- c(0.5,0.5)
#' StockEx@Frac_area_1 <- c(0.5,0.5)
#' StockEx@Prob_staying <- c(0.5,0.5)
#' ## Set the depletion level 
#' StockEx@D <- c(0.3, 0.4)
#' ## create Operation Model
#' OMex <- new("OM", Stock = StockEx, Fleet = FleetEx, 
#'                 Obs = ObsEx)
#' ## Set simulation options
#' OMex@nsim <- 10
#' OMex@nyears <- 25
#' OMex@proyears <- 5
#' ## Get SPiCT HCRs
#' MPname <- c(spict2DLMtool(),
#'             spict2DLMtool(pa=TRUE))
#' ## run MSE
#' MSEex <- runMSE(OMex,
#'                 MPs = MPname,
#'                 timelimit = 150,
#'                 CheckMPs = FALSE)
#' ## example plot of results
#' Pplot2(MSEex, traj="quant", quants=c(0.2, 0.8))
#' }
#'
spict2DLMtool <- function(fractileC = 0.5,
                          fractileFFmsy = 0.5,
                          pa = FALSE,
                          prob = 0.95,
                          fractileBBmsy = 0.5,
                          uncertaintyCap = FALSE,
                          lower = 0.8,
                          upper = 1.2,
                          interval = 1,
                          env = globalenv()){

    ## allowing for multiple generation of MPs
    argList <- list(fractileC, fractileFFmsy, pa, prob, fractileBBmsy,
                    uncertaintyCap, lower, upper)
    argLengths <- sapply(argList, length)
    maxi <- max(argLengths)
    maxl  <- which(argLengths == maxi)
    if(maxi>1){
        if(max(argLengths[(1:8)[-maxl]]) > 1)
            stop("Specified arguments have different lengths, they should have the same length or length = 1.")
    }
    argListCor <- argList
    argListCor[(1:8)[-maxl]] <- lapply(argList[(1:8)[-maxl]], function(x) rep(unlist(x), maxi))


    ## MP as function
    template  <- expression(paste0(
        'structure(function(x, Data, reps = 1,
          fractileC = ',a,',
          fractileFFmsy=',b,',
          pa=',c,',
          prob=',d,',
          fractileBBmsy=',e,',
          uncertaintyCap=',f,',
          lower=',g,',
          upper=',h,'){
            dependencies <- "Data@Year, Data@Cat, Data@Ind"
            time <- Data@Year
            Catch <- Data@Cat[x,]
            Index <- Data@Ind[x,]
            inp <- list(timeC=time, obsC=Catch, 
                        timeI=time, obsI=Index,
                        dteuler = 1 / 16,
                        do.sd.report=TRUE,
                        getReportCovariance = FALSE)
            inp <- check.inp(inp)
            inp$timepredi <- inp$timepredc + ',interval,'
            rep <- try(spict::fit.spict(inp),silent=TRUE)
            if(is(rep, "try-error") || rep$opt$convergence != 0) {
                TAC <- rep(NA, reps)
            } else {
                ## Reduction based on uncertainty in Fmsy. Default = median
                idx <- rep$inp$indpred[1]
                logFFmsy <- spict::get.par("logFFmsy", rep)[idx,]
                fi <- 1-fractileFFmsy
                fm <- exp( qnorm( fi, logFFmsy[2], logFFmsy[4] ) )
                fm5 <- exp( qnorm( 0.5, logFFmsy[2], logFFmsy[4] ) )
                red <- fm5 / fm            
                ## precautionary approach
                if(pa == 1){
                    Fmsy <- get.par("logFmsy", rep, exp=TRUE)[2]
                    Flast <- get.par("logF", rep, exp=TRUE)[rep$inp$indpred[1], 2]            
                    ffac <- (red + 1e-6) * Fmsy / Flast
                    bbmsyQ5 <- spict:::probdev(ffac, rep, bbmsyfrac=fractileBBmsy,
                                                     prob=prob, getFrac=TRUE)
                    if((0.5 - bbmsyQ5) > 0.001){
                        fy <- spict:::get.PAffac(rep, bbmsyfrac=fractileBBmsy,
                                                 prob=prob)
                        red <- fy * Flast / Fmsy
                    }
                }
                ## Uncertainty cap
                if(uncertaintyCap){
                    red[red < lower] <- lower
                    red[red > upper] <- upper
                }
                predcatch <- try(spict::pred.catch(rep, MSEmode = 1,
                                                   get.sd = TRUE, exp = FALSE, fmsyfac = red),
                                 silent=TRUE)
                if(is(predcatch, "try-error")) {
                    TAC <- rep(NA, reps)
                } else {
                    TACi <- exp(qnorm(fractileC, predcatch[2], predcatch[4]))
                    ## Reduction based on B/Bmsy. Default = median
                    ##idx <- rep$inp$indpred[1]
                    logBBmsy <- spict::get.par("logBpBmsy", rep, exp = TRUE)   ##[idx,]
                    predBBtrigger <- 2 * exp(qnorm(fractileBBmsy, logBBmsy[2], logBBmsy[4]))
                    TACi <- TACi * min(1, predBBtrigger)
                    ## hack to guarantee compatibility with other MPs (DLMtool takes median, thus rep no effect)
                    TAC <- rep(TACi, reps)
                }
            }
            res <- TACfilter(TAC)
            Rec <- new("Rec")
            Rec@TAC <- res
            return(Rec)
        },
        class="MP")'))


    nami <- rep(NA,maxi)
    for(I in 1:maxi){

        ## create MPs as functions
        subList <- lapply(argListCor, "[[", I)
        names(subList) <- letters[1:8]
        templati <- eval(parse(text=paste(parse(text = eval(template, subList)),collapse=" ")))

        ## save names of MPs
        if(argListCor[[1]][I] == 0.5){
            c1 <- ""
        }else{
            c1 <- paste0("_C",argListCor[[1]][I])
        }
        if(argListCor[[2]][I] == 0.5){
            c2 <- ""
        }else{
            c2 <- paste0("_FFmsy",argListCor[[2]][I])
        }
        if(argListCor[[3]][I] == FALSE){
            c3 <- ""
        }else{
            c3 <- paste0("_pa")
        }
        if(argListCor[[4]][I] == 0.95){
            c4 <- ""
        }else{
            c4 <- paste0("_P",argListCor[[4]][I])            
        }
        if(argListCor[[5]][I] == 0.5){
            c5 <- ""
        }else{
            c5 <- paste0("_BBmsy",argListCor[[5]][I])            
        }
        if(argListCor[[6]][I] == FALSE){
            c6 <- ""
        }else{
            c6 <- paste0("_uC")            
        }
        if(argListCor[[7]][I] == 0.8){
            c7 <- ""
        }else{
            c7 <- paste0("_lo",argListCor[[7]][I])            
        }
        if(argListCor[[8]][I] == 1.2){
            c8 <- ""
        }else{
            c8 <- paste0("_up",argListCor[[8]][I])            
        }
        ## put everythin together
        nami[I] <- paste0("spict",c1,c2,c3,c4,c5,c6,c7,c8)
        assign(value=templati, x=nami[I], envir=env)
    }

    ## allow for assigning names
    invisible(nami)        
}
