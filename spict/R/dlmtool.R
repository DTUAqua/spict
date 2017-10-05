#' @name spict2DLMtool
#' @title Get function to estimate TAC for the DLMtool package
#'
#' 
#' @details This function creates harvest control rules (HCRs) which can be incorporated into a
#' management strategy evaluation framework (DLMtool package). HCRs are saved with a
#' generic name to the global environment and name of HCR is returned if results of the
#' functin are assigned to an object. HCR runs a SPiCT assessment using catch and
#' relative biomass index observations and stock status estimates are used to set the TAC
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
#' @param fractileBBmsy The fractile of the distribution of B/Bmsy. Default is 0.5 (median).
#' @param uncertaintyCap Logical; If true TAC is bound between two values set in lower and upper. Default: FALSE.
#' @param lower lower bound of the uncertainty cap. Default is 0.8, used if uncertaintyCap = TRUE.
#' @param upper upper bound of the uncertainty cap. Default is 1.2, used if uncertaintyCap = TRUE.
#' @return A function which can estimate TAC recommendations based on SPiCT assessment,
#'   taking assessment uncertainty into account.
#' 
#' @export
#'
#' @examples
#' \dontrun{
#' library(DLMtool)
#' 
#' ## Put together an operating model from the available DLM toolkit examples
#' StockEx <- Herring
#' FleetEx <- Generic_IncE
#' ObsEx <- Precise_Unbiased
#' 
#' ## Remove changes in life history parameters
#' StockEx@Mgrad <- c(0,0)
#' StockEx@Kgrad <- c(0,0)
#' StockEx@Linfgrad <- c(0,0)
#' StockEx@Prob_staying <- c(1,1)
#' 
#' ## Set the depletion level 
#' StockEx@D <- c(0.3, 0.4)
#'
#' ## create Operation Model
#' OMex <- new("OM", Stock = StockEx, Fleet = FleetEx, 
#'                   Obs = ObsEx)
#' 
#' ## Set simulation options
#' OMex@nsim <- 10
#' OMex@nyears <- 25
#' OMex@proyears <- 3
#' 
#' ## Get SPiCT HCR
#' MPname <- spict2DLMtool(fractileC=0.3)
#'
#' ## run MSE
#' MSEex <- runMSE(OMex, MPs = MPname,
#'             interval = 1, reps = 1, timelimit = 150, CheckMPs = FALSE)
#' }
#'
#'
spict2DLMtool <- function(fractileC = 0.5,
                          fractileFFmsy = 0.5,
                          fractileBBmsy = 0.5,
                          uncertaintyCap = FALSE,
                          lower = 0.8,
                          upper = 1.2){

    ## allowing for multiple generation of MPs
    argList <- list(fractileC, fractileFFmsy, fractileBBmsy,
                    uncertaintyCap, lower, upper)
    argLengths <- sapply(argList, length)
    maxi <- max(argLengths)
    maxl  <- which(argLengths == maxi)
    if(maxi>1){
        if(max(argLengths[(1:6)[-maxl]]) > 1) stop("Specified arguments have different lengths, they should have the same length or length = 1.")
    }
    argListCor <- argList
    argListCor[(1:6)[-maxl]] <- lapply(argList[(1:6)[-maxl]], function(x) rep(unlist(x), maxi))
    uncertaintyCapPrint <- argListCor[[4]]
    uncertaintyCapPrint[which(uncertaintyCapPrint == TRUE)] <- "T"
    uncertaintyCapPrint[which(uncertaintyCapPrint == FALSE)] <- "F"                        


    ## MP as function
    template  <- expression(paste0(
        'structure(function(x, Data, reps = 1,
          fractileC = ',a,',
          fractileFFmsy=',b,',
          fractileBBmsy=',c,',
          uncertaintyCap=',d,',
          lower=',e,',
          upper=',f,'){
            dependencies <- "Data@Year, Data@Cat, Data@Ind"
            time <- Data@Year
            Catch <- Data@Cat[x,]
            Index <- Data@Ind[x,]
            inp <- list(timeC=time, obsC=Catch, 
                        timeI=time, obsI=Index,
                        ## timepredc = max(time) + 1,
                        dteuler = 1 / 16,
                        do.sd.report=TRUE,
                        getReportCovariance = FALSE)
            rep <- try(spict::fit.spict(inp))
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
                ## Uncertainty cap
                if(uncertaintyCap){
                   red[red < lower] <- lower
                   red[red > upper] <- upper
                }
                predcatch <- try(spict::pred.catch(rep, get.sd = TRUE, exp = FALSE, fmsyfac = red))
                if(is(predcatch, "try-error")) {
                    TAC <- rep(NA, reps)
                } else {
                    TACi <- exp(qnorm(fractileC, predcatch[2], predcatch[4]))

                    ## Reduction based on B/Bmsy. Default = median
                    idx <- rep$inp$indpred[1]
                    logBBmsy <- spict::get.par("logBBmsy", rep, exp = TRUE)[idx,]
                    predBBtrigger <- 2 * exp(qnorm(fractileBBmsy, logBBmsy[2], logBBmsy[4]))
                    TACi <- TACi * min(1, predBBtrigger)

                    ## hack to guarantee compatibility with other MPs (DLMtool takes median, thus rep no effect)
                    TAC <- rep(TACi, reps)
                }
            }
            res <- DLMtool:::TACfilter(TAC)
            return(res)
        },
        class="Output")'))

        
    nami <- rep(NA,maxi)
    for(I in 1:maxi){

        ## create MPs as functions
        subList <- lapply(argListCor, "[[", I)
        names(subList) <- letters[1:6]
        templati <- eval(parse(text=paste(parse(text = eval(template, subList)),collapse=" ")))

         ## save names of MPs
        nami[I] <- paste0("spict_C",argListCor[[1]][I],"_FFmsy",
                       argListCor[[2]][I],"_BBmsy",argListCor[[3]][I],"_uC",uncertaintyCapPrint[I])
        assign(value=templati, x=nami[I], envir=globalenv())
    }

    ## allow for assigning names
    invisible(nami)
}


