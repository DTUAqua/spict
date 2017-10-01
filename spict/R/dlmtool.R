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
#' @param percentileC The percentile of the catch distribution to be used for setting TAC. Default
#'   is median (0.5).
#' @param percentileFmsy The percentile of the distribution of Fmsy. Default is NA.
#' @param cap Logical; If true TAC is multiplied with 1 or the ratio of current biomass
#'   over Blim (0.5*Bmsy) if this ratio is smaller than 1. Default is FALSE.
#'
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
#' MPname <- spict2DLMtool(percentileC=0.3)
#'
#' ## run MSE
#' MSEex <- runMSE(OMex, MPs = MPname,
#'             interval = 1, reps = 1, timelimit = 150, CheckMPs = FALSE)
#' }
#'
#'

spict2DLMtool <- function(percentileC = 0.5, percentileFmsy = NA, cap = TRUE){

    maxl <- max(length(percentileC),length(percentileFmsy),length(cap))

    ## error messages if lengths are different and not 1
    if(length(percentileC) != maxl) percentileC <- rep(percentileC, maxl)
    if(length(percentileFmsy) != maxl) percentileFmsy <- rep(percentileFmsy, maxl)
    if(length(cap) != maxl) cap <- rep(cap, maxl)

    X  <- expression(paste0(
        'structure(function(x, Data, reps = 1,
          percentileC = ',xI,',
          percentileFmsy=',yI,', cap=',zI,'){
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
                if(!is.na(percentileFmsy) & !is.null(percentileFmsy)){
                    idx <- rep$inp$indpred[1]
                    logFFmsy <- spict::get.par("logFFmsy", rep)[idx,]
                    fi <- 1-percentileFmsy
                    fm <- exp( qnorm( fi, logFFmsy[2], logFFmsy[4] ) )
                    fm5 <- exp( qnorm( 0.5, logFFmsy[2], logFFmsy[4] ) )
                    red <- fm5 / fm
                } else {
                    red <- 1
                }
                    predcatch <- try(spict::pred.catch(rep, get.sd = TRUE, exp = FALSE, fmsyfac = red))
                if(is(predcatch, "try-error")) {
                    TAC <- rep(NA, reps)
                } else {
                    TACi <- exp(qnorm(percentileC, predcatch[2], predcatch[4]))
                    if(cap){
                        idx <- rep$inp$indpred[1]                    
                        blast <- spict::get.par("logB", rep, exp = TRUE)[idx,2]
                        bmsy <- spict::get.par("logBmsy", rep, exp = TRUE)[2]
                        btrigger <- 0.5 * bmsy
                        capi <- min(1, blast/btrigger)
                    } else {
                        capi <- 1
                    }
                    TACi <- TACi * capi
                    TAC <- rep(TACi, reps)
                }
            }
            rm(rep); gc()
            res <- DLMtool:::TACfilter(TAC)
            return(res)
        },
        class="Output")'))    

    wrapper <- function (...) {eval(parse(text=paste(...,collapse=" ")))}
    nami <- rep(NA,maxl)
    for(I in 1:maxl){
        ##X <- structure(X, class = "Output")
        XI <- eval(X, list(xI=percentileC[I], yI=percentileFmsy[I], zI=cap[I]))
        XI <- parse(text=XI)
        XI <- wrapper(XI)
        nami[I] <- paste0("spict_percentile_C",percentileC[I],"_Fmsy",
                       percentileFmsy[I],"_Cap",cap[I])
        assign(value=XI,x=nami[I],envir=globalenv())
        rm(XI)
    }
    rm(X);gc()
    invisible(nami)
}



