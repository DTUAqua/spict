#' @name spict_percentiles
#' @title SPiCT assessment with specified percentiles from the predicted catch distribution
#'   and/or the distribution of Fmsy
#' 
#' @details SPiCT assessment is done using catch and relative biomass index observations. 
#' Stock status estimates are used to set the TAC for the next year, equal to the specified
#' percentile of the distribution of predicted catches and/or the distribution of
#' the Fmsy reference level
#'
#' @param x A position in a data-limited mehods data object
#' @param DLM_data A data-limited methods data object (see DLMtool)
#' @param reps The number of stochastic samples of the TAC recommendation
#' @param percentileC The percentile of the catch distribution to be used for setting TAC
#' @param percentileFmsy The percentile of the distribution of Fmsy
#' @param cap Logical; If true TAC is multiplied with 1 or the ratio of current biomass
#'   over Blim (0.5*Bmsy) if this ratio is smaller than 1. Default is FALSE.
#'
#' @return A numeric vector of TAC recommendations (one TAC and a number of NAs corresponding to
#'   reps - 1).
#' @export
#'
#' @examples
#' \dontrun{
#' library(DLMtool)
#' 
#' ## Put together an operating model from the available DLM toolkit examples
#' stock <- Herring
#' Fleet.example <- Generic_IncE
#' Observation.example <- Precise_Unbiased
#' 
#' ## Remove changes in life history parameters
#' stock@Mgrad <- c(0,0)
#' stock@Kgrad <- c(0,0)
#' stock@Linfgrad <- c(0,0)
#' stock@Prob_staying <- c(1,1)
#' 
#' ## Set the depletion level 
#' stock@D <- c(0.3, 0.4)
#'
#' OM.example <- new("OM", Stock = stock, Fleet = Fleet.example, 
#'                   Obs = Observation.example)
#' OM.example@nsim <- 20
#' OM.example@nyears <- 25
#' OM.example@proyears <- 5
#'
#' MP.vec <- c("spict_percentiles")
#'
#' MSE.example <- runMSE(OM.example, MPs = MP.vec,
#'                       interval = 1, reps = 1, timelimit = 150, CheckMPs = FALSE)
#' }
#'

spict_percentiles <- structure(
    function(x, Data, reps = 1, percentileC = 0.50, percentileFmsy=NA, cap=FALSE){
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
                    blim <- 0.5 * bmsy
                    capi <- min(1, blast/blim)
                } else {
                    capi <- 1
                }
                TACi <- TACi * capi
                TAC <- c(TACi, rep(TACi, reps-1))                
            }
        }
        rm(rep); gc()
        res <- DLMtool:::TACfilter(TAC)
        return(res)
    },
    class = "Output"
)

