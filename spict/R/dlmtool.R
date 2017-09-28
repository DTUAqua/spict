#' @name spict_catchPercentile
#' @title SPiCT assessment with specified percentiles from the predicted catch distribution
#' 
#' @details SPiCT assessment is done using catch and relative biomass index observations. 
#' Stock status estimates are used to set the TAC for the next year, equal to the specified
#' percentile of the distribution of predicted catches.
#'
#' @param x A position in a data-limited mehods data object
#' @param DLM_data A data-limited methods data object (see DLMtool)
#' @param reps The number of stochastic samples of the TAC recommendation
#' @param percentile The percentile of the catch distribution to be used for setting TAC
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
#' MP.vec <- c("spict_catchPercentile")
#'
#' MSE.example <- runMSE(OM.example, MPs = MP.vec,
#'                       interval = 1, reps = 1, timelimit = 150, CheckMPs = FALSE)
#' }
#'
spict_catchPercentile <- structure(
    function(x, Data, reps = 1, percentile = 0.5, cap = FALSE){
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
        if(is(rep, "try-error") | rep$opt$convergence != 0) {
            TAC <- rep(NA, 1)
        } else {
            predcatch <- try(spict::pred.catch(rep, get.sd = TRUE, exp = FALSE, fmsyfac = 1))
            if(is(predcatch, "try-error")) {
                TAC <- rep(NA, 1)
            } else {
                TAC <- exp(qnorm(percentile, predcatch[2], predcatch[4]))
                if(cap){
                    ## cap
                    bio <- spict::get.par("logB", rep, exp = TRUE)
                    bioLast <- bio[nrow(bio),]
                    Bmsy <- spict::get.par("logBmsy", rep, exp = TRUE)
                    Blim <- 0.5 * Bmsy
                    cap <- min(1, bioLast[2]/Blim[2])
                } else {
                    cap <- 1
                }
                ## correct TAC
                TAC <- TAC * cap
            }
        }
        rm(rep); gc()
        tacTemp <- DLMtool:::TACfilter(TAC)
        res <- c(tacTemp, rep(NA, reps-1))
        return(res)
    },
    class = "Output"
)
