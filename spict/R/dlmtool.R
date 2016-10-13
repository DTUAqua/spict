#' @name SPiCT_Feq08Fmsy
#' @title SPiCT assessment with F equal 80\% Fmsy harvest control rule
#' 
#' @details SPiCT assessment is done using catch and relative biomass index observations. 
#' Stock status estimates are used to set the TAC for the next year, equal to the
#' catch that corresponds to fishing mortality equal to 80\% of Fmsy.
#'
#' @param x A position in a data-limited mehods data object
#' @param DLM_data A data-limited methods data object (see DLMtool)
#' @param reps The number of stochastic samples of the TAC recommendation
#'
#' @return A numeric vector of TAC recommendations
#' @export
#'
#' @examples
#' \dontrun{
#' library(DLMtool)
#' 
#' ## Put together an operating model from the available DLM toolkit examples
#' stock <- DLMdat[[6]] ## Herring
#' Fleet.example <- DLMdat[[22]] # Generic_IncE
#' Observation.example <- DLMdat[[34]] # Precise_Unbiased
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
#'                   Observation = Observation.example)
#'
#' MP.vec <- c("SPiCT_Feq08Fmsy")
#' 
#' MSE.example <- runMSE(OM.example, MPs = MP.vec, nsim = 200, proyears = 20,
#'                       interval = 1, reps = 100, timelimit = 150, CheckMPs = FALSE)
#' }
SPiCT_Feq08Fmsy <- structure(function(x, DLM_data, reps) {
  dependencies <- "DLM_data@Year DLM_data@Cat DLM_data@Ind"
  time <- DLM_data@Year
  Catch <- DLM_data@Cat[x, ]
  Index <- DLM_data@Ind[x, ]
  inp <- list(timeC=time, obsC=Catch, 
              timeI=time, obsI=Index,
              ## timepredc = max(time) + 1,
              dteuler = 1 / 4,
              do.sd.report=FALSE)
  rep <- try(spict::fit.spict(inp))
  if(is(rep, "try-error")) {
    TAC <- rep(NA, reps)
  } else {
    #predcatch <- spict::pred.catch(rep, get.sd = TRUE, exp = FALSE, fmsyfac = 0.8)
    predcatch <- spict::pred.catch(rep, get.sd = FALSE, exp = FALSE, fmsyfac = 0.8)
    ## The repetitions of TAC are log normally distributed with CV of 10% to match
    ## the implementation of DD management procedure. To accuratelly represent the
    ## uncertainty estimated by SPiCT one should use the estimated standard deviation
    ## reported by SPiCT, i.e. TAC <- exp(rnorm(reps, predcatch[2], predcatch[4]))
    ## remember to set do.sd.report = TRUE when creating inp
    TAC <- exp(rnorm(reps, predcatch[2], predcatch[2] * 0.1))
  }
  DLMtool:::TACfilter(TAC)
}, class = "DLM_output")


