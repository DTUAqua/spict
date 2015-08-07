#' @name retro
#' @title Conduct retrospective analysis 
#' @details A retrospective analysis consists of estimating the model with later data points removed sequentially one year at a time.
#' @param rep A valid result from fit.spict.
#' @param nretroyear Number of years of data to remove (this is also the total number of model runs).
#' @return A rep list with the added key retro containing the results of the retrospective analysis. Use plotspict.retro() to plot these results.
#' @examples
#' data(pol)
#' inp <- pol$albacore
#' rep <- fit.spict(inp)
#' rep <- retro(rep, nretroyear=6)
#' plotspict.retro(rep)
#' @export
retro <- function(rep, nretroyear=5){
    inp1 <- rep$inp
    inpall <- list()
    for(i in 1:nretroyear){
        inpall[[i]] <- list()
        inpall[[i]]$ini <- as.list(rep$par.fixed)
        inpall[[i]]$dtpredc <- 0
        indsC <- which(inp1$timeC <= inp1$timeC[inp1$nobsC]-i)
        inpall[[i]]$obsC <- inp1$obsC[indsC]
        inpall[[i]]$timeC <- inp1$timeC[indsC]
        inpall[[i]]$obsI <- list()
        inpall[[i]]$timeI <- list()
        for(j in 1:inp1$nindex){
            indsI <- which(inp1$timeI[[j]] <= inp1$timeI[[j]][inp1$nobsI[j]]-i)
            inpall[[i]]$obsI[[j]] <- inp1$obsI[[j]][indsI]
            inpall[[i]]$timeI[[j]] <- inp1$timeI[[j]][indsI]
        }
    }
    asd <- try(library(parallel))
    if(class(asd) == 'try-error'){
        rep$retro <- mclapply(inpall, fit.spict)
    } else {
        rep$retro <- lapply(inpall, fit.spict)
    }
    return(rep)
}
