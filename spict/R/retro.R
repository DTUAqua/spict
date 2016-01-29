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
        #inpall[[i]]$ini <- as.list(rep$par.fixed)
        #inpall[[i]]$dtpredc <- 0
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
    asd <- try(parallel::mclapply(inpall, fit.spict))
    if(class(asd) == 'try-error'){
        rep$retro <- asd
    } else {
        rep$retro <- lapply(inpall, fit.spict)
    }
    return(rep)
}
