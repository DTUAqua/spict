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


#' @name fit.meyermillar
#' @title Fit the model of Meyer & Millar (1999)
#' @details Same input structure as for fit.spict(). Fitting the model of Meyer & Millar requires the packages rjags and coda. It furthermore requires that priors are specified for K, r, q, sigma2 (process error variance) and tau2 (observation error variance). Following Meyer & Millar (1999) the priors are:
#' \itemize{
#'  \item{"K"}{ log-normal.}
#'  \item{"r"}{ log-normal.}
#'  \item{"q"}{ inverse-gamma.}
#'  \item{"tau2"}{ inverse-gamma.}
#'  \item{"sigma2"}{ inverse-gamma.}
#' }
#' See example for how to specify priors.
#' @param mminp Input list similar to the input to fit.spict()
#' @return List containing results
#' @examples
#' priors <- list()
#' priors$K <- c(5.042905, 3.76)
#' priors$r <- c(-1.38, 3.845)
#' priors$iq <- c(0.001, 0.0012)
#' priors$itau2 <- c(1.709, 0.00861342)
#' priors$isigma2 <- c(3.785518, 0.0102232)
#' priors$logPini <- -0.223
#' data(pol)
#' inp <- pol$albacore
#' inp$meyermillar$n.iter <- 10000
#' inp$meyermillar$burnin <- 1000
#' inp$meyermillar$thin <- 10
#' inp$meyermillar$n.chains <- 1
#' inp$meyermillar$priors <- priors
#' res <- fit.meyermillar(inp)
#' summary(res$jags)
#' @export
fit.meyermillar <- function(mminp){
    mminp <- check.inp(mminp)
    check.prior <- function(nm){
        if (!nm %in% names(mminp$meyermillar$priors)){
            stop(paste('Prior for', nm, 'is missing!'))
        } else {
            stopifnot(length(mminp$meyermillar$priors[[nm]]) > 1)
        }
    }
    # Check priors
    if (!"priors" %in% names(mminp$meyermillar)){
        stop('Cannot continue because no priors were specified for Meyer & Millar model, see example of fit.meyermillar for how to do that.')
    } else {
        check.prior('K')
        check.prior('r')
        check.prior('iq')
        check.prior('isigma2')
        check.prior('itau2')
        if (!'logPini' %in% names(mminp$meyermillar$priors)){
            stop(paste('Prior for logPini is missing!'))
        }
    }
    
    # Write tmeporary .bug file
    bugfn <- mminp$meyermillar$bugfn
    write.bug.file(mminp$meyermillar$priors, fn=bugfn)

    # Fit model
    res <- list()
    res$inp <- mminp
    #n.iter <- 225000
    res$jags <- fit.jags(mminp,
                         fn=bugfn,
                         n.iter=mminp$meyermillar$n.iter,
                         n.chains=mminp$meyermillar$n.chains,
                         burnin=mminp$meyermillar$burnin,
                         thin=mminp$meyermillar$thin)

    # Diagnostics
    alpha <- 0.05 # Significance level
    res$diag$convall <- numeric(mminp$meyermillar$n.chains)
    res$diag$pfail <- numeric(mminp$meyermillar$n.chains)
    # Heidelberger and Welch stationarity and half-width test
    hd <- coda::heidel.diag(res$jags, pvalue=alpha)
    res$diag$heidel.diag <- hd
    res$diag$hd.no.failed <- numeric(mminp$meyermillar$n.chains)
    # Geweke Z-score for convergence
    gd <- coda::geweke.diag(res$jags)
    res$diag$geweke.diag <- gd
    res$diag$gd.no.failed <- numeric(mminp$meyermillar$n.chains)
    # Raftery and Lewis check thinning, burnin, and sample size
    rd <- coda::raftery.diag(res$jags, q=0.025, r=0.02, s=0.9)
    res$diag$raftery.diag <- rd
    for (chainno in 1:mminp$meyermillar$n.chains){
        res$diag$hd.no.failed[chainno] <- sum(hd[[chainno]][, c(1, 4)]==0, na.rm=TRUE)
        res$diag$gd.no.failed[chainno] <- sum(pnorm(gd[[chainno]]$z) < alpha)
        # Determine overall convergence of chain
        cn <- chainno # Chain number
        ntests <- dim(res$diag$heidel.diag[[cn]])[1] + length(res$diag$geweke.diag[[cn]]$z)
        no.failed <- res$diag$hd.no.failed[chainno] + res$diag$gd.no.failed[chainno]
        res$diag$pfail[chainno] <- no.failed/ntests
        # convall = 0 means converged
        res$diag$convall[chainno] <- as.numeric(res$diag$pfail[chainno] > alpha) 
    }

    # Save a summary of results
    sum <- summary(res$jags)
    resmat <- cbind(sum$statistics[, 1:2], sum$quantiles[, c(1, 5)],
                    median=sum$quantiles[, 3])
    resmat <- resmat[-grep('P', rownames(resmat)), ]
    resmat <- cbind(resmat, CV=resmat[, 2]/resmat[, 1])
    res$resmat <- resmat
    
    # Clean up temporary files
    if (mminp$meyermillar$cleanup){
        unlink(bugfn)
    }
    return(res)
}


#' @name plotmm.priors
#' @title Plot priors of Meyer & Millar model
#' @param nm Name of prior
#' @param priorsin List of priors, typically inp$meyermillar$priors.
#' @param add If TRUE add to current plot.
#' @param ... Additional arguments to plot.
#' @return Nothing.
#' @export
plot.priors <- function(nm, priorsin, add=TRUE, ...){
    if (nm == 'K'){
        mu <- priorsin$K[1]
        sd <- 1/sqrt(priorsin$K[2])
        vec <- exp(seq(mu-3*sd, mu+3*sd, length=100))
        val <- dlnorm(vec, mu, sd)
    }
    if (nm == 'r'){    
        mu <- priorsin$r[1]
        sd <- 1/sqrt(priorsin$r[2])
        vec <- exp(seq(mu-3*sd, mu+3*sd, length=100))
        #val <- dnorm(log(vec), mu, sd)
        val <- dlnorm(vec, mu, sd)
    }
    if (nm == 'iq'){
        #require(pscl) # for densigamma
        qshape <- priorsin$iq[1]
        qrate <- priorsin$iq[2]
        vec <- seq(0.1, 100, length=1000)
        #val <- densigamma(vec, qshape, qrate)
        val <- dgamma(vec, shape=qshape, rate=qrate)
        #vec <- 1/vec
    }
    if (nm == 'isigma2'){
        #require(pscl) # for densigamma
        sigmashape <- priorsin$isigma2[1]
        sigmarate <- priorsin$isigma2[2]
        vec <- seq(1, 500, length=200)
        #vec <- seq(1e-6, 0.02, length=200)
        #val <- densigamma(vec, sigmashape, sigmarate)
        val <- dgamma(vec, sigmashape, sigmarate)
    }
    if (nm == 'itau2'){
        #require(pscl) # for densigamma
        taushape <- priorsin$itau2[1]
        taurate <- priorsin$itau2[2]
        vec <- seq(1, 500, length=200)
        #vec <- seq(1e-6, 0.05, length=200)
        #val <- densigamma(vec, taushape, taurate)
        val <- dgamma(vec, taushape, taurate)
    }
    if (add){
        lines(vec, val, ...)
    } else {
        plot(vec, val, typ='l', ...)
    }
}


#' @name write.bug.file
#' @title Write the BUGS code to a text file
#' @details The .bug file generated by this function contains code published in Meyer & Millar (1999).
#' @param priors List of priors, typically coming from inp$meyermillar$priors.
#' @param fn Filename of to put BUGS code in.
#' @references Meyer, R., & Millar, R. B. (1999). BUGS in Bayesian stock assessments. Canadian Journal of Fisheries and Aquatic Sciences, 56(6), 1078-1087.
#' @return Nothing.
write.bug.file <- function(priors, fn='tmp.bug'){
    cat('model{\n# prior distribution of K: lognormal with 10% and 90% quantile at 80 and 300\n', file=fn)
    cat(paste0('K ~ dlnorm(', priors$K[1], ',', priors$K[2], ')I(10,1000);\n'),
        file=fn, append=TRUE)
    cat('# prior distribution of r: lognormal with 10% and 90% quantile at 0.13 and 0.48\n',
        file=fn, append=TRUE)
    cat(paste0('r ~ dlnorm(', priors$r[1], ',', priors$r[2], ')I(0.01,1.2);\n'),
        file=fn, append=TRUE)
    cat('# prior distribution of q: instead of improper (prop. to 1/q) use just proper IG\n',
        file=fn, append=TRUE)
    cat(paste0('iq ~ dgamma(', priors$iq[1], ',', priors$iq[2], ')I(0.5,100);\nq <- 1/iq;\n'),
        file=fn, append=TRUE)
    cat('# prior distribution of sigma2: inv. gamma with 10% and 90% qu. at 0.04 and 0.08\n',
        file=fn, append=TRUE)
    cat(paste0('isigma2 ~ dgamma(', priors$isigma2[1], ',', priors$isigma2[2],
               ');\nsigma2 <- 1/isigma2;\n'), file=fn, append=TRUE)
    cat('# prior distribution of tau2: inv. gamma with 10% and 90% qu. at 0.05 and 0.15\n',
        file=fn, append=TRUE)
    cat(paste0('itau2 ~ dgamma(', priors$itau2[1], ',', priors$itau2[2],
               ');\ntau2 <- 1/itau2;\n'), file=fn, append=TRUE)
    cat('# (conditional) prior distribution of Ps (from state equations):\n',
        file=fn, append=TRUE)
    cat(paste0('Pmed[1] <- ', priors$logPini, ';\n'), file=fn, append=TRUE)
    cat('P[1] ~ dlnorm(Pmed[1],isigma2) #I(0.001,2.0)
for (t in 2:N) { 
    Pmed[t] <- log(max(P[t-1] + r*P[t-1]*(1-P[t-1]) - C[t-1]/K, 0.001));
    P[t] ~ dlnorm(Pmed[t],isigma2) #I(0.001,2.0) 
}
# sampling distribution:
for (t in 1:N) { 
    Imed[t] <- log(q*K*P[t]);
    I[t] ~ dlnorm(Imed[t],itau2);
}
# further management parameters and predictions:
MSY <- r*K/4;
EMSY <- r/(2*q);
Bmsy <- K/2;
Fmsy <- r/2;
Blast <- P[N]*K;
BBlast <- P[N]*2;
Ppred <- P[N] + r*P[N]*(1-P[N]) - C[N]/K;
Bpred <- Ppred*K;
BBpred <- Ppred*2;
}', file=fn, append=TRUE)
}


#' @name fit.jags
#' @title Fit the Meyer & Millar model using rjags
#' @param inp Input list containing data and settings.
#' @param fn Filename of containing BUGS code.
#' @param n.iter Number of iterations.
#' @param n.chains Number of chains.
#' @param burnin Number of burn-in iterations.
#' @param thin Thin chains by this value.
#' @return The raw output of rjags::coda.samples.
fit.jags <- function(inp, fn, n.iter=10000, n.chains=1, burnin=round(n.iter/2), thin=1000){
    N <- inp$nobsC
    datain <- list(C=inp$obsC,
                   I=inp$obsI[[1]],
                   N=N)
    initsin <- list(P=rep(0.5, N),
                    r=exp(inp$ini$logr),
                    K=exp(inp$ini$logK),
                    iq=1/exp(inp$ini$logq),
                    isigma2=1/exp(inp$ini$logsdb)^2,
                    itau2=1/exp(inp$ini$logsdi)^2)
    # Parameters
    pars <- c('P', 'r', 'K', 'iq', 'isigma2', 'itau2', 'MSY', 'EMSY', 'Bmsy', 'Fmsy',
              'Blast', 'BBlast', 'Ppred', 'Bpred', 'BBpred')
    # Number of iteration
    jags <- rjags::jags.model(file=fn, data=datain, inits=initsin, n.chains=n.chains)
    if (burnin > 0){
        update(jags, burnin) # Burnin
    }
    samples <- rjags::coda.samples(jags, pars, n.iter=n.iter, thin=thin)
    return(samples)
}


