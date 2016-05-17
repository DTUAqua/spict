
fit.meyermillar <- function(mminp){
    
    # Write tmeporary .bug file
    bugfn <- mminp$meyermillar$bugfn
    write.bug.file(fn=bugfn, mminp$meyermillar$priors)

    # Fit model
    niter <- 225000
    s <- fit.jags(inp, fn=bugfn, niter, n.chains=2, burnin=25000, thin=25)




    parnms <- c('K', 'r', 'iq', 'isigma2', 'itau2')
    nparnms <- length(parnms)

    sum <- summary(s)
    nms <- rownames(sum$quantiles)
    inds <- match(parnms, nms)

    # Diagnostics
    chainno <- 1
    # Heidelberger and Welch stationarity and half-width test
    hd <- heidel.diag(s)
    hd.no.failed <- sum(hd[[chainno]][, c(1, 4)]==0, na.rm=TRUE)
    # Geweke Z-score for convergence
    gd <- geweke.diag(s)
    gd.no.failed <- sum(pnorm(gd[[chainno]]$z) < 0.05)
    # Raftery and Lewis check thinning, burnin, and sample size
    stest <- fit.jags(inp, fn=bugfn, niter=100000, n.chains=1, burnin=2, thin=1)
    rd <- raftery.diag(stest, q=0.025, r=0.02, s=0.9)

}


plot.priors <- function(nm, priorsin, add=TRUE, ...){
    if (nm == 'K'){
        mu <- priorsin$K[1]
        sd <- 1/sqrt(priorsin$K[2])
        vec <- exp(seq(mu-3*sd, mu+3*sd, length=100))
        val <- dnorm(log(vec), mu, sd)
    }
    if (nm == 'r'){    
        mu <- priorsin$r[1]
        sd <- 1/sqrt(priorsin$r[2])
        vec <- exp(seq(mu-3*sd, mu+3*sd, length=100))
        val <- dnorm(log(vec), mu, sd)
    }
    if (nm == 'isigma2'){
        require(pscl) # for densigamma
        sigmashape <- priorsin$isigma2[1]
        sigmarate <- priorsin$isigma2[2]
        vec <- seq(1e-6, 0.02, length=200)
        val <- densigamma(vec, sigmashape, sigmarate)
    }
    if (nm == 'itau2'){
        require(pscl) # for densigamma
        taushape <- priorsin$itau2[1]
        taurate <- priorsin$itau2[2]
        vec <- seq(1e-6, 0.05, length=200)
        val <- densigamma(vec, taushape, taurate)
    }
    if (nm == 'iq'){
        require(pscl) # for densigamma
        qshape <- priorsin$iq[1]
        qrate <- priorsin$iq[2]
        vec <- seq(1e-6, 0.6, length=200)
        val <- densigamma(vec, qshape, qrate)
    }
    if (add){
        lines(vec, val, ...)
    } else {
        plot(vec, val, typ='l', ...)
    }
}



write.bug.file <- function(fn='tmp.bug', priors){
    cat('model{\n# prior distribution of K: lognormal with 10% and 90% quantile at 80 and 300\n', file=fn)
    cat(paste0('K ~ dlnorm(', priors$K[1], ',', priors$K[2], ')I(10,1000);\n'), file=fn, append=TRUE)
    cat('# prior distribution of r: lognormal with 10% and 90% quantile at 0.13 and 0.48\n', file=fn, append=TRUE)
    cat(paste0('r ~ dlnorm(', priors$r[1], ',', priors$r[2], ')I(0.01,1.2);\n'), file=fn, append=TRUE)
    cat('# prior distribution of q: instead of improper (prop. to 1/q) use just proper IG\n', file=fn, append=TRUE)
    cat(paste0('iq ~ dgamma(', priors$iq[1], ',', priors$iq[2], ')I(0.5,100);\nq <- 1/iq;\n'), file=fn, append=TRUE)
    cat('# prior distribution of sigma2: inv. gamma with 10% and 90% qu. at 0.04 and 0.08\n', file=fn, append=TRUE)
    cat(paste0('isigma2 ~ dgamma(', priors$isigma2[1], ',', priors$isigma2[2], ');\nsigma2 <- 1/isigma2;\n'), file=fn, append=TRUE)
    cat('# prior distribution of tau2: inv. gamma with 10% and 90% qu. at 0.05 and 0.15\n', file=fn, append=TRUE)
    cat(paste0('itau2 ~ dgamma(', priors$itau2[1], ',', priors$itau2[2], ');\ntau2 <- 1/itau2;\n'), file=fn, append=TRUE)
    cat('# (conditional) prior distribution of Ps (from state equations):
Pmed[1] <- 0;
P[1] ~ dlnorm(Pmed[1],isigma2) #I(0.001,2.0)
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
MSP <- r*K/4;
EMSP <- r/(2*q);
P1990 <- P[N]+r*P[N]*(1-P[N]) -C[N]/K;
B1990 <- P1990*K;
}', file=fn, append=TRUE)
}


fit.jags <- function(inp, fn, niter, n.chains=1, burnin=round(niter/2), thin=1000){
    N <- inp$nobsC
    datain <- list(C=inp$obsC, I=inp$obsI[[1]], N=N)
    initsin <- list(P=rep(0.5, N), r=exp(inp$ini$logr), K=exp(inp$ini$logK), iq=1/exp(inp$ini$logq),
                    isigma2=1/exp(inp$ini$logsdb)^2, itau2=1/exp(inp$ini$logsdi)^2)
    # Parameters
    pars <- c('P', 'r', 'K', 'iq', 'isigma2', 'itau2')
    # Number of iteration
    jags <- jags.model(file=fn, data=datain, inits=initsin, n.chains=n.chains)
    if (burnin > 0){
        update(jags, burnin) # Burnin
    }
    samples <- coda.samples(jags, pars, n.iter=niter, thin=thin)
    return(samples)
}


