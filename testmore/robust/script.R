## Tests vectorized robflagi and plotspict.priors for list priors

## Also tests using the logngamma prior
library(spict)


## Fleet=1: catches, Fleet > 1 : indices of exploitable biomass, Fleet=999 : effort
d<-read.csv("data.txt", skip=1, comment.char="#")

d$time_mid = (d$time_end + d$time_start)*0.5;
d$dt = (d$time_end - d$time_start)

eff = d[ d$fleet==999 , ]
d = subset( d, fleet!=999 )

dfleet=split(d,d$fleet)

getnamedEl<-function(x,name)x[[name]]

inp=list(obsC=dfleet[[1]]$obs,
         obsI=lapply(dfleet[-1],getnamedEl,name="obs"),
         timeC=dfleet[[1]]$time_start,
         timeI=lapply(dfleet[-1],getnamedEl,name="time_mid"),
         dtc=dfleet[[1]]$dt,
         stdevfacI=lapply(dfleet[-1],getnamedEl,name="stdevfac"),
         stdevfacC=dfleet[[1]]$stdevfac
         )



if(nrow(eff)>0){
    inp$obsE <- eff$obs
    inp$timeE <- eff$time_start
    inp$dte <- eff$dt
    inp$stdevfacE <- eff$stdevfac
}
## Numerical solver time step (probably don't need to change)
inp$dteuler <- 1/16

## Don't report covariances -- not needed and size of object is much smaller without
inp$getReportCovariance = FALSE
###########################
## Model configuration
###########################

## Priors
## (value, std.dev, 0/1=Off/On)
## logsdi may be given as a list of vectors, one for each index.
##inp$priors$logn <- c(log(2), 1, 1)
inp$priors$logalpha <- c(log(1), 2, 0)
inp$priors$logbeta <- c(log(1), 2, 0)
##inp$priors$logbkfrac <- c(log(1), 1, 1)
inp$priors$logsdi <- list( c(log(0.5), 1, 1) ,
                          c(log(0.5), 1, 1) )
inp$priors$logsdc <- c(log(0.1), 1, 1)
inp$robflagi <- c(0,1)
# Predict catches until this year
inp$timepredc <- max(inp$timeC,na.rm=TRUE) + 1

# Multiply last F with ffac and forecast using resulting F
inp$ffac <- 1.00

inp<-check.inp(inp)

inp$phases$logitpp <- -1

fit<-fit.spict(inp)

## Testing do.plot=1 only adds the first prior to existing plot
par(mfrow=c(2,1),mar=c(4,3,3,1))
plotspict.priors(fit,do.plot=1)
plotspict.priors(fit,do.plot=1)

## Should plot all active priors
plotspict.priors(fit)

## Gamma prior on logn with mode equal to Thorson's mean estimate and 90th percentile match
n.est <- 1.478
sdn <- 0.849
x90 <- qnorm(0.9,n.est,sdn)
sr <- modefrac2shaperate(log(n.est),log(x90))

inp$priors$logn <- c(log(2), 1, 0)
inp$priors$logngamma <- c(sr[1], sr[2], 1)

fit2<-fit.spict(inp)
plotspict.priors(fit2)
plot(fit2)

sink("res.out")
print(round(sumspict.parest(fit),3))
print(round(sumspict.parest(fit2),3))
sink()
