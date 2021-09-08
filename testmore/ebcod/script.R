library(spict)


inp<-list()

datQ <- read.csv("EBcod catch by Quarter.csv",sep=",",stringsAsFactors=F)
catchesQonly <- data.frame(Time=c(datQ$X,datQ$X+1/4,datQ$X+2/4,datQ$X+3/4),
                      catch=c(datQ$Q1,datQ$Q2,datQ$Q3,datQ$Q4)/1e3 )
catchesQonly<-catchesQonly[order(catchesQonly$Time),]

inp$timeC <- catchesQonly$Time
inp$obsC <- catchesQonly$catch


inp$catchunit <- "'000 t"

## BITS
lbi<-read.csv("ebcodindexsoft24hard.csv")

midL = read.table("midL.dat")$x
LWparam = read.csv("lw parameters_EBC.csv")

lbi$Year = as.numeric(sapply(lbi$X,substr,start=0,stop=4))
lbi$Quarter = as.numeric(sapply(lbi$X,substr,start=6,stop=6))
lbi$Time = lbi$Year+ (lbi$Quarter-1)/4+1/8

## Exploitable stock biomass correction here
esbcorr=read.table("ESBcorr.dat")
esbcorr = rbind( data.frame(length=1:8,corr=esbcorr[1,2]),
                 esbcorr,
                data.frame(length=91:120,corr=esbcorr[nrow(esbcorr),2]))

## get it at the right lengths
esb2 = approx(esbcorr$length,esbcorr$corr,xout=midL)

lbi$ESB = numeric( nrow(lbi) )
lbi$B = numeric( nrow(lbi) )
for(i in 1:nrow(lbi)){
    lw = LWparam[ LWparam$year==lbi$Year[i],3:4]
    for(l in 1:length(midL)){
        lbi[i,"ESB"] = lbi[i,"ESB"] + esb2$y[l]*lbi[i,l+1]*lw[1]*midL[l]^lw[2]
        lbi[i,"B"] = lbi[i,"B"] + lbi[i,l+1]*lw[1]*midL[l]^lw[2]
    }
}

mstd<-function(x) x/mean(x,na.rm=TRUE)

selQ1<-which(lbi$Quarter==1)
inp$timeI[[1]] <- lbi$Time[selQ1]
inp$obsI[[1]] <- mstd(lbi$ESB[selQ1])

selQ4<-which(lbi$Quarter==4)
inp$timeI[[2]] <- lbi$Time[selQ4]
inp$obsI[[2]] <-  mstd(lbi$ESB[selQ4])

inp$stdevfacI[[1]] = mstd(lbi$CVtotal[selQ1])
inp$stdevfacI[[2]] = mstd(lbi$CVtotal[selQ4])

## egg survey
eggsurv<-read.csv("SSB from egg production.csv")
eggsurv = eggsurv[ eggsurv$Year>=1991 ,]

inp$timeI[[3]] = eggsurv$Year + (eggsurv$Month-1)/12
inp$obsI[[3]] = eggsurv$SSB

inp$stdevfacI[[3]] = rep(1,length(inp$obsI[[3]]))

inp$stdevfacC = rep(1, length(inp$timeC))
inp$stdevfacC[ inp$timeC>=1993 & inp$timeC<1996  ] = 2
inp$stdevfacC[ inp$timeC>=2000 & inp$timeC<2008  ] = 2

inp$timepredc<-2019 ## forecast year start
inp$dtpredc<-1 ## one year ahead
inp$timepredi<-2020 ## forecast year end
inp$manstart<-2019 ## forecast year start

inporig<-inp

add.thorson.prior<-function(x, family="Pooled"){
    n.est <- 1.478
    sdn <- 0.849
    if(family=="Pleuronectiformes") { n.est <- 1.353; sdn <- 0.739 }
    if(family=="Gadiformes") { n.est <- 1.729; sdn <- 0.937 }
    if(family=="Perciformes") { n.est <- 1.064; sdn <- 0.590 }
    if(family=="Clupeiformes") { n.est <- 0.599; sdn <- 0.342 }
    if(family=="Scorpaeniformes") { n.est <- 1.970; sdn <- 1.074 }
    if(family=="Other") { n.est <- 1.431; sdn <- 0.805 }

    logn.est <- log( n.est)
    var.logn <- sdn^2 / n.est^2
    x$priors$logn <- c( logn.est, sqrt( var.logn) )
    x
}

AIC.spictcls<-function(x) 2*x$opt$objective + 2*length(x$opt$par)

#####################################################
##  C - Vg  --  Thorson prior
#####################################################

inp<-inporig

inp$seasontype=3

inp$timevaryinggrowth=TRUE ## only works if set before check.inp (!?)
inp$msytype="d" ## n is estimated to be < 1, i.e. stochastic ref points do not work

inp<-check.inp(inp)

inp$ini$logq[3]=0 ## helps convergence

## OBS OBS, horrible defaults of -18.42 ?!?
inp$ini$logsdm=-3
inp$ini$logpsi=-3


inp <- add.thorson.prior(inp,family="Gadiformes")
inp$priors$logn
inp$priors$logalpha=c(0,0,0)
inp$priors$logbeta=c(0,0,0)
inp$priors$logsdm=c(0,0,0)
inp$priors$logpsi=c(0,0,0)

res9 <- fit.spict(inp, dbg=0)

retr9<-retro(res9,5,mc.cores=1)

mvec2<-get.par("logMSYvec",res9,exp=TRUE)
##plot(res9$inp$time,mvec2[,2],type="l",lwd=2,ylim=c(0,120),xlab="Time",ylab="m")
for(i in 1:length(retr9$retro)){
    mvec2<-get.par("logMSYvec",retr9$retro[[i]],exp=TRUE)
    tt <- retr9$retro[[i]]$inp$time
    ##lines(tt,mvec2[,2])
}

cat(round(retr9$retro[[i]]$opt$objective,3),"\n", file="res.out")
