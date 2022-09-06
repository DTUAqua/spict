fn <- 'spict/DESCRIPTION'
date <- format(Sys.Date(), "%Y-%m-%d")
cat('Package: spict
Type: Package
Title: Stochastic surplus Production model in Continuous-Time (SPiCT)
Version: 1.3.7
Date:', date, '
Authors@R: c(person(given="Martin Waever",
                    family="Pedersen",
                    email="wpsgodd@gmail.com",
                    role=c("aut", "cre", "cph")),
             person(given="Casper Willestofte",
                    family="Berg",
                    email="cbe@aqua.dtu.dk",
                    role=c("aut"),
                    comment=c()),
             person(given="Tobias Karl",
                    family="Mildenberger",
                    email="tobm@aqua.dtu.dk",
                    role=c("aut"),
                    comment=c(ORCID="0000-0002-6631-7524")),
             person(given="Alexandros",
                    family="Kokkalis",
                    email="alko@aqua.dtu.dk",
                    role=c("aut"),
                    comment=c(ORCID="0000-0002-9662-6272")))
Description: Fits a surplus production model to fisheries catch and biomass index data.
License: GPL (>=3)
Depends:
    R (>= 3.0),
    TMB (>= 1.7.1)
LinkingTo: TMB, RcppEigen
Imports:
    ellipse
Suggests:
    parallel,
    mgcv,
    rjags,
    coda,
    knitr,
    rmarkdown
LazyData: true
Encoding: UTF-8
VignetteBuilder: knitr\n',
    file=fn)

# Add Git hub stuff
sha <- system('git rev-parse HEAD', intern=TRUE)
branch <- system('git rev-parse --abbrev-ref HEAD', intern=TRUE)

##cat(paste('GithubRepo: spict\n'), file=fn, append=TRUE)
##cat(paste('GithubRef:', branch, '\n'), file=fn, append=TRUE)
##cat(paste('GithubSHA1:', sha, '\n'), file=fn, append=TRUE)
pd <- utils::packageDescription('roxygen2')
v <- paste0(pd$Package, "_v", pd$Version)
v <- regmatches(v, regexpr("v(.*)$", v))
v <- substr(v, 2, nchar(v))

cat("RoxygenNote: ", v, "\n", file=fn, append=TRUE, sep = "")

##save('sha', file='spict/data/sha.rda')
