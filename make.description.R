fn <- 'spict/DESCRIPTION'
date <- format(Sys.Date(), "%Y-%m-%d")
cat('Package: spict
Type: Package
Title: Stochastic Surplus Production Model in Continuous-Time (SPiCT)
Version: 1.3
Date: ', date, '
Maintainer: Alexandros Kokkalis <alko@aqua.dtu.dk>
Authors@R: c(person(c("Martin", "W."), "Pedersen", email = "wpsgodd@gmail.com", role = "aut"),
             person(c("Casper", "W."), "Berg", email = "cbe@aqua.dtu.dk", role = "aut"),
             person("Tobias", "Mildenberger", email = "tobm@aqua.dtu.dk", role = "ctb"),
             person("Alexandros", "Kokkalis", email = "alko@aqua.dtu.dk", role = c("ctb", "cre")))
Author: Martin W. Pedersen [aut], 
    Casper W. Berg [aut],
    Tobias Mildenberger [ctb],
    Alexandros Kokkalis [ctb, cre]
Description: Fits a surplus production model to fisheries catch and biomass index data.
License: GPL (>=3)
Copyright: Martin Waever Pedersen
Depends:
    R (>= 3.0),
    TMB
Imports:
    graphics,
    methods,
    stats,
    utils
LinkingTo: TMB, RcppEigen
Suggests:
    parallel,
    mgcv,
    rjags,
    coda,
    knitr,
    rmarkdown,
    DLMtool
LazyData: true
VignetteBuilder: knitr
RoxygenNote: 6.0.1
BugReports: https://github.com/mawp/spict/issues\n', 
    file=fn)

# # Add Git hub stuff
# sha <- system('git rev-parse HEAD', intern=TRUE)
# branch <- system('git rev-parse --abbrev-ref HEAD', intern=TRUE)
# 
# cat(paste('GithubRepo: spict\n'), file=fn, append=TRUE)
# cat(paste('GithubRef:', branch, '\n'), file=fn, append=TRUE)
# cat(paste('GithubSHA1:', sha, '\n'), file=fn, append=TRUE)


#save('sha', file='spict/data/sha.rda')