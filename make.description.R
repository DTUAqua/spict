fn <- 'spict/DESCRIPTION'
date <- format(Sys.Date(), "%Y-%m-%d")
cat('Package: spict
Type: Package
Title: Stochastic suplus Production model in Continuous-Time (SPiCT)
Version: 1.2.2
Date:', date, '
Author: Martin Waever Pedersen
Maintainer: Martin Waever Pedersen <mawp@dtu.dk>
Description: Fits a surplus production model to fisheries catch and biomass index data.
License: GPL (>=3)
Copyright: Martin Waever Pedersen <mawp@dtu.dk>
Depends:
    R (>= 3.0),
    TMB
LinkingTo: TMB, RcppEigen
Suggests:
    ellipse,
    parallel,
    mgcv,
    rjags,
    coda,
    knitr,
    rmarkdown,
    DLMtool
LazyData: true
VignetteBuilder: knitr\n',
    file=fn)

# Add Git hub stuff
sha <- system('git rev-parse HEAD', intern=TRUE)
branch <- system('git rev-parse --abbrev-ref HEAD', intern=TRUE)

cat(paste('GithubRepo: spict\n'), file=fn, append=TRUE)
cat(paste('GithubRef:', branch, '\n'), file=fn, append=TRUE)
cat(paste('GithubSHA1:', sha, '\n'), file=fn, append=TRUE)


#save('sha', file='spict/data/sha.rda')
