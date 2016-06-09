
fn <- 'spict/DESCRIPTION'
cat('Package: spict
Type: Package
Title: Stochastic suplus Production model in Continuous-Time (SPiCT)
Version: 1.0
Date: 2016-06-09
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
    coda
LazyData: true\n',
    file=fn)

# Add Git hub stuff
sha <- system('git rev-parse HEAD', intern=TRUE)
branch <- system('git rev-parse --abbrev-ref HEAD', intern=TRUE)

cat(paste('GithubRepo: spict\n'), file=fn, append=TRUE)
cat(paste('GithubRef:', branch, '\n'), file=fn, append=TRUE)
cat(paste('GithubSHA1:', sha, '\n'), file=fn, append=TRUE)


#save('sha', file='spict/data/sha.rda')
