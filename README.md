spict
=====

Surplus Production in Continuous-Time

## Installing the package

The package requires [`TMB`](http://www.tmb-project.org) to be installed.

To install the package from GitHub (not on Windows) use

```
library(devtools)
install_github("mawp/spict/spict")
```

Windows
-------
1. Start 64 bit R and change working directory to the (cloned or unzipped) ```spict``` folder.

2. From R run: ```source("install_windows.R")```

This requires that Rtools is installed, which is probably the case if TMB is installed. When running install_windows.R remember to set your working directory to the spict directory containing install_windows.R.