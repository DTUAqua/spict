spict
=====

An R-package for fittng surplus production models in continuous-time to fisheries catch data and biomass indices (either scientific or commercial). Main advantages of spict are:

1. All estimated reference points (MSY, Fmsy, Bmsy) are reported with uncertainties.

2. The model can be used for short-term forecasting and management strategy evaluation.

3. The model is fully stochastic in that observation error is included in catch and index observations, and process error is included in fishing and stock dynamics.

4. The model is formulated in continuous-time and can therefore incorporate arbitrarily sampled data.

## Package requirements

The package requires [`TMB`](http://www.tmb-project.org) to be installed. TMB is now a part of CRAN and can therefore be installed using the install.packages() command. For more information about TMB click [`here`](https://github.com/kaskr/adcomp).

## Installing the spict package

To install spict from GitHub use

```
library(devtools)
install_github("mawp/spict/spict")            # master branch
install_github("mawp/spict/spict", ref="dev") # development branch
```

Windows
-------
The above procedure using install_github() should now work on Windows (make sure to remove spict before trying to reinstall). If it doesn't work the old, but tedious, procedure can be used:

1. Start 64 bit R and change working directory to the (cloned or unzipped) ```spict``` folder.

2. From R run: ```source("install_windows.R")```

This requires that Rtools is installed. Rtools can be obtained [`here`](https://cran.r-project.org/bin/windows/Rtools/). When running install_windows.R remember to set your working directory to the spict directory containing install_windows.R.