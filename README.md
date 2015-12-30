spict
=====

An R-package for fittng surplus production models in continuous-time to fisheries catch data and biomass indices (either scientific or commercial). Main advantages of spict are:

1. All estimated reference points (MSY, Fmsy, Bmsy) are reported with uncertainties.

2. The model can be used for short-term forecasting and management strategy evaluation.

3. The model is fully stochastic in that observation error is included in catch and index observations, and process error is included in fishing and stock dynamics.

4. The model is formulated in continuous-time and can therefore incorporate arbitrarily sampled data.

## Package requirements

The package requires [`TMB`](http://www.tmb-project.org) to be installed. TMB is available on CRAN and can therefore be installed using the install.packages() command. However, it appears that on Windows spict is incompatible with version 1.6.5 of TMB. Instead it is recommended to use TMB version 1.6.2, which can be obtained from [`github`](https://github.com/kaskr/adcomp/tree/v1.6.2). On Windows: download ZIP, unzip and use the install_windows.R installer to install TMB.

For more information about TMB click [`here`](https://github.com/kaskr/adcomp).

## Linux

The package can be installed from github:

```
library(devtools)
install_github("mawp/spict/spict")
```

Windows
-------
Important: on Windows the spict package has been tested and works with version 1.6.2 of TMB (but not 1.6.5). Get version 1.6.2 as described above.

1. Download the compressed spict source code by clicking the Download ZIP button. Alternatively install git and clone the spict repository.

2. Start 64 bit R and change working directory to the (cloned or unzipped) ```spict``` folder.

3. From R run: ```source("install_windows.R")```

This requires that Rtools is installed, which may be the case if TMB is installed. If not Rtools can be obtained [`here`](https://cran.r-project.org/bin/windows/Rtools/). When running install_windows.R remember to set your working directory to the spict directory containing install_windows.R. 