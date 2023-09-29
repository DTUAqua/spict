spict
=====

An R-package for fittng surplus production models in continuous-time to fisheries catch data and biomass indices (either scientific or commercial). Main advantages of spict are:

1. All estimated reference points (MSY, Fmsy, Bmsy) are reported with uncertainties.

2. The model can be used for short-term forecasting and management strategy evaluation.

3. The model is fully stochastic in that observation error is included in catch and index observations, and process error is included in fishing and stock dynamics.

4. The model is formulated in continuous-time and can therefore incorporate arbitrarily sampled data.

## Help files

A vignette for the package is available [`here`](https://github.com/DTUAqua/spict/raw/master/spict/inst/doc/spict_handbook.pdf), and serves as an introduction to the package and its functionality. The vignette also contains description of the more advanced features of the package.

A document with technical guidelines for using SPiCT is available [here](https://github.com/DTUAqua/spict/raw/master/spict/inst/doc/spict_guidelines.pdf). This is a living document that has a list of things to check before accepting an assessment and some options to deal with more dificult data sets.

The package also contains reasonable documentation in the form of help texts associated with each function (some may not be fully up-to-date). These can be accessed in the usual R manner by writing e.g. ```?check.inp```. A good place to start (in addition to reading the vignette) is to read ```?check.inp``` and ```?fit.spict```.

## Citation

The underlying model used in the package is described in a published [`paper`](http://onlinelibrary.wiley.com/doi/10.1111/faf.12174/full). A preprint of the paper is included in the package in the [`inst`](https://github.com/DTUAqua/spict/tree/master/spict/inst) folder and can be downloaded [here](https://github.com/DTUAqua/spict/raw/master/spict/inst/spict.pdf). To get citation information write `citation(spict)` in the command line.

## Package requirements

The package requires [`TMB`](http://www.tmb-project.org) to be installed. TMB is now a part of CRAN and can therefore be installed using ```install.packages("TMB", type="source")```. For more information about TMB click [`here`](https://github.com/kaskr/adcomp).

## Installing the spict package

To install spict from GitHub use

```
library(remotes)
install_github("DTUAqua/spict/spict")            # master branch
```

Windows
-------

Installing spict in Windows requires Rtools, available [`here`](https://cran.r-project.org/bin/windows/Rtools/). When installing Rtools, it is important to check the option ""**Add Rtools to the system PATH**".

Then using the above procedure using```install_github()``` should work (make sure to remove spict before trying to reinstall).

If it doesn't work the old, but tedious, procedure can be used:

1. Start 64 bit R and change working directory to the (cloned or unzipped) ```spict``` folder.

2. From R run: ```source("install_windows.R")```

When running install_windows.R remember to set your working directory to the spict directory containing install_windows.R.


Funding
-------

The development of `spict` was partially funded by the European Maritime and Fisheries Fund (EMFF) and the Danish Fisheries Agency via the RoMA project (Grant agreement number 33113-B-20-183).
