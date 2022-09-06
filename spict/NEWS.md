spict v1.3.7 (2022-09-06)
============

New features

* Added functionality for performing and plotting hindcast analysis 
  introducing the functions `hindcast()` and `plotspict.hindcast`.

Bug fixes:

* fixed retro when last year's observations where missing

spict v1.3.6 (2022-06-06)
============

Bug fixes:

* `sim.spict` saves the dtc (and dte) (closes #159)
* Fixed optimisation of absolute catch management scenario
* Fixed plotting of simulated data

spict v1.3.5 (2021-09-09)
============

New features:

* Use TMB simulation functionality to simulate data from fitted spict object.
  Set argument `use.tmb = TRUE` in `sim.spict()` function to use TMB
  `SIMULATE{}`.
* All functions using `parallel::mclapply()` have the argument `mc.cores` to
  specify the number of cores to be used. By default `mc.cores` is equal to 1.

Bug fixes

* Default labeling of management scenarios was not using increasing numbers, but
  reused the label "cutsomScenario_1". This is corrected and increasing numbers
  are now used, i.e. "customScenario_2", etc.
* In the case there is an abundance index after the last catch interval inside
  the intermediate period and the argument `intermediatePeriodCatch` in
  `manage()` is used, the intermediate period was calculated incorrectly, using
  the argument `inp$indpred` rather than `inp$indCpred`. This also led to an
  error message in `plotspict.catch()` with these management scenarios.

Minor changes:

* More informative error messages and user-friendly functionality of
  `sim.spict()`
* Three new arguments in the inp list (created by `check.inp`):
  `sim.random.effects` that allows to turn the simulation of random effects on
  and off, and `sim.fit` that allows to define whether the estimated parameters from the
  last fit or initial values should be used for simulation.
* The application of the retro function in the vignette was set to `mc.cores=1`
  to circumvent multithread MKL and parallel problems.


spict v1.3.4 (2021-02-21)
============

New features:

* Option to set `Blim` in the hockey-stick HCR to anything else than 0, by
  specifying two values in `breakpointB`, e.g. `breakpointB = c(0.3,0.5)`.
* New argument to set time for the evaluation of the hockey-stick HCR:
  `evalBreakpointB`.
* New plotting function: `plotspict.hcr` that visualises management scenarios
* See changes for each spict version using `news(package = "spict")`

Bug fixes

* Corrected biomass predictions for scenarios that use catch fractiles

Minor changes:

* Updated documentation, vignettes and guidelines



spict v1.3.3
============

Minor changes:

* Extra arguments `sumspict.manage` control what is returned (uncertainty, absolute states)
* `sumspict.manage` shows informative messages when scenarios are not comparable due to assumption differences
* Management scenarios targeting absolute catch (`cabs`) or fishing mortality (`fabs`). Absolute catch scenarios are based on finding the fishing mortalirt (`ffac`) that leads to the target catch using a simple optimisation.
* Management scenarios stored in the input object
* Legend in most plots have transparent backgrounds

Bug fixes:

* Fixed previously ignored argument `include.Ebinf` in `sumspict.manage`

spict v1.3.2
============

Minor changes:

* `TMB` version 1.7.1 is required
* Package `ellipse` is now imported
* Removed `plyr` dependency

Bug fixes:

* Fixed `test.spict` to work with the management updates

spict v1.3.1
============

New features:

* New function: `plotspict.retro.fixed`. Plots parameter estimates (and confidence intervals) in each of the retrospective runs

Minor changes:

* `plotspict.retro` returns Mohn's rho if `add.mohn` is TRUE
* `plotspict.retro` was polished, added legend, better spacing
* `mohns_rho` now works for the predicted index `Ipred`
* Retro plotting and Mohn's rho remove not converged retro runs (closes #49)
* Reduced the size of retro object and added an argument to control it (closes #97)
* Added argument checks in retro, plotspict_retro and mohns_rho

Bug fixes:

* Fixed Mohn's rho calculation that was ignoring the base run
* Fixed bug in Mohn's rho calculation that included forecast period
* Fixed shorten.inp for catch and effort

spict v1.3.0 (2020-10-23)
============

New features:

* The forecast & management section in SPiCT handbook was updated and describes
  the new management functionality of the package.
* Two new variables for controlling management-related settings are introduced:
 `inp$maninterval` and `inp$maneval`. The first is comprised of two numeric
 values which define the start and end of the management period, i.e. the period
 to which fishing mortality or catch restrictions are applied to and for which
 the catch (or TAC) is predicted for. `inp$maneval` defines the management
 evaluation time, i.e. the time point at which model states should be predicted
 for (corresponds to the 'p' quantities in the model, e.g. `logBpBmsy`). Thus, these
 variables replace the deprecated variables 'timepredc', 'dtpredc', and
 'manstart'. The new variables have a higher priority than the older ones, but
 old variables can still be used.
* With the new management variables come new default values for the management
  period, length of the period, and evaluation time: The management period
  starts now by default at the start of the next full new year after the last
  observation. This fixes a bug concerning default values for specific cases
  with seasonal catches (see below). Independent of annual or seasonal catches,
  the default duration of the management period is now 1 year. The default
  management evaluation time is at the end of the management period, i.e. by
  default at the end of the next year after the last observation. Please note
  that these changes might lead to different values for some predicted quantities
  (e.g. 'Predictions' in the summary output of a SPiCT assessment), because they
  affect the 'p' states in the model, e.g. `logBp`.
* The whole management functionality of the package is updated, including:
  - New management scenario labelled 'ices' in the `manage()` function
  - New function `add.man.scenario()` to add customised scenarios to `rep$man`
  - New function `get.TAC()` to estimate TAC of a specified scenario
  - New function `man.tac()` to estimate TAC of all scenarios in `rep$man`
  - New function `man.timeline()` to visualise management times in form of a timeline
  - New function `sumspict.manage()` to print a summary of all management scenarios
  - New function `check.man.time()` to check and correct model times in respect to management times
  - New function `man.select()` to select certain scenarios from `rep$man`
* New function `plot2()` allowing to plot a selection of 4 plots: biomass
  relative to Bmsy, Fishing mortality relative to Fmsy, catches, and the Kobe
  plot.
* New argument 'verbose' in all management-related functions, `fit.spict()` and
  `check.inp()` which allows to turn off console printing. This may be relevant
  for automatic reports and parallel scripts. By default `verbose = TRUE` and
  informative text messages and warnings are printed to the console.
* The new argument 'inp$reportmode' allows to change the number of ADreported
  quantities. So far, the default option (`0`) which reports all quantities, a
  minimal selection of quantities relevant for management (`1`), and the option
  to only report the predicted catch (`logCp` with `inp$reportmode == 2`) are
  implemented.
* A new family of ADreported quantities is introduced: The 'm' quantities, e.g.
  `logFmFmsy`, correspond to the start of the management period.
* New function `mohns_rho` that calculates Mohn's rho from the retrospective analysis
  of user selected quantities

Minor changes:

* The guidelines vignette was updated.
* New test scripts for all management-related functions.
* The manbase scenario was removed.
* The `manage()` function overwrites all scenarios in `rep$man`.
* The quantities 'perc.dF' and 'perc.dB' in the management summary output relate
  the fishing mortality or biomass, respectively, from the end of the management
  period to the start of the management period of each scenario. Thus, these
  quantities are omitted from the summary output if the timelines of the
  scenarios vary.
* DLMtool dependent functions were removed.

Bug fixes:

* The argument `inp$indlastobs` which is used to estimate all 'l' quantities,
  e.g. `logFlFmsy`, did not consider the length of the catch (or effort)
  interval. Fixing this bug might lead to different values for some quantities
  (e.g. 'States' in the summary output of a SPiCT assessment), because they
  affect the 'l' states in the model, e.g. `logBl`.
* The new default values of management-related variables explained above fix
  inconsistent default values for seasonal data, which might have lead
  'manstart' to be after 'timepredc'.
* Fixed bug in retro function when effort data are used.
* Allow to plot into same plotting device for all plotspict.* family plots, e.g.
  for changing the ylim or adding a horizontal reference line.
* All management scenarios allow for the definition of intermediate periods
  (period between last observation and start of the management period). By
  default the assumptions during the intermediate period are equal (continuing
  the F process) for all scenarios.
* The argument 'maxtime' in the function `shorten.inp()` did not account for the
  length of the catch (or effort) interval. Thus, a catch observation for the
  interval 1990-1991 was included although `maxtime = 1990`. After this bug fix,
  this observation is not included anymore as it ranges beyond 'maxtime'.


spict v1.2.8 (2019-12-02)
============

New features:

* Added new vignette: "Guidelines for the use of the stochastic production model
  in continuous time (SPiCT)"
* Added helper function `calc.bmsyk` that returns the Bmsy / K ratio of a fitted
  object
* Added helper function `calc.om` that return a matrix with the order of
  magnitude of the confidence interval range of F/Fmsy and B/Bmsy
* Added helper function `shorten.inp` that cuts the input time series to the
  specified range
* Added option for gamma distributed prior for the parameter determining the
  shape of the production (`n`) curve on the log scale
* Added option to report the states in terms of biomass and fishing mortality
  relative to the overall average over the whole time series (including the
  prediction period) `Brel` and `Frel`, respectively

Minor changes:

* Added sdfac argument in `manage` that sets the standard deviation factor keep
  current catch scenario
* Increased maximum number of iterations and function evaluations of the
  `nlminb` optimiser
* Increased tested code coverage

Bug fixes:

* Use non seasonal F in `manage`
* Use better default value for the standard deviation factor in `take.c`
* Remove inverse gamma priors which were not fully functional


spict v1.2.7 (2019-04-03)
============

New features:

* Added new management scenario: ICES MSY advice rule
* Added NEWS file

Minor changes:

* Updated documentation for check.inp and spict2DLMtool
* Makefile now takes the version number from DESCRIPTION file
