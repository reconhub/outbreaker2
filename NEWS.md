outbreaker2 1.1.4 (2022-12-16)
==================

### Bug fixes

* Added `negative_si` option to `create_config` to allow ancestries with negative serial intervals to be immediately discarded.
* Removed C++11 SystemRequirements specification (R 4.0+ defaults to C++11)



<br>
<br>

outbreaker2 1.1.3 (2022-05-18)
==================

### Bug fixes

* Email addresses fixed
* Adegenet added to Suggests



<br>
<br>

outbreaker2 1.1.2 (2021-01-25)
==================

### New features
* Custom likelihoods now support local likelihood calculations via the optional
  `i` argument; this will speed up likelihood calculations by avoiding global
  calculations where these are unnecessary

### Bug fixes

* Random index shuffling in cpp_move_swap_cases to ensure correct mixing
* Negative serial intervals now allowed in can_be_ances
* ID labels fix in plotting functions



<br>
<br>

outbreaker2 1.1.1 (2020-01-31)
==================

### Bug fixes

* Fixed uninitialized variable use
* Updated email addresses



<br>
<br>

outbreaker2 1.1.0 (2019-06-06)
==================

### Bug fixes

* The genetic likelihood of the [original *outbreaker* paper](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1003457) was found to contain a minor mistake in accounting for unobserved generations of infection. As of June 7th 2019, *outbreaker2* will use the correct genetic likelihood published [here](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1006930).

### New features

* *outbreaker2* now accepts `epicontacts` objects for contact data
* Contact tracing data can be directed or not by toggling `ctd_directed`
* Cyclical consensus trees can now be avoided by using `method = 'decyle'` in `summary.outbreaker_chains`
* Added support for case-labelling in `plot.outbreaker_chains`
* Added progress bar


<br>
<br>

outbreaker2 1.0.1 (2018-03-06)
==================

### Minor bug fixes

* Fixed convolution operator to keep values in log space with high kappa (thanks to [@gtonkinhill](https://github.com/gtonkinhill/))
* C++ garbage collection fix



<br>
<br>

outbreaker2 1.0.0 (2017-11-24)
==================

### First version of outbreaker2!

* This package re-implements and generalises the package
  [*outbreaker*](https://CRAN.R-project.org/package=outbreaker)

* It provides various functions to process input data, specify settings of the
  method, and even specify custom priors, likelihoods, or movement functions.

* Outputs of the main function `outbreaker` have a dedicated class extending the
  regular `data.frame`, with `print`, `summary` and `plot` methods.

* Its engine is written in C++, interfaced via Rcpp.

* Its C++ API is accessible via the function `get_cpp_api()`.

* Functionalities are documented in 4 vignettes.
