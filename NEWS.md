outbreaker2 1.1.0 (2018-11-15)
==================

### Bug fixes

* The genetic likelihood of the [original *outbreaker* paper](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1003457) was found to contain a minor mistake in accounting for unobserved generations of infection. As of February 19th 2019, the development version of *outbreaker2* hosted on GitHub will use the correct genetic likelihood. When the new likelihood passes through peer-review and is published, we will push these changes to the CRAN version. For a description of the changes and a derivation of the new likelihood, see the [*outbreaker2* README](https://github.com/reconhub/outbreaker2/blob/master/README.md).

### New features

* Added support for case-labelling in plot.outbreaker_chains
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
