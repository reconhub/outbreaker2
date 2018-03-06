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
