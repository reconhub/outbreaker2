
[![Travis-CI Build Status](https://travis-ci.org/reconhub/outbreaker2.svg?branch=master)](https://travis-ci.org/reconhub/outbreaker2)
[![Appveyor build status](https://ci.appveyor.com/api/projects/status/yj449x0yqhphvcrt/branch/master?svg=true)](https://ci.appveyor.com/project/thibautjombart/outbreaker2/branch/master)
[![Coverage Status](https://codecov.io/github/reconhub/outbreaker2/coverage.svg?branch=master)](https://codecov.io/github/reconhub/outbreaker2?branch=master)
[![CRAN Downloads](https://cranlogs.r-pkg.org/badges/outbreaker2)](https://cran.r-project.org/package=outbreaker2)
[![Downloads from Rstudio mirror](https://cranlogs.r-pkg.org/badges/grand-total/outbreaker2)](https://www.r-pkg.org:443/pkg/outbreaker2)
[![CRAN_Status_Badge](https://www.r-pkg.org/badges/version/outbreaker2)](https://cran.r-project.org/package=outbreaker2)

*outbreaker2: a framework for reconstructing disease outbreaks*
---------------------------------------------------------------

Welcome to the project page of *outbreaker2*, a Bayesian framework
 for integrating epidemiological and genetic data to reconstruct transmission
 trees of densely sampled outbreaks. It re-implements, generalises and replaces
 the model of [*outbreaker*](https://github.com/thibautjombart/outbreaker), and uses
 a modular approach which enables fine customisation of priors, likelihoods
 and parameter movements (see [customisation
 vignette](http://www.repidemicsconsortium.org/outbreaker2/articles/customisation.html)).

### NOTE: Correction to genetic likelihood

The genetic likelihood of the [original *outbreaker* paper](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1003457) was found to contain a minor mistake in accounting for unobserved generations of infection. As of June 7th 2019, *outbreaker2* will use the correct genetic likelihood published [here](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1006930).

The original genetic likelihood was:

$$\mu^{d(s_i,s_{\alpha_i})}(1 - \mu)^{(\kappa_i\times l(s_i, s_{\alpha_i})) - d(s_i,s_{\alpha_i})}$$

The corrected genetic likelihood is:

$$(\kappa_i\mu)^{d(s_i,s_{\alpha_i})}(1-\kappa_i\mu)^{l(s_i,s_{\alpha_i})-d(s_i,s_{\alpha_i})}$$

<br>

Installation
-------------

To install the stable version from CRAN:

```r
install.packages("outbreaker2")
```

To install the development version from github (requires Rtools on windows and
GSL headers on all platforms):


```r
devtools::install_github("reconhub/outbreaker2")
```

To add local copies of the vignettes, you will need to specify:

```r
devtools::install_github("reconhub/outbreaker2", build_vignettes = TRUE)
```

Then, to load the package, use:


```r
library("outbreaker2")
```



<br>

Documentation
-------------

*outbreaker2* is fully documented on a [dedicated
 website](http://www.repidemicsconsortium.org/outbreaker2/).

It also comes with the following vignettes:

- **`introduction`**: general introduction using a worked example.
- **`overview`**: brief overview of the package's content.
- **`customisation`**: customisation of priors, likelihoods, and movement functions.
- **`Rcpp_API`**: documentation for the Rcpp API.



<br>

Contributors
------------
- [Thibaut Jombart](https://github.com/thibautjombart)
- [Finlay Campbell](https://github.com/finlaycampbell)
- [Rich Fitzjohn](https://github.com/richfitz)
- [Gerry Tonkin-Hill](https://github.com/gtonkinhill)
- [Alexis Robert](https://github.com/alxsrobert)
- [Kristjan Eldjarn](https://github.com/kreldjarn)


See details of contributions
[here](https://github.com/reconhub/outbreaker2/graphs/contributors).

Contributions are welcome via **pull requests**.

Please note that this project is released with a [Contributor Code of
Conduct](https://github.com/reconhub/outbreaker2/blob/master/CONDUCT.md). By
participating in this project you agree to abide by its terms.

**Maintainer:** Finlay Campbell (finlaycampbell93@gmail.com)
