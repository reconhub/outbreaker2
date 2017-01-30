---
title: "outbreaker2: package overview"
author: "Thibaut Jombart"
date: "2017-01-30"
output:
   rmarkdown::html_vignette:
     toc: true
     toc_depth: 2
vignette: >
  %\VignetteIndexEntry{Overview}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---



*oubtreaker2* provides ...



<br>

# Installing the package

To install the current stable, CRAN version of the package, type:

```r
install.packages("outbreaker2")
```

To benefit from the latest features and bug fixes, install the development, *github* version of the package using:

```r
devtools::install_github("reconhub/outbreaker2")
```

Note that this requires the package *devtools* installed.



<br>

# Main visible functions of the package

The main functions of the package include:

- **`outbreaker`**: the returned object is an instance of the (S3) class *outbreaker_chains*.

- **`plot`**: this method (see `?plot.outbreaker_chains` for details) plots *outbreaker_chains* objects.

- **`summary`**: this method (see `?summary.outbreaker_chains` for details) provides summaries for the various outputs of *outbreaker2*, stored in an *outbreaker_chains* object.

- **`outbreaker_data`**: function processing input data.

- **`create_config`**: function creating default settings, also used for specifying customised settings for *outbreaker2*.

- **`custom_priors`**: function used for specifying customised functions to be used as priors in *outbreaker2*.

- **`custom_likelihoods`**: function used for specifying customised functions to be used for likelihood computation in *outbreaker2*.

- **`custom_moves`**: function used for specifying customised functions to be used for moving parameters and augmented data in *outbreaker2*.





<br>

# Main internal functions

*outbreaker2* uses many functions internally which are not visible to the user when loading the package. However, some of these functions will be useful when designing custom likelihoods or movement functions, or when contributing code. The most useful ones are C++ functions bound to R using Rcpp. The list of these functions is:


```r
env <- asNamespace("outbreaker2")
ls(envir = env, pattern = "^cpp")  
```

```
##  [1] "cpp_are_possible_ancestors" "cpp_find_descendents"      
##  [3] "cpp_find_local_cases"       "cpp_ll_all"                
##  [5] "cpp_ll_genetic"             "cpp_ll_reporting"          
##  [7] "cpp_ll_timing"              "cpp_ll_timing_infections"  
##  [9] "cpp_ll_timing_sampling"     "cpp_move_alpha"            
## [11] "cpp_move_kappa"             "cpp_move_mu"               
## [13] "cpp_move_pi"                "cpp_move_swap_cases"       
## [15] "cpp_move_t_inf"             "cpp_pick_possible_ancestor"
## [17] "cpp_prior_all"              "cpp_prior_mu"              
## [19] "cpp_prior_pi"               "cpp_sample1"               
## [21] "cpp_swap_cases"
```

See the vignette on Rcpp API for a detail of these functions.




