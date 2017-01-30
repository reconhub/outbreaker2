---
title: "outbreaker2: Rcpp API"
author: "Thibaut Jombart"
date: "2017-01-30"
output:
   rmarkdown::html_vignette:
     toc: true
     toc_depth: 2
vignette: >
  %\VignetteIndexEntry{Rcpp API}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---



# List of available functions

Here are all the C++ functions bound to R via Rcpp available in *outbreaker2*:


```r
env <- asNamespace("outbreaker2")
list_cpp_functions <- sort(ls(envir = env, pattern = "^cpp"))
list_cpp_functions
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



# Function calls

These functions take the following arguments:

```r
get_args <- function(f) args(get(f, envir = env))
get_formals <- function(f) formals(get(f, envir = env))

list_args <- lapply(list_cpp_functions, get_args)
names(list_args) <- list_cpp_functions
list_args
```

```
## $cpp_are_possible_ancestors
## function (t_inf, i) 
## NULL
## 
## $cpp_find_descendents
## function (alpha, i) 
## NULL
## 
## $cpp_find_local_cases
## function (alpha, i) 
## NULL
## 
## $cpp_ll_all
## function (data, param, i = NULL, custom_functions = NULL) 
## NULL
## 
## $cpp_ll_genetic
## function (data, param, i = NULL, custom_function = NULL) 
## NULL
## 
## $cpp_ll_reporting
## function (data, param, i = NULL, custom_function = NULL) 
## NULL
## 
## $cpp_ll_timing
## function (data, param, i = NULL, custom_functions = NULL) 
## NULL
## 
## $cpp_ll_timing_infections
## function (data, param, i = NULL, custom_function = NULL) 
## NULL
## 
## $cpp_ll_timing_sampling
## function (data, param, i = NULL, custom_function = NULL) 
## NULL
## 
## $cpp_move_alpha
## function (param, data, list_custom_ll = NULL) 
## NULL
## 
## $cpp_move_kappa
## function (param, data, config, list_custom_ll = NULL) 
## NULL
## 
## $cpp_move_mu
## function (param, data, config, custom_ll = NULL, custom_prior = NULL) 
## NULL
## 
## $cpp_move_pi
## function (param, data, config, custom_ll = NULL, custom_prior = NULL) 
## NULL
## 
## $cpp_move_swap_cases
## function (param, data, list_custom_ll = NULL) 
## NULL
## 
## $cpp_move_t_inf
## function (param, data, list_custom_ll = NULL) 
## NULL
## 
## $cpp_pick_possible_ancestor
## function (t_inf, i) 
## NULL
## 
## $cpp_prior_all
## function (param, config, custom_functions = NULL) 
## NULL
## 
## $cpp_prior_mu
## function (param, config, custom_function = NULL) 
## NULL
## 
## $cpp_prior_pi
## function (param, config, custom_function = NULL) 
## NULL
## 
## $cpp_sample1
## function (x) 
## NULL
## 
## $cpp_swap_cases
## function (param, i) 
## NULL
```




# Arguments

These arguments are:

```r
list_formals <- lapply(list_cpp_functions, get_formals)
names(list_formals) <- list_cpp_functions
args <- sort(unique(unlist(lapply(list_formals, names))))
```

- **`alpha`**: a vector of integers of length 'N' (number of cases), indicating infectors of each case, with values from 1 to N; missing values should be `NA`

- **`config`**: a list containing configuration settings as returned by `create_config`

- **`custom_function`**: a R function for a custom prior, with a single argument, which must be a list of parameters and augmented data with the class `outbreaker_param`; returned values must be **on the log scale**

- **`custom_functions`**: a list of R functions obeying the rules of `custom_function`, named according to the priors; currently available names are:

```
## [1] "mu" "pi"
```

- **`custom_ll`**: a R function for a custom likelihood, taking two arguments: `data` (see `data`), and `param` (see `param`)

- **`custom_prior`**: same as `custom_function`

- **`data`**: a valid 'outbreaker_data' list

- **`i`**: an integer scalar indicating the index of a case, from 1 to N (number of cases) 

- **`list_custom_ll`**: a list of R functions obeying the rules of `custom_ll`, named according to the computed likelihood component; available names are: 

```
## [1] "genetic"           "reporting"         "timing"           
## [4] "timing_infections" "timing_sampling"
```

- **`param`**: a list containing parameters and augmented data with the class `outbreaker_param`

- **`t_inf`**: a vector of integers of length N (number of cases), indicating infection dates of each case; missing values should be `NA`


- **`x`**: a vector of integers to be sampled from
