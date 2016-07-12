

## In all of these functions:
## --------------------------
##
## we return a function in which 'data' is enclosed this is substantially faster than having to
## pass the extra 'data' argument all the time, which does not change anyway; only things which
## effectively change are the parameters and augmented data, all stored in 'param', and the
## indices of cases for which we compute the likelihood 'i'.

## Arguments are:
##
## 'data' a list of named items containing input data as returned by \code{\link{outbreaker.data}}
## 'param' a list containing parameters as returned by \code{outbreaker.create.mcmc}
## 'i' an optional vector of integers, indicating subset of cases included in the likelihood computation; if NULL (default), all cases are used.




## This likelihood corresponds to the probability of observing infection dates of cases given the
## infection dates of their ancestors.

make.ll.timing.infections <- function(data) {
    if (data$N>1) {
        function(param, i=NULL) cpp.ll.timing.infections(data, param, i)
    } else {
        function(...) 0
    }
}




## This likelihood corresponds to the probability of reporting dates of cases given their
## infection dates.

make.ll.timing.sampling <- function(data) {
    if (data$N>1) {
        function(param, i=NULL) cpp.ll.timing.sampling(data, param, i)
    } else {
        function(...) 0
    }
}





## This likelihood corresponds to the probability of observing a number of mutations between cases
## and their ancestors. See src/likelihoods.cpp for details of the Rcpp implmentation.

make.ll.genetic <- function(data) {
    if (nrow(data$D)>1) {
        function(param, i=NULL) cpp.ll.genetic(data, param, i)
    } else {
        function(...) 0
    }
}





## This likelihood corresponds to the probability of a given number of unreported cases on an ancestry.

make.ll.reporting <- function(data) {
    if (nrow(data$D)>1) {
        function(param, i=NULL) cpp.ll.reporting(data, param, i)
    } else {
        function(...) 0
    }
}
