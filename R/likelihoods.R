

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
    if (data$N > 1) {
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
    if (nrow(data$D) > 1) {
        function(param, i=NULL) cpp.ll.genetic(data, param, i)
    } else {
        function(...) 0
    }
}





## This likelihood corresponds to the probability of a given number of unreported cases on an ancestry.

make.ll.reporting <- function(data) {
    if (data$N > 1) {
        function(param, i=NULL) cpp.ll.reporting(data, param, i)
    } else {
        function(...) 0
    }
}





## This likelihood corresponds to the sums of the separate timing likelihoods, which include:

##   - p(infection dates): see function cpp_ll_timing_infections
##   - p(collection dates): see function cpp_ll_timing_sampling

make.ll.timing <- function(data) {
    if (data$N > 1) {
        function(param, i=NULL) cpp.ll.timing(data, param, i)
    } else {
        function(...) 0
    }
}





##  This likelihood corresponds to the sums of the separate likelihoods, which include:

##   - p(infection dates): see function cpp_ll_timing_infections
##   - p(collection dates): see function cpp_ll_timing_sampling
##   - p(genetic diversity): see function cpp_ll_genetic
##   - p(missing cases): see function cpp_ll_reporting

make.ll.all <- function(data) {
    if (data$N > 1) {
        function(param, i=NULL) cpp.ll.all(data, param, i)
    } else {
        function(...) 0
    }
}





## This likelihood corresponds to the probability of observing contact between two individuals
## for a given ancestry

make.ll.contact <- function(data) {

    ## i will be the index of cases to be used, but it is useful to define it by default as all
    ## cases
    cases <- seq_len(data$N)

    if (nrow(data$contact)>0){
        function(param, i=cases) {

            ## discard cases with no ancestors as these cannot be informed by contact tracing data
            i <- i[!is.na(param$current.alpha[i])]

            ## Look up if contact is observed from the contact data matrix (data$contact)
            cij <- data$contact[cbind(i, param$current.alpha[i], deparse.level=0)]

            sum(log(param$current.eps^cij + param$current.eps*(cij-1)))

        }
    } else {
        function(...) 0
    }
}
