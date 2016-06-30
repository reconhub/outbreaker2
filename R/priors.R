
## In all of these functions:
## --------------------------
##
## we return a function in which prior parameters are enclosed; see 'likelihood.R' for more detail
## about the motivation (but basically, it's faster). These functions are later called by
## 'create.priors()' which will attach prior values and make a list of functions with enclosed
## parameters.

## Arguments are:
##
## 'config' a list of named items containing input data as returned by \code{\link{outbreaker.config}}


## We use an exponential prior for the mutation rate; the prior rate, which does not change in the
## MCMC, is enclosed in the returned function.

make.prior.mu <- function(config) {
    rate <- config$prior.mu
    function(param) {
        stats::dexp(param$current.mu, rate, log=TRUE)
    }
} 




## We use a beta prior for the reporting probability (which contrains it to lie between 0 and 1);
## the 2 shape parameters, which do not change in the MCMC, are enclosed in the returned function.

make.prior.pi <- function(config) {
    shape1 <- config$prior.pi[1]
    shape2 <- config$prior.pi[2]
    function(param) {
        stats::dbeta(param$current.pi, shape1, shape2, log=TRUE)
    }
}


## The prior of eps is uniform and must be constrained between 0.5 and 1 in order to positively 
## weight transmission networks supported by contact tracing data 

make.prior.eps <- function(config) {
  min <- config$prior.eps[1]
  max <- config$prior.eps[2]
  function(param) {
    stats::dunif(param$current.eps, min, max, log=TRUE)
  }
}