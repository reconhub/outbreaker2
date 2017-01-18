
## In all of these functions:
## --------------------------
##
## we return a function in which prior parameters are enclosed; see
## 'likelihood_R' for more detail about the motivation (but basically, it's
## faster). These functions are later called by 'create_priors()' which will
## attach prior values and make a list of functions with enclosed parameters.

## Arguments are:
##
## 'config' a list of named items containing input data as returned by
## \code{\link{outbreaker_c onfig}}


## We use an exponential prior for the mutation rate; the prior rate, which does
## not change in the MCMC, is enclosed in the returned function.

.prior_mu <- function(param, rate) {
    stats::dexp(param$mu, rate, log = TRUE)
}


bind_prior_mu <- function(config) {
    rate <- config$prior_mu
    function(param) .prior_mu(param, rate)
} 




## We use a beta prior for the reporting probability (which contrains it to lie
## between 0 and 1); the 2 shape parameters, which do not change in the MCMC,
## are enclosed in the returned function.

.prior_pi <- function(param, shape1, shape2) {
    stats::dbeta(param$pi, shape1, shape2, log = TRUE)
}

bind_prior_pi <- function(config) {
    shape1 <- config$prior_pi[1]
    shape2 <- config$prior_pi[2]
    
    function(param) .prior_pi(param, shape1, shape2)
                              
}


