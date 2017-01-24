
## ON THE PURPOSE OF THESE FUNCTIONS

## These functions are no longer used in outbreaker2, but were part of the
## original implementation, and are still used in testing procedures to ensure
## that Rcpp versions give identical results.






## We use an exponential prior for the mutation rate; the prior rate, which does
## not change in the MCMC, is enclosed in the returned function.

.prior_mu <- function(param, rate) {
    stats::dexp(param$mu, rate, log = TRUE)
}






## We use a beta prior for the reporting probability (which contrains it to lie
## between 0 and 1); the 2 shape parameters, which do not change in the MCMC,
## are enclosed in the returned function.

.prior_pi <- function(param, shape1, shape2) {
    stats::dbeta(param$pi, shape1, shape2, log = TRUE)
}

