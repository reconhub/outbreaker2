
## ## In all of these functions:
## ## --------------------------
## ##
## ## we return a function in which likelihood and prior functions are enclosed; see 'likelihood_R' for
## ## more detail about the motivation (but basically, it's faster). These functions are later called
## ## by 'create_posterior()' which will create a list of functions with enclosed
## ## loglikelihood and prior functions.
## ##




## ## The 'genetic posterior' includes the genetic likelihood and the priors for mutation rate 'mu'.

## make_post_genetic <- function(loglike, priors) {
##     ll_genetic <- loglike$genetic
##     prior_mu <- priors$mu
##     function(param) {
##         ll_genetic(param) + prior_mu(param)
##     }
## }



## ## The 'reporting posterior' includes the reporting likelihood and the priors for the reporting probability 'pi'.

## make_post_reporting <- function(loglike, priors) {
##     ll_reporting <- loglike$reporting
##     prior_pi <- priors$pi
##     function(param) {
##         ll_reporting(param) + prior_pi(param)
##     }
## }



## ## The global posterior includes the all likelihoods and priors

## make_post_all <- function(loglike, priors) {
##     ll_all <- loglike$all
##     prior_all <- priors$all
##     function(param) {
##         ll_all(param) + prior_all(param)
##     }
## }


