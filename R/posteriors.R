
## In all of these functions:
## --------------------------
##
## we return a function in which likelihood and prior functions are enclosed; see 'likelihood.R' for
## more detail about the motivation (but basically, it's faster). These functions are later called
## by 'create.posterior()' which will create a list of functions with enclosed
## loglikelihood and prior functions.
##




## The 'genetic posterior' includes the genetic likelihood and the priors for mutation rate 'mu'.

make.post.genetic <- function(loglike, priors){
    ll.genetic <- loglike$genetic
    prior.mu <- priors$mu
    function(param){
        ll.genetic(param) + prior.mu(param)
    }
}



## The 'reporting posterior' includes the reporting likelihood and the priors for the reporting probability 'pi'.

make.post.reporting <- function(loglike, priors){
    ll.reporting <- loglike$reporting
    prior.pi <- priors$pi
    function(param){
        ll.reporting(param) + prior.pi(param)
    }
}


