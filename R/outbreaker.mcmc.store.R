#' Stores MCMC samples for outbreaker
#'
#' This function creates stores MCMC samples for outbreaker: augmented data and parameter states, likelihood, priors and posterior.
#'
#' @author Thibaut Jombart (\email{t.jombart@@imperial.ac.uk})
#'
#' @inheritParams outbreaker
#'
#' @param param a list of output items as returned by \code{outbreaker.create.mcmc}
#'
#' @param densities a list containing lists of functions computing densities, named: 'loglike' (log-likelihoods), 'priors' and 'posteriors'
#'
#' @param step an integer indicating the MCMC iteration being stored
#'
#' @export
#'
outbreaker.mcmc.store <- function(param.current, param.store, densities, step) {
    ## UPDATE COUNTER
    counter <- param.current$counter

    ## STORE STEP
    param.store$step[counter] <- step

    ## STORE LIKELIHOOD, PRIOR, POSTERIOR
    param.store$like[counter] <- densities$loglike$all(param)
    param.store$prior[counter] <- densities$priors$all(param)
    param.store$post[counter] <- param.store$like[counter] + param.store$prior[counter]

    ## PARAMETERS AND AUGMENTED DATA
    param.store$mu[counter] <- param.current$mu
    param.store$alpha[[counter]] <- param.current$alpha
    param.store$t.inf[[counter]] <- param.current$t.inf
    param.store$kappa[[counter]] <- param.current$kappa
    param.store$pi[counter] <- param.current$pi

    return(param.store)
}

