#' Stores MCMC samples for outbreaker
#'
#' This function creates stores MCMC samples for outbreaker: augmented data and parameter states, likelihood, priors and posterior.
#'
#' @author Thibaut Jombart (\email{t_jombart@@imperial_ac_uk})
#'
#' @inheritParams outbreaker
#'
#' @param param a list of output items as returned by \code{outbreaker_create_mcmc}
#'
#' @param densities a list containing lists of functions computing densities, named: 'loglike' (log-likelihoods), 'priors' and 'posteriors'
#'
#' @param step an integer indicating the MCMC iteration being stored
#'
#' @export
#'
outbreaker_mcmc_store <- function(param_current, param_store, densities, step) {
    ## UPDATE COUNTER
    counter <- param_store$counter <- param_store$counter + 1

    ## STORE STEP
    param_store$step[counter] <- step

    ## STORE LIKELIHOOD, PRIOR, POSTERIOR
    param_store$like[counter] <- densities$loglike$all(param_current)
    param_store$prior[counter] <- densities$priors$all(param_current)
    param_store$post[counter] <- param_store$like[counter] + param_store$prior[counter]

    ## PARAMETERS AND AUGMENTED DATA
    param_store$mu[counter] <- param_current$mu
    param_store$pi[counter] <- param_current$pi
    param_store$alpha[[counter]] <- param_current$alpha
    param_store$t_inf[[counter]] <- param_current$t_inf
    param_store$kappa[[counter]] <- param_current$kappa

    return(param_store)
}

