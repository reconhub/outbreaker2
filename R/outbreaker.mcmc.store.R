#' Stores MCMC samples for outbreaker
#'
#' This function creates stores MCMC samples for outbreaker: augmented data and parameter states, likelihood, priors and posterior.
#'
#' @author Thibaut Jombart (\email{t.jombart@@imperial.ac.uk})
#'
#' @param param a list of output items as returned by \code{outbreaker.create.mcmc}
#' @param data a list of settings as returned by \code{outbreaker.data}
#' @param step an integer indicating the MCMC iteration being stored
#'
#' @export
#'
outbreaker.mcmc.store <- function(param, densities, step){
    ## UPDATE COUNTER
    counter <- param$counter <- param$counter + 1

    ## STORE STEP
    param$step[counter] <- step

    ## STORE LIKELIHOOD, PRIOR, POSTERIOR
    param$like[counter] <- densities$loglike$all(param)
    param$prior[counter] <- densities$priors$all(param)
    param$post[counter] <- param$like[counter] + param$prior[counter]

    ## PARAMETERS AND AUGMENTED DATA
    param$mu[counter] <- param$current.mu
    param$alpha[[counter]] <- param$current.alpha
    param$t.inf[[counter]] <- param$current.t.inf
    param$kappa[[counter]] <- param$current.kappa
    param$pi[counter] <- param$current.pi

    ## RETURN ##
    return(param)
} # end store.mcmc

