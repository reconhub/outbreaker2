#' Stores MCMC samples for outbreaker
#'
#' This function creates stores MCMC samples for outbreaker: augmented data and parameter states, likelihood, priors and posterior.
#'
#' @author Thibaut Jombart (\email{t.jombart@@imperial.ac.uk})
#'
#' @param chain a list of output items as returned by \code{outbreaker.mcmc.init}
#'
#' @param data a list of settings as returned by \code{outbreaker.data}
#'
#' @param config a list of output items as returned by \code{outbreaker.config}
#'
#' @export
#'
outbreaker.mcmc.store <- function(chain, data, config) {
    ## UPDATE COUNTER
    chain$counter <- chain$counter + 1

    ## STORE LIKELIHOOD, PRIOR, POSTERIOR
    counter <- chain$counter
    chain$like[counter] <- ll.all(data=data, config=config)
    chain$prior[counter] <- prior.all(chains)
    chain$post[counter] <- chain$like[counter] + chain$prior[counter]

    ## PARAMETERS AND AUGMENTED DATA
    chain$mu[counter] <- chain$current.mu
    chain$ances[[counter]] <- chain$current.ances
    chain$t.inf[[counter]] <- chain$current.t.inf

    ## RETURN ##
    return(chain)
} # end store.mcmc

