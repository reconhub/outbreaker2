#' Stores MCMC samples for outbreaker
#'
#' This function creates stores MCMC samples for outbreaker: augmented data and parameter states, likelihood, priors and posterior.
#'
#' @author Thibaut Jombart (\email{t.jombart@@imperial.ac.uk})
#'
#' @param param a list of output items as returned by \code{outbreaker.mcmc.init}
#'
#' @param data a list of settings as returned by \code{outbreaker.data}
#'
#' @export
#'
outbreaker.mcmc.store <- function(param, data)
{
    ## UPDATE COUNTER
    param$counter <- param$counter + 1

    ## STORE LIKELIHOOD, PRIOR, POSTERIOR
    counter <- param$counter
    param$like[counter] <- ll.all(data=data, param=param)
    param$prior[counter] <- prior.all(param)
    param$post[counter] <- param$like[counter] + param$prior[counter]

    ## PARAMETERS AND AUGMENTED DATA
    param$mu[counter] <- param$current.mu
    param$ances[[counter]] <- param$current.ances
    param$t.inf[[counter]] <- param$current.t.inf

    ## RETURN ##
    return(param)
} # end store.mcmc

