#' Initianilizes outputs for outbreaker
#'
#' This function creates initial outputs and parameter states for outbreaker.
#'
#' @author Thibaut Jombart (\email{t_jombart@@imperial_ac_uk})
#'
#' @param param a list of data items as returned by \code{outbreaker_create_mcmc}
#'
#' @param loglike a list of loglikelihood functions with enclosed data as returned by \code{create_loglike}
#'
#' @param priors a list of prior functions with enclosed parameters as returned by \code{create_priors}
#'
#' @export
#'
outbreaker_init_mcmc <- function(param_current, param_store, loglike, priors) {

    ## COMPUTE INITIAL LIKE/PRIOR/POST ##
    param_store$like[1] <- loglike$all(param_current)
    param_store$prior[1] <- priors$all(param_current)
    param_store$post[1] <- param_store$like[1] + param_store$prior[1]

    return(param_store)
}

