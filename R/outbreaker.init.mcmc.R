#' Initianilizes outputs for outbreaker
#'
#' This function creates initial outputs and parameter states for outbreaker.
#'
#' @author Thibaut Jombart (\email{t.jombart@@imperial.ac.uk})
#'
#' @param param a list of data items as returned by \code{outbreaker.create.mcmc}
#'
#' @param loglike a list of loglikelihood functions with enclosed data as returned by \code{create.loglike}
#'
#' @param priors a list of prior functions with enclosed parameters as returned by \code{create.priors}
#'
#' @export
#'
outbreaker.init.mcmc <- function(param.current, param.store, loglike, priors) {

    ## COMPUTE INITIAL LIKE/PRIOR/POST ##
    param.store$like[1] <- loglike$all(param.current)
    param.store$prior[1] <- priors$all(param.current)
    param.store$post[1] <- param.store$like[1] + param.store$prior[1]

    return(param.store)
}

