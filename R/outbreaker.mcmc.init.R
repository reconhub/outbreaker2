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
#' @export
#'
outbreaker.init.mcmc <- function(param, loglike){

    ## COMPUTE INITIAL LIKE/PRIOR/POST ##
    param$like[1] <- loglike$all(param)
    param$prior[1] <- prior.all(param=param)
    param$post[1] <- param$like[1] + param$prior[1]

    return(param)
}

