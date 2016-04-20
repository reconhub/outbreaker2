#' Priors of outbreaker2
#'
#' These functions compute various priors used in outbreaker.
#'
#' @author Thibaut Jombart \email{t.jombart@@imperial.ac.uk}
#'
#' @rdname priors
#'
#' @export
#'
#' @param param a list containing parameters as returned by \code{outbreaker.mcmc.init}
#' @param config a set of settings as returned by \code{\link{outbreaker.config}}
#'
#'
prior.mu <- function(param, config){
    return(stats::dexp(param$current.mu, 1000, log=TRUE))
} # end prior.mu



#' @rdname priors
#' @export
#'
prior.pi <- function(param, config){
    return(stats::dbeta(param$current.pi, 10, 1, log=TRUE))
}



#' @rdname priors
#' @export
#'
prior.all <- function(param, config){
    return(prior.mu(param=param) + prior.pi(param=param))
} # end prior.all

