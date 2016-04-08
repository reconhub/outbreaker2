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
#' @importFrom stats dexp dbeta
#'
prior.mu <- function(param){
    return(dexp(param$current.mu, 1000, log=TRUE))
} # end prior.mu



#' @rdname priors
#' @export
#'
prior.pi <- function(param){
    return(dbeta(param$current.pi, 10, 1, log=TRUE))
}



#' @rdname priors
#' @export
#'
prior.all <- function(param){
    return(prior.mu(param=param) + prior.pi(param=param))
} # end prior.all

