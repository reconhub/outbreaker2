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
#' @param mu a mutation rate
#'
prior.mu <- function(mu){
    return(dexp(mu, 1000, log=TRUE))
} # end prior.mu


#' @rdname priors
#' @export
#'
prior.all <- function(mu){
    return(prior.mu(mu))
} # end prior.all

