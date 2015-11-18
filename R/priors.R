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
prior.mu <- function(mu)
{
    return(dexp(mu, 1000, log=TRUE))
} # end prior.mu


#' @rdname priors
#' @export
#' @param chain a list of output items as returned by \code{outbreaker.mcmc.init}
#'
prior.all <- function(chain)
{
    return(prior.mu(chain$current.mu))
} # end prior.all

