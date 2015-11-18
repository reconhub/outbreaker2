#' Posteriors of outbreaker2
#'
#' These functions compute various posteriors used in outbreaker.
#'
#' @author Thibaut Jombart \email{t.jombart@@imperial.ac.uk}
#'
#' @rdname posteriors
#'
#' @param D a matrix of pairwise genetic distances
#' @param gen.length length of the genetic sequences
#' @param ances a vector of indices of ancestors
#' @param mu a mutation rate
#'
post.genetic <- function(D, gen.length, ances, mu)
{
    return(ll.genetic(D=D, gen.length=gen.length, ances=ances, mu=mu) +
           prior.mu(mu=mu))
} # end post.genetic





#' @rdname posteriors
#' @export
#' @param data a list of data items as returned by \code{outbreaker.data}
#' @param chain a list of output items as returned by \code{outbreaker.mcmc.init}
#'
post.all <- function(data=data, chain=chain)
{
    return(ll.all(data=data, chain=chain) + prior.all(chain))
}
