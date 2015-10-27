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
#' @param times a vector of integers indicating dates of infections of the cases
#' @param log.w a vector of log probabilities of time intervals (between infections), starting at p(T=0)
#'
post.genetic <- function(D, gen.length, ances, mu){
    return(ll.genetic(D=D, gen.length=gen.length, ances=ances, mu=mu) +
           prior.mu(mu))
} # end post.genetic





#' @rdname posteriors
#' @export
#'
post.all <- function(times, D, gen.length, log.w, ances, mu){
    return(ll.timing.infections(times=times, ances=ances, log.w=log.w) +
           ll.genetic(D=D, ances=ances, mu=mu, gen.length=gen.length) +
           prior.mu(mu))
}
