#' Movements of augmented data and parameters for outbreaker2
#'
#' These functions are used to move in the parameter space in the MCMC.
#'
#' @author Thibaut Jombart \email{t.jombart@@imperial.ac.uk}
#'
#' @rdname moves
#'
#' @param D a matrix of pairwise genetic distances
#' @param gen.length length of the genetic sequences
#' @param ances a vector of indices of ancestors
#' @param mu a mutation rate
#' @param t.inf a vector of integers indicating dates of infections of the cases
#' @param log.w a vector of log probabilities of time intervals (between infections), starting at p(T=0)
#' @param dev a deviation from the previous mu
#' @param lunif a logged random variate from Unif(0,1)
#' @importFrom stats rnorm
#'
move.mu <- function(D, gen.length, ances, mu, dev, lunif){ # assumes symmetric proposal
    new.mu <- mu + dev

    ## escape if new.mu<0
    if(new.mu<0) return(mu)

    ## compute log ratio
    logratio <- post.genetic(D=D, gen.length=gen.length, ances=ances, mu=new.mu) -
        post.genetic(D=D, gen.length=gen.length, ances=ances, mu=mu)

    ## accept/reject
    if(logratio >= lunif) return(new.mu)
    return(mu)
} # end move.mu



#' @rdname moves
#' @export
#'
move.t.inf <- function(t.inf, log.w, log.f, ances, sampling.times, lunif){ # assumes symmetric proposal
    ## propose new t.inf
    new.t.inf <- t.inf + sample(-1:1, size=length(t.inf), replace=TRUE, prob=c(.1,8,.1))

    ## compute log ratio
    logratio <- ll.timing(t.inf=new.t.inf, log.w=log.w, log.f=log.f, ances=ances, sampling.times=sampling.times) -
        ll.timing(t.inf=t.inf, log.w=log.w, log.f=log.f, ances=ances, sampling.times=sampling.times)

    ## accept/reject
    if(logratio >= lunif) return(new.t.inf)
    return(t.inf)
} # end move.mu




#' @rdname moves
#' @export
#'
move.ances <- function(t.inf, sampling.times, D, gen.length, log.w, log.f, ances, mu, lunif){

}
