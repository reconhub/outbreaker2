#' Internal functions for outbreaker2
#'
#' These functions compute various likelihoods used in outbreaker.
#'
#' @author Thibaut Jombart \email{t.jombart@@imperial.ac.uk}
#'
#' @rdname likelihoods
#'
#' @param times a vector of integers indicating dates of infections of the cases
#' @param ances a vector of indices of ancestors
#' @param log.w a vector of log probabilities of time intervals (between infections), starting at p(T=0)
#'
#' @export
ll.timing.infections <- function(times, log.w, ances){
    ## find indices in w (delay + 1 day)
    T <- (times-times[ances])+1

    ## avoid over-shooting
    T[T>length(log.w)] <- 0

    ## return
    return(sum(log.w[T], na.rm=TRUE))
} # end ll.timing.infections





#' @rdname likelihoods
#' @export
#'
#' @param sampling.time a vector of integers indicating dates of sampling/observation/reporting of the cases
#' @param log.f a vector of log probabilities of time intervals (from infection to collection), starting at p(T=0)
#'
ll.timing.sampling <- function(times, log.f, sampling.times){
    ## find indices in w (delay + 1 day)
    T <- (sampling.times-times)+1

    ## avoid over-shooting
    T[T>length(log.f)] <- 0

    ## return
    return(sum(log.f[T], na.rm=TRUE))
} # end ll.timing.sampling





#' @rdname likelihoods
#' @export
#'
ll.timing <- function(times, log.w, log.f, ances, sampling.times){
    return(ll.timing.infections(times=times, log.w=log.w, ances=ances) +
           ll.timing.sampling(times=times, log.f=log.f, sampling.times=sampling.times))
} # end ll.timing





#' @rdname likelihoods
#' @export
#'
#' @param D a matrix of pairwise genetic distances
#' @param mu a mutation rate
#' @param gen.length length of the genetic sequences
ll.genetic <- function(D, gen.length, ances, mu){
    nmut <- diag(D[,ances])
    return(sum(log(mu)*nmut + log(1-mu)*(gen.length-nmut), na.rm=TRUE))
} # end ll.genetic





#' @rdname likelihoods
#' @export
#'
ll.all <- function(times, D, gen.length, log.w, ances, mu){
    return(ll.timing.infections(times=times, ances=ances, log.w=log.w) +
           ll.genetic(D=D, ances=ances, mu=mu, gen.length=gen.length))
}
