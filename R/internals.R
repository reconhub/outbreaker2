#' Internal functions for outbreaker2
#'
#' These functions represent the internal machinary behind outbreaker2: priors, likelihood functions, movements of parameters and augmented data.
#'
#' @author Thibaut Jombart \email{t.jombart@@imperial.ac.uk}
#'
#' @rdname internals
#'
#' @param times a vector of integers indicating dates of infections of the cases
#' @param ances a vector of indices of ancestors
#' @param w a vector of log probabilities of time intervals, starting at p(T=0)
ll.timing.infections <- function(times, ances, log.w){
    ## find indices in w (delay + 1 day)
    T <- (times-times[ances])+1

    ## avoid over-shooting
    T[T>length(log.w)] <- 0

    ## return
    return(sum(log.w[T], na.rm=TRUE))
} # end ll.timing.infections



ll.genet <- function(){}

