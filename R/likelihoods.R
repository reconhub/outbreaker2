#' Internal functions for outbreaker2
#'
#' These functions compute various likelihoods used in outbreaker.
#'
#' @author Thibaut Jombart \email{t.jombart@@imperial.ac.uk}
#'
#' @name likelihoods
#' @export
#'
#' @param data a list of named items containing input data as returned by \code{\link{outbreaker.data}}
#' @param param a list containing parameters as returned by \code{outbreaker.mcmc.init}
#' @param i an optional vector of integers, indicating subset of cases included in the likelihood computation; if NULL (default), all cases are used.
ll.timing.infections <- function(data, param, i=NULL){
    ## check i
    i <- check.i(data=data, i=i)

    ## compute delays
    T <- param$current.t.inf[i] - param$current.t.inf[param$current.ances[i]]
    T <- T[!is.na(T)]

    ## avoid over-shooting
    if(any(T<1 | T>length(data$log.w.dens))) return(-Inf)

    ## return
    return(sum(data$log.w.dens[T], na.rm=TRUE))
} # end ll.timing.infections





#' @rdname likelihoods
#' @export
#'
ll.timing.sampling <- function(data, param, i=NULL){
    ## check i
    i <- check.i(data=data, i=i)

    ## compute delays
    T <- data$dates[i] - param$current.t.inf[i]
    T <- T[!is.na(T)]

    ## avoid over-shooting
    if(any(T<1 | T>length(data$log.f.dens))) return(-Inf)

    ## return
    return(sum(data$log.f.dens[T], na.rm=TRUE))
} # end ll.timing.sampling





#' @rdname likelihoods
#' @export
#'
ll.timing <- function(data, param, i=NULL){
    ## check i
    i <- check.i(data=data, i=i)

    ## compute likelihood
    return(ll.timing.infections(data=data, param=param, i=i) +
           ll.timing.sampling(data=data, param=param, i=i))
} # end ll.timing





#' @rdname likelihoods
#' @export
#'
ll.genetic <- function(data, param, i=NULL){
    ## check i
    i <- check.i(data=data, i=i)

    ## discard cases with no ancestors


    ## compute likelihood
    if(is.null(data$dna)) return(0)
    nmut <- diag(data$D[i, param$current.ances[i]])
    return(sum(log(param$current.mu)*nmut + log(1 - param$current.mu)*(data$L - nmut), na.rm=TRUE))
} # end ll.genetic





#' @rdname likelihoods
#' @export
#'
ll.all <- function(data, param, i=NULL){
    ## check i
    i <- check.i(data=data, i=i)

    ## compute likelihood
    return(ll.timing(data=data, param=param, i=i) +
           ll.genetic(data=data, param=param, i=i)
           )
} # end ll.all

