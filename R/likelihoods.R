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
#'
#' @param param a list containing parameters as returned by \code{outbreaker.mcmc.init}
#'
ll.timing.infections <- function(data, param){
    ## compute delays
    T <- param$current.t.inf - param$current.t.inf[param$current.ances]
    T <- T[!is.na(T)]

    ## avoid over-shooting
    if(any(T<1 | T>length(data$log.w.dens))) return(-Inf)

    ## return
    return(sum(data$log.w.dens[T], na.rm=TRUE))
} # end ll.timing.infections





#' @rdname likelihoods
#' @export
#'
ll.timing.sampling <- function(data, param){
    ## compute delays
    T <- data$dates - param$current.t.inf
    T <- T[!is.na(T)]

    ## avoid over-shooting
    if(any(T<1 | T>length(data$log.f.dens))) return(-Inf)

    ## return
    return(sum(data$log.f.dens[T], na.rm=TRUE))
} # end ll.timing.sampling





#' @rdname likelihoods
#' @export
#'
ll.timing <- function(data, param){
    return(ll.timing.infections(data=data, param=param) +
           ll.timing.sampling(data=data, param=param))
} # end ll.timing





#' @rdname likelihoods
#' @export
#'
ll.genetic <- function(data, param){
    if(is.null(data$dna)) return(0)
    nmut <- diag(data$D[, param$current.ances])
    return(sum(log(param$current.mu)*nmut + log(1 - param$current.mu)*(data$L - nmut), na.rm=TRUE))
} # end ll.genetic





#' @rdname likelihoods
#' @export
#'
ll.all <- function(data, param){
    return(ll.timing(data=data, param=param) +
           ll.genetic(data=data, param=param)
           )
} # end ll.all






#######################
## LOCAL LIKELIHOODS ##
#######################

#' @rdname likelihoods
#' @export
#'
ll.timing.infections.i <- function(data, param, i){
    ## check i
    check.i(data, i)

    ## escape if 'i' is imported
    if(is.na(param$current.ances[i])) return(0)

    ## compute delays
    T <- param$current.t.inf[i] - param$current.t.inf[param$current.ances[i]]
    T <- T[!is.na(T)]

    ## avoid over-shooting
    if(any(T<1 | T>length(data$log.w.dens))) return(-Inf)

    ## return
    return(sum(data$log.w.dens[T], na.rm=TRUE))
} # end ll.timing.infections.i





#' @rdname likelihoods
#' @export
#'
ll.timing.sampling.i <- function(data, param, i){
    ## check i
    check.i(data, i)

    ## escape if 'i' is imported
    if(is.na(param$current.ances[i])) return(0)

    ## compute delays
    T <- data$dates[i] - param$current.t.inf[i]
    T <- T[!is.na(T)]

    ## avoid over-shooting
    if(any(T<1 | T>length(data$log.f.dens))) return(-Inf)

    ## return
    return(sum(data$log.f.dens[T], na.rm=TRUE))
} # end ll.timing.sampling.i





#' @rdname likelihoods
#' @export
#'
ll.timing.i <- function(data, param, i){
    ## check i
    check.i(data, i)

    ## escape if 'i' is imported
    if(is.na(param$current.ances[i])) return(0)

    ## compute log-likelihood
    return(ll.timing.infections.i(data=data, param=param, i=i) +
           ll.timing.sampling.i(data=data, param=param, i=i))
} # end ll.timing.i





#' @rdname likelihoods
#' @export
#'
ll.genetic.i <- function(data, param, i){
    ## check i
    check.i(data, i)

    ## escape if 'i' is imported
    if(is.na(param$current.ances[i])) return(0)

    ## compute log-likelihood
    if(is.null(data$dna)) return(0)
    nmut <- data$D[i, param$current.ances[i]]
    return(sum(log(param$current.mu)*nmut + log(1 - param$current.mu)*(data$L - nmut), na.rm=TRUE))
} # end ll.genetic.i





#' @rdname likelihoods
#' @export
#'
ll.all.i <- function(data, param, i){
    ## check i
    check.i(data, i)

    ## escape if 'i' is imported
    if(is.na(param$current.ances[i])) return(0)

    ## compute log-likelihood
    return(ll.timing.i(data=data, param=param, i=i) +
           ll.genetic.i(data=data, param=param, i=i)
           )
} # end ll.all

