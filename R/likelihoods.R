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
    ## compute delays between infection dates of cases and of their ancestors
    T <- param$current.t.inf[i] - param$current.t.inf[param$current.alpha[i]]

    ## avoid over-shooting: delays outside the range of columns in pre-computed log-densities
    ## (data$log.w.dens) will give a likelihood of zero
    if(any(T<1 | T>ncol(data$log.w.dens), na.rm=TRUE)) return(-Inf)

    ## return
    sum(data$log.w.dens[cbind(param$current.kappa[i], T)], na.rm=TRUE)
}





#' @rdname likelihoods
#' @export
#'
ll.timing.sampling <- function(data, param, i=NULL){
    ## compute delays
    T <- data$dates[i] - param$current.t.inf[i]
    T <- T[!is.na(T)]

    ## avoid over-shooting
    if(any(T<1 | T>length(data$log.f.dens))) return(-Inf)

    ## return
    sum(data$log.f.dens[T], na.rm=TRUE)
}





#' @rdname likelihoods
#' @export
#'
ll.timing <- function(data, param, i=NULL){
    ll.timing.infections(data=data, param=param, i=i) +
           ll.timing.sampling(data=data, param=param, i=i)
}


## NOTE: if the data can be bound to the ll functions via closures, we
## can get about a 1 us (50%) speed up on the calls.  The same is not
## really possible for param of course.  Moving from $... to [["..."]]
## looks to save a small amount too.
##
## You _might_ get slightly better performance, especially as the
## number of things grows, by factoring as
##  log(mu / (1 - mu)) * sum(nmut) + length(nmut) * log(1 - mu) * L
## which limits to 2 operations rather than 2*n

#' @rdname likelihoods
#' @export
#'
ll.genetic <- function(data, param, i=NULL){
    ## discard cases with no ancestors to avoid subsetting data$D with 'NA'
    i <- i[!is.na(param$current.alpha[i])]

    ## likelihood is based on the number of mutations between a case and its ancestor;
    ## these are extracted from a pairwise genetic distance matrix (data$D)
    nmut <- data$D[cbind(i, param$current.alpha[i], deparse.level=0)]
    sum(log(param$current.mu)*nmut + log(1 - param$current.mu)*
               (data$L - nmut), na.rm=TRUE)
} # end ll.genetic





#' @rdname likelihoods
#' @export
#' @importFrom stats dgeom
ll.reporting <- function(data, param, i=NULL){
    sum(dgeom(param$current.kappa[i]-1, prob=param$current.pi, log=TRUE),na.rm=TRUE)
} # end ll.reporting





#' @rdname likelihoods
#' @export
#'
ll.all <- function(data, param, i=NULL){
    ## NOTE: this should probably be replaced by a function iterating over all likelihood functions
    ll.timing(data=data, param=param, i=i) +
        ll.genetic(data=data, param=param, i=i) +
            ll.reporting(data=data, param=param, i=i)

} # end ll.all

