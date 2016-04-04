#' Posteriors of outbreaker2
#'
#' These functions compute various posteriors used in outbreaker.
#'
#' @author Thibaut Jombart \email{t.jombart@@imperial.ac.uk}
#'
#' @rdname posteriors
#'
#' @param data a list of named items containing input data as returned by \code{\link{outbreaker.data}}
#'
#' @param param a list containing parameters as returned by \code{outbreaker.mcmc.init}
#'
post.genetic <- function(data, param){
    return(ll.genetic(data=data, param=param) +
           prior.mu(param=param))
} # end post.genetic




#' @rdname posteriors
#' @export
#'
post.reporting <- function(data, param){
    return(ll.reporting(data=data, param=param) +
           prior.pi(param=param))
} # end post.reporting




#' @rdname posteriors
#' @export
#'
post.all <- function(data, param){
    return(ll.all(data=data, param=param) + prior.all(param))
}
