#' Create likelihood functions for outbreaker
#'
#' This function creates a named list of functions, each of which are compute log-likelihood components of outbreaker.
#' If these functions are provided as input, they will be used. Otherwise, default functions will be used.
#'
#' User-provided functions will be checked to ensure the right arguments are present. All log-likelihood functions should have to following arguments:
#' \describe{
#' \item{data}{a list of named items containing input data as returned by \code{\link{outbreaker.data}}}
#' \item{param}{a list of parameters as returned by \code{outbreaker.mcmc.init}}
#' }
#'
#' Note that unlike \code{outbreaker.create.moves}, likelihood functions have no mandatory names. They only need to be compatible with the movement functions used in \code{outbreaker.create.moves}. See 'details' for defaults.
#'
#' @author Thibaut Jombart \email{t.jombart@@imperial.ac.uk}
#'
#' @param ... a named list (see details) of functions, each computing a log-likelihood component.
#' @param loglike a list of functions as returned by \code{outbreaker.create.loglike}
#'
#' @details
#' The likelihood functions used by default are:
#' \describe{
#' \item{timing.infections}{corresponds to the probability of delays between cases (generation time)}
#' \item{timing.sampling}{corresponds to the probability of delays between infection and reporting}
#' \item{timing}{corresponds to the probability of all delays}
#' \item{genetic}{corresponds to the probability of the genetic sequences given the transmission tree}
#' \item{all}{corresponds to the overal probability of the data}
#' }
#'
#'
#' @return a list of named functions
#'
outbreaker.create.loglike <- function(..., loglike=NULL){
    ## PROCESS ... ONLY IF NO MOVES IS PASSED
    if(is.null(loglike)){
        loglike <- list(...)
    }


    ## SET DEFAULTS ##
    defaults <- list(timing.infections = ll.timing.infections,
                     timing.sampling = ll.timing.sampling,
                     timing = ll.timing,
                     genetic = ll.genetic,
                     all = ll.all
                     )

    ## MODIFY DEFAULTS WITH ARGUMENTS ##
    loglike <- modify.defaults(defaults, loglike, strict=FALSE)


    ## CHECK FUNCTIONS ##
    check.function.args <- function(f){
        args <- names(formals(f))
        if(identical(sort(args), c("data", "param")) ||
           identical(sort(args), c("data", "i", "param"))) {
            return(TRUE)
        } else {
            return(FALSE)
        }
    }

    args.checks <- sapply(unlist(loglike), check.function.args)
    if(!all(args.checks)){
        culprits <- names(args.checks)[!args.checks]
        culprits <- paste(culprits, collapse=",")
        stop("problems in movements of: ", culprits, "\narguments shoud be: 'data', 'param' (and optionally 'i' for local likelihoods)")
    }

    ## RETURN ##
    return(loglike)

} # end outbreaker.create.loglike
