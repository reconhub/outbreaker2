#' Create posterior functions for outbreaker
#'
#' This function creates a named list of functions, each of which are compute log-posterior components of outbreaker.
#' If these functions are provided as input, they will be used. Otherwise, default functions will be used.
#'
#' User-provided functions will be checked to ensure the right arguments are present. All log-posterior functions should have to following arguments:
#' \describe{
#' \item{data}{a list of named items containing input data as returned by \code{\link{outbreaker.data}}}
#' \item{param}{a list of parameters as returned by \code{\link{outbreaker.mcmc.init}}} 
#' }
#'
#' Note that unlike \code{outbreaker.create.moves}, posterior functions have no mandatory names. They only need to be compatible with the movement functions used in \code{outbreaker.create.moves}. See 'details' for defaults.
#' 
#' @author Thibaut Jombart \email{t.jombart@@imperial.ac.uk}
#'
#' @param ... a named list (see details) of functions, each computing a log-posterior component.
#' @param posteriors a list of functions as returned by \code{outbreaker.create.posteriors}
#'
#' @details
#' The posterior functions used by default are:
#' \describe{
#' \item{genetic}{corresponds to the posterior distribution of the mutation rates; defaults to \code{\link{post.mu}}}
#' \item{all}{corresponds to the joint posterior; defaults to \code{\link{post.all}}}
#' }
#' 
#' 
#' @return a list of named functions
#' 
outbreaker.create.posteriors <- function(..., posteriors=NULL){
    ## PROCESS ... ONLY IF NO MOVES IS PASSED
    if(is.null(posteriors)){
        posteriors <- list(...)
    }
    
    ## SET DEFAULTS ##
    defaults <- list(genetic = post.genetic,
                     all = post.all
                     )

    ## MODIFY DEFAULTS WITH ARGUMENTS ##
    posteriors <- modify.defaults(defaults, posteriors, strict=FALSE)
    

    ## CHECK FUNCTIONS ##
    check.function.args <- function(f){
        args <- names(formals(f))
        if(!identical(sort(args), c("data", "param"))) {
            return(FALSE)
        } else {
            return(TRUE)
        }
    }

    args.checks <- sapply(unlist(posteriors), check.function.args)
    if(!all(args.checks)){
        culprits <- names(args.checks)[!args.checks]
        culprits <- paste(culprits, collapse=",")
        stop("problems in movements of: ", culprits, "\narguments shoud be: 'data', 'param'")
    }
    
    ## RETURN ##
    return(posteriors)
    
} # end outbreaker.create.posteriors
