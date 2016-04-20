#' Create prior functions for outbreaker
#'
#' This function creates a named list of functions, each of which are compute log-prior components of outbreaker.
#' If these functions are provided as input, they will be used. Otherwise, default functions will be used.
#'
#' User-provided functions will be checked to ensure the right arguments are present. All log-prior functions should have to following arguments:
#' \describe{
#' \item{param}{a list of parameters as returned by \code{outbreaker.mcmc.init}} 
#' }
#'
#' Note that unlike \code{outbreaker.create.moves}, prior functions have no mandatory names. They only need to be compatible with the movement functions used in \code{outbreaker.create.moves}. See 'details' for defaults.
#' 
#' @author Thibaut Jombart \email{t.jombart@@imperial.ac.uk}
#'
#' @param ... a named list (see details) of functions, each computing a log-prior component.
#' @param priors a list of functions as returned by \code{outbreaker.create.priors}
#'
#' @details
#' The prior functions used by default are:
#' \describe{
#' \item{mu}{corresponds to the prior of the mutation rates; defaults to \code{\link{prior.mu}}}
#' \item{all}{corresponds to the joint priors of all parameters; defaults to \code{\link{prior.all}}}
#' }
#' 
#' 
#' @return a list of named functions
#' 
outbreaker.create.priors <- function(..., priors=NULL){
    ## PROCESS ... ONLY IF NO MOVES IS PASSED
    if(is.null(priors)){
        priors <- list(...)
    }
    
    ## SET DEFAULTS ##
    defaults <- list(mu = prior.mu,
                     all = prior.all
                     )

    ## MODIFY DEFAULTS WITH ARGUMENTS ##
    priors <- modify.defaults(defaults, priors, strict=FALSE)
    

    ## CHECK FUNCTIONS ##
    check.function.args <- function(f){
        args <- names(formals(f))
        if(!identical(sort(args), c("param"))) {
            return(FALSE)
        } else {
            return(TRUE)
        }
    }

    args.checks <- sapply(unlist(priors), check.function.args)
    if(!all(args.checks)){
        culprits <- names(args.checks)[!args.checks]
        culprits <- paste(culprits, collapse=",")
        stop("problems in movements of: ", culprits, "\narguments shoud be: 'param'")
    }
    
    ## RETURN ##
    return(priors)
    
} # end outbreaker.create.priors
