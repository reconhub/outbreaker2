#' Create the parameter space explorer
#'
#' This function creates a named list of functions, each of which moves parameters of outbreaker.
#' If these functions are provided as input, they will be used. Otherwise, default functions will be used.
#'
#' User-provided functions will be checked to ensure the right arguments are present. All movement functions should have to following arguments:
#' \describe{
#' \item{data}{a list of named items containing input data as returned by \code{\link{outbreaker.data}}}
#' \item{config}{a set of settings as returned by \code{\link{outbreaker.config}}}
#' \item{param}{a list of parameters as returned by \code{outbreaker.mcmc.init}}
#' \item{rand}{a list of items as returned by \code{outbreaker.rand.vec}}
#' }
#'
#' @author Thibaut Jombart \email{t.jombart@@imperial.ac.uk}
#'
#' @param ... a named list of movement functions for parameters or augmented data; see details for available movements and the corresponding names.
#' @param moves a list of functions as returned by \code{outbreaker.create.moves}
#' @param config a list of settings as returned by \code{outbreaker.config}
#'
#' @details
#' The movement functions which are used by default:
#' \describe{
#' \item{move.mu}{a function to move mutation rates}
#' \item{move.t.inf}{a function to move dates of infection}
#' \item{move.ances}{a function to move ancestries, i.e. the transmission tree}
#' \item{move.swap.ances}{another function to move ancestries, relying on swapping ancestries (a->b becomes b->a)}
#' }
#'
#'
#' @return a list of named functions
#'
outbreaker.create.moves <- function(..., moves=NULL, config=outbreaker.config()){
    ## PROCESS ... ONLY IF NO MOVES IS PASSED
    if(is.null(moves)){
        moves <- list(...)
    }


    ## SET DEFAULTS ##
    defaults <- list(move.mu = move.mu,
                     move.t.inf = move.t.inf,
                     move.ances = move.ances,
                     move.swap.ances = move.swap.ances
                     )

    ## MODIFY DEFAULTS WITH ARGUMENTS ##
    moves <- modify.defaults(defaults, moves, strict=FALSE)


    ## REMOVE FUNCTIONS IF MOVEMENTS DISABLED ##
    if(!any(config$move.ances)) moves$move.ances <- NULL
    if(!any(config$move.t.inf)) moves$move.t.inf <- NULL
    if(!any(config$move.mu)) moves$move.mu <- NULL


    ## CHECK FUNCTIONS ##
    check.function.args <- function(f){
        if(is.null(f)) return(TRUE)
        args <- names(formals(f))
        if(!identical(sort(args), c("config","data","param", "rand"))) {
            return(FALSE)
        } else {
            return(TRUE)
        }
    }

    args.checks <- sapply(unlist(moves), check.function.args)
    if(!all(args.checks)){
        culprits <- names(args.checks)[!args.checks]
        culprits <- paste(culprits, collapse=",")
        stop("problems in movements of: ", culprits, "\narguments shoud be: 'data', 'config', 'param', 'rand'")
    }

    ## RETURN ##
    return(moves)

} # end outbreaker.create.moves
