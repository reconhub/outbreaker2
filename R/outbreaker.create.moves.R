#' Create the parameter space explorer
#'
#' This function creates a named list of functions, each of which moves parameters of outbreaker.
#' If these functions are provided as input, they will be used. Otherwise, default functions will be used.
#'
#' @author Thibaut Jombart \email{t.jombart@@imperial.ac.uk}
#'
#' @param ... a list of named movement functions; see details for available names.
#' @param moves a list of functions as returned by \code{outbreaker.create.moves}
#'
#' @return a list of named functions
#' 
outbreaker.create.moves <- function(..., moves=NULL){
    ## PROCESS ... ONLY IF NO MOVES IS PASSED
    if(is.null(moves)){
        moves <- list(...)
    }

    
    ## SET DEFAULTS ##
    defaults <- list(mu = list(move.mu),
                     t.inf = list(move.t.inf),
                     ances = list(move.ances, move.swap.ances),
                     )

    ## MODIFY DEFAULTS WITH ARGUMENTS ##
    moves <- modify.defaults(defaults, moves)
    

    
} # end outbreaker.create.moves
