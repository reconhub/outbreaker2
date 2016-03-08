
#' Internal functions for outbreaker2
#'
#' These functions are meant for internal use in outbreaker2. They are not exported, and their API might change in the future.
#'
#' \describe{
#' \item{modify.default}{modify default arguments using user-provided values}
#' \item{add.to.function.context}{add a set of objects to the environment of a function}
#' }
#' 
#' @rdname internals
#' 
#' @param defaults a list containing default arguments
#' @param x in \code{modify.defaults}, a list with user-provided arguments; in \code{add.to.context}, a function.
#' @param strict a logical indicating if errors shoul be returned when 'x' contains items not in 'defaults'
#'
#' @author Rich Fitzjohn, Thibaut Jombart
#' 
#' @importFrom utils "modifyList"
#' 
modify.defaults <- function(defaults, x, strict=TRUE){
    extra <- setdiff(names(x), names(defaults))
    if (strict && (length(extra) > 0L)){
        stop("Additional invalid options: ", paste(extra, collapse=", "))
    }
    modifyList(defaults, x)
} # end modify.defaults



#' @rdname internals
#'
#' @param objects a list of objects to be added to the environment of 'x'; all items need to be named.
#' 
add.to.context <- function(x, objects){
    ## escape if object is empty ##
    if(length(objects) < 1L) return(x)

    ## check that all objects are named ##
    if(any(names(objects == ""))) {
        warning("all objects to be added to the environment of 'x' must be named; only adding named objects")
        to.keep <- names(objects) != ""
        objects <- objects[to.keep]
    }

    ## add objects to environment ##
    for(i in seq_along(objects)){
        assign(x = names(objects)[i],
               value = objects[[i]],
               envir = environment(x))
    }

    ## return ##
    return(x)
} # add.to.context
