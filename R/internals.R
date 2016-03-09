
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
#' @param x in \code{modify.defaults}, a list with user-provided arguments; in \code{add.to.context}, a function or an environment.
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
#' @examples
#' f1 <- function(x) { foo(x)+1}
#'
add.to.context <- function(x, objects){
    ## get environment
    if(is.environment(x)){
        env <- x
    } else{
        env <- environment(x)
    }

    ## escape if object is empty
    if(length(objects) < 1L) return(invisible(NULL))

    ## check that all objects are named
    if(any(names(objects) == "")) {
        warning("all objects to be added to the environment of 'x' must be named; only adding named objects")
        to.keep <- names(objects) != ""
        objects <- objects[to.keep]
    }

    ## add objects to environment
    for(i in seq_along(objects)){
        ## recursive behaviour if object is a list
        if(is.list(objects[[i]])){
            add.to.context(env, objects[[i]])
        }
        assign(x = names(objects)[i],
               value = objects[[i]],
               envir = env)
    }

    return(invisible(NULL))
} # add.to.context








## check which ancestries can move (returns a TRUE/FALSE vector)
can.move.ances <- function(param, config){
    out <- !is.na(param$current.ances) & # non-imported case
        (param$current.t.inf > min(param$current.t.inf)) & # not the first date
            config$move.ances # add user-specification through move.ances
    return(out)
}


#' @rdname internals
## random selection of cases for which ancestries is moved
select.ances.to.move <- function(param, config){
    choices <- which(can.move.ances(param, config))
    n.to.move <- max(round(config$prop.ances.move * length(choices)),0)
    out <- sample(choices, n.to.move, replace=FALSE)
    return(out)
}


## find possible ancestors for a case 'i'
find.possible.ances <- function(t.inf, i){
    if(length(i)>1) stop("i has a length > 1")
    if(any(t.inf[i]==min(t.inf))) return(NA)
    return(sample(which(t.inf < t.inf[i[1]]), 1))
} # end find.possible.ances



## swaps ancestries in the tree
## x-> i becomes i->x
## plus all subsequent changes
swap.ances <- function(param, i){
    ## stop if 'i' out of range
    if(i>length(param$current.ances)) stop("trying to swap ancestry of case ",
                                           i, " while there are only ",
                                           length(param$current.ances), " cases")

    ## find ancestor of 'i'
    x <- param$current.ances[i]

    ## stop if case 'i' is imported
    if(is.na(x)){
        warning("trying to swap the ancestry of the imported case ", i)
        return(param)
    }

    ## find indices to swap
    to.be.x <- which(param$current.ances==i)
    to.be.i <- which(param$current.ances==x)

    ## swap 'i' and 'x' in ancestries
    param$current.ances[to.be.x] <- x
    param$current.ances[to.be.i] <- i

    ## the ancestor of 'i' is now has the ancestor of 'x'
    param$current.ances[i] <- param$current.ances[x]

    ## 'i' is now the ancestor of 'x'
    param$current.ances[x] <- i

    ## swap t.inf
    param$current.t.inf[c(x,i)] <- param$current.t.inf[c(i,x)]

    ## return
    return(param)
} # end swap.ancestries
