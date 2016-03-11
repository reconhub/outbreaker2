
#' Internal functions for outbreaker2
#'
#' These functions are meant for internal use in outbreaker2. They are not exported, and their API might change in the future.
#'
#' \describe{
#' \item{modify.default}{modify default arguments using user-provided values}
#' \item{add.to.function.context}{add a set of objects to the environment of a function}
#' \item{check.param}{a function performing various checks of the status of the chain; tries to pick up impossible values for augmented data and parameters, or -Inf log-likelihood or posterior values; this function is used only in 'safemode' (see \code{\link{outbreaker.config}})}
#' }
#'
#' @rdname internals
#'
#' @param defaults a list containing default arguments
#' @param x in \code{modify.defaults}, a list with user-provided arguments; in \code{add.to.context}, a function or an environment.
#' @param strict a logical indicating if errors shoul be returned when 'x' contains items not in 'defaults'
#' @param param a list of parameters as returned by \code{outbreaker.mcmc.init}
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





#' @rdname internals
#'
## checks are only sure for the 'current' state
##
look.for.trouble <- function(param, data){
    ## PREPARE OUTPUT ##
    out <- list(pass=TRUE, msg=NULL)


    ## LIEKLIHOOD / POSTERIOR / PRIOR
    ## look for NAs in loglike / post / prior
    if(any(is.na(param$post))){
        out$pass <- FALSE
        out$msg <- c(out$msg, "NA detected in posterior values (param$post)")
    }
    if(any(is.na(param$like))){
        out$pass <- FALSE
        out$msg <- c(out$msg, "NA detected in likelihood values (param$like)")
    }
    if(any(is.na(param$prior))){
        out$pass <- FALSE
        out$msg <- c(out$msg, "NA detected in prior values (param$prior)")
    }

    ## look for NAs in loglike / post / prior
    if(!all(is.finite(param$post))){
        out$pass <- FALSE
        out$msg <- c(out$msg, "non-finite posterior values detected (param$post)")
    }
    if(!all(is.finite(param$like))){
        out$pass <- FALSE
        out$msg <- c(out$msg, "non-finite likelihood values detected (param$like)")
    }
    if(!all(is.finite(param$prior))){
        out$pass <- FALSE
        out$msg <- c(out$msg, "non-finite prior values detected (param$prior)")
    }


    ## CHECKS ON MU ##
    ## check that mu > 0
    if(param$current.mu<0){
        out$pass <- FALSE
        out$msg <- c(out$msg, "mu has a negative value:", param$current.mu)
    }

    ## check if mu is NA
    if(is.na(param$current.mu)){
        out$pass <- FALSE
        out$msg <- c(out$msg, "mu is NA")
    }

    ## check if mu is finite
    if(!is.finite(param$current.mu)){
        out$pass <- FALSE
        out$msg <- c(out$msg, "mu is not finite and equals:", param$current.mu)
    }

    ## check if mu is numeric
    if(!is.numeric(param$current.mu)){
        out$pass <- FALSE
        out$msg <- c(out$msg, "mu is not numeric and equals:", param$current.mu)
    }


    ## ANCESTRIES ##
    ## look for new imported cases (should not happen)
    if(!identical(is.na(param$ances[[1]]), is.na(param$current.ances))){
        out$pass <- FALSE
        out$msg <- c(out$msg, "imported cases have changed")
    }

    ## look for negative ancestries
    if(any(param$current.ances<1,na.rm=TRUE)){
       out$pass <- FALSE
       out$msg <- c(out$msg, "some ancestries point to unknown cases (param$current.ances<1)")
    }

    ## look for ancestries greater than 'N'
    if(any(param$current.ances>length(param$ances[[1]]),na.rm=TRUE)){
       out$pass <- FALSE
       out$msg <- c(out$msg, "some ancestries point to unknown cases (param$current.ances>N)")
    }

    ## case infecting itself
    if(any(param$current.ances==1:length(param$current.ances),na.rm=TRUE)){
       out$pass <- FALSE
       out$msg <- c(out$msg, "auto-infections detected (param$current.ances[i]==i)")
    }


    ## INFECTION DATES ##
    ## check NA
    if(any(is.na(param$current.t.inf))){
        out$pass <- FALSE
        out$msg <- c(out$msg, "NA detected in infection dates (param$current.t.inf)")
    }

    ## check finite values
    if(any(!is.finite(param$current.t.inf))){
        out$pass <- FALSE
        out$msg <- c(out$msg, "some infection dates are not finite (param$current.t.inf)")
    }

    ## check that values are numeric
    if(any(!is.numeric(param$current.t.inf))){
        out$pass <- FALSE
        out$msg <- c(out$msg, "some infection dates are not numeric (param$current.t.inf)")
    }

    ## check that delays between infections are > 0
    if(any((param$current.t.inf - param$current.t.inf[param$current.ances]) < 1, na.rm=TRUE)){
        out$pass <- FALSE
        out$msg <- c(out$msg, "some delays between succesive infections are less than 1 (param$current.t.inf)")
    }

    ## check that delays to collection are > 0
    if(any((data$dates-param$current.t.inf) < 1, na.rm=TRUE)){
        out$pass <- FALSE
        out$msg <- c(out$msg, "some delays to collection are less than 1 (param$current.t.inf)")
    }

    ## SHAPE OUTPUT AND RETURN ##
    out$msg <- paste(out$msg, collapse="\n")
    return(out)
} # end check.param








##############################
## NON-DOCUMENTED FUNCTIONS ##
##############################

## check which ancestries can move (returns a TRUE/FALSE vector)
can.move.ances <- function(param, config){
    out <- !is.na(param$current.ances) & # non-imported case
        (param$current.t.inf > min(param$current.t.inf)) & # not the first date
            config$move.ances # add user-specification through move.ances
    return(out)
}


## check which ancestries can move (returns a TRUE/FALSE vector)
can.be.swapped <- function(param, config){
    out <- !is.na(param$current.ances) & # non-imported case
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


## which cases are possible ancestors for a case 'i'
are.possible.ances <- function(t.inf, i){
    if(length(i)>1) stop("i has a length > 1")
    if(any(t.inf[i]==min(t.inf))) return(NA)
    return(which(t.inf < t.inf[i[1]]))
} # end are.possible.ances


## choose one possible ancestor for a case 'i'
choose.possible.ances <- function(t.inf, i){
    return(sample(are.possible.ances(t.inf=t.inf, i=i), 1))
} # end choose.possible.ances



## swaps ancestries in the tree
## x-> i becomes i->x
## plus all subsequent changes
swap.ances <- function(param, config, i){
    ## stop if 'i' out of range
    if(i>length(param$current.ances)) stop("trying to swap ancestry of case ",
                                           i, " while there are only ",
                                           length(param$current.ances), " cases")
    ## find cases for which ancestries can move
    id.ok.to.swap <- which(can.be.swapped(param, config))

    ## find ancestor of 'i'
    x <- param$current.ances[i]

    ## stop if case 'i' is imported - this should not happen
    if(is.na(x)){
        warning("trying to swap the ancestry of the imported case ", i)
        return(param)
    }

    ## check that x can be swapped, stop if not
    if(!(x %in% id.ok.to.swap)){
        return(param)
    }

    ## find indices to swap
    to.be.x <- intersect(which(param$current.ances==i), id.ok.to.swap)
    to.be.i <- intersect(which(param$current.ances==x), id.ok.to.swap)

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
