#' Encloses argument in a function's environment
#'
#' This function takes a function \code{f} and a series of named arguments, and
#' returns a closure of \code{f} which will only rely on one single argument
#' 'param'. This is used to reduce the number of arguments passed around to
#' likelihood or movement functions. This is needed when using custom functions
#' via \code{\link{custom_priors}}, \code{\link{custom_likelihoods}}, or
#' \code{\link{custom_moves}}.
#' 
#' @param f The function to which arguments are bound.
#' 
#' @param ... Named arguments to bind to the function's environment.
#' 
#' @author Initial code by Rich FitzJohn (see 'references') with some
#' adaptations by Thibaut Jombart
#'
#' @references Initial code comes from the \code{partially_apply} function in
#' the 'rodeint' package (\url{/github.com/richfitz/rodeint/}).
#'
#' 
#' @export
#' 
bind_to_function <- function(f, ...) {

    ## We isolate the arguments of 'f' and identify those without defaults,
    ## which need to be provided through '...'. Arguments of 'f' which have a
    ## default value will be replaced with content of '...' if provided. The
    ## function returned is a closure with a single argument 'param', and with
    ## all non-default arguments in its environment.
    
    dots <- list(...)
    dots_names <- names(dots)
    f_args <- setdiff(names(formals(f)), "param")
    have_no_default <- vapply(formals(f)[f_args], is.symbol, logical(1))
    f_args_no_default <- names(have_no_default)[have_no_default]


    ## CHECKS ##
    if (is.primitive(f)) {
        stop("Cannot use with primitive functions")
    }

    ## Nothing to do if nothing provided
    if (length(dots) == 0) {
        if (length(f_args) > 1) {
            stop("'...' is empty but 'f' has more than one argument.")
        }
        return(f)
    }

    ## All objects passed in '...' need to be named
    if (is.null(dots_names) || !all(nzchar(dots_names))) {
        stop("All arguments provided through '...' need to be named.")
    }

    ## Name duplication is not allowed
    if (any(duplicated(dots_names))) {
        msg <- sprintf("Duplicated formal arguments: ",
             collapse(unique(dots_names[duplicated(dots_names)])))
        stop(msg)
    }

    ## ... cannot contain 'param'
    if ("param" %in% dots_names) {
        stop("'...' cannot contain an argument 'param'")
    }

    
    ## make sure all arguments of 'f' which don't have default values but
    ## 'param' are in '...'
    
    are_missing <- !f_args_no_default %in% dots_names
    if (any(are_missing)) {
        missing_args <- f_args_no_default[are_missing] 
        missing_args <- paste(missing_args, collapse = ", ")
        msg <- sprintf("Arguments of %s are missing from '...' and have no default: %s",
             deparse(substitute(f)),
             missing_args)
        stop(msg)
    }

    ## remove arguments that are not part of 'f'
    to_keep <- dots_names %in% f_args
    dots <- dots[to_keep]
    dots_names <- names(dots)
    

    ## Attach arguments to 'f'
    add_to_function_environment(f, dots)
}






## This function adds a list of objects to a function's environment

add_to_function_environment <- function(f, defaults) {
    e <- as.environment(defaults)
    parent.env(e) <- environment(f)
    ff <- formals(f)
    replace_formals(f, ff[c(setdiff(names(ff), names(defaults)))], e)
}






## This replaces forms, but preserves attributes except for srcref,
## which will be invalid for any nontrivial change (and will
## confusingly be printed with the wrong structure).

replace_formals <- function(fun, value, envir = environment(fun)) {
    old_attributes <- attributes(fun)
    formals(fun, envir = envir) <- value
    attributes(fun) <- old_attributes[names(old_attributes) != "srcref"]
    fun
}
