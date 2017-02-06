#'
#' Customise likelihood functions for outbreaker
#'
#' This function is used to specify customised likelihood functions for
#' outbreaker. Custom functions are specified as a named list or series of
#' comma-separated, named arguments, indicating which log-likelihood component
#' they compute. Values currently available are:
#'
#' \itemize{
#'
#' \item \code{genetic}: the genetic likelihood; by default, the function
#' \code{cpp_ll_genetic} is used.
#'
#' \item \code{timing_sampling}: the likelihood of sampling times; by default, the function
#' \code{cpp_ll_timing_sampling} is used.
#'
#' \item \code{timing_infections}: the likelihood of infection times; by default, the function
#' \code{cpp_ll_timing_infections} is used.
#'
#' \item \code{reporting}: the likelihood of the reporting process; by default,
#' the function \code{cpp_ll_reporting} is used.
#'
#' \item \code{contact}: the likelihood of the contact tracing data; by default,
#' the function \code{cpp_ll_contact} is used.
#' }
#'
#' All log-likelihood functions should have the following arguments, in this
#' order:
#'
#' \itemize{
#'
#' \item \code{data}: a list of named items containing input data as returned by
#' \code{\link{outbreaker_data}}
#'
#' \item \code{param}: a list of parameters with the class
#' \code{\link{create_param}}
#'
#' }
#'
#'
#' @return
#' A named list of functions with the class \code{custom_likelihood}, each
#'     implementing a customised log-likelihood components of
#'     outbreaker. Functions which are not customised will result in a NULL
#'     component.
#'
#' @author Thibaut Jombart (\email{thibautjombart@@gmail.com})
#'
#' @param ... a named list of functions, each computing a log-likelihood component.
#'
#' @return a list of named functions
#'
#' @export
#'


## USING CUSTOM LIKELIHOOD FUNCTIONS

## Likelihood functions in outbreaker2 are implemented using Rcpp. However,
## these functions can also be replaced by customized functions. These can be
## specified by the user, through the '...' argument of
## 'custom_likelihoods'. These functions must have 2 arguments:

## - data: a valid 'outbreaker_data' list

## - param: a list containing current parameter states, as returned by
## - create_param

custom_likelihoods <- function(...) {

    ll_functions <- list(...)

    if (length(ll_functions) == 1L && is.list(ll_functions[[1]])) {
        ll_functions <- ll_functions[[1]]
    }


    defaults <- list(genetic = NULL,
                     reporting = NULL,
                     timing_infections = NULL,
                     timing_sampling = NULL,
                     contact = NULL
                     )

    likelihoods <-  modify_defaults(defaults, ll_functions, FALSE)
    likelihoods_names <- names(likelihoods)



    ## check all likelihoods are functions

    function_or_null <- function(x) {
        is.null(x) || is.function(x)
    }

    is_ok <- vapply(likelihoods, function_or_null, logical(1))

    if (!all(is_ok)) {
        culprits <- likelihoods_names[!is_ok]
        msg <- paste0("The following likelihoods are not functions: ",
                      paste(culprits, collapse = ", "))
        stop(msg)
    }


    ## check they all have a single argument

    with_two_args <- function(x) {
        if(is.function(x)) {
            return (length(methods::formalArgs(x)) == 2L)
        }

        return(TRUE)
    }

    two_args <- vapply(likelihoods, with_two_args, logical(1))

    if (!all(two_args)) {
        culprits <- likelihoods_names[!two_args]
        msg <- paste0("The following likelihoods dont' have two arguments: ",
                      paste(culprits, collapse = ", "))
        stop(msg)
    }


    class(likelihoods) <- c("custom_likelihoods", "list")
    return(likelihoods)

}







#' @rdname custom_likelihoods
#'
#' @export
#'
#' @aliases print.custom_likelihoods
#'
#' @param x an \code{outbreaker_config} object as returned by \code{create_config}.
#'

print.custom_likelihoods <- function(x, ...) {
    cat("\n\n ///// outbreaker custom likelihoods ///\n")
    cat("\nclass:", class(x))
    cat("\nnumber of items:", length(x), "\n\n")

    is_custom <- !vapply(x, is.null, FALSE)


    names_default <- names(x)[!is_custom]
    if (length(names_default) > 0) {
        cat("/// custom likelihoods set to NULL (default used) //\n")
        print(x[!is_custom])
    }


    names_custom <- names(x)[is_custom]
    if (length(names_custom) > 0) {
        cat("/// custom likelihoods //\n")
        print(x[is_custom])
    }

    return(invisible(NULL))

}

