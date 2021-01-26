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
#' A named list of list(function, arity) pairs with the class
#'     \code{custom_likelihood}, each function implementing a customised
#'     log-likelihood component of outbreaker. Functions which are not
#'     customised will result in a list(NULL, 0) component.
#'
#' @author Thibaut Jombart (\email{thibautjombart@@gmail.com})
#'
#' @param ... a named list of functions, each computing a log-likelihood component.
#'
#' @return a list of named functions
#'
#' @seealso See \href{http://www.repidemicsconsortium.org/outbreaker2/articles/customisation.html#customizing-likelihood}{customization vignette} for detailed examples on how to customize likelihoods.
#'
#' @export
#' 
#' @examples
#' 
#' ## specify a null model by disabling all likelihood components
#' f_null <- function(data, param) {
#'   return(0.0)
#' }
#' 
#' null_model <- custom_likelihoods(genetic = f_null,
#'                                 timing_sampling = f_null,
#'                                 timing_infections = f_null,
#'                                 reporting = f_null,
#'                                 contact = f_null)
#'
#' null_config <- list(find_import = FALSE,
#'                     n_iter = 200,
#'                     sample_every = 1)
#'
#' ## load data
#' x <- fake_outbreak
#' data <- outbreaker_data(dates = x$sample, dna = x$dna, w_dens = x$w)
#' 
#' res_null <- outbreaker(data = data,
#'                        config = null_config,
#'                        likelihoods = null_model)
#'
#' ## visualise ancestries to see if all transmission trees have been explored
#' plot(res_null, type = "alpha")


## USING CUSTOM LIKELIHOOD FUNCTIONS

## Likelihood functions in outbreaker2 are implemented using Rcpp. However,
## these functions can also be replaced by customized functions. These can be
## specified by the user, through the '...' argument of
## 'custom_likelihoods'. These functions must have at least 2 arguments:

## - data: a valid 'outbreaker_data' list

## - param: a list containing current parameter states, as returned by
## - create_param

## - [i=NULL]: (optional) a list of the cases for which the loglikelihoods
## - should be calculated. Needs to default to `NULL` in which case the
## - loglikelihood of the entire tree is calculated.

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
    list_function_or_null <- function(x) {
        is.list(x) && (is.null(x[[1]]) || is.function(x[[1]]))
    }

    # Ensure that custom_likelihoods(l) == custom_likelihoods(custom_likelihoods(l))
    is_list_function_or_null <- vapply(likelihoods, list_function_or_null, logical(1))
    is_function_or_null <- vapply(likelihoods, function_or_null, logical(1))

    if (!all(is_function_or_null) & !all(is_list_function_or_null)) {
        culprits <- likelihoods_names[!is_function_or_null]
        msg <- paste0("The following likelihoods are not functions: ",
                      paste(culprits, collapse = ", "))
        stop(msg)
    }

    # If the arity of the likelihood functions is three, the last argumen should
    # be the (1-based) indices of the cases we're currently perturbing. This
    # allows us to calculate the local likelihood delta, rather than having to
    # calculate the likelihood of the entire tree twice for every single
    # perturbation we make.
    if (!all(is_list_function_or_null)) {
        likelihoods <- lapply(likelihoods, function(x) { if (is.null(x)) return(list(x, 0)); list(x, length(methods::formalArgs(x))) })
    }

    arity_two_or_three <- function(x) {
        if (is.function(x[[1]])) {
            return (x[[2]] == 2L | x[[2]] == 3L)
        }
        return(T)
    }

    legal_arity <- vapply(likelihoods, arity_two_or_three, logical(1))

    if (!all(legal_arity)) {
        culprits <- likelihood_names[!legal_arity]
        msg <- paste0("The following likelihoods do not have arity two or three: ",
                      paste(culprits, collapse=", "))
        stop(msg)
    }

    names(likelihoods) <- likelihoods_names
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

