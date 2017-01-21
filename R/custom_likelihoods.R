## #' Create likelihood functions for outbreaker
## #'
## #' This function creates a named list of functions, each of which are compute
## #' log-likelihood components of outbreaker.  If these functions are provided as
## #' input, they will be used. Otherwise, default functions will be used.
## #'
## #' User-provided functions will be checked to ensure the right arguments are
## #' present. All log-likelihood functions should have to following arguments:
## #'  
## #' \describe{
## #'
## #' \item{data}{a list of named items containing input data as returned by
## #' \code{\link{outbreaker_data}}}
## #'
## #' \item{param}{a list of parameters as returned by \code{create_mcmc}}
## #' }
## #'
## #' Note that unlike \code{outbreaker_custom_moves}, likelihood functions have no
## #' mandatory names. They only need to be compatible with the movement functions
## #' used in \code{outbreaker_custom_moves}. See 'details' for defaults.
## #'
## #' @author Thibaut Jombart \email{t_jombart@@imperial_ac_uk}
## #'
## #' @param ... a named list (see details) of functions, each computing a log-likelihood component.
## #'
## #' @param loglike a list of functions as returned by \code{outbreaker_custom_likelihood}
## #'
## #' @details
## #' The likelihood functions used by default are:
## #' \describe{
## #'
## #' \item{timing_infections}{corresponds to the probability of delays between
## #' cases (generation time)}
## #'
## #' \item{timing_sampling}{corresponds to the probability of delays between
## #' infection and reporting}
## #'
## #' \item{timing}{corresponds to the probability of all delays}
## #'
## #' \item{genetic}{corresponds to the probability of the genetic sequences given
## #' the transmission tree}
## #'
## #' \item{all}{corresponds to the overal probability of the data}
## #'
## #' }
## #'
## #'
## #' @return a list of named functions
## #'




## USING CUSTOM LIKELIHOOD FUNCTIONS

## Likelihood functions in outbreaker2 are implemented using Rcpp. However,
## these functions can also be replaced by customized functions. These can be
## specified by the user, through the '...' argument of
## 'custom_likelihoods'. These functions must have 2 arguments:

## - data: a valid 'outbreaker_data' list

## - param: a list containing current parameter states, as returned by create_mcmc



custom_likelihoods <- function(...) {
    
    ll_functions <- list(...)
    
    if (length(ll_functions) == 1L && is.list(ll_functions[[1]])) {
        ll_functions <- ll_functions[[1]]
    }
    

    defaults <- list(genetic = NULL,
                     reporting = NULL,
                     timing_infections = NULL,
                     timing_sampling = NULL
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
    

    class(likelihoods) <- c("outbreaker_likelihoods", "list")
    return(likelihoods)

}
