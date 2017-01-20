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
## #' \item{param}{a list of parameters as returned by \code{outbreaker_create_mcmc}}
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



## This function creates a named list of log-likelihood functions with enclosed
## data

custom_likelihood <- function(...) {

    ## These are all the functions generating various log-likelihood functions;
    ## we list them by alphabetic order

    ## default_functions <- list(genetic = bind_ll_genetic,
    ##                           reporting = bind_ll_reporting,
    ##                           timing_infections = bind_ll_timing_infections,
    ##                           timing_sampling = bind_ll_timing_sampling,
    ##                           timing = bind_ll_timing,
    ##                           all = bind_ll_all
    ##                           )

    ## lapply(default_functions, function(f) f(data))
  
    ll_functions <- list(...)
    
    if (length(ll_functions) == 1L && is.list(ll_functions[[1]])) {
        ll_functions <- ll_functions[[1]]
    }
    

    defaults <- list(genetic = NULL,
                     reporting = NULL,
                     timing_infections = NULL,
                     timing_sampling = NULL
                     )

    ## out <-  modify_defaults(defaults, ll_functions, FALSE)

    out <- defaults

    return(out)
}
