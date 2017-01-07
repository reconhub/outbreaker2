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
## #' \code{\link{outbreaker.data}}}
## #'
## #' \item{param}{a list of parameters as returned by \code{outbreaker.create.mcmc}}
## #' }
## #'
## #' Note that unlike \code{outbreaker.create.moves}, likelihood functions have no
## #' mandatory names. They only need to be compatible with the movement functions
## #' used in \code{outbreaker.create.moves}. See 'details' for defaults.
## #'
## #' @author Thibaut Jombart \email{t.jombart@@imperial.ac.uk}
## #'
## #' @param ... a named list (see details) of functions, each computing a log-likelihood component.
## #'
## #' @param loglike a list of functions as returned by \code{outbreaker.create.loglike}
## #'
## #' @details
## #' The likelihood functions used by default are:
## #' \describe{
## #'
## #' \item{timing.infections}{corresponds to the probability of delays between
## #' cases (generation time)}
## #'
## #' \item{timing.sampling}{corresponds to the probability of delays between
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

create.loglike <- function(data) {

    ## These are all the functions generating various log-likelihood functions;
    ## we list them by alphabetic order

    default.functions <- list(genetic = make.ll.genetic,
                              reporting = make.ll.reporting,
                              timing.infections = make.ll.timing.infections,
                              timing.sampling = make.ll.timing.sampling,
                              timing = make.ll.timing,
                              all = make.ll.all
                              )

    lapply(default.functions, function(f) f(data))

}
