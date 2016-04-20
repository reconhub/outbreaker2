## #' Create likelihood functions for outbreaker
## #'
## #' This function creates a named list of functions, each of which are compute log-likelihood components of outbreaker.
## #' If these functions are provided as input, they will be used. Otherwise, default functions will be used.
## #'
## #' User-provided functions will be checked to ensure the right arguments are present. All log-likelihood functions should have to following arguments:
## #' \describe{
## #' \item{data}{a list of named items containing input data as returned by \code{\link{outbreaker.data}}}
## #' \item{param}{a list of parameters as returned by \code{outbreaker.mcmc.init}}
## #' }
## #'
## #' Note that unlike \code{outbreaker.create.moves}, likelihood functions have no mandatory names. They only need to be compatible with the movement functions used in \code{outbreaker.create.moves}. See 'details' for defaults.
## #'
## #' @author Thibaut Jombart \email{t.jombart@@imperial.ac.uk}
## #'
## #' @param ... a named list (see details) of functions, each computing a log-likelihood component.
## #' @param loglike a list of functions as returned by \code{outbreaker.create.loglike}
## #'
## #' @details
## #' The likelihood functions used by default are:
## #' \describe{
## #' \item{timing.infections}{corresponds to the probability of delays between cases (generation time)}
## #' \item{timing.sampling}{corresponds to the probability of delays between infection and reporting}
## #' \item{timing}{corresponds to the probability of all delays}
## #' \item{genetic}{corresponds to the probability of the genetic sequences given the transmission tree}
## #' \item{all}{corresponds to the overal probability of the data}
## #' }
## #'
## #'
## #' @return a list of named functions
## #'



##
## This function creates a named list of log-likelihood functions with enclosed data
##
create.loglike <- function(data){

    ## These are all the functions generating various log-likelihood functions;
    ## we list them by alphabetic order

    default.functions <- list(genetic = make.ll.genetic,
                              reporting = make.ll.reporting,
                              timing.infections = make.ll.timing.infections,
                              timing.sampling = make.ll.timing.sampling
                              )
    ll.function.names <- names(default.functions)
    out <- lapply(default.functions, function(f) f(data))


    ## We need a function summing all log-likelihoods - useful as a shortcut for several
    ## movements of parameters and augmented data.

    out$all <- function(param, i=NULL){
        sum(vapply(out[ll.function.names], function(f) f(param, i)), na.rm=TRUE)
    }


    ## We need a function computing likelihood relating to timing, which includes:
    ## - p(sampling dates | infections dates)
    ## - p(infection dates | ancestral infection dates)

    out$timing <- function(param, i=NULL){
        out$timing.infections(param, i) + out$timing.sampling(param, i)
    }


    return(out)
}
